"""Define firetasks for running psi4 calculations.

Unlike Gaussian, which runs as an external process and requires a separate
"write input file" / "run program" / "parse output file" split (mispr.gaussian.
firetasks.write_inputs.WriteInput, run_calc.RunGaussianCustodian,
parse_outputs.ProcessRun), psi4 is driven entirely through its Python API and keeps
its results in memory. This module therefore combines all three steps into a single
Firetask, RunPsi4, which runs the calculation and immediately builds a result
dictionary with the same schema as mispr.gaussian.utilities.gout.process_run (see
mispr.gaussian.utilities.dbdoc._cleanup_gout for the reference schema) and stores it
under fw_spec["gaussian_output"], the same fw_spec key used by the Gaussian firetasks.
This lets any downstream Firetask written for Gaussian runs (e.g. BDEtoDB, ESPtoDB,
IPEAtoDB, BindingEnergytoDB) consume psi4 results unmodified.
"""

import os
import json
import logging

from timeit import default_timer as timer

import numpy as np
import psi4

from pymatgen.core.structure import Molecule

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.gaussian.defaults import JOB_TYPES
from mispr.gaussian.utilities.misc import recursive_signature_remove
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.gaussian.utilities.db_utilities import get_db

__author__ = "Ruiqi"
__status__ = "Development"

logger = logging.getLogger(__name__)

DEFAULT_KEY = "gout_key"

# default psi4 job parameters used if not specified by the caller
DEFAULT_FUNCTIONAL = "b3lyp"
DEFAULT_BASIS_SET = "6-31g(d)"
DEFAULT_MEMORY = "4 GB"
DEFAULT_NUM_THREADS = 4

def _json_default(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# default PCM block used when a solvent is requested but no explicit pcm_string
# is provided
DEFAULT_PCM_TEMPLATE = """
    Units = Angstrom
    Medium {{
    SolverType = CPCM
    Solvent = {solvent}
    }}
    Cavity {{
    RadiiSet = bondi
    Type = GePol
    Scaling = True
    Area = 0.3
    Mode = Implicit
    }}
"""


GAS_CONSTANT = 8.314462618  # J/(mol K)
AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23  # J/K
PLANCK = 6.62607015e-34  # J s
AMU_TO_KG = 1.66053906660e-27
HARTREE_TO_J = 4.3597447222071e-18
STANDARD_T = 298.15  # K
STANDARD_P = 101325.0  # Pa


def _atomic_thermo_corrections(mass_amu, multiplicity, t=STANDARD_T, p=STANDARD_P):
    """
    Analytic ideal-gas thermochemistry corrections for a single free atom.

    psi4's frequency/Hessian driver is not well-behaved for one-atom systems (no
    vibrational degrees of freedom), so single-atom fragments -- which are common in
    BDE calculations, since breaking any X-H bond produces a lone H atom -- are
    handled by computing the (purely translational + electronic) corrections
    analytically instead of calling psi4.frequencies().

    Returns:
        dict: "Zero-point correction", "Enthalpy", and "Gibbs Free Energy"
            corrections, in Hartree, using the same additive convention as
            mispr.gaussian (i.e. to be added to the raw electronic energy).
    """
    mass_kg = mass_amu * AMU_TO_KG
    q_trans = (2 * np.pi * mass_kg * BOLTZMANN * t / PLANCK**2) ** 1.5 * (
        BOLTZMANN * t / p
    )
    s_trans = GAS_CONSTANT * (np.log(q_trans) + 2.5)  # J/(mol K)
    s_elec = GAS_CONSTANT * np.log(multiplicity)  # J/(mol K)
    u_corr = 1.5 * GAS_CONSTANT * t  # J/mol; no ZPE, no rotational dof
    h_corr = u_corr + GAS_CONSTANT * t  # J/mol
    g_corr = h_corr - t * (s_trans + s_elec)  # J/mol

    def to_hartree(x_per_mol):
        return x_per_mol / AVOGADRO / HARTREE_TO_J

    return {
        "Zero-point correction": 0.0,
        "Enthalpy": to_hartree(h_corr),
        "Gibbs Free Energy": to_hartree(g_corr),
    }


def _mol_to_psi4_geometry(mol, charge, multiplicity):
    """Build a psi4.core.Molecule from a pymatgen Molecule."""
    lines = [f"{charge} {multiplicity}"]
    for site in mol:
        lines.append(f"{site.specie.symbol} {site.x:.10f} {site.y:.10f} {site.z:.10f}")
    lines.append("units angstrom")
    lines.append("symmetry c1")
    lines.append("no_reorient")
    lines.append("no_com")
    return psi4.geometry("\n".join(lines))


def _psi4_geometry_to_mol(psi4_mol):
    """Convert an (already updated) psi4.core.Molecule back to a pymatgen Molecule."""
    psi4_mol.update_geometry()
    bohr_to_angstrom = 0.52917721067
    coords = psi4_mol.geometry().np * bohr_to_angstrom
    species = [psi4_mol.symbol(i) for i in range(psi4_mol.natom())]
    return Molecule(species, coords)


@explicit_serialize
class RunPsi4(FiretaskBase):
    """
    Run a psi4 calculation (single point, geometry optimization, and/or frequency
    analysis) and store the result under fw_spec["gaussian_output"] using the same
    dictionary schema mispr uses for Gaussian runs.

    Args:
        molecule (Molecule, optional): pymatgen Molecule to run the calculation on;
            required unless ``prev_calc_key`` is provided.
        prev_calc_key (str, optional): Key of a previous psi4/Gaussian-compatible run
            in fw_spec["gaussian_output"] whose output molecule should be used as the
            input structure for this calculation (mirrors chaining an opt firework
            into a freq firework); ignored if ``molecule`` is provided.
        charge (int, optional): Molecular charge; defaults to the charge on
            ``molecule`` if it carries one, otherwise 0.
        multiplicity (int, optional): Spin multiplicity; defaults to the value on
            ``molecule`` if it carries one, otherwise 1.
        functional (str, optional): DFT functional (or "hf"/post-HF method name);
            defaults to "b3lyp".
        basis_set (str, optional): Basis set; defaults to "6-31g(d)".
        route_parameters (dict, optional): Gaussian-style route parameters used only
            to decide which psi4 driver call to make; recognized keys (case
            insensitive) are "opt" (calls psi4.optimize), "freq" (calls
            psi4.frequencies, implies "opt" was already done on the input geometry),
            and anything else falls back to a single point energy
            (psi4.energy). Defaults to a single point calculation.
        solvent (dict, optional): Implicit solvent (PCM) options; supported keys are
            "solvent" (a psi4 PCM solvent name, e.g. "water") and "pcm_string" (a full
            custom PCM input block, overrides "solvent" if given). If not provided, the
            calculation is run in the gas phase.
        memory (str, optional): psi4 memory setting; defaults to "4 GB".
        num_threads (int, optional): Number of threads for psi4; defaults to 4.
        db (str or dict, optional): Database credentials; path to db.json or a dict.
        save_to_db (bool, optional): Whether to insert the run into the runs collection.
        save_to_file (bool, optional): Whether to save the run to a json file.
        filename (str, optional): Name of the json file to save the run to; defaults
            to "run.json".
        gout_key (str, optional): Key to store this run under in
            fw_spec["gaussian_output"]; the run is always additionally stored under
            the default key "gout_key" (as with Gaussian runs).
    """

    required_params = []
    optional_params = [
        "molecule",
        "prev_calc_key",
        "charge",
        "multiplicity",
        "functional",
        "basis_set",
        "route_parameters",
        "solvent",
        "memory",
        "num_threads",
        "db",
        "save_to_db",
        "save_to_file",
        "filename",
        "gout_key",
        "tag",
    ]

    def _get_molecule(self, fw_spec):
        mol = self.get("molecule")
        if mol is not None:
            if not isinstance(mol, Molecule):
                raise TypeError("molecule must be a pymatgen Molecule object")
            return mol

        prev_calc_key = self.get("prev_calc_key")
        if prev_calc_key:
            gout_dict = fw_spec.get("gaussian_output", {}).get(prev_calc_key)
            if not gout_dict:
                raise KeyError(
                    f"No previous run found under fw_spec['gaussian_output']"
                    f"['{prev_calc_key}']"
                )
            return Molecule.from_dict(gout_dict["output"]["output"]["molecule"])

        raise KeyError(
            "No molecule present; provide 'molecule' or 'prev_calc_key', or check "
            "fw_spec"
        )

    def run_task(self, fw_spec):
        working_dir = os.getcwd()
        mol = self._get_molecule(fw_spec)

        charge = self.get("charge", getattr(mol, "charge", 0) or 0)
        multiplicity = self.get(
            "multiplicity", getattr(mol, "spin_multiplicity", 1) or 1
        )

        functional = self.get("functional", DEFAULT_FUNCTIONAL)
        basis_set = self.get("basis_set", DEFAULT_BASIS_SET)
        route_parameters = self.get("route_parameters") or {}
        job_types = sorted(
            k.lower() for k in route_parameters if k.lower() in JOB_TYPES
        )
        if not job_types:
            job_types = ["sp"]

        solvent = self.get("solvent")

        psi4.core.clean()
        psi4.set_memory(self.get("memory", DEFAULT_MEMORY))
        psi4.set_num_threads(self.get("num_threads", DEFAULT_NUM_THREADS))
        psi4.core.set_output_file(
            os.path.join(working_dir, "psi4_output.dat"), False
        )

        psi4_mol = _mol_to_psi4_geometry(mol, charge, multiplicity)

        # psi4 does not auto-select an open-shell reference from the molecule's
        # multiplicity; radicals/odd-electron fragments (the norm after breaking a
        # bond) need an explicit UHF/UKS reference, otherwise psi4 defaults to RHF
        # and errors out ("RHF reference is only for singlets")
        is_open_shell = multiplicity != 1
        if is_open_shell:
            is_dft = functional.lower() not in ("hf", "scf")
            psi4.set_options({"reference": "uks" if is_dft else "uhf"})
        else:
            psi4.set_options({"reference": "rks" if functional.lower() not in ("hf", "scf") else "rhf"})

        is_pcm = bool(solvent)
        if is_pcm:
            psi4.set_options({"pcm": True, "pcm_scf_type": "total"})
            pcm_string = solvent.get("pcm_string") or DEFAULT_PCM_TEMPLATE.format(
                solvent=solvent.get("solvent", "water")
            )
            psi4.pcm_helper(pcm_string)
        else:
            psi4.set_options({"pcm": False})

        method = f"{functional}/{basis_set}"

        st = timer()
        corrections = {}
        has_completed = True
        error_msg = None
        try:
            if "opt" in job_types:
                energy, wfn = psi4.optimize(method, molecule=psi4_mol, return_wfn=True)
                final_mol = _psi4_geometry_to_mol(wfn.molecule())
            elif "freq" in job_types and mol.num_sites == 1:
                # a single atom has no vibrational degrees of freedom; psi4's
                # Hessian/frequency driver is not safe to call on 1-atom systems
                # (has been observed to segfault), so compute the analytic
                # ideal-gas atomic thermochemistry instead
                energy, wfn = psi4.energy(method, molecule=psi4_mol, return_wfn=True)
                final_mol = _psi4_geometry_to_mol(wfn.molecule())
                corrections.update(
                    _atomic_thermo_corrections(mol.species[0].atomic_mass, multiplicity)
                )
            elif "freq" in job_types:
                energy, wfn = psi4.frequencies(
                    method, molecule=psi4_mol, return_wfn=True
                )
                final_mol = _psi4_geometry_to_mol(wfn.molecule())
                # ZPVE is reported by psi4 as a correction on its own, while
                # ENTHALPY/GIBBS FREE ENERGY are reported as totals (electronic
                # energy + correction); mispr expects the correction alone in all
                # three cases, to be added back to final_energy downstream (see
                # BDEtoDB/IPEAtoDB)
                corrections["Zero-point correction"] = float(psi4.variable("ZPVE"))
                for label, qcvar in [
                    ("Enthalpy", "ENTHALPY"),
                    ("Gibbs Free Energy", "GIBBS FREE ENERGY"),
                ]:
                    total = psi4.variable(qcvar)
                    corrections[label] = float(total) - float(energy)
            else:
                energy, wfn = psi4.energy(method, molecule=psi4_mol, return_wfn=True)
                final_mol = _psi4_geometry_to_mol(wfn.molecule())

            dipole = None
            try:
                dipole = list(psi4.variable("CURRENT DIPOLE"))
            except Exception:
                pass

        except Exception as e:
            has_completed = False
            error_msg = str(e)
            energy = None
            final_mol = mol
            dipole = None
        run_time = timer() - st
        fw_spec["run_time"] = run_time

        output_block = {
            "final_energy": energy,
            "molecule": final_mol.as_dict(),
        }
        if corrections:
            output_block["corrections"] = corrections
        if dipole:
            output_block["dipole_moment"] = dipole
        if error_msg:
            output_block["error_message"] = error_msg

        gout_dict = {
            "input": {
                "functional": functional,
                "basis_set": basis_set,
                "route_parameters": route_parameters,
                "charge": charge,
                "spin_multiplicity": multiplicity,
                "molecule": mol.as_dict(),
            },
            "output": {
                "output": output_block,
                "has_gaussian_completed": has_completed,
                "is_pcm": is_pcm,
            },
            "functional": functional,
            "basis": basis_set,
            "phase": "solution" if is_pcm else "gas",
            "type": ";".join(job_types),
            **get_chem_schema(final_mol),
            "gauss_version": f"psi4-{psi4.__version__}",
        }
        gout_dict = {
            i: j
            for i, j in gout_dict.items()
            if i not in ["sites", "@module", "@class", "charge", "spin_multiplicity"]
        }
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]
        gout_dict["wall_time (s)"] = run_time
        gout_dict = json.loads(json.dumps(gout_dict, default=_json_default))
        gout_dict = recursive_signature_remove(gout_dict)

        if not has_completed:
            raise ValueError(f"psi4 did not complete normally: {error_msg}")

        run_list = {}
        db = self.get("db")
        if self.get("save_to_db"):
            runs_db = get_db(db)
            run_id = runs_db.insert_run(gout_dict)
            run_list["run_id_list"] = run_id
            logger.info("Saved parsed psi4 output to db")

        if self.get("save_to_file"):
            filename = self.get("filename", "run")
            file_path = os.path.join(working_dir, f"{filename}.json")
            with open(file_path, "w") as f:
                f.write(json.dumps(gout_dict, default=str))
            run_list["run_loc_list"] = file_path
            logger.info("Saved parsed psi4 output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
        mod_dict = {"_set": set_dict}
        if run_list:
            mod_dict.update({"_push": run_list})
        return FWAction(mod_spec=mod_dict, propagate=True)
