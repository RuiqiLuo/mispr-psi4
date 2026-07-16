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
import resp

from copy import deepcopy

from pymatgen.core.structure import Molecule, IMolecule

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.gaussian.defaults import JOB_TYPES
from mispr.gaussian.utilities.misc import recursive_signature_remove
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.gaussian.utilities.db_utilities import get_db

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_8"
__version__ = "0.0.5"

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


def _mol_to_zmatrix_block(mol):
    """
    Build a psi4-compatible Z-matrix block (atoms + internal coordinates, numeric
    values inlined) from a pymatgen Molecule.

    pymatgen's GaussianInput.get_zmatrix() already does the cartesian -> internal
    coordinate conversion, but emits Gaussian's "named variable" style (e.g.
    "H 1 B1" with "B1=0.957776" defined separately below); psi4's Z-matrix parser
    expects the numeric values inline instead, so the two are merged here.
    """
    # local import to avoid a module-level dependency for the (default) cartesian
    # path, which is the common case
    from pymatgen.io.gaussian import GaussianInput

    zmatrix_text = GaussianInput(mol, charge=0, spin_multiplicity=1).get_zmatrix()
    coord_lines, _, variable_lines = zmatrix_text.partition("\n\n")

    variables = {}
    for line in variable_lines.strip().splitlines():
        name, value = line.split("=")
        variables[name.strip()] = value.strip()

    resolved_lines = []
    for line in coord_lines.strip().splitlines():
        tokens = line.split()
        resolved_lines.append(
            " ".join(variables.get(tok, tok) for tok in tokens)
        )
    return resolved_lines


def _mol_to_psi4_geometry(mol, charge, multiplicity, cart_coords=True, ghost_indices=None):
    """
    Build a psi4.core.Molecule from a pymatgen Molecule.

    ghost_indices (set of int, optional): site indices to write as ghost atoms
        (prefixed "Gh(...)" -- basis functions present, but no nucleus/electrons),
        used for counterpoise (BSSE) correction; only supported with cart_coords,
        since ghost atoms in an internal-coordinate Z-matrix are not a common/
        well-defined combination.
    """
    ghost_indices = ghost_indices or set()
    if ghost_indices and not cart_coords:
        raise ValueError("ghost_indices is only supported with cart_coords=True")
    lines = [f"{charge} {multiplicity}"]
    if cart_coords:
        for i, site in enumerate(mol):
            symbol = site.specie.symbol
            if i in ghost_indices:
                symbol = f"Gh({symbol})"
            lines.append(f"{symbol} {site.x:.10f} {site.y:.10f} {site.z:.10f}")
    else:
        lines += _mol_to_zmatrix_block(mol)
    lines.append("units angstrom")
    lines.append("symmetry c1")
    lines.append("no_reorient")
    lines.append("no_com")
    return psi4.geometry("\n".join(lines))


def _psi4_geometry_to_mol(psi4_mol, charge, multiplicity, has_ghost=False):
    """
    Convert an (already updated) psi4.core.Molecule back to a pymatgen Molecule.

    charge/multiplicity must be passed explicitly and set via set_charge_and_spin:
    a plain Molecule(species, coords) call defaults to charge 0, with a
    multiplicity guessed from the (charge-0) electron count -- silently wrong for
    anything but a neutral species, and this molecule is what downstream steps
    (e.g. a chained frequency calculation reading it back via prev_calc_key) use to
    determine what charge/multiplicity to run at.

    has_ghost (bool, optional): if True, skip set_charge_and_spin's electron-count
        consistency check -- pymatgen's Molecule has no concept of "ghost" atoms
        (basis functions only, no real nucleus/electrons), so it would count a
        ghost atom's full atomic number towards the electron total and (correctly,
        from its point of view, but wrongly for this case) reject the requested
        charge/multiplicity as inconsistent. The resulting molecule is only used
        for bookkeeping in a counterpoise-correction calculation, not further
        chemistry, so an unchecked charge/multiplicity is fine here.
    """
    psi4_mol.update_geometry()
    bohr_to_angstrom = 0.52917721067
    coords = psi4_mol.geometry().np * bohr_to_angstrom
    # psi4 returns element symbols in all caps for multi-letter elements (e.g. "LI"
    # instead of "Li"); pymatgen requires standard capitalization to recognize them
    species = [psi4_mol.symbol(i).capitalize() for i in range(psi4_mol.natom())]
    result = Molecule(species, coords)
    if has_ghost:
        result._charge = charge
        result._spin_multiplicity = multiplicity
    else:
        result.set_charge_and_spin(charge, multiplicity)
    return result


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
        "oxidation_states",
        "cart_coords",
        "ghost_indices",
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

    def _charge_from_oxidation_states(self, mol):
        """
        Calculate the charge of a molecule/cluster from the oxidation state of its
        individual elements (e.g. {"Li": 1, "O": -2}); mirrors
        mispr.gaussian.firetasks.write_inputs.WriteInput._update_charge, since this
        is plain pymatgen bookkeeping with no dependency on the QM engine used.
        """
        mol_copy = deepcopy(mol)
        mol_copy.add_oxidation_state_by_element(self["oxidation_states"])
        mol_copy.set_charge_and_spin(super(IMolecule, mol_copy).charge)
        return int(mol_copy.charge)

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

        # set by mispr.gaussian.firetasks.geo_transformation.ProcessMoleculeInput
        # (e.g. when it ran with operation_type="link_molecules" to combine two
        # previously-computed molecules into one, ahead of this Firetask in the
        # same Firework); reused unmodified since it's plain molecule bookkeeping,
        # not tied to any QM engine
        if fw_spec.get("prev_calc_molecule"):
            return fw_spec["prev_calc_molecule"]

        raise KeyError(
            "No molecule present; provide 'molecule' or 'prev_calc_key', or check "
            "fw_spec"
        )

    def run_task(self, fw_spec):
        working_dir = os.getcwd()
        mol = self._get_molecule(fw_spec)

        if self.get("oxidation_states"):
            charge = self._charge_from_oxidation_states(mol)
        else:
            charge = self.get("charge", getattr(mol, "charge", 0) or 0)

        if self.get("multiplicity") is not None:
            multiplicity = self.get("multiplicity")
        else:
            # the default multiplicity must be recomputed for whatever "charge"
            # ends up being (e.g. when it comes from oxidation_states, which can
            # differ from mol's own charge) -- mol's own spin_multiplicity
            # attribute reflects mol's original charge, not this one, and using
            # it here would silently pair a charge with an inconsistent
            # multiplicity (caught downstream by psi4 as a validation error)
            mol_for_mult = deepcopy(mol)
            mol_for_mult.set_charge_and_spin(charge)
            multiplicity = mol_for_mult.spin_multiplicity

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

        ghost_indices = self.get("ghost_indices")
        psi4_mol = _mol_to_psi4_geometry(
            mol,
            charge,
            multiplicity,
            cart_coords=self.get("cart_coords", True),
            ghost_indices=set(ghost_indices) if ghost_indices else None,
        )

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
            if "opt" in job_types and mol.num_sites == 1:
                # a single atom has no internal coordinates to optimize (no
                # bond/angle for OPTKING to vary), and psi4.optimize() crashes
                # trying anyway ("not enough values to unpack") -- it's already
                # at its one and only possible geometry, so just take the
                # energy there (same bare-nucleus check as the single-atom
                # "freq" branch below, since a monomer here could equally be one)
                n_electrons = mol.species[0].Z - charge
                if n_electrons <= 0:
                    energy = 0.0
                    final_mol = mol
                else:
                    energy, wfn = psi4.energy(method, molecule=psi4_mol, return_wfn=True)
                    final_mol = _psi4_geometry_to_mol(wfn.molecule(), charge, multiplicity, has_ghost=bool(ghost_indices))
            elif "opt" in job_types:
                energy, wfn = psi4.optimize(method, molecule=psi4_mol, return_wfn=True)
                final_mol = _psi4_geometry_to_mol(wfn.molecule(), charge, multiplicity, has_ghost=bool(ghost_indices))
            elif "freq" in job_types and mol.num_sites == 1:
                # a single atom has no vibrational degrees of freedom; psi4's
                # Hessian/frequency driver is not safe to call on 1-atom systems
                # (has been observed to segfault), so compute the analytic
                # ideal-gas atomic thermochemistry instead
                n_electrons = mol.species[0].Z - charge
                if n_electrons <= 0:
                    # a bare nucleus (e.g. H+, which BDE's charge-state enumeration
                    # assigns to a single-atom fragment when it gives that fragment
                    # both bonding electrons) has no electronic structure for psi4 to
                    # solve; RHF/UHF divide by the electron count internally and
                    # crash with "float division by zero" if asked to try, so skip
                    # psi4 entirely -- there is no electronic energy to compute
                    energy = 0.0
                    final_mol = mol
                    corrections.update(
                        {
                            "Zero-point correction": 0.0,
                            "Enthalpy": 0.0,
                            "Gibbs Free Energy": 0.0,
                        }
                    )
                else:
                    energy, wfn = psi4.energy(method, molecule=psi4_mol, return_wfn=True)
                    final_mol = _psi4_geometry_to_mol(wfn.molecule(), charge, multiplicity, has_ghost=bool(ghost_indices))
                    corrections.update(
                        _atomic_thermo_corrections(mol.species[0].atomic_mass, multiplicity)
                    )
            elif "freq" in job_types:
                energy, wfn = psi4.frequencies(
                    method, molecule=psi4_mol, return_wfn=True
                )
                final_mol = _psi4_geometry_to_mol(wfn.molecule(), charge, multiplicity, has_ghost=bool(ghost_indices))
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
                final_mol = _psi4_geometry_to_mol(wfn.molecule(), charge, multiplicity, has_ghost=bool(ghost_indices))

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
            # get_chem_schema builds a SMILES/InChI/formula schema by treating every
            # site as a real atom; that's meaningless (and errors out on the
            # charge/electron-count check) for a ghost-atom counterpoise-correction
            # calculation, where some "atoms" contribute basis functions only, no
            # real nucleus/electrons -- skip it in that case, since only the energy
            # is actually needed downstream
            **(get_chem_schema(final_mol) if not ghost_indices else {}),
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
    
@explicit_serialize
class ESP(FiretaskBase):

    optional_params = [
        "molecule", "prev_calc_key",
        "charge", "multiplicity", "oxidation_states",
        "cart_coords",
        "method_esp", "basis_esp",
        "resp_options",
        "memory", "num_threads",
        "db", "save_to_db", "save_to_file", "filename",
        "gout_key", "tag",
    ]

    def run_task(self, fw_spec):
        working_dir = os.getcwd()
        mol = self.get("molecule")
        if mol is not None:
            pass
        elif self.get("prev_calc_key"):
            prev_calc_key = self.get("prev_calc_key")
            mol = Molecule.from_dict(
                fw_spec.get("gaussian_output", {}).get(prev_calc_key)["output"]["output"]["molecule"]
            )
        else:
            raise KeyError("Don't have 'molecule' and 'prev_calc_key' ")

        charge = self.get("charge", getattr(mol, "charge", 0) or 0)
        multiplicity = self.get(
            "multiplicity", getattr(mol, "spin_multiplicity", 1) or 1
        )

        psi4.core.clean()
        psi4.set_memory(self.get("memory", DEFAULT_MEMORY))
        psi4.set_num_threads(self.get("num_threads", DEFAULT_NUM_THREADS))
        psi4.core.set_output_file(
            os.path.join(working_dir, "psi4_output.dat"), False
        )

        psi4_mol = _mol_to_psi4_geometry(
            mol, charge, multiplicity, cart_coords=self.get("cart_coords", True)
        )

        method_esp = self.get("method_esp", "scf")
        basis_esp = self.get("basis_esp", "6-31g*")

        default_resp_options = {
            "VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0],
            "VDW_POINT_DENSITY": 1.0,
            "RESP_A": 0.0005,
            "RESP_B": 0.1,
            "METHOD_ESP": method_esp,
            "BASIS_ESP": basis_esp,
        }
        # ** means inserting a list into another list
        # if there don't have any input files, it will use the default_resp_options
        resp_options = {**default_resp_options, **self.get("resp_options", {})}

        has_completed = True
        error_msg = None
        try:
            energy = psi4.energy(
                name=f"{method_esp}/{basis_esp}", molecule=psi4_mol, return_wfn=False
            )  # return_wfn = False means only return one value
            resp_charges = resp.resp([psi4_mol], resp_options)
        except Exception as e:
            has_completed = False
            error_msg = str(e)
            energy = None
            resp_charges = None

        output_block = {
            "final_energy": energy,
            "molecule": mol.as_dict(),
        }
        if resp_charges is not None:
            output_block["ESP_charges"] = resp_charges[0].tolist()
            output_block["RESP_charges"] = resp_charges[1].tolist()
        if error_msg:
            output_block["error_message"] = error_msg

        gout_dict = {
            "input": {
                "functional": method_esp,
                "basis_set": basis_esp,
                "charge": charge,
                "spin_multiplicity": multiplicity,
                "molecule": mol.as_dict(),
            },
            "output": {
                "output": output_block,
                "has_gaussian_completed": has_completed,
                "is_pcm": False,
            },
            "functional": method_esp,
            "basis": basis_esp,
            "phase": "gas",
            "type": "esp",
            **get_chem_schema(mol),
            "gauss_version": f"psi4-{psi4.__version__}",
        }
        gout_dict = {
            i: j
            for i, j in gout_dict.items()
            if i not in ["sites", "@module", "@class", "charge", "spin_multiplicity"]
        }
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]
        gout_dict = json.loads(json.dumps(gout_dict, default=_json_default))
        gout_dict = recursive_signature_remove(gout_dict)

        if not has_completed:
            raise ValueError(f"psi4 ESP calculation did not complete normally: {error_msg}")

        run_list = {}
        db = self.get("db")
        if self.get("save_to_db"):
            runs_db = get_db(db)
            run_id = runs_db.insert_run(gout_dict)
            run_list["run_id_list"] = run_id
            logger.info("Saved parsed psi4 ESP output to db")

        if self.get("save_to_file"):
            filename = self.get("filename", "run")
            file_path = os.path.join(os.getcwd(), f"{filename}.json")
            with open(file_path, "w") as f:
                f.write(json.dumps(gout_dict, default=str))
            run_list["run_loc_list"] = file_path
            logger.info("Saved parsed psi4 ESP output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
        mod_dict = {"_set": set_dict}
        if run_list:
            mod_dict.update({"_push": run_list})
        return FWAction(mod_spec=mod_dict, propagate=True)

        
