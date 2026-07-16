"""Define the bond dissociation energy workflow, using psi4 instead of Gaussian.

Mirrors ``mispr.gaussian.workflows.base.bde.get_bde``. The bond-breaking/fragment
enumeration (``BreakMolFW``) and the final energy bookkeeping (``BDEtoDB``) are pure
cheminformatics/arithmetic that don't care which QM engine produced the underlying
energies, so both are reused as-is; only the optimization/frequency Fireworks are
psi4-specific (via ``mispr.psi4.workflows.base.core.common_fw``).
"""

import os

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.utilities.metadata import get_job_name, get_mol_formula
from mispr.gaussian.firetasks.parse_outputs import BDEtoDB
from mispr.psi4.fireworks.break_mol import BreakMolFW
from mispr.psi4.workflows.base.core import common_fw, WORKFLOW_KWARGS

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_8"
__version__ = "0.0.5"

def get_bde(
    mol_operation_type,
    mol,
    ref_charge=0,
    fragment_charges=None,
    bonds=None,
    open_rings=False,
    db=None,
    name="bde_calculation_psi4",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    skips=False,
    visualize=True,
    **kwargs,
):
    """
    Define a dynamic workflow for calculating the bond dissociation energy using psi4:
    optimize + frequency the principle molecule, break it into fragments (trying
    several charge splits per bond), optimize + frequency each fragment, then
    compute the BDE for every bond/charge-split combination. Mirrors
    ``mispr.gaussian.workflows.base.bde.get_bde``, which has the same overall
    structure and arguments; see that function's docstring for more background on
    the fragment charge-enumeration behavior.

    ``opt_gaussian_inputs``/``freq_gaussian_inputs`` use the same keys as the Gaussian
    workflow (e.g. {"functional": "B3LYP", "basis_set": "6-31G(d)",
    "route_parameters": {"Opt": None}}); Gaussian-only keys such as
    "link0_parameters" are accepted but ignored.

    Args:
        mol_operation_type (str): Type of molecule operation; see ``process_mol``
            in ``mispr/gaussian/utilities/mol.py`` for supported operations.
        mol (Molecule, GaussianOutput, str, dict): Source of the principle
            molecule to be processed; should match ``mol_operation_type``.
        ref_charge (int, optional): Charge on the principle molecule; defaults
            to 0.
        fragment_charges (list, optional): Additional charges to try for each
            fragment, beyond the defaults derived from ``ref_charge`` (0, +1,
            -1 for a neutral principle molecule); see
            ``mispr.gaussian.firetasks.geo_transformation.BreakMolecule`` for
            the full charge-enumeration rules. Defaults to ``None``.
        bonds (list, optional): List of bonds (as atom index pairs) to break;
            if ``None``, every bond found in the molecule is broken in turn.
        open_rings (bool, optional): Whether to open rings instead of skipping
            ring bonds that can't be simply split; defaults to ``False``.
        db (str or dict, optional): Database credentials; path to db.json or a
            dict; if ``None``, read from the configuration files.
        name (str, optional): Name of the workflow; defaults to
            "bde_calculation_psi4".
        working_dir (str, optional): Working directory for input/output files;
            defaults to the current working directory.
        opt_gaussian_inputs (dict, optional): Parameters for the optimization
            step; defaults to B3LYP/6-31G(d).
        freq_gaussian_inputs (dict, optional): Parameters for the frequency
            step; defaults to B3LYP/6-31G(d).
        skips (bool or list, optional): If ``True``, skips both "opt" and "freq"
            for the principle molecule; can also be a list like ``["opt"]`` to
            skip only specific jobs. Defaults to ``False`` (skip nothing).
        visualize (bool, optional): Whether ``BDEtoDB`` should generate a
            visualization of the computed BDEs; defaults to ``True``.
        kwargs (keyword arguments): Additional kwargs passed to
            ``BreakMolFW``/``BDEtoDB``/``Workflow`` (e.g. ``tag``).

    Returns:
        Workflow
    """
    fws = []
    working_dir = working_dir or os.getcwd()
    gout_key = "ref_mol"

    opt_gaussian_inputs = opt_gaussian_inputs or {
        "functional": "b3lyp",
        "basis_set": "6-31g(d)",
        "route_parameters": {"Opt": None},
    }
    freq_gaussian_inputs = freq_gaussian_inputs or {
        "functional": "b3lyp",
        "basis_set": "6-31g(d)",
        "route_parameters": {"Freq": None},
    }

    processed_mol = process_mol(mol_operation_type, mol, db=db, charge=ref_charge)
    label = get_mol_formula(processed_mol)

    if skips:
        skips = ["opt", "freq"] if skips is True else skips
    else:
        skips = None

    _, _, opt_freq_fws = common_fw(
        mol=processed_mol,
        working_dir=working_dir,
        dir_structure=["principle_mol"],
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name=label,
        skips=skips,
        gout_key=gout_key,
        tag=kwargs.get("tag", "unknown"),
    )
    fws += opt_freq_fws

    break_fw = BreakMolFW(
        mol=gout_key,
        mol_operation_type="get_from_run_dict",
        from_fw_spec=True,
        bonds=bonds,
        open_rings=open_rings,
        ref_charge=ref_charge,
        fragment_charges=fragment_charges,
        db=db,
        calc_frags=True,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        name=get_job_name(label, "breaking"),
        parents=fws[:],
        working_dir=os.path.join(working_dir, label, "fragments"),
        **kwargs,
    )
    fws.append(break_fw)

    fw_analysis = Firework(
        BDEtoDB(
            principle_mol_key=gout_key,
            db=db,
            visualize=visualize,
            **{
                i: j
                for i, j in kwargs.items()
                if i in BDEtoDB.required_params + BDEtoDB.optional_params
            },
        ),
        parents=fws[:],
        name="{}-{}".format(label, "bde_analysis"),
        spec={
            "_launch_dir": os.path.join(working_dir, label, "analysis"),
            "_allow_fizzled_parents": True,
        },
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(label, name),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
    )
