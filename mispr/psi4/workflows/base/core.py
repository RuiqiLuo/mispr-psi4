"""Define a list of common Fireworks used in psi4 workflows."""

import os
import logging

from fireworks import Workflow

from mispr.psi4.fireworks.core import Psi4FW

__author__ = "Ruiqi"
__status__ = "Development"

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = Workflow.__init__.__code__.co_varnames


def common_fw(
    mol,
    working_dir,
    opt_gaussian_inputs,
    freq_gaussian_inputs,
    gout_key=None,
    db=None,
    mol_name=None,
    dir_head=None,
    dir_structure=None,
    skips=None,
    tag="unknown",
    **kwargs,
):
    """
    Define a list of Fireworks commonly used in psi4 workflows: optimize a molecule,
    then run a frequency analysis on the optimized structure. Mirrors
    ``mispr.gaussian.workflows.base.core.common_fw``, minus the Gaussian-specific
    molecule-processing/labeling machinery (callers are expected to already have a
    concrete pymatgen ``Molecule``).

    Args:
        mol (Molecule): pymatgen Molecule to run the calculations on.
        working_dir (str): Path of the working directory.
        opt_gaussian_inputs (dict): Parameters for the optimization step (see
            ``mispr.psi4.fireworks.core.Psi4FW``).
        freq_gaussian_inputs (dict): Parameters for the frequency step.
        gout_key (str, optional): Unique key for the run dict; defaults to "mol".
        db (str or dict, optional): Database credentials.
        mol_name (str, optional): Name of the molecule, used for labeling.
        dir_head (str, optional): Head of the directory where the workflow will run;
            defaults to ``mol_name``.
        dir_structure (list, optional): Additional subfolders appended under
            ``dir_head``.
        skips (list, optional): List of jobs to skip; e.g. ["opt"] or ["freq"].
        tag (str, optional): Tag stored in the db documents.
        kwargs: Additional kwargs; unrecognized ones are ignored (accepted for
            interface parity with the Gaussian ``common_fw``).

    Returns:
        Molecule, str, list: The input molecule, a label, and a list of Fireworks.
    """
    fws = []
    if not gout_key:
        gout_key = "mol"

    label = mol_name or "mol"
    dir_head = dir_head or label
    dir_structure = dir_structure or []
    working_dir = os.path.join(working_dir, dir_head, *dir_structure)

    skips = skips or []

    if not skips:
        opt_fw = Psi4FW(
            molecule=mol,
            gaussian_input_params=opt_gaussian_inputs,
            db=db,
            name=f"{label}_optimization",
            working_dir=os.path.join(working_dir, "Optimization"),
            gout_key=gout_key + "_opt",
            tag=tag,
        )
        fws.append(opt_fw)

        freq_fw = Psi4FW(
            prev_calc_key=gout_key + "_opt",
            gaussian_input_params=freq_gaussian_inputs,
            db=db,
            name=f"{label}_frequency",
            parents=opt_fw,
            working_dir=os.path.join(working_dir, "Frequency"),
            gout_key=gout_key,
            tag=tag,
        )
        fws.append(freq_fw)

    elif len(skips) == 1 and skips[0].lower() == "opt":
        freq_fw = Psi4FW(
            molecule=mol,
            gaussian_input_params=freq_gaussian_inputs,
            db=db,
            name=f"{label}_frequency",
            working_dir=os.path.join(working_dir, "Frequency"),
            gout_key=gout_key,
            tag=tag,
        )
        fws.append(freq_fw)

    elif len(skips) == 1 and skips[0].lower() == "freq":
        opt_fw = Psi4FW(
            molecule=mol,
            gaussian_input_params=opt_gaussian_inputs,
            db=db,
            name=f"{label}_optimization",
            working_dir=os.path.join(working_dir, "Optimization"),
            gout_key=gout_key + "_opt",
            tag=tag,
        )
        fws.append(opt_fw)

    return mol, label, fws
