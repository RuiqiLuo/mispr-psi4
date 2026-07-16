"""Define the electrostatic partial charge (ESP/RESP) workflow, using psi4
instead of Gaussian.

Mirrors ``mispr.gaussian.workflows.base.esp.get_esp_charges``: optimize +
frequency the molecule, run an ESP single-point + RESP charge fit on the
optimized geometry, then save the result.
"""

import os


from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.firetasks.parse_outputs import ESPtoDB
from mispr.psi4.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.psi4.firetasks.run_calc import ESP
from fireworks import Firework, Workflow

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_13"
__version__ = "0.0.5"

def get_esp_charges(
    mol_operation_type,
    mol,
    db=None,
    name="esp_charges_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    esp_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    **kwargs,
):
    """
    Define a workflow for calculating the electrostatic partial charges.

    * **Firework 1**: Optimize the molecule.
    * **Firework 2**: Run a frequency analysis.
    * **Firework 3**: Run an ESP calculation.
    * **Firework 4**: Create ESP document/json file.

    Args:
        mol_operation_type (str): The type of molecule operation. See ``process_mol``
            defined in ``mispr/gaussian/utilities/mol.py`` for supported operations.
        mol (Molecule, GaussianOutput, str, dict): Source of the molecule to be
            processed. Should match the ``mol_operation_type``.
        db (str or dict, optional): Database credentials; could be provided as the path
            to the "db.json" file or in the form of a dictionary; if none is provided,
            attempts to get it from the configuration files.
        name (str, optional): Name of the workflow. Defaults to "esp_charges_calculation".
        working_dir (str, optional): Path of the working directory where any required
            input files can be found and output will be created. Defaults to the current
            working directory.
        opt_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the optimization step; e.g.:

            .. code-block:: python

                {
                    "functional": "B3LYP",
                    "basis_set": "6-31G(d)",
                    "route_parameters": {"Opt": None},
                    "link0_parameters": {
                        "%chk": "checkpoint.chk",
                        "%mem": "45GB",
                        "%NProcShared": "24"}
                }

            The above default parameters will be used if not specified.
        freq_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the frequency step; default parameters will be used if not specified.
        esp_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the ESP step; default parameters will be used if not specified.
        solvent_gaussian_inputs (str, optional): Gaussian input parameters corresponding
            to the implicit solvent model to be used in the ESP calculations, if any;
            e.g.:

            .. code-block:: python

                "(Solvent=TetraHydroFuran)"

            These parameters should only be specified here and not included in the main
            gaussian_inputs dictionary for each job (i.e. ``opt_gaussian_inputs``,
            ``freq_gaussian_inputs``, etc.). Defaults to None.
        solvent_properties (dict, optional): Additional input parameters to be used in
            the ESP calculations and relevant to the solvent model, if any; e.g.,
            {"EPS":12}. Defaults to None.
        cart_coords (bool, optional): Uses cartesian coordinates in writing Gaussian
            input files if set to ``True``, otherwise uses z-matrix. Defaults to ``True``.
        oxidation_states (dict, optional): Dictionary of oxidation states that can be
            used in setting the charge and spin multiplicity of the molecule; e.g.:
            {"Li":1, "O":-2}. Defaults to None.
        skips (list, optional): List of jobs to skip; e.g.: ["opt", "freq"]; defaults
            to None.
        kwargs (keyword arguments): Additional kwargs to be passed to the workflow.

    Returns:
        tuple:
            - Workflow
            - str: Label of the molecule (e.g. "H2O", "water", etc.).
    """

    # Firework 1: Optimize the molecule
    # Firework 2: Run a frequency analysis
    working_dir = working_dir or os.getcwd()
    # transforming to get a moledule
    processed_mol = process_mol(mol_operation_type, mol, db=db)

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

    solvent = None
    tag = kwargs.get("tag", "unknown")

    _, label, fws = common_fw(
        mol=processed_mol,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name="mol",
        gout_key="mol",
        skips=skips,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        solvent=solvent,
        tag=tag,
    )

    # Firework 3: Run an ESP calculation
    # spec must carry "tag" itself -- unlike the opt/freq Fireworks built by
    # common_fw (which get it via Psi4FW's spec.update({"tag": tag, ...})), a plain
    # Firework doesn't inherit spec from its parents, so ESPtoDB would otherwise
    # KeyError on fw_spec["tag"] downstream
    esp_fw = Firework(
        tasks=[ESP(prev_calc_key="mol", gout_key="mol_esp", db=db)],
        parents=fws[:],
        name=f"{label}_esp",
        spec={"tag": tag},
    )
    fws.append(esp_fw)

    # Firework 4: Create ESP document/json file
    fw_analysis = Firework(
        ESPtoDB(
            db=db,
            keys=["mol", "mol_esp"],
            **{
                i: j
                for i, j in kwargs.items()
                if i in ESPtoDB.required_params + ESPtoDB.optional_params
            },
        ),
        parents=fws[:],
        name=f"{label}_esp_analysis",
        spec={"tag": tag},
    )
    fws.append(fw_analysis)

    return (
        Workflow(
            fws,
            name=name,
            **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
        ),
        label,
    )

