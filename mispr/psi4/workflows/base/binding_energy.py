"""Define the binding energy workflow, using psi4 instead of Gaussian.

Mirrors ``mispr.gaussian.workflows.base.binding_energy.get_binding_energies``.

* **Fireworks 1-4**: Optimize + frequency for each of the two molecules.
* **Firework 5**: Combine the two optimized molecules at ``index`` and optimize the
  resulting complex (``LinkedMolPsi4FW`` -- linking and optimizing must happen in
  the same Firework, since the linked molecule is only passed via
  ``fw_spec["prev_calc_molecule"]``, which does not survive a Firework boundary).
* **Firework 6**: Frequency calculation on the optimized complex.
* **Fireworks 6a/6b** (optional, ``counterpoise_correction=True``): Single-point
  energies of each monomer, frozen at its position within the optimized complex,
  with the other monomer present as ghost atoms (basis functions only) -- used
  for a Boys-Bernardi counterpoise (BSSE) correction.
* **Firework 7**: Compute the binding energy (and, if requested, the BSSE
  correction) and save it (``BindingEnergytoDB``, reused unmodified from the
  Gaussian firetasks -- pure energy bookkeeping, not tied to any QM engine).
"""

import os

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.firetasks.parse_outputs import BindingEnergytoDB
from mispr.psi4.fireworks.core import Psi4FW, LinkedMolPsi4FW
from mispr.psi4.firetasks.parse_outputs import CounterpoiseToDB
from mispr.psi4.workflows.base.core import common_fw, WORKFLOW_KWARGS

__author__ = "Ruiqi Luo"
__status__ = "Development"


def _to_psi4_solvent(solvent_gaussian_inputs):
    """
    Translate the Gaussian-style solvent string (e.g. "(Solvent=Water)") into the
    dict RunPsi4 expects (e.g. {"solvent": "water"}); returns None for gas phase.
    """
    if not solvent_gaussian_inputs:
        return None
    solvent_inputs = [
        i.lower() for i in solvent_gaussian_inputs.strip("()").split(",")
    ]
    solvent_name = next(
        (s.split("=")[1] for s in solvent_inputs if "solvent" in s), "water"
    )
    return {"solvent": solvent_name}


def get_binding_energies(
    mol_operation_type,
    mol,
    index,
    bond_order=1,
    db=None,
    name="binding_energy_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    counterpoise_correction=False,
    **kwargs
):
    """
    Define a workflow for calculating the binding energy between two molecules
    using psi4. See ``mispr.gaussian.workflows.base.binding_energy.get_binding_energies``
    for the full description of arguments; this is a psi4-backed counterpart.

    ``opt_gaussian_inputs``/``freq_gaussian_inputs`` use the same keys as the
    Gaussian workflow (e.g. {"functional": "B3LYP", "basis_set": "6-31G(d)",
    "route_parameters": {"Opt": None}}); Gaussian-only keys are accepted but
    ignored. ``solvent_gaussian_inputs``/``solvent_properties`` are kept only for
    ``BindingEnergytoDB``'s db metadata bookkeeping (same as the Gaussian
    workflow); they are translated internally into the ``solvent`` dict RunPsi4
    expects to actually drive the PCM calculation.

    counterpoise_correction (bool, optional): If ``True``, additionally compute
        a Boys-Bernardi counterpoise (BSSE) correction (two extra single-point
        calculations) and save it as "be_eV_cp_corrected" alongside the regular
        "be_eV"; see the module docstring for what the two numbers mean.
        Defaults to ``False`` (skip it, keeping the original, cheaper workflow).

    Returns:
        Workflow
    """
    working_dir = working_dir or os.getcwd()
    tag = kwargs.get("tag", "unknown")
    solvent = _to_psi4_solvent(solvent_gaussian_inputs)

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

    skips = skips or [None, None]

    # Fireworks 1-4: optimize + frequency for each of the two molecules
    processed_mol_1 = process_mol(mol_operation_type[0], mol[0], db=db)
    processed_mol_2 = process_mol(mol_operation_type[1], mol[1], db=db)

    _, label_1, fws_1 = common_fw(
        mol=processed_mol_1,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name="mol_1",
        gout_key="mol_1",
        skips=skips[0],
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        solvent=solvent,
        tag=tag,
    )

    _, label_2, fws_2 = common_fw(
        mol=processed_mol_2,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name="mol_2",
        gout_key="mol_2",
        skips=skips[1],
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        solvent=solvent,
        tag=tag,
    )

    # Firework 5: combine the two optimized molecules and optimize the complex
    # (linking and optimizing must be one Firework -- see module docstring)
    linked_opt_fw = LinkedMolPsi4FW(
        gout_keys=["mol_1", "mol_2"],
        index=index,
        bond_order=bond_order,
        db=db,
        name="linked_mol_optimization",
        parents=fws_1 + fws_2,
        working_dir=os.path.join(working_dir, "linked_mol", "Optimization"),
        gaussian_input_params=opt_gaussian_inputs,
        gout_key="mol_linked_opt",
        cart_coords=cart_coords,
        solvent=solvent,
        tag=tag,
    )

    # Firework 6: frequency calculation on the optimized complex
    linked_freq_fw = Psi4FW(
        prev_calc_key="mol_linked_opt",
        gaussian_input_params=freq_gaussian_inputs,
        db=db,
        name="linked_mol_frequency",
        parents=linked_opt_fw,
        working_dir=os.path.join(working_dir, "linked_mol", "Frequency"),
        gout_key="mol_linked",
        cart_coords=cart_coords,
        solvent=solvent,
        tag=tag,
    )

    fws = fws_1 + fws_2 + [linked_opt_fw, linked_freq_fw]

    # Fireworks 6a/6b (optional): counterpoise (BSSE) correction -- two extra
    # single-point energies at the complex's optimized geometry, each keeping
    # one monomer "real" and ghosting out the other
    cp_fws = []
    if counterpoise_correction:
        n1 = len(processed_mol_1)
        n2 = len(processed_mol_2)
        charge_1, mult_1 = processed_mol_1.charge, processed_mol_1.spin_multiplicity
        charge_2, mult_2 = processed_mol_2.charge, processed_mol_2.spin_multiplicity

        mono1_ghost_fw = Psi4FW(
            prev_calc_key="mol_linked_opt",
            gaussian_input_params={
                "functional": opt_gaussian_inputs["functional"],
                "basis_set": opt_gaussian_inputs["basis_set"],
            },
            charge=charge_1,
            multiplicity=mult_1,
            ghost_indices=list(range(n1, n1 + n2)),
            db=db,
            name="mono1_ghost",
            parents=linked_opt_fw,
            working_dir=os.path.join(working_dir, "linked_mol", "counterpoise"),
            gout_key="mono1_ghost",
            cart_coords=cart_coords,
            tag=tag,
        )
        mono2_ghost_fw = Psi4FW(
            prev_calc_key="mol_linked_opt",
            gaussian_input_params={
                "functional": opt_gaussian_inputs["functional"],
                "basis_set": opt_gaussian_inputs["basis_set"],
            },
            charge=charge_2,
            multiplicity=mult_2,
            ghost_indices=list(range(0, n1)),
            db=db,
            name="mono2_ghost",
            parents=linked_opt_fw,
            working_dir=os.path.join(working_dir, "linked_mol", "counterpoise"),
            gout_key="mono2_ghost",
            cart_coords=cart_coords,
            tag=tag,
        )
        cp_fw = Firework(
            CounterpoiseToDB(
                mol_linked_key="mol_linked_opt",
                mono1_ghost_key="mono1_ghost",
                mono2_ghost_key="mono2_ghost",
            ),
            parents=[mono1_ghost_fw, mono2_ghost_fw],
            name="counterpoise_correction",
            spec={"_launch_dir": os.path.join(working_dir, "linked_mol", "counterpoise")},
        )
        cp_fws = [mono1_ghost_fw, mono2_ghost_fw, cp_fw]
        fws += cp_fws

    # Firework 7: gather energies, compute binding energy, save to db
    fw_analysis = Firework(
        BindingEnergytoDB(
            index=index,
            db=db,
            keys=["mol_1", "mol_2", "mol_linked"],
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in BindingEnergytoDB.required_params + BindingEnergytoDB.optional_params
            },
        ),
        parents=fws[:],
        name="binding_energy_analysis",
        spec={"_launch_dir": os.path.join(working_dir, "analysis")},
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(name, "_".join([label_1, label_2])),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
    )
