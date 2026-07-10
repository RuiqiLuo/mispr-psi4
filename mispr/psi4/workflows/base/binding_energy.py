
'''
Firework 1 & 2: Perform optimization on two molecules, respectively
Firework 3 & 4: do frequency calculation for the two molecules
Firework 5: combine the two molecules and optimize them
Firework 6: do the frequency calculation after optimizing
Firework 7: gather all energy which from the three processes. Get the binding energy, and save them in the dataset.
'''

__author__ = "Ruiqi Luo"
__Date__ = "2026_7_10"
__version__ = "0.0.5"

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
    **kwargs
):