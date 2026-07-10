"""
Minimal example: run a BDE (bond dissociation energy) workflow with the psi4
backend, end to end, using your real local config (config/my_launchpad.yaml and
config/db.json).

Usage:
    python run_bde_test.py

After it finishes, check the result yourself with mongosh:
    mongosh mongodb://localhost:27017/gaussian
    > db.bde.find({tag: "learning_test"}).pretty()

And clean it up when you're done:
    mongosh mongodb://localhost:27017/gaussian
    > db.bde.deleteMany({tag: "learning_test"})
    > db.molecules.deleteMany({inchi: "InChI=1S/H2O2/c1-2/h1-2H"})
"""

import os

from pymatgen.core.structure import Molecule
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from mispr.psi4.workflows.base.bde import get_bde

RUN_DIR = os.path.join(os.path.dirname(__file__), "run_output")
os.makedirs(RUN_DIR, exist_ok=True)
os.chdir(RUN_DIR)

# a simple test molecule: hydrogen peroxide
h2o2 = Molecule(
    ["O", "O", "H", "H"],
    [
        [0.0000, 0.7375, -0.0528],
        [0.0000, -0.7375, -0.0528],
        [0.8190, 0.8170, 0.4220],
        [-0.8190, -0.8170, 0.4220],
    ],
)

# cheap method/basis, just to keep this example fast
opt_params = {"functional": "hf", "basis_set": "sto-3g", "route_parameters": {"Opt": None}}
freq_params = {"functional": "hf", "basis_set": "sto-3g", "route_parameters": {"Freq": None}}

wf = get_bde(
    mol_operation_type="get_from_mol",
    mol=h2o2,
    working_dir=RUN_DIR,
    opt_gaussian_inputs=opt_params,
    freq_gaussian_inputs=freq_params,
    visualize=False,
    save_to_db=True,
    tag="learning_test",
)

lp = LaunchPad.from_file(
    os.path.join(os.path.dirname(__file__), "..", "config", "my_launchpad.yaml")
)
lp.add_wf(wf)

print("Workflow submitted. Running...")
rapidfire(lp, m_dir=RUN_DIR)

print("Done. Check MongoDB (see docstring at the top of this file for the commands).")
