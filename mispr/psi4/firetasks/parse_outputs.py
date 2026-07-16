"""Define psi4-specific result-processing firetasks.

Unlike the firetasks in mispr.gaussian.firetasks.parse_outputs (pure energy
bookkeeping, reused unmodified by the psi4 workflows since they don't care which
engine produced the underlying numbers), the counterpoise (BSSE) correction below
depends on psi4's ghost-atom support (mispr.psi4.firetasks.run_calc.RunPsi4's
ghost_indices), so it lives here instead.
"""

import logging

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.gaussian.utilities.misc import pass_gout_dict

__author__ = "Ruiqi Luo"
__status__ = "Development"

logger = logging.getLogger(__name__)

HARTREE_TO_EV = 27.2114


@explicit_serialize
class CounterpoiseToDB(FiretaskBase):
    """
    Compute the Boys-Bernardi counterpoise (BSSE) correction for a binding energy
    and add it to fw_spec, for BindingEnergytoDB (or a similar downstream analysis
    task) to pick up.

    The correction is the standard 2-term counterpoise-corrected interaction
    energy, evaluated at the (already-optimized) complex geometry:

        be_eV_cp_corrected = E(complex) - E(mol_1, ghost mol_2 present)
                                         - E(mol_2, ghost mol_1 present)

    Note this uses frozen, complex-geometry monomers (via ghost-atom single point
    calculations), not the separately-optimized isolated monomer geometries used
    by the "raw" binding energy (mispr.gaussian.firetasks.parse_outputs.
    BindingEnergytoDB's be_eV) -- the two numbers are complementary, not
    interchangeable: be_eV includes monomer relaxation energy but has some BSSE
    contamination; be_eV_cp_corrected is BSSE-free but excludes monomer
    relaxation energy. Both are reported.

    Args:
        mol_linked_key (str): Key of the optimized complex's run in
            fw_spec["gaussian_output"].
        mono1_ghost_key (str): Key of the "mol_1 real, mol_2 ghost" single-point
            run (at the complex geometry) in fw_spec["gaussian_output"].
        mono2_ghost_key (str): Key of the "mol_2 real, mol_1 ghost" single-point
            run (at the complex geometry) in fw_spec["gaussian_output"].
    """

    required_params = ["mol_linked_key", "mono1_ghost_key", "mono2_ghost_key"]
    optional_params = []

    def run_task(self, fw_spec):
        """Compute be_eV_cp_corrected from the three referenced runs and pass it
        forward via update_spec for BindingEnergytoDB to pick up."""
        complex_gout = pass_gout_dict(fw_spec, self["mol_linked_key"])
        mono1_ghost_gout = pass_gout_dict(fw_spec, self["mono1_ghost_key"])
        mono2_ghost_gout = pass_gout_dict(fw_spec, self["mono2_ghost_key"])

        e_complex = complex_gout["output"]["output"]["final_energy"]
        e_mono1_ghost = mono1_ghost_gout["output"]["output"]["final_energy"]
        e_mono2_ghost = mono2_ghost_gout["output"]["output"]["final_energy"]

        be_ev_cp_corrected = (
            e_complex - e_mono1_ghost - e_mono2_ghost
        ) * HARTREE_TO_EV

        logger.info(
            f"counterpoise correction complete: be_eV_cp_corrected = "
            f"{be_ev_cp_corrected}"
        )

        return FWAction(
            update_spec={"be_eV_cp_corrected": be_ev_cp_corrected},
            propagate=True,
        )
