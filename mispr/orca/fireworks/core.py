"""Define common Fireworks used in ORCA workflows."""

import os
import logging

from fireworks import Firework

from mispr.orca.firetasks.run_calc import RunOrca

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


class OrcaFW(Firework):
    """
    Run a single ORCA calculation (single point, optimization, or frequency
    analysis, depending on ``gaussian_input_params["route_parameters"]``) and
    store the result.

    Reuses the "gaussian_input_params" naming convention from the Gaussian
    workflows (e.g. {"functional": "B3LYP", "basis_set": "6-31G(d)",
    "route_parameters": {"Opt": None}}) so that the same input dictionaries used
    to configure Gaussian jobs can be passed to an ORCA job unmodified
    (irrelevant Gaussian-only keys, e.g. "link0_parameters", are simply
    ignored). Note that basis set names must be ones ORCA recognizes.
    """

    def __init__(
        self,
        molecule=None,
        prev_calc_key=None,
        db=None,
        name="orca_calc",
        parents=None,
        working_dir=None,
        gaussian_input_params=None,
        tag="unknown",
        gout_key=None,
        **kwargs,
    ):
        """
        Args:
            molecule (Molecule, optional): pymatgen Molecule to run the
                calculation on; required unless ``prev_calc_key`` is given.
            prev_calc_key (str, optional): Key of a previous run in
                fw_spec["gaussian_output"] to use as the input structure (e.g.
                to chain a frequency calculation onto a prior optimization).
            db (str or dict, optional): Database credentials.
            name (str, optional): Name of the Firework.
            parents (Firework or [Firework], optional): Parent FWs this FW
                depends on.
            working_dir (str, optional): Working directory for the calculation.
            gaussian_input_params (dict, optional): Dictionary of parameters
                (see class docstring) controlling the calculation.
            tag (str, optional): Tag stored in the db documents for easy
                retrieval.
            gout_key (str, optional): Key to store this run under in
                fw_spec["gaussian_output"].
            kwargs: other kwargs passed to Firework.__init__ and RunOrca (e.g.
                ``orca_cmd``, ``num_cores``, ``memory``).
        """
        working_dir = working_dir or os.getcwd()
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        gaussian_input_params = gaussian_input_params or {}
        task_kwargs = {
            i: j
            for i, j in {**gaussian_input_params, **kwargs}.items()
            if i in RunOrca.optional_params
        }

        task = RunOrca(
            molecule=molecule,
            prev_calc_key=prev_calc_key,
            gout_key=gout_key,
            db=db,
            **task_kwargs,
        )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(OrcaFW, self).__init__(
            [task],
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS},
        )
