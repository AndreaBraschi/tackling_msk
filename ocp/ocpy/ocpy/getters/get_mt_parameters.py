"""
get_mt_parameters.py

Returns the muscle-tendon parameters for a list of muscles from an
OpenSim model.

Original author: Antoine Falisse (12/19/2018)

Uses the OpenSim Python bindings (opensim-core).
"""

import numpy as np
import opensim as osim


def get_mt_parameters(model: osim.Model, muscle_names: list[str]) -> np.ndarray:
    """
    Retrieve muscle-tendon parameters from an OpenSim model.

    Parameters
    ----------
    model : opensim.Model
        Loaded (and optionally initialised) OpenSim model.
    muscle_names : list[str]
        Ordered list of muscle names to retrieve parameters for.

    Returns
    -------
    params : np.ndarray, shape (5, num_muscles)
        Row 0: max isometric force         (FMo)
        Row 1: optimal fibre length        (lMo)
        Row 2: tendon slack length         (lTs)
        Row 3: pennation angle at lMo      (alphao)  [radians]
        Row 4: max contraction velocity    (vMmax = vMaxRel * lMo)
    """
    num_muscles = len(muscle_names)
    params = np.zeros((5, num_muscles))

    muscles = model.getMuscles()

    for i, name in enumerate(muscle_names):
        muscle = muscles.get(name)
        params[0, i] = muscle.getMaxIsometricForce()
        params[1, i] = muscle.getOptimalFiberLength()
        params[2, i] = muscle.getTendonSlackLength()
        params[3, i] = muscle.getPennationAngleAtOptimalFiberLength()
        params[4, i] = muscle.getMaxContractionVelocity() * params[1, i]

    return params
