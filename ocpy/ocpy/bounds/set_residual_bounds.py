import numpy as np
from typing import Dict

def set_residual_bounds(cfg) -> Dict[str, Dict[str, np.ndarray]]:
    # This function assign the user defined residual forces.


    bounds: Dict[str, Dict[str, np.ndarray]] = {}


    residuals = cfg["bounds"]["residuals"]
    keys = residuals.keys()

    for i, key in enumerate(list(keys)):
        items = residuals[key]

        residual_value: float = items[0]
        num_coords: int = items[1]

        lower = -residual_value * np.ones(num_coords)
        upper =  residual_value * np.ones(num_coords)
        scaling = np.maximum(np.abs(lower), np.abs(upper))

        bounds[key] = {
            "lower": lower / scaling,
            "upper": upper / scaling,
        }

    return bounds