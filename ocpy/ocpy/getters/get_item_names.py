"""
get_item_names.py

Retrieves the string name of every item in an OpenSim ComponentSet
(e.g. CoordinateSet, ForceSet, MuscleSet).

Uses the OpenSim Python bindings (opensim-core).
"""

import opensim as osim


def get_item_names(component_set) -> list[str]:
    """
    Return the names of all items in an OpenSim ComponentSet.

    Parameters
    ----------
    component_set : opensim.ComponentSet (or any Set with getSize / get)
        Any OpenSim set-type object.

    Returns
    -------
    list[str]
        Ordered list of item name strings.
    """
    num_items = component_set.getSize()
    item_names = []

    for i in range(num_items):
        item = component_set.get(i)
        item_names.append(str(item.getName()))

    return item_names
