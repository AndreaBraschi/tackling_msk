import numpy as np
from numpy import array, ndarray


def find_vector_between_two_points(p0: ndarray, p1: ndarray) -> ndarray:
    return p1 - p0

def norm(vector: ndarray) -> ndarray:

    arr_quad = np.power(vector, 2)
    arr_sum_quad = arr_quad.sum()
    return np.sqrt(arr_sum_quad)


def unit_vector(vector: ndarray) -> ndarray:
    vector_norm = norm(vector)

    def divide(vector_dimension, vector_norm):
        return vector_dimension / vector_norm

    vfunc = np.vectorize(divide)
    res = vfunc([vector[0], vector[1], vector[2]], vector_norm)

    return res

def cross_product(vec1: ndarray, vec2: ndarray) -> ndarray:
    vec1_x, vec1_y, vec1_z = vec1[0], vec1[1], vec1[2]
    vec2_x, vec2_y, vec2_z = vec2[0], vec2[1], vec2[2]

    res_x = vec1_y * vec2_z - vec1_z * vec2_y
    res_y = vec1_z * vec2_x - vec1_x * vec2_z
    res_z = vec1_x * vec2_y - vec1_y * vec2_x

    if type(res_x) != array and type(res_y) != array and type(res_z) != array:
        res_x = res_x[None]
        res_y = res_y[None]
        res_z = res_z[None]

    res = np.concatenate([res_x, res_y, res_z], axis=0)
    return res

def dot_product(vec1: ndarray, vec2: ndarray) -> ndarray:

    """

    Args:
        - vec1: vector that is going to be projected
        - vec2: vector that is going to get projected on

    """

    if len(vec1) != len(vec2):
        raise ValueError(f"Vectors {vec1} and {vec2} do not have the same length")

    if len(vec1.shape) < 2:
        vec1 = np.expand_dims(vec1, axis=-1)

    if len(vec2.shape) < 2:
        vec2 = np.expand_dims(vec2, axis=-1)

    p = vec1.T @ vec2
    return p