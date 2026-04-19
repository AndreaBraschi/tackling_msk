from numpy import (
    array, float64, abs, eye, all
)
from numpy.typing import NDArray
from numpy.linalg import inv
from typing import Tuple

def compute_coefficients(x_mat: NDArray[float64], y: NDArray[float64]) -> NDArray[float64]:
    # here we compute the coefficients of whatever polynomial by solving a linear system of equations. The idea is
    # to represent y = X @ coefficients.
    # --> coefficients = X{^-1} @ y

    x_inv: NDArray[float64] = inv(x_mat)

    # check that it's close to an identity matrix
    tol: float = 1e-6
    shape: int = x_mat.shape[0]
    identity: NDArray[float64] = eye(shape)
    diff = abs(x_inv @ x_mat - identity)

    if all(diff > tol):
        raise ValueError("x_inv and x_mat product is not close enough to identity matrix.")

    coefficients: NDArray[float64] = x_inv @ y
    return coefficients

def fit_quadratic(x: NDArray[float64], y: NDArray[float64]) -> NDArray[float64]:

    expected_shape: Tuple[int, ...] = (3, )

    # check shape of arrays
    if x.shape != expected_shape:
        raise ValueError("x must have 3 points for quadratic fitting")

    elif y.shape != expected_shape:
        raise ValueError("y must have 3 points for quadratic fitting")

    # create a matrix that contains x{^2}, x and 1 for each element of the original x array
    x_mat: NDArray[float64] = array([
        [x[0] ** 2, x[0], 1],
        [x[1] ** 2, x[1], 1],
        [x[2] ** 2, x[2], 1],
    ])

    coefficients: NDArray[float64] = compute_coefficients(x_mat, y)
    return coefficients

def fit_cubic(x: NDArray[float64], y: NDArray[float64]) -> NDArray[float64]:
    min_number_of_points: int = 4
    expected_shape: Tuple[int, ...] = (min_number_of_points,)

    # check shape of arrays
    if x.shape != expected_shape:
        raise ValueError(f"x must have {min_number_of_points} points for quadratic fitting. Current array has: {x.shape}")

    elif y.shape != expected_shape:
        raise ValueError(f"y must have {min_number_of_points} points for quadratic fitting. Current array has: {y.shape}")

    # create a matrix that contains x{^2}, x and 1 for each element of the original x array
    x_mat: NDArray[float64] = array([
        [x[0] ** 3, x[0] ** 2, x[0], 1],
        [x[1] ** 3, x[1] ** 2, x[1], 1],
        [x[2] ** 3, x[2] ** 2, x[2], 1],
        [x[3] ** 3, x[3] ** 2, x[3], 1],
    ])

    coefficients: NDArray[float64] = compute_coefficients(x_mat, y)
    return coefficients

def quadratic(x: NDArray[float64], coefficients: NDArray[float64]) -> NDArray[float64]:
    a, b, c = coefficients
    y: NDArray[float64] = (a * x ** 2) + (b * x) + c
    return y

def cubic(x: NDArray[float64], coefficients: NDArray[float64]) -> NDArray[float64]:
    a, b, c, d = coefficients
    y: NDArray[float64] = (a * x ** 3) + (b * x ** 2) + (c * x) + d
    return y