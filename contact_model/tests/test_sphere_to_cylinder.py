from numpy import ndarray, array, zeros, ones, linspace, float64
from contact_model.contact_forces.smooth_forces import smooth_hunt_crossley
from osim_utils.read import readStoFile
from utils.geometry import compute_effective_radius
from utils.polynomials import fit_cubic, cubic
import matplotlib.pyplot as plt



# cylinder radius
cylinder_r: float = 0.1702085

# sphere radius
sphere_r: float = 0.025

# define contact model parameters
k: float = 2300832.0
c: float = 2.5

# pre-allocate force array
T: int = 121
time: ndarray = linspace(0, 1, T)
force_arr: ndarray = zeros((T, ), dtype=float64)
x_points: ndarray = array([0, 0.3, 0.5, 1])
y_points: ndarray = array([-0.01, 0.003, 0.007, 0.000025])

coefficients: ndarray = fit_cubic(x_points, y_points)
x: ndarray = cubic(time, coefficients)
x_dot: ndarray = ones((T, ), dtype=float64) * 5

# ----------------------------------------------------------------------------------------------------------------------
# Double loop: 1) loop through time, 2) at each time point, loop through the coordinate set, read the coordinate value
# from the IK output and assign that value to the model coordinate
# ----------------------------------------------------------------------------------------------------------------------
for n in range(T):
    print(f"Frame number: {n}")
    print(f"Penetration: {float(x[n])}")
    print(f"Velocity: {float(x_dot[n])}")

    R: float = compute_effective_radius(sphere_r, cylinder_r)
    force = smooth_hunt_crossley(float(x[n]), float(x_dot[n]), R, k, c)
    force_arr[n] = force

    print(f"Scalar force: {force}")
    print("\n")


plt.figure()
plt.title("Scalar Force over time")
plt.plot(time, force_arr)
plt.xlabel("Time (s)")
plt.ylabel("Force (N)")
plt.show()

plt.figure()
plt.title("Force vs. Indentation")
plt.plot(x, force_arr)
plt.xlabel("Penetration (m)")
plt.ylabel("Force (N)")
plt.show()