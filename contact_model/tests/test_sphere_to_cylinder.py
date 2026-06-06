from osim_utils.read import readStoFile
from ..utils.geometry import compute_effective_radius
from ..utils.polynomials import fit_cubic, cubic
from ..contact_forces.smooth_forces import *
from numpy import sqrt as np_sqrt
from numpy import ndarray, array, tanh, zeros, ones, linspace, float64
import matplotlib.pyplot as plt


# cylinder parameters: radius and height
cylinder_r: float = 0.1702085
cylinder_h: float = 0.53859999999999997

# here we create 2 numpy arrays representing the center of the top and bottom surfaces of the cylinder, respectively, in
# the cylinder local frame, whose origin is located at the cylinder geometrical center.
cylinder_top_arr: ndarray = array([[0.0], [cylinder_h / 2], [0.0]], dtype=float64)
cylinder_bottom_arr: ndarray = cylinder_top_arr * (-1)

# sphere radius
sphere_r: float = 0.025

# define contact model parameters
k: float = 1000000.0
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
# x_dot: ndarray = quadratic(time, coefficients[:-1])

# plot x and x_dot
fig, axes = plt.subplots(2, 1)
axes[0].plot(time, x)
axes[0].set_title('indentation')

axes[1].plot(time, x_dot)
axes[1].set_title('rate of change')
fig.tight_layout()
plt.show()

b: float = 50
bc: ndarray = b * x
tanh_bc: ndarray = tanh(bc)

bc_dot: ndarray = b * (x_dot + 2 / (3 * c))
tanh_bc_dot: ndarray = tanh(bc_dot)

f_p = compute_fp(np_sqrt(x ** 2 + 1e-8))
f_v = compute_fv(x_dot, c)

f_p_smooth = f_p * (0.5 + 0.5 * tanh(b * x))
f_v_smooth = f_v * (0.5 + 0.5 * tanh(b * (x_dot + 2 / (3 * c))))

fig, axes = plt.subplots(2, 2)
axes[0, 0].plot(x, bc)
axes[0, 0].set_title('bc')

axes[0, 1].plot(x, tanh_bc)
axes[0, 1].set_title('tanh_bc')

axes[1, 0].plot(x, f_p)
axes[1, 0].set_title('fp')

axes[1, 1].plot(x, f_p_smooth)
axes[1, 1].set_title('f_p_smooth')

fig.tight_layout()
plt.show()

fig, axes = plt.subplots(2, 2)

axes[0, 0].plot(time, bc_dot)
axes[0, 0].set_title('50_bc_dot')

axes[0, 1].plot(time, tanh_bc_dot)
axes[0, 1].set_title('tanh_bc_dot')

axes[1, 0].plot(time, f_v)
axes[1, 0].set_title('f_v')

axes[1, 1].plot(f_v_smooth)
axes[1, 1].set_title('f_v_smooth')

fig.tight_layout()
plt.show()
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