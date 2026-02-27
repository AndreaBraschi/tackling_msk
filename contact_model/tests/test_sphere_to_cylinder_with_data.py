import os
from numpy import ndarray, array, concatenate, zeros, ones, linspace, float64
from pandas import Series, DataFrame
from opensim import Model, Coordinate, Body
from contact_model.sphere_to_cylinder import compute_x_and_x_dot
from contact_model.contact_forces.smooth_forces import smooth_hunt_crossley
from osim_utils.read import readStoFile
from utils.geometry import compute_effective_radius
from utils.polynomials import fit_cubic, cubic
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt



# ----------------------------------------------------------------------------------------------------------------------
# Load .osim Model and get CoordinateSet and BodySet
# ----------------------------------------------------------------------------------------------------------------------
model = Model(r"osim_models\modifiedWrapping_pb.osim")

s = model.initSystem()
coordinateSet = model.getCoordinateSet()
bodySet = model.getBodySet()

# ----------------------------------------------------------------------------------------------------------------------
# Read .sto file
# ----------------------------------------------------------------------------------------------------------------------
velData = readStoFile(r"files/P8_BodyKinematics_vel_global.sto")

stateFileName: str = f"files/P8_StatesReporter_states.sto"
stateData: DataFrame = readStoFile(stateFileName)

time = stateData["time"]
T: int = len(time)   # period
coordinates = stateData.iloc[:, 1::2]
speeds = stateData.iloc[:, 2::2]
t_flags = ["tx", "ty", "tz", "t1", "t2", "t3"]

# sphere location in local clavicle reference frame
sphere_loc: ndarray = array([[-0.05], [0.015], [0.1]])

# get cylinder velocity
cylinderVelx: Series = velData["punching_bag_X"].values[:, None]
cylinderVely: Series = velData["punching_bag_Y"].values[:, None]
cylinderVelz: Series = velData["punching_bag_Z"].values[:, None]

cylinderVel: ndarray = concatenate([cylinderVelx, cylinderVely, cylinderVelz], axis=-1)

# cylinder radius
cylinder_r: float = 0.1702085

# sphere radius
sphere_r: float = 0.025

# define contact model parameters
k: float = 2300832.0
c: float = 2.5

# pre-allocate force array and indentation array
force_arr: ndarray = zeros((T, ), dtype=float64)
x_arr: ndarray = zeros((T, ), dtype=float64)
# ----------------------------------------------------------------------------------------------------------------------
# Double loop: 1) loop through time, 2) at each time point, loop through the coordinate set, read the coordinate value
# from the IK output and assign that value to the model coordinate
# ----------------------------------------------------------------------------------------------------------------------
for n in range(T):
    print(f"Frame number: {n}")
    for i, coord in enumerate(coordinateSet):
        coordName: str = coordinates.columns[i].split("/")[-2]
        coordValue: float = float(coordinates.iloc[n, i])
        speedValue: float = float(speeds.iloc[n, i])
        coordinate: Coordinate = model.getCoordinateSet().get(i)
        coordinate.setValue(s, coordValue)
        coordinate.setSpeedValue(s, speedValue)

    model.realizeVelocity(s)
    clavicleBody: Body = bodySet.get("rclavicle")

    d, vel = compute_x_and_x_dot(model, sphere_loc, cylinderVel[n][:, None], clavicleBody, sphere_r, cylinder_r, s)

    x_arr[n] = n

    print(f"Penetration: {d}")
    print(f"Velocity: {vel}")

    R: float = compute_effective_radius(sphere_r, cylinder_r)
    force = smooth_hunt_crossley(d, vel, R, k, c)
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
plt.plot(x_arr, force_arr)
plt.xlabel("Penetration (m)")
plt.ylabel("Force (N)")
plt.show()