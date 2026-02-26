import os
from numpy import ndarray, array, concatenate, zeros, linspace, float64
from pandas import Series, DataFrame
from opensim import Model, Coordinate, Body, Point, Vec3
from contact_model.sphere_to_cylinder import compute_x_and_x_dot
from contact_model.contact_forces.smooth_forces import smooth_hunt_crossley
from osim_utils.read import readStoFile
from utils.numerical_methods import central_diff
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------
# Declare global paths
# ----------------------------------------------------------------------------------------------------------------------
subject: str = "P8"
trial = "P08R0002"
ik_dir = f"X:/Health/ResearchProjects/DCazzola/RE-FH1285-c-spine-injuries-in-rugby/Darios data/Phase 3/OpenSim/P8/ik/{trial}"
bodyKinematics_dir: str = f"{ik_dir}/bodyKinematics"
# ----------------------------------------------------------------------------------------------------------------------
# Read .sto file
# ----------------------------------------------------------------------------------------------------------------------
posData = readStoFile(os.path.join(bodyKinematics_dir, f"{subject}_BodyKinematics_pos_global.sto"))
velData = readStoFile(os.path.join(bodyKinematics_dir, f"{subject}_BodyKinematics_vel_global.sto"))

ik_file = f"{trial}_IK.mot"
ikData = readStoFile(os.path.join(ik_dir, ik_file))

stateFileName: str = f"{subject}_StatesReporter_states.sto"
stateData: DataFrame = readStoFile(os.path.join(bodyKinematics_dir, stateFileName))

time = stateData["time"]
T: int = len(time)   # period
coordinates = stateData.iloc[:, 1::2]
speeds = stateData.iloc[:, 2::2]
t_flags = ["tx", "ty", "tz", "t1", "t2", "t3"]

# ----------------------------------------------------------------------------------------------------------------------
# Load .osim Model and get CoordinateSet and BodySet
# ----------------------------------------------------------------------------------------------------------------------
model = Model(r"C:\Users\ab3758\Downloads\modifiedWrapping_pb.osim")

s = model.initSystem()
coordinateSet = model.getCoordinateSet()
bodySet = model.getBodySet()
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
k: float = 1533888.0
c: float = 2.0

# pre-allocate force array
force_arr: ndarray = zeros((T, ), dtype=float64)

# ----------------------------------------------------------------------------------------------------------------------
# Double loop: 1) loop through time, 2) at each time point, loop through the coordinate set, read the coordinate value
# from the IK output and assign that value to the model coordinate
# ----------------------------------------------------------------------------------------------------------------------
for n in range(50, 60):
    print(f"Frame number: {n}, Time: {time[n]}")
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

    print(f"Penetration: {d}")
    print(f"Velocity: {vel}")

    force = smooth_hunt_crossley(d, vel, sphere_r, k, c)
    force_arr[n] = force

    print(f"Scalar force: {force}")
    print("\n")


plt.figure()
plt.title("Scalar Force")
plt.plot(time[50:60], force_arr[50:60])
plt.xlabel("Time (s)")
plt.ylabel("Force (N)")
plt.show()