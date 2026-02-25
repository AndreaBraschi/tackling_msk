import os
from numpy import ndarray, array, concatenate
from pandas import Series, DataFrame
from opensim import Model, Coordinate, Body, Point, Vec3
from contact_model.sphere_to_cylinder import compute_x_and_x_dot
from osim_utils.read import readStoFile

# ----------------------------------------------------------------------------------------------------------------------
# Declare global paths
# ----------------------------------------------------------------------------------------------------------------------
trial = "P08R0002"
ik_dir = f"X:/Health/ResearchProjects/DCazzola/RE-FH1285-c-spine-injuries-in-rugby/Darios data/Phase 3/OpenSim/P8/ik/{trial}"

# ----------------------------------------------------------------------------------------------------------------------
# Read .sto file
# ----------------------------------------------------------------------------------------------------------------------
posData = readStoFile(os.path.join(ik_dir, "optimized_scale_and_markers_BodyKinematics_pos_global.sto"))
velData = readStoFile(os.path.join(ik_dir, "optimized_scale_and_markers_BodyKinematics_vel_global.sto"))

ik_file = f"{trial}_IK.mot"
ikData = readStoFile(os.path.join(ik_dir, ik_file))

stateFileName: str = "model_StatesReporter_states.sto"
stateData: DataFrame = readStoFile(os.path.join(ik_dir, stateFileName))

time = stateData["time"]
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
cylinderVelx: Series = velData["punching_bag_2_X"].values[:, None]
cylinderVely: Series = velData["punching_bag_2_Y"].values[:, None]
cylinderVelz: Series = velData["punching_bag_2_Z"].values[:, None]

cylinderVel: ndarray = concatenate([cylinderVelx, cylinderVely, cylinderVelz], axis=-1)

# ----------------------------------------------------------------------------------------------------------------------
# Double loop: 1) loop through time, 2) at each time point, loop through the coordinate set, read the coordinate value
# from the IK output and assign that value to the model coordinate
# ----------------------------------------------------------------------------------------------------------------------
for n in range(50, 120):
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

    d, vel = compute_x_and_x_dot(model, sphere_loc, cylinderVel[n][:, None], clavicleBody, s)
    print(f"Penetration: {d}")
    print(f"Velocity: {vel}")
    print("\n")
