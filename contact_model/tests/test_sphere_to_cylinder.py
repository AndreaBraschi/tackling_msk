import os
import math
from numpy import ndarray, array, concatenate
from pandas import Series, DataFrame
from opensim import Model
from contact_model.sphere_to_cylinder import compute_penetration
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

time = ikData["time"]
ikCoords = ikData.iloc[:, 1:]
t_flags = ["tx", "ty", "tz", "t1", "t2", "t3"]

# ----------------------------------------------------------------------------------------------------------------------
# Load .osim Model and get CoordinateSet and BodySet
# ----------------------------------------------------------------------------------------------------------------------
model = Model("X:\Health\ResearchProjects\DCazzola\RE-FH1285-c-spine-injuries-in-rugby\Darios data\Phase 3\OpenSim\P8/scaled_model_markers.osim")
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
        coordName = ikCoords.columns[i]
        coordValue = ikCoords.iloc[n, i] * (math.pi / 180)
        if len(coordName.split("_")) > 1:
            if coordName.split("_")[1] in t_flags:
                coordValue /= (math.pi / 180)

        model.getCoordinateSet().get(i).setValue(s, coordValue)

    model.realizePosition(s)

    cylinderBody = bodySet.get("punching_bag_2")
    clavicleBody = bodySet.get("rclavicle")

    d = compute_penetration(sphere_loc, cylinderBody, cylinderVel[n][:, None], clavicleBody, s)
    print("\n")