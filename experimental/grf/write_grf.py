import os
import json
from numpy import array, arange, einsum, concatenate
from osim_utils.rotations import QUALYSIS_TO_OPENSIM
from osim_utils.write import writeMotFromDataFrame
import matplotlib.pyplot as plt
from pandas import DataFrame
from numpy import float64
from numpy.typing import NDArray as ndarray


def write_grf(json_filepath: str, output_dir: str, hz: float) -> int:
    """

    This function reads the Qualisys force data that was exported as a .json file, converts them to OpenSim convention
    it finally writes the transformed forces and moments as a .mot file to be used in other OpenSim analysis.

    :param json_filepath: path to json file containing the experimental force data
    :param mot_filepath: path mot file that will be written
    :return:
    """

    force_indices = [0, 1, 2]
    moment_indices = [3, 4, 5]

    with open(json_filepath, "r") as f:
        data = json.load(f)

    fp_one: ndarray[float64] = array(data["ForcePlates"][0]['Parts'][0]['Values'])
    fp_two: ndarray[float64] = array(data["ForcePlates"][1]['Parts'][0]['Values'])

    # qualisys raw data
    force_one_q = fp_one[:, force_indices]
    force_two_q  = fp_two[:, force_indices]

    moment_one_q  = fp_one[:, moment_indices]
    moment_two_q  = fp_two[:, moment_indices]


    # rotate force and moment vectors, following OpenSim convention: x: ap, y: vertical, z: medio-lat (pointing left)
    # for this dataset, we have x: medio-lat (pointing left), y: ap, z: vertical
    R: ndarray = QUALYSIS_TO_OPENSIM
    force_one = einsum('ij, tj -> ti', R, force_one_q)
    force_two  = einsum('ij, tj -> ti', R, force_two_q)

    moment_one  = einsum('ij, tj -> ti', R, moment_one_q)
    moment_two  = einsum('ij, tj -> ti', R, moment_two_q)


    # read information about the data: sampling frequency, total period, etc.
    init_frame = data["Timebase"]["Range"]["Start"] - 1
    final_frame = len(force_one_q)

    dt = 1 / hz
    init_time = init_frame / hz
    end_time = final_frame / hz

    time: ndarray[float64] = arange(init_time, end_time, dt, dtype=float64)

    filename = os.path.basename(os.path.normpath(json_filepath)).split(".")[0]

    # -------------------- figures -------------------- #
    fig_dir = os.path.join(output_dir, "figures")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # plot forces
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(time, force_one_q[:, 0], color='coral', label="right")
    axes[0].plot(time, force_two_q[:, 0], color='seagreen', label="left")
    axes[0].set_title('Qualysis x')

    axes[1].plot(time, force_one_q[:, 1], color='coral', label="right")
    axes[1].plot(time, force_two_q[:, 1], color='seagreen', label="left")
    axes[1].set_title('Qualysis y')

    axes[2].plot(time, force_one_q[:, 2], color='coral', label="right")
    axes[2].plot(time, force_two_q[:, 2], color='seagreen', label="left")
    axes[2].set_title('Qualysis z')

    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(fig_dir, filename, "qualisys_force.png"))
    # plt.show()
    plt.close()


    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(time, force_one[:, 0], color='coral', label='right')
    axes[0].plot(time, force_two[:, 0], color='seagreen', label='left')
    axes[0].set_title('OpenSim x')
    plt.legend()

    axes[1].plot(time, force_one[:, 1], color='coral', label='right')
    axes[1].plot(time, force_two[:, 1], color='seagreen', label='left')
    axes[1].set_title('OpenSim y')
    plt.legend()

    axes[2].plot(time, force_one[:, 2], color='coral', label='right')
    axes[2].plot(time, force_two[:, 2], color='seagreen', label='left')
    axes[2].set_title('OpenSim z')
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, filename, "opensim_force.png"))
    #plt.show()
    plt.close()


    # plot moments
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(time, moment_one_q[:, 0], color='coral', label="right")
    axes[0].plot(time, moment_two_q[:, 0], color='seagreen', label="left")
    axes[0].set_title('Qualysis x')
    plt.legend()

    axes[1].plot(time, moment_one_q[:, 1], color='coral', label="right")
    axes[1].plot(time, moment_two_q[:, 1], color='seagreen', label="left")
    axes[1].set_title('Qualysis y')
    plt.legend()

    axes[2].plot(time, moment_one_q[:, 2], color='coral', label="right")
    axes[2].plot(time, moment_two_q[:, 2], color='seagreen', label="left")
    axes[2].set_title('Qualysis z')
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, filename, "qualisys_moment.png"))
    plt.close()


    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(time, moment_one[:, 0], color='coral', label="right")
    axes[0].plot(time, moment_two[:, 0], color='seagreen', label="left")
    axes[0].set_title('Qualysis x')
    plt.legend()

    axes[1].plot(time, moment_one[:, 1], color='coral', label="right")
    axes[1].plot(time, moment_two[:, 1], color='seagreen', label="left")
    axes[1].set_title('Qualysis y')
    plt.legend()

    axes[2].plot(time, moment_one[:, 2], color='coral', label="right")
    axes[2].plot(time, moment_two[:, 2], color='seagreen', label="left")
    axes[2].set_title('Qualysis z')
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, filename, "opensim_moment.png"))
    plt.close()


    # ----------- write to .mot file ---------- #
    # list of labels
    labels = ["time",
              "force_lx", "force_ly", "force_lz", "moment_lx", "moment_ly", "moment_lz",
              "force_rx", "force_ry", "force_rz", "moment_rx", "moment_ry", "moment_rz"]

    arr: ndarray[float64] = concatenate([time[:, None], force_two, moment_two, force_one, moment_one], axis=-1)

    df = DataFrame(data=arr, columns=labels)

    mot_filepath = os.path.join(output_dir, filename + ".mot")

    writeMotFromDataFrame(df, mot_filepath, inDegrees="no")

    return 0