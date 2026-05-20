from opensim import Model
from osim_utils.tools.run_ik import run_ik
from osim_utils.read import readTrc

from numpy import ndarray
from pandas import DataFrame

def run_single_ik(model_filepath: str, mot_filepath: str, setup_filepath: str, trc_filepath: str, output_dir: str) -> int:

    # load model
    model = Model(model_filepath)

    # read from the trc file the initial and final times
    trc_df: DataFrame = readTrc(trc_filepath)
    time: ndarray = trc_df['Time'].values

    init: float = float(time[0])
    end: float = float(time[-1])

    run_ik(model, trc_filepath, mot_filepath, output_dir, init, end, setup_filepath)

    return 0


