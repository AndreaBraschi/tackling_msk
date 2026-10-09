from scipy import signal
import numpy as np

def butterworth(time: np.ndarray, data: np.ndarray, dt: float, order: int, cutoff: int, filter_type: str):

    fs = 1 / np.mean(np.diff(time))
    adj_cutoff = cutoff / (0.5 * fs)
    b, a = signal.butter(order, adj_cutoff, btype=filter_type)
    res = signal.filtfilt(b, a, data, padlen=10)

    return res.T