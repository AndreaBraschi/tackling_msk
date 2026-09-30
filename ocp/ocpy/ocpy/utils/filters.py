from scipy import signal
from numpy import ndarray

def butterworth(data: ndarray, dt: float, order: int, cutoff: int, filter_type: str):

    adj_cutoff = 2 * cutoff * dt
    b, a = signal.butter(order, adj_cutoff, btype=filter_type)
    res = signal.filtfilt(b, a, data, padlen=10)

    return res.T