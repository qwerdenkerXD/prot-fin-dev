from scipy import signal
import numpy as np
from tools import *
from os import environ as env
import pandas as pd

# parameters for the STFT
WINDOW_SIZE = int(env.get("WINDOW_SIZE", 30))
WINDOW_TYPE = env.get("WINDOW_TYPE", ["boxcar", "triang", "blackman", "hamming", "hann", "bartlett", "flattop", "parzen", "bohman", "blackmanharris", "nuttall", "barthann", "cosine", "exponential", "tukey", "taylor", "lanczos"]\
                                     [0])
OVERLAP = int(env.get("OVERLAP", 15))
N_PEAKS = int(env.get("N_PEAKS", 0))  # 0 means all

AMPL_QUANTILES = {
    10: {
        95: [1.576, 1.0551, 0.36524, 0.58613, 0.35894, 0.59],
        5: [0.844, 0.54312, 0.0482, 0.1085, 0.04832, 0.025],
    },
    20: {
        95: [1.47, 0.96363, 0.25633, 0.45731, 0.25985, 0.37825, 0.2688, 0.33932, 0.25623, 0.32337, 0.2875],
        5: [0.9265, 0.58593, 0.03295, 0.12298, 0.0336, 0.0614, 0.03466, 0.05002, 0.0334, 0.04578, 0.009],
    },
    30: {
        95: [1.426, 0.92711, 0.20942, 0.41527, 0.20886, 0.31861, 0.21161, 0.28917, 0.22348, 0.27425, 0.21433, 0.25526, 0.20979, 0.24952, 0.2072, 0.278],
        5: [0.96233, 0.60741, 0.02684, 0.13735, 0.02703, 0.0618, 0.0273, 0.04408, 0.02884, 0.03842, 0.02779, 0.03539, 0.02723, 0.03403, 0.02677, 0.00933],
    },
    40: {
        95: [1.40075, 0.90677, 0.18135, 0.39088, 0.18179, 0.29322, 0.18058, 0.25347, 0.18294, 0.23882, 0.19231, 0.23472, 0.1918, 0.22034, 0.18369, 0.21159, 0.18169, 0.20922, 0.18049, 0.20535, 0.20225],
        5: [0.985, 0.62111, 0.02311, 0.14885, 0.0235, 0.06486, 0.02337, 0.04227, 0.0236, 0.03498, 0.02493, 0.03231, 0.02477, 0.02989, 0.02375, 0.02861, 0.02354, 0.02804, 0.02329, 0.02736, 0.00625],
    },
    50: {
        95: [1.3842, 0.89365, 0.16329, 0.37428, 0.16251, 0.27711, 0.16205, 0.23588, 0.16137, 0.21476, 0.16336, 0.20636, 0.17051, 0.20539, 0.17526, 0.19843, 0.16774, 0.18882, 0.16347, 0.18405, 0.16241, 0.18309, 0.16234, 0.18014, 0.15974, 0.2026],
        5: [1.0006, 0.63076, 0.02058, 0.15777, 0.02102, 0.06872, 0.02094, 0.04221, 0.02089, 0.03309, 0.02107, 0.0294, 0.02206, 0.02803, 0.02252, 0.02664, 0.02165, 0.02532, 0.02116, 0.02457, 0.02104, 0.02429, 0.021, 0.02382, 0.02053, 0.0066],
    },
    60: {
        95: [1.37183, 0.88416, 0.15074, 0.36246, 0.14746, 0.26519, 0.14853, 0.22382, 0.14754, 0.20151, 0.1473, 0.18865, 0.14897, 0.18342, 0.15454, 0.18351, 0.1599, 0.18209, 0.15699, 0.17321, 0.1515, 0.16778, 0.14879, 0.16478, 0.14824, 0.16449, 0.14877, 0.16251, 0.14613, 0.16061, 0.16483],
        5: [1.0125, 0.63852, 0.01874, 0.16428, 0.019, 0.07258, 0.01915, 0.04301, 0.01908, 0.03235, 0.01904, 0.02775, 0.01922, 0.02562, 0.01995, 0.02488, 0.02062, 0.0242, 0.02023, 0.02301, 0.01955, 0.02226, 0.01922, 0.02182, 0.01917, 0.02168, 0.01917, 0.02139, 0.01879, 0.02104, 0.00517],
    },
    70: {
        95: [1.362, 0.87674, 0.14103, 0.35387, 0.13535, 0.2555, 0.13765, 0.2146, 0.13711, 0.19186, 0.1364, 0.17814, 0.13628, 0.16964, 0.13774, 0.1663, 0.14221, 0.16682, 0.14723, 0.16737, 0.14869, 0.162, 0.14257, 0.156, 0.13926, 0.15234, 0.13746, 0.15041, 0.13721, 0.15035, 0.13804, 0.14923, 0.13593, 0.147, 0.13481, 0.16614],
        5: [1.02114, 0.64475, 0.01728, 0.16884, 0.01742, 0.07655, 0.01777, 0.04418, 0.01767, 0.03215, 0.01765, 0.0269, 0.01762, 0.02423, 0.01777, 0.02291, 0.01833, 0.02247, 0.01899, 0.02215, 0.01899, 0.02138, 0.01837, 0.02056, 0.01798, 0.02008, 0.01778, 0.01982, 0.01774, 0.01974, 0.01778, 0.01955, 0.01752, 0.01923, 0.0173, 0.00529],
    },
    80: {
        95: [1.35438, 0.87097, 0.13288, 0.34738, 0.12608, 0.2473, 0.12854, 0.20713, 0.12867, 0.18442, 0.1279, 0.17026, 0.12757, 0.16102, 0.12746, 0.15508, 0.12876, 0.15292, 0.13241, 0.15354, 0.13695, 0.15457, 0.13977, 0.15327, 0.13611, 0.14683, 0.13196, 0.143, 0.12966, 0.14046, 0.12845, 0.13916, 0.12836, 0.13922, 0.12922, 0.13867, 0.12762, 0.13672, 0.12605, 0.1359, 0.14263],
        5: [1.02737, 0.64944, 0.01609, 0.17235, 0.01619, 0.08037, 0.01657, 0.04552, 0.01656, 0.03229, 0.0165, 0.02641, 0.01647, 0.0234, 0.01645, 0.0217, 0.01661, 0.02083, 0.01703, 0.02058, 0.01764, 0.02045, 0.01781, 0.02002, 0.0175, 0.01929, 0.01698, 0.01879, 0.01672, 0.01844, 0.01661, 0.01827, 0.01658, 0.01822, 0.01662, 0.01812, 0.01643, 0.01784, 0.01619, 0.01769, 0.0045],
    },
    90: {
        95: [1.348, 0.86621, 0.12572, 0.34209, 0.11924, 0.24054, 0.12052, 0.20078, 0.12149, 0.17829, 0.12099, 0.16394, 0.12034, 0.1545, 0.12019, 0.14786, 0.1202, 0.14354, 0.12132, 0.14214, 0.12442, 0.14272, 0.12836, 0.14401, 0.13115, 0.14451, 0.13114, 0.13993, 0.1262, 0.13556, 0.12363, 0.13282, 0.12184, 0.13097, 0.12103, 0.13004, 0.12102, 0.1302, 0.12184, 0.13007, 0.12067, 0.1283, 0.11918, 0.12703, 0.1188, 0.14389],
        5: [1.03189, 0.65308, 0.01509, 0.17548, 0.01525, 0.08371, 0.01551, 0.04708, 0.01567, 0.03259, 0.01559, 0.02616, 0.01555, 0.02285, 0.01551, 0.02096, 0.01551, 0.01981, 0.01565, 0.01921, 0.016, 0.01904, 0.01654, 0.01902, 0.01679, 0.01882, 0.01666, 0.0183, 0.01621, 0.01774, 0.01592, 0.01739, 0.01569, 0.01717, 0.01564, 0.01703, 0.01563, 0.017, 0.01567, 0.01694, 0.01556, 0.01672, 0.01535, 0.01651, 0.01524, 0.00456],
    },
    100: {
        95: [1.3426, 0.86222, 0.11938, 0.33751, 0.11407, 0.23522, 0.11333, 0.19521, 0.11517, 0.17309, 0.11514, 0.15879, 0.11452, 0.14891, 0.11415, 0.14232, 0.11395, 0.13737, 0.11408, 0.13416, 0.11502, 0.13323, 0.11767, 0.13371, 0.12109, 0.13524, 0.12365, 0.13602, 0.12545, 0.13445, 0.1218, 0.12945, 0.11853, 0.12677, 0.11669, 0.12454, 0.11529, 0.12322, 0.11475, 0.1225, 0.11478, 0.1227, 0.11552, 0.12281, 0.11475, 0.12123, 0.11336, 0.11996, 0.11244, 0.11974, 0.1275],
        5: [1.0351, 0.65587, 0.01423, 0.17838, 0.01452, 0.08649, 0.01457, 0.04882, 0.01487, 0.03302, 0.0148, 0.02606, 0.01475, 0.02248, 0.01474, 0.02042, 0.01471, 0.01914, 0.01472, 0.01832, 0.01485, 0.01789, 0.01514, 0.01777, 0.01559, 0.01781, 0.01591, 0.01772, 0.01591, 0.01741, 0.01562, 0.01691, 0.01524, 0.01653, 0.01501, 0.01627, 0.01489, 0.01612, 0.01482, 0.01601, 0.01482, 0.01598, 0.01487, 0.01597, 0.01478, 0.01578, 0.01458, 0.0156, 0.01442, 0.01554, 0.004],
    }
}


def create_constellation(
        aa_vec: np.ndarray,
        window_size=WINDOW_SIZE,
        n_peaks=N_PEAKS,
        window=WINDOW_TYPE,
        overlap=OVERLAP
        ) -> ConstellationMap:
    """
    The function carries out a windowed fast fourier transformation,
    often called short time FFT (STFT), on the given vector and creates a list
    of frequencies, found by the STFT, and the index of the window they are found

    ...

    Parameter
    ---------
    aa_vec : np.ndarray
        An array of floats representing an amino acid sequence
    window_size : int
        size of a window, related to amino acids -> size 2 means two amino acids
    n_peaks : int, optional
        the number of frequency peaks that are selected from the STFT result,
        if 0 then all are selected (defaults to 0)
    window : str
        The window type for the STFT, read https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window
        for supported ones (defaults to "boxcar")
    overlap : int
        the overlap between two windows during doing the STFT (defaults to half the window size)

    Returns
    -------
     A ConstellationMap, a list of coordinates, meaning the pairs of
     window index and prominent frequency peak of the window
    """

    # adjust window size and overlap if invalid
    if overlap >= window_size:
        overlap = window_size - 1

    if len(aa_vec) < window_size:
        return []

    # executing the STFT
    stft_result = signal.stft(
        aa_vec,
        nperseg=window_size,
        noverlap=overlap,
        window=window
    )

    return stft_to_constellation(*stft_result, n_peaks)


def stft_to_constellation(
        frequencies: np.ndarray,
        window_indexes: np.ndarray,
        stft: np.ndarray,
        n_peaks: int
        ) -> ConstellationMap:
    constellation_map: ConstellationMap = []

    # find and collect the most prominent frequencies from STFT per window
    for amplitudes in stft.T:

        # get rid of complex values to make them comparable
        spectrum: np.ndarray = abs(amplitudes)
        positive = np.argwhere(spectrum > AMPL_QUANTILES[WINDOW_SIZE][95]).flatten()
        negative = np.argwhere(spectrum < AMPL_QUANTILES[WINDOW_SIZE][5]).flatten()

        for i, filt in enumerate((negative, positive)):
            # find peaks
            peaks: List[int] = find_peaks(spectrum[filt], n_peaks)

            constellation_map.append(tuple((int(filt[peak]), float(spectrum[filt][peak]), i) for peak in peaks))

    return constellation_map


def find_peaks(spectrum: np.ndarray, n_peaks: int) -> List[int]:
    # prominence=0 includes all peaks, but weights their prominence as well
    peaks, props = signal.find_peaks(spectrum, prominence=0)

    # Only want the most prominent peaks
    peaks: List[Tuple[int, int]] = sorted(zip(props["prominences"], peaks), reverse=True)
    if n_peaks:
        peaks = peaks[:n_peaks]

    return [int(p[1]) for p in peaks]
