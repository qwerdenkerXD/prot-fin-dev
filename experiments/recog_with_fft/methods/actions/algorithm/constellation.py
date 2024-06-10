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

try:
    AMPL_QUANTILES = {
        i: pd.read_csv(f"../results/_v0.4-exp-uniref_sampling/sample_WINSIZE_{i}.csv", sep=",", usecols=["5%%-Quantile", "95%%-Quantile"]).rename(columns={"5%%-Quantile": 5, "95%%-Quantile": 95})
        for i in range(10, 101, 10)
    }
except Exception as e:
    AMPL_QUANTILES = {
        i: pd.DataFrame({
            5: [float("+inf")] * (WINDOW_SIZE // 2 + 1),
            95: [float("-inf")] * (WINDOW_SIZE // 2 + 1)
        })
        for i in range(10, 101, 10)
    }
    eprint("not filtering amplitudes")


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
