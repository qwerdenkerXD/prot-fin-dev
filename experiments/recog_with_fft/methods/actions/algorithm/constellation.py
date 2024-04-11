from scipy import signal
import numpy as np
from tools import *


def create_constellation(
        aa_vec: np.ndarray,
        window_size: int,
        n_peaks=0,
        window="boxcar",
        **kwargs
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
    **kwargs:
        overlap : int
            the overlap between two windows during doing the STFT (defaults to half the window size)

    Returns
    -------
     A ConstellationMap, a list of coordinates, meaning the pairs of
     window index and prominent frequency peak of the window
    """

    # extract overlap from keyword arguments
    overlap: int = kwargs.get("overlap", window_size // 2)

    # adjust window size and overlap if invalid
    if len(aa_vec) < window_size:
        window_size = len(aa_vec)
    if overlap >= window_size:
        overlap = window_size - 1

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
    for window_idx, amplitudes in zip(window_indexes, stft.T):

        # get rid of complex values to make them comparable
        spectrum: np.ndarray = abs(amplitudes)

        # find peaks
        peaks: List[int] = find_peaks(spectrum, n_peaks)

        for peak in peaks:
            frequency = frequencies[peak]
            constellation_map.append([window_idx, frequency, spectrum[peak]])

    return constellation_map


def find_peaks(spectrum: np.ndarray, n_peaks: int) -> List[int]:
    # prominence=0 includes all peaks, but weights their prominence as well
    peaks, props = signal.find_peaks(spectrum, prominence=0)

    # Only want the most prominent peaks
    peaks: List[Tuple[int, int]] = sorted(zip(props["prominences"], peaks), reverse=True)
    if n_peaks:
        peaks = peaks[:n_peaks]

    return [p[1] for p in peaks]
