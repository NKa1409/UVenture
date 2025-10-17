import copy
import datetime
import sys
import numpy as np
import scipy
from scipy.sparse import diags

def baseline_als(y, lam=1e6, p=0.01, niter=10, window_min_vals=4):
    """Asymmetric Least Squares (ALS) baseline correction with correct matrix sizing."""
    # Parameters:
    # y: input data (1D array)
    # lam: smoothness parameter (larger values give smoother baseline)
    # p: asymmetry parameter (0 < p < 1, larger values give more weight to the left side)
    # niter: number of iterations for convergence
    # Returns:
    if not isinstance(window_min_vals, int):
        print("window_min_vals must be an integer!")
        try:
            window_min_vals = int(window_min_vals)
        except ValueError:
            print("Setting to integer did not work. Setting to default value of 4!")
            window_min_vals = 4
    # Get rolling window minimum values
    y_mins = []
    for i in range(len(y)):
        if i < window_min_vals:
            y_mins.append( sorted(y[:i + window_min_vals])[ int(window_min_vals*0.8)-1 ] )
        elif i > len(y)-(window_min_vals+1):
            y_mins.append( sorted(y[i-window_min_vals:])[ int(window_min_vals*0.8)-1 ] )
        else:
            y_mins.append( sorted(y[i - window_min_vals:i + window_min_vals])[ int(window_min_vals*1.6)-1 ] )

    y_original = y.copy()     
    y = np.array(y_mins)
    L = len(y)
    # Fix: Adjust differentiation matrix to match length L
    D = scipy.sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L))
    DTD = (D.T @ D).tocsc()  # Precompute once
    w = np.ones(L)
    W_diag = w.copy()
    for _ in range(niter):
        W_diag[:] = w
        # Fix: Use the same size for D as y
        Z = diags(W_diag) + lam * DTD
        # Solve using np.linalg.solve
        baseline = scipy.sparse.linalg.spsolve(Z, w * y)
        # Update weights
        w = p * (y > baseline) + (1 - p) * (y < baseline)
    for i in range(len(baseline)):
        if baseline[i] <= 0:
            baseline[i] = 0
    for i in range(len(baseline)):
        if y_original[i] - baseline[i] < 0:
            baseline[i] = y_original[i]
    return baseline


def get_window_size_by_frequency(intensityvals, timevals, min_width=3, max_width=45):
    if min_width < 0:
        print("min_width must be greater than 0!  Setting to default value of 3!")
        min_width = 3
    if max_width < min_width:
        print("max_width must be greater than min_width!  Setting to default value of 45!")
        max_width = 45
    # Estimate a window size based on the frequency of noise oscillations in multiple windows across the whole series
    window_size = 0
    dom_freqs = []
    for i in range(0, len(intensityvals), int(len(intensityvals)/20)):
        if i == 0:
            continue
        elif i + int(len(intensityvals)/20) >= len(intensityvals):
            break
        else:
            window = intensityvals[i:i + int(len(intensityvals)/20)]
            if sum(window) == 0:
                continue
            window_mean = np.mean(window)
            window_std = np.std(window)
            if window_std > 0.1 * window_mean:
                window_size += 1
            fs = len(timevals) / (max(timevals) - min(timevals))
            N = len(window)
            yf = scipy.fft.fft(window)
            xf = scipy.fft.fftfreq(N, 1 / fs)
            idxs = np.where(xf >= 0)
            freqs = xf[idxs]
            mags = np.abs(yf[idxs])
            if len(mags) < 6:
                continue
            dom_freqs.append(freqs[np.argmax(mags[4:])])
    if len(dom_freqs) > 0:
        #print("dom_freqs: " + str(dom_freqs))
        if not sum(dom_freqs) <= 0.000001:
            avg_dom_freq = sum(dom_freqs) / len(dom_freqs) # 1/s
            window_size = int(1 / avg_dom_freq) # 1 oscillation every x seconds
            window_size = window_size * 2 # Increase the window size by a factor of 2 to get better smoothing results --> a window size of 1 oscillation is not enough!
            window_size = int(window_size * fs) # in measurements --> Window size is so that the window contains 1 oscillation and depending on the sampling rate x measurements
        else:
            window_size = (max_width-min_width) / 2
        if window_size < min_width:
            window_size = min_width
        elif window_size > max_width:
            window_size = max_width
    else:
        window_size = (max_width-min_width) / 2
    return int(window_size)


def do_smoothing_without_effecting_peaks(y, window_size=5):
    if window_size < 3:
        print("window_size must be greater than 3!  Setting to default value of 5!")
        window_size = 5
    s_intensityvals = scipy.signal.savgol_filter(y, window_length=window_size, polyorder=2, mode="nearest")
    peaks, _ = scipy.signal.find_peaks(s_intensityvals,
                                       prominence=np.std(s_intensityvals)*3, width=(5, 30))
    timeseries_lists = []
    for peak in range(len(peaks)+1):
        if len(peaks) == 0:
            break
        if peak == 0:
            left = 0
            right = int(peaks[peak] - window_size)
            if right <= 0:
                #print("right <= 0: " + str(right))
                continue
        elif peak == len(peaks):
            left = int(peaks[peak-1] + window_size)
            right = len(y)
            if left > right:
                #print("left > right: " + str(left) + " > " + str(right))
                continue
        else:
            left = int(peaks[peak-1] + window_size)
            right = int(peaks[peak] - window_size)
            if left > right:
                #print("left > right: " + str(left) + " > " + str(right))
                continue
        smooth_vals = scipy.signal.savgol_filter(y[left:right], window_length=window_size, polyorder=2, mode="nearest")
        indices_list = list(range(left, right, 1))
        for i in range(len(smooth_vals)):
            timeseries_lists.append( (indices_list[i], smooth_vals[i]) )
    new_int_vals = copy.deepcopy(y)
    for index,value in timeseries_lists:
        new_int_vals[index] = value
    return new_int_vals


def get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=2, max_width_seconds=40):
    if min_width_seconds < 0:
        print("min_width_seconds must be greater than 0!")
        print("Setting to default value of 2 seconds!")
        min_width_seconds = 2
    if max_width_seconds < min_width_seconds:
        print("max_width_seconds must be greater than min_width_seconds!")
        print("Setting to default value of 40 seconds!")
        max_width_seconds = 40
    starttime = datetime.datetime.now()
    # Do first smoothing
    window_size = get_window_size_by_frequency(intensities, times, min_width=2, max_width=30)
    intensityvals = do_smoothing_without_effecting_peaks(intensities, window_size=window_size)
    # Apply baseline correction
    baseline = baseline_als(intensityvals, lam=1e4, p=0.05, niter=100, window_min_vals=window_size)
    corrected_intensity = intensityvals - baseline
    smoothed_intensity = do_smoothing_without_effecting_peaks(corrected_intensity, window_size=window_size)
    # Find peaks with adaptive height and width detection
    peaks, properties = scipy.signal.find_peaks(smoothed_intensity,
                                prominence=np.std(smoothed_intensity), width=(min_width_seconds, max_width_seconds))
    return peaks, properties, smoothed_intensity
    



if __name__ == "__main__":
    print("This is a module with functions for peak detection and baseline correction. It is not meant to be run directly.")
    sys.exit(0)