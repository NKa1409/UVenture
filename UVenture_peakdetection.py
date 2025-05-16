import ast
import copy
import datetime
import math
import operator
import pathlib
import sys
import traceback
import PIL
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import matplotlib
import scipy.fft
import UVenture
import matplotlib.gridspec
import numpy as np
import psutil
import scipy
import MS_functions
import pyteomics
import flask
import time
import threading
import werkzeug
import os
from pyteomics import mzml
import similaritymeasures

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import shutil



def create_gaussian_function(xcenter, xsigma, spacing):
    def gaussian(x, mu, sig):
        return ( 1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2) )
    x_vals = []
    y_vals = []




class PeakFinding:
    def __init__(self, ms_file, **kwargs):
        self.ms_file = ms_file
        self.peak_properties_list = []

        default_kwargs = {"mass_deviation": 11,
                          "minimum_required_max_intensity": 100000,
                          "minimum_peak_width_seconds": 0.1,
                          "maximum_peak_width_seconds": 10,
                          "minimum_peak_prominence": 50000,
                          "peak_finding_log_filepath": str(self.ms_file.parentfolder + "/" + "peak_finding_log.txt"),
                          "round_for_masses_summary": 4}
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["peak_finding_log_filepath"], exist_ok=True)

        print("Starting peak finding...")
        ordered_mass_intensity_dict = self.get_all_available_masses_in_list()
        print("Finished getting all available masses in list!")
        print(ordered_mass_intensity_dict)

        while len(ordered_mass_intensity_dict) > 0:
            curr_mass = list(ordered_mass_intensity_dict.keys())[0]
            mass_evaluation = self.search_peaks_for_mass(curr_mass)
            print(mass_evaluation)
            for peak in mass_evaluation:
                self.peak_properties_list.append(peak)
            del ordered_mass_intensity_dict[curr_mass]
        print("Finished peak finding!")
        self.peak_properties_list = sorted([peak for peak in self.peak_properties_list], key=lambda x: x[2])
        print("Peak properties list: " + str(self.peak_properties_list))

    def search_peaks_for_mass(self, mass):
        print("Evaluating peak for mass: " + str(mass))
        xic = MS_functions.get_xic(self.ms_file.rawdata, mass, ((self.kwargs["mass_deviation"] * mass) / 1000000), requested_filter_mode="Full scan")
        times = xic[0]
        intensities = xic[1]
        window, order = 5, 3
        intensities = scipy.signal.savgol_filter(intensities, window, order, mode="nearest")
        peak_properties = scipy.signal.find_peaks(intensities, height=max(intensities) / 100, distance=2, prominence=max(intensities) / 50, width=(2, 20))
        identified_peaks = []
        for element in range(len(peak_properties[1]["peak_heights"])):
            if peak_properties[1]["peak_heights"][element] < self.kwargs["minimum_required_max_intensity"] or \
                    peak_properties[1]["prominences"][element] < self.kwargs["minimum_peak_prominence"]:
                continue
            one_peak = []
            one_peak.append(mass)
            one_peak.append(times[int(peak_properties[1]["left_ips"][element] + (peak_properties[1]["widths"][element] / 2))])
            one_peak.append(peak_properties[1]["peak_heights"][element])
            one_peak.append(times[int(peak_properties[1]["right_ips"][element])] - times[int(peak_properties[1]["left_ips"][element])])
            one_peak.append(peak_properties[1]["prominences"][element])
            identified_peaks.append(one_peak)
        # [[mass, time, height, duration, prominence], [....], ...]
        return identified_peaks

    def combine_and_sum_dicts(self, dict1, dict2):
        return {k: dict1.get(k, 0) + dict2.get(k, 0) for k in set(dict1) | set(dict2)}

    def get_best_approx_for_ppm_spacing_within_peak(self, mass_list, worst_expected_ppm_deviation=20):
        mass_list = sorted(mass_list)
        ppm_spacing_list = [(((mass_list[i + 1] - mass_list[i]) / mass_list[i]) * 1000000) for i in range(len(mass_list) - 1)]
        ppm_spacing_list = [ppm for ppm in ppm_spacing_list if ppm < worst_expected_ppm_deviation]
        avg_ppm = sum(ppm_spacing_list) / len(ppm_spacing_list)
        return avg_ppm

    def summarize_mass_intensity_dict_with_deviation(self, dictio, deviation=11):
        print("=============================================================")
        print("=============================================================")
        start_time = datetime.datetime.now()
        print("Start time: " + str(datetime.datetime.now()))
        print("summarizing dict according to new method")
        dictio = {k: v for k, v in dictio.items() if v >= 1}
        # sort the dictio by its keys
        dictio = dict(sorted(dictio.items(), key=lambda item: item[0]))
        print("old length of start dictio:")
        print(len(dictio))

        old_masses_list = list(dictio.keys())
        old_abundances_list = list(dictio.values())
        new_masses_list = []
        new_abundances_list = []

        best_approx_ppm_spacing_within_peak = self.get_best_approx_for_ppm_spacing_within_peak(old_masses_list)
        print("best approx for ppm spacing within peak: " + str(best_approx_ppm_spacing_within_peak))

        while len(old_masses_list) > 0:
            try:
                remove_all_lower = False
                remove_all_upper = False
                index_of_highest_abundance = old_abundances_list.index(max(old_abundances_list))
                curr_mass = old_masses_list[index_of_highest_abundance]
                mass_lower_border = curr_mass - ((deviation * curr_mass) / 1000000)
                mass_upper_border = curr_mass + ((deviation * curr_mass) / 1000000)
                iteration_step_lower = 0
                integration_step_upper = 0
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                expected_min_spacing_between_measurement_points = (best_approx_ppm_spacing_within_peak * curr_mass) / 1000000
                try:
                    while (mass_lower_border < old_masses_list[index_of_highest_abundance - iteration_step_lower]) and \
                            (last_existing_mass_within_border - old_masses_list[index_of_highest_abundance - iteration_step_lower] <= 2.2 * expected_min_spacing_between_measurement_points) and \
                            (old_abundances_list[index_of_highest_abundance - iteration_step_lower] <= (1.1 * last_abundance)):
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance - iteration_step_lower]
                        last_abundance = old_abundances_list[index_of_highest_abundance - iteration_step_lower]
                        iteration_step_lower += 1
                except IndexError:
                    remove_all_lower = True
                    iteration_step_lower = 0
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                try:
                    while (mass_upper_border > old_masses_list[index_of_highest_abundance + integration_step_upper]) and \
                            (old_masses_list[index_of_highest_abundance + integration_step_upper] - last_existing_mass_within_border <= 2.2 * expected_min_spacing_between_measurement_points) and \
                            (old_abundances_list[index_of_highest_abundance + integration_step_upper] <= (1.1 * last_abundance)):
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance + integration_step_upper]
                        last_abundance = old_abundances_list[index_of_highest_abundance + integration_step_upper]
                        integration_step_upper += 1
                except IndexError:
                    remove_all_upper = True
                    integration_step_upper = 0

                if remove_all_lower == True and remove_all_upper == True:
                    break
                if remove_all_upper:
                    summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, len(old_abundances_list), 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, len(old_abundances_list), 1)]) / summed_intensity
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, len(old_masses_list), 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, len(old_abundances_list), 1)]
                    continue
                if remove_all_lower:
                    summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance + integration_step_upper, 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance + integration_step_upper, 1)]) / summed_intensity
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance + integration_step_upper, 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance + integration_step_upper, 1)]
                    continue

                if iteration_step_lower == 0 and integration_step_upper == 0:
                    summed_intensity = old_abundances_list[index_of_highest_abundance]
                    weighted_mass_average = old_masses_list[index_of_highest_abundance]
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance]
                    continue

                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + integration_step_upper, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + integration_step_upper, 1)]) / summed_intensity
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + integration_step_upper, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + integration_step_upper, 1)]
            except Exception as e:
                print("EXCEPTION IN summarize_mass_intensity_dict_with_deviation()!!!")
                print(traceback.format_exc())
                print(e)
                break
        outdict = dict(zip(new_masses_list, new_abundances_list))
        print("new length of summarized dictio:")
        print(len(outdict))
        print("End time: " + str(datetime.datetime.now()))
        print("Time taken = " + str(datetime.datetime.now() - start_time))
        return outdict

    def get_all_available_masses_in_list(self):
        ordered_mass_intensity_dict = {}
        for entry in self.ms_file.all_ms_spectra:
            if not self.ms_file.all_modes[self.ms_file.all_ms_spectra.index(entry)] == "Full scan":
                continue
            if self.ms_file.all_ms_spectra.index(entry) > 300:
                break
            print(str(self.ms_file.all_ms_spectra.index(entry)) + " / " + str(len(self.ms_file.all_ms_spectra)))
            curr_save = entry
            curr_save = {round(key, self.kwargs["round_for_masses_summary"]): value for key, value in curr_save.items() if value > self.kwargs["minimum_required_max_intensity"] / 10}
            ordered_mass_intensity_dict = self.combine_and_sum_dicts(ordered_mass_intensity_dict, curr_save)
            ordered_mass_intensity_dict = {round(key, self.kwargs["round_for_masses_summary"]): value for key, value in ordered_mass_intensity_dict.items() if value > self.kwargs["minimum_required_max_intensity"] / 10}
        print("ORDERED MASS INTENSITY DICT____FIRST STEP:")
        print(len(ordered_mass_intensity_dict))
        print(ordered_mass_intensity_dict)

        ordered_mass_intensity_dict = dict(sorted(ordered_mass_intensity_dict.items(), key=lambda x: x[1], reverse=True))
        ordered_mass_intensity_dict = self.summarize_mass_intensity_dict_with_deviation(ordered_mass_intensity_dict, self.kwargs["mass_deviation"])
        print("ORDERED MASS INTENSITY DICT____SECOND STEP:")
        print(len(ordered_mass_intensity_dict))
        print(ordered_mass_intensity_dict)

        ordered_mass_intensity_dict = dict(sorted(ordered_mass_intensity_dict.items(), key=lambda x: x[1], reverse=True))
        ordered_mass_intensity_dict = {key: value for key, value in ordered_mass_intensity_dict.items() if value > self.kwargs["minimum_required_max_intensity"]}
        print("ORDERED MASS INTENSITY DICT____LAST STEP:")
        print(len(ordered_mass_intensity_dict))
        print(ordered_mass_intensity_dict)
        return ordered_mass_intensity_dict


def baseline_als(y, lam=1e6, p=0.01, niter=10, window_min_vals=4):
    """Asymmetric Least Squares (ALS) baseline correction with correct matrix sizing."""
    # Parameters:
    # y: input data (1D array)
    # lam: smoothness parameter (larger values give smoother baseline)
    # p: asymmetry parameter (0 < p < 1, larger values give more weight to the left side)
    # niter: number of iterations for convergence
    # Returns:
    from scipy.sparse import diags
    import numpy as np
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
    D = diags([1, -2, 1], [0, -1, -2], shape=(L, L)).toarray()
    w = np.ones(L)
    for _ in range(niter):
        W = np.diag(w)
        # Fix: Use the same size for D as y
        Z = W + lam * (D.T @ D)
        # Solve using np.linalg.solve
        baseline = np.linalg.solve(Z, w * y)
        # Update weights (asymmetry)
        #w = np.where(y < baseline, 1.0, p)  # Fix: Use np.where for element-wise comparison
        w = p * (y > baseline) + (1 - p) * (y < baseline)
    for i in range(len(baseline)):
        if baseline[i] <= 0:
            baseline[i] = 0
    for i in range(len(baseline)):
        if y_original[i] - baseline[i] < 0:
            baseline[i] = y_original[i]
    return baseline


def get_window_size_by_frequency(intensityvals, timevals, min_width=3, max_width=30):
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
            dom_freqs.append(freqs[np.argmax(mags[4:])])
    if len(dom_freqs) > 0:
        #print("dom_freqs: " + str(dom_freqs))
        avg_dom_freq = sum(dom_freqs) / len(dom_freqs) # 1/s
        window_size = int(1 / avg_dom_freq) # 1 oscillation every x seconds
        window_size = window_size * 2 # Increase the window size by a factor of 2 to get better smoothing results --> a window size of 1 oscillation is not enough!
        window_size = int(window_size * fs) # in measurements --> Window size is so that the window contains 1 oscillation and depending on the sampling rate x measurements
        if window_size < min_width:
            window_size = min_width
        elif window_size > max_width:
            window_size = max_width
    else:
        window_size = (max_width-min_width) / 2
    return window_size


def do_smoothing_without_effecting_peaks(y, window_size=5):
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


def get_one_mass_range(ms_file, mass, mass_window, threshold_area, threshold_intensity, mass_deviation=30):
    print("Mass: " + str(mass))
    starttime = datetime.datetime.now()
    xic = MS_functions.get_xic(ms_file.rawdata, mass=mass, mass_deviation=mass_window *1.5, requested_filter_mode="Full scan")
    print("XIC created: " + str(datetime.datetime.now() - starttime))
    times = xic[0]
    intensities = xic[1]
    true_indices_of_entries = xic[2]
    # Do first smoothing
    window_size = get_window_size_by_frequency(intensities, times, min_width=2, max_width=30)
    intensityvals = do_smoothing_without_effecting_peaks(intensities, window_size=window_size)
    # Apply baseline correction
    baseline = baseline_als(intensityvals, lam=1e4, p=0.05, niter=100, window_min_vals=window_size)
    corrected_intensity = intensityvals - baseline
    smoothed_intensity = do_smoothing_without_effecting_peaks(corrected_intensity, window_size=window_size)
    print("Smoothing and bg subst done: " + str(datetime.datetime.now() - starttime))
    # Find peaks with adaptive height and width detection
    min_width_seconds = int( (len(times)/(max(times)-min(times)) ) * 1)
    max_width_seconds = int( (len(times)/(max(times)-min(times)) ) * 50)
    min_width_measurements = 2
    max_width_measurements = 30
    min_width = min_width_seconds if min_width_seconds > min_width_measurements else min_width_measurements
    max_width = max_width_seconds if max_width_seconds > max_width_measurements else max_width_measurements
    peaks, properties = scipy.signal.find_peaks(smoothed_intensity,
                                prominence=np.std(smoothed_intensity), width=(min_width, max_width))
    #Calculate peak area
    peak_areas = []
    for peak in range(len(peaks)):
        left = int(properties["left_ips"][peak])
        right = int(properties["right_ips"][peak])
        print(peak)
        area = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
        if area <= threshold_area:
            continue
        spec = UVenture.Spec(ms_file, true_indices_of_entries[peaks[peak]], spec_requested_filter_mode="Full scan", mass_deviation=mass_deviation)
        interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(mass_window*1.5) and k < mass+(mass_window*1.5)}
        print("interesting_range: " + str(interesting_range))
        interesting_range = {k: v for k, v in interesting_range.items() if v > threshold_intensity}
        for k, v in interesting_range.items():
            print("Entry: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
            peak_areas.append((k, times[peaks[peak]], v, area))
    print(peak_areas)
    return peak_areas



if __name__ == "__main__":
    fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(5, 15), dpi=300, layout="tight")

    print("Start")
    mzml_filename = "V://005 Mitarbeiter-aktuell/Karbach/UVenture-main/webserver_save/mzml_files/05Dec2024_PAMPain2_neg_20241008_ISP_OH_O3_NO.mzML"
    ms_file = UVenture.MS_File(mzml_filename)
    min_mz, max_mz = ms_file.mz_range

    mass_range = 1
    threshold_area = 40000
    threshold_intensity = 10000
    identified_peaks = []
    for i in range(int(185), int(183), -mass_range):
        peak_areas = get_one_mass_range(ms_file, i, mass_range, threshold_area, threshold_intensity)
        identified_peaks.extend(peak_areas)
    peak_df = pd.DataFrame(identified_peaks, columns=["mass", "rt", "height", "area"])
    print(peak_df)
    peak_df["unique_descriptor"] = peak_df["mass"].round(4).astype(str) + peak_df["rt"].round(5).astype(str)
    peak_df.drop_duplicates(subset=["unique_descriptor"], inplace=True)
    print(peak_df)
    peak_df["peak_bins"] = pd.cut(peak_df["rt"], bins=100, labels=False)
    print(peak_df)
    peak_df = peak_df.sort_values(by=["rt"], ascending=[True])
    final_peaks = {}
    for i in identified_peaks:
        descriptor = str(round(i[0], 3)) + str(round(i[1], 1))
        if descriptor not in list(final_peaks.keys()):
            final_peaks[descriptor] = i
    final_peaks = list(final_peaks.values())
    final_peaks = sorted(final_peaks, key=operator.itemgetter(0), reverse=True)
    final_peaks = sorted(final_peaks, key=operator.itemgetter(2), reverse=True)
    final_peaks = sorted(final_peaks, key=operator.itemgetter(1))
    for v in final_peaks:
        print("Final peak: " + str(v))
    
    rt_dict = {}
    for i in final_peaks:
        if round(i[1], 0) not in list(rt_dict.keys()):
            rt_dict[round(i[1], 0)] = []
        rt_dict[round(i[1], 0)].append(v)
    for k, v in rt_dict.items():
        masses_at_rt = v



    #xic = MS_functions.get_xic(ms_file.rawdata, mass=183.10248, mass_deviation=2, requested_filter_mode="Full scan")
    xic = MS_functions.get_xic(ms_file.rawdata, mass=186.097, mass_deviation=1, requested_filter_mode="Full scan")
    #xic = MS_functions.get_xic(ms_file.rawdata, mass=138.0196, mass_deviation=.001, requested_filter_mode="Full scan")
    timevals = xic[0]
    intensityvals = xic[1]
    true_indices_of_entries = xic[2]
    df = pd.DataFrame({"time": timevals, "intensity": intensityvals})
    ax[0].plot(df["time"], df["intensity"], label="Original Chromatogram", alpha=0.5, color="red")

    window_size = get_window_size_by_frequency(intensityvals, timevals, min_width=2, max_width=30)
    print("window_size: " + str(window_size))
    intensityvals = do_smoothing_without_effecting_peaks(intensityvals, window_size=window_size)
    ax[0].plot(df["time"], intensityvals, label="Smoothed Signal", color="blue")
    
    # Apply baseline correction
    baseline = baseline_als(intensityvals, lam=1e4, p=0.05, niter=100, window_min_vals=window_size)
    corrected_intensity = intensityvals - baseline
    ax[1].plot(df["time"], intensityvals, label="Smoothed Signal", color="blue")
    ax[1].plot(df["time"], corrected_intensity, label="Corrected Signal", color="green")
    ax[1].plot(df["time"], baseline, label="Estimated Baseline", linestyle="dashed", color="red")

    smoothed_intensity = do_smoothing_without_effecting_peaks(corrected_intensity, window_size=window_size)
    ax[2].plot(df["time"], smoothed_intensity, label="Smooth", color="blue")
    
    from scipy.signal import find_peaks
    # Find peaks with adaptive height and width detection
    min_width_seconds = int( (len(timevals)/(max(timevals)-min(timevals)) ) * 1)
    max_width_seconds = int( (len(timevals)/(max(timevals)-min(timevals)) ) * 50)
    min_width_measurements = 2
    max_width_measurements = 30
    min_width = min_width_seconds if min_width_seconds > min_width_measurements else min_width_measurements
    max_width = max_width_seconds if max_width_seconds > max_width_measurements else max_width_measurements
    peaks, properties = find_peaks(smoothed_intensity,
                                   prominence=np.std(smoothed_intensity), width=(min_width, max_width))
    ax[2].plot(df["time"], smoothed_intensity, label="Smooth", color="blue")
    ax[2].plot(df["time"][peaks], smoothed_intensity[peaks], "rx", label="Detected Peaks")
    
    #print(peaks)
    #curr_peak = peaks[2]
    #spec = UVenture.Spec(ms_file, true_indices_of_entries[curr_peak], spec_requested_filter_mode="Full scan")
    #interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > 186.1 - 1 and k < 186.1 + 1}
    #print("interesting_range: " + str(interesting_range))
    #x = ax[4].bar(list(interesting_range.keys()), list(interesting_range.values()), label="MS Spectrum", color="blue", width=0.01)
    #ax[4].bar_label(x, label_type="edge")
    #ax[4].text(1,1, "RT" + str(spec.rt) + " min", fontsize=8, ha='right', va='top')
    #Calculate peak area
    peak_areas = []
    for peak in range(len(peaks)):
        #print("Peak " + str(peak) + ": " + str(properties))
        left = int(properties["left_ips"][peak])
        right = int(properties["right_ips"][peak])
        area = np.trapz(smoothed_intensity[left:right], dx=(timevals[1] - timevals[0]))
        peak_areas.append(area)
        print("Peak area for peak " + str(peak) + ": " + str(area))
        ax[3].plot(df["time"][left:right], smoothed_intensity[left:right], label="Peak Area", alpha=0.5, color="orange")
        ax[3].fill_between(df["time"][left:right], smoothed_intensity[left:right], alpha=0.5, color="orange")
    ax[3].plot(df["time"], smoothed_intensity, label="Smooth", color="blue")
    ax[3].plot(df["time"][peaks], smoothed_intensity[peaks], "rx", label="Detected Peaks")




    ax[3].set_xlabel("Time (s)")
    ax[3].set_ylabel("Intensity")
    ax[3].legend()
    plt.tight_layout()
    plt.savefig("test" + str(datetime.datetime.now().strftime("%Y%m%d%H%M%S")) + ".png", dpi=300)
    print("End")
    

