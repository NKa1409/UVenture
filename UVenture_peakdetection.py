import ast
import copy
import datetime
import math
import pathlib
import sys
import traceback
import PIL
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import matplotlib
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


def baseline_als(y, lam=1e6, p=0.01, niter=10):
    """Asymmetric Least Squares (ALS) baseline correction with correct matrix sizing."""
    from scipy.sparse import diags
    import numpy as np
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
        w = p * (y > baseline) + (1 - p) * (y < baseline)
    for i in range(len(baseline)):
        if baseline[i] <= 0:
            baseline[i] = 0
    return baseline


if __name__ == "__main__":
    print("Start")
    mzml_filename = "U://MyFolder//MONOTONS//GithubUVenture//UVenture-20240520_Tested1//webserver_save//mzml_files//05Dec2024_PAMPain2_neg_20241008_ISP_OH_O3_NO.mzML"
    ms_file = UVenture.MS_File(mzml_filename)
    print(ms_file)
    xic = MS_functions.get_xic(ms_file.rawdata, mass=138.0191, mass_deviation=.001, requested_filter_mode="Full scan")
    timevals = xic[0]
    intensityvals = xic[1]
    df = pd.DataFrame({"time": timevals, "intensity": intensityvals})
    print(df)

    # Apply baseline correction
    baseline = baseline_als(intensityvals, lam=1e4, p=0.05)
    corrected_intensity = intensityvals - baseline

    from scipy.signal import savgol_filter
    smoothed_intensity = savgol_filter(corrected_intensity, window_length=11, polyorder=2)

    from scipy.signal import find_peaks
    # Find peaks with adaptive height and width detection
    peaks, properties = find_peaks(smoothed_intensity, height=np.mean(smoothed_intensity) + 2 * np.std(smoothed_intensity),
                                   prominence=0.5 * np.std(smoothed_intensity), width=(5, 50))


    # Plot baseline-corrected chromatogram
    plt.figure(figsize=(10, 5))
    plt.plot(timevals, intensityvals, label="Original Chromatogram", alpha=0.5, color="black")
    plt.plot(timevals, baseline, label="Estimated Baseline", linestyle="dashed", color="red")
    plt.plot(timevals, corrected_intensity, label="Corrected Signal", color="green")
    plt.plot(timevals, smoothed_intensity, label="Smooth", color="blue")
    #plt.plot(timevals[peaks], smoothed_intensity[peaks], "rx", label="Detected Peaks")
    plt.xlabel("Time (min)")
    plt.ylabel("Intensity")
    plt.title("Baseline Correction")
    plt.legend()
    plt.savefig("test.png")

