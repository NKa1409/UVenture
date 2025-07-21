# MS_functions.py
import sys
import traceback
import similaritymeasures
import scipy
import numpy as np
from pyteomics import mzml
import pandas as pd



def read_mzml_file(mzml_filename):
    f = mzml.read(mzml_filename)
    f = list(f)
    return f


def get_xic(f, mass, mass_deviation, requested_filter_mode="Full scan"):
    # The mass deviation is defined as the requested mass +1x the mass deviation and -1x the mass deviation.
    # If the requested mass is 1000 and the mass deviation is 0.005, the range is 999.995 to 1000.005.
    if (not requested_filter_mode == "Full scan") and (not requested_filter_mode == "AIF") and (not requested_filter_mode == "MS/MS"):
        requested_filter_mode = "Full scan"
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    intensity_list = []
    index_with_wrong_ms_level_list = []
    for entry in f:
        filter = entry["scanList"]["scan"][0]["filter string"]
        filter_mode = ""
        if (" d " in filter) and ("@hcd" in filter):
            filter_mode = "MS/MS"
        if (not " d " in filter) and (not "hcd" in filter):
            filter_mode = "Full scan"
        if (not " d " in filter) and ("hcd" in filter):
            filter_mode = "AIF"
        if not filter_mode == requested_filter_mode:
            index_with_wrong_ms_level_list.append(f.index(entry))
        mass_indices = (i for i in range(len((entry["m/z array"]))) if mass-mass_deviation <= entry["m/z array"][i] <= mass+mass_deviation)
        try:
            curr_sum_int = 0
            for mass_index in mass_indices:
                curr_sum_int = curr_sum_int + entry["intensity array"][mass_index]
            intensity_list.append(curr_sum_int)
        except:
            intensity_list.append(0)
    original_index_list = list(range(len(f)))
    if not requested_filter_mode == "all":
        for index in sorted(index_with_wrong_ms_level_list, reverse=True):
            del rt_list[index]
            del intensity_list[index]
            del original_index_list[index]
    return [rt_list, intensity_list, original_index_list]


def compare_peak_shape_similarity(xic1, xic2, peak_rt, peakwidth=10, debug_output=False):
    index = xic1[0].index(min(xic1[0], key=lambda x: abs(peak_rt - x)))
    intensity1_at_peak_rt = xic1[1][index]
    if intensity1_at_peak_rt <= 0:
        intensity1_at_peak_rt = 0.000001
    intensity2_at_peak_rt = xic2[1][index]
    if intensity2_at_peak_rt <= 0:
        intensity2_at_peak_rt = 0.000001
    try:
        peakintensity1 = xic1[1][(index - peakwidth):(index + peakwidth)]
        peak1_rt = xic1[0][int(index - peakwidth):int(index + peakwidth)]
        peakintensity2 = xic2[1][int(index - peakwidth):int(index + peakwidth)]
        peak2_rt = xic2[0][int(index - peakwidth):int(index + peakwidth)]

        neighbour_list1 = xic1[1][int(index - 3 * peakwidth):int(index + 3 * peakwidth)]
        neighbour_list1 = [neighbour_list1[i] for i in range(len(neighbour_list1)) if (i < len(neighbour_list1) / 3) or (i > ((len(neighbour_list1) / 3) + (len(neighbour_list1) / 2)))]
        average_surrounding1 = sum(neighbour_list1) / len(neighbour_list1)
        try:
            peakintensity1 = [((i - average_surrounding1) / (xic1[1][index] - average_surrounding1)) for i in peakintensity1]
        except:
            max_peakint1 = max(peakintensity1)
            if max_peakint1 == 0:
                max_peakint1 = 0.001
            peakintensity1 = [((i - average_surrounding1) / (max_peakint1)) for i in peakintensity1]

        neighbour_list2 = xic2[1][int(index - 3 * peakwidth):int(index + 3 * peakwidth)]
        neighbour_list2 = [neighbour_list2[i] for i in range(len(neighbour_list2)) if (i < len(neighbour_list2) / 3) or (i > ((len(neighbour_list2) / 3) + (len(neighbour_list2) / 2)))]
        average_surrounding2 = sum(neighbour_list2) / len(neighbour_list2)
        try:
            peakintensity2 = [((i - average_surrounding2) / (xic2[1][index] - average_surrounding2)) for i in peakintensity2]
        except:
            max_peakint2 = max(peakintensity2)
            if max_peakint2 == 0:
                max_peakint2 = 0.0001
            peakintensity2 = [[((i - average_surrounding2) / (max_peakint2)) for i in peakintensity2]]
    except:
        if index - peakwidth <= 2:
            peakwidth = index - 2
        if index + peakwidth >= len(xic1[1]):
            peakwidth = (len(xic1[1]) - index - 2)
        peakintensity1 = xic1[1][(index - peakwidth):(index + peakwidth)]
        peakintensity1 = [i / intensity1_at_peak_rt for i in peakintensity1]
        peak1_rt = xic1[0][int(index - peakwidth):int(index + peakwidth)]
        peakintensity2 = xic2[1][int(index - peakwidth):int(index + peakwidth)]
        peakintensity2 = [i / intensity2_at_peak_rt for i in peakintensity2]
        peak2_rt = xic2[0][int(index - peakwidth):int(index + peakwidth)]

    try:
        P = np.array([peak1_rt, peakintensity1]).T
        Q = np.array([peak2_rt, peakintensity2]).T
        area = similaritymeasures.area_between_two_curves(P, Q)
    except:
        area = -1

    if not (isinstance(area, float) or isinstance(area, int)):
        try:
            area = float(area)
        except:
            area = -1
    if area == np.nan or (str(area).lower() == "nan"):
        area = -1
        if debug_output == True:
            print("area was nan. Chaning area to: " + str(area))

    return area, peak1_rt, peakintensity1, peak2_rt, peakintensity2


def get_mode_of_spec(filter_string):
    if " d " in filter_string and "@hcd" in filter_string:
        ms_ms_masses = []
        filter_parsed = filter_string.split(" ")
        filter_parsed = [x for x in filter_parsed if "hcd" in x]
        for element in filter_parsed:
            ms_ms_masses.append(round(float(element.split("@")[0]), 2))
        return "MS/MS"
    elif " d " not in filter_string and "hcd" not in filter_string:
        return "Full scan"
    elif " d " not in filter_string and "hcd" in filter_string:
        return "AIF"


def min_deviation_between_list_elements(input_list):
    # Sort the list in ascending order
    sorted_list = sorted(input_list)
    # Initialize the minimum deviation with a large value
    min_deviation = float('inf')
    # Initialize the pair of values where the minimum deviation occurs
    min_deviation_values = None
    # Iterate over the sorted list and compare adjacent elements
    for i in range(len(sorted_list) - 1):
        deviation = abs(sorted_list[i + 1] - sorted_list[i])
        if deviation < min_deviation:
            min_deviation = deviation
            min_deviation_values = (sorted_list[i], sorted_list[i + 1])
    return [min_deviation, min_deviation_values]


def get_best_approx_for_ppm_spacing_within_peak(mass_list, worst_expected_ppm_deviation=20):
    mass_list = sorted(mass_list)
    ppm_spacing_list = [(((mass_list[i+1] - mass_list[i]) / mass_list[i]) * 1000000) for i in range(len(mass_list)-1)]
    ppm_spacing_list = [ppm for ppm in ppm_spacing_list if ppm < worst_expected_ppm_deviation]
    avg_ppm = sum(ppm_spacing_list) / len(ppm_spacing_list)
    return avg_ppm

def get_peaks_in_xy_series(x, y, sg_window=10, sg_order=3):
        if len(x) != len(y):
            print("Error in MS_functions.get_peaks_in_xy_series(): x and y have different lengths!")
            return None
        y = scipy.signal.savgol_filter(y, sg_window, sg_order, mode="nearest")
        peak_properties = scipy.signal.find_peaks(y, height=max(y)/100, distance=2, prominence=max(y)/100, width=(2, len(y)/10))
        identified_peaks = []
        for element in range(len(peak_properties[1]["peak_heights"])):
            one_peak = []
            one_peak.append(x[int(peak_properties[1]["left_ips"][element] + (peak_properties[1]["widths"][element]/2))])
            one_peak.append(peak_properties[1]["peak_heights"][element])
            one_peak.append(x[int(peak_properties[1]["left_ips"][element])])
            one_peak.append(x[int(peak_properties[1]["right_ips"][element])])
            one_peak.append(peak_properties[1]["prominences"][element])
            identified_peaks.append(one_peak)
        #[[time, height, lefttime, righttime, prominence], [time, height, left, right, prominence], ...]
        return identified_peaks


def summarize_mass_intensity_dict(dictio, deviation=11, debug_output=True):
        dictio = {k: v for k, v in dictio.items() if v > 0}
        if debug_output == True:
            print("summarizing dict according to new method. old length of start dictio:" + str(len(dictio)))
        #sort the dictio by its keys
        dictio = dict(sorted(dictio.items(), key=lambda item: item[0]))
        old_masses_list = list(dictio.keys())
        old_abundances_list = list(dictio.values())
        new_masses_list = []
        new_abundances_list = []

        best_approx_ppm_spacing_within_peak = get_best_approx_for_ppm_spacing_within_peak(old_masses_list)
        if debug_output == True:
            print("best approx for ppm spacing within peak: " + str(best_approx_ppm_spacing_within_peak))

        while len(old_masses_list) > 0:
            try:
                remove_all_lower = False
                remove_all_upper = False
                index_of_highest_abundance = old_abundances_list.index(max(old_abundances_list))
                curr_mass = old_masses_list[index_of_highest_abundance]
                mass_lower_border = curr_mass - 4*((deviation*curr_mass)/1000000)
                mass_upper_border = curr_mass + 4*((deviation*curr_mass)/1000000)
                iteration_step_lower = 0
                integration_step_upper = 0

                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                expected_min_spacing_between_measurement_points = (best_approx_ppm_spacing_within_peak * curr_mass) / 1000000
                try:
                    while (mass_lower_border < old_masses_list[index_of_highest_abundance-(iteration_step_lower+1)]) and \
                            (last_existing_mass_within_border - old_masses_list[index_of_highest_abundance-(iteration_step_lower+1)] <= 3.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance-(iteration_step_lower+1)] <= (1.1 * last_abundance)):

                        iteration_step_lower += 1
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance - iteration_step_lower]
                        last_abundance = old_abundances_list[index_of_highest_abundance - iteration_step_lower]
                except IndexError:
                    remove_all_lower = True
                    iteration_step_lower = 0

                
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                try:
                    while (mass_upper_border > old_masses_list[index_of_highest_abundance+integration_step_upper+1]) and \
                            (old_masses_list[index_of_highest_abundance+integration_step_upper+1] - last_existing_mass_within_border <= 3.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance+integration_step_upper+1] <= (1.1 * last_abundance)):
                        integration_step_upper += 1
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance + integration_step_upper]
                        last_abundance = old_abundances_list[index_of_highest_abundance + integration_step_upper]
                except IndexError:
                    remove_all_upper = True
                    integration_step_upper = 0
                
                if remove_all_lower == True and remove_all_upper == True:
                    break
                if remove_all_upper:
                    summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]) / summed_intensity
                    if summed_intensity >= 1:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_masses_list), 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]
                    continue
                if remove_all_lower:
                    summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper+1, 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]) / summed_intensity
                    if summed_intensity >= 1:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]
                    continue

                if iteration_step_lower == 0 and integration_step_upper == 0:
                    summed_intensity = old_abundances_list[index_of_highest_abundance]
                    weighted_mass_average = old_masses_list[index_of_highest_abundance]
                    if summed_intensity >= 1:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance ]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance ]
                    continue

                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]) / summed_intensity
                if summed_intensity >= 1:
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]
            except Exception as e:
                print("EXCEPTION IN summarize_mass_intensity_dict()!!!")
                print(traceback.format_exc())
                break
        outdict = dict(zip(new_masses_list, new_abundances_list))
        if debug_output == True:
            print("New length of summarized dictio: " + str(len(outdict)))
        return outdict



def summarize_mass_intensity_dict_for_isotopo_simulation(dictio, deviation=30, debug_output=True):
    dictio = dict(sorted(dictio.items(), key=lambda item: item[0]))
    old_masses_list = list(dictio.keys())
    old_abundances_list = list(dictio.values())
    new_masses_list = []
    new_abundances_list = []

    while len(old_masses_list) > 0:
        try:
            remove_all_lower = False
            remove_all_upper = False
            index_of_highest_abundance = old_abundances_list.index(max(old_abundances_list))
            curr_mass = old_masses_list[index_of_highest_abundance]
            mass_lower_border = curr_mass - ((deviation * curr_mass) / 1000000)
            mass_upper_border = curr_mass + ((deviation * curr_mass) / 1000000)

            iteration_step_lower = 0
            iteration_step_upper = 0

            try:
                while (mass_lower_border < old_masses_list[index_of_highest_abundance - (iteration_step_lower + 1)]):
                    iteration_step_lower += 1
            except IndexError:
                remove_all_lower = True
                iteration_step_lower = 0

            try:
                while (mass_upper_border > old_masses_list[index_of_highest_abundance + iteration_step_upper + 1]):
                    iteration_step_upper += 1
            except IndexError:
                remove_all_upper = True
                iteration_step_upper = 0

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
                summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)]) / summed_intensity
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)]
                continue

            if iteration_step_lower == 0 and iteration_step_upper == 0:
                summed_intensity = old_abundances_list[index_of_highest_abundance]
                weighted_mass_average = old_masses_list[index_of_highest_abundance]
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance]
                continue

            summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)])
            weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)]) / summed_intensity
            new_masses_list.append(weighted_mass_average)
            new_abundances_list.append(summed_intensity)
            old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)]
            old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)]
        except Exception as e:
            print("EXCEPTION IN summarize_mass_intensity_dict_for_isotopo_simulation()!!!")
            print(traceback.format_exc())
            break
    outdict = dict(zip(new_masses_list, new_abundances_list))
    return outdict


def read_summary_to_df(filepath, sep="\t", colnames=("MASS", "FORMULA", "INDEX", "RT", "INTENSITY", "SCORE", "FRAGS", "F_SCORES", "F_INTENSITIES", "NLS", "NLS_DEV"), header=None):
    df = pd.read_csv(filepath, sep=sep, header=header)
    new_colnames = []
    for entry in range(len(df.columns)):
        try:
            new_colnames.append(colnames[entry])
        except Exception as e:
            print("Error setting columnames in read_summary_to_df(): " + str(e))
            print(traceback.format_exc())
            new_colnames.append(entry)
    df.columns = new_colnames
    return df






if __name__ == "__main__":
    print("This module is not meant to be run directly. Please import it in your script.")
    sys.exit()

