import statistics
import traceback

import pyteomics
from pyteomics import mzml
import copy
import matplotlib.pyplot as plt
import scipy
import time
from pyteomics import mass
import math
import ast
import os
import datetime
import PIL
import matplotlib
import pathlib
import multiprocessing
import sys
import time
import psutil
import inspect
import MS_functions
import multiprocessing
import sys
import time
import psutil
import inspect


class MyMultiprocessing:
    def __init__(self, only_multithread=False):
        self.only_multithread = only_multithread
        self.functions = []
        self.args = []
        self.pool = None
        self.pool_list = []
        self.pool_list_starmap = []
        self.ret_values = []
        self.max_cores = psutil.cpu_count() - 1
        self.timeout = 0.02
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        self.caller = mod.__name__
        self.caller_file_location = mod.__file__
        if self.caller == "__main__":
            pass
        else:
            print("You need to put the class into the following if-statement: \n"
                  "if __name__ == '__main__':\n"
                  "    mp = MultiprocessingClass.MyMultiprocessing()\n"
                  "    mp.add_function(fn, arglist)\n"
                  "    mp.run_pool()\n"
                  "    ret = mp.get_return_values()")
            time.sleep(1)
            sys.exit()
        if self.only_multithread == False:
            self.pool = multiprocessing.Pool(processes=self.max_cores)
        else:
            from multiprocessing.pool import ThreadPool
            self.pool = ThreadPool(processes=self.max_cores)
            print("onlyMultithread")

    def add_function(self, function, arg_list):
        if hasattr(function, "__self__"):
            print("Function is bound method!\n"
                  "The function you want to add is not static. If the function is inside a class, consider adding the "
                  "'@staticmethod' decorator.\n"
                  "Otherwise use a function that is not inside the class "
                  "(that is not relying on the 'self' argument). \n"
                  "Only functions that do not have the 'self' argument are valid functions to use with this "
                  "multiprocessing module.\n"
                  "The given function was therefore NOT added!")
            return None
        self.functions.append(function)

        if arg_list is None:
            print("arg_list must not be empty!\n"
                  "At least one argument is needed in the function for the class to work properly!!!")
        if arg_list == []:
            print("arg_list must not be empty!\n"
                  "At least one argument is needed in the function for the class to work properly!!!")
        if arg_list == ():
            print("arg_list must not be empty!\n"
                  "At least one argument is needed in the function for the class to work properly!!!")
        if isinstance(arg_list, list):
            pass
        elif isinstance(arg_list, int):
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        elif isinstance(arg_list, str):
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        elif isinstance(arg_list, float):
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        else:
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        self.args.append(arg_list)

    def get_return_values(self):
        start_time = time.time()
        while len(self.ret_values) < len(self.functions):
            for pool in list(self.pool_list):
                try:
                    self.ret_values.append(pool.get(timeout=self.timeout))
                    self.pool_list.remove(pool)
                except:
                    pass
            for pool in list(self.pool_list_starmap):
                try:
                    self.ret_values.append(pool.get(timeout=self.timeout))
                    self.pool_list_starmap.remove(pool)
                except:
                    pass
            if start_time < (time.time() - self.timeout):
                break
        return self.ret_values

    def run_pool(self, wait_for_finish=True):
        for function in self.functions:
            self.pool_list.append(self.pool.apply_async(func=function, args=tuple(self.args[len(self.pool_list)])))
            self.pool_list_starmap.append(
                self.pool.starmap_async(function, tuple(self.args[len(self.pool_list_starmap)])))
        if wait_for_finish:
            self.pool.close()
            self.pool.join()
            start_time = time.time()
            while len(self.ret_values) < len(self.functions):
                for pool in list(self.pool_list):
                    try:
                        self.ret_values.append(pool.get(timeout=self.timeout))
                        self.pool_list.remove(pool)
                    except:
                        pass
                for pool in list(self.pool_list_starmap):
                    try:
                        self.ret_values.append(pool.get(timeout=self.timeout))
                        self.pool_list_starmap.remove(pool)
                    except:
                        pass
                if start_time < (time.time() - self.timeout):
                    break
        else:
            pass


def do_bg_substraction(y, frequency_in_seconds=1, prominence=400000, distance=3, min_height=0, width=(0, 0),
                       sg_windows=(6, 50, 2), sg_orders=(5, 5, 1), not_including_peak_width_multiplier=3):
    if min_height == 0:
        min_height = max(y) / 10
    else:
        pass
    approx_capture_duration_per_element = frequency_in_seconds
    if width == (0, 0):
        width = [int(5 / approx_capture_duration_per_element), int(20 / approx_capture_duration_per_element)]
    else:
        width = width

    peaks, peak_properties = scipy.signal.find_peaks(scipy.signal.savgol_filter(y, sg_windows[0], sg_orders[0]),
                                                     height=min_height, prominence=prominence, distance=distance,
                                                     width=width)

    background = copy.deepcopy(y)
    for peak in range(len(peaks)):
        for index in range(int(peaks[peak] - peak_properties["widths"][peak] * not_including_peak_width_multiplier),
                           int(peaks[peak] + peak_properties["widths"][peak] * not_including_peak_width_multiplier), 1):

            try:
                m = (background[
                         int(peaks[peak] + peak_properties["widths"][peak] * not_including_peak_width_multiplier)] -
                     background[
                         int(peaks[peak] - peak_properties["widths"][peak] * not_including_peak_width_multiplier)]) / (
                            int(
                                peaks[peak] + peak_properties["widths"][
                                    peak] * not_including_peak_width_multiplier) - int(
                        peaks[peak] - peak_properties["widths"][peak] * not_including_peak_width_multiplier))
            except:
                m = 1
            try:
                background[index] = background[int(
                    peaks[peak] - peak_properties["widths"][peak] * not_including_peak_width_multiplier)] + (m * (
                        index - int(
                    peaks[peak] - peak_properties["widths"][peak] * not_including_peak_width_multiplier)))
            except:
                break
    background = scipy.signal.savgol_filter(background, sg_windows[1], sg_orders[1])
    bg_subst_series = []
    for entry in range(len(y)):
        bg_subst_series.append(y[entry] - background[entry])
    bg_subst_series = scipy.signal.savgol_filter(bg_subst_series, sg_windows[2], sg_orders[2])
    return [bg_subst_series, background, peaks, peak_properties]


def check_for_isotopologue_peaks_in_mass_spectrum(spectrum,
                                                  formula_to_check,
                                                  max_ppm_deviation=20,
                                                  charge_of_formula_to_check=-1,
                                                  noise_divisor_for_isotopologue_calculation=2,
                                                  multiplier_isotopologue_influence_of_ppm_deviation_on_score=1,
                                                  isotopologue_score_e_function_exponent=0.2,
                                                  minimum_assumed_noise=1000,
                                                  add_value_to_score_if_isotopo_was_found=100,
                                                  stop_isotopologue_search_if_score_lower_than=-200,
                                                  max_ppm_deviation_change_for_isotopologue=1.3,
                                                  summarized_mass_intensity_dict=None):
    round_masses_to_decimalplace = 6  # not needed to adjust. BUT LEAVE THIS IN THE FUNCTION!! Leftover variable from previous method...

    ordered_mass_intensities_dict = {}
    for entry in range(len(spectrum[0])):
        ordered_mass_intensities_dict[spectrum[0][entry]] = spectrum[1][entry]
    ordered_mass_intensities_dict = dict(
        sorted(ordered_mass_intensities_dict.items(), key=lambda item: item[1], reverse=True))
    ordered_mass_intensities_dict = {m: i for m, i in ordered_mass_intensities_dict.items() if i > 1}
    all_masses_intensities_dict = copy.deepcopy(ordered_mass_intensities_dict)

    isotope_pattern_dict = MS_functions.simulate_isotope_pattern_of_formula(formula_to_check)
    isotope_pattern_list = []
    for item in list(isotope_pattern_dict.items()):
        isotope_pattern_list.append([item[0], item[1]])
    isotope_pattern_dict = {}
    for element in isotope_pattern_list:
        try:
            isotope_pattern_dict[round(element[0], round_masses_to_decimalplace)] = isotope_pattern_dict[
                                                                                        round(element[0],
                                                                                              round_masses_to_decimalplace)] + \
                                                                                    element[1]
        except:
            isotope_pattern_dict[round(element[0], round_masses_to_decimalplace)] = element[1]
    isotope_pattern_dict = dict(sorted(isotope_pattern_dict.items(), key=lambda item: item[1], reverse=True))

    isotope_pattern_dict = {abs((m - 0.000548 * charge_of_formula_to_check) / charge_of_formula_to_check): a for m, a in
                            isotope_pattern_dict.items()}

    mass = min(spectrum[0], key=lambda x: abs(list(isotope_pattern_dict.items())[0][0] - x))
    within_deviation_mass_list = [masse for masse in list(all_masses_intensities_dict.keys()) if ((abs(
        list(isotope_pattern_dict.items())[0][0] - masse) / list(isotope_pattern_dict.items())[0][
                                                                                                       0]) * 1000000) < max_ppm_deviation]
    intensity_of_requested_mass = 0
    for masse in within_deviation_mass_list:
        intensity_of_requested_mass = intensity_of_requested_mass + all_masses_intensities_dict[masse]

    summarized_masses_intensities_dict = summarized_mass_intensity_dict

    isotopes_found = {}
    previous_score = 99999
    previous_deviation_list = []
    startitem = list(isotope_pattern_dict.items())[0]
    isotopo_mass = round(startitem[0], round_masses_to_decimalplace)
    measured_mass_with_minimal_deviation = min(list(summarized_masses_intensities_dict.keys()),
                                               key=lambda x: abs(isotopo_mass - x))
    initial_deviation = (abs(isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000

    for isotopologue in list(isotope_pattern_dict.items()):
        isotopo_mass = round(isotopologue[0], round_masses_to_decimalplace)
        within_deviation_mass_list = [masse for masse in list(summarized_masses_intensities_dict.keys()) if abs((((
                                                                                                                   isotopo_mass - masse) / isotopo_mass) * 1000000) - initial_deviation) < max_ppm_deviation_change_for_isotopologue]
        measured_mass_with_minimal_deviation = min(list(summarized_masses_intensities_dict.keys()), key=lambda x: abs(
            (((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
        deviation = (abs(isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000

        if (deviation < max_ppm_deviation) and (
                abs(deviation - initial_deviation) <= max_ppm_deviation_change_for_isotopologue):
            previous_deviation_list.append(copy.deepcopy(deviation))
            theoretical_intensity = (intensity_of_requested_mass / (list(isotope_pattern_dict.items())[0][1])) * \
                                    isotopologue[1]
            noise = sorted(list(summarized_masses_intensities_dict.values()))[int(len(list(
                summarized_masses_intensities_dict.items())) / noise_divisor_for_isotopologue_calculation)] + minimum_assumed_noise
            measured_intensity = 0
            for masse in within_deviation_mass_list:
                measured_intensity = measured_intensity + summarized_masses_intensities_dict[masse]
            if measured_intensity <= 0:
                break
            isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)] = [True, theoretical_intensity,
                                                                                    measured_intensity]
            score = measured_intensity / theoretical_intensity  # je naeher an 1 desto besser; wenn <1: weniger gemessen als da sein sollte; wenn >1: mehr gemessen als da sein sollte.
            if score > 1:
                score = -(-1.5 + (1 / (0.7 + (2.718281828459045 ** (-isotopologue_score_e_function_exponent * score)))))
            else:
                score = score * 1
            # score: je naeher an 1 desto besser; wenn negativ: weniger gemessen als theoretisch da; wenn positiv: mehr gemessen als theoretisch da.
            score = score * 100
            score = score + 0.001
            score = abs(score)
            score = score + add_value_to_score_if_isotopo_was_found
            if measured_intensity < noise:
                score_multiplier_by_intensity_and_noise = abs(1 / (math.log(measured_intensity / noise) - 1))
            else:
                score_multiplier_by_intensity_and_noise = abs(math.log(measured_intensity / noise) + 1)

            if isotopologue == list(isotope_pattern_dict.items())[0]:
                isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)].append(0)
                previous_score = 9999999
                continue
            if theoretical_intensity <= (noise / 1.5):
                score = score * score_multiplier_by_intensity_and_noise
                score = score - abs(
                    deviation * multiplier_isotopologue_influence_of_ppm_deviation_on_score * score_multiplier_by_intensity_and_noise)
                isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)].append(score)
                if (score - abs(stop_isotopologue_search_if_score_lower_than) >= previous_score):
                    isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)][-1] = \
                        isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)][-1] * 0.1
                break
            else:
                score = score * score_multiplier_by_intensity_and_noise
                score = score - abs(
                    deviation * multiplier_isotopologue_influence_of_ppm_deviation_on_score * score_multiplier_by_intensity_and_noise)
                isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)].append(score)
            if (score < stop_isotopologue_search_if_score_lower_than):
                break
            if (score - abs(stop_isotopologue_search_if_score_lower_than) >= previous_score):
                isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)][-1] = \
                    isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)][-1] * 0.1
                break
            previous_score = score
        else:
            theoretical_intensity = (intensity_of_requested_mass / (list(isotope_pattern_dict.items())[0][1])) * \
                                    isotopologue[1]
            noise = sorted(list(all_masses_intensities_dict.values()))[
                int(len(list(all_masses_intensities_dict.items())) / 2)]
            measured_intensity = 0
            isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)] = [False, theoretical_intensity,
                                                                                    measured_intensity]
            negative_score = theoretical_intensity / noise
            negative_score = -10 * negative_score
            isotopes_found[round(isotopologue[0], round_masses_to_decimalplace)].append(negative_score)
            break
    return [isotopes_found, formula_to_check, isotope_pattern_dict, intensity_of_requested_mass]


def get_spectra_associated_with_mass(f, mass, mass_deviation=0.005,
                                     requested_xic_mode="Full scan",
                                     sg_windows_bg_removal=(10, 75, 2), width_bg_removal=(1, 15),
                                     min_distance_of_chromatographic_peaks_in_seconds=3,
                                     prominence_for_peak_detection=40000, min_relative_height_for_peak_detection=0.1,
                                     remove_if_combined_fragment_score_lower_than=-10, charge_of_measured_mass=-1,
                                     remove_fragments_with_score_lower_than=-10, max_nl_deviation_ppm=7.2,
                                     divisor_for_fragment_noise=3,
                                     max_ppm_deviation=20, ppm_deviation_scoring_multiplier=3,
                                     neighbour_window=10,
                                     include_peak_if_height_multiple_of_std_dev=2,
                                     noise_divisor_for_isotopologue_calculation=2,
                                     multiplier_isotopologue_influence_of_ppm_deviation_on_score=1,
                                     isotopologue_score_e_function_exponent=0.2,
                                     summary_filepath="",
                                     nl_formula_cache_dict=None,
                                     nl_formula_cache_dict_max_mass=230,
                                     fragment_formula_cache_dict=None,
                                     fragment_formula_cache_dict_max_mass=200,
                                     reject_predicted_formula_if_score_lower_than=10,
                                     minimum_assumed_noise=1000,
                                     minimal_required_peak_height_to_background_ratio=0.8,
                                     add_value_to_score_if_isotopo_was_found=100,
                                     max_ppm_deviation_find_formula=25,
                                     stop_isotopologue_search_if_score_lower_than=-200,
                                     min_peakarea=100000):
    def set_size(w, h, ax=None):
        """ w, h: width, height in inches """
        if not ax: ax = plt.gca()
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        figw = float(w) / (r - l)
        figh = float(h) / (t - b)
        ax.figure.set_size_inches(figw, figh)

    current_folder = output_folder + str(round(mass, 4)) + "//"
    if not os.path.isdir(current_folder):
        os.makedirs(current_folder)

    xic = MS_functions.get_xic(f, mass, mass_deviation, requested_filter_mode=requested_xic_mode)
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    approx_capture_duration_per_element = (max(rt_list) * 60) / len(xic[1])

    bg_subst_int = do_bg_substraction(xic[1], approx_capture_duration_per_element, sg_windows=sg_windows_bg_removal,
                                      width=width_bg_removal,
                                      distance=((
                                                        min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element) + 1))

    peaks = scipy.signal.find_peaks(bg_subst_int[0],
                                    height=max(bg_subst_int[0]) * min_relative_height_for_peak_detection,
                                    prominence=prominence_for_peak_detection,
                                    distance=((
                                                      min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element) + 1),
                                    width=(2, 30))

    # dynamic peak finding.
    # get std dev of xic in direct neighbourhood of the peak (excluding the peak itself)
    peaks_original = copy.deepcopy(peaks)
    intensities_in_neighbourhood = []
    new_peak_list = []
    for peak_index in peaks[0]:
        try:
            intensities_in_neighbourhood_1 = bg_subst_int[0][
                                             peak_index - neighbour_window:int(peak_index - neighbour_window / 2)]
            intensities_in_neighbourhood_2 = bg_subst_int[0][
                                             int(peak_index + neighbour_window / 2):peak_index + neighbour_window]
            for element in intensities_in_neighbourhood_1:
                intensities_in_neighbourhood.append(element)
            for element in intensities_in_neighbourhood_2:
                intensities_in_neighbourhood.append(element)
            mean = sum(intensities_in_neighbourhood) / len(intensities_in_neighbourhood)
            variance = sum([((x - mean) ** 2) for x in intensities_in_neighbourhood]) / len(
                intensities_in_neighbourhood)
            res = variance ** 0.5
        except:
            res = 0
        if (bg_subst_int[0][peak_index] >= (include_peak_if_height_multiple_of_std_dev * res)) and (
                bg_subst_int[0][peak_index] >= minimal_required_peak_height_to_background_ratio * bg_subst_int[1][
            peak_index]):
            new_peak_list.append(peak_index)
    peaks = list(peaks)
    # calculate area of peak and delete peaks below the threshhold limit.
    peak_areas = []
    for peak in range(len(peaks_original[0])):
        width_of_peak = list(peaks_original[1]["widths"])[peak]
        index_of_peak = peaks_original[0][peak]
        area = 0
        for element in range(int(index_of_peak - (width_of_peak / 2)), int(index_of_peak + (width_of_peak + 2)), 1):
            try:
                area = area + xic[1][element]
            except Exception as e:
                print("Error in calculating area of peak: " + str(e))
                print(traceback.format_exc())
                break
        peak_areas.append(area)
    new_peak_list2 = [peaks[0][p_i] for p_i in range(len(peak_areas)) if peak_areas[p_i] >= min_peakarea]
    new_peak_list = list(set(new_peak_list).intersection(new_peak_list2))
    new_peak_list.sort()

    keep_list_indices = [l_i for l_i in list(range(len(peaks[0]))) if peaks[0][l_i] in new_peak_list]
    peaks[0] = copy.deepcopy(new_peak_list)
    for key, value in peaks[1].items():
        peaks[1][key] = [value[v] for v in range(len(value)) if v in keep_list_indices]
    fig, ax = plt.subplots()
    for peak_index in peaks[0]:
        ax.scatter(xic[0][peak_index], bg_subst_int[0][peak_index], color="red", marker="x", s=100)
    ax.plot(xic[0], xic[1], label="XIC measured")
    ax.plot(xic[0], bg_subst_int[1], label="Background")
    ax.plot(xic[0], bg_subst_int[0], label="XIC background_substracted")
    ax.set_title("xic_" + str(round(mass, 4)) + "+-" + str(mass_deviation) + "_" + str(requested_xic_mode))
    ax.set_xlabel("retention time / min")
    ax.set_ylabel("intensity / a.u.")
    plt.legend()
    matplotlib.rcParams.update({'figure.autolayout': True})
    image_filepath = current_folder + "xic_" + str(round(mass, 4)) + "+-" + str(mass_deviation) + "_" + str(
        requested_xic_mode) + ".png"
    plt.savefig(image_filepath, bbox_inches='tight')
    metadata = PIL.PngImagePlugin.PngInfo()
    metadata.add_text("xic_retention_time", str(xic[0]))
    metadata.add_text("xic_intensity_list", str(xic[1]))
    metadata.add_text("xic_original_index_list", str(xic[2]))
    metadata.add_text("xic_background", str(bg_subst_int[1]))
    metadata.add_text("xic_background_subst_intensity", str(bg_subst_int[0]))
    target_image = PIL.Image.open(image_filepath)
    target_image.save(image_filepath, pnginfo=metadata)

    log_f = open(current_folder + "identified_peaks_for_mass.txt", "w")
    log_f.write("Number of peaks\t" + str(len(peaks[0])) + "\n")
    log_f.write("\n")
    for peak_index in peaks[0]:
        log_f.write("ID\t" + str(list(peaks[0]).index(peak_index)) + "\t Original peak index\t" + str(
            xic[2][peak_index]) + "\t RT \t" + str(round(rt_list[xic[2][peak_index]], 2)) + "\t Intensity \t" + str(
            round(xic[1][peak_index], 1)) + "\n")
    log_f.close()

    if len(peaks[0]) == 0:
        print("ERROR: NO PEAKS WERE FOUND!!!!!!!!!!")

    best_formula_for_each_index_list = []

    results_list = []
    for peak_index in peaks[0]:
        try:
            original_peak_index = xic[2][peak_index]
            specific_folder = current_folder + str(original_peak_index) + "//"
            if not os.path.isdir(specific_folder):
                os.makedirs(specific_folder)

            full_scan_spec = MS_functions.get_mass_spectrum(f, original_peak_index, mode="as_index",
                                                            requested_filter_mode="Full scan",
                                                            save_folder=specific_folder)
            spec = MS_functions.get_mass_spectrum(f, original_peak_index, mode="as_index",
                                                  requested_filter_mode="MS/MS",
                                                  requested_ms_ms_mass=mass, save_folder=specific_folder)

            mass_spectra_list = []
            spec = MS_functions.get_mass_spectrum(f, original_peak_index, mode="as_index", requested_filter_mode="AIF",
                                                  save_folder=specific_folder)
            aif_processed = process_aif_spectrum(spec,
                                                 mass,
                                                 save_folder=specific_folder,
                                                 charge_of_measured_mass=charge_of_measured_mass,
                                                 remove_fragments_with_score_lower_than=remove_fragments_with_score_lower_than,
                                                 max_nl_deviation_ppm=max_nl_deviation_ppm,
                                                 divisor_for_fragment_noise=divisor_for_fragment_noise,
                                                 max_ppm_deviation=max_ppm_deviation,
                                                 ppm_deviation_scoring_multiplier=ppm_deviation_scoring_multiplier,
                                                 noise_divisor_for_isotopologue_calculation=noise_divisor_for_isotopologue_calculation,
                                                 multiplier_isotopologue_influence_of_ppm_deviation_on_score=multiplier_isotopologue_influence_of_ppm_deviation_on_score,
                                                 isotopologue_score_e_function_exponent=isotopologue_score_e_function_exponent,
                                                 full_scan_spec=full_scan_spec,
                                                 nl_formula_cache_dict=nl_formula_cache_dict,
                                                 nl_formula_cache_dict_max_mass=nl_formula_cache_dict_max_mass,
                                                 fragment_formula_cache_dict=fragment_formula_cache_dict,
                                                 fragment_formula_cache_dict_max_mass=fragment_formula_cache_dict_max_mass,
                                                 reject_predicted_formula_if_score_lower_than=reject_predicted_formula_if_score_lower_than,
                                                 minimum_assumed_noise=minimum_assumed_noise,
                                                 add_value_to_score_if_isotopo_was_found=add_value_to_score_if_isotopo_was_found,
                                                 max_ppm_deviation_find_formula=max_ppm_deviation_find_formula,
                                                 stop_isotopologue_search_if_score_lower_than=stop_isotopologue_search_if_score_lower_than)
            if aif_processed == None:
                print("Peak is not further processed, as process_aif_spectrum returned None!")
                continue
        except Exception as e_process_aif:
            print("ERROR IN PROCESSING AIF SPECTRUM! Error code:")
            print(str(e_process_aif))
            print(traceback.format_exc())
            continue
        try:
            true_fragment_list = aif_processed[0]
            molecule_predictions_list = aif_processed[5]
        except Exception as e:
            print("Error while processing AIF spectrum! Error code: ")
            print(str(e))
            print(traceback.format_exc())
            continue

        fig, ax = plt.subplots()
        ax.bar(spec[0], spec[1], color="gray", label="Spectrum", alpha=0.2)
        ax.bar(aif_processed[3], aif_processed[4], color="red", width=1.5, label="Molecule")
        delete_index_list = []
        combined_score_list = []
        for fragment in range(len(true_fragment_list)):
            combined_score = true_fragment_list[fragment][0][2] + true_fragment_list[fragment][1][2] - abs(
                true_fragment_list[fragment][2][1])
            if combined_score < remove_if_combined_fragment_score_lower_than:
                delete_index_list.append(fragment)
            else:
                combined_score_list.append(combined_score)
        for i in sorted(delete_index_list, reverse=True):
            del true_fragment_list[i]
        sorted_index_list = [combined_score_list.index(i) for i in sorted(combined_score_list, reverse=True)]
        true_fragment_list = [true_fragment_list[i] for i in sorted_index_list]
        fragment_mass_list = []
        fragment_intensity_list = []
        nl_mass_list = []
        nl_formula_list = []
        fragment_formula_list = []
        molecule_formula_list = []
        scoring_list = []
        fragment_score_list = []
        nl_deviation_list = []
        molecule_score_list = []
        text_str = ""
        text_str = text_str + "{:<10}|{:<10}|{:<10}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format("score",
                                                                                                          "Frag mass",
                                                                                                          "intensity",
                                                                                                          "Frag formula",
                                                                                                          "NL formula",
                                                                                                          "Molec formula",
                                                                                                          "score F",
                                                                                                          "dev. NL",
                                                                                                          "score M")
        text_str = text_str + "{:_<10}|{:_<10}|{:_<10}|{:_<15}_{:_<15}___{:_<15}|{:_>10}|{:_>10}|{:_>10}\n".format("",
                                                                                                                   "",
                                                                                                                   "",
                                                                                                                   "",
                                                                                                                   "",
                                                                                                                   "",
                                                                                                                   "",
                                                                                                                   "",
                                                                                                                   "")
        for molecule in molecule_predictions_list:
            text_str = text_str + \
                       "{:<10}|{:<10}|{:<10}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format(
                           str(round(float(molecule[1]), 2)),
                           0,
                           "100%",
                           "",
                           "",
                           molecule[0],
                           0,
                           0,
                           str(round(float(molecule[1]), 2)))
        print(text_str)
        for fragment in true_fragment_list:
            combined_score = fragment[0][2] + fragment[1][2] - abs(fragment[2][1])
            text_str = text_str + \
                       "{:<10}|{:<10}|{:<10}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format(
                           round(combined_score, 2), str(round(fragment[1][3], 4)),
                           str(round(((fragment[1][1] / aif_processed[4]) * 100), 2)) + "%",
                           str("".join([str(a) + str(n) for a, n in fragment[1][0].items()])),
                           str("".join([str(a) + str(n) for a, n in fragment[2][0].items()])),
                           str("".join([str(a) + str(n) for a, n in fragment[0][0].items()])), round(fragment[1][2], 2),
                           round(abs(fragment[2][1]), 2), round(fragment[0][2], 2))
            fragment_mass_list.append(fragment[1][3])
            fragment_intensity_list.append(fragment[1][1])
            nl_mass_list.append(fragment[2][2])
            nl_formula_list.append(fragment[2][0])
            fragment_formula_list.append(fragment[1][0])
            molecule_formula_list.append(fragment[0][0])
            fragment_score_list.append(fragment[1][2])
            nl_deviation_list.append(fragment[2][1])
            molecule_score_list.append(fragment[0][2])
            combined_score = fragment[0][2] + fragment[1][2] - abs(fragment[2][1])
            scoring_list.append(combined_score)
        molecule_formula_dict = {}
        for formula in range(len(molecule_formula_list)):
            formula_string = str("".join([str(a) + str(n) for a, n in molecule_formula_list[formula].items()]))
            if formula_string not in list(molecule_formula_dict.keys()):
                molecule_formula_dict[formula_string] = [[fragment_mass_list[formula]],
                                                         [fragment_intensity_list[formula]]]
            else:
                molecule_formula_dict[formula_string][0].append(fragment_mass_list[formula])
                molecule_formula_dict[formula_string][1].append(fragment_intensity_list[formula])
        for molecule_formula_item in list(molecule_formula_dict.items()):
            ax.bar(molecule_formula_item[1][0], molecule_formula_item[1][1], width=1.4,
                   label="Fragments for: " + str(molecule_formula_item[0]))
        text_str = text_str.strip()
        props = dict(boxstyle='round', facecolor='grey', alpha=0.05)  # bbox features
        text = ax.text(1.03, 0.98, text_str, fontfamily="monospace", transform=ax.transAxes, fontsize=10,
                       verticalalignment="top", bbox=props)
        text_width = text.get_window_extent(renderer=fig.canvas.get_renderer()).width
        text_height = text.get_window_extent(renderer=fig.canvas.get_renderer()).height
        set_size(4 + (text_width / fig.dpi), 4 + (text_height / fig.dpi), ax)

        filter_mode = spec[6]
        index = spec[2]
        filter = spec[5]
        ms_ms_masses = spec[7]
        ax.legend(loc="upper left")
        ax.set_xlabel("masses / Da")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("Mode:" + str(filter_mode) + "; Index: " + str(index) + "; RT: " + str(
            round(rt_list[index], 2)) + " min" + ";\nFilter: " + str(filter) + ";\nMS/MS masses: " + str(ms_ms_masses))
        matplotlib.rcParams.update({'figure.autolayout': True})
        image_filepath = specific_folder + "Mass_spectrum_FRAGMENTS_index" + str(index) + "_" + str(
            filter_mode) + "_" + str(ms_ms_masses) + ".png"
        plt.savefig(image_filepath, bbox_inches='tight')
        metadata = PIL.PngImagePlugin.PngInfo()
        metadata.add_text("masses", str(spec[0]))
        metadata.add_text("fragment_masses", str(fragment_mass_list))
        metadata.add_text("intensities", str(spec[1]))
        metadata.add_text("fragment_intensities", str(fragment_intensity_list))
        metadata.add_text("molecule_mass", str(aif_processed[3]))
        metadata.add_text("molecule_intensities", str(aif_processed[4]))
        metadata.add_text("index", str(index))
        metadata.add_text("orig_file_entry_for_index", str(f[index]))
        metadata.add_text("mode", str(filter_mode))
        metadata.add_text("filter", str(filter))
        metadata.add_text("msms_masses", str(ms_ms_masses))
        target_image = PIL.Image.open(image_filepath)
        target_image.save(image_filepath, pnginfo=metadata)

        total_combined_score_per_formula_prediction = {}
        for entry in range(len(molecule_formula_list)):
            formula_string = "".join([str(a) + str(n) for a, n in molecule_formula_list[entry].items()])
            try:
                total_combined_score_per_formula_prediction[formula_string] = \
                    total_combined_score_per_formula_prediction[formula_string] + scoring_list[entry]
            except:
                total_combined_score_per_formula_prediction[formula_string] = scoring_list[entry]
        fragments_for_best_formula = []
        score_of_fragments = []
        nl_for_best_formula = []
        score_of_nl = []
        frag_intensity_for_best_formula = []
        score_of_molecule_list = []

        for entry in range(len(molecule_formula_list)):
            formula_string = "".join([str(a) + str(n) for a, n in molecule_formula_list[entry].items()])
            fragment_string = "".join([str(a) + str(n) for a, n in fragment_formula_list[entry].items()])
            nl_string = "".join([str(a) + str(n) for a, n in nl_formula_list[entry].items()])
            if formula_string == max(total_combined_score_per_formula_prediction,
                                     key=total_combined_score_per_formula_prediction.get):
                fragments_for_best_formula.append(fragment_string)
                nl_for_best_formula.append(nl_string)
                score_of_fragments.append(fragment_score_list[entry])
                score_of_nl.append(nl_deviation_list[entry])
                frag_intensity_for_best_formula.append(fragment_intensity_list[entry])
                score_of_molecule_list.append(molecule_score_list[entry])

        best_formula = []
        try:
            if len(molecule_formula_list) == 0:
                for molecule in molecule_predictions_list:
                    best_formula.append(molecule[0])
                    best_formula.append(molecule[1])
                    best_formula.append(mass)
                    best_formula.append(index)
                    best_formula.append(aif_processed[4])  # intensity of molecule
                    best_formula.append(fragments_for_best_formula)
                    best_formula.append(frag_intensity_for_best_formula)
                    best_formula.append(nl_for_best_formula)
                    best_formula.append(score_of_fragments)
                    best_formula.append(score_of_nl)
                    best_formula.append(score_of_molecule_list)

            else:
                best_formula.append(max(total_combined_score_per_formula_prediction,
                                        key=total_combined_score_per_formula_prediction.get))
                best_formula.append(total_combined_score_per_formula_prediction[
                                        max(total_combined_score_per_formula_prediction,
                                            key=total_combined_score_per_formula_prediction.get)])
                best_formula.append(mass)
                best_formula.append(index)
                best_formula.append(aif_processed[4])  # intensity of molecule
                best_formula.append(fragments_for_best_formula)
                best_formula.append(frag_intensity_for_best_formula)
                best_formula.append(nl_for_best_formula)
                best_formula.append(score_of_fragments)
                best_formula.append(score_of_nl)
                best_formula.append(score_of_molecule_list)
        except Exception as e_best_formula:
            print("Error in best_formula.append part of get_spectra_associated_with_mass. Error: ")
            print(str(e_best_formula))
            print(traceback.format_exc())
            continue

        if len(best_formula) >= 1:
            with open(specific_folder + "BEST_FORMULA_APPROXIMATION.txt", "w") as file:
                file.write("{:<30}".format("Best Formula: ") + str(best_formula[0]) + "\n")
                file.write("{:<30}".format("Score: ") + str(round(best_formula[1], 2)) + "\n")
                file.write("{:<30}".format("Intensity of Molecule: ") + str(round(best_formula[4], 1)) + "\n")
                file.write("_____________________________________________________________\n")
                file.write("{:<10}".format("Mass: ") + str(round(best_formula[2], 4)) + "\n")
                file.write("{:<10}".format("Index: ") + str(best_formula[3]) + "\n")
                file.write("_____________________________________________________________\n")
                file.write("===============FRAGMENTS===============\n")
                file.write(
                    "{:<20}{:<20}{:<20}{:<20}{:<20}\n".format("Fragment", "Score", "NL", "NL Deviation", "Intensity"))
                for fragment_index in range(len(best_formula[5])):
                    file.write("{:<20}".format(str(best_formula[5][fragment_index])))
                    file.write("{:<20}".format(str(round(best_formula[8][fragment_index], 2))))
                    file.write("{:<20}".format(str(best_formula[7][fragment_index])))
                    file.write("{:<20}".format(str(round(best_formula[9][fragment_index], 2))))
                    file.write("{:<20}\n".format(str(round(best_formula[6][fragment_index], 2))))

            if not summary_filepath == "":
                abs_path = os.path.abspath(summary_filepath)
                os.makedirs(os.path.dirname(abs_path), exist_ok=True)
                with open(summary_filepath, "a") as summary_file:
                    # mass \t formula \t index \t rt \t intensity \t score \t [fragments] \t [fragments_score] \t [fragments_intensity] \t [nl] \t [nl_deviation]
                    summary_file.write(str(round(best_formula[2], 4)) + "\t")
                    summary_file.write(str(best_formula[0]) + "\t")
                    summary_file.write(str(best_formula[3]) + "\t")
                    summary_file.write(str(round(rt_list[index], 2)) + "\t")
                    summary_file.write(str(round(best_formula[4], 1)) + "\t")
                    summary_file.write(str(round(best_formula[1], 4)) + "\t")
                    summary_file.write(str(best_formula[5]) + "\t")
                    summary_file.write(str(best_formula[8]) + "\t")
                    summary_file.write(str(best_formula[6]) + "\t")
                    summary_file.write(str(best_formula[7]) + "\t")
                    summary_file.write(str(best_formula[9]) + "\n")

        best_formula_for_each_index_list.append(best_formula)
        results_list.append(aif_processed)

    for best_formula in best_formula_for_each_index_list:
        if len(best_formula) == 0:
            continue
        with open(current_folder + "BEST_FORMULA_APPROXIMATIONS.txt", "a") as file:
            file.write(
                "=======================================================================================================\n")
            file.write(
                "=======================================================================================================\n")
            file.write("{:<30}".format("Best Formula: ") + str(best_formula[0]) + "\n")
            file.write("{:<30}".format("Score: ") + str(round(best_formula[1], 2)) + "\n")
            file.write("{:<30}".format("Intensity of Molecule: ") + str(round(best_formula[4], 2)) + "\n")
            file.write("_____________________________________________________________\n")
            file.write("{:<10}".format("Mass: ") + str(round(best_formula[2], 4)) + "\n")
            file.write("{:<10}".format("Index: ") + str(best_formula[3]) + "\n")
            file.write("_____________________________________________________________\n")
            file.write("===============FRAGMENTS===============\n")
            file.write(
                "{:<20}{:<20}{:<20}{:<20}{:<20}\n".format("Fragment", "Score", "NL", "NL Deviation", "Intensity"))
            for fragment_index in range(len(best_formula[5])):
                file.write("{:<20}".format(str(best_formula[5][fragment_index])))
                file.write("{:<20}".format(str(round(best_formula[8][fragment_index], 2))))
                file.write("{:<20}".format(str(best_formula[7][fragment_index])))
                file.write("{:<20}".format(str(round(best_formula[9][fragment_index], 2))))
                file.write("{:<20}\n".format(str(round(best_formula[6][fragment_index], 2))))
    return [best_formula_for_each_index_list, results_list]


def predict_formula_of_mass_with_mass_spectrum(spectrum, mass,
                                               ppm_deviation_scoring_multiplier=3,
                                               save_folder="",
                                               charge_of_measured_mass=-1,
                                               max_ppm_deviation=20,
                                               noise_divisor_for_isotopologue_calculation=2,
                                               multiplier_isotopologue_influence_of_ppm_deviation_on_score=1,
                                               isotopologue_score_e_function_exponent=0.2,
                                               fragment_formula_cache_dict=None,
                                               fragment_formula_cache_dict_max_mass=200,
                                               reject_predicted_formula_if_score_lower_than=10,
                                               minimum_assumed_noise=1000,
                                               add_value_to_score_if_isotopo_was_found=100,
                                               stop_isotopologue_search_if_score_lower_than=-200,
                                               include_likelyhood_of_formula=True,
                                               summarized_mass_intensity_dict=None):
    # create formula prediction for the desired molecule
    log_f_filepath = ""
    if not save_folder == "":
        log_f_filepath = save_folder + "predictFormulaForMassWithSpectrum_mass_" + str(round(mass, 4)) + ".txt"
        log_f = open(log_f_filepath, "w")
        log_f.write("Predicting formula for mass of: " + str(round(mass, 4)) + "\n")
        log_f.write("max_deviation_ppm \t" + str(max_ppm_deviation) + "\n")
        log_f.write("ppm_deviation_scoring_multiplier \t" + str(ppm_deviation_scoring_multiplier) + "\n")

        log_f.close()

    best_formulas_molecule = MS_functions.get_formula_from_cache(
        "U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//", mass, max_ppm_deviation)

    # delete all formula predictions of molecule from dict where deviation is too high
    # best_formulas_molecule = [mass, {formula1: ppm_deviation1, formula2: ppm_deviation2, .....}]
    max_molecule_mass_deviation = max_ppm_deviation
    for formula in list(best_formulas_molecule[1].keys()):
        try:
            if abs(best_formulas_molecule[1][formula]) >= max_molecule_mass_deviation:
                del best_formulas_molecule[1][formula]
                continue
            if ((charge_of_measured_mass < 0) and (("Na" in formula) or ("K" in formula))):
                del best_formulas_molecule[1][formula]
                continue
        except Exception as e:
            print("Error with formulaprediction exclusion: " + str(e))
            pass

    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("================================================\n")
        log_f.write("OUTPUT OF SIMPLE FORMULA PREDICTION (ONLY MASS DEVIATION): \n")
        log_f.write("{:>20} \t {:>20} \n".format("Formula", "ppm_deviation"))
        for item in list(best_formulas_molecule[1].items()):
            log_f.write("{:>20} \t {:>20} \n".format(str(item[0]), str(round(item[1], 2))))
        log_f.close()

    # create list of all possible formulas in dict format
    # gives a list with [[formula_dict, ppm_deviation], ...]
    molecule_formulas = []
    for entry in best_formulas_molecule[1]:
        formula_dict_molecule = MS_functions.get_formula_to_dict(entry)
        molecule_formulas.append([formula_dict_molecule, best_formulas_molecule[1][entry]])

    # simulate isotope pattern to check if the formula prediction might be correct.
    # isotopo_check = [{simulated_isotopologue_mass1: [bool_isotopolog_found, theoretical_intensity_of_isotopolo, measured_intensity_of_isotopolo, score_of_isotopologue], ......},
    #                 formula_prediction_that_was_checked,
    #                 {simulated_isotopologue_mass1: simulated_isotopologue_abundance1, ....},
    #                 intensity_of_requested_mass]
    formula_score_dict = {}
    # formula_score_dict = {formula_prediction1: combined_score1, formula_prediction2: combined_score2, ......}
    for entry in molecule_formulas:
        formula_prediction = entry[0]
        isotope_check = check_for_isotopologue_peaks_in_mass_spectrum(spectrum,
                                                                      formula_prediction,
                                                                      max_ppm_deviation=max_ppm_deviation,
                                                                      charge_of_formula_to_check=charge_of_measured_mass,
                                                                      noise_divisor_for_isotopologue_calculation=noise_divisor_for_isotopologue_calculation,
                                                                      multiplier_isotopologue_influence_of_ppm_deviation_on_score=multiplier_isotopologue_influence_of_ppm_deviation_on_score,
                                                                      isotopologue_score_e_function_exponent=isotopologue_score_e_function_exponent,
                                                                      minimum_assumed_noise=minimum_assumed_noise,
                                                                      add_value_to_score_if_isotopo_was_found=add_value_to_score_if_isotopo_was_found,
                                                                      stop_isotopologue_search_if_score_lower_than=stop_isotopologue_search_if_score_lower_than,
                                                                      summarized_mass_intensity_dict=summarized_mass_intensity_dict)

        formula_score = 0
        for element in list(isotope_check[0].items()):
            try:
                formula_score = formula_score + element[1][3]
            except Exception as e:
                print("Error in calculating complete formula score! Code: " + str(e) + ", however it will continue!!!")
                continue
        formula_score = formula_score - abs(ppm_deviation_scoring_multiplier * entry[1])
        formula_score_dict[str(formula_prediction)] = formula_score

        if not save_folder == "":
            log_f = open(log_f_filepath, "a")
            log_f.write("\n")
            log_f.write("================================================\n")
            log_f.write("OUTPUT OF FORMULA PREDICTION WITH MASS SPECTRUM (ISOTOPOLOGUES)\n")
            log_f.write("FORMULA: " + str("".join([str(a) + str(n) for a, n in formula_prediction.items()])) + "\n")
            log_f.write("COMBINED SCORE: " + str(round(formula_score, 2)) + "\n")
            log_f.write("PPM DEVIATION OF MOLECULE FORMULA: " + str(round(entry[1], 2)) + "\n")
            log_f.write(
                "{:>25} \t {:>25} \t {:>25} \t {:>25} \t {:>25} \n".format("Mass_of_isotopologue", "Isotopo_found?",
                                                                           "measured_intensity", "required_intensity",
                                                                           "score_of_isotopologue"))
            for item in list(isotope_check[0].items()):
                try:
                    log_f.write(
                        "{:>25} \t {:>25} \t {:>25} \t {:>25} \t {:>25} \n".format(str(round(item[0], 4)),
                                                                                   str(item[1][0]),
                                                                                   str(round(item[1][2], 4)),
                                                                                   str(round(item[1][1], 4)),
                                                                                   str(round(item[1][3], 4))))
                except Exception as e:
                    continue
            log_f.close()

    if include_likelyhood_of_formula:
        for f_dict_str, s in formula_score_dict.items():
            try:
                f_dict = ast.literal_eval(f_dict_str)
                if ("C" in list(f_dict.keys())) and ("N" in list(f_dict.keys())):
                    if int(f_dict["N"]) > 2 and (int(f_dict["C"]) / int(f_dict["N"]) <= 4):
                        formula_score_dict[f_dict_str] = float(s) - ((float(f_dict["N"]) - 2) ** 2) * 50

                if ("C" in list(f_dict.keys())) and ("H" in list(f_dict.keys())):
                    if (float(f_dict["H"]) / float(f_dict["C"]) >= 2):
                        formula_score_dict[f_dict_str] = float(s) - (
                                ((float(f_dict["H"]) / float(f_dict["C"])) - 2) ** 3) * 50

                dbe = MS_functions.calc_dbe(f_dict)
                if dbe < 0:
                    score_substract = (abs(dbe + 2) * 50) ** 3
                elif (dbe - f_dict.get("O", 0)) > 0:
                    score_substract = abs(dbe - f_dict.get("O", 0) - 1) * 50
                else:
                    score_substract = 0
                formula_score_dict[f_dict_str] = float(s) - score_substract


            except Exception as e:
                print("ERROR in formula score likelyhood: " + str(e))

    for formula in list(formula_score_dict.items()):
        if formula[1] == 0:
            del formula_score_dict[formula[0]]
    formula_score_dict = dict(sorted(formula_score_dict.items(), key=lambda item: item[1], reverse=True))
    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("\n")
        log_f.write("===============================================================\n")
        log_f.write("===============================================================\n")
        log_f.write("Formula_score_dict\n")
        for item in list(formula_score_dict.items()):
            log_f.write("{:>30} \t {:>30} \n".format(str(item[0]), str(item[1])))
        log_f.close()
    formula_score_dict = {f: s for f, s in formula_score_dict.items() if
                          s >= reject_predicted_formula_if_score_lower_than}
    return formula_score_dict


def process_aif_spectrum(aif_spectrum, mass, save_folder="",
                         charge_of_measured_mass=-1,
                         remove_fragments_with_score_lower_than=-20,
                         max_nl_deviation_ppm=7.2,
                         divisor_for_fragment_noise=3,
                         max_ppm_deviation=30,
                         ppm_deviation_scoring_multiplier=3,
                         noise_divisor_for_isotopologue_calculation=2,
                         multiplier_isotopologue_influence_of_ppm_deviation_on_score=1,
                         isotopologue_score_e_function_exponent=0.2,
                         full_scan_spec=None,
                         nl_formula_cache_dict=None,
                         nl_formula_cache_dict_max_mass=230,
                         fragment_formula_cache_dict=None,
                         fragment_formula_cache_dict_max_mass=200,
                         reject_predicted_formula_if_score_lower_than=10,
                         minimum_assumed_noise=1000,
                         add_value_to_score_if_isotopo_was_found=100,
                         max_ppm_deviation_find_formula=25,
                         stop_isotopologue_search_if_score_lower_than=-200,
                         include_likelyhood_of_formula=True,
                         whole_file=None,
                         peak_rt=None,
                         peak_identification_width=3,
                         current_folder=None):
    # get new (measured) mass
    mass_old = copy.deepcopy(mass)
    print("Old mass: " + str(mass_old))
    if not full_scan_spec is None:
        dict_for_mass_determination = MS_functions.summarize_mass_intensity_dict_with_deviation(
            dict(zip(full_scan_spec[0], full_scan_spec[1])), deviation=max_ppm_deviation)
        mass = min(list(dict_for_mass_determination.keys()), key=lambda x: abs(mass - x))
    else:
        dict_for_mass_determination = MS_functions.summarize_mass_intensity_dict_with_deviation(
            dict(zip(aif_spectrum[0], aif_spectrum[1])), deviation=max_ppm_deviation)
        mass = min(list(dict_for_mass_determination.keys()), key=lambda x: abs(mass - x))
    if abs(((mass_old - mass) / mass) * 1000000) >= max_ppm_deviation_find_formula:
        print(
            "The deviation of the given mass and the closest mass exceeds the given max_ppm_deviation_find_formula of " + str(
                max_ppm_deviation_find_formula) + ". \nReturning False!!!")
        return False
    print("The given mass was adjusted from: " + str(mass_old) + " m/z to the new mass: " + str(mass) + " m/z. ")
    mass = mass + (charge_of_measured_mass * 0.000548)
    print("The corrected mass is: " + str(
        mass) + " m/z. The corrected mass is corrected for the electron. Continuing with the corrected mass!")
    if not full_scan_spec is None:
        intensity_of_requested_mass = sum([full_scan_spec[1][i] for i in range(len(full_scan_spec[0])) if
                                           abs(((full_scan_spec[0][i] - mass) / mass) * 1000000) <= max_ppm_deviation])
    else:
        intensity_of_requested_mass = sum([aif_spectrum[1][i] for i in range(len(aif_spectrum[0])) if
                                           abs(((aif_spectrum[0][i] - mass) / mass) * 1000000) <= max_ppm_deviation])
    print("Intensity of the requested mass is: " + str(intensity_of_requested_mass))

    log_f_filepath = ""
    if not save_folder == "":
        log_f_filepath = save_folder + "process_aifSpectrum_forMass_" + str(mass) + ".txt"
        log_f = open(log_f_filepath, "w")
        log_f.write(
            "The given mass was adjusted from: " + str(mass_old) + " m/z to the new mass: " + str(mass) + " m/z\n")
        log_f.write("Intensity of the requested mass is: " + str(intensity_of_requested_mass) + "\n")
        log_f.close()

    if full_scan_spec is None:
        formula_predict_spec = aif_spectrum
    else:
        formula_predict_spec = full_scan_spec
    prediction = predict_formula_of_mass_with_mass_spectrum(formula_predict_spec,
                                                            mass,
                                                            max_ppm_deviation=max_ppm_deviation,
                                                            ppm_deviation_scoring_multiplier=ppm_deviation_scoring_multiplier,
                                                            charge_of_measured_mass=charge_of_measured_mass,
                                                            save_folder=save_folder,
                                                            noise_divisor_for_isotopologue_calculation=noise_divisor_for_isotopologue_calculation,
                                                            multiplier_isotopologue_influence_of_ppm_deviation_on_score=multiplier_isotopologue_influence_of_ppm_deviation_on_score,
                                                            isotopologue_score_e_function_exponent=isotopologue_score_e_function_exponent,
                                                            reject_predicted_formula_if_score_lower_than=reject_predicted_formula_if_score_lower_than,
                                                            minimum_assumed_noise=minimum_assumed_noise,
                                                            add_value_to_score_if_isotopo_was_found=add_value_to_score_if_isotopo_was_found,
                                                            stop_isotopologue_search_if_score_lower_than=stop_isotopologue_search_if_score_lower_than,
                                                            fragment_formula_cache_dict=fragment_formula_cache_dict,
                                                            fragment_formula_cache_dict_max_mass=fragment_formula_cache_dict_max_mass,
                                                            include_likelyhood_of_formula=include_likelyhood_of_formula,
                                                            summarized_mass_intensity_dict=dict_for_mass_determination)

    # prediction = {{formula1}: score1, {formula2}: score2.......}
    for f_dict, s in list(prediction.items()):
        try:
            f_dict_dict = ast.literal_eval(f_dict)
            if ("C" in list(f_dict_dict.keys())) and ("N" in list(f_dict_dict.keys())):
                if (float(f_dict_dict["C"]) / 2 <= float(f_dict_dict["N"])) and int(f_dict_dict["N"]) >= 4:
                    del prediction[f_dict]
                    continue
        except Exception as e:
            print("Error in formula exclusion2: " + str(e))

    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("============================================================\n")
        log_f.write("PREDICTIONS FOR MOLECULE FORMULA (mass=" + str(round(mass, 4)) + ")\n")
        for item in list(prediction.items()):
            log_f.write(str(item[0]) + "\t" + str(item[1]) + "\n")
        log_f.write("============================================================\n \n")
        log_f.close()

    print({m: a for m, a in dict(zip(full_scan_spec[0], full_scan_spec[1])).items() if abs(m - mass) < 0.001})

    # create a dict with ordered intensities (high intensity to low intensity) {mass: intensity, ....} --> ordered_mass_intensities_dict
    # and remove measured masses that are below the noise.
    ordered_mass_intensities_dict = {}
    for entry in range(len(aif_spectrum[0])):
        ordered_mass_intensities_dict[aif_spectrum[0][entry]] = aif_spectrum[1][entry]
    ordered_mass_intensities_dict = dict(
        sorted(ordered_mass_intensities_dict.items(), key=lambda item: item[1], reverse=True))
    ordered_mass_intensities_dict = {k: v for k, v in ordered_mass_intensities_dict.items() if k <= mass}
    noise = intensity_of_requested_mass / divisor_for_fragment_noise
    ordered_mass_intensities_dict = {k: v for k, v in ordered_mass_intensities_dict.items() if v >= noise}
    ordered_mass_intensities_dict = MS_functions.summarize_mass_intensity_dict_with_deviation(
        ordered_mass_intensities_dict, max_ppm_deviation)
    print(ordered_mass_intensities_dict)
    print(len(ordered_mass_intensities_dict))
    if not full_scan_spec is None:
        full_scan_masses = full_scan_spec[0]
        full_scan_intensities = full_scan_spec[1]
        for element in list(ordered_mass_intensities_dict.items()):
            deviation_list = [abs(((element[0] - m) / m) * 1000000) for m in full_scan_masses]
            intensities_within_deviation_list = [full_scan_intensities[m] for m in range(len(full_scan_masses)) if
                                                 deviation_list[m] < max_ppm_deviation]
            total_intensity = sum(intensities_within_deviation_list)
            ordered_mass_intensities_dict[element[0]] = ordered_mass_intensities_dict[element[0]] - total_intensity
        ordered_mass_intensities_dict = {m: (i) for m, i in ordered_mass_intensities_dict.items() if i >= noise}

    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("============================================================\n")
        log_f.write(
            "Create a dictionary with all masses that are present in the current spectrum and their intensity {mass: intensity}...\n")
        log_f.write("Order the dictionary (highest to lowest intensity)...\n")
        log_f.write(
            "Remove all masses that are lower than the mass of the desired molecule (only possible fragments remain in the spectrum)...")
        log_f.write(
            "Delete all intensities that are below the noise-level (noise = intensity_of_requested_mass / divisor_for_fragment_noise)\n")
        log_f.write("         " + str(round(intensity_of_requested_mass, 2)) + " / " + str(
            divisor_for_fragment_noise) + " = " + str(round(noise, 2)) + "\n")
        log_f.write("{:>20} \t {:>20}\n".format("measured_mass", "intensity_of_mass"))
        for item in list(ordered_mass_intensities_dict.items()):
            log_f.write("{:>20} \t {:>20}\n".format(str(round(item[0], 4)), str(round(item[1], 1))))
        log_f.write("============================================================\n \n")
        log_f.close()

    # create a fragment peak list with every possible fragment (all fragments mass must be < mass of the molecule (mass))
    fragment_peak_list = []
    for element in ordered_mass_intensities_dict.keys():
        mass_of_fragment = element + (charge_of_measured_mass * 0.000548)
        intensity = ordered_mass_intensities_dict[element]
        fragment_peak_list.append([(mass - mass_of_fragment), mass_of_fragment, intensity])

    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("============================================================\n")
        log_f.write(
            "Create a fragment_peak_list where all possible fragments are listed according to [[neutral_loss, fragment_mz, intensity], ....]\n")
        log_f.write("{:>20} \t {:>20} \t {:>20}\n".format("neutral_loss_mass", "fragment_mass", "intensity"))
        for entry in fragment_peak_list:
            log_f.write("{:>20} \t {:>20} \t {:>20}\n".format(str(round(entry[0], 4)), str(round(entry[1], 4)),
                                                              str(round(entry[2], 1))))
        log_f.write("============================================================\n \n")
        log_f.close()
    # make fromula prediction with help of spectrum for all fragments.
    molecule_xic = MS_functions.get_xic(whole_file, mass, ((max_ppm_deviation * mass) / 1000000),
                                        requested_filter_mode="Full scan")
    rt_index = molecule_xic[0].index(min(molecule_xic[0], key=lambda x: abs(peak_rt - x)))
    molecule_normalized_intensity_list_surrounding = [molecule_xic[1][i] for i in range(len(molecule_xic[1])) if
                                                      rt_index - peak_identification_width <= i <= rt_index + peak_identification_width]
    molecule_normalized_intensity_list_surrounding = [round(((i - min(
        molecule_normalized_intensity_list_surrounding)) / (max(molecule_normalized_intensity_list_surrounding) - min(
        molecule_normalized_intensity_list_surrounding))), 1) for i in molecule_normalized_intensity_list_surrounding]

    dict_for_mass_determination = MS_functions.summarize_mass_intensity_dict_with_deviation(
        dict(zip(aif_spectrum[0], aif_spectrum[1])), deviation=max_ppm_deviation)

    for fragment in range(len(fragment_peak_list)):
        try:
            print("in for loop predict formula of fragments")
            predicted_formulas = predict_formula_of_mass_with_mass_spectrum(aif_spectrum,
                                                                            fragment_peak_list[fragment][1],
                                                                            max_ppm_deviation=max_ppm_deviation,
                                                                            ppm_deviation_scoring_multiplier=ppm_deviation_scoring_multiplier,
                                                                            charge_of_measured_mass=charge_of_measured_mass,
                                                                            noise_divisor_for_isotopologue_calculation=noise_divisor_for_isotopologue_calculation,
                                                                            multiplier_isotopologue_influence_of_ppm_deviation_on_score=multiplier_isotopologue_influence_of_ppm_deviation_on_score,
                                                                            isotopologue_score_e_function_exponent=isotopologue_score_e_function_exponent,
                                                                            fragment_formula_cache_dict=fragment_formula_cache_dict,
                                                                            fragment_formula_cache_dict_max_mass=fragment_formula_cache_dict_max_mass,
                                                                            reject_predicted_formula_if_score_lower_than=reject_predicted_formula_if_score_lower_than,
                                                                            minimum_assumed_noise=minimum_assumed_noise,
                                                                            add_value_to_score_if_isotopo_was_found=add_value_to_score_if_isotopo_was_found,
                                                                            stop_isotopologue_search_if_score_lower_than=stop_isotopologue_search_if_score_lower_than,
                                                                            include_likelyhood_of_formula=include_likelyhood_of_formula,
                                                                            summarized_mass_intensity_dict=dict_for_mass_determination)
            # {"{'C': 1, 'H': 2, 'P': 1, 'Cl': 1}": -254.60289862716849}
            predicted_formulas = {f: s for f, s in list(predicted_formulas.items()) if
                                  s > remove_fragments_with_score_lower_than}

            predicted_formulas = {f: s for f, s in list(predicted_formulas.items())}

            fragment_peak_list[fragment].append(predicted_formulas)
        except Exception as e:
            fragment_peak_list[fragment].append({})
            print("ERROR IN FRAGMENT FORMULA PREDICTION!!!")
            print(str(e))
            print(traceback.format_exc())
    fragment_peak_list = [[nl_m, f_m, i, d] for nl_m, f_m, i, d in fragment_peak_list if len(d) >= 1]
    # predict formulas based on the neutral loss mass. Does not include spectrum, so no isotopic information are used.
    for fragment in range(len(fragment_peak_list)):
        try:
            formulas_nl = MS_functions.get_formula_from_cache(
                "U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//", fragment_peak_list[fragment][0],
                max_ppm_deviation)
            formulas_nl = formulas_nl[1]
            formulas_nl = [[f, dev] for f, dev in list(formulas_nl.items()) if (abs(dev) < max_nl_deviation_ppm) or ()]
            fragment_peak_list[fragment].append(formulas_nl)
        except Exception as e1:
            print("Error in process_aif_spectrum, Formula prediction of NL: ")
            print(str(e1))
            print(traceback.format_exc())
    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("============================================================\n")
        log_f.write("Create formula predictions of fragment_peak_list...\n")
        log_f.write("Delete all fragments where no formula approximation for the fragment could be determined...\n")
        log_f.write(
            "Predict formula approximations for the neutral loss of each fragment ion (in relation to the desired molecule)...\n")
        log_f.write("Deleting all neutral loss formula approximations where the ppm deviation is higher than " + str(
            max_nl_deviation_ppm) + " ppm...\n")
        header = "{:>50} | {:>50} | {:>50} | {:>50} | {:>50}\n".format("neutral_loss_mass",
                                                                       "neutral_loss_formulas with score",
                                                                       "fragment_mass", "fragment_formulas with score",
                                                                       "intensity")
        log_f.write(header)
        for entry in fragment_peak_list:
            line = "{:>50} | {:>50} | {:>50} | {:>50} | {:>50}\n".format(str(round(entry[0], 4)), str(entry[4]),
                                                                         str(round(entry[1], 4)), str(entry[3]),
                                                                         str(round(entry[2], 1)))
            log_f.write(line)
        log_f.write("============================================================\n \n")
        log_f.close()

    # calculate sum of formulas of neutral loss and fragment. See if this matches with one of the predicted formulas of the requested mass.
    molecule_xic = MS_functions.get_xic(whole_file, mass, ((max_ppm_deviation * mass) / 1000000),
                                        requested_filter_mode="Full scan")
    rt_index = molecule_xic[0].index(min(molecule_xic[0], key=lambda x: abs(peak_rt - x)))
    molecule_normalized_intensity_list_surrounding = [molecule_xic[1][i] for i in range(len(molecule_xic[1])) if
                                                      rt_index - peak_identification_width <= i <= rt_index + peak_identification_width]
    molecule_normalized_intensity_list_surrounding = [round(((i - min(
        molecule_normalized_intensity_list_surrounding)) / (max(molecule_normalized_intensity_list_surrounding) - min(
        molecule_normalized_intensity_list_surrounding))), 1) for i in molecule_normalized_intensity_list_surrounding]
    max_molecule_index = molecule_normalized_intensity_list_surrounding.index(max(molecule_normalized_intensity_list_surrounding))


    true_fragment_list = []
    formula_predictions_list = []
    for molecule_formula in list(prediction.items()):
        molecule_formula_f = molecule_formula[0]
        molecule_formula_f = ast.literal_eval(molecule_formula_f)
        molecule_formula_f = "".join([str(a) + str(n) for a, n in molecule_formula_f.items()])
        formula_predictions_list.append([molecule_formula_f, molecule_formula[1]])
    for fragment in range(len(fragment_peak_list)):
        fragment_mass = round(fragment_peak_list[fragment][1], 4)
        nl_mass = round(fragment_peak_list[fragment][0], 4)

        for nl in fragment_peak_list[fragment][4]:
            nl_deviation = nl[1]
            nl = nl[0]
            nl = MS_functions.get_formula_to_dict(nl)
            for fragment_formula in list(fragment_peak_list[fragment][3].items()):
                fragment_formula_score = fragment_formula[1]
                fragment_formula = ast.literal_eval(fragment_formula[0])
                molecule_formula_nl_f = {key: nl.get(key, 0) + fragment_formula.get(key, 0) for key in
                                         set(nl) | set(fragment_formula)}
                for molecule_formula in list(prediction.items()):
                    molecule_formula_f = molecule_formula[0]
                    molecule_formula_f = ast.literal_eval(molecule_formula_f)
                    if molecule_formula_nl_f == molecule_formula_f:
                        print("TRUETRUEdjlgajkölsdgjklöasdjklgösjdklöfgajksöldgjöl")

                        frag_xic = MS_functions.get_xic(whole_file, fragment_peak_list[fragment][1], (
                                (max_ppm_deviation * fragment_peak_list[fragment][1]) / 1000000),
                                                        requested_filter_mode="AIF")

                        fig, ax = plt.subplots()
                        ax.plot(frag_xic[0], frag_xic[1], label="XIC of fragment")
                        ax.set_title("XIC of Fragment with mass: " + str(round(fragment_mass, 4)) + ". Peak at RT: " + str(round(peak_rt, 2)))
                        ax.set_xlabel("retention time / min")
                        ax.set_ylabel("intensity / a.u.")
                        plt.legend()
                        matplotlib.rcParams.update({'figure.autolayout': True})
                        image_filepath = current_folder + "fragments//"
                        if not os.path.isdir(image_filepath):
                            os.makedirs(image_filepath)
                        image_filepath = image_filepath + "xic_" + str(round(fragment_mass, 4)) + ".png"
                        plt.savefig(image_filepath, bbox_inches='tight')
                        rt_index = frag_xic[0].index(min(frag_xic[0], key=lambda x: abs(peak_rt - x)))
                        noise_frag = [frag_xic[1][i] for i in range(len(frag_xic[1])) if (rt_index - (2*peak_identification_width) <= i <= rt_index-peak_identification_width) or (rt_index + (peak_identification_width) <= i <= rt_index+2*peak_identification_width)]
                        std_dev_frag = statistics.stdev(noise_frag)
                        frag_intensity_list_surrounding = [frag_xic[1][i] for i in range(len(frag_xic[1])) if
                                                           rt_index - peak_identification_width <= i <= rt_index + peak_identification_width]
                        std_dev_multiplier = (sum(frag_intensity_list_surrounding)/len(frag_intensity_list_surrounding)) / std_dev_frag
                        print(std_dev_multiplier)
                        frag_normalized_intensity_list_surrounding = [round(((i - min( frag_intensity_list_surrounding)) / ( max(frag_intensity_list_surrounding) - min(
                                                                                     frag_intensity_list_surrounding))),
                                                                            1) for i in
                                                                      frag_intensity_list_surrounding]
                        max_frag_index = frag_normalized_intensity_list_surrounding.index(max(frag_normalized_intensity_list_surrounding))
                        delta_max_intensities = abs(max_frag_index - max_molecule_index)
                        print(delta_max_intensities)
                        try:
                            index_delta_multiplier = (1/(delta_max_intensities/peak_identification_width)) - 1
                        except:
                            index_delta_multiplier = peak_identification_width
                        print(index_delta_multiplier)
                        peak_height_multiplier = (max(frag_intensity_list_surrounding) / (sum(frag_intensity_list_surrounding) / len(frag_intensity_list_surrounding))) - 1
                        if peak_height_multiplier < 0.75:
                            peak_height_multiplier = (-1/(peak_height_multiplier**2)) * 20
                        else:
                            peak_height_multiplier = (peak_height_multiplier+0.3)**4
                        print(peak_height_multiplier)
                        if (peak_identification_width >= 4) and (peak_height_multiplier > 5):
                            if index_delta_multiplier == 0:
                                mysimilarity = -200
                            elif index_delta_multiplier < 1:
                                mysimilarity = -(1/(index_delta_multiplier+0.0001))**3
                            else:
                                mysimilarity = index_delta_multiplier**2
                            mysimilarity = mysimilarity*10 + (peak_height_multiplier)
                        else:
                            mysimilarity = peak_height_multiplier*3
                        print(mysimilarity)

                        if mysimilarity < -100:
                            continue

                        add_to_score = mysimilarity
                        #print(fragment_formula_score)
                        #fragment_formula_score = fragment_formula_score + add_to_score
                        print(fragment_formula_score)

                        true_fragment_list.append(
                            [[molecule_formula_f, intensity_of_requested_mass, molecule_formula[1], mass],
                             [fragment_formula, fragment_peak_list[fragment][2], fragment_formula_score, fragment_mass],
                             [nl, nl_deviation, nl_mass]])

    # remove duplica entries with the same neutral loss and the same frag formula.
    nl_list = [true_fragment_list[i][2][0] for i in range(len(true_fragment_list))]
    frag_formula_list = [true_fragment_list[i][1][0] for i in range(len(true_fragment_list))]
    myset_of_unique_combination = []
    delete_index_list = []
    for entry in range(len(list(zip(frag_formula_list, nl_list)))):
        if not list(zip(frag_formula_list, nl_list))[entry] in myset_of_unique_combination:
            myset_of_unique_combination.append(list(zip(frag_formula_list, nl_list))[entry])
        else:
            delete_index_list.append(entry)
    true_fragment_list = [true_fragment_list[i] for i in range(len(true_fragment_list)) if not i in delete_index_list]

    # removing entries, where two different formula have the same neutral loss. Then only the best formula prediction is kept and the other is discarded.
    # shortened_true_fragment_list = []
    # nl_list = []
    # formula_list = []
    # combined_score_list = []
    # print(true_fragment_list)
    # for entry in range(len(true_fragment_list)):
    #     nl_list.append(true_fragment_list[entry][2][0])
    #     formula_list.append(true_fragment_list[entry][0][0])
    #     combined_score_list.append(true_fragment_list[entry][0][2] + true_fragment_list[entry][1][2])
    #
    # nl_list = [nl for _, nl in sorted(zip(combined_score_list, nl_list), reverse=False)]
    # nl_string_list = []
    # for entry in nl_list:
    #     nl_string_list.append("".join([str(a) + str(n) for a, n in entry.items()]))
    #
    # formula_list = [f for _, f in sorted(zip(combined_score_list, formula_list), reverse=False)]
    # formula_string_list = []
    # for entry in formula_list:
    #     formula_string_list.append("".join([str(a) + str(n) for a, n in entry.items()]))
    #
    # mydict = dict(zip(nl_string_list, formula_string_list))
    # mydict_keys_list = []
    # mydict_values_list = []
    # for item in list(mydict.items()):
    #     mydict_keys_list.append(MS_functions.get_formula_to_dict(item[0]))
    #     mydict_values_list.append(MS_functions.get_formula_to_dict(item[1]))
    #
    # for entry in range(len(true_fragment_list)):
    #     if (true_fragment_list[entry][2][0] in mydict_keys_list) and (
    #             true_fragment_list[entry][0][0] in mydict_values_list):
    #         shortened_true_fragment_list.append(true_fragment_list[entry])

    # true_fragment_list = copy.deepcopy(shortened_true_fragment_list)

    if not save_folder == "":
        log_f = open(log_f_filepath, "a")
        log_f.write("============================================================\n")
        log_f.write(
            "Checking if ANY neutral_loss_formula_prediction that was found fits to ANY fragment_formula_prediction to give ANY formula_prediction of the molecule...\n")
        header = "{:>50} | {:>50} | {:>50}\n".format("Molecule", "Fragment", "Neutral_loss")
        log_f.write(header)
        for entry in true_fragment_list:
            line = "{:>50} | {:>50} | {:>50}\n".format(str(entry[0]), str(entry[1]), str(entry[2]))
            log_f.write(line)
        log_f.write("============================================================\n \n")
        log_f.close()

    return [true_fragment_list, fragment_peak_list, prediction, mass, intensity_of_requested_mass,
            formula_predictions_list]


def find_masses_with_peak(f, range_to_search=(50, 300),
                          sg_windows_bg_removal=(10, 75, 2), width_bg_removal=(4, 20),
                          min_distance_of_chromatographic_peaks_in_seconds=10,
                          prominence_for_peak_detection=0.05, min_relative_height_for_peak_detection=0.1,
                          neighbour_window=20,
                          include_peak_if_height_multiple_of_std_dev=10):
    valid_masses = []
    stepsize = 0.01
    total_masses_range = [m * stepsize for m in
                          range(int(range_to_search[0] / stepsize), int(range_to_search[1] / stepsize))]
    total_one_decimal_range = [m for m in total_masses_range if ((m % 1) <= 0.35 or (m % 1) >= 0.7)]
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    xic = MS_functions.get_xic(f, 138.0191, mass_deviation=0.001, requested_filter_mode="Full scan")
    # plt.ion()
    # fig, ax = plt.subplots()
    # line1, = ax.plot(xic[0], xic[1])
    for mass in total_one_decimal_range:
        xic = MS_functions.get_xic(f, mass, mass_deviation=stepsize, requested_filter_mode="Full scan")
        approx_capture_duration_per_element = (max(rt_list) * 60) / len(xic[1])
        bg_subst_int = do_bg_substraction(xic[1], approx_capture_duration_per_element, sg_windows=sg_windows_bg_removal,
                                          width=width_bg_removal,
                                          distance=min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element)
        prominence_for_peak_detection_curr = prominence_for_peak_detection * (max(bg_subst_int[0]) + 1000)
        peaks = scipy.signal.find_peaks(bg_subst_int[0],
                                        height=max(bg_subst_int[0]) * min_relative_height_for_peak_detection + 10000,
                                        prominence=prominence_for_peak_detection_curr,
                                        distance=min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element,
                                        width=width_bg_removal)
        if len(peaks[0]) >= 1:
            # dynamic peak finding.
            # get std dev of xic in direct neighbourhood of the peak (excluding the peak itself)
            std_dev_neighbourhood = 0
            intensities_in_neighbourhood = []
            new_peak_list = []
            for peak_index in peaks[0]:
                try:
                    intensities_in_neighbourhood_1 = bg_subst_int[0][
                                                     peak_index - neighbour_window:int(
                                                         peak_index - neighbour_window / 2)]
                    intensities_in_neighbourhood_2 = bg_subst_int[0][
                                                     int(peak_index + neighbour_window / 2):peak_index + neighbour_window]
                    for element in intensities_in_neighbourhood_1:
                        intensities_in_neighbourhood.append(element)
                    for element in intensities_in_neighbourhood_2:
                        intensities_in_neighbourhood.append(element)
                    mean = sum(intensities_in_neighbourhood) / len(intensities_in_neighbourhood)
                    variance = sum([((x - mean) ** 2) for x in intensities_in_neighbourhood]) / len(
                        intensities_in_neighbourhood)
                    res = variance ** 0.5
                except:
                    res = 0
                if (bg_subst_int[0][peak_index] >= (include_peak_if_height_multiple_of_std_dev * res)):
                    new_peak_list.append(peak_index)
            peaks = list(peaks)
            keep_list_indices = [l_i for l_i in list(range(len(peaks[0]))) if peaks[0][l_i] in new_peak_list]
            peaks[0] = copy.deepcopy(new_peak_list)
            for key, value in peaks[1].items():
                peaks[1][key] = [value[v] for v in range(len(value)) if v in keep_list_indices]
        if not len(peaks[0]) == 0 and len(peaks[0]) <= 4:
            valid_masses.append(mass)
            # line1.set_ydata(xic[1])
            # ax.relim()
            # ax.autoscale_view()
            # fig.canvas.draw()
            # fig.canvas.flush_events()

    new_stepsize = 0.001
    new_masslist = []
    for mass in valid_masses:
        values_added_to_individual_mass = [offset * new_stepsize for offset in
                                           range(-round(stepsize / new_stepsize), round(stepsize / new_stepsize), 1)]
        new_masses_for_mass = [mass + offset for offset in values_added_to_individual_mass]
        for element in new_masses_for_mass:
            new_masslist.append(element)

    valid_masses = []
    stepsize = 0.001
    for mass in new_masslist:
        xic = MS_functions.get_xic(f, mass, mass_deviation=new_stepsize, requested_filter_mode="Full scan")
        approx_capture_duration_per_element = (max(rt_list) * 60) / len(xic[1])
        bg_subst_int = do_bg_substraction(xic[1], approx_capture_duration_per_element, sg_windows=sg_windows_bg_removal,
                                          width=width_bg_removal,
                                          distance=min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element)
        prominence_for_peak_detection_curr = prominence_for_peak_detection * (max(bg_subst_int[0]) + 50000)
        peaks = scipy.signal.find_peaks(bg_subst_int[0],
                                        height=max(bg_subst_int[0]) * min_relative_height_for_peak_detection + 30000,
                                        prominence=prominence_for_peak_detection_curr,
                                        distance=min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element,
                                        width=width_bg_removal)
        if len(peaks[0]) >= 1:
            # dynamic peak finding.
            # get std dev of xic in direct neighbourhood of the peak (excluding the peak itself)
            std_dev_neighbourhood = 0
            intensities_in_neighbourhood = []
            new_peak_list = []
            for peak_index in peaks[0]:
                try:
                    intensities_in_neighbourhood_1 = bg_subst_int[0][
                                                     peak_index - neighbour_window:int(
                                                         peak_index - neighbour_window / 2)]
                    intensities_in_neighbourhood_2 = bg_subst_int[0][
                                                     int(peak_index + neighbour_window / 2):peak_index + neighbour_window]
                    for element in intensities_in_neighbourhood_1:
                        intensities_in_neighbourhood.append(element)
                    for element in intensities_in_neighbourhood_2:
                        intensities_in_neighbourhood.append(element)
                    mean = sum(intensities_in_neighbourhood) / len(intensities_in_neighbourhood)
                    variance = sum([((x - mean) ** 2) for x in intensities_in_neighbourhood]) / len(
                        intensities_in_neighbourhood)
                    res = variance ** 0.5
                except:
                    res = 0
                if (bg_subst_int[0][peak_index] >= (include_peak_if_height_multiple_of_std_dev * res)):
                    new_peak_list.append(peak_index)
            peaks = list(peaks)
            keep_list_indices = [l_i for l_i in list(range(len(peaks[0]))) if peaks[0][l_i] in new_peak_list]
            peaks[0] = copy.deepcopy(new_peak_list)
            for key, value in peaks[1].items():
                peaks[1][key] = [value[v] for v in range(len(value)) if v in keep_list_indices]
        if not len(peaks[0]) == 0 and len(peaks[0]) <= 4:
            valid_masses.append(mass)
            # line1.set_ydata(xic[1])
            # ax.relim()
            # ax.autoscale_view()
            # fig.canvas.draw()
            # fig.canvas.flush_events()

    new_stepsize = 0.0003
    new_masslist = []
    for mass in valid_masses:
        values_added_to_individual_mass = [offset * new_stepsize for offset in
                                           range(-round(stepsize / new_stepsize), round(stepsize / new_stepsize), 1)]
        new_masses_for_mass = [mass + offset for offset in values_added_to_individual_mass]
        for element in new_masses_for_mass:
            new_masslist.append(element)
    stepsize = 0.001
    valid_masses = []
    for mass in new_masslist:
        xic = MS_functions.get_xic(f, mass, mass_deviation=new_stepsize, requested_filter_mode="Full scan")
        approx_capture_duration_per_element = (max(rt_list) * 60) / len(xic[1])
        bg_subst_int = do_bg_substraction(xic[1], approx_capture_duration_per_element, sg_windows=sg_windows_bg_removal,
                                          width=width_bg_removal,
                                          distance=min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element)
        prominence_for_peak_detection_curr = prominence_for_peak_detection * (max(bg_subst_int[0]) + 50000)
        peaks = scipy.signal.find_peaks(bg_subst_int[0],
                                        height=max(bg_subst_int[0]) * min_relative_height_for_peak_detection + 30000,
                                        prominence=prominence_for_peak_detection_curr,
                                        distance=min_distance_of_chromatographic_peaks_in_seconds / approx_capture_duration_per_element,
                                        width=width_bg_removal)
        if len(peaks[0]) >= 1:
            # dynamic peak finding.
            # get std dev of xic in direct neighbourhood of the peak (excluding the peak itself)
            std_dev_neighbourhood = 0
            intensities_in_neighbourhood = []
            new_peak_list = []
            for peak_index in peaks[0]:
                try:
                    intensities_in_neighbourhood_1 = bg_subst_int[0][
                                                     peak_index - neighbour_window:int(
                                                         peak_index - neighbour_window / 2)]
                    intensities_in_neighbourhood_2 = bg_subst_int[0][
                                                     int(peak_index + neighbour_window / 2):peak_index + neighbour_window]
                    for element in intensities_in_neighbourhood_1:
                        intensities_in_neighbourhood.append(element)
                    for element in intensities_in_neighbourhood_2:
                        intensities_in_neighbourhood.append(element)
                    mean = sum(intensities_in_neighbourhood) / len(intensities_in_neighbourhood)
                    variance = sum([((x - mean) ** 2) for x in intensities_in_neighbourhood]) / len(
                        intensities_in_neighbourhood)
                    res = variance ** 0.5
                except:
                    res = 0
                if (bg_subst_int[0][peak_index] >= (include_peak_if_height_multiple_of_std_dev * res)):
                    new_peak_list.append(peak_index)
            peaks = list(peaks)
            keep_list_indices = [l_i for l_i in list(range(len(peaks[0]))) if peaks[0][l_i] in new_peak_list]
            peaks[0] = copy.deepcopy(new_peak_list)
            for key, value in peaks[1].items():
                peaks[1][key] = [value[v] for v in range(len(value)) if v in keep_list_indices]
        if not len(peaks[0]) == 0 and len(peaks[0]) <= 4:
            valid_masses.append(mass)
            print("FOUND MASS: " + str(mass))
            # line1.set_ydata(xic[1])
            # ax.relim()
            # ax.autoscale_view()
            # fig.canvas.draw()
            # fig.canvas.flush_events()
    valid_masses = [(round(m * 5, 3) / 5) for m in valid_masses]
    valid_masses = list(set(valid_masses))
    return valid_masses


def one_analysis(mass, rt, f, mass_deviation=0.005,
                 requested_xic_mode="Full scan",
                 remove_if_combined_fragment_score_lower_than=-10, charge_of_measured_mass=-1,
                 remove_fragments_with_score_lower_than=-10, max_nl_deviation_ppm=7.2,
                 divisor_for_fragment_noise=3,
                 max_ppm_deviation=20, ppm_deviation_scoring_multiplier=3,
                 noise_divisor_for_isotopologue_calculation=2,
                 multiplier_isotopologue_influence_of_ppm_deviation_on_score=1,
                 isotopologue_score_e_function_exponent=0.2,
                 summary_filepath="",
                 additional_information_about_peak=[],
                 nl_formula_cache_dict=None,
                 nl_formula_cache_dict_max_mass=230,
                 fragment_formula_cache_dict=None,
                 fragment_formula_cache_dict_max_mass=200,
                 reject_predicted_formula_if_score_lower_than=10,
                 minimum_assumed_noise=1000,
                 add_value_to_score_if_isotopo_was_found=100,
                 max_ppm_deviation_find_formula=25,
                 stop_isotopologue_search_if_score_lower_than=-200,
                 include_likelyhood_of_formula=True,
                 peak_identification_width=8):
    def set_size(w, h, ax=None):
        """ w, h: width, height in inches """
        if not ax: ax = plt.gca()
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        figw = float(w) / (r - l)
        figh = float(h) / (t - b)
        ax.figure.set_size_inches(figw, figh)

    current_folder = output_folder + str(round(mass, 4)) + "_" + str(rt) + "//"
    if not os.path.isdir(current_folder):
        os.makedirs(current_folder)

    xic = MS_functions.get_xic(f, mass, mass_deviation, requested_filter_mode=requested_xic_mode)
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    peak_index = rt_list.index(min(rt_list, key=lambda x: abs(rt - x)))
    xic_peak_index = xic[0].index(min(xic[0], key=lambda x: abs(rt - x)))

    fig, ax = plt.subplots()
    ax.scatter(rt_list[peak_index], xic[1][xic_peak_index], color="red", marker="x", s=100)
    ax.plot(xic[0], xic[1], label="XIC measured")
    ax.set_title("xic_" + str(round(mass, 4)) + "+-" + str(mass_deviation) + "_" + str(requested_xic_mode))
    ax.set_xlabel("retention time / min")
    ax.set_ylabel("intensity / a.u.")
    plt.legend()
    matplotlib.rcParams.update({'figure.autolayout': True})
    image_filepath = current_folder + "xic_" + str(round(mass, 4)) + "+-" + str(mass_deviation) + "_" + str(
        requested_xic_mode) + ".png"
    plt.savefig(image_filepath, bbox_inches='tight')
    metadata = PIL.PngImagePlugin.PngInfo()
    metadata.add_text("xic_retention_time", str(xic[0]))
    metadata.add_text("xic_intensity_list", str(xic[1]))
    metadata.add_text("xic_original_index_list", str(xic[2]))
    target_image = PIL.Image.open(image_filepath)
    target_image.save(image_filepath, pnginfo=metadata)

    try:
        full_scan_spec = MS_functions.get_mass_spectrum(f, peak_index, mode="as_index",
                                                        requested_filter_mode="Full scan",
                                                        save_folder=current_folder)
        print({m: a for m, a in dict(zip(full_scan_spec[0], full_scan_spec[1])).items() if abs(m - mass) < 0.001})
        MS_functions.get_mass_spectrum(f, peak_index, mode="as_index", requested_filter_mode="MS/MS",
                                       requested_ms_ms_mass=mass, save_folder=current_folder)

        spec = MS_functions.get_mass_spectrum(f, peak_index, mode="as_index", requested_filter_mode="AIF",
                                              save_folder=current_folder)
        aif_processed = process_aif_spectrum(spec,
                                             mass,
                                             save_folder=current_folder,
                                             charge_of_measured_mass=charge_of_measured_mass,
                                             remove_fragments_with_score_lower_than=remove_fragments_with_score_lower_than,
                                             max_nl_deviation_ppm=max_nl_deviation_ppm,
                                             divisor_for_fragment_noise=divisor_for_fragment_noise,
                                             max_ppm_deviation=max_ppm_deviation,
                                             ppm_deviation_scoring_multiplier=ppm_deviation_scoring_multiplier,
                                             noise_divisor_for_isotopologue_calculation=noise_divisor_for_isotopologue_calculation,
                                             multiplier_isotopologue_influence_of_ppm_deviation_on_score=multiplier_isotopologue_influence_of_ppm_deviation_on_score,
                                             isotopologue_score_e_function_exponent=isotopologue_score_e_function_exponent,
                                             full_scan_spec=full_scan_spec,
                                             nl_formula_cache_dict=nl_formula_cache_dict,
                                             nl_formula_cache_dict_max_mass=nl_formula_cache_dict_max_mass,
                                             fragment_formula_cache_dict=fragment_formula_cache_dict,
                                             fragment_formula_cache_dict_max_mass=fragment_formula_cache_dict_max_mass,
                                             reject_predicted_formula_if_score_lower_than=reject_predicted_formula_if_score_lower_than,
                                             minimum_assumed_noise=minimum_assumed_noise,
                                             add_value_to_score_if_isotopo_was_found=add_value_to_score_if_isotopo_was_found,
                                             max_ppm_deviation_find_formula=max_ppm_deviation_find_formula,
                                             stop_isotopologue_search_if_score_lower_than=stop_isotopologue_search_if_score_lower_than,
                                             include_likelyhood_of_formula=include_likelyhood_of_formula,
                                             whole_file=f,
                                             peak_rt=rt,
                                             peak_identification_width=peak_identification_width,
                                             current_folder=current_folder)
        if aif_processed == None:
            print("Peak is not further processed, as process_aif_spectrum returned None!")
            return
    except Exception as e_process_aif:
        print("ERROR IN PROCESSING AIF SPECTRUM! Error code:")
        print(str(e_process_aif))
        print(traceback.format_exc())
        return
    try:
        true_fragment_list = aif_processed[0]
        molecule_predictions_list = aif_processed[5]
    except Exception as e:
        print("Error while processing AIF spectrum! Error code: ")
        print(str(e))
        print(traceback.format_exc())
        return

    fig, ax = plt.subplots()
    ax.bar(spec[0], spec[1], color="gray", label="Spectrum", alpha=0.2)
    ax.bar(aif_processed[3], aif_processed[4], color="red", width=1.5, label="Molecule")
    delete_index_list = []
    combined_score_list = []
    for fragment in range(len(true_fragment_list)):
        combined_score = true_fragment_list[fragment][0][2] + true_fragment_list[fragment][1][2] - abs(
            true_fragment_list[fragment][2][1])
        if combined_score < remove_if_combined_fragment_score_lower_than:
            delete_index_list.append(fragment)
        else:
            combined_score_list.append(combined_score)
    for i in sorted(delete_index_list, reverse=True):
        del true_fragment_list[i]
    sorted_index_list = [combined_score_list.index(i) for i in sorted(combined_score_list, reverse=True)]
    true_fragment_list = [true_fragment_list[i] for i in sorted_index_list]
    fragment_mass_list = []
    fragment_intensity_list = []
    nl_mass_list = []
    nl_formula_list = []
    fragment_formula_list = []
    molecule_formula_list = []
    scoring_list = []
    fragment_score_list = []
    nl_deviation_list = []
    molecule_score_list = []
    text_str = ""
    text_str = text_str + "{:<10}|{:<10}|{:<10}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format("score",
                                                                                                      "Frag mass",
                                                                                                      "intensity",
                                                                                                      "Frag formula",
                                                                                                      "NL formula",
                                                                                                      "Molec formula",
                                                                                                      "score F",
                                                                                                      "dev. NL",
                                                                                                      "score M")
    text_str = text_str + "{:_<10}|{:_<10}|{:_<10}|{:_<15}_{:_<15}___{:_<15}|{:_>10}|{:_>10}|{:_>10}\n".format("",
                                                                                                               "",
                                                                                                               "",
                                                                                                               "",
                                                                                                               "",
                                                                                                               "",
                                                                                                               "",
                                                                                                               "",
                                                                                                               "")
    for molecule in molecule_predictions_list:
        text_str = text_str + \
                   "{:<10}|{:<10}|{:<10}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format(
                       str(round(float(molecule[1]), 2)),
                       0,
                       "100%",
                       "",
                       "",
                       molecule[0],
                       0,
                       0,
                       str(round(float(molecule[1]), 2)))
    print(text_str)
    for fragment in true_fragment_list:
        combined_score = fragment[0][2] + fragment[1][2] - abs(fragment[2][1])
        text_str = text_str + \
                   "{:<10}|{:<10}|{:<10}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format(
                       round(combined_score, 2), str(round(fragment[1][3], 4)),
                       str(round(((fragment[1][1] / aif_processed[4]) * 100), 2)) + "%",
                       str("".join([str(a) + str(n) for a, n in fragment[1][0].items()])),
                       str("".join([str(a) + str(n) for a, n in fragment[2][0].items()])),
                       str("".join([str(a) + str(n) for a, n in fragment[0][0].items()])), round(fragment[1][2], 2),
                       round(abs(fragment[2][1]), 2), round(fragment[0][2], 2))
        fragment_mass_list.append(fragment[1][3])
        fragment_intensity_list.append(fragment[1][1])
        nl_mass_list.append(fragment[2][2])
        nl_formula_list.append(fragment[2][0])
        fragment_formula_list.append(fragment[1][0])
        molecule_formula_list.append(fragment[0][0])
        fragment_score_list.append(fragment[1][2])
        nl_deviation_list.append(fragment[2][1])
        molecule_score_list.append(fragment[0][2])
        combined_score = fragment[0][2] + fragment[1][2] - abs(fragment[2][1])
        scoring_list.append(combined_score)
    molecule_formula_dict = {}
    for formula in range(len(molecule_formula_list)):
        formula_string = str("".join([str(a) + str(n) for a, n in molecule_formula_list[formula].items()]))
        if formula_string not in list(molecule_formula_dict.keys()):
            molecule_formula_dict[formula_string] = [[fragment_mass_list[formula]],
                                                     [fragment_intensity_list[formula]]]
        else:
            molecule_formula_dict[formula_string][0].append(fragment_mass_list[formula])
            molecule_formula_dict[formula_string][1].append(fragment_intensity_list[formula])
    for molecule_formula_item in list(molecule_formula_dict.items()):
        ax.bar(molecule_formula_item[1][0], molecule_formula_item[1][1], width=1.4,
               label="Fragments for: " + str(molecule_formula_item[0]))
    text_str = text_str.strip()
    props = dict(boxstyle='round', facecolor='grey', alpha=0.05)  # bbox features
    text = ax.text(1.03, 0.98, text_str, fontfamily="monospace", transform=ax.transAxes, fontsize=10,
                   verticalalignment="top", bbox=props)
    text_width = text.get_window_extent(renderer=fig.canvas.get_renderer()).width
    text_height = text.get_window_extent(renderer=fig.canvas.get_renderer()).height
    set_size(4 + (text_width / fig.dpi), 4 + (text_height / fig.dpi), ax)

    filter_mode = spec[6]
    index = spec[2]
    filter = spec[5]
    ms_ms_masses = spec[7]
    ax.legend(loc="upper left")
    ax.set_xlabel("masses / Da")
    ax.set_ylabel("intensity / a.u.")
    ax.set_title("Mode:" + str(filter_mode) + "; Index: " + str(index) + "; RT: " + str(
        round(rt_list[index], 2)) + " min" + ";\nFilter: " + str(filter) + ";\nMS/MS masses: " + str(ms_ms_masses))
    matplotlib.rcParams.update({'figure.autolayout': True})
    image_filepath = current_folder + "Mass_spectrum_FRAGMENTS_index" + str(index) + "_" + str(
        filter_mode) + "_" + str(ms_ms_masses) + ".png"
    plt.savefig(image_filepath, bbox_inches='tight')
    metadata = PIL.PngImagePlugin.PngInfo()
    metadata.add_text("masses", str(spec[0]))
    metadata.add_text("fragment_masses", str(fragment_mass_list))
    metadata.add_text("intensities", str(spec[1]))
    metadata.add_text("fragment_intensities", str(fragment_intensity_list))
    metadata.add_text("molecule_mass", str(aif_processed[3]))
    metadata.add_text("molecule_intensities", str(aif_processed[4]))
    metadata.add_text("index", str(index))
    metadata.add_text("orig_file_entry_for_index", str(f[index]))
    metadata.add_text("mode", str(filter_mode))
    metadata.add_text("filter", str(filter))
    metadata.add_text("msms_masses", str(ms_ms_masses))
    target_image = PIL.Image.open(image_filepath)
    target_image.save(image_filepath, pnginfo=metadata)

    total_combined_score_per_formula_prediction = {}
    for entry in range(len(molecule_formula_list)):
        formula_string = "".join([str(a) + str(n) for a, n in molecule_formula_list[entry].items()])
        try:
            total_combined_score_per_formula_prediction[formula_string] = total_combined_score_per_formula_prediction[
                                                                              formula_string] + scoring_list[entry]
        except:
            total_combined_score_per_formula_prediction[formula_string] = scoring_list[entry]
    fragments_for_best_formula = []
    score_of_fragments = []
    nl_for_best_formula = []
    score_of_nl = []
    frag_intensity_for_best_formula = []
    score_of_molecule_list = []

    for entry in range(len(molecule_formula_list)):
        formula_string = "".join([str(a) + str(n) for a, n in molecule_formula_list[entry].items()])
        fragment_string = "".join([str(a) + str(n) for a, n in fragment_formula_list[entry].items()])
        nl_string = "".join([str(a) + str(n) for a, n in nl_formula_list[entry].items()])
        if formula_string == max(total_combined_score_per_formula_prediction,
                                 key=total_combined_score_per_formula_prediction.get):
            fragments_for_best_formula.append(fragment_string)
            nl_for_best_formula.append(nl_string)
            score_of_fragments.append(fragment_score_list[entry])
            score_of_nl.append(nl_deviation_list[entry])
            frag_intensity_for_best_formula.append(fragment_intensity_list[entry])
            score_of_molecule_list.append(molecule_score_list[entry])

    best_formula = []
    try:
        if len(molecule_formula_list) == 0:
            for molecule in molecule_predictions_list:
                best_formula.append(molecule[0])
                best_formula.append(molecule[1])
                best_formula.append(mass)
                best_formula.append(index)
                best_formula.append(aif_processed[4])  # intensity of molecule
                best_formula.append(fragments_for_best_formula)
                best_formula.append(frag_intensity_for_best_formula)
                best_formula.append(nl_for_best_formula)
                best_formula.append(score_of_fragments)
                best_formula.append(score_of_nl)
                best_formula.append(score_of_molecule_list)

        else:
            best_formula.append(
                max(total_combined_score_per_formula_prediction, key=total_combined_score_per_formula_prediction.get))
            best_formula.append(total_combined_score_per_formula_prediction[
                                    max(total_combined_score_per_formula_prediction,
                                        key=total_combined_score_per_formula_prediction.get)])
            best_formula.append(mass)
            best_formula.append(index)
            best_formula.append(aif_processed[4])  # intensity of molecule
            best_formula.append(fragments_for_best_formula)
            best_formula.append(frag_intensity_for_best_formula)
            best_formula.append(nl_for_best_formula)
            best_formula.append(score_of_fragments)
            best_formula.append(score_of_nl)
            best_formula.append(score_of_molecule_list)
    except Exception as e_best_formula:
        print("Error in best_formula.append part of get_spectra_associated_with_mass. Error: ")
        print(str(e_best_formula))
        print(traceback.format_exc())
        return

    if len(best_formula) >= 1:
        with open(current_folder + "BEST_FORMULA_APPROXIMATION.txt", "w") as file:
            file.write("{:<30}".format("Best Formula: ") + str(best_formula[0]) + "\n")
            file.write("{:<30}".format("Score: ") + str(round(best_formula[1], 2)) + "\n")
            file.write("{:<30}".format("Intensity of Molecule: ") + str(round(best_formula[4], 1)) + "\n")
            file.write("_____________________________________________________________\n")
            file.write("{:<10}".format("Mass: ") + str(round(best_formula[2], 4)) + "\n")
            file.write("{:<10}".format("Index: ") + str(best_formula[3]) + "\n")
            file.write("_____________________________________________________________\n")
            file.write("===============FRAGMENTS===============\n")
            file.write(
                "{:<20}{:<20}{:<20}{:<20}{:<20}\n".format("Fragment", "Score", "NL", "NL Deviation", "Intensity"))
            for fragment_index in range(len(best_formula[5])):
                file.write("{:<20}".format(str(best_formula[5][fragment_index])))
                file.write("{:<20}".format(str(round(best_formula[8][fragment_index], 2))))
                file.write("{:<20}".format(str(best_formula[7][fragment_index])))
                file.write("{:<20}".format(str(round(best_formula[9][fragment_index], 2))))
                file.write("{:<20}\n".format(str(round(best_formula[6][fragment_index], 2))))

        if not summary_filepath == "":
            abs_path = os.path.abspath(summary_filepath)
            os.makedirs(os.path.dirname(abs_path), exist_ok=True)
            with open(summary_filepath, "a") as summary_file:
                # mass \t formula \t index \t rt \t intensity \t score \t [fragments] \t [fragments_score] \t [fragments_intensity] \t [nl] \t [nl_deviation]
                summary_file.write(str(round(best_formula[2], 4)) + "\t")
                summary_file.write(str(best_formula[0]) + "\t")
                summary_file.write(str(best_formula[3]) + "\t")
                summary_file.write(str(round(rt_list[index], 2)) + "\t")
                summary_file.write(str(round(best_formula[4], 1)) + "\t")
                summary_file.write(str(round(best_formula[1], 4)) + "\t")
                summary_file.write(str(best_formula[5]) + "\t")
                summary_file.write(str(best_formula[8]) + "\t")
                summary_file.write(str(best_formula[6]) + "\t")
                summary_file.write(str(best_formula[7]) + "\t")
                summary_file.write(str(best_formula[9]) + "\t")
                for element in additional_information_about_peak:
                    summary_file.write(str(element) + "\t")
                summary_file.write("\n")

    return [aif_processed, best_formula]


if __name__ == "__main__":

    mzml_filename = "U://MyFolder//MONOTONS//Filtermessungen//ACROSS Filter//Neue Messungen//Kalibrierung//Cal2.mzML"
    f = MS_functions.read_mzml_file(mzml_filename)
    output_folder = ".".join(mzml_filename.split(".")[:-1]) + "_outputfolder//"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    results_list = []
    mass_list = []
    peak_rt_list = []
    b_subst_area_list = []
    additional_info_lists = []
    with open(
            "U://MyFolder//MONOTONS//Filtermessungen//ACROSS Filter//Neue Messungen//Kalibrierung//Masses.csv",
            "r") as mass_file:
        lines = mass_file.readlines()
        for line in range(len(lines)):
            if line == 0:
                continue
            lines[line] = lines[line].replace("\n", "")
            entries = lines[line].split(";")
            peak_rt_list.append(float(entries[4]) * 60)
            mass_list.append(float(entries[3]))
            b_subst_area_list.append(float(entries[3]))
            additional_info_lists.append(entries)
    print(len(peak_rt_list))
    print(peak_rt_list)
    print(len(mass_list))
    print(mass_list)

    peak_dict = dict(zip(mass_list, peak_rt_list))
    peak_dict = {m: rt for m, rt in peak_dict.items() if rt < 1000}
    mass_list = list(peak_dict.keys())
    peak_rt_list = list(peak_dict.values())
    print(len(peak_dict))

    #mass_list = [154.0872]
    #peak_rt_list = [126.88]

    for entry in range(len(mass_list)):
        print("====================================================")
        print("Calculating for mass: " + str(mass_list[entry]) + " at RT: " + str(peak_rt_list[entry]))
        print("====================================================")
        additional_info_list = additional_info_lists[entry]
        try:
            result = one_analysis(mass_list[entry], peak_rt_list[entry], f, mass_deviation=0.001,
                                  # in u; for formula prediction -> where only the mass is considered
                                  requested_xic_mode="Full scan",
                                  # how should the xic be constructed?
                                  remove_if_combined_fragment_score_lower_than=20,
                                  # for creation of the plot. Remove the fragment, if the combined score is lower than -10
                                  charge_of_measured_mass=-1,
                                  # what is the charge of the measured mass (in neg. ion mode == -1, in pos ion mode == 1)
                                  remove_fragments_with_score_lower_than=125,
                                  # process aif spectrum; Fragment formula predictions where the score (isotope + mass deviation) is lower than this value, are removed and not considered for neutral loss formula prediction! There is no harm if this value is a little bit lower (-30 or -50 even).
                                  max_nl_deviation_ppm=7,  # orig 5
                                  # process aif spectrum; how big can the mass deviation in ppm be for a formula to be considered in further calculations?
                                  divisor_for_fragment_noise=300,  # orig 100
                                  # process aif spectrum; noise = intensity_of_requested_mass / divisor_for_fragment_noise, if the intensity of a peak is below the noise, the peak will not be considered.
                                  max_ppm_deviation=9,  # orig 15
                                  # Do not consider any formulas where the deviation to the measured mass is higher than this value.
                                  ppm_deviation_scoring_multiplier=14,  # orig 10
                                  # how big should the influence of the mass ppm deviation be on the final score of a formula prediction with isotope composition and mass deviation.
                                  noise_divisor_for_isotopologue_calculation=5,  # orig 1
                                  # 1 < variable < inf how big should the noise level be? Higher variable value -> lower noise level assumed -> more masses are considered in isotopologue calculation
                                  multiplier_isotopologue_influence_of_ppm_deviation_on_score=30,
                                  # the higher the variable is, the more influence does the ppm deviation has on the score of an isotopologue
                                  isotopologue_score_e_function_exponent=0.4,
                                  # the higher the value is, the more rapid the score decreases if the measured intensity is higher than the predicted intensity.
                                  summary_filepath=output_folder + "SUMMARY.txt",
                                  # where should the summary be saved?
                                  additional_information_about_peak=additional_info_list,
                                  # pass additional values that were standing in the mzmine file. They will also be saved to the summary file.
                                  nl_formula_cache_dict=None,  # neutral loss cache dict
                                  nl_formula_cache_dict_max_mass=1,
                                  # what is the mass that the cache dict was created with?
                                  fragment_formula_cache_dict=None,  # cache dict for fragments
                                  fragment_formula_cache_dict_max_mass=1,
                                  # what is the mass that the cache dict was created with?
                                  reject_predicted_formula_if_score_lower_than=10,
                                  # below what score should a formula be excluded from any further calculations?
                                  minimum_assumed_noise=7500,
                                  # what is the minimum noise that is considered for all calculations?
                                  add_value_to_score_if_isotopo_was_found=250,
                                  # score to be added to the combined score for each isotopologue that was found with its exact mass
                                  max_ppm_deviation_find_formula=25,
                                  # what is the maximum ppm deviation that may be applied to the given masses list to find a peak? Here you can set a high value.
                                  include_likelyhood_of_formula=True,
                                  # should the likelyhood of a formula be considered (somehow prevents the discovery of unknown unknowns)
                                  stop_isotopologue_search_if_score_lower_than=-80)  # stops the search for new isotopologues if a isotopologue with this score has been identified.

            results_list.append(result)
            plt.close("all")
            print(result)
        except ZeroDivisionError:
            pass

    print(results_list)
    with open(output_folder + "all_results_file.txt", "w") as results_file:
        for element in results_list:
            results_file.write(str(element))
            results_file.write("\n")
