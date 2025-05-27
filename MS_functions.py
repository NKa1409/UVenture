import ast
import copy
import sys
import traceback
import similaritymeasures
import scipy
import matplotlib
matplotlib.use('Agg')
import gzip
import numpy as np




def read_mzml_file(mzml_filename):
    from pyteomics import mzml
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


def get_formula_from_cache(cache_folder, mass, max_ppm_dev, max_unlikeliness_of_formula=10, min_unlikelyness_of_formula=0, file_basename="formulas_for_mass_", filetype=".txt.gz", max_mass=1200, remove_unlogical_formulas=False, debug_output=True):
    if debug_output == True:
        print("Starting formula generation from cache...")
    if mass >= max_mass:
        #on the fly prediction
        result = get_one_formula(mass, max_dev=0.005, charge_of_measured_mass=0, dbe_range=(-10, 150), c_oxidation_state_range=(-4, 4),
                    heteroatom_ratio_range=(0, 99999), aromaticity_eq_range=(-20, 100), h_c_ratio=(0, 4),
                    max_h=100, max_c=30, max_n=30, max_o=20, max_s=3, max_p=2, max_cl=4, max_br=2, max_i=0, print_progress=True, return_cache_dict=False)
        pred_dict = result[1]

        #lookup in cache files
        files_to_search = [cache_folder + file_basename + str(int(mass)) + filetype,
                           cache_folder + file_basename + str(int(mass) - 1) + filetype,
                           cache_folder + file_basename + str(int(mass) + 1) + filetype]
        all_dict = {}
        for file in files_to_search:
            try:
                with open(file, "r") as f:
                    lines = f.readlines()
                    formulas = []
                    deviations = []
                    for line_index in range(len(lines)):
                        lines[line_index] = lines[line_index].strip()
                        if lines[line_index] == "":
                            continue
                        formulas.append(lines[line_index].split(",")[0])
                        deviations.append((((float(lines[line_index].split(",")[-1])) - mass) / mass) * 1000000)
                    mydict = dict(zip(formulas, deviations))
                all_dict.update(mydict)
            except Exception as e:
                print("Error with formula determination! " + str(e))
                print(traceback.format_exc())
                if debug_output == True:
                    print(all_dict)
                continue
        all_dict = {f: dev for f, dev in all_dict.items() if abs(dev) <= max_ppm_dev}
        all_dict = {k: v for k, v in sorted(all_dict.items(), key=lambda item: abs(item[1]), reverse=False)}
        all_dict.update(pred_dict)
    else:
        files_to_search = [cache_folder + file_basename + str(int(mass)) + filetype,
                           cache_folder + file_basename + str(int(mass)-1) + filetype,
                           cache_folder + file_basename + str(int(mass)+1) + filetype]
        all_dict = {}
        for file in files_to_search:
            try:
                with gzip.open(file, "rt") as f:
                    lines = f.readlines()
                    formulas = []
                    deviations = []
                    for line_index in range(len(lines)):
                        lines[line_index] = lines[line_index].strip()
                        if lines[line_index] == "":
                            continue
                        formulas.append(lines[line_index].split(",")[0])
                        deviations.append( (((float(lines[line_index].split(",")[-1])) - mass) / mass) * 1000000 )
                    mydict = dict(zip(formulas, deviations))
                all_dict.update(mydict)
            except Exception as e:
                print("Error with formula determination! " + str(e))
                print(traceback.format_exc())
                if debug_output == True:
                    print(all_dict)
                continue
        all_dict = {f: dev for f, dev in all_dict.items() if abs(dev) <= max_ppm_dev}
        all_dict = {k: v for k, v in sorted(all_dict.items(), key=lambda item: abs(item[1]), reverse=False)}
    formula_dict_list = [get_formula_to_dict(f) for f, dev in all_dict.items()]
    formula_dict_list = [dict(t) for t in {tuple(d.items()) for d in formula_dict_list}]
    formula_dict_list = [get_formula_string_from_dict(f) for f in formula_dict_list]
    all_dict = {f: dev for f, dev in all_dict.items() if f in formula_dict_list}

    if remove_unlogical_formulas == True:
        for formula in list(all_dict.keys()):
            curr_formula_dict = get_formula_to_dict(formula)
            f_likelyness = calculate_likelyhood_of_formula_dict(curr_formula_dict)
            f_unlikelyness = -1 * f_likelyness
            if f_unlikelyness >= max_unlikeliness_of_formula or f_unlikelyness < min_unlikelyness_of_formula:
                if debug_output == True:
                    print("Unlikelyness of current formula: " + str(f_unlikelyness) + "Formula: " + str(formula))
                del all_dict[formula]

    return [mass, all_dict]


def get_one_formula(measured_mass, max_dev=0.005, charge_of_measured_mass=-1, dbe_range=(-1, 15), c_oxidation_state_range=(-4, 4),
                    heteroatom_ratio_range=(0, 0.75), aromaticity_eq_range=(-20, 100), h_c_ratio=(0.3, 3),
                    max_h=100, max_c=30, max_n=30, max_o=20, max_s=5, max_p=3, max_cl=5, max_br=1, max_i=1, print_progress=True, return_cache_dict=False):
    mass_dict = {"H": 1.007825,
                 "B": 11.009305,
                 "C": 12.00000,
                 "N": 14.003074,
                 "O": 15.994915,
                 "F": 18.998403,
                 "Cl": 34.968853,
                 "P": 30.973763,
                 "S": 31.972072,
                 "Br": 80.916290,
                 "I": 126.904477,
                 "Fe": 55.934939,
                 "Si": 27.976928,
                 "As": 74.921596,
                 "Se": 79.916521}
    cache_dict = {}

    if not charge_of_measured_mass == 0:
        mass = abs(measured_mass * charge_of_measured_mass)
    else:
        mass = measured_mass

    mass = mass + (charge_of_measured_mass * 0.0005485)

    h_atoms = [0, max_h]
    c_atoms = [0, max_c]
    n_atoms = [0, max_n]
    o_atoms = [0, max_o]
    s_atoms = [0, max_s]
    p_atoms = [0, max_p]
    cl_atoms = [0, max_cl]
    br_atoms = [0, max_br]
    i_atoms = [0, max_i]
    dbe_range = dbe_range
    max_deviation_mz = max_dev
    formulas = {}
    if c_atoms[-1] > int(mass / mass_dict["C"]) + 1:
        c_atoms[-1] = int(mass / mass_dict["C"])
    if n_atoms[-1] > int(mass / mass_dict["N"]) + 1:
        n_atoms[-1] = int(mass / mass_dict["N"])
    if o_atoms[-1] > int(mass / mass_dict["O"]) + 1:
        o_atoms[-1] = int(mass / mass_dict["O"])
    if s_atoms[-1] > int(mass / mass_dict["S"]) + 1:
        s_atoms[-1] = int(mass / mass_dict["S"])
    if p_atoms[-1] > int(mass / mass_dict["P"]) + 1:
        p_atoms[-1] = int(mass / mass_dict["P"])
    if cl_atoms[-1] > int(mass / mass_dict["Cl"]) + 1:
        cl_atoms[-1] = int(mass / mass_dict["Cl"])
    if br_atoms[-1] > int(mass / mass_dict["Br"]) + 1:
        br_atoms[-1] = round(mass / mass_dict["Br"])
    if i_atoms[-1] > int(mass / mass_dict["I"]) + 1:
        i_atoms[-1] = round(mass / mass_dict["I"])
    if h_atoms[-1] > int(mass / mass_dict["H"]) + 1:
        h_atoms[-1] = int(mass / mass_dict["H"])
    if 0.25 <= (mass - int(mass)) <= 0.75:
        return [mass, formulas]
    for c_count in range(c_atoms[0], c_atoms[-1] + 1, 1):
        if print_progress:
            print(str(c_count) + " / " + str(c_atoms[-1]))
        for n_count in range(n_atoms[0], n_atoms[-1] + 1, 1):
            for br_count in range(br_atoms[0], br_atoms[1]+1, 1):
                for i_count in range(i_atoms[0], i_atoms[1]+1, 1):
                    for o_count in range(o_atoms[0], o_atoms[-1] + 1, 1):
                        for s_count in range(s_atoms[0], s_atoms[-1] + 1, 1):
                            for p_count in range(p_atoms[0], p_atoms[-1] + 1, 1):
                                for cl_count in range(cl_atoms[0], cl_atoms[-1] + 1, 1):
                                    for h_count in range(h_atoms[0], 1 + int(2 * (c_count + (2 / 3 * n_count) + (1 / 2 * s_count) + (1 / 2 + o_count)) + 2) + 1, 1):
                                        halogen_count = 0
                                        heavy_atom_count = 0
                                        if c_count == 0:
                                            current_heteroatom_ratio = 999
                                            current_cos = 0
                                            current_h_c_ratio = 999
                                        else:
                                            current_heteroatom_ratio = (n_count + o_count + s_count + p_count + cl_count + br_count + i_count) / (c_count + (n_count + o_count + s_count + p_count + cl_count + br_count + i_count))
                                            current_cos = 2 * (o_count / c_count) - (h_count / c_count)  # see https://doi.org/10.1016/j.atmosenv.2018.06.036
                                            current_h_c_ratio = h_count / c_count

                                        dbe = (c_count + 1) - ((h_count - n_count + halogen_count) / 2)
                                        try:
                                            aromaticity_equivalent = ((c_count - (h_count - c_count)) / (
                                                dbe)) + 1  # see https://doi.org/10.1002/rcm.7038
                                        except ZeroDivisionError:
                                            aromaticity_equivalent = 0
                                        if (dbe_range[0] <= dbe <= dbe_range[-1] and
                                                (c_oxidation_state_range[0] <= current_cos <= c_oxidation_state_range[-1]) and
                                                (heteroatom_ratio_range[0] <= current_heteroatom_ratio <= heteroatom_ratio_range[-1]) and
                                                (aromaticity_eq_range[0] <= aromaticity_equivalent <= aromaticity_eq_range[-1]) and
                                                (h_c_ratio[0] <= current_h_c_ratio <= h_c_ratio[-1])):
                                            calc_mass = h_count * mass_dict["H"] + \
                                                        c_count * mass_dict["C"] + \
                                                        n_count * mass_dict["N"] + \
                                                        o_count * mass_dict["O"] + \
                                                        s_count * mass_dict["S"] + \
                                                        p_count * mass_dict["P"] + \
                                                        cl_count * mass_dict["Cl"] + \
                                                        br_count * mass_dict["Br"] + \
                                                        i_count * mass_dict["I"]
                                            if calc_mass >= mass+2:
                                                continue
                                            if return_cache_dict:
                                                print_str = ""
                                                if c_count > 0:
                                                    print_str = print_str + "C" + str(c_count)
                                                if h_count > 0:
                                                    print_str = print_str + "H" + str(h_count)
                                                if n_count > 0:
                                                    print_str = print_str + "N" + str(n_count)
                                                if o_count > 0:
                                                    print_str = print_str + "O" + str(o_count)
                                                if s_count > 0:
                                                    print_str = print_str + "S" + str(s_count)
                                                if p_count > 0:
                                                    print_str = print_str + "P" + str(p_count)
                                                if cl_count > 0:
                                                    print_str = print_str + "Cl" + str(cl_count)
                                                if br_count > 0:
                                                    print_str = print_str + "Br" + str(br_count)
                                                if i_count > 0:
                                                    print_str = print_str + "I" + str(i_count)
                                                cache_dict[print_str] = calc_mass
                                            if calc_mass - max_deviation_mz <= mass <= calc_mass + max_deviation_mz:
                                                # calculate score
                                                score = 0
                                                ppm_deviation = ((calc_mass - mass) / mass) * 1000000
                                                if ppm_deviation == 0:
                                                    score = 99999
                                                else:
                                                    try:
                                                        score = (1 / ppm_deviation)
                                                    except ZeroDivisionError:
                                                        score = 999999
                                                print_str = ""
                                                if c_count > 0:
                                                    print_str = print_str + "C" + str(c_count)
                                                if h_count > 0:
                                                    print_str = print_str + "H" + str(h_count)
                                                if n_count > 0:
                                                    print_str = print_str + "N" + str(n_count)
                                                if o_count > 0:
                                                    print_str = print_str + "O" + str(o_count)
                                                if s_count > 0:
                                                    print_str = print_str + "S" + str(s_count)
                                                if p_count > 0:
                                                    print_str = print_str + "P" + str(p_count)
                                                if cl_count > 0:
                                                    print_str = print_str + "Cl" + str(cl_count)
                                                if br_count > 0:
                                                    print_str = print_str + "Br" + str(br_count)
                                                if i_count > 0:
                                                    print_str = print_str + "I" + str(i_count)
                                                formulas[print_str] = ppm_deviation
                                        else:
                                            pass

    formulas = {k: v for k, v in sorted(formulas.items(), key=lambda item: abs(item[1]), reverse=False)}
    if len(formulas) == 0:
        formulas = {"": 99999999}
    if return_cache_dict:
        cache_dict = {k: v for k, v in sorted(cache_dict.items(), key=lambda item: abs(item[1]), reverse=False)}
        return cache_dict
    return [mass, formulas]


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


def calculate_likelyhood_of_formula_dict(f_dict, charge_of_measured_mass=None, score_subst_rdbe_non_integer=10, score_subst_senior_rule=10, debug_output=False):
    formula_likelyness = 0
    try:
        # X/C ratio check from #https://pmc.ncbi.nlm.nih.gov/articles/PMC1851972/ The numbers are representing 99.7% of all available molecular formulas.
        if "C" in list(f_dict.keys()):
            if (f_dict.get("H", 0) / f_dict.get("C", 0) < 0.2) or (f_dict.get("H", 0) / f_dict.get("C", 0) > 3.1):
                formula_likelyness = formula_likelyness -40
            if (f_dict.get("F", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("F", 0) / f_dict.get("C", 0) > 1.5):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("Cl", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("Cl", 0) / f_dict.get("C", 0) > 0.8):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("Br", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("Br", 0) / f_dict.get("C", 0) > 0.8):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("N", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("N", 0) / f_dict.get("C", 0) > 1.3):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("O", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("O", 0) / f_dict.get("C", 0) > 3):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("P", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("P", 0) / f_dict.get("C", 0) > 0.3):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("S", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("S", 0) / f_dict.get("C", 0) > 0.8):
                formula_likelyness = formula_likelyness - 40
            if (f_dict.get("Si", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("Si", 0) / f_dict.get("C", 0) > 0.5):
                formula_likelyness = formula_likelyness - 40

        heteroatom_count = float(f_dict.get("O", 0)) + float(f_dict.get("N", 0)) + float(f_dict.get("S", 0)) + float(f_dict.get("P", 0))
        if heteroatom_count > 5 and f_dict.get("C", 0) <= int(heteroatom_count / 4):
            formula_likelyness = formula_likelyness - float((heteroatom_count / 4) - f_dict.get("C", 0)) * 50

        if ("C" in list(f_dict.keys())) and ("N" in list(f_dict.keys())):
            if int(f_dict["N"]) > 2 and (int(f_dict["C"]) / int(f_dict["N"]) <= 4):
                formula_likelyness = formula_likelyness - ((float(f_dict["N"]) - 2) ** 2) * 50

        if ("C" in list(f_dict.keys())) and ("H" in list(f_dict.keys())):
            if (float(f_dict["H"]) / float(f_dict["C"]) >= 2):
                formula_likelyness = formula_likelyness - (((float(f_dict["H"]) / float(f_dict["C"])) - 2) ** 3) * 50
                formula_likelyness = formula_likelyness - ((float(f_dict.get("N", 0)) ** 2) * 50)
                formula_likelyness = formula_likelyness - ((float(f_dict.get("N", 0)) ** 2) * 50)

        if ("P" in list(f_dict.keys())):
            if (float(f_dict.get("P", 0)) * 3.1) >= float(f_dict.get("O", 0)):
                formula_likelyness = formula_likelyness - (((float(f_dict.get("P", 0)) * 3) / float(f_dict.get("O", 1))) ** 4) * 250

        if 2 < f_dict.get("C", 0) < heteroatom_count:
            formula_likelyness = formula_likelyness - ((heteroatom_count - float(f_dict.get("C", 0))) ** 2) * 50

        dbe = calc_dbe(f_dict)

        if dbe < 0:
            score_substract = (abs(dbe + 2) * 50) ** 3
        elif (dbe - f_dict.get("O", 0)) > 7:
            score_substract = abs(dbe - f_dict.get("O", 0) - 7) * 50
        else:
            score_substract = 0
        if not charge_of_measured_mass == None:
            if charge_of_measured_mass < 0:
                dbe = dbe + (charge_of_measured_mass * 0.5)
        if dbe - 1 > (get_mass_of_most_abundant_isotopologue_formula(f_dict) * (62 / 1000)):  # Senior Rule
            score_substract += score_subst_senior_rule
            if debug_output == True:
                print("Senior rule of formula: "+ str(f_dict) + " not fulfilled.")
        if not dbe - int(dbe) == 0:
            score_substract += score_subst_rdbe_non_integer
            if debug_output == True:
                print("RDBE rule of formula: " + str(f_dict) + " not fulfilled.")

        # https://pmc.ncbi.nlm.nih.gov/articles/PMC1851972/
        # Element ratios

        formula_likelyness = formula_likelyness - score_substract
    except Exception as e:
        print("ERROR in formula score likelyhood: " + str(e))
        print(traceback.format_exc())
        formula_likelyness = 0

    return formula_likelyness


def get_mass_of_most_abundant_isotopologue_formula(formula):
    if len(formula) == 0:
        return None
    import pyteomics.mass
    isotopologues = pyteomics.mass.mass.isotopologues(formula)
    isotopologues_list = []
    while True:
        try:
            isotopologues_list.append(next(isotopologues))
        except:
            break
    abundance_dict = {}
    for entry in isotopologues_list:
        abundance_dict[pyteomics.mass.mass.calculate_mass(entry)] = pyteomics.mass.mass.isotopic_composition_abundance(entry)
    abundance_dict = {m: a for m, a in abundance_dict.items() if a > 0.0001}
    return list(abundance_dict.keys())[0]


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


def simulate_isotope_pattern_of_formula(formula, mass_resolution_ppm=10, debug_output=False):
    if len(formula) == 0:
        return None
    import pyteomics.mass
    isotopologues = pyteomics.mass.mass.isotopologues(formula)
    isotopologues_list = []
    while True:
        try:
            isotopologues_list.append(next(isotopologues))
        except:
            break
    abundance_dict = {}
    for entry in isotopologues_list:
        abundance_dict[pyteomics.mass.mass.calculate_mass(entry)] = pyteomics.mass.mass.isotopic_composition_abundance(entry)
    abundance_dict = {m: a for m, a in abundance_dict.items() if a > 0.0001}
    if debug_output == True:
        print("Abundance dict for formula " + str(formula) + ":  " + str(abundance_dict))
    try:
        abundance_dict = summarize_mass_intensity_dict_for_isotopo_simulation(abundance_dict, deviation=mass_resolution_ppm, debug_output=debug_output)
        abundance_dict = dict(sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True))
    except Exception as e:
        print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
        print(traceback.format_exc())
        abundance_dict = dict(sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True))
    
    return abundance_dict


def get_formula_to_dict(formula_string):
    if formula_string == "":
        return {}
    if not formula_string[-1].isnumeric():
        formula_string += "1"
    new_string = ""
    for position in range(len(formula_string)):
        if (not formula_string[position].isnumeric()) and (not formula_string[position+1].isnumeric()) and (formula_string[position+1].isupper()):
            new_string += formula_string[position] + "1"
        else:
            new_string += formula_string[position]
    formula_string = new_string

    atom_dict = {}
    for entry in range(len(formula_string)):
        atom_name = ""
        atom_count = "0"
        iterator = 0
        if not formula_string[entry - 1].isnumeric():
            continue
        if formula_string[entry].isnumeric():
            continue
        if (not formula_string[entry].isnumeric()) and (formula_string[entry + 1].isnumeric()):
            atom_name = formula_string[entry]
            try:
                while formula_string[entry + 1 + iterator].isnumeric():
                    atom_count = atom_count + formula_string[entry + 1 + iterator]
                    iterator = iterator + 1
            except:
                pass
        elif (not formula_string[entry].isnumeric()) and (not formula_string[entry + 1].isnumeric()):
            atom_name = formula_string[entry] + formula_string[entry + 1]
            try:
                while formula_string[entry + 2 + iterator].isnumeric():
                    atom_count = atom_count + formula_string[entry + 2 + iterator]
                    iterator = iterator + 1
            except:
                pass
        if int(atom_count) >= 1:
            atom_dict[atom_name] = int(atom_count)
    return atom_dict


def get_formula_string_from_dict(formula_dict):
    if isinstance(formula_dict, str):
        try:
            eventual_formula_dict = ast.literal_eval(formula_dict)
            if isinstance(eventual_formula_dict, dict):
                formula_dict = eventual_formula_dict
            else:
                pass
        except:
            pass
    if isinstance(formula_dict, str):
        try:
            formula_dict = get_formula_to_dict(formula_dict)
        except:
            pass
    if len(formula_dict) == 0:
        return ""
    return "".join([str(a) + str(n) for a, n in formula_dict.items()])


def read_summary_to_df(filepath, sep="\t", colnames=("MASS", "FORMULA", "INDEX", "RT", "INTENSITY", "SCORE", "FRAGS", "F_SCORES", "F_INTENSITIES", "NLS", "NLS_DEV"), header=None):
    import pandas as pd
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


def calc_dbe(formula):
    try:
        if isinstance(formula, str):
            formula = get_formula_to_dict(formula)
        if isinstance(formula, dict):
            pass
        dbe = ((formula.get("Si", 0) * 2) + (formula.get("C", 0)*2) + 2 + formula.get("N", 0) + formula.get("P", 0) - formula.get("H", 0) - formula.get("F", 0) - formula.get("Cl", 0) - formula.get("Br", 0) - formula.get("I", 0) ) / 2
        return dbe
    except Exception as e:
        print("Error calculating DBE in MS_functions: " + str(e))
        print(traceback.format_exc())
        return 0



if __name__ == "__main__":
    pattern = simulate_isotope_pattern_of_formula("C3H11N1O1S2Br1")
    print(pattern)
    sys.exit()


################################################################################################################

    dictio = {221.94504457713: 0.4669847312210661,
              219.94709107713: 0.43961483934980194,
              223.94084047713: 0.03914717156374999,
              222.94839941493: 0.015152339909224929,
              220.95044591493001: 0.014264263967581771,
              222.94443233713: 0.00834003332386434,
              220.94647883713: 0.006942017675804853,
              222.94207947053: 0.002016631029850584,
              220.94412597053: 0.0016060440154495154,
              224.94419531493: 0.001241661217539753,
              221.95133745757: 0.001128654627464162,
              223.94929095756999: 0.0010979277901017523,
              225.93663637713001: 0.0008560657853599821,
              224.94022823713: 0.00030214086542117006,
              221.95380075273002: 0.0001542784033691751,
              223.95175425273: 0.00015007828112318058,
              224.93787537052998: 0.00013980129449213708}

    print(summarize_mass_intensity_dict(dictio))

    sys.exit()



    # start_time1 = time.time()
    # cache_dict = get_one_formula(150, return_cache_dict=True, charge_of_measured_mass=0)
    # end_time1 = time.time()
    # print(len(cache_dict))
    # start_time2 = time.time()
    # formulas = get_formula_from_cache_dict(cache_dict, 138.0191, 0.005)
    # end_time2 = time.time()
    # print(len(formulas))
    # print(formulas)
    # print("Time1: " + str(end_time1-start_time1))
    # print("Time2: " + str(end_time2 - start_time2))
    # formula = get_one_formula(138.0196, max_dev=0.005, charge_of_measured_mass=-1)
    # print(formula[1])
    # input()

    mzml_filename = "C://Users//Admin//Desktop//UVenture//UVenture//webserver_save//mzml_files//F12_HRAIF4_3.mzML"
    f = read_mzml_file(mzml_filename)
    starttime = time.time()
    xic = get_xic(f, mass=138.0191, mass_deviation=.001, requested_filter_mode="Full scan", starttime=starttime)
    print("Time: " + str(time.time()-starttime))
    fig = matplotlib.figure.Figure()
    ax = fig.subplots()
    ax.plot(xic[0], xic[1])
    print("Time: " + str(time.time()-starttime))
    plt.show()
    
    input()
#######################################################################################################################
    #returns {"C": 3, "H": 4, "N": 2, "F": 1}
    formula_dict = get_formula_to_dict("C3H4N2F1")

    #returns "C3H4N2F1"
    formula_string = get_formula_string_from_dict(formula_dict)
#######################################################################################################################


    # returns{mass_with_highest_abundance_for_formula: abundance, mass_with_second_highest_abundance: abundance, mass_with_third_highest_..........}
    pattern = simulate_isotope_pattern_of_formula("C3H7")
#######################################################################################################################

    # returns [rt_list, intensity_list, list_of_original_indices] -> as only full scan intensities are plotted,
    # only those indices where a full scan was measured are plotted.
    # So the original index of each datapoint is also given
    xic = get_xic(f, mass=138.0191, mass_deviation=.001, requested_filter_mode="Full scan")
    fig = matplotlib.figure.Figure()
    ax = fig.subplots()
    ax.plot(xic[0], xic[1])
    plt.show()
#######################################################################################################################

    #returns [masses_list, intensities_list, index_of_spectrum, all_information_contained_in_file_for_index, ms_level, filter, filter_mode, measured_msms_masses]

    # as_mass first creates the bg_substracted xic for the requested_mass, looks if it finds any peaks for the requested_mass
    # and then returns the spectrum that is closest to the identified peak and has the requested filter mode.

    # as_index returns the spectrum that is closest to the requested index and has the requested filter mode.

    # as_rt returns the spectrum that is closest to the requested retention time and has the requested filter mode

    # requested_filter_mode = "Full scan", "AIF", "MS/MS", "whatever" --> "whatever" only works for as_rt and as_index.
    # it returns the exact spectrum that is first found. No filtering based on requested filter mode.

    # if save_folder="" no logfile and no plot of the spectrum will be saved
    # if save_folder="path_to_folder//" a logfile and a plot of the spectrum will be saved there.

    spec = get_mass_spectrum(f, 138.0191, mode="as_mass", mass_deviation=0.005, requested_filter_mode="Full scan")
    fig = matplotlib.figure.Figure()
    ax = fig.subplots()
    ax.bar(spec[0], spec[1])
    plt.show()
#######################################################################################################################

    #returns [mass, {formula_approx1: ppm_deviation1, formula_approx2: ppm_deviation2, .....}]
    formula_approximations = get_one_formula(138.0191, max_dev=0.003, charge_of_measured_mass=-1, dbe_range=(-1, 15), c_oxidation_state_range=(-4, 4),
                    heteroatom_ratio_range=(0, 0.75), aromaticity_eq_range=(-20, 100), h_c_ratio=(0.3, 3),
                    max_h=100, max_c=100, max_n=100, max_o=100, max_s=0, max_p=0, max_cl=0, print_progress=False)
#######################################################################################################################





