import copy
import traceback
from functools import lru_cache
import matplotlib.pyplot as plt
import scipy
import PIL
import matplotlib
import pathlib
import time




def read_mzml_file(mzml_filename):
    from pyteomics import mzml
    f = mzml.read(mzml_filename)
    f = list(f)
    return f


def get_xic_fast(f, mass, mass_deviation):
    def is_within_deviation(m):
        return mass-mass_deviation <= m <= mass+mass_deviation
    rt_list = [element["scanList"]["scan"][0]["scan time"] for element in f]
    intensity_list = []
    for entry in f:
        try:
            intensity_list.append(sum([entry["intensity array"][i] for i in range(len(entry["m/z array"])) if is_within_deviation(entry["m/z array"][i])]))
        except:
            intensity_list.append(0)
    return [rt_list, intensity_list]

def get_xic(f, mass, mass_deviation, requested_filter_mode="Full scan"):
    def is_within_deviation(m):
        return mass-mass_deviation <= m <= mass+mass_deviation
    if (not requested_filter_mode == "Full scan") and (not requested_filter_mode == "AIF") and (not requested_filter_mode == "MS/MS"):
        requested_filter_mode = "Full scan"
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    intensity_list = []
    counter = 0
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
        counter = counter + 1
        #print(str(counter) + "/" + str(len(f)))
        y = (i for i,v in enumerate(entry["m/z array"]) if is_within_deviation(v))
        mass_indices = []
        while True:
            try:
                mass_indices.append(next(y))
            except:
                break
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


def get_mass_spectrum(f, value, mode="as_index", mass_deviation=0.005, requested_filter_mode="Full scan",
                      requested_ms_ms_mass="", save_folder=""):
    # FTMS - p APCI corona d Full ms2 138.0197@hcd29.67 [50.0000-160.0000]   --> MSMS
    # FTMS - p APCI corona Full ms [50.0000-500.0000]   --> Full scan
    # FTMS - p APCI corona Full ms2 225.0000@hcd29.67 [50.0000-400.0000]   --> AIF

    # requested_filter_mode="Full scan" / "AIF" / "MS/MS" / "whatever"
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
                               int(peaks[peak] + peak_properties["widths"][peak] * not_including_peak_width_multiplier),
                               1):

                try:
                    m = (background[
                             int(peaks[peak] + peak_properties["widths"][peak] * not_including_peak_width_multiplier)] -
                         background[
                             int(peaks[peak] - peak_properties["widths"][
                                 peak] * not_including_peak_width_multiplier)]) / (int(
                        peaks[peak] + peak_properties["widths"][peak] * not_including_peak_width_multiplier) - int(
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

    def is_within_deviation(curr_mass, requested_mass, mass_deviation):
        return requested_mass - mass_deviation <= curr_mass <= requested_mass + mass_deviation

    if requested_filter_mode == "AIF":
        requested_ms_ms_mass = ""
    if requested_filter_mode == "Full scan":
        requested_ms_ms_mass = ""
    if (mode == "as_mass") and (requested_ms_ms_mass == "") and (requested_filter_mode == "MS/MS"):
        requested_ms_ms_mass = value
    ms_ms_masses = []
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    log_f_filename = ""
    if not save_folder == "":
        log_f_filename = save_folder + "Mass_spectrum_index" + str(value) + "_AdditionalInfo" + ".txt"
        log_f = open(log_f_filename, "w")
        log_f.write("given_value\t" + str(value) + "\n")
        log_f.write("requested_filter\t" + str(requested_filter_mode) + "\n")
        log_f.write("mode\t" + str(mode) + "\n")
        log_f.write("requested_ms_ms_mass\t" + str(requested_ms_ms_mass) + "\n")
        log_f.write("STARTING....................................................." + "\n")
        log_f.close()
    index = 0
    start_index = 0
    # mode == as_index, as_rt, as_mass
    if mode == "as_rt":
        change_index_value = 1
        change_index_direction = "+"
        index = 0

        abs_err_list = []
        for element in rt_list:
            abs_err_list.append(abs(element - value))
        index = abs_err_list.index(min(abs_err_list))
        start_index = copy.deepcopy(index)

        filter_mode = ""
        while True:
            if requested_filter_mode == "whatever":
                break
            filter = f[index]["scanList"]["scan"][0]["filter string"]
            filter_mode = ""
            ms_ms_masses = []
            if (" d " in filter) and ("@hcd" in filter):
                filter_mode = "MS/MS"
                ms_ms_masses = []
                filter_parsed = filter.split(" ")
                filter_parsed = [x for x in filter_parsed if "hcd" in x]
                for element in filter_parsed:
                    ms_ms_masses.append(round(float(element.split("@")[0]), 2))
            if (not " d " in filter) and (not "hcd" in filter):
                filter_mode = "Full scan"
            if (not " d " in filter) and ("hcd" in filter):
                filter_mode = "AIF"

            if (index >= len(f)-3) or (index <= 3):
                if not save_folder == "":
                    log_f = open(log_f_filename, "a")
                    log_f.write(
                        "INFO:\t" + "NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning..." + "\n")
                    log_f.close()
                print("NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...")
                return
            if not filter_mode == requested_filter_mode:
                if change_index_direction == "+":
                    index = index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    index = index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1
                continue
            if filter_mode == requested_filter_mode and requested_ms_ms_mass == "":
                break
            if (filter_mode == requested_filter_mode) and (round(requested_ms_ms_mass, 2) in ms_ms_masses):
                if not save_folder == "":
                    log_f = open(log_f_filename, "a")
                    log_f.write("INFO:\t" + "FOUND MSMSSPECTRUM" + "\n")
                    log_f.close()
                print("FOUND MSMSSPECTRUM")
                break
            else:
                if change_index_direction == "+":
                    index = index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    index = index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1

    if mode == "as_index":
        print("Value: " + str(value))
        change_index_value = 1
        change_index_direction = "+"
        index = value
        filter = f[index]["scanList"]["scan"][0]["filter string"]
        filter_mode = ""
        ms_ms_masses = []
        if (" d " in filter) and ("@hcd" in filter):
            filter_mode = "MS/MS"
            filter_parsed = filter.split(" ")
            filter_parsed = [x for x in filter_parsed if "hcd" in x]
            for element in filter_parsed:
                ms_ms_masses.append(round(float(element.split("@")[0]), 2))
        if (not " d " in filter) and (not "hcd" in filter):
            filter_mode = "Full scan"
        if (not " d " in filter) and ("hcd" in filter):
            filter_mode = "AIF"

        if (not filter_mode == requested_filter_mode) and (not requested_filter_mode == "whatever"):
            if not save_folder == "":
                log_f = open(log_f_filename, "a")
                log_f.write(
                    "INFO:\t" + "Filter mode of the provided index does not match the requested filter mode! Adapting index to nearest spectrum that has requested filter mode!" + "\n")
                log_f.close()
            print(
                "Filter mode of the provided index does not match the requested filter mode! Adapting index to nearest spectrum that has requested filter mode!")

        start_index = copy.deepcopy(index)

        while True:
            if requested_filter_mode == "whatever":
                break
            filter = f[index]["scanList"]["scan"][0]["filter string"]
            filter_mode = ""
            ms_ms_masses = []
            if (" d " in filter) and ("@hcd" in filter):
                filter_mode = "MS/MS"
                ms_ms_masses = []
                filter_parsed = filter.split(" ")
                filter_parsed = [x for x in filter_parsed if "hcd" in x]
                for element in filter_parsed:
                    ms_ms_masses.append(round(float(element.split("@")[0]), 2))
            if (not " d " in filter) and (not "hcd" in filter):
                filter_mode = "Full scan"
            if (not " d " in filter) and ("hcd" in filter):
                filter_mode = "AIF"

            if (index >= len(f)-3) or (index <= 3):
                if not save_folder == "":
                    log_f = open(log_f_filename, "a")
                    log_f.write(
                        "INFO:\t" + "NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning..." + "\n")
                    log_f.close()
                print("NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...")
                return
            if not filter_mode == requested_filter_mode:
                if change_index_direction == "+":
                    index = index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    index = index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1
                continue
            if filter_mode == requested_filter_mode and requested_ms_ms_mass == "":
                break
            if (filter_mode == requested_filter_mode) and (round(requested_ms_ms_mass, 2) in ms_ms_masses):
                if not save_folder == "":
                    log_f = open(log_f_filename, "a")
                    log_f.write("INFO:\t" + "FOUND MSMSSPECTRUM" + "\n")
                    log_f.close()
                print("FOUND MSMSSPECTRUM")
                break
            else:
                if change_index_direction == "+":
                    index = index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    index = index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1

    if mode == "as_mass":
        change_index_value = 1
        change_index_direction = "+"
        xic = get_xic(f, value, mass_deviation, requested_filter_mode="Full scan")
        dur = (60 * max(rt_list)) / len(xic[1])
        bg_subst_int = do_bg_substraction(xic[1], dur)
        peaks = scipy.signal.find_peaks(bg_subst_int[0], height=max(bg_subst_int[0]) / 10)
        if not save_folder == "":
            fig, ax = plt.subplots()
            ax.plot(xic[0], bg_subst_int[0], label="BG subst. intensity")
            ax.plot(xic[0], xic[1], label="Intensity")
            ax.plot(xic[0], bg_subst_int[1], label="Background")
            plt.legend()
            ax.set_xlabel("RT / min")
            ax.set_ylabel("intensity / a.u.")
            ax.set_title("XIC of mass " + str(value))
            matplotlib.rcParams.update({'figure.autolayout': True})
            image_filepath = save_folder + "XIC_forMass" + str(index) + ".png"
            plt.savefig(image_filepath, bbox_inches='tight')

        if len(peaks[0]) >= 2:
            print(
                "ATTENTION: There are more than 1 peak for this mass!!! The program will take the highest value it can find in the xic!!")
            print("Press enter to continue...")
            input()
        index = xic[1].index(max(xic[1]))
        index = xic[2][index]
        filter_mode = ""

        start_index = copy.deepcopy(index)

        while True:
            if requested_filter_mode == "whatever":
                break
            filter = f[index]["scanList"]["scan"][0]["filter string"]
            filter_mode = ""
            ms_ms_masses = []
            if (" d " in filter) and ("@hcd" in filter):
                filter_mode = "MS/MS"
                ms_ms_masses = []
                filter_parsed = filter.split(" ")
                filter_parsed = [x for x in filter_parsed if "hcd" in x]
                for element in filter_parsed:
                    ms_ms_masses.append(round(float(element.split("@")[0]), 2))
            if (not " d " in filter) and (not "hcd" in filter):
                filter_mode = "Full scan"
            if (not " d " in filter) and ("hcd" in filter):
                filter_mode = "AIF"

            if (index >= len(f)-3) or (index <= 3):
                if not save_folder == "":
                    log_f = open(log_f_filename, "a")
                    log_f.write(
                        "INFO:\t" + "NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning..." + "\n")
                    log_f.close()
                print("NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...")
                return
            if not filter_mode == requested_filter_mode:
                if change_index_direction == "+":
                    index = index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    index = index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1
                continue
            if filter_mode == requested_filter_mode and requested_ms_ms_mass == "":
                break
            if (filter_mode == requested_filter_mode) and (round(requested_ms_ms_mass, 2) in ms_ms_masses):
                if not save_folder == "":
                    log_f = open(log_f_filename, "a")
                    log_f.write("INFO:\t" + "FOUND MSMSSPECTRUM" + "\n")
                    log_f.close()
                print("FOUND MSMSSPECTRUM")
                break
            else:
                if change_index_direction == "+":
                    index = index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    index = index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1

    masses = list(f[index]["m/z array"])
    intensities = list(f[index]["intensity array"])
    filter = f[index]["scanList"]["scan"][0]["filter string"]
    filter_mode = ""
    if (" d " in filter) and ("@hcd" in filter):
        filter_mode = "MS/MS"
    if (not " d " in filter) and (not "hcd" in filter):
        filter_mode = "Full scan"
    if (not " d " in filter) and ("hcd" in filter):
        filter_mode = "AIF"
    print(filter_mode)
    print("MS level of the mass spectrum: " + str(f[index]["ms level"]))
    print(filter)
    ms_ms_masses = []
    filter_parsed = filter.split(" ")
    filter_parsed = [x for x in filter_parsed if "hcd" in x]
    for element in filter_parsed:
        ms_ms_masses.append((float(element.split("@")[0])))
    print(ms_ms_masses)
    print(index)
    print(rt_list[index])
    print(f[index]["scanList"]["scan"][0]["scan time"])
    if not save_folder == "":
        fig, ax = plt.subplots()
        ax.bar(masses, intensities, label="measured ions")
        try:
            ax.get_legend().remove()
        except:
            pass
        ax.set_xlabel("masses / Da")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("Mode:" + str(filter_mode) + "; Index: " + str(index) + "; RT: " + str(
            round(rt_list[index], 2)) + ";\nFilter: " + str(filter) + ";\nMS/MS masses: " + str(ms_ms_masses))
        matplotlib.rcParams.update({'figure.autolayout': True})
        image_filepath = save_folder + "Mass_spectrum_index" + str(index) + "_" + str(
            filter_mode.replace("/", "")) + "_" + str(ms_ms_masses) + ".png"
        plt.savefig(image_filepath, bbox_inches='tight')
        metadata = PIL.PngImagePlugin.PngInfo()
        metadata.add_text("masses", str(masses))
        metadata.add_text("intensities", str(intensities))
        metadata.add_text("index", str(index))
        metadata.add_text("orig_file_entry_for_index", str(f[index]))
        metadata.add_text("mode", str(filter_mode))
        metadata.add_text("filter", str(filter))
        metadata.add_text("msms_masses", str(ms_ms_masses))
        target_image = PIL.Image.open(image_filepath)
        target_image.save(image_filepath, pnginfo=metadata)

        old_log_f = open(log_f_filename, "r")
        lines = old_log_f.readlines()
        old_log_f.close()
        pathlib.Path(log_f_filename).unlink()
        log_f_filename = save_folder + "Mass_spectrum_requestedValue" + str(value) + "_trueIndex" + str(
            index) + "_" + str(filter_mode.replace("/", "")) + "_AdditionalInfo" + ".txt"
        log_f = open(log_f_filename, "w")
        for entry in lines:
            log_f.write(entry)
        log_f.write("===================================================================================\n")
        log_f.write("Deviation_to_initial_index\t" + str(abs(index - start_index)) + "\n")
        log_f.write("Index\t" + str(index) + "\n")
        log_f.write("Filter mode\t" + str(filter_mode) + "\n")
        log_f.write("Filter\t" + str(filter) + "\n")
        log_f.write("MS_level\t" + str(f[index]["ms level"]) + "\n")
        log_f.write("MS/MS_masses\t" + str(ms_ms_masses) + "\n")
        log_f.write("===================================================================================\n")
        log_f.write("Masses\t" + str(masses) + "\n")
        log_f.write("Intensities\t" + str(intensities) + "\n")
        log_f.write("Header_of_spectrum\t" + str(f[index]) + "\n")
        log_f.close()

    return [masses, intensities, index, f[index], f[index]["ms level"], filter, filter_mode, ms_ms_masses]


def get_formula_from_cache(cache_folder, mass, max_ppm_dev, file_basename="formulas_for_mass_", filetype=".txt", max_mass=219):
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
                # print(traceback.format_exc())
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
                with open(file, "r") as f:
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
                #print(traceback.format_exc())
                print(all_dict)
                continue
        all_dict = {f: dev for f, dev in all_dict.items() if abs(dev) <= max_ppm_dev}
        all_dict = {k: v for k, v in sorted(all_dict.items(), key=lambda item: abs(item[1]), reverse=False)}

    formula_dict_list = [get_formula_to_dict(f) for f, dev in all_dict.items()]
    formula_dict_list = [dict(t) for t in {tuple(d.items()) for d in formula_dict_list}]
    formula_dict_list = [get_formula_string_from_dict(f) for f in formula_dict_list]
    all_dict = {f: dev for f, dev in all_dict.items() if f in formula_dict_list}
    for formula in list(all_dict.keys()):
        curr_formula_dict = get_formula_to_dict(formula)
        try:
            if ( ("H" in list(curr_formula_dict.keys())) and ("C" in list(curr_formula_dict.keys())) and (not "N" in list(curr_formula_dict.keys())) and (not "O" in list(curr_formula_dict.keys())) ):
                if curr_formula_dict["H"] > curr_formula_dict["C"] * 3:
                    del all_dict[formula]
                    continue
            if ( ("H" in list(curr_formula_dict.keys())) and ("C" in list(curr_formula_dict.keys())) and ("N" in list(curr_formula_dict.keys())) and (not "O" in list(curr_formula_dict.keys())) ):
                if curr_formula_dict["H"] > curr_formula_dict["C"] * 3 + curr_formula_dict["N"]*2:
                    del all_dict[formula]
                    continue
            if ( ("H" in list(curr_formula_dict.keys())) and ("C" in list(curr_formula_dict.keys())) and ("N" in list(curr_formula_dict.keys())) and ("O" in list(curr_formula_dict.keys())) ):
                if curr_formula_dict["H"] > curr_formula_dict["C"] * 3 + curr_formula_dict["N"]*2 + curr_formula_dict["O"]:
                    del all_dict[formula]
                    continue
            if ( ("H" in list(curr_formula_dict.keys())) and ("C" in list(curr_formula_dict.keys())) and (not "N" in list(curr_formula_dict.keys())) and ("O" in list(curr_formula_dict.keys())) ):
                if curr_formula_dict["H"] > curr_formula_dict["C"] * 3 + curr_formula_dict["O"]:
                    del all_dict[formula]
                    continue
            if ("C" in list(curr_formula_dict.keys()) and "O" in list(curr_formula_dict.keys()) and (not "N" in list(curr_formula_dict.keys())) and (not "S" in list(curr_formula_dict.keys())) ):
                if curr_formula_dict["C"]*4 < (curr_formula_dict["O"]):
                    del all_dict[formula]
                    continue
            if ( ("C" in list(curr_formula_dict.keys())) and ("O" in list(curr_formula_dict.keys())) and ("S" in list(curr_formula_dict.keys())) and (not "N" in list(curr_formula_dict.keys()))):
                if curr_formula_dict["C"]*3 + curr_formula_dict["S"]*4 < (curr_formula_dict["O"]):
                    del all_dict[formula]
                    continue
            if (("C" in list(curr_formula_dict.keys())) and ("O" in list(curr_formula_dict.keys())) and ("S" in list(curr_formula_dict.keys())) and ("N" in list(curr_formula_dict.keys()))):
                if curr_formula_dict["C"] * 3 + curr_formula_dict["S"] * 4 + curr_formula_dict["N"]*2 < (curr_formula_dict["O"]):
                    del all_dict[formula]
                    continue
            if ("P" in list(curr_formula_dict.keys())) and (not "O" in list(curr_formula_dict.keys())):
                del all_dict[formula]
                continue
            if ("P" in list(curr_formula_dict.keys()) and "O" in list(curr_formula_dict.keys())):
                if (curr_formula_dict["P"] > curr_formula_dict["O"]/3):
                    del all_dict[formula]
                    continue
        except:
            pass
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


def summarize_mass_intensity_dict_with_deviation(dictio, deviation=20, noise=5000):
    print(len(dictio))
    dictio = {k: v for k, v in dictio.items() if v >= 1}
    old_masses_list = list(dictio.keys())
    old_abundances_list = list(dictio.values())
    new_masses_list = []
    new_abundances_list = []
    for entry in range(len(old_masses_list)):
        try:
            mass_dev_list = [abs(((m - old_masses_list[entry]) / old_masses_list[entry]) * 1000000) for m in old_masses_list]
            curr_mass_plus_deviation_list = [old_masses_list[m] for m in range(len(old_masses_list)) if
                                             mass_dev_list[m] <= deviation]
            curr_abundance_plus_deviation_list = [old_abundances_list[m] for m in range(len(old_masses_list)) if
                                                  mass_dev_list[m] <= deviation]

            if not curr_mass_plus_deviation_list[curr_abundance_plus_deviation_list.index(max(curr_abundance_plus_deviation_list))] in new_masses_list:
                new_masses_list.append(curr_mass_plus_deviation_list[curr_abundance_plus_deviation_list.index( max(curr_abundance_plus_deviation_list))])
                new_abundances_list.append(sum(curr_abundance_plus_deviation_list))

        except Exception as e:
            print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
            print(e)
            break
    outdict = dict(zip(new_masses_list, new_abundances_list))

    new_masses_list = []
    new_abundances_list = []
    for entry in range(len(outdict)):
        try:
            mass_dev_list = [abs(((m - list(outdict.keys())[entry]) / list(outdict.keys())[entry]) * 1000000) for m in old_masses_list]
            curr_mass_plus_deviation_list = [old_masses_list[m] for m in range(len(old_masses_list)) if mass_dev_list[m] <= deviation]
            curr_abundance_plus_deviation_list = [old_abundances_list[m] for m in range(len(old_masses_list)) if mass_dev_list[m] <= deviation]
            cumm_mass = 0
            cumm_abundance = sum(curr_abundance_plus_deviation_list)
            for entry in range(len(curr_mass_plus_deviation_list)):
                cumm_mass = cumm_mass + (curr_mass_plus_deviation_list[entry] * curr_abundance_plus_deviation_list[entry])
            cumm_mass = cumm_mass / cumm_abundance

            new_masses_list.append(cumm_mass)
            new_abundances_list.append(cumm_abundance)

        except Exception as e:
            print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
            print(e)
            break
    outdict = dict(zip(new_masses_list, new_abundances_list))
    return outdict


def simulate_isotope_pattern_of_formula(formula, mass_resolution_ppm=10):
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
    masses_list = []
    abundances_list = []
    for entry in list(abundance_dict.items()):
        masses_list.append(entry[0])
        abundances_list.append(entry[1])
    new_masses_list = []
    new_abundances_list = []
    old_masses_list = copy.deepcopy(masses_list)
    old_abundances_list = copy.deepcopy(abundances_list)
    for entry in range(len(masses_list)):
        try:
            mass_dev_list = [abs(((m-old_masses_list[entry])/old_masses_list[entry])*1000000) for m in old_masses_list]
            curr_mass_plus_deviation_list = [old_masses_list[m] for m in range(len(old_masses_list)) if mass_dev_list[m] <= mass_resolution_ppm]
            curr_abundance_plus_deviation_list = [old_abundances_list[m] for m in range(len(old_masses_list)) if mass_dev_list[m] <= mass_resolution_ppm]

            if not curr_mass_plus_deviation_list[curr_abundance_plus_deviation_list.index(max(curr_abundance_plus_deviation_list))] in new_masses_list:
                new_masses_list.append(curr_mass_plus_deviation_list[curr_abundance_plus_deviation_list.index(max(curr_abundance_plus_deviation_list))])
                new_abundances_list.append(sum(curr_abundance_plus_deviation_list))

            # indices_list = [masses_list.index(m) for m in curr_mass_plus_deviation_list]
            #
            # for i in range(len(masses_list)-1, -1, -1):
            #     if i in indices_list:
            #         del masses_list[i]
            #         del abundances_list[i]
            #     print(masses_list)
            # print(curr_mass_plus_deviation_list)
        except Exception as e:
            print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
            print(e)
            break
    sorted_masses_according_to_abundance = [m for _,m in sorted(zip(new_abundances_list, new_masses_list), reverse=True)]
    sorted_abundance_list = sorted(new_abundances_list, reverse=True)
    abundance_dict = dict(zip(sorted_masses_according_to_abundance, sorted_abundance_list))
    return abundance_dict


def get_formula_to_dict(formula_string):
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
    return "".join([str(a) + str(n) for a, n in formula_dict.items()])


def read_summary_to_df(filepath, sep="\t", colnames=("MASS", "FORMULA", "INDEX", "RT", "INTENSITY", "SCORE", "FRAGS", "F_SCORES", "F_INTENSITIES", "NLS", "NLS_DEV"), header=None):
    import pandas as pd
    df = pd.read_csv(filepath, sep=sep, header=header)
    new_colnames = []
    for entry in range(len(df.columns)):
        try:
            new_colnames.append(colnames[entry])
        except Exception as e:
            print("Error setting columnames: " + str(e))
            new_colnames.append(entry)
    df.columns = new_colnames
    return df


def calc_dbe(formula):
    try:
        if isinstance(formula, str):
            formula = get_formula_to_dict(formula)
        if isinstance(formula, dict):
            pass
        dbe = ((formula.get("C", 0)*2) + 2 + formula.get("N", 0) - formula.get("H", 0) - formula.get("F", 0) - formula.get("Cl", 0) - formula.get("Br", 0) - formula.get("I", 0) ) / 2
        return dbe
    except Exception as e:
        print("Error calculating DBE in MS_functions: " + str(e))
        print(traceback.format_exc())
        return 0



if __name__ == "__main__":


    print(get_formula_from_cache("U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//", 220.0194, 20))



    input()
    dictio = {384.93646240234375: 132925.58, 384.9373779296875: 72943.84, 384.9382629394531: 23663.477,
              384.93914794921875: 1052.3583, 388.4117431640625: 7730.6597, 388.41265869140625: 12910.038,
              388.4135437011719: 15305.983, 388.4144592285156: 13459.421, 388.4153747558594: 8582.9,
              396.9542541503906: 5247.326, 396.9551696777344: 11537.47, 396.95611572265625: 15755.13,
              396.9570617675781: 14744.662, 396.9580078125: 9259.735, 396.96923828125: 10694.229,
              396.9701843261719: 15516.561, 396.97113037109375: 16795.674, 396.9720458984375: 13634.491,
              396.9729919433594: 7864.8403}
    newdictio = summarize_mass_intensity_dict_with_deviation(dictio, 2)
    print(dictio)
    print(newdictio)

    input()

    possible_combinations = coin_change_problem(138, coins=[1, 12, 15, 16])
    print(possible_combinations)


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

    mzml_filename = "PaulineF7_APCI_neg_ALL_EXPERIMENTS.mzML"
    f = read_mzml_file(mzml_filename)
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
    fig, ax = plt.subplots()
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
    fig, ax = plt.subplots()
    ax.bar(spec[0], spec[1])
    plt.show()
#######################################################################################################################

    #returns [mass, {formula_approx1: ppm_deviation1, formula_approx2: ppm_deviation2, .....}]
    formula_approximations = get_one_formula(138.0191, max_dev=0.003, charge_of_measured_mass=-1, dbe_range=(-1, 15), c_oxidation_state_range=(-4, 4),
                    heteroatom_ratio_range=(0, 0.75), aromaticity_eq_range=(-20, 100), h_c_ratio=(0.3, 3),
                    max_h=100, max_c=100, max_n=100, max_o=100, max_s=0, max_p=0, max_cl=0, print_progress=False)
#######################################################################################################################





