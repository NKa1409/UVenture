import copy
import datetime
import os
import pandas as pd
from matplotlib import pyplot as plt
import scipy.fft
import UVenture
import numpy as np
import scipy
import math

import MS_functions
import peakdetection_funcs

def get_rounding_precision(number):
    if number == 0:
        return 0  # Handle the edge case for 0
    # Get the order of magnitude of the number
    magnitude = math.floor(math.log10(abs(number)))
    # Calculate the precision for rounding
    precision = -magnitude
    return precision

def remove_isotopo_signals(peak_df, mass_deviation_isotopo=2, height_deviation_isotopo=0.5):
    # Convert mass_deviation_isotopo from ppm to absolute value
    mass_deviation_ppm = copy.deepcopy(mass_deviation_isotopo)
    # Remove isotopo signals
    isotopo_signal_rows = []
    for i in peak_df["peak_bins"].unique():
        subgroup_df = peak_df[peak_df["peak_bins"] == i]
        subgroup_df = subgroup_df.sort_values(by=["height"], ascending=[True])
        for index, row in subgroup_df.iterrows():
            curr_row = row
            curr_mass = row["mass"]
            curr_height = row["height"]
            #Start checking if an isotopo signal exists for curr_row
            for index2, row2 in subgroup_df.iterrows():
                if index == index2:
                    continue
                mass_deviation_isotopo = (mass_deviation_ppm / 1000000) * curr_mass
                if (curr_mass+1.00336-mass_deviation_isotopo) <= row2["mass"] <= (curr_mass+1.00336+mass_deviation_isotopo):
                    height_ratio = row2["height"] / curr_height
                    theo_height_ratio = (((curr_mass*0.85)*0.7 ) / 12) * 0.01082 # first multiplicator is for hydrogen atoms and second is for heteroatoms
                    if theo_height_ratio * (1 - height_deviation_isotopo) <= height_ratio <= theo_height_ratio * (1 + height_deviation_isotopo):
                        print("Found isotopo signal: " + str(row) + " " + str(row2))
                        isotopo_signal_rows.append(row2)
    print("Isotopo signal rows: " + str(isotopo_signal_rows))
    indices_to_delete = peak_df.index[peak_df.apply(tuple, axis=1).isin([tuple(row) for row in isotopo_signal_rows])]
    indices_to_delete = list(set(indices_to_delete))
    print("Indices to delete: " + str(indices_to_delete))
    peak_df.drop(indices_to_delete, inplace=True)
    print("Peak df after deleting isotopo signals: " + str(peak_df))
    return peak_df

def remove_duplicates_from_peak_df(peak_df, rt_bins=100, mass_deviation_isotopo=0.0007):
    peak_df = peak_df.drop("area", axis=1)
    rt_round_value = max(ms_file.rt_range) / (rt_bins*2)
    rt_round_value = get_rounding_precision(rt_round_value)
    mz_round_value = get_rounding_precision(mass_deviation_isotopo)
    peak_df["unique_descriptor"] = peak_df["mass"].round(mz_round_value).astype(str) + peak_df["rt"].round(rt_round_value).astype(str)
    peak_df.drop_duplicates(subset=["unique_descriptor"], inplace=True)
    peak_df["peak_bins"] = pd.cut(peak_df["rt"], bins=rt_bins, labels=False)
    peak_df = peak_df.sort_values(by=["rt"], ascending=[True])
    return peak_df

def remove_duplicates_from_peak_df_2(peak_df, rt_bins=100, mass_deviation_ppm=15):
    rt_threshold = 2
    def calculate_mass_threshold(mass):
        return mass * (mass_deviation_ppm / 1000000)
    peak_df = peak_df.sort_values(by=["mass", "rt"]).reset_index(drop=True)
    keep_mask = [True] * len(peak_df)

    for i, row_i in peak_df.iterrows():
        if not keep_mask[i]:
            continue
        for j in range(i + 1, len(peak_df)):  # Compare only subsequent rows
            if not keep_mask[j]:
                continue  # Skip rows already marked as duplicates
            row_j = peak_df.iloc[j]
            mass_diff = abs(row_i["mass"] - row_j["mass"])
            rt_diff = abs(row_i["rt"] - row_j["rt"])
            if rt_diff <= rt_threshold and mass_diff <= calculate_mass_threshold(row_i["mass"]):
                keep_mask[j] = False
    print("Keep mask: " + str(keep_mask))
    peak_df = peak_df[keep_mask].reset_index(drop=True)
    
    peak_df["peak_bins"] = pd.cut(peak_df["rt"], bins=rt_bins, labels=False)
    peak_df = peak_df.sort_values(by=["rt"], ascending=[True])
    return peak_df

def get_all_possible_peaks_2(ms_file, mass_range=1, 
                             threshold_area=40000, threshold_intensity=10000,
                             rt_bins=400, mass_deviation_isotopo=2, height_deviation_isotopo=0.5,
                             min_peak_width=4, max_peak_width=40,
                             peaklist_filename="", mzrt_filename=""):
    if not isinstance(mass_range, int):
        print("mass_range must be an integer")
        print("Setting mass_range to 1")
        mass_range = 1
    min_mz, max_mz = ms_file.mz_range
    min_mz = int(min_mz)
    max_mz = int(max_mz)

    peak_df = pd.DataFrame(columns=["mass", "rt", "height", "area"])
    for mass in range(max_mz, min_mz-1, -mass_range):
        xic = MS_functions.get_xic(ms_file.rawdata, mass=mass, mass_deviation=mass_range, requested_filter_mode="Full scan")
        times = np.array(xic[0])
        intensities = np.array(xic[1])
        true_indices_of_entries = np.array(xic[2])
        peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=min_peak_width, max_width_seconds=max_peak_width)
        #Calculate peak area
        peak_areas = []
        for peak in range(len(peaks)):
            left = int(properties["left_ips"][peak])
            right = int(properties["right_ips"][peak])
            area = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
            if area <= threshold_area:
                continue
            spec = UVenture.Spec(ms_file, true_indices_of_entries[peaks[peak]], spec_requested_filter_mode="Full scan", mass_deviation=mass_range)
            interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(mass_range*1.1) and k < mass+(mass_range*1.1)}
            interesting_range = {k: v for k, v in interesting_range.items() if v > threshold_intensity}
            for k, v in interesting_range.items():
                print("Mass: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                peak_areas.append((k, times[peaks[peak]], v, area))
        
        # Add the peak areas to the DataFrame
        if len(peak_areas) > 0:
            new_peak_df = pd.DataFrame(peak_areas, columns=["mass", "rt", "height", "area"])
            peak_df = pd.concat([peak_df, new_peak_df], ignore_index=True)

        if len(peak_df) <=2:
            continue

        # Remove duplicates
        #print("Peak df before removing duplicates: " + str(peak_df))
        peak_df = remove_duplicates_from_peak_df_2(peak_df, rt_bins=rt_bins, mass_deviation_ppm=7)
        print("Peak df after removing duplicates: " + str(peak_df))

        # Remove isotopo signals
        peak_df = remove_isotopo_signals(peak_df, mass_deviation_isotopo=mass_deviation_isotopo, height_deviation_isotopo=height_deviation_isotopo)

        new_peak_df = peak_df.copy()
        for index, row in peak_df.iterrows():
            if not row["mass"] >= mass + 5:
                continue
            new_peak_df = new_peak_df.drop(index)
            print("Mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " Height: " + str(row["height"]))
            xic_mass_deviation = (mass_deviation_isotopo * row["mass"]) / 1000000
            xic_mass_deviation = xic_mass_deviation * 3
            xic = MS_functions.get_xic(ms_file.rawdata, mass=row["mass"], mass_deviation=xic_mass_deviation, requested_filter_mode="Full scan")
            times = np.array(xic[0])
            intensities = np.array(xic[1])
            true_indices_of_entries = np.array(xic[2])
            peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=min_peak_width, max_width_seconds=max_peak_width)
                
            #Check if a peak exists at row["rt"]
            if len(peaks) == 0:
                print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                continue
            else:
                peak_found = False
                for peak in range(len(peaks)):
                    left = int(properties["left_ips"][peak])
                    right = int(properties["right_ips"][peak])
                    if times[left] <= row["rt"] <= times[right]:
                        # A peak was found at the specified RT
                        # Peak will be added to the peaklist
                        print("Peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                        new_row = copy.deepcopy(row)
                        new_row["area"] = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
                        peak_found = True
                        if not peaklist_filename == "":
                            # Save the peak to the peaklist file
                            with open(peaklist_filename, "a") as f:
                                f.write(f"{new_row['mass']}\t{new_row['rt']}\t{new_row['area']}\t{new_row["height"]}\t{times[left]}\t{times[right]}\n")
                        if not mzrt_filename == "":
                            # Save the peak to the mzrt file to be processed directly
                            with open(mzrt_filename, "a") as f:
                                file = os.path.normpath(ms_file.filename)
                                file = file.split(os.sep)[-1]
                                print(file)
                                f.write(f"{file}\t{new_row['mass']}\t{new_row['rt']}\n")
                        break
                if not peak_found:
                    print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    continue
        peak_df = new_peak_df.copy()
    return True


def perform_peak_detection_of_complete_file(ms_file, mass_range=1, 
                                            threshold_area=40000, threshold_intensity=10000, 
                                            rt_bins=400, mass_deviation_isotopo=0.0007, height_deviation_isotopo=0.5,
                                            peaklist_filename="", mzrt_filename=""):
    identified_peaks = get_all_possible_peaks(ms_file, mass_range=mass_range, threshold_area=threshold_area, threshold_intensity=threshold_intensity)
    # identified_peaks = [[mass, height, area, rt], ...]

    # Convert to DataFrame and remove duplicates
    peak_df = pd.DataFrame(identified_peaks, columns=["mass", "rt", "height", "area"])
    # Remove duplicates
    peak_df = remove_duplicates_from_peak_df(peak_df, rt_bins=rt_bins, mass_deviation_isotopo=mass_deviation_isotopo)
    

    # Remove isotopo signals
    peak_df = remove_isotopo_signals(peak_df, mass_deviation_isotopo=mass_deviation_isotopo, height_deviation_isotopo=height_deviation_isotopo)
    

    keep_rows = []
    for index, row in peak_df.iterrows():
        print("Mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " Height: " + str(row["height"]))
        xic = MS_functions.get_xic(ms_file.rawdata, mass=row["mass"], mass_deviation=mass_deviation_isotopo, requested_filter_mode="Full scan")
        times = np.array(xic[0])
        intensities = np.array(xic[1])
        true_indices_of_entries = np.array(xic[2])
        peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=2, max_width_seconds=30)
            
        #Check if a peak exists at row["rt"]
        if len(peaks) == 0:
            print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
            continue
        else:
            peak_found = False
            for peak in range(len(peaks)):
                left = int(properties["left_ips"][peak])
                right = int(properties["right_ips"][peak])
                if times[left] <= row["rt"] <= times[right]:
                    # A peak was found at the specified RT
                    # Peak will be added to the peaklist
                    print("Peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    new_row = copy.deepcopy(row)
                    new_row["area"] = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
                    keep_rows.append(new_row)
                    peak_found = True
                    if not peaklist_filename == "":
                        # Save the peak to the peaklist file
                        with open(peaklist_filename, "a") as f:
                            f.write(f"{new_row['mass']}\t{new_row['rt']}\t{new_row['area']}\t{new_row["height"]}\t{times[left]}\t{times[right]}\n")
                    if not mzrt_filename == "":
                        # Save the peak to the mzrt file to be processed directly
                        with open(mzrt_filename, "a") as f:
                            f.write(f"{new_row['mass']}\t{new_row['rt']}\n")
                    break
            if not peak_found:
                print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                continue
    print("Keep rows: " + str(keep_rows))
    keep_rows = pd.DataFrame(keep_rows)
    print(keep_rows)



if __name__ == "__main__":

    print("Start")
    mzml_filename = "C://Users//Admin//Desktop//UVenture//UVenture-20250516_StartPeakDetection//webserver_save//mzml_files//05Dec2024_PAMPain2_neg_20241008_ISP_OH_O3_NO.mzML"
    ms_file = UVenture.MS_File(mzml_filename)
    min_mz, max_mz = ms_file.mz_range

    mass_range = 1
    threshold_area = 40000
    threshold_intensity = 10000
    identified_peaks = []
    for mass in range(int(725), int(723), -mass_range):
        starttime = datetime.datetime.now()
        xic = MS_functions.get_xic(ms_file.rawdata, mass=mass, mass_deviation=mass_range, requested_filter_mode="Full scan")
        print("XIC done: " + str(datetime.datetime.now() - starttime))
        times = xic[0]
        times = np.array(times)
        intensities = xic[1]
        intensities = np.array(intensities)
        true_indices_of_entries = xic[2]
        true_indices_of_entries = np.array(true_indices_of_entries)
        peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities)
        #Calculate peak area
        peak_areas = []
        for peak in range(len(peaks)):
            left = int(properties["left_ips"][peak])
            right = int(properties["right_ips"][peak])
            print(peak)
            area = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
            if area <= threshold_area:
                continue
            spec = UVenture.Spec(ms_file, true_indices_of_entries[peaks[peak]], spec_requested_filter_mode="Full scan", mass_deviation=mass_range)
            interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(mass_range*1.5) and k < mass+(mass_range*1.5)}
            print("interesting_range: " + str(interesting_range))
            interesting_range = {k: v for k, v in interesting_range.items() if v > threshold_intensity}
            for k, v in interesting_range.items():
                print("Entry: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                peak_areas.append((k, times[peaks[peak]], v, area))
        print(peak_areas)
        identified_peaks.extend(peak_areas)

    peak_df = pd.DataFrame(identified_peaks, columns=["mass", "rt", "height", "area"])
    print(peak_df)
    peak_df = peak_df.drop("area", axis=1)
    peak_df["unique_descriptor"] = peak_df["mass"].round(4).astype(str) + peak_df["rt"].round(5).astype(str)
    peak_df.drop_duplicates(subset=["unique_descriptor"], inplace=True)
    print(peak_df)
    peak_df["peak_bins"] = pd.cut(peak_df["rt"], bins=100, labels=False)
    print(peak_df)
    peak_df = peak_df.sort_values(by=["rt"], ascending=[True])
    mass_deviation_isotopo = 0.0007
    height_deviation_isotopo = 0.5
    isotopo_signal_rows = []
    for i in peak_df["peak_bins"].unique():
        print("Peak bin: " + str(i))
        subgroup_df = peak_df[peak_df["peak_bins"] == i]
        subgroup_df = subgroup_df.sort_values(by=["height"], ascending=[True])
        for index, row in subgroup_df.iterrows():
            print("Mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " Height: " + str(row["height"]))
            curr_row = row
            curr_mass = row["mass"]
            curr_height = row["height"]
            #Start checking if an isotopo signal exists for curr_row
            for index2, row2 in subgroup_df.iterrows():
                if index == index2:
                    continue
                if (curr_mass+1.00336-mass_deviation_isotopo) <= row2["mass"] <= (curr_mass+1.00336+mass_deviation_isotopo):
                    height_ratio = row2["height"] / curr_height
                    theo_height_ratio = (((curr_mass*0.85)*0.7 ) / 12) * 0.01082 # first multiplicator is for hydrogen atoms and second is for heteroatoms
                    if theo_height_ratio * (1 - height_deviation_isotopo) <= height_ratio <= theo_height_ratio * (1 + height_deviation_isotopo):
                        print("Found isotopo signal: " + str(row) + " " + str(row2))
                        isotopo_signal_rows.append(row2)
    print("Isotopo signal rows: " + str(isotopo_signal_rows))
    indices_to_delete = peak_df.index[peak_df.apply(tuple, axis=1).isin([tuple(row) for row in isotopo_signal_rows])]
    indices_to_delete = list(set(indices_to_delete))
    print("Indices to delete: " + str(indices_to_delete))
    peak_df.drop(indices_to_delete, inplace=True)
    print("Peak df after deleting isotopo signals: " + str(peak_df))
    peak_df = peak_df.drop("unique_descriptor", axis=1)

    mass_deviation = 0.001
    keep_rows = []
    counter = 0
    for index, row in peak_df.iterrows():
        if counter == 20:
            break
        counter += 1
        print("Mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " Height: " + str(row["height"]))
        xic = MS_functions.get_xic(ms_file.rawdata, mass=row["mass"], mass_deviation=mass_deviation, requested_filter_mode="Full scan")
        times = xic[0]
        times = np.array(times)
        intensities = xic[1]
        intensities = np.array(intensities)
        true_indices_of_entries = xic[2]
        true_indices_of_entries = np.array(true_indices_of_entries)
        # Do first smoothing
        window_size = peakdetection_funcs.get_window_size_by_frequency(intensities, times, min_width=2, max_width=30)
        intensityvals = peakdetection_funcs.do_smoothing_without_effecting_peaks(intensities, window_size=window_size)
        # Apply baseline correction
        baseline = peakdetection_funcs.baseline_als(intensityvals, lam=1e4, p=0.05, niter=100, window_min_vals=window_size)
        corrected_intensity = intensityvals - baseline
        smoothed_intensity = peakdetection_funcs.do_smoothing_without_effecting_peaks(corrected_intensity, window_size=window_size)
        # Find peaks with adaptive height and width detection
        min_width_seconds = int( (len(times)/(max(times)-min(times)) ) * 1)
        max_width_seconds = int( (len(times)/(max(times)-min(times)) ) * 50)
        min_width_measurements = 2
        max_width_measurements = 30
        min_width = min_width_seconds if min_width_seconds > min_width_measurements else min_width_measurements
        max_width = max_width_seconds if max_width_seconds > max_width_measurements else max_width_measurements
        peaks, properties = scipy.signal.find_peaks(smoothed_intensity,
                                    prominence=np.std(smoothed_intensity), width=(min_width, max_width))
        #Check if a peak exists at row["rt"]
        if len(peaks) == 0:
            print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
            continue
        else:
            peak_found = False
            for peak in range(len(peaks)):
                left = int(properties["left_ips"][peak])
                right = int(properties["right_ips"][peak])
                if times[left] <= row["rt"] <= times[right]:
                    print("Peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    new_row = copy.deepcopy(row)
                    new_row["area"] = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
                    keep_rows.append(new_row)
                    peak_found = True
                    break
            if not peak_found:
                print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                continue
    print("Keep rows: " + str(keep_rows))
    keep_rows = pd.DataFrame(keep_rows)
    print(keep_rows)

    print("End")
    

