import copy
import os
import sys
import pandas as pd
import numpy as np

import UVenture.class_Spec as class_Spec
import UVenture.peakdetection_funcs as peakdetection_funcs
import UVenture.MS_functions as MS_functions


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
    return peak_df

def remove_duplicates_from_peak_df(peak_df, rt_bins=100, mass_deviation_ppm=15):
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

def get_all_possible_peaks(ms_file, settings_dict, mass_range=1, 
                             threshold_area=40000, threshold_intensity=10000,
                             rt_bins=400, mass_deviation_isotopo=2, height_deviation_isotopo=0.5,
                             min_peak_width=4, max_peak_width=40,
                             peaklist_filename="", mzrt_filename=""):
    if not isinstance(mass_range, int):
        print("mass_range must be an integer. Setting mass_range to 1")
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
            spec = class_Spec.Spec(ms_file, true_indices_of_entries[peaks[peak]], spec_requested_filter_mode="Full scan", mass_deviation=mass_range)
            interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(mass_range*1.1) and k < mass+(mass_range*1.1)}
            interesting_range = {k: v for k, v in interesting_range.items() if v > threshold_intensity}
            for k, v in interesting_range.items():
                print("Mass: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                peak_areas.append((k, times[peaks[peak]], v, area))
        
        # Add the peak areas to the DataFrame
        if len(peak_areas) > 0:
            new_peak_df = pd.DataFrame(peak_areas, columns=["mass", "rt", "height", "area"])
            if len(peak_df) >= 1:
                peak_df = pd.concat([peak_df, new_peak_df], ignore_index=True)
            elif len(peak_df) == 0:
                peak_df = new_peak_df

        if len(peak_df) <= 2:
            continue

        # Remove duplicates
        #print("Peak df before removing duplicates: " + str(peak_df))
        peak_df = remove_duplicates_from_peak_df(peak_df, rt_bins=rt_bins, mass_deviation_ppm=7)
        # Remove isotopo signals
        peak_df = remove_isotopo_signals(peak_df, mass_deviation_isotopo=mass_deviation_isotopo, height_deviation_isotopo=height_deviation_isotopo)
        print("Peak df after removing duplicates and isotopo signals: " + str(peak_df))

        

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
                                f.write(f"{new_row['mass']}\t{new_row['rt']}\t{round(new_row['area'], 2)}\t{round(new_row["height"], 2)}\t{times[left]}\t{times[right]}\n")
                            print("Peak saved to peaklist: " + str(peaklist_filename) + ".")
                        if not mzrt_filename == "":
                            # Save the peak to the mzrt file to be processed directly
                            with open(mzrt_filename, "a") as f:
                                file = os.path.normpath(ms_file.filename)
                                file = file.split(os.sep)[-1]
                                print("Taskstorage file: " + str(file))
                                f.write(f"{file}\t{new_row['mass']}\t{new_row['rt']}\t{settings_dict}\n")
                            print("Peak saved to mzrt file: " + str(mzrt_filename) + ".")
                        break
                if not peak_found:
                    print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    continue
        peak_df = new_peak_df.copy()
    return True


if __name__ == "__main__":

    print("This module is not meant to be run directly. Please use it as part of the UVenture package.")
    sys.exit()