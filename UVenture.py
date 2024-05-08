import ast
import copy
import datetime
import math
import pathlib
import sys
import traceback
import PIL
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
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

from plotly.subplots import make_subplots
import plotly.graph_objects as go



class MS_File:
    def __init__(self, filename, **kwargs):
        if "parentfolder_msfile" in kwargs:
            parentfolder = kwargs["parentfolder_msfile"]
        else:
            parentfolder = str(".".join(filename.split(".")[:-1]) + "/")
        default_kwargs = {"parentfolder_msfile":parentfolder,
                          "logfile_filepath":parentfolder + "MSfile_logfile.txt" }
        kwargs = {**default_kwargs, **kwargs}
        os.makedirs(kwargs["parentfolder_msfile"], exist_ok=True)

        self.parentfolder = kwargs["parentfolder_msfile"]
        self.filename = filename
        self.file = mzml.read(self.filename)
        self.rawdata = list(self.file)
        self.method_duration = self.rawdata[-1]["scanList"]["scan"][0]["scan time"]
        self.rt_list = [element["scanList"]["scan"][0]["scan time"] for element in self.rawdata]
        self.all_ms_spectra = []
        self.all_filters = []
        self.all_modes = []
        self.tic = []
        self.available_modes = []

        for index in range(len(self.rawdata)):
            masses = list(self.rawdata[index]["m/z array"])
            intensities = list(self.rawdata[index]["intensity array"])
            mass_intensity_dict = dict(zip(masses, intensities))
            self.all_ms_spectra.append(mass_intensity_dict)
            self.tic.append(self.rawdata[index]["total ion current"])
            self.all_filters.append(self.rawdata[index]["scanList"]["scan"][0]["filter string"])
        for filter_string in self.all_filters:
            if " d " in filter_string and "@hcd" in filter_string:
                self.available_modes.append("MS/MS")
            elif " d " not in filter_string and "hcd" not in filter_string:
                self.available_modes.append("Full scan")
            elif " d " not in filter_string and "hcd" in filter_string:
                self.available_modes.append("AIF")
        self.all_modes = self.available_modes
        self.available_modes = list(set(self.available_modes))
        

class Spec:
    def __init__(self, ms_file, index, **kwargs):
        self.original_index = index
        self.ms_file = ms_file

        if "spec_requested_filter_mode" in kwargs:
            requested_mode = kwargs["spec_requested_filter_mode"]
        else:
            requested_mode = "Full scan"
        if "spec_requested_ms_ms_mass" in kwargs:
            requested_ms_ms_mass = kwargs["spec_requested_ms_ms_mass"]
        else:
            requested_ms_ms_mass = ""
        self.index = self.search_for_required_spec_close_to_rt(index, requested_mode=requested_mode, requested_ms_ms_mass=requested_ms_ms_mass)

        if "spec_subfolder" in kwargs:
            spec_folder = self.ms_file.parentfolder + "/spectra/" + kwargs["spec_subfolder"]
        else:
            spec_folder = str(self.ms_file.parentfolder + "/spectra/" + "massspec_requestIndex_" + str(self.index) + "/")
        if "absolute_spec_folder" in kwargs:
            spec_folder = kwargs["absolute_spec_folder"]
        if "spec_subfolder" in kwargs and "absolute_spec_folder" in kwargs:
            print("ERROR: You can only use either 'spec_subfolder=' or 'absolute_spec_folder=', not both!")
            print("Using absolute_spec_folder now...")
            print(kwargs["absolute_spec_folder"])

        self.rt = self.ms_file.rt_list[self.index]
        self.masses = []
        self.summarized_masses = []
        self.intensities = []
        self.summarized_intensities = []
        self.mass_intensity_dict = {}
        self.summarized_mass_intensity_dict = {}
        self.filter = None
        self.filter_mode = None
        self.ms_ms_masses = None
        self.spec_rawdata = None
        self.ms_level = None
        
        default_kwargs = {"spec_folder":spec_folder,
                          "spec_requested_filter_mode":requested_mode,
                          "spec_requested_ms_ms_mass":requested_ms_ms_mass,
                          "spec_log_filepath":spec_folder + "spectrum_creation_log_" + str(self.index) + ".txt",
                          "mass_deviation": 11,

                          "spec_save_matplotlib_plot":False,
                          "spec_save_go_plot":False }
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["spec_folder"], exist_ok=True)

        

        self.make_spec_log_entry("")
        self.make_spec_log_entry("=============================================================")
        self.make_spec_log_entry("=============================================================")
        self.make_spec_log_entry("=============================================================")
        self.make_spec_log_entry("INFO:\t" + "Creating spectrum object for index: " + str(self.index))
        self.make_spec_log_entry("INFO:\t" + "The original index was: " + str(self.original_index) + " This results in a deviation of " + str(abs(self.original_index - self.index)) + " indices.")
        self.make_spec_log_entry("INFO:\t" + "This translates to a deviation of " + str(round(self.ms_file.rt_list[self.original_index] - self.ms_file.rt_list[self.index], 2)) + " seconds in retention time.")
        self.make_spec_log_entry("INFO:\t" + "Spectrum is at retention time: " + str(self.rt))
        self.make_spec_log_entry("INFO:\t" + "Requested filter mode: " + str(self.kwargs["spec_requested_filter_mode"]))
        self.make_spec_log_entry("INFO:\t" + "Requested MS/MS mass: " + str(self.kwargs["spec_requested_ms_ms_mass"]))
        self.make_spec_log_entry("INFO:\t" + "Saving plot?: " + str(self.kwargs["spec_save_matplotlib_plot"]))

        print("Getting mass spec...")
        self.mass_intensity_dict = self.get_mass_spec(self.index)
        print("Got mass spec")
        print("Summarizing mass intensity dict...")
        self.summarized_mass_intensity_dict = self.summarize_mass_intensity_dict_with_deviation(dict(zip(self.masses, self.intensities)), deviation=self.kwargs["mass_deviation"])
        self.summarized_masses = list(self.summarized_mass_intensity_dict.keys())
        self.summarized_intensities = list(self.summarized_mass_intensity_dict.values())
        print("Summarized mass intensity dict shortened from: " + str(len(self.mass_intensity_dict)) + " to: " + str(len(self.summarized_mass_intensity_dict)) + " elements.")

        if self.kwargs["spec_save_go_plot"] == True:
            self.create_go_plot()
        if self.kwargs["spec_save_matplotlib_plot"] == True:
            self.create_matplotlib_plot()

        self.make_spec_log_entry("INFO:\t" + "Found MS/MS masses: " + str(self.ms_ms_masses))
        self.make_spec_log_entry("INFO:\t" + "Filter mode: " + str(self.filter_mode))
        self.make_spec_log_entry("INFO:\t" + "Filter: " + str(self.filter))
        self.make_spec_log_entry("INFO:\t" + "Length of summarized mass intensity dict: " + str(len(self.summarized_mass_intensity_dict)) + " elements.")
        self.make_spec_log_entry("INFO:\t" + "TIC: " + str(self.ms_file.tic[self.index]))
    
    def min_deviation(self, input_list):
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

    def get_best_approx_for_ppm_spacing_within_peak(self, mass_list, worst_expected_ppm_deviation=20):
        mass_list = sorted(mass_list)
        ppm_spacing_list = [(((mass_list[i+1] - mass_list[i]) / mass_list[i]) * 1000000) for i in range(len(mass_list)-1)]
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
        #sort the dictio by its keys
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
                mass_lower_border = curr_mass - ((deviation*curr_mass)/1000000)
                mass_upper_border = curr_mass + ((deviation*curr_mass)/1000000)
                iteration_step_lower = 0
                integration_step_upper = 0
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                expected_min_spacing_between_measurement_points = (best_approx_ppm_spacing_within_peak * curr_mass) / 1000000
                try:
                    while (mass_lower_border < old_masses_list[index_of_highest_abundance-iteration_step_lower]) and \
                            (last_existing_mass_within_border - old_masses_list[index_of_highest_abundance-iteration_step_lower] <= 2.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance-iteration_step_lower] <= (1.1 * last_abundance)):
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance-iteration_step_lower]
                        last_abundance = old_abundances_list[index_of_highest_abundance-iteration_step_lower]
                        iteration_step_lower += 1
                except IndexError:
                    remove_all_lower = True
                    iteration_step_lower = 0

                
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                try:
                    while (mass_upper_border > old_masses_list[index_of_highest_abundance+integration_step_upper]) and \
                            (old_masses_list[index_of_highest_abundance+integration_step_upper] - last_existing_mass_within_border <= 2.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance+integration_step_upper] <= (1.1 * last_abundance)):
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance+integration_step_upper]
                        last_abundance = old_abundances_list[index_of_highest_abundance+integration_step_upper]
                        integration_step_upper += 1
                except IndexError:
                    remove_all_upper = True
                    integration_step_upper = 0
                
                if remove_all_lower == True and remove_all_upper == True:
                    break
                if remove_all_upper:
                    summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]) / summed_intensity
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_masses_list), 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]
                    continue
                if remove_all_lower:
                    summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper, 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper, 1)]) / summed_intensity
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper, 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper, 1)]
                    continue
                    
                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)]) / summed_intensity
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)]
            except Exception as e:
                print("EXCEPTION IN summarize_mass_intensity_dict_with_deviation()!!!")
                print(traceback.format_exc())
                print(e)
                break
        outdict = dict(zip(new_masses_list, new_abundances_list))
        print("new length of summarized dictio:")
        print(len(outdict))
        print("minimal deviation between elements: " + str(self.min_deviation(new_masses_list)))
        print("End time: " + str(datetime.datetime.now()))
        print("Time taken = " + str(datetime.datetime.now() - start_time))
        return outdict

    def get_mode_of_spec(self, filter_string):
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

    def make_spec_log_entry(self, log_entry):
        #get the dirname of the spec_log_filepath
        directory_logfile = os.path.dirname(self.kwargs["spec_log_filepath"])
        #create the directory if it does not exist
        os.makedirs(directory_logfile, exist_ok=True)
                
        log_f = open(self.kwargs["spec_log_filepath"], "a")
        log_f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t" + log_entry + "\n")
        log_f.close()
        return True

    def create_go_plot(self):
        # Sort the data by intensities
        data = sorted(zip(self.summarized_masses, self.summarized_intensities), key=lambda x: x[1], reverse=True)
        # Select the top 10
        top_10_data = data[:10]
        # Create the figure
        go_fig = go.Figure(data=go.Bar(x=self.summarized_masses, y=self.summarized_intensities, marker=dict(color='black', opacity=1), width=0.1))
        # Add data labels for the top 10 values
        for mass, intensity in top_10_data:
            go_fig.add_annotation(x=mass, y=intensity, text=str(mass), showarrow=False, font=dict(size=12, color="Black"), bgcolor="White", opacity=0.8, textangle=-90)
        go_fig.update_layout(plot_bgcolor='white')
        go_fig.write_html(self.kwargs["spec_folder"] + "Mass_spectrum_index" + str(self.index) + "_" + str(self.filter_mode.replace("/", "")) + "_" + str(self.ms_ms_masses) + ".html")

    def create_matplotlib_plot(self):
        self.make_spec_log_entry("INFO:\t" + "Saving mass spectrum plot...")
        fig, ax = plt.subplots()
        ax.bar(self.summarized_masses, self.summarized_intensities, label="measured ions", width=0.2, color="blue")
        # Sort the data by intensities
        data = sorted(zip(self.summarized_masses, self.summarized_intensities), key=lambda x: x[1], reverse=True)
        # Select the top 4
        top_10_data = data[:4]
        # Add labels for the top 4 values
        for mass, intensity in top_10_data:
            ax.text(mass, (intensity), str(round(mass, 4)), ha='center', va='bottom', rotation=0)
        try:
            ax.get_legend().remove()
        except:
            pass
        ax.set_xlabel("masses / Da")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("Mode:" + str(self.filter_mode) + "; Index: " + str(self.index) + "; RT: " + str(
            round(self.ms_file.rt_list[self.index], 2)) + ";\nFilter: " + str(self.filter) + ";\nMS/MS masses: " + str(self.ms_ms_masses))
        matplotlib.rcParams.update({'figure.autolayout': True})
        image_filepath = self.kwargs["spec_folder"] + "Mass_spectrum_index" + str(self.index) + "_" + str(self.filter_mode.replace("/", "")) + "_" + str(self.ms_ms_masses) + ".png"
        plt.savefig(image_filepath, bbox_inches='tight', dpi=1000)
        plt.close(fig)
        metadata = PIL.PngImagePlugin.PngInfo()
        metadata.add_text("masses", str(self.summarized_masses))
        metadata.add_text("intensities", str(self.summarized_intensities))
        metadata.add_text("index", str(self.index))
        metadata.add_text("orig_file_entry_for_index", str(self.ms_file.rawdata[self.index]))
        metadata.add_text("mode", str(self.filter_mode))
        metadata.add_text("filter", str(self.filter))
        metadata.add_text("msms_masses", str(self.ms_ms_masses))
        target_image = PIL.Image.open(image_filepath)
        target_image.save(image_filepath, pnginfo=metadata)
        print("Saved image")
        self.make_spec_log_entry("INFO:\t" + "Finished saving mass spectrum plot...")

    def search_for_required_spec_close_to_rt(self, index, requested_mode="whatever", requested_ms_ms_mass=""):
        curr_index = copy.deepcopy(index)
        change_index_value = 1
        change_index_direction = "+"
        start_index = copy.deepcopy(index)
        print("Starting ms spectrum search at index: " + str(curr_index))

        while True and not requested_mode == "whatever":
            filter = self.ms_file.rawdata[curr_index]["scanList"]["scan"][0]["filter string"]
            filter_mode = self.get_mode_of_spec(filter)
            if filter_mode == "MS/MS":
                ms_ms_masses = []
                filter_parsed = filter.split(" ")
                filter_parsed = [x for x in filter_parsed if "hcd" in x]
                for element in filter_parsed:
                    ms_ms_masses.append(round(float(element.split("@")[0]), 1))

            if (curr_index >= len(self.ms_file.rawdata)-3) or (curr_index <= 3):
                print("NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...")
                print("Returning the start index: " + str(start_index))
                print("Mode of the start index: " + str(self.get_mode_of_spec(self.ms_file.rawdata[start_index]["scanList"]["scan"][0]["filter string"])))
                return start_index
            
            if not filter_mode == requested_mode:
                if change_index_direction == "+":
                    curr_index = curr_index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    curr_index = curr_index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1
                continue

            elif filter_mode == requested_mode and requested_ms_ms_mass == "":
                print("Found requested spectrum")
                break

            elif (filter_mode == requested_mode) and (round(requested_ms_ms_mass, 1) in ms_ms_masses):
                print("FOUND MSMSSPECTRUM")
                break

            else:
                if change_index_direction == "+":
                    curr_index = curr_index + change_index_value
                    change_index_direction = "-"
                    change_index_value = change_index_value + 1
                else:
                    curr_index = curr_index - change_index_value
                    change_index_direction = "+"
                    change_index_value = change_index_value + 1
        return curr_index

    def get_mass_spec(self, index):
        masses = list(self.ms_file.rawdata[index]["m/z array"])
        intensities = list(self.ms_file.rawdata[index]["intensity array"])
        filter = self.ms_file.rawdata[index]["scanList"]["scan"][0]["filter string"]
        filter_mode = self.get_mode_of_spec(filter)
        print(filter_mode)
        print("MS level of the mass spectrum: " + str(self.ms_file.rawdata[index]["ms level"]))
        print(filter)
        ms_ms_masses = []
        filter_parsed = filter.split(" ")
        filter_parsed = [x for x in filter_parsed if "hcd" in x]
        for element in filter_parsed:
            ms_ms_masses.append((float(element.split("@")[0])))

        self.masses = masses
        self.intensities = intensities
        self.filter = filter
        self.filter_mode = filter_mode
        self.ms_ms_masses = ms_ms_masses
        self.spec_rawdata = self.ms_file.rawdata[index]
        self.ms_level = self.ms_file.rawdata[index]["ms level"]
        return dict(zip(masses, intensities))


class Prediction:
    def __init__(self, ms_file, mass, spec, **kwargs):
        self.ms_file = ms_file
        self.mass = min(list(spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(mass - x))
        self.spec = spec
        self.rt = self.spec.rt
        self.formula_score_dict = {}
        self.best_formula_prediction = None
        self.score_of_best_formula = None
        self.intensity_of_ion = None
        self.identified_peaks_for_mass = None
        self.peak_found = None # True or False

        if "prediction_subfolder" in kwargs:
            pred_folder = self.ms_file.parentfolder + "/predictions/" + kwargs["prediction_subfolder"]
        else:
            pred_folder = str(self.ms_file.parentfolder + "/predictions/" + str(self.spec.index) + "_" + str(round(self.mass, 4)) + "/")
        if "absolute_pred_folder" in kwargs:
            pred_folder = kwargs["absolute_pred_folder"]
        if "prediction_subfolder" in kwargs and "absolute_pred_folder" in kwargs:
            print("ERROR: You can only use either 'prediction_subfolder=' or 'absolute_pred_folder=', not both!")
            print("Using absolute_pred_folder now...")
            print(kwargs["absolute_pred_folder"])
        default_kwargs = {"pred_folder":pred_folder,  
                          "pred_log_filepath": pred_folder + "prediction_log_for_mass_" + str(round(self.mass, 4)) + ".txt",
                          "pred_save_matplotlib_plot_of_isotopologues":False,
                          "pred_save_go_plot_of_isotopologues":False,
                          "pred_save_xic_plot":False,
                          "pred_save_detailed_log":True,
                          "pred_formula_cache_folder_path":"U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//",

                          "mass_deviation":11,
                          "charge_of_measured_mass":-1,

                          "pred_max_ppm_deviation_change_for_isotopologue":2,
                          "pred_minimum_assumed_noise":10000,
                          "pred_noise_divisor_for_isotopologue_calculation":5,
                          "pred_add_value_to_score_if_isotopo_was_found":250,
                          "pred_isotopologue_score_e_function_exponent":0.4,
                          "pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score":30,
                          "pred_stop_isotopologue_search_if_score_lower_than":-80,
                          "pred_ppm_deviation_score_multiplier":14,
                          "pred_include_likelyhood_of_formula":True,
                          "pred_reject_formula_if_score_lower_than":10,
                          "pred_return_if_no_peak_is_found": True }
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["pred_folder"], exist_ok=True)

        self.intensity_of_ion = sum([self.spec.intensities[i] for i in range(len(self.spec.masses)) if abs(((self.spec.masses[i] - mass)/mass)*1000000) <= self.kwargs["mass_deviation"]])
        self.make_op_log_entry("")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("=============================================================")   
        self.make_op_log_entry("INFO:\t" + "Creating prediction object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.spec.index))    
        self.make_op_log_entry("INFO:\t" + "Intensity of the molecular ion: " + str(self.intensity_of_ion))
        self.make_op_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))

        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, ((self.kwargs["mass_deviation"]*self.mass)/1000000), requested_filter_mode=self.spec.filter_mode)
        self.identified_peaks_for_mass = self.get_peak_properties(self.xic[0], self.xic[1])
        self.prediction_spec_within_peak_range = self.check_if_provided_spectrum_within_peak_range(self.spec.index, self.identified_peaks_for_mass)
        if not self.prediction_spec_within_peak_range:
            self.peak_found = False
            self.make_op_log_entry("WARNING:\t" + "No peak is detected where a prediction should be made!")
            self.make_op_log_entry("WARNING:\t" + "RT: " + str(self.spec.rt) + " Index: " + str(self.spec.index) + " Mass: " + str(self.mass))
            print("No peak is detected where a prediction should be made!")
            print("Looking at peak: RT: " + str(self.spec.rt) + " Index: " + str(self.spec.index) + " Mass: " + str(self.mass))
            if self.kwargs["pred_return_if_no_peak_is_found"]:
                print("Returning...")
                self.make_op_log_entry("Returning...")
                return
            print("However, the program will continue...")
            self.make_op_log_entry("However, the program will continue...")
        else:
            self.peak_found = True
            self.make_op_log_entry("INFO:\t" + "Peak is detected where a prediction should be made!")

        _, possible_formulas = MS_functions.get_formula_from_cache(self.kwargs["pred_formula_cache_folder_path"], self.mass, self.kwargs["mass_deviation"])
        if self.kwargs["charge_of_measured_mass"] < 0:
            self.make_op_log_entry("INFO:\t" + "Charge of the measured mass is negative. Removing formulas with Li, Na and K...")
            possible_formulas = {key: value for key, value in possible_formulas.items() if "Na" not in key and "K" not in key and "Li" not in key}
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("OUTPUT OF SIMPLE FORMULA PREDICTION (ONLY MASS DEVIAITON):")
        for formula, deviation in possible_formulas.items():
            self.make_op_log_entry("{:>20} \t {:>20}".format(str(formula), str(round(deviation, 2))))
        
        molecule_formulas = [[MS_functions.get_formula_to_dict(formula), deviation] for formula, deviation in possible_formulas.items()]
        self.mass_possible_formulas = molecule_formulas


        self.formula_score_dict = {}
        for entry in molecule_formulas:
            self.make_op_log_entry("=============================================================")
            isotope_check = self.check_for_isotope_pattern(entry[0])
            formula_score = sum(e_s[3] for e_s in list(isotope_check.values())) - abs(self.kwargs["pred_ppm_deviation_score_multiplier"] * entry[1])
            self.make_op_log_entry("INFO:\t" + "Formula: " + str("".join([str(a) + str(n) for a, n in entry[0].items()])) + " \t Score: " + str(round(formula_score, 2)))
            self.formula_score_dict[str(entry[0])] = formula_score
            if self.kwargs["pred_save_detailed_log"] == True:
                self.make_op_log_entry("OUTPUT OF FORMULA PREDICTION WITH MASS SPECTRUM (ISOTOPOLOGUES)")
                self.make_op_log_entry("Mass of ion: " + str(self.mass))
                self.make_op_log_entry("PPM DEVIATION OF MOLECULE FORMULA: " + str(round(entry[1], 2)))
                self.make_op_log_entry("{:>25} \t {:>25} \t {:>25} \t {:>25} \t {:>25}".format("Mass_of_isotopologue", "Isotopo_found?", "measured_intensity", "required_intensity", "score_of_isotopologue"))
                for item in list(isotope_check.items()):
                    try:
                        self.make_op_log_entry("{:>25} \t {:>25} \t {:>25} \t {:>25} \t {:>25}".format(str(round(item[0], 4)), str(item[1][0]), str(round(item[1][2], 4)), str(round(item[1][1], 4)), str(round(item[1][3], 4))))
                    except Exception as e:
                        continue
        
        if self.kwargs["pred_include_likelyhood_of_formula"] == True:
            self.formula_score_dict = self.add_likelyhood_of_formula_to_score(self.formula_score_dict)

        self.formula_score_dict = {f: s for f, s in self.formula_score_dict.items() if s > self.kwargs["pred_reject_formula_if_score_lower_than"]}
        self.formula_score_dict = dict(sorted(self.formula_score_dict.items(), key=lambda x: x[1], reverse=True))

        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICITIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>30} \t {:>20}".format(str(item[0]), str(round(item[1], 2))))

        if self.kwargs["pred_save_xic_plot"] == True:
            self.xic_plot_filepath = self.kwargs["pred_folder"] + "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["mass_deviation"]) + "_" + str( self.spec.filter_mode) + ".png"
            self.save_xic_plot(self.xic[0], self.xic[1], self.xic_plot_filepath)
        if self.kwargs["pred_save_matplotlib_plot_of_isotopologues"] == True and len(list(self.formula_score_dict.keys())) > 0:
            self.matplotlib_plot_filepath = self.kwargs["pred_folder"] + "isotopo_matplotlib_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.formula_score_dict.keys())[0]).items()])) + ".png"
            self.make_plot_of_isotopologues(ast.literal_eval(list(self.formula_score_dict.keys())[0]), self.matplotlib_plot_filepath)
        if self.kwargs["pred_save_go_plot_of_isotopologues"] == True and len(list(self.formula_score_dict.keys())) > 0:
            self.go_plot_filepath = self.kwargs["pred_folder"] + "isotopo_go_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.formula_score_dict.keys())[0]).items()])) + ".html"
            self.make_go_plot_of_isotopologues(ast.literal_eval(list(self.formula_score_dict.keys())[0]), self.go_plot_filepath)

    def check_if_provided_spectrum_within_peak_range(self, index, peak_properties):
        for peak in peak_properties:
            if (self.ms_file.rt_list[index] >= peak[2]) and (self.ms_file.rt_list[index] <= peak[3]):
                return True
        return False

    def add_likelyhood_of_formula_to_score(self, formula_score_dict):
        #formula score dict given in the form of {"{'C': 2, 'H': 4, 'O': 1}": 100, "{'C': 3, 'H': 6, 'O': 1}": 200, ...}
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("INFO:\t" + "Adding likelyhood of formula to the score...")
        self.make_op_log_entry("INFO:\t" + "Substract points for N > 2 and C/N <= 4")
        self.make_op_log_entry("INFO:\t" + "Substract points for H/C >= 2")
        self.make_op_log_entry("INFO:\t" + "Substract points for dbe < 0")
        self.make_op_log_entry("INFO:\t" + "Substract points for dbe - O > 7")
        for formula_dict_str, score in formula_score_dict.items():
            try:
                f_dict = ast.literal_eval(formula_dict_str)
                if ("C" in list(f_dict.keys())) and ("N" in list(f_dict.keys())):
                    if int(f_dict["N"]) > 2 and (int(f_dict["C"]) / int(f_dict["N"]) <= 4):
                        formula_score_dict[formula_dict_str] = float(score) - ((float(f_dict["N"]) - 2) ** 2) * 50
                if ("C" in list(f_dict.keys())) and ("H" in list(f_dict.keys())):
                    if (float(f_dict["H"]) / float(f_dict["C"]) >= 2):
                        formula_score_dict[formula_dict_str] = float(score) - (((float(f_dict["H"]) / float(f_dict["C"])) - 2) ** 3) * 50
                dbe = MS_functions.calc_dbe(f_dict)
                if dbe < 0:
                    score_substract = (abs(dbe + 2) * 50) ** 3
                elif (dbe - f_dict.get("O", 0)) > 7:
                    score_substract = abs(dbe - f_dict.get("O", 0) - 7) * 50
                else:
                    score_substract = 0
                formula_score_dict[formula_dict_str] = float(score) - score_substract
            except Exception as e:
                print("ERROR in formula score likelyhood: " + str(e))    
        self.make_op_log_entry("Finished adding likelyhood of formula to the score...")
        return formula_score_dict

    def get_peak_properties(self, times, intensities):
        xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, ((self.kwargs["mass_deviation"]*self.mass)/1000000), requested_filter_mode=self.spec.filter_mode)
        window, order = 5, 3
        print(len(xic[1]))
        intensities = scipy.signal.savgol_filter(intensities, window, order, mode="nearest")
        print(len(intensities))
        peak_properties = scipy.signal.find_peaks(intensities, height=max(intensities)/100, distance=2, prominence=max(intensities)/50, width=(2, 20))
        identified_peaks = []
        for element in range(len(peak_properties[1]["peak_heights"])):
            one_peak = []
            one_peak.append(times[int(peak_properties[1]["left_ips"][element] + (peak_properties[1]["widths"][element]/2))])
            one_peak.append(peak_properties[1]["peak_heights"][element])
            one_peak.append(times[int(peak_properties[1]["left_ips"][element])])
            one_peak.append(times[int(peak_properties[1]["right_ips"][element])])
            one_peak.append(peak_properties[1]["prominences"][element])
            identified_peaks.append(one_peak)
        self.identified_peaks_for_mass = identified_peaks
        print(self.identified_peaks_for_mass)
        #[[time, height, lefttime, righttime, prominence], [time, height, left, right, prominence], ...]
        return identified_peaks

    def make_go_plot_of_isotopologues(self, formula_to_simulate, filepath, include_actual_spectral_data=True):
        isotope_simulation = MS_functions.simulate_isotope_pattern_of_formula(formula_to_simulate)
        fig = make_subplots(rows=1, cols=len(isotope_simulation), shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0, vertical_spacing=0)
        intensities_to_include = []
        for entry in list(isotope_simulation.keys()):
            summed_intensity = sum([self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if abs(((self.spec.summarized_masses[i] - entry) / entry)*1000000) <= (self.kwargs["mass_deviation"])])
            intensities_to_include.append(summed_intensity)
        intensities_to_include = [(i/sum(intensities_to_include)) for i in intensities_to_include]

        for entry in range(len(intensities_to_include)):
            curr_x_values = [self.spec.summarized_masses[i] for i in range(len(self.spec.summarized_masses)) if list(isotope_simulation.keys())[entry]-0.1 <= self.spec.summarized_masses[i] <= list(isotope_simulation.keys())[entry]+0.1] 
            curr_y_values = [self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if list(isotope_simulation.keys())[entry]-0.1 <= self.spec.summarized_masses[i] <= list(isotope_simulation.keys())[entry]+0.1]  
            curr_y_values = [(i/sum(curr_y_values)) for i in curr_y_values]
            if include_actual_spectral_data:
                fig.add_trace(go.Bar(x=curr_x_values, y=curr_y_values, name="all ions", marker=dict(color='green', opacity=1), width=0.001), row=1, col=entry+1)
            fig.add_trace(go.Bar(x=[list(isotope_simulation.keys())[entry]], y=[intensities_to_include[entry]], name="measured ions", marker=dict(color='blue', opacity=0.65), width=0.01), row=1, col=entry+1)
            fig.add_trace(go.Bar(x=[list(isotope_simulation.keys())[entry]], y=[-1 * list(isotope_simulation.values())[entry]], name="simulated intensity", marker=dict(color='red', opacity=0.65), width=0.01), row=1, col=entry+1)
            fig.update_xaxes(range=[(list(isotope_simulation.keys())[entry]) - 0.1, (list(isotope_simulation.keys())[entry]) + 0.1], row=1, col=entry+1)
            fig.add_shape(type="line", x0=0, x1=400, y0=-0.01, y1=0.01, line=dict(color="black", width=10), row=1, col=entry+1)
        fig.update_yaxes(range=[-1, 1], showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor="Gray", )
        fig.update_layout(shapes=[dict(type="rect", xref="paper", yref="paper", x0=0, y0=0, x1=1, y1=1, line=dict(color="Black", width=4))],
                          barmode='overlay', title_text="Isotopologues Plot", xaxis_title="masses / Da", yaxis_title="intensity / a.u.", plot_bgcolor='white')
        fig.write_html(filepath)
        return True
    
    def check_for_isotope_pattern(self, formula_dict):
        ordered_mass_intensity_dict = {self.spec.masses[i]: self.spec.intensities[i] for i in range(len(self.spec.masses))}
        ordered_mass_intensity_dict = dict(sorted(ordered_mass_intensity_dict.items(), key=lambda x: x[1], reverse=True))
        ordered_mass_intensity_dict = {key: value for key, value in ordered_mass_intensity_dict.items() if value > 1}
        
        isotope_pattern_dict = MS_functions.simulate_isotope_pattern_of_formula(formula_dict)
        isotope_pattern_dict = dict(sorted(isotope_pattern_dict.items(), key=lambda x: x[1], reverse=True))
        isotope_pattern_dict = {abs((m - 0.000548 * self.kwargs["charge_of_measured_mass"]) / self.kwargs["charge_of_measured_mass"]): a for m, a in isotope_pattern_dict.items()}
        self.make_op_log_entry("Simulated isotope pattern (already corrected for electron mass): " + str(isotope_pattern_dict))

        isotopes_found = {}
        previous_deviation_list = []
        startitem = list(isotope_pattern_dict.items())[0]
        isotopo_mass = startitem[0]
        measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(isotopo_mass - x))
        initial_deviation = (abs(isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000


        for isotopo in list(isotope_pattern_dict.items()):
            isotopo_mass = isotopo[0]
            within_deviation_mass_list = [masse for masse in list(self.spec.summarized_mass_intensity_dict.keys()) if abs((((isotopo_mass - masse) / isotopo_mass) * 1000000) - initial_deviation) < self.kwargs["pred_max_ppm_deviation_change_for_isotopologue"]]
            measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
            deviation = (abs(isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000
            theoretical_intensity = (self.intensity_of_ion / (list(isotope_pattern_dict.items())[0][1])) * isotopo[1]
            if theoretical_intensity <= 0.01:
                theoretical_intensity = 0.01
            noise = sorted(list(self.spec.summarized_mass_intensity_dict.values()))[int(len(list(self.spec.summarized_mass_intensity_dict.items())) / self.kwargs["pred_noise_divisor_for_isotopologue_calculation"])] + self.kwargs["pred_minimum_assumed_noise"]
            measured_intensity = 0

            if (deviation < self.kwargs["mass_deviation"]) and (abs(deviation - initial_deviation) <= self.kwargs["pred_max_ppm_deviation_change_for_isotopologue"]):
                previous_deviation_list.append(copy.deepcopy(deviation))
                for masse in within_deviation_mass_list:
                    measured_intensity = measured_intensity + self.spec.summarized_mass_intensity_dict[masse]
                if measured_intensity <= 1:
                    break
                isotopes_found[isotopo[0]] = [True, theoretical_intensity, measured_intensity]
                score = measured_intensity / theoretical_intensity  # je naeher an 1 desto besser; wenn <1: weniger gemessen als da sein sollte; wenn >1: mehr gemessen als da sein sollte.
                if score > 1:
                    score = -(-1.5 + (1 / (0.7 + (2.718281828459045 ** (-self.kwargs["pred_isotopologue_score_e_function_exponent"] * score)))))
                # score: je naeher an 1 desto besser; wenn negativ: weniger gemessen als theoretisch da; wenn positiv: mehr gemessen als theoretisch da.
                score = (abs(score * 100) + 0.001) + self.kwargs["pred_add_value_to_score_if_isotopo_was_found"]
                if measured_intensity < noise:
                    score_multiplier_by_intensity_and_noise = abs(1 / (math.log(measured_intensity / noise) - 1))
                else:
                    score_multiplier_by_intensity_and_noise = abs(math.log(measured_intensity / noise) + 1)

                if isotopo == list(isotope_pattern_dict.items())[0]:
                    isotopes_found[isotopo[0]].append(0)
                    previous_score = 9999999
                    continue

                if theoretical_intensity <= (noise / 1.5):
                    score = (score * score_multiplier_by_intensity_and_noise) - abs(deviation * self.kwargs["pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    isotopes_found[isotopo[0]].append(score)
                    if (score - abs(self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
                        isotopes_found[isotopo[0]][-1] = isotopes_found[isotopo[0]][-1] * 0.1
                    break

                else:
                    score = (score * score_multiplier_by_intensity_and_noise) - abs(deviation * self.kwargs["pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    isotopes_found[isotopo[0]].append(score)
                
                if (score < self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]):
                    break

                if (score - abs(self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
                    isotopes_found[isotopo[0]][-1] = isotopes_found[isotopo[0]][-1] * 0.1
                    break
                previous_score = score
            else:
                isotopes_found[isotopo[0]] = [False, theoretical_intensity, measured_intensity]
                negative_score = theoretical_intensity / noise
                negative_score = -10 * negative_score
                isotopes_found[isotopo[0]].append(negative_score)
                break
        return isotopes_found
    
    def make_op_log_entry(self, log_entry, error=False):
        if self.kwargs["pred_save_detailed_log"] == False and error == False:
            return False
        #get the dirname of the spec_log_filepath
        directory_logfile = os.path.dirname(self.kwargs["pred_log_filepath"])
        #create the directory if it does not exist
        os.makedirs(directory_logfile, exist_ok=True)
        log_f = open(self.kwargs["pred_log_filepath"], "a")
        log_f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t" + log_entry + "\n")
        log_f.close()
        return True

    def make_plot_of_isotopologues(self, formula_to_simulate, filepath, include_actual_spectral_data=True):
        print("Making plot of isotopologues for formula: " + str(formula_to_simulate))
        isotope_simulation = MS_functions.simulate_isotope_pattern_of_formula(formula_to_simulate)
        print("Isotope simulation: " + str(isotope_simulation))
        self.make_op_log_entry("INFO:\t" + "Simulated isotope pattern will be plotted now. Including actual spectral data: " + str(include_actual_spectral_data))
        self.make_op_log_entry("INFO:\t" + "Simulated isotope pattern: " + str(isotope_simulation))
        fig = plt.figure()
        gs = fig.add_gridspec(1, len(isotope_simulation), hspace=0, wspace=0)
        ax = gs.subplots(sharex="col", sharey="row")
        intensities_to_include = []
        summarized_intensities = [i for m, i in self.spec.summarized_mass_intensity_dict.items()]
        summarized_masses = [m for m, i in self.spec.summarized_mass_intensity_dict.items()]
        for entry in list(isotope_simulation.keys()):
            summed_intensity = sum([summarized_intensities[i] for i in range(len(summarized_masses)) if abs(((summarized_masses[i] - entry) / entry)*1000000) <= (self.kwargs["mass_deviation"])])
            intensities_to_include.append(summed_intensity)
        intensities_to_include = [(i/sum(intensities_to_include)) for i in intensities_to_include]
        
        for entry in range(len(intensities_to_include)):
            if include_actual_spectral_data:
                cmap = matplotlib.colormaps.get_cmap('Greens')
                colors = cmap((summarized_intensities - min(summarized_intensities)) / (max(summarized_intensities) - min(summarized_intensities)) * 0.3 + 0.7)                
                ax[entry].bar(summarized_masses, summarized_intensities, color=colors, label="all ions", width=0.0005, alpha=1)
            ax[entry].bar(list(isotope_simulation.keys())[entry], intensities_to_include[entry], label="measured ions", color="blue", width=0.005, alpha=0.5)
            ax[entry].bar(list(isotope_simulation.keys())[entry], -1 * list(isotope_simulation.values())[entry], label="simulated intensity", color="red", width=0.005, alpha=0.5)
            ax[entry].axhline(0, color='black', linewidth=1)
            ax[entry].set_xlim(list(isotope_simulation.keys())[entry] - 0.1, list(isotope_simulation.keys())[entry] + 0.1)
            ax[entry].set_ylim(-1, 1)
            if entry == 0:  # only set the y-label for the first (leftmost) subplot
                ax[entry].set_ylabel("intensity / a.u.")
            if entry == len(intensities_to_include) // 2:  # only set the x-label for the middle subplot
                ax[entry].set_xlabel("masses / Da")
            ax[entry].set_title(str(round(list(isotope_simulation.keys())[entry], 4)), rotation='vertical')
        for a in fig.get_axes():
            a.label_outer()
        #set the title of the whole figure so that it will be displayed above the subplots
        fig.suptitle("Isotopologues Plot for Formula: " + str(formula_to_simulate) + 
                     "\n At RT: " + str(round(self.spec.rt, 2)) + " and Mass: " + str(round(self.mass, 4)))
        plt.legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        plt.savefig(filepath, bbox_inches='tight', dpi=1000)
        plt.close(fig)
        return True

    def save_xic_plot(self, times, intensities, xic_plot_filepath):
        fig, ax = plt.subplots()
        if self.identified_peaks_for_mass is None:
            self.identified_peaks_for_mass = self.get_peak_properties(times, intensities)
        identified_peak_times = [entry[0] for entry in self.identified_peaks_for_mass]
        plot_heights = [max(intensities)/4 for i in range(len(identified_peak_times))]
        ax.scatter(identified_peak_times, plot_heights, color="green", label="identified peak", s=20, alpha=0.5)
        ax.bar(self.spec.rt, max(intensities), color="red", label="observed time", alpha=0.5, width=10)
        ax.plot(times, intensities, color="blue", label="XIC")
        ax.set_xlabel("retention time / s")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("Extracted Ion Chromatogram")
        ax.legend()
        plt.savefig(xic_plot_filepath, bbox_inches='tight', dpi=1000)
        plt.close(fig)
        return True


class OneAnalysis:
    def __init__(self, ms_file, mass, rt, **kwargs):
        self.mass = mass
        self.rt = rt
        self.ms_file = ms_file
        self.peak_index = self.ms_file.rt_list.index(min(self.ms_file.rt_list, key=lambda x: abs(self.rt - x)))

        default_kwargs = {"one_analysis_folder":str(self.ms_file.parentfolder + "/" + str(self.mass) + "_" + str(self.rt) + "/"),
                          "oa_log_filepath":str(self.ms_file.parentfolder + "/" + str(self.mass) + "_" + str(self.rt) + "/" + "oa_log.txt"),
                          "mass_deviation":11,
                          "charge_of_measured_mass":-1,

                          "oa_save_xic_plot":True,
                          "oa_mass_deviation_xic":(50*self.mass)/1000000,
                          "oa_xic_requested_filter_mode":"Full scan",

                          "oa_spec_acquisition_save_matplotlib_plot":True,
                          "oa_spec_acquisition_save_go_plot":True,
                          "oa_spec_acquisition_save_additional_info":True,

                          "oa_molecular_ion_pred_save_matplotlib_plot_of_isotopologues":True,
                          "oa_molecular_ion_pred_save_go_plot_of_isotopologues":True,
                          "oa_molecular_ion_pred_save_xic_plot":True,
                          "oa_molecular_ion_pred_save_detailed_log":True,
                          "oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found":False,

                          "oa_molecular_ion_reject_formula_if_score_lower_than":10,

                          "oa_molecular_ion_pred_before_spec_acquisition_save_matplotlib_plot":False,
                          "oa_molecular_ion_pred_before_spec_acquisition_save_go_plot":False,
                          "oa_molecular_ion_pred_after_spec_acquisition_save_matplotlib_plot":False,
                          "oa_molecular_ion_pred_after_spec_acquisition_save_go_plot":False,

                          "oa_multiplespec_pred_save_matplotlib_plot_of_isotopologues":False,
                          "oa_multiplespec_pred_save_go_plot_of_isotopologues":False,
                          "oa_multiplespec_pred_save_xic_plot":False,
                          "oa_multiplespec_pred_save_detailed_log":False,

                          "oa_fragments_pred_save_matplotlib_plot_of_isotopologues":True,
                          "oa_fragments_pred_save_go_plot_of_isotopologues":True,
                          "oa_fragments_spec_save_matplotlib_plot":True,
                          "oa_fragments_spec_save_go_plot":False,
                          "oa_fragments_pred_save_xic_plot":True,
                          "oa_fragments_pred_save_detailed_log":True,
                          "oa_fragments_only_calc_prediction_if_peak_is_found":True,

                          "oa_include_frag_intensity_noise_multiplier":0.01,
                          "oa_reject_formula_if_score_lower_than": 10,
                          "oa_make_good_fragment_formula_prediction":False,
                          "pred_formula_cache_folder_path":"C://Users//Admin//Desktop//UVenture//Formula_Predictions//Formula_Predictions//" }
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["one_analysis_folder"], exist_ok=True)

        print()
        print("OneAnalysis kwargs: " + str(self.kwargs))
        print()

        self.make_oa_log_entry("")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("INFO:\t" + "Creating OneAnalysis object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.peak_index))
        self.make_oa_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))
        self.make_oa_log_entry("INFO:\t" + "Assuming a charge of the measured mass of: " + str(self.kwargs["charge_of_measured_mass"]))
        self.make_oa_log_entry("INFO:\t" + "Available MS modes: " + str(self.ms_file.available_modes))

        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, self.kwargs["oa_mass_deviation_xic"], requested_filter_mode=self.kwargs["oa_xic_requested_filter_mode"])
        if self.kwargs["oa_save_xic_plot"] == True:
            print("Saving XIC plot...")
            self.xic_plot_filepath = self.kwargs["one_analysis_folder"] + "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["oa_mass_deviation_xic"]) + "_" + str(
                                    self.kwargs["oa_xic_requested_filter_mode"]) + ".png"
            self.plot_and_save_xic(self.xic_plot_filepath, self.peak_index, self.xic)
        
        print(kwargs)

        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue
        
        print(additional_kwargs)
        print()
        print(kwargs)

        if "Full scan" in self.ms_file.available_modes:
            self.full_scan_spec = Spec(self.ms_file, 
                                             self.peak_index, 
                                             spec_requested_filter_mode="Full scan", 
                                             absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                             mass_deviation=self.kwargs["mass_deviation"],
                                             spec_save_matplotlib_plot=self.kwargs["oa_spec_acquisition_save_matplotlib_plot"],
                                             spec_save_go_plot=self.kwargs["oa_spec_acquisition_save_go_plot"],
                                             **additional_kwargs)
            print("Full scan spec found!")
            self.best_molecular_ion_spec = self.full_scan_spec
            self.best_frag_spec = self.full_scan_spec
        else:
            self.full_scan_spec = None
            self.best_frag_spec = None
            self.make_oa_log_entry("WARNING:\t" + "No full scan spectrum found!")
            print("No full scan spectrum found!")
            self.best_molecular_ion_spec = None
        
        if "AIF" in self.ms_file.available_modes:
            self.aif_spec = Spec(self.ms_file, 
                                       self.peak_index, 
                                       spec_requested_filter_mode="AIF", 
                                       absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                       mass_deviation=self.kwargs["mass_deviation"],
                                       spec_save_matplotlib_plot=self.kwargs["oa_spec_acquisition_save_matplotlib_plot"],
                                       spec_save_go_plot=self.kwargs["oa_spec_acquisition_save_go_plot"],
                                       **additional_kwargs)
            print("AIF spec found!")
            self.best_frag_spec = self.aif_spec
            if self.best_molecular_ion_spec is None:
                self.best_molecular_ion_spec = self.aif_spec
        else:
            self.aif_spec = None
            self.make_oa_log_entry("WARNING:\t" + "No AIF spectrum found!")
            print("No AIF spectrum found!")
        
        if "MS/MS" in self.ms_file.available_modes:
            self.ms_ms_spec = Spec(self.ms_file, 
                                         self.peak_index, 
                                         spec_requested_filter_mode="MS/MS", 
                                         absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                         mass_deviation=self.kwargs["mass_deviation"], 
                                         spec_requested_ms_ms_mass=self.mass,
                                         spec_save_matplotlib_plot=self.kwargs["oa_spec_acquisition_save_matplotlib_plot"],
                                         spec_save_go_plot=self.kwargs["oa_spec_acquisition_save_go_plot"],
                                         **additional_kwargs)
            print("MS/MS spec found!")
            if self.best_molecular_ion_spec is None:
                self.best_molecular_ion_spec = self.ms_ms_spec
            if self.best_molecular_ion_spec is None:
                self.best_molecular_ion_spec = self.ms_ms_spec
        else:
            self.ms_ms_spec = None
            self.make_oa_log_entry("INFO:\t" + "No MS/MS spectrum found!")
            print("No MS/MS spectrum found!")
        
        if self.best_molecular_ion_spec is None:
            self.make_oa_log_entry("ERROR:\t" + "No spectrum found to use as molecular ion prediction spectrum!")
            self.make_oa_log_entry("ERROR:\t" + "Stopping prediction...")
            print("No spectrum found to use as molecular ion prediction spectrum!")
            print("Stopping prediction...")
            return
        if self.best_frag_spec is None:
            self.make_oa_log_entry("ERROR:\t" + "No spectrum found to use as fragment ion prediction spectrum!")
            self.make_oa_log_entry("ERROR:\t" + "Continuing prediction without fragment ion prediction...")
            print("No spectrum found to use as fragment ion prediction spectrum!")
            print("Continuing prediction without fragment ion prediction...")

        self.mass_old = self.mass
        self.mass = self.adjust_mass_to_closest_measured_mass(self.best_molecular_ion_spec)
        self.make_oa_log_entry("INFO:\t" + "Adjusting the mass of the molecular ion...")
        self.make_oa_log_entry("INFO:\t" + "New mass: " + str(self.mass) + "  Old mass: " + str(self.mass_old))
        self.mass = self.mass + (self.kwargs["charge_of_measured_mass"] * 0.000548)
        self.make_oa_log_entry("INFO:\t" + "Charge of the measured mass: " + str(self.kwargs["charge_of_measured_mass"]))
        self.make_oa_log_entry("INFO:\t" + "New mass after charge correction: " + str(self.mass))
        self.intensity_of_molecular_ion = sum([self.best_molecular_ion_spec.summarized_intensities[i] for i in range(len(self.best_molecular_ion_spec.summarized_intensities)) if abs(((self.best_molecular_ion_spec.summarized_masses[i] - self.mass)/self.mass)*1000000) <= self.kwargs["mass_deviation"]])
        self.make_oa_log_entry("INFO:\t" + "Intensity of the molecular ion: " + str(self.intensity_of_molecular_ion))

        self.molecular_ion_prediction = self.get_molecular_ion_prediction(self.best_molecular_ion_spec)
        if self.molecular_ion_prediction.peak_found == False and self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"] == True:
            self.make_oa_log_entry("INFO:\t" + "No molecular ion prediction found. Stopping prediction...")
            return
        
        self.make_oa_log_entry("INFO:\t" + "Finished prediction of molecular ion...")
        self.make_oa_log_entry("INFO:\t" + "Summarized molecular ion formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))
        self.make_oa_log_entry("INFO:\t" + "Predicted molecular ion: " + str(self.best_molecular_ion_prediction))
        self.make_oa_log_entry("INFO:\t" + "Predicted molecular ion score: " + str(self.score_of_best_molecular_ion_prediction))
        self.make_oa_log_entry("INFO:\t" + "Ion intensity: " + str(self.molecular_ion_prediction.intensity_of_ion))
        print("Finished prediction of molecular ion! Molecular ion prediction: " + str(self.best_molecular_ion_prediction))

        self.fragment_predictions, self.fragment_predictions_formula_score_dicts = self.get_fragment_predictions(make_good_fragment_formula_prediction=self.kwargs["oa_make_good_fragment_formula_prediction"])

        self.make_oa_log_entry("INFO:\t" + "Finished prediction of fragment ions...")
        self.make_oa_log_entry("INFO:\t" + "Fragment predictions: " + str(self.fragment_predictions))
        self.make_oa_log_entry("INFO:\t" + "Fragment predictions formula score dicts: " + str(self.fragment_predictions_formula_score_dicts))

        self.summary_list_of_one_analysis = self.make_true_fragment_list()
        print("Summary list of one analysis: " + str(self.summary_list_of_one_analysis))
        
        self.make_oa_log_entry("INFO:\t" + "Finished creating summary list of one analysis...")
        self.make_oa_log_entry("INFO:\t" + "Summary list of one analysis: " + str(self.summary_list_of_one_analysis))
        
        self.append_oa_summary_to_raw_file_summary(self.summary_list_of_one_analysis)

    def append_oa_summary_to_raw_file_summary(self, summary_list):
        with open(self.ms_file.parentfolder + "/" + "SUMMARY.txt", "a") as oa_summary:
            oa_summary.write(str(summary_list) + "\n")

    def make_true_fragment_list(self):
        print()
        print()
        print("MAKING TRUE FRAGMENT LIST")
        true_fragment_list = []
        
        molecular_ion_pred_object = self.molecular_ion_prediction
        molecular_ion_mass = self.mass
        molecular_ion_best_approx = self.best_molecular_ion_prediction
        molecular_ion_best_approx_dict = MS_functions.get_formula_to_dict(molecular_ion_best_approx)
        molecular_ion_score = self.score_of_best_molecular_ion_prediction
        molecular_ion_intensity = self.intensity_of_molecular_ion
        molecular_ion_formula_score_dict = self.summarized_molecular_ion_formula_score_dict
        mi_list = [molecular_ion_mass, molecular_ion_best_approx, molecular_ion_score, molecular_ion_intensity, molecular_ion_formula_score_dict]
        true_fragment_list.append(mi_list)
        print("Molecular ion list: " + str(mi_list))
        print("Starting for loop...")
        for f_mass, f_score_dict in self.fragment_predictions_formula_score_dicts.items():
            print()
            f_intensity = self.fragment_predictions[f_mass].intensity_of_ion
            neutral_loss_mass = molecular_ion_mass - f_mass
            neutral_loss_formula_predictions = (MS_functions.get_formula_from_cache(self.kwargs["pred_formula_cache_folder_path"], neutral_loss_mass, self.kwargs["mass_deviation"]))[1]
            # neutral_loss_formula_prediction = [mass, alldict] ---> alldict = {formula: score, "C1H3O2": score, ...}
            if len(neutral_loss_formula_predictions) == 0:
                continue
            found_pair = False
            print(neutral_loss_formula_predictions)
            for f_formula, f_score in f_score_dict.items():
                f_formula = ast.literal_eval(f_formula)
                for nl_formula, nl_deviation in neutral_loss_formula_predictions.items():
                    nl_deviation = neutral_loss_formula_predictions[nl_formula]
                    nl_formula = MS_functions.get_formula_to_dict(nl_formula)
                    print(nl_formula)
                    print(f_formula)
                    print(molecular_ion_best_approx_dict)
                    if (self.combine_and_sum_dicts(nl_formula, f_formula) == molecular_ion_best_approx_dict):
                        print("TRUETRUETRUEskjaskjhlgfaivwzbevwuief")
                        true_fragment_list.append([f_mass, f_formula, f_score, f_intensity, neutral_loss_mass, nl_formula, nl_deviation])
                        found_pair = True
                        break
                if found_pair == True:
                    break
        return true_fragment_list
        
    def adjust_mass_to_closest_measured_mass(self, spec, masse=0):
        if masse == 0:
            masse = self.mass
        mass = min(list(spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(masse - x))
        return mass

    def combine_and_sum_dicts(self, dict1, dict2):
        return {k: dict1.get(k, 0) + dict2.get(k, 0) for k in set(dict1) | set(dict2)}

    def get_molecular_ion_prediction(self, best_molecular_ion_spec=None):
        if best_molecular_ion_spec is None:
            print("No best molecular ion prediction spec provided. Returning...")
            return None
        available_specs = []
        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["absolute_pred_folder", "pred_formula_cache_folder_path", "pred_save_matplotlib_plot_of_isotopologues", "pred_save_go_plot_of_isotopologues", "pred_save_xic_plot", "pred_save_detailed_log", "pred_return_if_no_peak_is_found"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue

        self.molecular_ion_prediction = Prediction(self.ms_file, self.mass, best_molecular_ion_spec, 
                                                   absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/",
                                                   pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                                   pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_matplotlib_plot_of_isotopologues"],
                                                   pred_save_go_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_go_plot_of_isotopologues"],
                                                   pred_save_xic_plot=self.kwargs["oa_molecular_ion_pred_save_xic_plot"],
                                                   pred_save_detailed_log=self.kwargs["oa_molecular_ion_pred_save_detailed_log"],
                                                   pred_return_if_no_peak_is_found=self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"],
                                                   **additional_kwargs)
        #the object that is returned, has a bool variable self.peak_found = False if no peak was found or =True if a peak was found
        available_specs.append(best_molecular_ion_spec)
        #try to get the two full scan spectra next to the provided spectrum
        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue

        try:
            spec_before = Spec(self.ms_file, 
                               self.peak_index - len(self.ms_file.available_modes), 
                               spec_requested_filter_mode=best_molecular_ion_spec.filter_mode, 
                               absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                               mass_deviation=self.kwargs["mass_deviation"],
                               spec_save_matplotlib_plot=self.kwargs["oa_molecular_ion_pred_before_spec_acquisition_save_matplotlib_plot"],
                               spec_save_go_plot=self.kwargs["oa_molecular_ion_pred_before_spec_acquisition_save_go_plot"],
                               **additional_kwargs)
            print("Spec before found!")
            available_specs.append(spec_before)
        except Exception as e:
            spec_before = None
            print("Error in getting spectrum & prediction before the actual analyzed spectrum: " + str(e))
            print(traceback.format_exc())
        try:
            spec_after = Spec(self.ms_file, 
                              self.peak_index + len(self.ms_file.available_modes), 
                              requested_filter_mode=best_molecular_ion_spec.filter_mode, 
                              absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                              mass_deviation=self.kwargs["mass_deviation"],
                              spec_save_matplotlib_plot=self.kwargs["oa_molecular_ion_pred_after_spec_acquisition_save_matplotlib_plot"],
                              spec_save_go_plot=self.kwargs["oa_molecular_ion_pred_after_spec_acquisition_save_go_plot"],
                              **additional_kwargs)
            print("Spec after found!")
            available_specs.append(spec_after)
        except Exception as e:
            spec_after = None
            print("Error in getting spectrum & prediction after the actual analyzed spectrum: " + str(e))
            print(traceback.format_exc())
        
        #add the scores of all the available prediction formula_score_dicts and create a summarized formula_score_dict
        self.summarized_molecular_ion_formula_score_dict, all_formula_score_dicts = self.get_formula_score_dict_with_multiple_specs(available_specs, self.mass, return_if_no_peak_found=self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"])
        
        self.summarized_molecular_ion_formula_score_dict = {f: s for f, s in self.summarized_molecular_ion_formula_score_dict.items() if s > self.kwargs["oa_molecular_ion_reject_formula_if_score_lower_than"]}
        self.summarized_molecular_ion_formula_score_dict = dict(sorted(self.summarized_molecular_ion_formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        self.make_oa_log_entry("INFO:\t" + "Summarized formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))
        print("Summarized formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))

        self.molecular_ion_prediction.make_op_log_entry("=============================================================")
        try:
            self.best_molecular_ion_prediction = str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0]).items()]))  
            self.score_of_best_molecular_ion_prediction = list(self.summarized_molecular_ion_formula_score_dict.values())[0]
        except Exception as e:
            self.best_molecular_ion_prediction = "No formula found"
            self.score_of_best_molecular_ion_prediction = -9999
            self.molecular_ion_prediction.make_op_log_entry("ERROR:\t" + "No formula found!")
            print("Error in getting best molecular ion prediction: " + str(e))
            print(traceback.format_exc())
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Formula predicted with " + str(len(available_specs)) + " spectra.")
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "All formula score dicts: ")
        for i, f_score_dict in enumerate(all_formula_score_dicts):
            self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Formula score dict " + str(i) + ": " + str(f_score_dict))
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Best formula prediction: " + str(self.best_molecular_ion_prediction) + " Score: " + str(self.score_of_best_molecular_ion_prediction))
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Summarized formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))
        self.molecular_ion_prediction.make_op_log_entry("=============================================================")
        return self.molecular_ion_prediction
    
    def get_formula_score_dict_with_multiple_specs(self, specs, mass, return_if_no_peak_found=True):
        predictions = []
        for spec in specs:
            try:
                additional_kwargs = copy.deepcopy(self.kwargs)
                pop_keys = ["absolute_pred_folder", "pred_formula_cache_folder_path", "pred_save_matplotlib_plot_of_isotopologues", "pred_save_go_plot_of_isotopologues", "pred_save_xic_plot", "pred_save_detailed_log", "pred_return_if_no_peak_is_found"]
                for key in pop_keys:
                    try:
                        additional_kwargs.pop(key)
                    except KeyError:
                        continue

                pred = Prediction(self.ms_file, mass, spec,
                                  absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/",
                                  pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                  pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_multiplespec_pred_save_matplotlib_plot_of_isotopologues"],
                                  pred_save_go_plot_of_isotopologues=self.kwargs["oa_multiplespec_pred_save_go_plot_of_isotopologues"],
                                  pred_save_xic_plot=self.kwargs["oa_multiplespec_pred_save_xic_plot"],
                                  pred_save_detailed_log=self.kwargs["oa_multiplespec_pred_save_detailed_log"],
                                  pred_return_if_no_peak_is_found=return_if_no_peak_found,
                                  **additional_kwargs)
                predictions.append(pred)
            except Exception as e:
                predictions.append(None)
                print("Error in getting prediction with multiple specs: " + str(e))
                print(traceback.format_exc())
                continue
        print(predictions)
        summarized_formula_score_dict = {}
        for pred in predictions:
            if pred is not None:
                summarized_formula_score_dict = self.combine_and_sum_dicts(summarized_formula_score_dict, pred.formula_score_dict)
                print("Formula score dict after adding prediction: " + str(summarized_formula_score_dict))
        summarized_formula_score_dict = dict(sorted(summarized_formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        all_formula_score_dicts = [pred.formula_score_dict for pred in predictions if pred is not None]
        return summarized_formula_score_dict, all_formula_score_dicts
        
    def get_fragment_predictions(self, make_good_fragment_formula_prediction=False):
        if self.best_frag_spec is None:
            self.make_oa_log_entry("INFO:\t" + "No best fragment prediction spec provided. Returning...")
            print("No best fragment prediction spec provided. Returning...")
            return None
        if make_good_fragment_formula_prediction:
            additional_kwargs = copy.deepcopy(self.kwargs)
            pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
            for key in pop_keys:
                try:
                    additional_kwargs.pop(key)
                except KeyError:
                    continue

            available_specs = []
            available_specs.append(self.aif_spec)
            try:
                spec_before = Spec(self.ms_file, 
                                   self.peak_index - len(self.ms_file.available_modes), 
                                   spec_requested_filter_mode=self.best_frag_spec.filter_mode, 
                                   absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                   mass_deviation=self.kwargs["mass_deviation"],
                                   spec_save_matplotlib_plot=self.kwargs["oa_fragments_spec_save_matplotlib_plot"],
                                   spec_save_go_plot=self.kwargs["oa_fragments_spec_save_go_plot"],
                                   **additional_kwargs)
                print("Spec before found!")
                available_specs.append(spec_before)
            except Exception as e:
                spec_before = None
                print("Error in getting spectrum & prediction before the actual analyzed spectrum: " + str(e))
                print(traceback.format_exc())
            try:
                spec_after = Spec(self.ms_file, 
                                  self.peak_index + len(self.ms_file.available_modes), 
                                  spec_requested_filter_mode=self.best_frag_spec.filter_mode, 
                                  absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                  mass_deviation=self.kwargs["mass_deviation"],
                                  spec_save_matplotlib_plot=self.kwargs["oa_fragments_spec_save_matplotlib_plot"],
                                  spec_save_go_plot=self.kwargs["oa_fragments_spec_save_go_plot"],
                                  **additional_kwargs)
                print("Spec after found!")
                available_specs.append(spec_after)
            except Exception as e:
                spec_after = None
                print("Error in getting spectrum & prediction after the actual analyzed spectrum: " + str(e))
                print(traceback.format_exc())
        print("Starting prediction of fragment ions...")
        self.possible_fragment_masses = [m for m in self.best_frag_spec.summarized_masses if m < self.mass and self.best_frag_spec.summarized_intensities[self.best_frag_spec.summarized_masses.index(m)] > (self.intensity_of_molecular_ion * self.kwargs["oa_include_frag_intensity_noise_multiplier"])]
        print("Possible fragment masses: " + str(self.possible_fragment_masses))
        self.make_oa_log_entry("INFO:\t" + "Possible fragment masses: " + str(self.possible_fragment_masses))
        self.fragment_predictions = {}
        self.fragment_predictions_formula_score_dicts = {}
        for frag_mass in self.possible_fragment_masses:
            additional_kwargs = copy.deepcopy(self.kwargs)
            pop_keys = ["absolute_pred_folder", "pred_formula_cache_folder_path", "pred_save_matplotlib_plot_of_isotopologues", "pred_save_go_plot_of_isotopologues", "pred_save_xic_plot", "pred_save_detailed_log", "pred_return_if_no_peak_is_found"]
            for key in pop_keys:
                try:
                    additional_kwargs.pop(key)
                except KeyError:
                    continue

            curr_prediction = Prediction(self.ms_file, frag_mass, self.best_frag_spec, 
                                                              absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/fragments/",
                                                              pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                                              pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_fragments_pred_save_matplotlib_plot_of_isotopologues"],
                                                              pred_save_go_plot_of_isotopologues=self.kwargs["oa_fragments_pred_save_go_plot_of_isotopologues"],
                                                              pred_save_xic_plot=self.kwargs["oa_fragments_pred_save_xic_plot"],
                                                              pred_save_detailed_log=self.kwargs["oa_fragments_pred_save_detailed_log"],
                                                              pred_return_if_no_peak_is_found=self.kwargs["oa_fragments_only_calc_prediction_if_peak_is_found"],
                                                              **additional_kwargs)
            if not curr_prediction.prediction_spec_within_peak_range and \
                self.kwargs["oa_fragments_only_calc_prediction_if_peak_is_found"] == True:
                continue
            self.fragment_predictions[frag_mass] = curr_prediction

            if make_good_fragment_formula_prediction:
                self.fragment_predictions_formula_score_dicts[frag_mass], _ = self.get_formula_score_dict_with_multiple_specs(available_specs, frag_mass, return_if_no_peak_found=False)
            else:
                self.fragment_predictions_formula_score_dicts[frag_mass] = self.fragment_predictions[frag_mass].formula_score_dict
            
            self.make_oa_log_entry("INFO:\t" + "Finished prediction of fragment ion with mass: " + str(frag_mass) + " Best formula prediction: " + str(self.fragment_predictions[frag_mass].best_formula_prediction) + " Score: " + str(self.fragment_predictions[frag_mass].score_of_best_formula))
        return self.fragment_predictions, self.fragment_predictions_formula_score_dicts

    def plot_and_save_xic(self, save_filepath, peak_index, xic):
        fig, ax = plt.subplots()
        xic_peak_index = xic[2].index(min(xic[2], key=lambda x: abs(peak_index - x)))
        ax.scatter(self.ms_file.rt_list[peak_index], xic[1][xic_peak_index], color="red", marker="x", s=100)
        ax.plot(xic[0], xic[1], label="XIC measured")
        ax.set_title("xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["oa_mass_deviation_xic"]) + "_" + str(self.kwargs["oa_xic_requested_filter_mode"]))
        ax.set_xlabel("retention time / seconds")
        ax.set_ylabel("intensity / a.u.")
        plt.legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        plt.savefig(save_filepath, bbox_inches='tight', dpi=1000)
        plt.close(fig)
        metadata = PIL.PngImagePlugin.PngInfo()
        metadata.add_text("xic_retention_time", str(xic[0]))
        metadata.add_text("xic_intensity_list", str(xic[1]))
        metadata.add_text("xic_original_index_list", str(xic[2]))
        target_image = PIL.Image.open(save_filepath)
        target_image.save(save_filepath, pnginfo=metadata)

    def make_oa_log_entry(self, log_entry):
        #get the dirname of the spec_log_filepath
        directory_logfile = os.path.dirname(self.kwargs["oa_log_filepath"])
        #create the directory if it does not exist
        os.makedirs(directory_logfile, exist_ok=True)
        log_f = open(self.kwargs["oa_log_filepath"], "a")
        log_f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t" + log_entry + "\n")
        log_f.close()
        return True


class PeakFinding:
    def __init__(self, ms_file, **kwargs):
        self.ms_file = ms_file
        self.peak_properties_list = []



        default_kwargs = {"mass_deviation":11,
                          "minimum_required_max_intensity":100000,
                          "minimum_peak_width_seconds":0.1,
                          "maximum_peak_width_seconds":10,
                          "minimum_peak_prominence":50000,
                          "peak_finding_log_filepath":str(self.ms_file.parentfolder + "/" + "peak_finding_log.txt"),
                          "round_for_masses_summary": 4 }
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
        xic = MS_functions.get_xic(self.ms_file.rawdata, mass, ((self.kwargs["mass_deviation"]*mass)/1000000), requested_filter_mode="Full scan")
        times = xic[0]
        intensities = xic[1]
        window, order = 5, 3
        intensities = scipy.signal.savgol_filter(intensities, window, order, mode="nearest")
        peak_properties = scipy.signal.find_peaks(intensities, height=max(intensities)/100, distance=2, prominence=max(intensities)/50, width=(2, 20))
        identified_peaks = []
        for element in range(len(peak_properties[1]["peak_heights"])):
            if peak_properties[1]["peak_heights"][element] < self.kwargs["minimum_required_max_intensity"] or \
                peak_properties[1]["prominences"][element] < self.kwargs["minimum_peak_prominence"]:
                continue
            one_peak = []
            one_peak.append(mass)
            one_peak.append(times[int(peak_properties[1]["left_ips"][element] + (peak_properties[1]["widths"][element]/2))])
            one_peak.append(peak_properties[1]["peak_heights"][element])
            one_peak.append(times[int(peak_properties[1]["right_ips"][element])] - times[int(peak_properties[1]["left_ips"][element])])
            one_peak.append(peak_properties[1]["prominences"][element])
            identified_peaks.append(one_peak)
        #[[mass, time, height, duration, prominence], [....], ...]
        return identified_peaks

    def combine_and_sum_dicts(self, dict1, dict2):
        return {k: dict1.get(k, 0) + dict2.get(k, 0) for k in set(dict1) | set(dict2)}

    def get_best_approx_for_ppm_spacing_within_peak(self, mass_list, worst_expected_ppm_deviation=20):
        mass_list = sorted(mass_list)
        ppm_spacing_list = [(((mass_list[i+1] - mass_list[i]) / mass_list[i]) * 1000000) for i in range(len(mass_list)-1)]
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
        #sort the dictio by its keys
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
                mass_lower_border = curr_mass - ((deviation*curr_mass)/1000000)
                mass_upper_border = curr_mass + ((deviation*curr_mass)/1000000)
                iteration_step_lower = 0
                integration_step_upper = 0
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                expected_min_spacing_between_measurement_points = (best_approx_ppm_spacing_within_peak * curr_mass) / 1000000
                try:
                    while (mass_lower_border < old_masses_list[index_of_highest_abundance-iteration_step_lower]) and \
                            (last_existing_mass_within_border - old_masses_list[index_of_highest_abundance-iteration_step_lower] <= 2.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance-iteration_step_lower] <= (1.1 * last_abundance)):
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance-iteration_step_lower]
                        last_abundance = old_abundances_list[index_of_highest_abundance-iteration_step_lower]
                        iteration_step_lower += 1
                except IndexError:
                    remove_all_lower = True
                    iteration_step_lower = 0
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                try:
                    while (mass_upper_border > old_masses_list[index_of_highest_abundance+integration_step_upper]) and \
                            (old_masses_list[index_of_highest_abundance+integration_step_upper] - last_existing_mass_within_border <= 2.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance+integration_step_upper] <= (1.1 * last_abundance)):
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance+integration_step_upper]
                        last_abundance = old_abundances_list[index_of_highest_abundance+integration_step_upper]
                        integration_step_upper += 1
                except IndexError:
                    remove_all_upper = True
                    integration_step_upper = 0
                
                if remove_all_lower == True and remove_all_upper == True:
                    break
                if remove_all_upper:
                    summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]) / summed_intensity
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_masses_list), 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]
                    continue
                if remove_all_lower:
                    summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper, 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper, 1)]) / summed_intensity
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper, 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper, 1)]
                    continue
                    
                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)]) / summed_intensity
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper, 1)]
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
            curr_save = {round(key, self.kwargs["round_for_masses_summary"]): value for key, value in curr_save.items() if value > self.kwargs["minimum_required_max_intensity"]/10}
            ordered_mass_intensity_dict = self.combine_and_sum_dicts(ordered_mass_intensity_dict, curr_save)
            ordered_mass_intensity_dict = {round(key, self.kwargs["round_for_masses_summary"]): value for key, value in ordered_mass_intensity_dict.items() if value > self.kwargs["minimum_required_max_intensity"]/10}
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
    
        
            

            
if __name__ == "__main__":
    settings_dict_filepath = "static/settings.txt"
    settings_dict = {}
    with open(settings_dict_filepath, "r") as f:
        file_contents_raw = f.read()
        lines = file_contents_raw.split("\n")
        for line in lines:
            line = line.strip()
            if not "=" in line:
                continue
            if line[0] == "#":
                continue
            key, value = line.split("=")
            try:
                value = float(value)
            except:
                try:
                    value = int(value)
                except:
                    try:
                        value = ast.literal_eval(value)
                    except:
                        value = str(value)
            settings_dict[key] = value  
    for key in settings_dict:
        print(key + ": " + str(settings_dict[key]))

    settings_dict["pred_formula_cache_folder_path"] = "C://Users//Admin//Desktop//UVenture//Formula_Predictions//"

    mzml_filename = "C://Users//Admin//Desktop//UVenture//Evaluation_folder//F12_HRAIF4_3.mzML"

    ms_file = MS_File(mzml_filename)

    #find_peaks = PeakFinding(ms_file)
    #myspec = Spec(ms_file, 400, requested_filter_mode="Full scan", save_plot=True, unique_spec_folder="test/")

    #print(myspec.index)
    #print(myspec.masses)
    #print(myspec.intensities)
    #print(myspec.filter)
    #print(myspec.filter_mode)
    #print(myspec.ms_ms_masses)
    #print(myspec.spec_rawdata)
    #print(myspec.ms_level)

    #print("=====================================================================================================")
    #print("=====================================================================================================")

    
    #myprediction = Prediction(ms_file, 117.0554, myspec, save_plot_of_isotopologues=True, charge_of_measured_mass=-1, formula_cache_folder_path="C://Users//Admin//Desktop//UVenture//Formula_Predictions//Formula_Predictions//")

    myanalysis = OneAnalysis(ms_file, 117.0554, ms_file.rt_list[79], **settings_dict)

    print(myanalysis.molecular_ion_prediction.best_formula_prediction)
    print(myanalysis.molecular_ion_prediction.score_of_best_formula)
    print(myanalysis.molecular_ion_prediction.intensity_of_ion)


    #print(myprediction.mass_possible_formulas)
    #print(myprediction.formulas_score_dict)
    #print(myprediction.mass)
    #print(myprediction.intensity_of_ion)


    







