import ast
import copy
import datetime
import math
import pathlib
import sys
import traceback
import PIL
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import matplotlib
matplotlib.use('Agg')
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

DPI = 300


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
        self.make_spec_log_entry("Summarized_mass_intensity_dict: ")
        self.make_spec_log_entry(str(self.summarized_mass_intensity_dict))

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
                    if summed_intensity >= 500:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_masses_list), 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]
                    continue
                if remove_all_lower:
                    summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper+1, 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]) / summed_intensity
                    if summed_intensity >= 200:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]
                    continue

                if iteration_step_lower == 0 and integration_step_upper == 0:
                    summed_intensity = old_abundances_list[index_of_highest_abundance]
                    weighted_mass_average = old_masses_list[index_of_highest_abundance]
                    if summed_intensity >= 500:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance ]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance ]
                    continue

                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]) / summed_intensity
                if summed_intensity >= 500:
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]
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
        fig = matplotlib.figure.Figure(layout="tight")
        ax = fig.subplots()
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
        ax.set_title("Mode:" + str(self.filter_mode) + "; Index: " + str(self.index) + "; RT: " + str(round(self.ms_file.rt_list[self.index], 2)) + ";\nFilter: " + str(self.filter) + ";\nMS/MS masses: " + str(self.ms_ms_masses))
        matplotlib.rcParams.update({'figure.autolayout': True})
        image_filepath = self.kwargs["spec_folder"] + "Mass_spectrum_index" + str(self.index) + "_" + str(self.filter_mode.replace("/", "")) + "_" + str(self.ms_ms_masses) + ".png"
        fig.savefig(image_filepath, bbox_inches='tight', dpi=DPI)
        fig.clf()
        fig.clear()
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
    def __init__(self, ms_file, mass, spec, spec_before=None, spec_after=None, **kwargs):
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
        self.simulated_isotopologue_pattern_for_best_formula = None

        self.spec_before = spec_before
        self.spec_after = spec_after

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
                          "pred_multiplier_isotopologue_influence_of_too_little_intensity_on_score": 100,
                          "pred_stop_isotopologue_search_if_score_lower_than":-80,
                          "pred_ppm_deviation_score_multiplier":14,
                          "pred_include_likelyhood_of_formula":True,
                          "pred_reject_formula_if_score_lower_than":10,
                          "pred_return_if_no_peak_is_found": True,
                          "pred_include_peak_matching_of_isotopologues": True,
                          "pred_formula_likelihood_substract_score_rdbe_non_integer": 1000,
                          "pred_formula_likelihood_substract_score_rdbe_higher_than_senior_rule": 10000,

                          "oa_fragments_do_peak_computation_if_area_higher_than": 4}
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["pred_folder"], exist_ok=True)

        self.intensity_of_ion = sum([self.spec.intensities[i] for i in range(len(self.spec.masses)) if abs(((self.spec.masses[i] - self.mass)/self.mass)*1000000) <= self.kwargs["mass_deviation"]])
        self.make_op_log_entry("")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("=============================================================")   
        self.make_op_log_entry("INFO:\t" + "Creating prediction object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.spec.index))    
        self.make_op_log_entry("INFO:\t" + "Intensity of the molecular ion: " + str(self.intensity_of_ion))
        self.make_op_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))

        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, ((self.kwargs["mass_deviation"]*self.mass)/1000000), requested_filter_mode=self.spec.filter_mode)
        self.make_op_log_entry("INFO:\t" + "XIC calculated.")
        self.identified_peaks_for_mass = self.get_peak_properties(self.xic[0], self.xic[1])
        self.make_op_log_entry("INFO:\t" + "Identified peaks in XIC.")
        self.prediction_spec_within_peak_range = self.check_if_provided_spectrum_within_peak_range(self.spec.index, self.identified_peaks_for_mass)
        self.make_op_log_entry("INFO:\t" + "Checked if the provided spectrum is within the peak range of the XIC. Result: " + str(self.prediction_spec_within_peak_range))
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


        neutral_mass = self.mass + (0.000548 * self.kwargs["charge_of_measured_mass"])
        _, possible_formulas = MS_functions.get_formula_from_cache(self.kwargs["pred_formula_cache_folder_path"], neutral_mass, self.kwargs["mass_deviation"])
        self.make_op_log_entry("INFO:\t" + "Retrieved formulas from cache.")
        if self.kwargs["charge_of_measured_mass"] < 0:
            self.make_op_log_entry("INFO:\t" + "Charge of the measured mass is negative. Removing formulas with Li, Na and K...")
            possible_formulas = {key: value for key, value in possible_formulas.items() if "Na" not in key and "K" not in key and "Li" not in key and "Mg" not in key and "Ca" not in key and "Be" not in key}


        atoms_to_keep_in_prediction = ["C", "H", "N", "O", "S", "Se", "Cl", "Ni", "Cu", "Zn", "Br"]
        possible_formulas_dict_list = [MS_functions.get_formula_to_dict(formula_string=f) for f, dev in possible_formulas.items()]
        possible_deviations_list = [dev for f, dev in possible_formulas.items()]
        new_possible_formulas = {}
        for i in range(len(possible_deviations_list)):
            if set(list(possible_formulas_dict_list[i].keys())) - set(atoms_to_keep_in_prediction):
                continue
            else:
                new_possible_formulas[MS_functions.get_formula_string_from_dict(possible_formulas_dict_list[i])] = possible_deviations_list[i]
        possible_formulas = new_possible_formulas

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

            ppm_dev_subtract = abs(self.kwargs["pred_ppm_deviation_score_multiplier"] * ((self.mass * 0.1/100) * (entry[1] ** 2).real))

            formula_score = sum(e_s[3] for e_s in list(isotope_check.values())) - ppm_dev_subtract
            self.make_op_log_entry("INFO:\t" + "Formula: " + str("".join([str(a) + str(n) for a, n in entry[0].items()])) + " \t Score: " + str(round(formula_score, 2)))
            self.formula_score_dict[str(entry[0])] = formula_score
            if self.kwargs["pred_save_detailed_log"] == True:
                self.make_op_log_entry("OUTPUT OF FORMULA PREDICTION WITH MASS SPECTRUM (ISOTOPOLOGUES)")
                self.make_op_log_entry("Mass of ion: " + str(self.mass))
                self.make_op_log_entry("PPM DEVIATION OF MOLECULE FORMULA: " + str(round(entry[1], 2)))
                self.make_op_log_entry("{:>25} \t {:>15} \t {:>17} \t {:>20} \t {:>20} \t {:>25}".format("Mass_of_isotopologue", "Measured_mass", "Isotopo_found?", "measured_intensity", "required_intensity", "score_of_isotopologue"))
                for item in list(isotope_check.items()):
                    try:
                        self.make_op_log_entry("{:>25} \t {:>15} \t {:>17} \t {:>20} \t {:>20} \t {:>25}".format(str(round(item[0], 5)), str(round(item[1][4], 5)), str(item[1][0]), str(round(item[1][2], 2)), str(round(item[1][1], 2)), str(round(item[1][3], 2))))
                    except Exception as e:
                        continue
        
        self.formula_score_dict = dict(sorted(self.formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICITIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>30} \t {:>20}".format(str(item[0]), str(round(item[1], 2))))
        
        if self.kwargs["pred_include_likelyhood_of_formula"] == True:
            self.formula_score_dict = self.add_likelyhood_of_formula_to_score(self.formula_score_dict)


        self.formula_score_dict = {f: s for f, s in self.formula_score_dict.items() if s > self.kwargs["pred_reject_formula_if_score_lower_than"]}
        self.formula_score_dict = dict(sorted(self.formula_score_dict.items(), key=lambda x: x[1], reverse=True))

        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICITIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>30} \t {:>20}".format(str(item[0]), str(round(item[1], 2))))

        if self.kwargs["pred_include_peak_matching_of_isotopologues"]:
            self.formula_score_dict = self.get_peak_matching_of_isotopologues(self.formula_score_dict)
            self.formula_score_dict = {f: s for f, s in self.formula_score_dict.items() if s > self.kwargs["pred_reject_formula_if_score_lower_than"]}
            self.formula_score_dict = dict(sorted(self.formula_score_dict.items(), key=lambda x: x[1], reverse=True))

        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICITIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>30} \t {:>20}".format(str(item[0]), str(round(item[1], 2))))


        if self.kwargs["pred_save_xic_plot"] == True:
            self.make_op_log_entry("INFO:\t" + "Start saving XIC plot...")
            self.xic_plot_filepath = self.kwargs["pred_folder"] + "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["mass_deviation"]) + "_" + str( self.spec.filter_mode) + ".png"
            self.save_xic_plot(self.xic[0], self.xic[1], self.xic_plot_filepath)
            self.make_op_log_entry("INFO:\t" + "Finished saving XIC plot...")
        if self.kwargs["pred_save_matplotlib_plot_of_isotopologues"] == True and len(list(self.formula_score_dict.keys())) > 0:
            self.make_op_log_entry("INFO:\t" + "Start saving matplotlib plot of isotopologues...")
            self.matplotlib_plot_filepath = self.kwargs["pred_folder"] + "isotopo_matplotlib_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.formula_score_dict.keys())[0]).items()])) + ".png"
            self.make_plot_of_isotopologues(ast.literal_eval(list(self.formula_score_dict.keys())[0]), self.matplotlib_plot_filepath)
            self.make_op_log_entry("INFO:\t" + "Finished saving matplotlib plot of isotopologues...")
        if self.kwargs["pred_save_go_plot_of_isotopologues"] == True and len(list(self.formula_score_dict.keys())) > 0:
            self.make_op_log_entry("INFO:\t" + "Start saving go plot of isotopologues...")
            self.go_plot_filepath = self.kwargs["pred_folder"] + "isotopo_go_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.formula_score_dict.keys())[0]).items()])) + ".html"
            self.make_go_plot_of_isotopologues(ast.literal_eval(list(self.formula_score_dict.keys())[0]), self.go_plot_filepath)
            self.make_op_log_entry("INFO:\t" + "Finished saving go plot of isotopologues...")
        try:
            if len(list(self.formula_score_dict.keys())) >= 1:
                self.make_op_log_entry("INFO:\t" + "Setting best formula prediction...")
                self.best_formula_prediction = ast.literal_eval(list(self.formula_score_dict.keys())[0])
                self.score_of_best_formula = self.formula_score_dict[(list(self.formula_score_dict.keys())[0])]
                self.simulated_isotopologue_pattern_for_best_formula = MS_functions.simulate_isotope_pattern_of_formula(self.best_formula_prediction, mass_resolution_ppm=self.kwargs["mass_deviation"])
                #self.simulated_isotopologue_pattern_for_best_formula = {mass: abundance, ...}
                self.make_op_log_entry("INFO:\t" + "Best formula prediction set. " + str(self.best_formula_prediction) + "   " + str(self.score_of_best_formula))
            else:
                self.make_op_log_entry("INFO:\t" + "No best formula was found!")
        except Exception as e:
            print("ERROR: Problem with best formula prediction setting!")
            print(str(e))
            print(self.formula_score_dict)
            print(traceback.format_exc())
        self.make_op_log_entry("INFO:\t" + "Finished creating prediction object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.spec.index))

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
                #X/C ratio check from #https://pmc.ncbi.nlm.nih.gov/articles/PMC1851972/ The numbers are representing 99.7% of all available molecular formulas.
                if "C" in list(f_dict.keys()):
                    if (f_dict.get("H", 0)/f_dict.get("C", 0) < 0.2) or (f_dict.get("H", 0)/f_dict.get("C", 0) > 3.1):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("F", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("F", 0)/f_dict.get("C", 0) > 1.5):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("Cl", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("Cl", 0)/f_dict.get("C", 0) > 0.8):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("Br", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("Br", 0)/f_dict.get("C", 0) > 0.8):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("N", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("N", 0)/f_dict.get("C", 0) > 1.3):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("O", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("O", 0)/f_dict.get("C", 0) > 3):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("P", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("P", 0)/f_dict.get("C", 0) > 0.3):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 4000
                    if (f_dict.get("S", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("S", 0)/f_dict.get("C", 0) > 0.8):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000
                    if (f_dict.get("Si", 0)/f_dict.get("C", 0) < 0) or (f_dict.get("Si", 0)/f_dict.get("C", 0) > 0.5):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - 1000



                heteroatom_count = float(f_dict.get("O", 0)) + float(f_dict.get("N", 0)) + float(f_dict.get("S", 0)) + float(f_dict.get("P", 0))
                if heteroatom_count > 5 and f_dict.get("C", 0) <= int(heteroatom_count/4):
                    formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - float((heteroatom_count/4) - f_dict.get("C", 0))*50
                if ("C" in list(f_dict.keys())) and ("N" in list(f_dict.keys())):
                    if int(f_dict["N"]) > 2 and (int(f_dict["C"]) / int(f_dict["N"]) <= 4):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - ((float(f_dict["N"]) - 2) ** 2) * 50
                if ("C" in list(f_dict.keys())) and ("H" in list(f_dict.keys())):
                    if (float(f_dict["H"]) / float(f_dict["C"]) >= 2):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - (((float(f_dict["H"]) / float(f_dict["C"])) - 2) ** 3) * 50
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - ((float(f_dict.get("N", 0))**2) * 50)
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - ((float(f_dict.get("N", 0)) ** 2) * 50)
                if ("P" in list(f_dict.keys())):
                    if (float(f_dict.get("P", 0)) * 3.1) >= float(f_dict.get("O", 0)):
                        formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - (((float(f_dict.get("P", 0)) * 3) / float(f_dict.get("O", 1))) ** 4) * 250
                if 2 < f_dict.get("C", 0) < heteroatom_count:
                    formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - ((heteroatom_count - float(f_dict.get("C", 0)))**2) * 50
                dbe = MS_functions.calc_dbe(f_dict)

                if dbe < 0:
                    score_substract = (abs(dbe + 2) * 50) ** 3
                elif (dbe - f_dict.get("O", 0)) > 7:
                    score_substract = abs(dbe - f_dict.get("O", 0) - 7) * 50
                else:
                    score_substract = 0
                if self.kwargs["charge_of_measured_mass"] < 0:
                    dbe = dbe + (self.kwargs["charge_of_measured_mass"] * 0.5)
                    if dbe - 1 > (self.mass * (62/1000)): #Senior Rule
                        score_substract += self.kwargs["pred_formula_likelihood_substract_score_rdbe_higher_than_senior_rule"]
                    if not dbe - int(dbe) == 0:
                        score_substract += self.kwargs["pred_formula_likelihood_substract_score_rdbe_non_integer"]

                #https://pmc.ncbi.nlm.nih.gov/articles/PMC1851972/
                #Element ratios

                formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) - score_substract
            except Exception as e:
                print("ERROR in formula score likelyhood: " + str(e))    
        self.make_op_log_entry("Finished adding likelyhood of formula to the score...")
        return formula_score_dict

    def get_peak_properties(self, times, intensities):
        self.make_op_log_entry("INFO:\t" + "Getting peak properties...")
        window, order = 10, 3
        intensities = scipy.signal.savgol_filter(intensities, window, order, mode="nearest")
        print(len(intensities))
        self.make_op_log_entry("INFO:\t" + "Intensities smoothed with Savitzky-Golay filter.")
        peak_properties = scipy.signal.find_peaks(intensities, height=max(intensities)/100, distance=2, prominence=max(intensities)/100, width=(2, 20))
        identified_peaks = []
        self.make_op_log_entry("INFO:\t" + "Peaks identified.")
        for element in range(len(peak_properties[1]["peak_heights"])):
            one_peak = []
            one_peak.append(times[int(peak_properties[1]["left_ips"][element] + (peak_properties[1]["widths"][element]/2))])
            one_peak.append(peak_properties[1]["peak_heights"][element])
            one_peak.append(times[int(peak_properties[1]["left_ips"][element])])
            one_peak.append(times[int(peak_properties[1]["right_ips"][element])])
            one_peak.append(peak_properties[1]["prominences"][element])
            identified_peaks.append(one_peak)
        self.make_op_log_entry("INFO:\t" + "Peak properties extracted.")
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

    def get_surrounding_spectra(self):
        index_before = None
        for i in range(self.spec.index - 1, -1, -1):
            if self.ms_file.all_modes[i] == self.spec.filter_mode:
                index_before = i
                break
        index_after = None
        for i in range(self.spec.index + 1, len(self.ms_file.all_modes), 1):
            if self.ms_file.all_modes[i] == self.spec.filter_mode:
                index_after = i
                break
        if index_before == None:
            print("Getting spectra before and after:")
            print(self.spec.index)
            print(self.ms_file.all_modes[self.spec.index])
            print(self.spec.filter_mode)
            print(self.ms_file.all_modes)
        print("getting surrounding spectra:")
        print("Start index = " + str(self.spec.index))
        print("before_index = " + str(index_before))
        print("after_index = " + str(index_after))
        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue
        self.spec_before = Spec(self.ms_file,
                                index_before,
                                spec_requested_filter_mode=self.spec.filter_mode,
                                absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/",
                                mass_deviation=self.kwargs["mass_deviation"],
                                spec_save_matplotlib_plot=False,
                                spec_save_go_plot=False,
                                **additional_kwargs)
        self.spec_after = Spec(self.ms_file,
                               index_after,
                               spec_requested_filter_mode=self.spec.filter_mode,
                               absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/",
                               mass_deviation=self.kwargs["mass_deviation"],
                               spec_save_matplotlib_plot=False,
                               spec_save_go_plot=False,
                               **additional_kwargs)
        return self.spec_before, self.spec_after

    def get_intensity_of_mass_in_spec(self, mass, spec, initial_deviation):
        within_deviation_mass_list = [masse for masse in list(spec.summarized_mass_intensity_dict.keys()) if abs((((mass - masse) / mass) * 1000000) - initial_deviation) < self.kwargs["pred_max_ppm_deviation_change_for_isotopologue"]]
        measured_intensity = 0
        for masse in within_deviation_mass_list:
            measured_intensity = measured_intensity + spec.summarized_mass_intensity_dict[masse]
        return measured_intensity

    def check_for_isotope_pattern(self, formula_dict):

        ordered_mass_intensity_dict = {self.spec.masses[i]: self.spec.intensities[i] for i in range(len(self.spec.masses))}
        ordered_mass_intensity_dict = dict(sorted(ordered_mass_intensity_dict.items(), key=lambda x: x[1], reverse=True))
        ordered_mass_intensity_dict = {key: value for key, value in ordered_mass_intensity_dict.items() if value > 1}
        
        isotope_pattern_dict = MS_functions.simulate_isotope_pattern_of_formula(formula_dict)
        isotope_pattern_dict = dict(sorted(isotope_pattern_dict.items(), key=lambda x: x[1], reverse=True))
        self.make_op_log_entry("Simulated isotope pattern (not corrected for charge and electron mass (equals neutral charge)): " + str(isotope_pattern_dict))
        isotope_pattern_dict = {abs((m - (0.000548 * self.kwargs["charge_of_measured_mass"])) / self.kwargs["charge_of_measured_mass"]): a for m, a in isotope_pattern_dict.items()}
        self.make_op_log_entry("Simulated isotope pattern (already corrected for electron mass): " + str(isotope_pattern_dict))

        isotope_pattern_deviation_list = [abs(self.mass - m) for m in list(isotope_pattern_dict.keys())]
        if isotope_pattern_deviation_list.index(min(isotope_pattern_deviation_list)) == 0:
            pass
        else:
            self.make_op_log_entry("The mass of the given ion is not the same as the most intense peak for this specific molecule. Adapting the isotope_pattern_dict... (swapping first and second isotope)")
            temp_isotope_pattern_dict_list = list(zip(    list(isotope_pattern_dict.keys()), list(isotope_pattern_dict.values())    ))
            if len(temp_isotope_pattern_dict_list) >= 2:
                temp_isotope_pattern_dict_list[1], temp_isotope_pattern_dict_list[0] = temp_isotope_pattern_dict_list[0], temp_isotope_pattern_dict_list[1]
                isotope_pattern_dict = dict(temp_isotope_pattern_dict_list)
            else:
                pass
            self.make_op_log_entry("Adapted isotope_pattern_dict: " + str(isotope_pattern_dict))

        isotopes_found = {}
        previous_deviation_list = []
        startitem = list(isotope_pattern_dict.items())[0]
        isotopo_mass = startitem[0]
        measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(isotopo_mass - x))
        initial_deviation = ((isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000
        print("Initial deviation for isotope check: " + str(initial_deviation))

        if (not self.spec_before is None) and (not self.spec_after is None):
            self.intensity_of_ion = 0
            self.intensity_of_ion += self.get_intensity_of_mass_in_spec(isotopo_mass, spec=self.spec, initial_deviation=initial_deviation)
            self.intensity_of_ion += self.get_intensity_of_mass_in_spec(isotopo_mass, spec=self.spec_before, initial_deviation=initial_deviation)
            self.intensity_of_ion += self.get_intensity_of_mass_in_spec(isotopo_mass, spec=self.spec_after, initial_deviation=initial_deviation)
        else:
            self.intensity_of_ion = sum([self.spec.intensities[i] for i in range(len(self.spec.masses)) if abs(((self.spec.masses[i] - self.mass) / self.mass) * 1000000) <= self.kwargs["mass_deviation"]])

        for isotopo in list(isotope_pattern_dict.items()):
            break_the_isotopo_prediction = False
            isotopo_mass = isotopo[0]
            within_deviation_mass_list = [masse for masse in list(self.spec.summarized_mass_intensity_dict.keys()) if abs((((isotopo_mass - masse) / isotopo_mass) * 1000000) - initial_deviation) < self.kwargs["pred_max_ppm_deviation_change_for_isotopologue"]]
            measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
            deviation = ((isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000
            theoretical_intensity = (self.intensity_of_ion / (list(isotope_pattern_dict.items())[0][1])) * isotopo[1]
            if theoretical_intensity <= 0.01:
                theoretical_intensity = 0.01
            noise = sorted(list(self.spec.summarized_mass_intensity_dict.values()))[int(len(list(self.spec.summarized_mass_intensity_dict.items())) / self.kwargs["pred_noise_divisor_for_isotopologue_calculation"])] + self.kwargs["pred_minimum_assumed_noise"]

            if (not self.spec_before is None) and (not self.spec_after is None):
                measured_intensity_with_surrounding_spectra = 0
                measured_intensity_with_surrounding_spectra += self.get_intensity_of_mass_in_spec(isotopo_mass, spec=self.spec, initial_deviation=initial_deviation)
                measured_intensity_with_surrounding_spectra += self.get_intensity_of_mass_in_spec(isotopo_mass, spec=self.spec_before, initial_deviation=initial_deviation)
                measured_intensity_with_surrounding_spectra += self.get_intensity_of_mass_in_spec(isotopo_mass, spec=self.spec_after, initial_deviation=initial_deviation)
                #print("Old ratio: " + str(intensity_of_ion_with_surrounding_spectra/measured_intensity_with_surrounding_spectra) + "  New ratio: " + str(self.intensity_of_ion/measured_intensity))
                measured_intensity = measured_intensity_with_surrounding_spectra
            else:
                measured_intensity = 0
                for masse in within_deviation_mass_list:
                    measured_intensity = measured_intensity + self.spec.summarized_mass_intensity_dict[masse]

            print("Measured intensity for ion: " + str(measured_intensity))

            if (deviation < self.kwargs["mass_deviation"]) and (abs(deviation - initial_deviation) <= self.kwargs["pred_max_ppm_deviation_change_for_isotopologue"]) and (abs(measured_mass_with_minimal_deviation-isotopo_mass) <= (self.kwargs["mass_deviation"]/1000000)*150 ):
                print("In if statement for isotope check...")
                previous_deviation_list.append(copy.deepcopy(deviation))
                if measured_intensity <= 1:
                    print("Breaking the isotope check, as measured intensity was <= 1")
                    break_the_isotopo_prediction = True
                    measured_intensity = 0.01
                isotopes_found[isotopo[0]] = [True, theoretical_intensity, measured_intensity]
                score = measured_intensity / theoretical_intensity  # je naeher an 1 desto besser; wenn <1: weniger gemessen als da sein sollte; wenn >1: mehr gemessen als da sein sollte.
                if score > 1:
                    score = -(-1.6 + (1 / (0.65 + (2.718281828459045 ** (-self.kwargs["pred_isotopologue_score_e_function_exponent"] * (score-1 ))))))
                # score: je naeher an 1 desto besser; wenn negativ: weniger gemessen als theoretisch da; wenn positiv: mehr gemessen als theoretisch da.
                if score < 1:
                    score = score ** (0.7/(2-self.kwargs["pred_isotopologue_score_e_function_exponent"]))
                    #score = score - (score * 0.15)
                    #score = score ** (1/(1-self.kwargs["pred_isotopologue_score_e_function_exponent"]))

                score = score.real

                score = score - 0.5

                score = ((score * 100) + 0.001)
                print("Score1: " + str(score))

                #if theoretical_intensity < noise:
                #    try:
                #        second_highest_abundance = (list(isotope_pattern_dict.items())[1][1])
                #        curr_abundance = isotopo[1]
                #    except:
                #        second_highest_abundance = 1
                #        curr_abundance = 1
                #    score_multiplier_by_intensity_and_noise = abs(1 / (math.log(theoretical_intensity / noise) - 1)) * ((1 - math.exp(-((curr_abundance / second_highest_abundance) * 6)) )+0.2)
                #else:
                #    try:
                #        second_highest_abundance = (list(isotope_pattern_dict.items())[1][1])
                #        curr_abundance = isotopo[1]
                #    except:
                #        second_highest_abundance = 1
                #        curr_abundance = 1
                #    #score_multiplier_by_intensity_and_noise = abs(math.log(measured_intensity / noise) + 1) * ((1 - math.exp(-((curr_abundance / second_highest_abundance) * 6)) )+0.2)
                #    score_multiplier_by_intensity_and_noise = abs(math.log(theoretical_intensity / noise) + 1) * ((1 - math.exp(-((curr_abundance / second_highest_abundance) * 6))) + 0.2)

                score_multiplier_by_intensity_and_noise = theoretical_intensity / (theoretical_intensity + noise)

                if isotopo == list(isotope_pattern_dict.items())[0]:
                    isotopes_found[isotopo[0]].append(0)
                    isotopes_found[isotopo[0]].append(measured_mass_with_minimal_deviation)
                    previous_score = 9999999
                    continue

                score = score * score_multiplier_by_intensity_and_noise

                print("Score2: " + str(score))
                print("score multiplier: " + str(score_multiplier_by_intensity_and_noise))

                if measured_intensity <= 10:
                    pass
                else:
                    score = score + self.kwargs["pred_add_value_to_score_if_isotopo_was_found"]
                print("Score3: " + str(score))

                if theoretical_intensity <= (noise / 1.5):
                    ppm_influence = ((2 ** ((self.kwargs["pred_isotopologue_score_e_function_exponent"] * 2) * (abs(deviation) - (self.kwargs["mass_deviation"] / 3)))) * 2)
                    ppm_influence = abs(ppm_influence * self.kwargs["pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    if ppm_influence > 30:
                        ppm_influence = 30
                    score = (score) - ppm_influence
                    print("Score4: " + str(score))
                    isotopes_found[isotopo[0]].append(score)
                    if (score - abs(self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
                        isotopes_found[isotopo[0]][-1] = isotopes_found[isotopo[0]][-1] * 0.1
                    isotopes_found[isotopo[0]].append(measured_mass_with_minimal_deviation)
                    break
                else:
                    ppm_influence = ((2 ** ((self.kwargs["pred_isotopologue_score_e_function_exponent"] * 2) * (abs(deviation) - (self.kwargs["mass_deviation"] / 3)))) * 2)
                    ppm_influence = abs(ppm_influence * self.kwargs["pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    if ppm_influence > 50:
                        ppm_influence = 50
                    score = (score) - ppm_influence
                    print("Score4: " + str(score))
                    isotopes_found[isotopo[0]].append(score)
                    isotopes_found[isotopo[0]].append(measured_mass_with_minimal_deviation)
                
                if (score < self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]):
                    break

                if (score - abs(self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
                    isotopes_found[isotopo[0]][-1] = isotopes_found[isotopo[0]][-1] * 0.1
                    break
                previous_score = score
                if break_the_isotopo_prediction:
                    break
            else:
                isotopes_found[isotopo[0]] = [False, theoretical_intensity, measured_intensity]
                negative_score = 50 * (theoretical_intensity / (theoretical_intensity + noise))
                try:
                    second_highest_abundance = (list(isotope_pattern_dict.items())[1][1])
                    curr_abundance = isotopo[1]
                except:
                    second_highest_abundance = 1
                    curr_abundance = 1
                negative_score = negative_score * ((1 - math.exp(-((curr_abundance / second_highest_abundance) * 6))) + 0.2)
                negative_score = -1 * negative_score
                isotopes_found[isotopo[0]].append(negative_score)
                isotopes_found[isotopo[0]].append(measured_mass_with_minimal_deviation)
                break
        self.intensity_of_ion = sum([self.spec.intensities[i] for i in range(len(self.spec.masses)) if abs(((self.spec.masses[i] - self.mass) / self.mass) * 1000000) <= self.kwargs["mass_deviation"]])
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
        fig = matplotlib.figure.Figure()
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
        for a in ax:
            a.legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        fig.savefig(filepath, bbox_inches='tight', dpi=DPI)
        fig.clf()
        fig.clear()
        return True

    def get_peak_matching_of_isotopologues(self, formula_score_dict):
        # formula score dict given in the form of {"{'C': 2, 'H': 4, 'O': 1}": 200, "{'C': 3, 'H': 6, 'O': 1}": 100, ...}
        noise = sorted(list(self.spec.summarized_mass_intensity_dict.values()))[int(len(list(self.spec.summarized_mass_intensity_dict.items())) / self.kwargs["pred_noise_divisor_for_isotopologue_calculation"])] + self.kwargs["pred_minimum_assumed_noise"]
        intensity_of_mi = self.intensity_of_ion

        new_formula_score_dict = copy.deepcopy(formula_score_dict)

        mi_xic = self.xic
        self.make_op_log_entry("Starting to calculate the matching of the areas for every prediction in the formula_score_dict.")
        self.make_op_log_entry("Good matching mass traces for isotopologues will increase the score. Bad matching will lead to a decreasing score.")
        for formula, score in formula_score_dict.items():
            print("Starting isotopo matching for formula: " + str(formula))
            formula_dict = ast.literal_eval(formula)
            simulated_isotopo_abundances_dict = MS_functions.simulate_isotope_pattern_of_formula(formula_dict)

            startitem = list(simulated_isotopo_abundances_dict.items())[0]
            isotopo_mass = startitem[0]
            measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(isotopo_mass - x))
            initial_deviation = ((isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000
            old_area = 0
            change_score = 0
            for isotopo_masse, isotopo_abundance in simulated_isotopo_abundances_dict.items():
                if list(simulated_isotopo_abundances_dict.keys()).index(isotopo_masse) == 0:
                    continue
                mass_deviation_for_xic = ((self.kwargs["mass_deviation"] * isotopo_masse)/1000000)
                isotopo_xic = MS_functions.get_xic(self.ms_file.rawdata, isotopo_masse, mass_deviation_for_xic, requested_filter_mode=self.spec.filter_mode)
                peak1_rt, peakintensity1, peak2_rt, peakintensity2 = self.isolate_and_prepare_peak(mi_xic, isotopo_xic, self.rt)
                break_formula_evaluation = False
                try:
                    P = np.array([peak1_rt, peakintensity1]).T
                    Q = np.array([peak2_rt, peakintensity2]).T
                    area = similaritymeasures.area_between_two_curves(P, Q)
                    if not (isinstance(area, float) or isinstance(area, int)):
                        try:
                            area = float(area)
                        except:
                            area = self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"] * 30
                    if not 0 < area < self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"] * 30:
                        print("Area was not with good value! Setting area to: self.kwargs['oa_fragments_do_peak_computation_if_area_higher_than'] * 30")
                        area = self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"] * 30
                    if area == np.nan or (str(area).lower() == "nan"):
                        print("area was nan. Chaning area..." + str(area))
                        area = self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"] * 30
                        break_formula_evaluation = True
                except:
                    continue
                print("Area for isotopologue with mass: " + str(isotopo_masse) + "   ; area = " + str(area))
                theoretical_intensity_of_isotopo_peak = (intensity_of_mi / (list(simulated_isotopo_abundances_dict.items())[0][1])) * isotopo_abundance
                within_deviation_mass_list = [masse for masse in list(self.spec.summarized_mass_intensity_dict.keys()) if abs((((isotopo_masse - masse) / isotopo_masse) * 1000000) - initial_deviation) < self.kwargs["pred_max_ppm_deviation_change_for_isotopologue"]]
                measured_intensity = 0
                for masse in within_deviation_mass_list:
                    measured_intensity = measured_intensity + self.spec.summarized_mass_intensity_dict[masse]
                isotopo_int_percent_of_mi_int = theoretical_intensity_of_isotopo_peak / intensity_of_mi
                area_threshhold = self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"]
                if area < area_threshhold: #the peaks match good
                    area = area + 0.03
                    delta_score = 1 / (1-((area_threshhold-area)/area_threshhold)) # the closer to 1 the better is the matching of the peaks.
                    delta_score = (delta_score ** 0.5)
                    delta_score = delta_score.real
                else: #the peaks do not match good
                    delta_score = (1 + (area-area_threshhold))
                    delta_score = (delta_score ** 0.5).real
                    delta_score = -1 * delta_score
                rel_isotopo_abundance = (isotopo_abundance / (list(simulated_isotopo_abundances_dict.items())[0][1] ))
                if delta_score < 0:
                    delta_score = delta_score * 1.3
                    break_formula_evaluation = True
                if break_formula_evaluation:
                    rel_isotopo_abundance = (sum(list(simulated_isotopo_abundances_dict.values())[list(simulated_isotopo_abundances_dict.values()).index(isotopo_abundance):]) / (list(simulated_isotopo_abundances_dict.items())[0][1] ))
                isotopo_xic_matching_multiplier = 1
                change_score = change_score + (delta_score * (score*rel_isotopo_abundance) * isotopo_xic_matching_multiplier)
                if rel_isotopo_abundance <= 0.001 or break_formula_evaluation:
                    break
            print("Change score for formula: " + str(formula) + "    ; Change score: " + str(change_score) + "    ; New score = " + str(score + change_score))
            new_score = score + change_score
            new_formula_score_dict[formula] = new_score
            self.make_op_log_entry("Change score for formula: " + str(formula) + "    ; Change score: " + str(change_score) + "    ; New score = " + str(score + change_score))
        return new_formula_score_dict

    def isolate_and_prepare_peak(self, xic1, xic2, peak_rt, peakwidth=10):
        index = xic1[0].index(min(xic1[0], key=lambda x: abs(peak_rt - x)))
        intensity1_at_peak_rt = xic1[1][index]
        if intensity1_at_peak_rt <= 0:
            intensity1_at_peak_rt = 0.0001
        intensity2_at_peak_rt = xic2[1][index]
        if intensity2_at_peak_rt <= 0:
            intensity2_at_peak_rt = 0.0001
        try:
            peakintensity1 = xic1[1][(index - peakwidth):(index + peakwidth)]
            peak1_rt = xic1[0][int(index - peakwidth):int(index + peakwidth)]
            peakintensity2 = xic2[1][int(index - peakwidth):int(index + peakwidth)]
            peak2_rt = xic2[0][int(index - peakwidth):int(index + peakwidth)]

            neighbour_list1 = xic1[1][int(index - 3 * peakwidth):int(index + 3 * peakwidth)]
            neighbour_list1 = [neighbour_list1[i] for i in range(len(neighbour_list1)) if (i < len(neighbour_list1) / 3) or (i > ((len(neighbour_list1) / 3) + (len(neighbour_list1) / 2)))]
            average_surrounding1 = sum(neighbour_list1) / len(neighbour_list1)
            try:
                peakintensity1 = [((i-average_surrounding1) / (xic1[1][index]-average_surrounding1)) for i in peakintensity1]
            except:
                max_peakint1 = max(peakintensity1)
                if max_peakint1 == 0:
                    max_peakint1 = 0.001
                peakintensity1 = [((i-average_surrounding1) / (max_peakint1)) for i in peakintensity1]


            neighbour_list2 = xic2[1][int(index - 3*peakwidth):int(index + 3*peakwidth)]
            neighbour_list2 = [neighbour_list2[i] for i in range(len(neighbour_list2)) if (i < len(neighbour_list2)/3) or (i > ((len(neighbour_list2)/3) + (len(neighbour_list2)/2))) ]
            average_surrounding2 = sum(neighbour_list2) / len(neighbour_list2)
            try:
                peakintensity2 = [((i-average_surrounding2) / (xic2[1][index]-average_surrounding2)) for i in peakintensity2]
            except:
                max_peakint2 = max(peakintensity2)
                if max_peakint2 == 0:
                    max_peakint2 = 0.0001
                peakintensity2 = [[((i-average_surrounding2) / (max_peakint2)) for i in peakintensity2]]
        except:
            if index - peakwidth <= 2:
                peakwidth = index - 2
            if index + peakwidth >= len(xic1[1]):
                peakwidth = (len(xic1[1]) - index - 2)
            peakintensity1 = xic1[1][(index - peakwidth):(index + peakwidth)]
            peakintensity1 = [i/intensity1_at_peak_rt for i in peakintensity1]
            peak1_rt = xic1[0][int(index - peakwidth):int(index + peakwidth)]
            peakintensity2 = xic2[1][int(index - peakwidth):int(index + peakwidth)]
            peakintensity2 = [i / intensity2_at_peak_rt for i in peakintensity2]
            peak2_rt = xic2[0][int(index - peakwidth):int(index + peakwidth)]

        return peak1_rt, peakintensity1, peak2_rt, peakintensity2

    def save_xic_plot(self, times, intensities, xic_plot_filepath):
        fig = matplotlib.figure.Figure()
        ax = fig.subplots()
        if self.identified_peaks_for_mass is None:
            self.identified_peaks_for_mass = self.get_peak_properties(times, intensities)
        identified_peak_times = [entry[0] for entry in self.identified_peaks_for_mass]
        plot_heights = [max(intensities)/4 for i in range(len(identified_peak_times))]
        ax.scatter(identified_peak_times, plot_heights, color="green", label="identified peak", s=20, alpha=0.5)
        ax.bar(self.spec.rt, max(intensities), color="red", label="observed time", alpha=0.5, width=10)
        ax.plot(times, intensities, color="blue", label="XIC")
        ax.set_xlabel("retention time / s")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("Extracted Ion Chromatogram (XIC) for mass: " + str(round(self.mass, 4)) + " at RT: " + str(round(self.rt, 2)) + " s. Filtermode: " + str(self.spec.filter_mode) , wrap=True)
        ax.legend()
        fig.savefig(xic_plot_filepath, bbox_inches='tight', dpi=DPI)
        fig.clf()
        fig.clear()
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

                          "oa_fragments_absolute_max_number_of_fragment_masses": 100,
                          "oa_fragments_pred_save_matplotlib_plot_of_isotopologues":True,
                          "oa_fragments_pred_save_go_plot_of_isotopologues":True,
                          "oa_fragments_spec_save_matplotlib_plot":True,
                          "oa_fragments_spec_save_go_plot":False,
                          "oa_fragments_pred_save_xic_plot":True,
                          "oa_fragments_pred_save_detailed_log":True,
                          "oa_fragments_only_calc_prediction_if_peak_is_found":True,
                          "oa_fragments_do_good_peak_comparison_with_area_between_curves": True,
                          "oa_fragments_do_peak_computation_if_area_higher_than": 3,
                          "oa_nl_deviation_threshold_influence_on_abs_score": 6,

                          "oa_include_frag_intensity_noise_multiplier":0.01,
                          "oa_reject_formula_if_score_lower_than": 10,
                          "oa_make_good_fragment_formula_prediction":False,
                          "pred_formula_cache_folder_path":"U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//",
                          "oa_zip_folder_when_finished": True}
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["one_analysis_folder"], exist_ok=True)

        print()
        print("OneAnalysis kwargs: " + str(self.kwargs))
        print()
        print("Make good fragment prediction: " + str(self.kwargs["oa_make_good_fragment_formula_prediction"]))
        self.kwargs["oa_make_good_fragment_formula_prediction"] = bool(self.kwargs["oa_make_good_fragment_formula_prediction"])
        print("Make good fragment prediction: " + str(self.kwargs["oa_make_good_fragment_formula_prediction"]))

        self.make_oa_log_entry("")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("INFO:\t" + "Creating OneAnalysis object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.peak_index))
        self.make_oa_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))
        self.make_oa_log_entry("INFO:\t" + "Assuming a charge of the measured mass of: " + str(self.kwargs["charge_of_measured_mass"]))
        self.make_oa_log_entry("INFO:\t" + "Available MS modes: " + str(self.ms_file.available_modes))

        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, self.kwargs["oa_mass_deviation_xic"], requested_filter_mode=self.kwargs["oa_xic_requested_filter_mode"])
        self.make_oa_log_entry("INFO:\t" + "XIC calculated. Continuing...")
        if self.kwargs["oa_save_xic_plot"] == True:
            print("Saving XIC plot...")
            self.xic_plot_filepath = self.kwargs["one_analysis_folder"] + "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["oa_mass_deviation_xic"]) + "_" + str(
                                    self.kwargs["oa_xic_requested_filter_mode"]) + ".png"
            self.plot_and_save_xic(self.xic_plot_filepath, self.peak_index, self.xic)
        
        self.make_oa_log_entry("INFO:\t" + "Extracted Ion Chromatogram (XIC) plot saved at: " + str(self.xic_plot_filepath))

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
            self.mi_spec_before, self.mi_spec_after = self.get_surrounding_best_mi_spectra()
            self.best_frag_spec = self.full_scan_spec
            self.frag_spec_before, self.frag_spec_after = self.get_surrounding_best_fragment_spectra()
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
                self.mi_spec_before, self.mi_spec_after = self.get_surrounding_best_mi_spectra()
            self.frag_spec_before, self.frag_spec_after = self.get_surrounding_best_fragment_spectra()
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
        
        self.make_oa_log_entry("INFO:\t" + "Finished searching for spectra...")

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

        self.make_oa_log_entry("INFO:\t" + "Best molecular ion spectrum at index: " + str(self.best_molecular_ion_spec.index))
        self.make_oa_log_entry("INFO:\t" + "Best fragment ion spectrum at index: " + str(self.best_frag_spec.index))
        
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
            self.make_oa_log_entry("INFO:\t" + "No molecular ion prediction peak found. Stopping prediction...")
            return

        if len(list(self.summarized_molecular_ion_formula_score_dict.keys())) == 0:
            print("No molecular ion prediction could be found! Returning....")
            self.make_oa_log_entry("INFO:\t" + "No molecular ion prediction could be found! Returning.....")
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

        self.true_fragment_list1, self.matching_fragments_list = self.make_true_fragment_list()
        if len(self.matching_fragments_list) == 0:
            try:
                mi_formula_prediction = self.best_molecular_ion_prediction
            except:
                mi_formula_prediction = "None"
            self.matching_fragments_list.append([self.mass, self.score_of_best_molecular_ion_prediction, self.intensity_of_molecular_ion, mi_formula_prediction, 0, 0, 0, "None", "None", 0])
        print("Summary list of one analysis: " + str(self.true_fragment_list1))
        self.make_oa_log_entry("INFO:\t" + "Finished creating matching fragments list...")
        self.make_oa_log_entry("INFO:\t" + "All matching fragments: " + str(self.matching_fragments_list))

        self.true_fragment_list = self.process_matching_fragments_list(self.matching_fragments_list)
        print("Summary list of one analysis: " + str(self.true_fragment_list))
        
        self.make_oa_log_entry("INFO:\t" + "Finished creating summary list of one analysis...")
        self.make_oa_log_entry("INFO:\t" + "Summary list 1 of one analysis: " + str(self.true_fragment_list1))
        self.make_oa_log_entry("INFO:\t" + "Summary list 2 of one analysis: " + str(self.true_fragment_list))

        write_to_summary_file_list = [self.rt]
        write_to_summary_file_list.extend(self.true_fragment_list)
        
        self.append_oa_summary_to_raw_file_summary(write_to_summary_file_list)

        self.make_oa_log_entry("INFO:\t" + "Finished appending summary to raw file summary...")
        self.make_oa_log_entry("INFO:\t" + "Creating summary plot and txt file...")
        save_filepath = self.kwargs["one_analysis_folder"] + "oa_summary_plot.png"
        try:
            self.create_summary_plot(self.true_fragment_list, self.best_frag_spec, self.fragment_predictions, save_filepath)
        except Exception as e:
            self.make_oa_log_entry("ERROR:\t" + "Error while creating summary plot: " + str(e))
            print("Error while creating summary plot: " + str(e))
            print(traceback.format_exc())
        self.make_oa_log_entry("INFO:\t" + "Summary plot saved at: " + str(save_filepath))
        oa_txt_save_filepath = self.kwargs["one_analysis_folder"] + "BEST_FORMULA_PREDICTION.txt"
        self.create_oa_summary_txtfile(self.true_fragment_list, self.fragment_predictions, oa_txt_save_filepath)
        oa_plt_save_filepath = self.kwargs["one_analysis_folder"] + "FRAGMENT_SUMMARY.png"
        self.plot_summary_xic_of_mi_and_fragments(oa_plt_save_filepath)

        if self.kwargs["oa_zip_folder_when_finished"] == True:
            shutil.make_archive(self.ms_file.parentfolder + "/" + str(self.mass) + "_" + str(self.rt), "zip", self.kwargs["one_analysis_folder"])
            shutil.rmtree(self.kwargs["one_analysis_folder"])



    def get_surrounding_best_fragment_spectra(self):
        index_before = None
        for i in range(self.best_frag_spec.index - 1, -1, -1):
            if self.ms_file.all_modes[i] == self.best_frag_spec.filter_mode:
                index_before = i
                break
        index_after = None
        for i in range(self.best_frag_spec.index + 1, len(self.ms_file.all_modes), 1):
            if self.ms_file.all_modes[i] == self.best_frag_spec.filter_mode:
                index_after = i
                break
        if index_before == None:
            print("Getting spectra before and after:")
            print(self.best_frag_spec.index)
            print(self.ms_file.all_modes[self.best_frag_spec.index])
            print(self.best_frag_spec.filter_mode)
            print(self.ms_file.all_modes)
        print("getting surrounding spectra:")
        print("Start index = " + str(self.best_frag_spec.index))
        print("before_index = " + str(index_before))
        print("after_index = " + str(index_after))
        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue
        self.spec_before = Spec(self.ms_file,
                                index_before,
                                spec_requested_filter_mode=self.best_frag_spec.filter_mode,
                                absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/",
                                mass_deviation=self.kwargs["mass_deviation"],
                                spec_save_matplotlib_plot=False,
                                spec_save_go_plot=False,
                                **additional_kwargs)
        self.spec_after = Spec(self.ms_file,
                               index_after,
                               spec_requested_filter_mode=self.best_frag_spec.filter_mode,
                               absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/",
                               mass_deviation=self.kwargs["mass_deviation"],
                               spec_save_matplotlib_plot=False,
                               spec_save_go_plot=False,
                               **additional_kwargs)
        return self.spec_before, self.spec_after

    def get_surrounding_best_mi_spectra(self):
        index_before = None
        for i in range(self.best_molecular_ion_spec.index - 1, -1, -1):
            if self.ms_file.all_modes[i] == self.best_molecular_ion_spec.filter_mode:
                index_before = i
                break
        index_after = None
        for i in range(self.best_molecular_ion_spec.index + 1, len(self.ms_file.all_modes), 1):
            if self.ms_file.all_modes[i] == self.best_molecular_ion_spec.filter_mode:
                index_after = i
                break
        if index_before == None:
            print("Getting spectra before and after:")
            print(self.best_molecular_ion_spec.index)
            print(self.ms_file.all_modes[self.best_molecular_ion_spec.index])
            print(self.best_molecular_ion_spec.filter_mode)
            print(self.ms_file.all_modes)
        print("getting surrounding spectra:")
        print("Start index = " + str(self.best_molecular_ion_spec.index))
        print("before_index = " + str(index_before))
        print("after_index = " + str(index_after))
        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue
        self.spec_before = Spec(self.ms_file,
                                index_before,
                                spec_requested_filter_mode=self.best_molecular_ion_spec.filter_mode,
                                absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/",
                                mass_deviation=self.kwargs["mass_deviation"],
                                spec_save_matplotlib_plot=False,
                                spec_save_go_plot=False,
                                **additional_kwargs)
        self.spec_after = Spec(self.ms_file,
                               index_after,
                               spec_requested_filter_mode=self.best_molecular_ion_spec.filter_mode,
                               absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/",
                               mass_deviation=self.kwargs["mass_deviation"],
                               spec_save_matplotlib_plot=False,
                               spec_save_go_plot=False,
                               **additional_kwargs)
        return self.spec_before, self.spec_after

    def process_matching_fragments_list(self, matching_fragments_list):
        #matching_fragments_list.append([molecular_ion_mass, mi_score, molecular_ion_intensity, mi_formula_str, f_mass, f_score, f_intensity, f_formula_str, nl_formula_str, nl_deviation])
        mi_predictions_combined_scores_dict = {}
        for entry in matching_fragments_list:
            mi_mass = entry[0]
            mi_score = entry[1]
            mi_intensity = entry[2]
            mi_formula_str = entry[3]
            f_mass = entry[4]
            f_score = entry[5]
            f_intensity = entry[6]
            f_formula_str = entry[7]
            nl_formula_str = entry[8]
            nl_deviation = entry[9]
            if not mi_formula_str in list(mi_predictions_combined_scores_dict.keys()):
                mi_predictions_combined_scores_dict[mi_formula_str] = ( mi_score * (1 - ( 1 / (1+math.exp( -0.5*(nl_deviation-self.kwargs["oa_nl_deviation_threshold_influence_on_abs_score"]) ) ) ) + 0.045 ) )
            else:
                mi_predictions_combined_scores_dict[mi_formula_str] = mi_predictions_combined_scores_dict[mi_formula_str] + ( mi_score * (1 - ( 1 / (1+math.exp( -0.5*(nl_deviation-self.kwargs["oa_nl_deviation_threshold_influence_on_abs_score"]) ) ) ) + 0.045) )
        best_formula_str = max(mi_predictions_combined_scores_dict, key=mi_predictions_combined_scores_dict.get)
        true_frag_list = []
        for entry in matching_fragments_list:
            mi_mass = entry[0]
            mi_score = entry[1]
            mi_intensity = entry[2]
            mi_formula_str = entry[3]
            f_mass = entry[4]
            f_score = entry[5]
            f_intensity = entry[6]
            f_formula_str = entry[7]
            nl_formula_str = entry[8]
            nl_deviation = entry[9]
            if mi_formula_str == best_formula_str:
                true_frag_list.append([mi_mass, mi_score, mi_intensity, mi_formula_str, f_mass, f_score, f_intensity, f_formula_str, nl_formula_str, nl_deviation])
        return true_frag_list

    def create_oa_summary_txtfile(self, true_fragment_list, fragment_predictions, save_filepath):
        #matching_fragments_list.append([molecular_ion_mass, mi_score, molecular_ion_intensity, mi_formula_str, f_mass, f_score, f_intensity, f_formula_str, nl_formula_str, nl_deviation])
        #in true_fragment list all molecular ion formulas are the same!
        with open(save_filepath, "w") as f:
            f.write("===========================================\n")
            f.write("BEST FORMULA APPROXIMATION\n")
            f.write("===========================================\n")
            f.write("Mass: \t" + str(round(self.molecular_ion_prediction.mass, 4)) + " u\n")
            f.write("Retention time: \t" + str(round(self.molecular_ion_prediction.rt, 2)) + " s\n")
            f.write("Formula: \t" + str(true_fragment_list[0][3]) + "\n")
            f.write("Combined Score: \t" + str(sum([frag[1] for frag in true_fragment_list])  ) + "\n")
            f.write("___________________________________________\n")
            f.write("Peak found: \t" + str(self.molecular_ion_prediction.peak_found) + "\n")
            f.write("Intensity: \t" + str(self.molecular_ion_prediction.intensity_of_ion) + "\n")
            f.write("Spec filter: \t" + str(self.molecular_ion_prediction.spec.filter) + "\n")
            f.write("\n")
            f.write("=================FRAGMENTS=================\n")
            f.write("{:<15}|{:<15}|{:<15}|{:<15}|{:<15}|{:<15}\n".format("Mass", "Formula", "Score", "Intensity", "NL Formula", "NL Deviation"))
            f.write("{:_<15}|{:_<15}|{:_<15}|{:_<15}|{:_<15}|{:_<15}\n".format("", "", "", "", "", ""))
            for frag in true_fragment_list:
                mass = frag[4]
                formula = frag[7]
                score = frag[5]
                intensity = frag[6]
                nl_formula = frag[8]
                nl_deviation = frag[9]
                f.write("{:<15}|{:<15}|{:<15}|{:<15}|{:<15}|{:<15}\n".format(str(round(mass, 4)), str(formula), str(round(score, 1)), str(round(intensity, 1)), str(nl_formula), str(round(nl_deviation, 1))))
            f.write("Created at: " + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    def create_summary_plot(self, true_fragment_list, best_frag_spec, fragment_predictions, save_filepath):
        #matching_fragments_list.append([molecular_ion_mass, mi_score, molecular_ion_intensity, mi_formula_str, f_mass, f_score, f_intensity, f_formula_str, nl_formula_str, nl_deviation])
        #in true_fragment list all molecular ion formulas are the same!
        print()
        print("Fragment Predictions for summary plot!")
        print(fragment_predictions)
        print()

        ncols = 2
        nrows = 2
        frag_masses = [frag[4] for frag in true_fragment_list if frag[4] >= 1]
        nrows = nrows + int( len(frag_masses)/2 + 0.5 )
        #fig, ax = plt.subplots(ncols=ncols, nrows=nrows, dpi=1000, layout="tight", figsize=(5*ncols, 5*nrows))
        fig = matplotlib.figure.Figure(constrained_layout=True, figsize=(5*ncols, 5*nrows), dpi=DPI)
        gs = matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, figure=fig)
        ax = [[None, None],
              [None      ]]
        ax[0][0] = fig.add_subplot(gs[0, 0])
        ax[0][1] = fig.add_subplot(gs[0, 1])
        ax[1][0] = fig.add_subplot(gs[1, :])
        ax[1][0].axis("off")
        frag_axs = []
        for i in range(0, len(frag_masses), 1):
            row_index = int(i/2 + 2)
            col_index = int(i%2)
            frag_axs.append(fig.add_subplot(gs[row_index, col_index]))

        

        fig.suptitle("Summary Plot for Mass: " + str(round(self.mass, 4)) + " at RT: " + str(round(self.rt, 2)), wrap=True, fontsize=16, fontweight="bold")

        ax[0][0].bar(best_frag_spec.summarized_masses, best_frag_spec.summarized_intensities, color="gray", label="Spectrum", alpha=0.2)
        intensity_of_molecular_ion_in_frag_spec = sum([best_frag_spec.summarized_intensities[i] for i in range(len(best_frag_spec.summarized_intensities)) if abs(((best_frag_spec.summarized_masses[i] - self.mass)/self.mass)*1000000) <= self.kwargs["mass_deviation"]])
        ax[0][0].bar(self.mass, intensity_of_molecular_ion_in_frag_spec, color="red", width=1.5, label="Molecule")
        ax[0][0].text(self.mass, (intensity_of_molecular_ion_in_frag_spec), str(true_fragment_list[0][3]), ha='center', va='bottom', rotation=0)
        for fragment in true_fragment_list:
            ax[0][0].bar(fragment[4], fragment[6], color="blue", label="Fragment")
            ax[0][0].text(fragment[4], (fragment[6]), str(fragment[7]), ha='center', va='bottom', rotation=0)
        ax[0][0].set_title("Fragment spectrum\n" + str(best_frag_spec.filter) + "\nIndex: " + str(best_frag_spec.index) + " RT: " + str(best_frag_spec.rt), wrap=True)
        ax[0][0].set_xlabel("ion mass / u")
        ax[0][0].set_ylabel("intensity / a.u.")
        ax[0][0].legend()

        text_str = ""
        text_str = text_str + "{:<12}|{:<15}|{:<12}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format("Comb. score", "Frag mass", "Intensity", "Frag formula", "NL formula", "Molec formula", "score F", "dev. NL", "score M")
        text_str = text_str + "{:_<12}|{:_<15}|{:_<12}|{:_<15}_{:_<15}___{:_<15}|{:_>10}|{:_>10}|{:_>10}\n".format("", "", "", "", "", "", "", "", "")
        for fragment in self.matching_fragments_list:
            comb_mi_score = sum([frag[1] for frag in self.matching_fragments_list if frag[3] == fragment[3]])
            text_str = text_str + "{:<12}|{:<15}|{:<12}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format(str(round(comb_mi_score, 1)), 
                                                                                                              str(round(fragment[4], 5)),
                                                                                                            str(round(fragment[6], 1)),
                                                                                                            str(fragment[7]),
                                                                                                            str(fragment[8]),
                                                                                                            str(fragment[3]),
                                                                                                            str(round(fragment[5], 1)),
                                                                                                            str(round(fragment[9], 3)),
                                                                                                            str(round(fragment[1], 1)))
        text_str = text_str + "================================================================================\n" + "All fragment predictions for spec: \n"

        text_str_all = "{:<10}|{:<15}|{:<12}\n".format("Mass", "Formula", "Score")
        text_str_all = text_str_all + "{:_<10}|{:_<15}|{:_<12}\n".format("", "", "")
        for pred_mass, pred in list(fragment_predictions.items()):
            try:
                try:
                    pred_formula = "".join([str(a) + str(n) for a, n in pred.best_formula_prediction.items()])
                except AttributeError:
                    pred_formula = "None"
                pred_score = pred.score_of_best_formula
                if pred_score is None:
                    pred_score = 0
                text_str_all = text_str_all + "{:<10}|{:<15}|{:<12}\n".format(str(round(pred_mass, 4)), str(pred_formula), str(round(pred_score, 1)))
            except Exception as e:
                print("Error creating text_str_all for summary plot.")
                print(e)
                print(traceback.format_exc())
        
        text_str = text_str + text_str_all

        props = dict(boxstyle='round', facecolor='grey', alpha=0.05)  # bbox features
        text = ax[1][0].text(0.02, 0.98, text_str, fontfamily="monospace", transform=ax[1][0].transAxes, fontsize=8,
                    verticalalignment="top", bbox=props)

        xic = self.xic
        peak_index = self.peak_index
        xic_peak_index = xic[2].index(min(xic[2], key=lambda x: abs(peak_index - x)))
        ax[0][1].scatter(self.ms_file.rt_list[peak_index], xic[1][xic_peak_index], color="red", marker="x", s=100)
        ax[0][1].plot(xic[0], xic[1], label="XIC measured")
        ax[0][1].set_title("xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["oa_mass_deviation_xic"]) + "_" + str(self.kwargs["oa_xic_requested_filter_mode"]), wrap=True)
        ax[0][1].set_xlabel("retention time / seconds")
        ax[0][1].set_ylabel("intensity / a.u.")

        for i, ax in enumerate(frag_axs):
            curr_xic = MS_functions.get_xic(self.ms_file.rawdata, frag_masses[i], self.kwargs["oa_mass_deviation_xic"], requested_filter_mode=self.best_frag_spec.filter_mode)
            ax.bar(self.ms_file.rt_list[peak_index], max(curr_xic[1]), color="red", label="observed time", alpha=0.5, width=10)
            ax.plot(curr_xic[0], curr_xic[1], color="blue", label="XIC " + str(round(frag_masses[i], 4)))
            ax.set_xlabel("retention time / s")
            ax.set_ylabel("intensity / a.u.")
            ax.set_title("Extracted Ion Chromatogram for fragment: " + str(round(frag_masses[i], 4)) + " u at RT: " + str(round(self.rt, 2)) + " s with mode: " + str(self.best_frag_spec.filter_mode), wrap=True)
            ax.legend()




        matplotlib.rcParams.update({'figure.autolayout': True})
        fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
        fig.clf()
        fig.clear()

    def append_oa_summary_to_raw_file_summary(self, summary_list):
        with open(self.ms_file.parentfolder + "/" + "SUMMARY.txt", "a") as oa_summary:
            for entry in summary_list:
                oa_summary.write(str(entry) + "\t")
            oa_summary.write("\n")

    def make_true_fragment_list(self):
        print()
        print()
        print("MAKING TRUE FRAGMENT LIST")
        matching_fragments_list = []
        true_fragment_list = []

        molecular_ion_mass = self.mass
        molecular_ion_best_approx = self.best_molecular_ion_prediction
        molecular_ion_best_approx_dict = MS_functions.get_formula_to_dict(molecular_ion_best_approx)
        molecular_ion_score = self.score_of_best_molecular_ion_prediction
        molecular_ion_intensity = self.intensity_of_molecular_ion
        molecular_ion_formula_score_dict = self.summarized_molecular_ion_formula_score_dict # {"{"C": 1, "H": 2, ....}": score, "{}": score2}
        mi_list = [molecular_ion_mass, molecular_ion_best_approx_dict, molecular_ion_score, molecular_ion_intensity, molecular_ion_formula_score_dict]
        true_fragment_list.append(mi_list)
        print("Molecular ion list: " + str(mi_list))
        print("Starting for loop...")
        for molecular_ion_pred_dict_str, mi_score in molecular_ion_formula_score_dict.items():
            molecular_ion_pred_dict = ast.literal_eval(molecular_ion_pred_dict_str)

            for f_mass, f_score_dict in self.fragment_predictions_formula_score_dicts.items():
                print()
                f_intensity = self.fragment_predictions[f_mass].intensity_of_ion
                neutral_loss_mass = (molecular_ion_mass - f_mass) - (self.kwargs["charge_of_measured_mass"]*0.0005)
                print(neutral_loss_mass)
                neutral_loss_formula_predictions = (MS_functions.get_formula_from_cache(self.kwargs["pred_formula_cache_folder_path"], neutral_loss_mass, self.kwargs["mass_deviation"]))[1]
                # neutral_loss_formula_prediction = [mass, alldict] ---> alldict = {formula: score, "C1H3O2": score, ...}
                print(neutral_loss_formula_predictions)
                if len(neutral_loss_formula_predictions) == 0:
                    continue
                found_pair = False
                for f_formula, f_score in f_score_dict.items():
                    f_formula_dict = ast.literal_eval(f_formula)
                    for nl_formula, nl_deviation in neutral_loss_formula_predictions.items():
                        nl_formula_dict = MS_functions.get_formula_to_dict(nl_formula)
                        print(nl_formula_dict)
                        print(f_formula_dict)
                        print(molecular_ion_pred_dict)
                        if (self.combine_and_sum_dicts(nl_formula_dict, f_formula_dict) == molecular_ion_pred_dict):
                            print("TRUETRUETRUEskjaskjhlgfaivwzbevwuief")
                            true_fragment_list.append([f_mass, f_formula_dict, f_score, f_intensity, neutral_loss_mass, nl_formula_dict, nl_deviation])
                            mi_formula_str = "".join([str(a) + str(n) for a, n in molecular_ion_pred_dict.items()])
                            f_formula_str = "".join([str(a) + str(n) for a, n in f_formula_dict.items()])
                            nl_formula_str = "".join([str(a) + str(n) for a, n in nl_formula_dict.items()])
                            matching_fragments_list.append([molecular_ion_mass, mi_score, molecular_ion_intensity, mi_formula_str, f_mass, f_score, f_intensity, f_formula_str, nl_formula_str, nl_deviation])
                            found_pair = True
                            #break
                    #if found_pair:
                        #break
        return true_fragment_list, matching_fragments_list
        
    def adjust_mass_to_closest_measured_mass(self, spec, masse=0):
        if masse == 0:
            masse = self.mass
        mass = min(list(spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(masse - x))
        return mass

    def isolate_and_prepare_peak(self, xic1, xic2, peak_rt, peakwidth=10):
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
                average_surrounding1 = sorted(peakintensity1)[1]
                if isinstance(average_surrounding1, int):
                    average_surrounding1 = float(average_surrounding1)
                if (not isinstance(average_surrounding1, float)) and (not isinstance(average_surrounding1, int)):
                    print("average_surrounding1 is not a float or int!!! --> " + str(average_surrounding1))
                    average_surrounding1 = 0.0001
            except Exception as e5:
                print("Error during isolating and preparing peak (1, 2): " + str(e5))
                average_surrounding1 = 0.0001
            try:
                peakintensity1 = [((i-average_surrounding1) / (intensity1_at_peak_rt-average_surrounding1)) for i in peakintensity1]
            except Exception as e4:
                print("Error during isolating and preparing peak (1, 1): " + str(e4))
                max_peakint1 = max(peakintensity1)
                if max_peakint1 == 0:
                    max_peakint1 = 0.001
                peakintensity1 = [((i-average_surrounding1) / (max_peakint1)) for i in peakintensity1]


            neighbour_list2 = xic2[1][int(index - 3*peakwidth):int(index + 3*peakwidth)]
            neighbour_list2 = [neighbour_list2[i] for i in range(len(neighbour_list2)) if (i < len(neighbour_list2)/3) or (i > ((len(neighbour_list2)/3) + (len(neighbour_list2)/2))) ]
            average_surrounding2 = sum(neighbour_list2) / len(neighbour_list2)
            try:
                average_surrounding2 = sorted(peakintensity2)[1]
                if (not isinstance(average_surrounding2, float)) and (not isinstance(average_surrounding2, int)):
                    print("average_surrounding2 is not a float or int!!! --> " + str(average_surrounding2))
                    average_surrounding2 = 0.0001
            except Exception as e3:
                print("Error during isolating and preparing peak (2, 2): " + str(e3))
                average_surrounding2 = 0.0001
            try:
                #print("Peakintensity2 before normalization: " + str(peakintensity2))
                peakintensity2 = [(float(i-average_surrounding2) / float(intensity2_at_peak_rt-average_surrounding2)) for i in peakintensity2]
            except Exception as e2:
                print("Error during isolating and preparing peak (2, 1): " + str(e2))
                max_peakint2 = max(peakintensity2)
                if max_peakint2 == 0:
                    max_peakint2 = 0.0001
                peakintensity2 = [((i-average_surrounding2) / (max_peakint2)) for i in peakintensity2]
                print("Peakintensity1: " + str(peakintensity1))
                print("Peakintensity2: " + str(peakintensity2))
                print("Peak RT 1: " + str(peak1_rt))
                print("Peak RT 2: " + str(peak2_rt))
        except Exception as e1:
            print("Error during isolating and preparing peak: " + str(e1))
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
            print("Peakintensity1: " + str(peakintensity1))
            print("Peakintensity2: " + str(peakintensity2))
            print("Peak RT 1: " + str(peak1_rt))
            print("Peak RT 2: " + str(peak2_rt))

        return peak1_rt, peakintensity1, peak2_rt, peakintensity2

    def get_peak_similarity(self, x1, x2, y1, y2):
        P = np.array([x1, y1]).T
        Q = np.array([x2, y2]).T
        area = similaritymeasures.area_between_two_curves(P, Q)
        return area

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

        self.make_oa_log_entry("INFO:\t" + "Starting first prediction of molecular ion...")

        self.molecular_ion_prediction = Prediction(self.ms_file, self.mass, best_molecular_ion_spec, spec_before=self.mi_spec_before, spec_after=self.mi_spec_after,
                                                   absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/",
                                                   pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                                   pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_matplotlib_plot_of_isotopologues"],
                                                   pred_save_go_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_go_plot_of_isotopologues"],
                                                   pred_save_xic_plot=self.kwargs["oa_molecular_ion_pred_save_xic_plot"],
                                                   pred_save_detailed_log=self.kwargs["oa_molecular_ion_pred_save_detailed_log"],
                                                   pred_return_if_no_peak_is_found=self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"],
                                                   **additional_kwargs)
        self.make_oa_log_entry("INFO:\t" + "Finished first prediction of molecular ion...")
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
        self.make_oa_log_entry("INFO:\t" + "Trying to other spectra...")
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
                              spec_requested_filter_mode=best_molecular_ion_spec.filter_mode, 
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
        self.make_oa_log_entry("INFO:\t" + "Finished searching for other spectra...")
        self.make_oa_log_entry("INFO:\t" + "Starting prediction of molecular ion with multiple spectra...")
        #add the scores of all the available prediction formula_score_dicts and create a summarized formula_score_dict
        self.summarized_molecular_ion_formula_score_dict, all_formula_score_dicts, best_pred = self.get_formula_score_dict_with_multiple_specs(available_specs, self.mass, return_if_no_peak_found=self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"])
        self.make_oa_log_entry("INFO:\t" + "Finished prediction of molecular ion with multiple spectra...")

        if len(self.summarized_molecular_ion_formula_score_dict) == 0:
            print("No molecular ion formula could be predicted!!!")
        elif not self.molecular_ion_prediction.best_formula_prediction == ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0]):
            additional_kwargs = copy.deepcopy(self.kwargs)
            pop_keys = ["absolute_pred_folder", "pred_formula_cache_folder_path", "pred_save_matplotlib_plot_of_isotopologues", "pred_save_go_plot_of_isotopologues", "pred_save_xic_plot", "pred_save_detailed_log", "pred_return_if_no_peak_is_found"]
            for key in pop_keys:
                try:
                    additional_kwargs.pop(key)
                except KeyError:
                    continue
            for curr_spec in available_specs:
                Prediction(self.ms_file, self.mass, curr_spec,
                                    absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/",
                                    pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                    pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_matplotlib_plot_of_isotopologues"],
                                    pred_save_go_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_go_plot_of_isotopologues"],
                                    pred_save_xic_plot=self.kwargs["oa_molecular_ion_pred_save_xic_plot"],
                                    pred_save_detailed_log=self.kwargs["oa_molecular_ion_pred_save_detailed_log"],
                                    pred_return_if_no_peak_is_found=self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"],
                                    **additional_kwargs)
            self.make_oa_log_entry("INFO:\t" + "Old prediction did not match with the prediction of multiple spectra!")
            self.make_oa_log_entry("INFO:\t" + "Plots are created for every spectrum next to the original index. Plots will be available.")
            with open(self.kwargs["one_analysis_folder"] + "predictions/BEST_FORMULA_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0]).items()])) + ".txt", "a") as txt_file:
                txt_file.write("Best prediction according to three spectra which are located next to the original peak: " + str(ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0])))
                txt_file.write("\n")
                txt_file.write("Summarized formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))
                txt_file.write("\n")

        if self.molecular_ion_prediction.best_formula_prediction is None:
            self.make_oa_log_entry("ERROR:\t" + "No molecular formula found in the first spectrum! Trying to find prediction for the otehr spectra...")
            print("No formula found! Searching for prediction success next to original spectrum...")
            if best_pred is not None:
                self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Best prediction was adjusted to another spectrum. This prediction is not working!")
                self.molecular_ion_prediction = best_pred
                self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Best prediction was adjusted to this spectrum. Continuing with prediction...")
                self.make_oa_log_entry("INFO:\t" + "Best prediction was adjusted to another spectrum. Continuing with prediction...")
                print("Best prediction was adjusted to another spectrum. Continuing with prediction...")

        self.summarized_molecular_ion_formula_score_dict = {f: s for f, s in self.summarized_molecular_ion_formula_score_dict.items() if s > self.kwargs["oa_molecular_ion_reject_formula_if_score_lower_than"]}
        self.summarized_molecular_ion_formula_score_dict = dict(sorted(self.summarized_molecular_ion_formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        self.make_oa_log_entry("INFO:\t" + "Summarized formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))
        print("Summarized formula score dict: " + str(self.summarized_molecular_ion_formula_score_dict))



        self.molecular_ion_prediction.make_op_log_entry("=============================================================")
        try:
            self.best_molecular_ion_prediction = str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0]).items()]))  
            self.score_of_best_molecular_ion_prediction = list(self.summarized_molecular_ion_formula_score_dict.values())[0]
            self.make_oa_log_entry("INFO:\t" + "New best_molecular_ion_prediction: " + str(self.best_molecular_ion_prediction))
            self.make_oa_log_entry("INFO:\t" + "New score_of_best_molecular_ion_prediction: " + str(self.score_of_best_molecular_ion_prediction))
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
        self.make_oa_log_entry("INFO:\t" + "Molecular Ion Prediction finished. \nBest formula prediction: " + str(self.best_molecular_ion_prediction) + " single score: " + str(self.score_of_best_molecular_ion_prediction))
        return self.molecular_ion_prediction
    
    def get_formula_score_dict_with_multiple_specs(self, specs, mass, return_if_no_peak_found=True, absolute_pred_subfolder_path="predictions/"):
        self.make_oa_log_entry("INFO:\t" + "Starting prediction with multiple specs...")
        predictions = []
        for spec in specs:
            self.make_oa_log_entry("INFO:\t" + "Starting prediction with spec with index: " + str(spec.index))
            try:
                additional_kwargs = copy.deepcopy(self.kwargs)
                pop_keys = ["absolute_pred_folder", "pred_formula_cache_folder_path", "pred_save_matplotlib_plot_of_isotopologues", "pred_save_go_plot_of_isotopologues", "pred_save_xic_plot", "pred_save_detailed_log", "pred_return_if_no_peak_is_found"]
                for key in pop_keys:
                    try:
                        additional_kwargs.pop(key)
                    except KeyError:
                        continue
                self.make_oa_log_entry("INFO:\t" + "Starting prediction just now.")
                pred = Prediction(self.ms_file, mass, spec,
                                  absolute_pred_folder=self.kwargs["one_analysis_folder"] + absolute_pred_subfolder_path,
                                  pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                  pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_multiplespec_pred_save_matplotlib_plot_of_isotopologues"],
                                  pred_save_go_plot_of_isotopologues=self.kwargs["oa_multiplespec_pred_save_go_plot_of_isotopologues"],
                                  pred_save_xic_plot=self.kwargs["oa_multiplespec_pred_save_xic_plot"],
                                  pred_save_detailed_log=self.kwargs["oa_multiplespec_pred_save_detailed_log"],
                                  pred_return_if_no_peak_is_found=return_if_no_peak_found,
                                  **additional_kwargs)
                predictions.append(pred)
                self.make_oa_log_entry("INFO:\t" + "Finished prediction with spec with index: " + str(spec.index))
            except Exception as e:
                predictions.append(None)
                print("Error in getting prediction with multiple specs: " + str(e))
                print(traceback.format_exc())
                continue
        print(predictions)
        self.make_oa_log_entry("INFO:\t" + "Finished prediction with multiple specs...")
        self.make_oa_log_entry("INFO:\t" + "Starting summarizing formula score dicts...")
        summarized_formula_score_dict = {}
        for pred in predictions:
            if pred is not None:
                summarized_formula_score_dict = self.combine_and_sum_dicts(summarized_formula_score_dict, pred.formula_score_dict)
                print("Formula score dict after adding prediction: " + str(summarized_formula_score_dict))
        summarized_formula_score_dict = dict(sorted(summarized_formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        all_formula_score_dicts = [pred.formula_score_dict for pred in predictions if pred is not None]
        self.make_oa_log_entry("INFO:\t" + "Finished summarizing formula score dicts...")
        self.make_oa_log_entry("INFO:\t" + "Searching for best prediction...")
        best_prediction = predictions[0]
        best_score = -9999999
        for entry in predictions:
            if entry is not None:
                if entry.best_formula_prediction is None:
                    continue
                if entry.score_of_best_formula > best_score:
                    best_score = entry.score_of_best_formula
                    best_prediction = entry
        self.make_oa_log_entry("INFO:\t" + "Best prediction found: " + str(best_prediction))
        return summarized_formula_score_dict, all_formula_score_dicts, best_prediction
        
    def get_fragment_predictions(self, make_good_fragment_formula_prediction=False):
        if self.best_frag_spec is None:
            self.make_oa_log_entry("INFO:\t" + "No best fragment prediction spec provided. Returning...")
            print("No best fragment prediction spec provided. Returning...")
            return None
        if make_good_fragment_formula_prediction:
            self.make_oa_log_entry("INFO:\t" + "Starting prediction of fragment ions with good formula prediction...")
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
            self.make_oa_log_entry("INFO:\t" + "Finished searching for other spectra...")
        print("Starting prediction of fragment ions...")

        self.possible_fragment_masses = [m for m in self.best_frag_spec.summarized_masses if m < self.mass-0.35 and self.best_frag_spec.summarized_intensities[self.best_frag_spec.summarized_masses.index(m)] > (self.intensity_of_molecular_ion * self.kwargs["oa_include_frag_intensity_noise_multiplier"])]
        print(self.possible_fragment_masses)
        self.possible_fragment_masses = sorted(self.possible_fragment_masses, key=lambda m: self.best_frag_spec.summarized_intensities[self.best_frag_spec.summarized_masses.index(m)], reverse=True)
        print(self.possible_fragment_masses)
        self.possible_fragment_masses = [self.possible_fragment_masses[m] for m in range(len(self.possible_fragment_masses)) if m <= self.kwargs["oa_fragments_absolute_max_number_of_fragment_masses"]]
        print("Possible fragment masses: " + str(self.possible_fragment_masses))
        self.make_oa_log_entry("INFO:\t" + "Possible fragment masses: " + str(self.possible_fragment_masses))
        self.fragment_predictions = {}
        self.fragment_predictions_formula_score_dicts = {}

        while len(self.possible_fragment_masses) >= 1:
            frag_mass = self.possible_fragment_masses[0]
            self.possible_fragment_masses.pop(self.possible_fragment_masses.index(frag_mass))

            try:
                if self.kwargs["oa_fragments_do_good_peak_comparison_with_area_between_curves"]:
                    os.makedirs(self.kwargs["one_analysis_folder"] + "predictions/fragments/peak_matching/", exist_ok=True)
                    fragment_xic = MS_functions.get_xic(self.ms_file.rawdata, frag_mass, self.kwargs["oa_mass_deviation_xic"], requested_filter_mode=self.best_frag_spec.filter_mode)
                    peak1_rt, peakintensity1, peak2_rt, peakintensity2 = self.isolate_and_prepare_peak(xic1=self.xic, xic2=fragment_xic, peak_rt=self.rt)
                    print("Fragment mass:   " + str(frag_mass))
                    try:
                        area_between_curves = self.get_peak_similarity(peak1_rt, peak2_rt, peakintensity1, peakintensity2)
                    except Exception as e2:
                        print("Error in calculating area between the two curves. " + str(e2))
                        print("Peakintensity1: " + str(peakintensity1))
                        print("Peakintensity2: " + str(peakintensity2))
                        print("Peak RT 1: " + str(peak1_rt))
                        print("Peak RT 2: " + str(peak2_rt))
                        continue
                    if area_between_curves > self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"]:
                        print("No good fragment peak shape was detected. Continuing with the next fragment...")
                        self.make_oa_log_entry("INFO:\t" + "Fragment peak with mass: " + str(frag_mass) + "  -> Does not have a good fragment peak shape. Area between curves too high: " + str(area_between_curves) + " Continuing....")
                        fig = matplotlib.figure.Figure()
                        ax = fig.subplots()
                        ax.plot(peak1_rt, peakintensity1, color="blue", label="Normalized molecular ion peak")
                        ax.plot(peak2_rt, peakintensity2, color="red", label="Normalized fragment peak")
                        ax.text(min(peak1_rt), max(peakintensity1), "Area between curves: " + str(round(area_between_curves, 2)))
                        ax.set_title("Peak matching evaluation" + str(round(self.mass, 4)) + " / " + str(round(frag_mass, 4)))
                        ax.set_xlabel("retention time / seconds")
                        ax.set_ylabel("normalized intensity / a.u.")
                        ax.legend()
                        matplotlib.rcParams.update({'figure.autolayout': True})
                        save_filepath = self.kwargs["one_analysis_folder"]  + "predictions/fragments/peak_matching/" + "FALSE_Fragmass_" + str(round(frag_mass, 4)) + ".png"
                        fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
                        fig.clf()
                        fig.clear()
                        continue
                    else:
                        fig = matplotlib.figure.Figure()
                        ax = fig.subplots()
                        ax.plot(peak1_rt, peakintensity1, color="blue", label="Normalized molecular ion peak")
                        ax.plot(peak2_rt, peakintensity2, color="red", label="Normalized fragment peak")
                        ax.text(min(peak1_rt), max(peakintensity1), "Area between curves: " + str(round(area_between_curves, 2)))
                        ax.set_title("Peak matching evaluation" + str(round(self.mass, 4)) + " / " + str(round(frag_mass, 4)))
                        ax.set_xlabel("retention time / seconds")
                        ax.set_ylabel("normalized intensity / a.u.")
                        ax.legend()
                        matplotlib.rcParams.update({'figure.autolayout': True})
                        save_filepath = self.kwargs["one_analysis_folder"] + "predictions/fragments/peak_matching/" + "TRUE_Fragmass_" + str(round(frag_mass, 4)) + ".png"
                        fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
                        fig.clf()
                        fig.clear()
            except Exception as e:
                print("Error in making peak matching plot in get_fragment_predictions!" + str(e))
                self.make_oa_log_entry("ERROR:\t" + "Error making peak matching plot!")
                self.make_oa_log_entry("ERROR:\t" + str(traceback.format_exc()))

            try:
                additional_kwargs = copy.deepcopy(self.kwargs)
                pop_keys = ["absolute_pred_folder", "pred_formula_cache_folder_path", "pred_save_matplotlib_plot_of_isotopologues", "pred_save_go_plot_of_isotopologues", "pred_save_xic_plot", "pred_save_detailed_log", "pred_return_if_no_peak_is_found"]
                for key in pop_keys:
                    try:
                        additional_kwargs.pop(key)
                    except KeyError:
                        continue

                curr_prediction = Prediction(self.ms_file, frag_mass, self.best_frag_spec, spec_before=self.frag_spec_before, spec_after=self.frag_spec_after,
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
                    print("Remaining length of fragment prediction list: " + str(len(self.possible_fragment_masses)))
                    continue
                self.fragment_predictions[frag_mass] = curr_prediction

                if make_good_fragment_formula_prediction:
                    self.fragment_predictions_formula_score_dicts[frag_mass], _, curr_prediction = self.get_formula_score_dict_with_multiple_specs(available_specs, frag_mass, return_if_no_peak_found=False)
                    print("Good fragment prediction finished.")
                    print(self.fragment_predictions_formula_score_dicts[frag_mass])
                else:
                    self.fragment_predictions_formula_score_dicts[frag_mass] = self.fragment_predictions[frag_mass].formula_score_dict

                self.make_oa_log_entry("INFO:\t" + "Finished prediction of fragment ion with mass: " + str(frag_mass) + " Best formula prediction: " + str(self.fragment_predictions[frag_mass].best_formula_prediction) + " Score: " + str(self.fragment_predictions[frag_mass].score_of_best_formula))
                simulated_isotopo_mass_abundance_dict = curr_prediction.simulated_isotopologue_pattern_for_best_formula
                print(self.possible_fragment_masses)
                if simulated_isotopo_mass_abundance_dict is None:
                    continue
                for entry in list(simulated_isotopo_mass_abundance_dict.keys()):
                    print("Reducing the number of possible fragments by looking at simulated isotope pattern... \t old length --> \t new length")
                    print(entry)
                    print(len(self.possible_fragment_masses))
                    self.possible_fragment_masses = [m for m in self.possible_fragment_masses if (abs( ((entry - m)/m)*1000000 ) >= self.kwargs["mass_deviation"])
                                                                                                    or ((simulated_isotopo_mass_abundance_dict[entry]*curr_prediction.intensity_of_ion) + self.kwargs["pred_minimum_assumed_noise"] <= self.best_frag_spec.summarized_mass_intensity_dict[m]) ]
                    print(len(self.possible_fragment_masses))
                    print()
            except Exception as e:
                print("Error in making fragment peak prediction in get_fragment_predictions!" + str(e))
                self.make_oa_log_entry("ERROR:\t" + "Error making fragment prediction!")
                self.make_oa_log_entry("ERROR:\t" + str(traceback.format_exc()))

        return self.fragment_predictions, self.fragment_predictions_formula_score_dicts

    def plot_and_save_xic(self, save_filepath, peak_index, xic):
        fig = matplotlib.figure.Figure()
        ax = fig.subplots()
        xic_peak_index = xic[2].index(min(xic[2], key=lambda x: abs(peak_index - x)))
        ax.scatter(self.ms_file.rt_list[peak_index], xic[1][xic_peak_index], color="red", marker="x", s=100)
        ax.plot(xic[0], xic[1], label="XIC measured")
        ax.set_title("xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["oa_mass_deviation_xic"]) + "_" + str(self.kwargs["oa_xic_requested_filter_mode"]))
        ax.set_xlabel("retention time / seconds")
        ax.set_ylabel("intensity / a.u.")
        ax.legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
        fig.clf()
        fig.clear()
        metadata = PIL.PngImagePlugin.PngInfo()
        metadata.add_text("xic_retention_time", str(xic[0]))
        metadata.add_text("xic_intensity_list", str(xic[1]))
        metadata.add_text("xic_original_index_list", str(xic[2]))
        target_image = PIL.Image.open(save_filepath)
        target_image.save(save_filepath, pnginfo=metadata)

    def plot_summary_xic_of_mi_and_fragments(self, save_filepath):
        #e.g. [173.00860345677103, 29845.72513229523, 1214042584.515625, 'C6H5O6', 111.00881936704444, 2496.3943129810064, 21546167.083251953, 'C5H3O3', 'C1H2O3', 1.8695119725597358]
        #self.true_fragment_list = [[mass_mi, score_mi, int_mi, formula_mi, mass_frag1, score_frag1, int_frag1, formula_frag1, nl_frag1, ppmdev_nl_frag1], [mass_mi, ..., mass_frag2, ...]]
        mi_mass = self.true_fragment_list[0][0]
        frag_masses = [frag[4] for frag in self.true_fragment_list]
        frag_pred = [frag[7] for frag in self.true_fragment_list]
        frag_intensities = [frag[6] for frag in self.true_fragment_list]
        mi_xic = MS_functions.get_xic(self.ms_file.rawdata, mi_mass, ((self.kwargs["mass_deviation"]*mi_mass)/1000000), requested_filter_mode="Full scan")
        fig = matplotlib.figure.Figure()
        ax = fig.subplots(nrows=1, ncols=2)
        ax[0].plot(mi_xic[0], mi_xic[1], label="Molecular Ion")
        ax[1].plot(mi_xic[0], [inten/self.true_fragment_list[0][2] for inten in mi_xic[1]], label="Molecular Ion")
        for i in range(len(frag_masses)):
            try:
                xic = MS_functions.get_xic(self.ms_file.rawdata, frag_masses[i], ((self.kwargs["mass_deviation"]*frag_masses[i])/1000000), requested_filter_mode="AIF")
                ax[1].plot(xic[0], [inten / frag_intensities[i] for inten in xic[1]], label="Fragment " + str(frag_pred[i]), alpha=0.4)
                ax[0].plot(xic[0], xic[1], label="Fragment " + str(frag_pred[i]), alpha=0.4)
            except Exception as e:
                print("Error in creating summary XIC of fragments. " + str(e))
        fig.suptitle("Summary for mass " + str(round(self.mass, 4)))
        ax[1].set_ylim([0, 1.2])
        ax[0].set_xlabel("retention time / seconds")
        ax[0].set_ylabel("intensity / a.u.")
        ax[0].legend()
        ax[1].legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
        fig.clf()
        fig.clear()

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

                if iteration_step_lower == 0 and integration_step_upper == 0:
                    summed_intensity = old_abundances_list[index_of_highest_abundance]
                    weighted_mass_average = old_masses_list[index_of_highest_abundance]
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance ]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance ]
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

    mzml_filename = "C://Users//Admin//Desktop//UVenture//UVenture//webserver_save//mzml_files//Cal2.mzML"

    ms_file = MS_File(mzml_filename)

    xic1 = MS_functions.get_xic(ms_file.rawdata, 138.0191, 0.001, requested_filter_mode="Full scan")
    xic2 = MS_functions.get_xic(ms_file.rawdata, 108.0207, 0.001, requested_filter_mode="AIF")
    peak1_rt, peakintensity1, peak2_rt, peakintensity2 = OneAnalysis.isolate_and_prepare_peak(xic1=xic1, xic2=xic2, peak_rt=259)

    settings_dict.pop("spec_requested_filter_mode")

    myspec = Spec(ms_file, 345, spec_requested_filter_mode="AIF", **settings_dict)
    input()

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


    







