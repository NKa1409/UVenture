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
        default_kwargs = {"parentfolder":str(".".join(filename.split(".")[:-1]) + "/"),
                          "peak_infos":[],
                          "logfile_filepath": str(".".join(filename.split(".")[:-1]) + "/" + "logfile.txt") }
        kwargs = {**default_kwargs, **kwargs}
        self.filename = filename
        self.parentfolder = kwargs["parentfolder"]
        os.makedirs(self.parentfolder, exist_ok=True)
        self.peak_infos = kwargs["peak_infos"]

        self.file = None
        self.rawdata = None
        self.all_ms_spectra = []
        self.all_filters = []
        self.method_duration = None
        self.available_modes = []
        self.tic = []
        self.open_mzml_file()
        self.rt_list = [element["scanList"]["scan"][0]["scan time"] for element in self.rawdata]
        self.spectra_object_dict = {}
        self.analytes_list = []
    
    def __call__(self):
        return self.rawdata
    
    def __repr__(self):
        return str(self.filename)
    
    def open_mzml_file(self):
        self.file = mzml.read(self.filename)
        self.rawdata = list(self.file)
        for index in range(len(self.rawdata)):
            masses = list(self.rawdata[index]["m/z array"])
            intensities = list(self.rawdata[index]["intensity array"])
            self.all_ms_spectra.append({"masses":masses, "intensities":intensities})
            self.tic.append(self.rawdata[index]["total ion current"])
            self.all_filters.append(self.rawdata[index]["scanList"]["scan"][0]["filter string"])
        for filter_string in self.all_filters:
            if " d " in filter_string and "@hcd" in filter_string:
                self.available_modes.append("MS/MS")
            elif " d " not in filter_string and "hcd" not in filter_string:
                self.available_modes.append("Full scan")
            elif " d " not in filter_string and "hcd" in filter_string:
                self.available_modes.append("AIF")
        self.available_modes = list(set(self.available_modes))
        self.method_duration = self.rawdata[-1]["scanList"]["scan"][0]["scan time"]
        return self

    def do_one_analysis(self, mass, rt, **kwargs):
        oa_object = OneAnalysis(self, mass, rt, **kwargs)
        self.analytes_list.append(oa_object)
        return True
        

class Spec:
    def __init__(self, ms_file, index, **kwargs):
        self.ms_file = ms_file
        self.index = index
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

        if "spec_subfolder" in kwargs:
            spec_folder = self.ms_file.parentfolder + "/spectra/" + kwargs["spec_subfolder"]
        else:
            spec_folder = str(self.ms_file.parentfolder + "/spectra/" + "massspec_" + str(self.index) + "/")
        if "absolute_spec_folder" in kwargs:
            spec_folder = kwargs["absolute_spec_folder"]
        if "spec_subfolder" in kwargs and "absolute_spec_folder" in kwargs:
            print("ERROR: You can only use either 'spec_subfolder=' or 'absolute_spec_folder=', not both!")
            print("Using absolute_spec_folder now...")
            print(kwargs["absolute_spec_folder"])
        default_kwargs = {"spec_folder":spec_folder,
                          "requested_filter_mode":"Full scan",
                          "save_plot":True,
                          "requested_ms_ms_mass":"" ,
                          "spec_log_filepath":spec_folder + "spectrum_creation_log_" + str(self.index) + ".txt",
                          "mass_deviation": 11,
                          "save_additional_spectrum_info":False,
                          "additional_spectrum_info_filepath":spec_folder + "additional_spectrum_info_" + str(self.index) + ".txt"}
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["spec_folder"], exist_ok=True)

        self.make_spec_log_entry("=============================================================")
        self.make_spec_log_entry("INFO:\t" + "Searching spectrum for retention time: " + str(self.rt))
        self.make_spec_log_entry("INFO:\t" + "Requested filter mode: " + str(self.kwargs["requested_filter_mode"]))
        if self.kwargs["requested_filter_mode"] == "MS/MS":
            self.make_spec_log_entry("INFO:\t" + "Requested MS/MS mass: " + str(self.kwargs["requested_ms_ms_mass"]))
        self.make_spec_log_entry("INFO:\t" + "Saving plot?: " + str(self.kwargs["save_plot"]))

        self.mass_intensity_dict = self.get_mass_spec(self.index, requested_mode=self.kwargs["requested_filter_mode"], requested_ms_ms_mass=self.kwargs["requested_ms_ms_mass"])
        if self.index in list(self.ms_file.spectra_object_dict.keys()):
            print("Summarized spectra already in dict")
            self.summarized_mass_intensity_dict = self.ms_file.spectra_object_dict[self.index].summarized_mass_intensity_dict
        else:
            print("Summarizing mass intensity dict...")
            self.summarized_mass_intensity_dict = self.summarize_mass_intensity_dict_with_deviation(dict(zip(self.masses, self.intensities)), deviation=self.kwargs["mass_deviation"])
            print("Finished summarizing mass intensity dict")
            

        self.summarized_masses = list(self.summarized_mass_intensity_dict.keys())
        self.summarized_intensities = list(self.summarized_mass_intensity_dict.values())
        print(len(self.masses))
        print(len(self.intensities))
        print(len(self.summarized_masses))
        print(len(self.summarized_intensities))



        if self.kwargs["save_plot"] == True:
            self.create_go_plot()
            self.create_matplotlib_plot()

        self.kwargs["additional_spectrum_info_filepath"] = self.kwargs["spec_folder"] + "additional_spectrum_info_" + str(self.index) + ".txt"
        if self.kwargs["save_additional_spectrum_info"] == True:
            log_f = open(self.kwargs["additional_spectrum_info_filepath"], "a")
            log_f.write(str(self.filter) + "\n")
            log_f.write(str(self.index) + "\n")
            log_f.write(str(self.filter_mode) + "\n")
            log_f.write(str(self.summarized_mass_intensity_dict) + "\n")
            log_f.write(str(self.spec_rawdata) + "\n")
            log_f.close()
        
        self.ms_file.spectra_object_dict[self.index] = self
    
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
            
    def summarize_mass_intensity_dict_with_deviation_old_and_slow(self, dictio, deviation=11):
        print("=============================================================")
        print("=============================================================")
        start_time = datetime.datetime.now()
        print("Start time: " + str(datetime.datetime.now()))
        print("summarizing dict according to old method")
        dictio = {k: v for k, v in dictio.items() if v >= 1}
        print("old length of start dictio:")
        print(len(dictio))
        print("minimal deviation between elements: " + str(self.min_deviation(list(dictio.keys()))))
        old_masses_list = list(dictio.keys())
        old_abundances_list = list(dictio.values())
        new_masses_list = []
        new_abundances_list = []
        self.make_spec_log_entry("=============================================================")
        self.make_spec_log_entry("INFO:\t" + "Summarizing mass intensity dict with deviation: " + str(deviation))
        for entry in range(len(old_masses_list)):
            try:
                mass_dev_list = [abs(((m - old_masses_list[entry]) / old_masses_list[entry]) * 1000000) for m in old_masses_list]
                curr_mass_plus_deviation_list = [old_masses_list[m] for m in range(len(old_masses_list)) if mass_dev_list[m] <= deviation]
                curr_abundance_plus_deviation_list = [old_abundances_list[m] for m in range(len(old_masses_list)) if mass_dev_list[m] <= deviation]

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

                cumm_abundance = sum(curr_abundance_plus_deviation_list)
                cumm_mass = sum([a*b for a,b in zip(curr_mass_plus_deviation_list, curr_abundance_plus_deviation_list)]) / cumm_abundance
                new_masses_list.append(cumm_mass)
                new_abundances_list.append(cumm_abundance)

            except Exception as e:
                print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
                print(e)
                break
        outdict = dict(zip(new_masses_list, new_abundances_list))
        print("new length of dictio according to old method")
        print(len(outdict))
        print("minimal deviation between elements: " + str(self.min_deviation(new_masses_list)))
        print("End time:" + str(datetime.datetime.now()))
        print("Time taken = " + str(datetime.datetime.now() - start_time))
        self.make_spec_log_entry("Finished calculating summarized mass intensity dict!!")
        self.make_spec_log_entry("=============================================================")
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
        image_filepath = self.kwargs["spec_folder"] + "Mass_spectrum_index" + str(self.index) + "_" + str(
            self.filter_mode.replace("/", "")) + "_" + str(self.ms_ms_masses) + ".png"
        plt.savefig(image_filepath, bbox_inches='tight', dpi=1000)
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
                self.make_spec_log_entry("WARNING:\t" + "NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...", error=True)
                print("NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...")
                return
            
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
                self.make_spec_log_entry("INFO:\t" + "FOUND REQUESTED SPECTRUM!")
                self.make_spec_log_entry("INFO:\t" + "Deviation to initial index: " + str(abs(curr_index - start_index)))
                self.make_spec_log_entry("INFO:\t" + "RT Difference to requested spectrum: " + str(round(self.ms_file.rt_list[curr_index] - self.ms_file.rt_list[start_index], 2)) + " seconds")
                print("Found requested spectrum")
                break

            elif (filter_mode == requested_mode) and (round(requested_ms_ms_mass, 1) in ms_ms_masses):
                self.make_spec_log_entry("INFO:\t" + "FOUND MS/MS SPECTRUM")
                self.make_spec_log_entry("INFO:\t" + "Deviation to initial index: " + str(abs(curr_index - start_index)))
                self.make_spec_log_entry("INFO:\t" + "RT Difference to requested spectrum: " + str(round(self.ms_file.rt_list[curr_index] - self.ms_file.rt_list[start_index], 2) + " seconds"))
                self.make_spec_log_entry("INFO:\t" + "All available MS/MS masses: " + str(ms_ms_masses))
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

    def get_mass_spec(self, index, requested_mode="whatever", requested_ms_ms_mass=""):
        
        index = self.search_for_required_spec_close_to_rt(index, requested_mode=requested_mode, requested_ms_ms_mass=requested_ms_ms_mass)
        self.index = index
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
        masse = min(list(spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(mass - x))
        self.mass = masse
        self.spec = spec
        self.rt = self.spec.rt
        self.formulas_score_dict = {}
        self.best_formula_prediction = None
        self.score_of_best_formula = None
        self.intensity_of_ion = None
        self.identified_peaks_for_mass = None

        if "prediction_subfolder" in kwargs:
            pred_folder = self.ms_file.parentfolder + "/predictions/" + kwargs["prediction_subfolder"]
        else:
            pred_folder = str(self.ms_file.parentfolder + "/predictions/" + str(self.spec.index) + "_" + str(round(self.mass, 2)) + "/")
        if "absolute_pred_folder" in kwargs:
            pred_folder = kwargs["absolute_pred_folder"]
        if "prediction_subfolder" in kwargs and "absolute_pred_folder" in kwargs:
            print("ERROR: You can only use either 'prediction_subfolder=' or 'absolute_pred_folder=', not both!")
            print("Using absolute_pred_folder now...")
            print(kwargs["absolute_pred_folder"])
        
        default_kwargs = {"pred_folder":pred_folder,  
                          "mass_deviation":11,
                          "op_log_filepath": pred_folder + "prediction_log_for_mass_" + str(round(self.mass, 4)) + ".txt",
                          "charge_of_measured_mass":-1,
                          "formula_cache_folder_path":"U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//",
                          "save_detailed_log":True,
                          "max_ppm_deviation_change_for_isotopologue":2,
                          "minimum_assumed_noise":10000,
                          "noise_divisor_for_isotopologue_calculation":5,
                          "add_value_to_score_if_isotopo_was_found":250,
                          "isotopologue_score_e_function_exponent":0.4,
                          "multiplier_isotopologue_influence_of_ppm_deviation_on_score":30,
                          "stop_isotopologue_search_if_score_lower_than":-80,
                          "ppm_deviation_score_multiplier":14,
                          "include_likelyhood_of_formula":True,
                          "reject_formula_if_score_lower_than":10,
                          "save_plot_of_isotopologues":False,
                          "filepath_of_plot_of_isotopologues":pred_folder + "isotopologues_plot_" + str(round(self.mass, 4)) + ".png" }
        
        self.kwargs = {**default_kwargs, **kwargs}
        os.makedirs(self.kwargs["pred_folder"], exist_ok=True)

        self.intensity_of_ion = sum([self.spec.intensities[i] for i in range(len(self.spec.masses)) if abs(((self.spec.masses[i] - mass)/mass)*1000000) <= self.kwargs["mass_deviation"]])

        self.make_op_log_entry("INFO:\t" + "Calculating formula of the ion with mass: " + str(self.mass) + " at retention time: " + str(self.spec.rt) + " at index: " + str(self.spec.index))
        self.make_op_log_entry("INFO:\t" + "Intensity of the molecular ion: " + str(self.intensity_of_ion))
        self.make_op_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))

        _, possible_formulas = MS_functions.get_formula_from_cache(self.kwargs["formula_cache_folder_path"], self.mass, self.kwargs["mass_deviation"])
        if self.kwargs["charge_of_measured_mass"] < 0:
            self.make_op_log_entry("INFO:\t" + "Charge of the measured mass is negative. Removing formulas with Li, Na and K...")
            possible_formulas = {key: value for key, value in possible_formulas.items() if "Na" not in key and "K" not in key and "Li" not in key}
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("OUTPUT OF SIMPLE FORMULA PREDICTION (ONLY MASS DEVIAITON):")
        for formula, deviation in possible_formulas.items():
            self.make_op_log_entry("{:>20} \t {:>20}".format(str(formula), str(round(deviation, 2))))
        
        molecule_formulas = [[MS_functions.get_formula_to_dict(formula), deviation] for formula, deviation in possible_formulas.items()]
        self.mass_possible_formulas = molecule_formulas

        formula_score_dict = {}
        for entry in molecule_formulas:
            self.make_op_log_entry("")
            self.make_op_log_entry("=============================================================")
            isotope_check = self.check_for_isotope_pattern(entry[0])
            formula_score = sum(e_s[3] for e_s in list(isotope_check.values())) - abs(self.kwargs["ppm_deviation_score_multiplier"] * entry[1])
            self.make_op_log_entry("INFO:\t" + "Formula: " + str("".join([str(a) + str(n) for a, n in entry[0].items()])) + " \t Score: " + str(round(formula_score, 2)))
            formula_score_dict[str(entry[0])] = formula_score
            if self.kwargs["save_detailed_log"] == True:
                log_f = open(self.kwargs["op_log_filepath"], "a")
                log_f.write("OUTPUT OF FORMULA PREDICTION WITH MASS SPECTRUM (ISOTOPOLOGUES)\n")
                log_f.write("Mass of ion: " + str(self.mass) + "\n")
                log_f.write("PPM DEVIATION OF MOLECULE FORMULA: " + str(round(entry[1], 2)) + "\n")
                log_f.write("{:>25} \t {:>25} \t {:>25} \t {:>25} \t {:>25} \n".format("Mass_of_isotopologue", 
                                                                                       "Isotopo_found?",
                                                                                       "measured_intensity",
                                                                                       "required_intensity",
                                                                                       "score_of_isotopologue"))
                for item in list(isotope_check.items()):
                    try:
                        log_f.write("{:>25} \t {:>25} \t {:>25} \t {:>25} \t {:>25} \n".format(str(round(item[0], 4)),
                                                                                                str(item[1][0]),
                                                                                                str(round(item[1][2], 4)),
                                                                                                str(round(item[1][1], 4)),
                                                                                                str(round(item[1][3], 4))))
                    except Exception as e:
                        continue
                log_f.close()
        if self.kwargs["include_likelyhood_of_formula"] == True:
            self.make_op_log_entry("")
            self.make_op_log_entry("=============================================================")
            self.make_op_log_entry("INFO:\t" + "Adding likelyhood of formula to the score...")
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
                    elif (dbe - f_dict.get("O", 0)) > 0:
                        score_substract = abs(dbe - f_dict.get("O", 0) - 1) * 50
                    else:
                        score_substract = 0
                    formula_score_dict[formula_dict_str] = float(score) - score_substract
                except Exception as e:
                    print("ERROR in formula score likelyhood: " + str(e))    
        self.make_op_log_entry("Finished adding likelyhood of formula to the score...")
        formula_score_dict = {f: s for f, s in formula_score_dict.items() if s > self.kwargs["reject_formula_if_score_lower_than"]}
        formula_score_dict = dict(sorted(formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        self.formulas_score_dict = formula_score_dict
        self.make_op_log_entry("")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICITIONARY:")
        for item in list(formula_score_dict.items()):
            self.make_op_log_entry("{:>30} \t {:>20}".format(str(item[0]), str(round(item[1], 2))))
        

        self.formulas_score_dict = formula_score_dict
        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, ((self.kwargs["mass_deviation"]*self.mass)/1000000), requested_filter_mode=self.spec.filter_mode)
        self.identified_peaks_for_mass = self.get_peak_properties(self.xic[0], self.xic[1])
        if self.kwargs["save_xic_plot"] == True:
            self.xic_plot_filepath = self.kwargs["pred_folder"] + "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["mass_deviation"]) + "_" + str( self.spec.filter_mode) + ".png"
            self.save_xic_plot(self.xic[0], self.xic[1], self.xic_plot_filepath)
        if self.kwargs["save_plot_of_isotopologues"] == True and len(list(formula_score_dict.keys())) > 0:
                
            print(list(formula_score_dict.keys())[0])
            print(ast.literal_eval(list(formula_score_dict.keys())[0]))
            self.make_plot_of_isotopologues(ast.literal_eval(list(formula_score_dict.keys())[0]))
            self.make_go_plot_of_isotopologues(ast.literal_eval(list(formula_score_dict.keys())[0]))

    def get_peak_properties(self, times, intensities):
        xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, ((self.kwargs["mass_deviation"]*self.mass)/1000000), requested_filter_mode=self.spec.filter_mode)
        window, order = 5, 3
        print(len(xic[1]))
        intensities = scipy.signal.savgol_filter(intensities, window, order, mode="nearest")
        print(len(intensities))
        peak_properties = scipy.signal.find_peaks(intensities, height=max(intensities)/100, distance=2, prominence=max(intensities)/50, width=(2, 20))
        print(peak_properties)
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

    def make_go_plot_of_isotopologues(self, formula_to_simulate, include_actual_spectral_data=True):
        isotope_simulation = MS_functions.simulate_isotope_pattern_of_formula(formula_to_simulate)
        fig = make_subplots(rows=1, cols=len(isotope_simulation), shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0, vertical_spacing=0)
        intensities_to_include = []
        for entry in list(isotope_simulation.keys()):
            summed_intensity = sum([self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if abs(((self.spec.summarized_masses[i] - entry) / entry)*1000000) <= (self.kwargs["mass_deviation"])])
            intensities_to_include.append(summed_intensity)
        intensities_to_include = [(i/(intensities_to_include[0]*list(isotope_simulation.keys())[0])) for i in intensities_to_include]

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
        filepath = self.kwargs["pred_folder"] + "isotopologues_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in formula_to_simulate.items()])) + ".html"
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
            within_deviation_mass_list = [masse for masse in list(self.spec.summarized_mass_intensity_dict.keys()) if abs((((isotopo_mass - masse) / isotopo_mass) * 1000000) - initial_deviation) < self.kwargs["max_ppm_deviation_change_for_isotopologue"]]
            measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
            deviation = (abs(isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000
            theoretical_intensity = (self.intensity_of_ion / (list(isotope_pattern_dict.items())[0][1])) * isotopo[1]
            noise = sorted(list(self.spec.summarized_mass_intensity_dict.values()))[int(len(list(self.spec.summarized_mass_intensity_dict.items())) / self.kwargs["noise_divisor_for_isotopologue_calculation"])] + self.kwargs["minimum_assumed_noise"]
            measured_intensity = 0

            if (deviation < self.kwargs["mass_deviation"]) and (abs(deviation - initial_deviation) <= self.kwargs["max_ppm_deviation_change_for_isotopologue"]):
                previous_deviation_list.append(copy.deepcopy(deviation))
                for masse in within_deviation_mass_list:
                    measured_intensity = measured_intensity + self.spec.summarized_mass_intensity_dict[masse]
                if measured_intensity <= 1:
                    break
                isotopes_found[isotopo[0]] = [True, theoretical_intensity, measured_intensity]
                score = measured_intensity / theoretical_intensity  # je naeher an 1 desto besser; wenn <1: weniger gemessen als da sein sollte; wenn >1: mehr gemessen als da sein sollte.
                if score > 1:
                    score = -(-1.5 + (1 / (0.7 + (2.718281828459045 ** (-self.kwargs["isotopologue_score_e_function_exponent"] * score)))))
                # score: je naeher an 1 desto besser; wenn negativ: weniger gemessen als theoretisch da; wenn positiv: mehr gemessen als theoretisch da.
                score = (abs(score * 100) + 0.001) + self.kwargs["add_value_to_score_if_isotopo_was_found"]
                if measured_intensity < noise:
                    score_multiplier_by_intensity_and_noise = abs(1 / (math.log(measured_intensity / noise) - 1))
                else:
                    score_multiplier_by_intensity_and_noise = abs(math.log(measured_intensity / noise) + 1)

                if isotopo == list(isotope_pattern_dict.items())[0]:
                    isotopes_found[isotopo[0]].append(0)
                    previous_score = 9999999
                    continue

                if theoretical_intensity <= (noise / 1.5):
                    score = (score * score_multiplier_by_intensity_and_noise) - abs(deviation * self.kwargs["multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    isotopes_found[isotopo[0]].append(score)
                    if (score - abs(self.kwargs["stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
                        isotopes_found[isotopo[0]][-1] = isotopes_found[isotopo[0]][-1] * 0.1
                    break

                else:
                    score = (score * score_multiplier_by_intensity_and_noise) - abs(deviation * self.kwargs["multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    isotopes_found[isotopo[0]].append(score)
                
                if (score < self.kwargs["stop_isotopologue_search_if_score_lower_than"]):
                    break

                if (score - abs(self.kwargs["stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
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
        if self.kwargs["save_detailed_log"] == False and error == False:
            return False
        #get the dirname of the spec_log_filepath
        directory_logfile = os.path.dirname(self.kwargs["op_log_filepath"])
        #create the directory if it does not exist
        os.makedirs(directory_logfile, exist_ok=True)
        log_f = open(self.kwargs["op_log_filepath"], "a")
        log_f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t" + log_entry + "\n")
        log_f.close()
        return True

    def make_plot_of_isotopologues(self, formula_to_simulate, include_actual_spectral_data=True):
        isotope_simulation = MS_functions.simulate_isotope_pattern_of_formula(formula_to_simulate)
        print(formula_to_simulate)
        print(isotope_simulation)
        self.make_op_log_entry("INFO:\t" + "Simulated isotope pattern will be plotted now...")
        if include_actual_spectral_data == True:
            self.make_op_log_entry("INFO:\t" + "Actual spectral data will be included!")
        self.make_op_log_entry("INFO:\t" + "Simulated isotope pattern: " + str(isotope_simulation))
        fig = plt.figure()
        gs = fig.add_gridspec(1, len(isotope_simulation), hspace=0, wspace=0)
        ax = gs.subplots(sharex="col", sharey="row")
        intensities_to_include = []
        for entry in list(isotope_simulation.keys()):
            summed_intensity = sum([self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if abs(((self.spec.summarized_masses[i] - entry) / entry)*1000000) <= (self.kwargs["mass_deviation"])])
            intensities_to_include.append(summed_intensity)
        intensities_to_include = [(i/(intensities_to_include[0]*list(isotope_simulation.keys())[0])) for i in intensities_to_include]
        
        for entry in range(len(intensities_to_include)):
            if include_actual_spectral_data:
                cmap = matplotlib.colormaps.get_cmap('Greens')
                colors = cmap((self.spec.summarized_intensities - min(self.spec.summarized_intensities)) / (max(self.spec.summarized_intensities) - min(self.spec.summarized_intensities)) * 0.3 + 0.7)                
                ax[entry].bar(self.spec.summarized_masses, self.spec.summarized_intensities, color=colors, label="all ions", width=0.0005, alpha=1)
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
        plt.legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        filepath = self.kwargs["pred_folder"] + "isotopologues_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in formula_to_simulate.items()])) + ".png"
        plt.savefig(filepath, bbox_inches='tight', dpi=1000)
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
        return True


class OneAnalysis:
    def __init__(self, ms_file, mass, rt, **kwargs):
        self.mass = mass
        self.rt = rt
        self.ms_file = ms_file
        self.peak_index = self.ms_file.rt_list.index(min(self.ms_file.rt_list, key=lambda x: abs(self.rt - x)))

        default_kwargs = {"peak_infos":[],  
                          "one_analysis_folder":str(self.ms_file.parentfolder + "/" + str(self.mass) + "_" + str(self.rt) + "/"),
                          "mass_deviation":11,
                          "mass_deviation_xic":(50*self.mass)/1000000,
                          "xic_requested_filter_mode":"Full scan",
                          "oa_log_filepath":str(self.ms_file.parentfolder + "/" + str(self.mass) + "_" + str(self.rt) + "/" + "oa_log.txt"),
                          "charge_of_measured_mass":-1,
                          "formula_cache_folder_path":"C://Users//Admin//Desktop//UVenture//Formula_Predictions//Formula_Predictions//",
                          "include_frag_intensity_noise_multiplier":0.01}
        self.kwargs = {**default_kwargs, **kwargs}

        os.makedirs(self.kwargs["one_analysis_folder"], exist_ok=True)

        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, self.kwargs["mass_deviation_xic"], requested_filter_mode=self.kwargs["xic_requested_filter_mode"])
        self.xic_plot_filepath = self.kwargs["one_analysis_folder"] + "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["mass_deviation_xic"]) + "_" + str(
                                self.kwargs["xic_requested_filter_mode"]) + ".png"
        self.plot_and_save_xic(self.xic_plot_filepath, self.peak_index, self.xic)
        
        if "Full scan" in self.ms_file.available_modes:
            self.full_scan_spec = Spec(self.ms_file, 
                                             self.peak_index, 
                                             requested_filter_mode="Full scan", 
                                             absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                             mass_deviation=self.kwargs["mass_deviation"])
        if "AIF" in self.ms_file.available_modes:
            self.aif_spec = Spec(self.ms_file, 
                                       self.peak_index, 
                                       requested_filter_mode="AIF", 
                                       absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                       mass_deviation=self.kwargs["mass_deviation"])
        if "MS/MS" in self.ms_file.available_modes:
            self.ms_ms_spec = Spec(self.ms_file, 
                                         self.peak_index, 
                                         requested_filter_mode="MS/MS", 
                                         absolute_spec_folder=self.kwargs["one_analysis_folder"] + "spectra/", 
                                         mass_deviation=self.kwargs["mass_deviation"], 
                                         requested_ms_ms_mass=self.mass)

        self.mass_old = self.mass
        self.mass = self.adjust_mass_to_closest_measured_mass(self.full_scan_spec)
        self.make_oa_log_entry("INFO:\t" + "Adjusting the mass of the molecular ion...")
        self.make_oa_log_entry("INFO:\t" + "New mass: " + str(self.mass) + "  Old mass: " + str(self.mass_old))
        self.mass = self.mass + (self.kwargs["charge_of_measured_mass"] * 0.000548)
        self.make_oa_log_entry("INFO:\t" + "Charge of the measured mass: " + str(self.kwargs["charge_of_measured_mass"]))
        self.make_oa_log_entry("INFO:\t" + "New mass after charge correction: " + str(self.mass))
        self.intensity_of_molecular_ion = sum([self.full_scan_spec.summarized_intensities[i] for i in range(len(self.full_scan_spec.summarized_intensities)) if abs(((self.full_scan_spec.summarized_masses[i] - self.mass)/self.mass)*1000000) <= self.kwargs["mass_deviation"]])
        self.make_oa_log_entry("INFO:\t" + "Intensity of the molecular ion: " + str(self.intensity_of_molecular_ion))

        self.molecular_ion_prediction = Prediction(self.ms_file, self.mass, self.full_scan_spec, 
                                                   absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/",
                                                   formula_cache_folder_path=self.kwargs["formula_cache_folder_path"],
                                                   save_plot_of_isotopologues=True,
                                                   save_xic_plot=True)
        self.make_oa_log_entry("INFO:\t" + "Finished prediction of molecular ion...")
        self.make_oa_log_entry("INFO:\t" + "Predicted molecular ion: " + str(self.molecular_ion_prediction.mass))
        self.make_oa_log_entry("INFO:\t" + "Ion intensity: " + str(self.molecular_ion_prediction.intensity_of_ion))
        self.possible_fragment_masses = [m for m in self.aif_spec.summarized_masses if m < self.mass and self.aif_spec.summarized_intensities[self.aif_spec.summarized_masses.index(m)] > (self.intensity_of_molecular_ion * self.kwargs["include_frag_intensity_noise_multiplier"])]
        print(self.possible_fragment_masses)
        self.fragment_predictions = {}
        for frag_mass in self.possible_fragment_masses:
            self.fragment_predictions[frag_mass] = Prediction(self.ms_file, frag_mass, self.aif_spec, 
                                                              absolute_pred_folder=self.kwargs["one_analysis_folder"] + "predictions/fragments/",
                                                              formula_cache_folder_path=self.kwargs["formula_cache_folder_path"],
                                                              save_plot_of_isotopologues=True,
                                                              save_xic_plot=True)
            plt.close("all")
            self.make_oa_log_entry("INFO:\t" + "Finished prediction of fragment ion with mass: " + str(frag_mass))
            self.make_oa_log_entry("INFO:\t" + "Fragment ion intensity: " + str(self.fragment_predictions[frag_mass].intensity_of_ion))
        

    def adjust_mass_to_closest_measured_mass(self, spec, masse=0):
        if masse == 0:
            masse = self.mass
        mass = min(list(spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(masse - x))
        return mass

    def plot_and_save_xic(self, save_filepath, peak_index, xic):
        fig, ax = plt.subplots()
        xic_peak_index = xic[2].index(min(xic[2], key=lambda x: abs(peak_index - x)))
        ax.scatter(self.ms_file.rt_list[peak_index], xic[1][xic_peak_index], color="red", marker="x", s=100)
        ax.plot(xic[0], xic[1], label="XIC measured")
        ax.set_title("xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["mass_deviation_xic"]) + "_" + str(self.kwargs["xic_requested_filter_mode"]))
        ax.set_xlabel("retention time / seconds")
        ax.set_ylabel("intensity / a.u.")
        plt.legend()
        matplotlib.rcParams.update({'figure.autolayout': True})
        plt.savefig(save_filepath, bbox_inches='tight', dpi=1000)
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




        
        
            
            
if __name__ == "__main__":
    mzml_filename = "C://Users//Admin//Desktop//UVenture//Evaluation_folder//F12_HRAIF4_3.mzML"

    ms_file = MS_File(mzml_filename)
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

    myanalysis = OneAnalysis(ms_file, 117.0554, ms_file.rt_list[79], mass_deviation=11, mass_deviation_xic=0.001, charge_of_measured_mass=-1, formula_cache_folder_path="C://Users//Admin//Desktop//UVenture//Formula_Predictions//Formula_Predictions//")

    print(myanalysis.molecular_ion_prediction.best_formula_prediction)
    print(myanalysis.molecular_ion_prediction.score_of_best_formula)
    print(myanalysis.molecular_ion_prediction.intensity_of_ion)


    #print(myprediction.mass_possible_formulas)
    #print(myprediction.formulas_score_dict)
    #print(myprediction.mass)
    #print(myprediction.intensity_of_ion)


    







