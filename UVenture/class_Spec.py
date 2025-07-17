import copy
import datetime
import os
import numpy as np

import UVenture.plotting as plotting
import UVenture.MS_functions as MS_functions


class Spec:
    def __init__(self, ms_file, index, debug_output=False, **kwargs):
        self.debug_output = debug_output
        self.original_index = index
        self.ms_file = ms_file

        requested_mode = kwargs.get("spec_requested_filter_mode", "Full scan")
        requested_ms_ms_mass = kwargs.get("spec_requested_ms_ms_mass", "")
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
        self.summarized_masses = []
        self.summarized_intensities = []
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

        if self.debug_output:
            print("Getting mass spec...")
        self.summarized_mass_intensity_dict = MS_functions.summarize_mass_intensity_dict(self.get_mass_spec(self.index), deviation=self.kwargs["mass_deviation"], debug_output=False)
        self.summarized_masses = list(self.summarized_mass_intensity_dict.keys())
        self.summarized_intensities = list(self.summarized_mass_intensity_dict.values())
        if self.debug_output:
            print("Summarized mass intensity dict shortened to: " + str(len(self.summarized_mass_intensity_dict)) + " elements.")
        self.make_spec_log_entry("Summarized_mass_intensity_dict: ")
        log_summarized_mass_intensity_dict = {float(dictkey) if isinstance(dictkey, np.float64) else dictkey : float(dictvalue) if isinstance(dictvalue, np.float32) else dictvalue for dictkey, dictvalue in self.summarized_mass_intensity_dict.items()}
        self.make_spec_log_entry(str(log_summarized_mass_intensity_dict))

        # Create a mass spectrum plot using go
        if self.kwargs["spec_save_go_plot"] == True:
            plot_filepath = self.kwargs["spec_folder"] + "Mass_spectrum_index" + str(self.index) + "_" + str(self.filter_mode.replace("/", "")) + "_" + str(self.ms_ms_masses) + ".html"
            self.create_barchart_massspec_with_go(self.summarized_masses, self.summarized_intensities, plot_filepath=plot_filepath)

        # Create a mass spectrum plot using matplotlib
        if self.kwargs["spec_save_matplotlib_plot"] == True:
            self.make_spec_log_entry("INFO:\t" + "Saving mass spectrum plot...")
            title = "Mode:" + str(self.filter_mode) + "; Index: " + str(self.index) + "; RT: " + str(round(self.ms_file.rt_list[self.index], 2)) + ";\nFilter: " + str(filter)
            image_filepath = self.kwargs["spec_folder"] + "Mass_spectrum_index" + str(self.index) + "_" + str(self.filter_mode.replace("/", "")) + ".png"
            plotting.create_barchart_massspec(self.summarized_masses, self.summarized_intensities, title=title, image_filepath=image_filepath)

        self.make_spec_log_entry("INFO:\t" + "Found MS/MS masses: " + str(self.ms_ms_masses))
        self.make_spec_log_entry("INFO:\t" + "Filter mode: " + str(self.filter_mode))
        self.make_spec_log_entry("INFO:\t" + "Filter: " + str(self.filter))
        self.make_spec_log_entry("INFO:\t" + "Length of summarized mass intensity dict: " + str(len(self.summarized_mass_intensity_dict)) + " elements.")
        self.make_spec_log_entry("INFO:\t" + "TIC: " + str(self.ms_file.tic[self.index]))

    def make_spec_log_entry(self, log_entry):
        #get the dirname of the spec_log_filepath
        directory_logfile = os.path.dirname(self.kwargs["spec_log_filepath"])
        #create the directory if it does not exist
        os.makedirs(directory_logfile, exist_ok=True)
        log_f = open(self.kwargs["spec_log_filepath"], "a")
        log_f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t" + log_entry + "\n")
        log_f.close()
        return True

    def search_for_required_spec_close_to_rt(self, index, requested_mode="whatever", requested_ms_ms_mass=""):
        curr_index = copy.deepcopy(index)
        change_index_value = 1
        change_index_direction = "+"
        start_index = copy.deepcopy(index)
        if self.debug_output:
            print("Starting ms spectrum search at index: " + str(curr_index))

        while True and not requested_mode == "whatever":
            filter = self.ms_file.rawdata[curr_index]["scanList"]["scan"][0]["filter string"]
            filter_mode = MS_functions.get_mode_of_spec(filter)
            if filter_mode == "MS/MS":
                ms_ms_masses = []
                filter_parsed = filter.split(" ")
                filter_parsed = [x for x in filter_parsed if "hcd" in x]
                for element in filter_parsed:
                    ms_ms_masses.append(round(float(element.split("@")[0]), 1))

            if (curr_index >= len(self.ms_file.rawdata)-3) or (curr_index <= 3):
                if self.debug_output:
                    print("NOTHING WAS FOUND IN THE WHOLE CHROMATOGRAM! Index at boundaries. Returning...")
                    print("Returning the start index: " + str(start_index))
                    print("Mode of the start index: " + str(MS_functions.get_mode_of_spec(self.ms_file.rawdata[start_index]["scanList"]["scan"][0]["filter string"])))
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
                if self.debug_output:
                    print("Found requested spectrum")
                break

            elif (filter_mode == requested_mode) and (round(requested_ms_ms_mass, 1) in ms_ms_masses):
                if self.debug_output:
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
        filter_mode = MS_functions.get_mode_of_spec(filter)
        if self.debug_output:
            print("Filter: " + filter)
            print("MS level of the mass spectrum: " + str(self.ms_file.rawdata[index]["ms level"]))
        ms_ms_masses = []
        filter_parsed = filter.split(" ")
        filter_parsed = [x for x in filter_parsed if "hcd" in x]
        for element in filter_parsed:
            ms_ms_masses.append((float(element.split("@")[0])))

        self.filter = filter
        self.filter_mode = filter_mode
        self.ms_ms_masses = ms_ms_masses
        self.spec_rawdata = self.ms_file.rawdata[index]
        self.ms_level = self.ms_file.rawdata[index]["ms level"]
        return dict(zip(masses, intensities))