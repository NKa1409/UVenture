import ast
import copy
import datetime
import math
import traceback
import scipy
import os
import shutil
import numpy as np

import UVenture.formula_cache as formula_cache
import UVenture.chemical_formula_parsing as chemical_formula_parsing
import UVenture.MS_functions as MS_functions
import UVenture.plotting as plotting
import UVenture.class_Spec as class_Spec
import UVenture.class_Prediction as class_Prediction


class OneAnalysis:
    def __init__(self, ms_file, mass, rt, **kwargs):
        self.debug_output = True
        self.mass = mass
        self.rt = rt
        self.ms_file = ms_file
        self.unit_of_rt = "sec"  # default unit of retention time

        self.peak_index = self.ms_file.rt_list.index(min(self.ms_file.rt_list, key=lambda x: abs(self.rt - x)))

        default_kwargs = {"one_analysis_folder":os.path.join(self.ms_file.parentfolder, str(self.mass) + "_" + str(self.rt)),
                          "oa_log_filepath":os.path.join(self.ms_file.parentfolder, str(self.mass) + "_" + str(self.rt), "oa_log.txt"),
                          "mass_deviation":11,
                          "charge_of_measured_mass":-1,

                          "oa_save_xic_plot":True,
                          "oa_xic_requested_filter_mode":"Full scan",

                          "oa_stop_if_oa_given_ion_is_a_fragment": True,

                          "oa_spec_acquisition_save_matplotlib_plot":True,
                          "oa_spec_acquisition_save_go_plot":True,
                          "oa_spec_acquisition_save_additional_info":True,

                          "oa_molecular_ion_pred_save_matplotlib_plot_of_isotopologues":True,
                          "oa_molecular_ion_pred_save_go_plot_of_isotopologues":True,
                          "oa_molecular_ion_pred_save_xic_plot":True,
                          "oa_molecular_ion_pred_save_detailed_log":True,
                          "oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found":False,

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

        if self.debug_output:
            print()
            print("OneAnalysis kwargs: " + str(self.kwargs))
            print()
            print("Make good fragment prediction: " + str(self.kwargs["oa_make_good_fragment_formula_prediction"]))
        self.kwargs["oa_make_good_fragment_formula_prediction"] = bool(self.kwargs["oa_make_good_fragment_formula_prediction"])

        additional_kwargs = copy.deepcopy(self.kwargs)
        pop_keys = ["spec_requested_filter_mode", "absolute_spec_folder", "mass_deviation", "spec_save_matplotlib_plot", "spec_save_go_plot"]
        for key in pop_keys:
            try:
                additional_kwargs.pop(key)
            except KeyError:
                continue
        
        if self.debug_output:
            print("Additional Kwargs: " + str(additional_kwargs))


        self.make_oa_log_entry("")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("=============================================================")
        self.make_oa_log_entry("INFO:\t" + "Creating OneAnalysis object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.peak_index))
        self.make_oa_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))
        self.make_oa_log_entry("INFO:\t" + "Assuming a charge of the measured mass of: " + str(self.kwargs["charge_of_measured_mass"]))
        self.make_oa_log_entry("INFO:\t" + "Available MS modes: " + str(self.ms_file.available_modes))

        self.xic = MS_functions.get_xic(self.ms_file.rawdata, self.mass, round(((self.kwargs["mass_deviation"]*self.mass) / 1000000), 4), requested_filter_mode=self.kwargs["oa_xic_requested_filter_mode"])
        self.make_oa_log_entry("INFO:\t" + "XIC calculated. Continuing...")

        self.make_oa_log_entry("INFO:\t" + "Adjusting retention time....")
        new_rt = self.adjust_retention_time_to_peak_maximum(original_rt=self.rt, max_rt_shift=20, xic=self.xic)
        self.make_oa_log_entry("INFO:\t" + "Retention time adjusted to: " + str(new_rt))
        self.rt = new_rt
        self.peak_index = self.ms_file.rt_list.index(min(self.ms_file.rt_list, key=lambda x: abs(self.rt - x)))

        if self.kwargs["oa_save_xic_plot"] == True:
            print("Saving XIC plot...")
            self.xic_plot_filepath = os.path.join(self.kwargs["one_analysis_folder"], "xic_" + str(round(self.mass, 4)) + "+-" + str( round(((self.kwargs["mass_deviation"]*self.mass) / 1000000), 4) ) + "_" + str(
                                    self.kwargs["oa_xic_requested_filter_mode"].replace(" ", "_")) + ".png")
            title = "Extracted Ion Chromatogram (XIC) for mass: " + str(round(self.mass, 4)) + " at retention time: " + str(round(self.rt, 4)) + " sec"
            plotting.create_xic(self.xic[0], self.xic[1], title=title, filepath=self.xic_plot_filepath, retention_time=self.rt)
        
        self.make_oa_log_entry("INFO:\t" + "Extracted Ion Chromatogram (XIC) plot saved at: " + str(self.xic_plot_filepath))


        if "Full scan" in self.ms_file.available_modes:
            self.full_scan_spec = class_Spec.Spec(self.ms_file, 
                                             self.peak_index, 
                                             spec_requested_filter_mode="Full scan", 
                                             absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"), 
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
            self.aif_spec = class_Spec.Spec(self.ms_file, 
                                       self.peak_index, 
                                       spec_requested_filter_mode="AIF", 
                                       absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"), 
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
            self.ms_ms_spec = class_Spec.Spec(self.ms_file, 
                                         self.peak_index, 
                                         spec_requested_filter_mode="MS/MS", 
                                         absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"), 
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

        intensity_ratio, int_in_mi, int_in_frag, type_of_ion = self.get_most_likely_type_of_ion(self.mass, self.best_molecular_ion_spec, self.best_frag_spec)
        self.make_oa_log_entry("INFO:\t" + "Type of the analyzed ion is most likely to be: " + str(type_of_ion))
        self.make_oa_log_entry("INFO:\t" + "Intensity of the ion in best molecular ion spec and best frag spec: " + str(int_in_mi) + " / " + str(int_in_frag))
        if type_of_ion == "fragment_ion" and self.kwargs["oa_stop_if_oa_given_ion_is_a_fragment"] == True:
            self.make_oa_log_entry("INFO:\t" + "Stopping the prediction for the molecular ion, as the given mass is most likely to be a fragment.")
            if self.kwargs["oa_zip_folder_when_finished"] == True:
                shutil.make_archive(os.path.join(self.ms_file.parentfolder, str(round(self.mass, 5)) + "_" + str(round(self.rt, 3))), "zip", self.kwargs["one_analysis_folder"])
                shutil.rmtree(self.kwargs["one_analysis_folder"])
            return


        self.molecular_ion_prediction = self.get_molecular_ion_prediction(self.best_molecular_ion_spec)
        if self.molecular_ion_prediction.peak_found == False and self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"] == True:
            self.make_oa_log_entry("INFO:\t" + "No molecular ion prediction peak found. Stopping prediction...")
            if self.kwargs["oa_zip_folder_when_finished"] == True:
                shutil.make_archive(os.path.join(self.ms_file.parentfolder, str(round(self.mass, 5)) + "_" + str(round(self.rt, 3))), "zip", self.kwargs["one_analysis_folder"])
                shutil.rmtree(self.kwargs["one_analysis_folder"])
            return

        if len(list(self.summarized_molecular_ion_formula_score_dict.keys())) == 0:
            print("No molecular ion prediction could be found! Returning....")
            self.make_oa_log_entry("INFO:\t" + "No molecular ion prediction could be found! Returning.....")
            if self.kwargs["oa_zip_folder_when_finished"] == True:
                shutil.make_archive(os.path.join(self.ms_file.parentfolder, str(round(self.mass, 5)) + "_" + str(round(self.rt, 3))), "zip", self.kwargs["one_analysis_folder"])
                shutil.rmtree(self.kwargs["one_analysis_folder"])
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

        for i in range(len(write_to_summary_file_list)):
            if isinstance(write_to_summary_file_list[i], np.float64) or isinstance(write_to_summary_file_list[i], np.float32):
                write_to_summary_file_list[i] = float(write_to_summary_file_list[i])
            elif isinstance(write_to_summary_file_list[i], list):
                for k in range(len(write_to_summary_file_list[i])):
                    if isinstance(write_to_summary_file_list[i][k], np.float64) or isinstance(write_to_summary_file_list[i][k], np.float32):
                        write_to_summary_file_list[i][k] = float(write_to_summary_file_list[i][k])
        self.append_oa_summary_to_raw_file_summary(write_to_summary_file_list)

        self.make_oa_log_entry("INFO:\t" + "Finished appending summary to raw file summary...")
        self.make_oa_log_entry("INFO:\t" + "Creating summary plot and txt file...")
        try:
            save_filepath = os.path.join(self.kwargs["one_analysis_folder"], "oa_summary_plot.png")
            plotting.create_oa_summary_plot(self, self.true_fragment_list, self.best_frag_spec, self.fragment_predictions, save_filepath)
        except Exception as e:
            self.make_oa_log_entry("ERROR:\t" + "Error while creating summary plot: " + str(e))
            print("Error while creating summary plot: " + str(e))
            print(traceback.format_exc())
        self.make_oa_log_entry("INFO:\t" + "Summary plot saved at: " + str(save_filepath))
        oa_txt_save_filepath = os.path.join(self.kwargs["one_analysis_folder"], "BEST_FORMULA_PREDICTION.txt")
        self.create_oa_summary_txtfile(self.true_fragment_list, self.fragment_predictions, oa_txt_save_filepath)

        if self.kwargs["oa_zip_folder_when_finished"] == True:
            shutil.make_archive(os.path.join(self.ms_file.parentfolder, str(round(self.mass, 5)) + "_" + str(round(self.rt, 3))), "zip", self.kwargs["one_analysis_folder"])
            shutil.rmtree(self.kwargs["one_analysis_folder"])

    def adjust_retention_time_to_peak_maximum(self, original_rt, max_rt_shift, xic):
        try:
            index_window = int(((max_rt_shift * (max(xic[0]) / len(xic[0])))+1) / 2)
            print(index_window)
            peak_index = xic[0].index(min(xic[0], key=lambda x: abs(original_rt - x)))
            rt_window, int_window = xic[0][peak_index-index_window:peak_index+index_window], xic[1][peak_index-index_window:peak_index+index_window]
            print(rt_window)
            print(int_window)
            peak_avg_int_dict = {}
            for i in range(len(int_window)):
                try:
                    tripple_summed_int = int_window[i-1] + int_window[i] + int_window[i+1]
                except IndexError:
                    continue
                peak_avg_int_dict[rt_window[i]] = tripple_summed_int
            avg_peak_maximum_rt = max(peak_avg_int_dict, key=peak_avg_int_dict.get)
            avg_peak_maximum_index = rt_window.index(avg_peak_maximum_rt)
            avg_ints = [peak_avg_int_dict[i] for i in peak_avg_int_dict.keys()]
            avg_rts = [i for i in peak_avg_int_dict.keys()]

            try:
                smooth_intensities = scipy.signal.savgol_filter(avg_ints, int(len(avg_ints) / 5), 2)
            except Exception as e:
                print("Error while smoothing intensities: " + str(e))
                print(traceback.format_exc())
                smooth_intensities = list(peak_avg_int_dict.values())
            new_indices_of_identified_peaks = [i for i in range(1, len(smooth_intensities) - 1) if smooth_intensities[i - 1] < smooth_intensities[i] > smooth_intensities[i + 1]]
            if len(new_indices_of_identified_peaks) >= 2:
                peak_index_in_smooth = min(new_indices_of_identified_peaks, key=lambda x: abs(peak_index - x))
            new_rt = avg_rts[peak_index_in_smooth]
            if new_rt == avg_peak_maximum_rt:
                print("Peak was adjusted to the average peak maximum retention time: " + str(new_rt))
                print("It is likely a good adjustment.")
                return new_rt
            else:
                print("Peak was adjusted to the average peak maximum retention time: " + str(new_rt))
                print("It is likely not a good adjustment as multiple peaks were found in the observation window.")
                print("Adjusted the peak retention time to the identified peak that is closest to the original peak index.")
                print("Original peak retention time: " + str(original_rt))
                print("Adjusted peak retention time: " + str(new_rt))
                return new_rt

        except:
            return original_rt


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
        self.spec_before = class_Spec.Spec(self.ms_file,
                                index_before,
                                spec_requested_filter_mode=self.best_frag_spec.filter_mode,
                                absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"),
                                mass_deviation=self.kwargs["mass_deviation"],
                                spec_save_matplotlib_plot=False,
                                spec_save_go_plot=False,
                                **additional_kwargs)
        self.spec_after = class_Spec.Spec(self.ms_file,
                               index_after,
                               spec_requested_filter_mode=self.best_frag_spec.filter_mode,
                               absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"),
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
        self.spec_before = class_Spec.Spec(self.ms_file,
                                index_before,
                                spec_requested_filter_mode=self.best_molecular_ion_spec.filter_mode,
                                absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"),
                                mass_deviation=self.kwargs["mass_deviation"],
                                spec_save_matplotlib_plot=False,
                                spec_save_go_plot=False,
                                **additional_kwargs)
        self.spec_after = class_Spec.Spec(self.ms_file,
                               index_after,
                               spec_requested_filter_mode=self.best_molecular_ion_spec.filter_mode,
                               absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"),
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

    def append_oa_summary_to_raw_file_summary(self, summary_list):
        self.make_oa_log_entry("INFO:\t" + "Appending summary to raw file summary. Summary list: " + str(summary_list))
        with open(os.path.join(self.ms_file.parentfolder, "SUMMARY.txt"), "a") as oa_summary:
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
        molecular_ion_best_approx_dict = chemical_formula_parsing.get_formula_to_dict(molecular_ion_best_approx)
        molecular_ion_score = self.score_of_best_molecular_ion_prediction
        molecular_ion_intensity = self.intensity_of_molecular_ion
        molecular_ion_formula_score_dict = self.summarized_molecular_ion_formula_score_dict # {"{"C": 1, "H": 2, ....}": score, "{}": score2}
        self.make_oa_log_entry("INFO:\t" + "Molecular ion formula score dict: " + str(molecular_ion_formula_score_dict))
        self.make_oa_log_entry("INFO:\t" + "Fragment ion formula score dict: " + str(self.fragment_predictions_formula_score_dicts))
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
                neutral_loss_formula_predictions = (formula_cache.get_formula_from_cache(self.kwargs["pred_formula_cache_folder_path"], neutral_loss_mass, self.kwargs["mass_deviation"]))[1]
                # neutral_loss_formula_prediction = [mass, alldict] ---> alldict = {formula: score, "C1H3O2": score, ...}
                print(neutral_loss_formula_predictions)
                if len(neutral_loss_formula_predictions) == 0:
                    continue
                found_pair = False
                for f_formula, f_score in f_score_dict.items():
                    f_formula_dict = ast.literal_eval(f_formula)
                    for nl_formula, nl_deviation in neutral_loss_formula_predictions.items():
                        nl_formula_dict = chemical_formula_parsing.get_formula_to_dict(nl_formula)
                        print(nl_formula_dict)
                        print(f_formula_dict)
                        print(molecular_ion_pred_dict)
                        if (self.combine_and_sum_dicts(nl_formula_dict, f_formula_dict) == molecular_ion_pred_dict):
                            print("TRUETRUETRUEskjaskjhlgfaivwzbevwuief")
                            self.make_oa_log_entry("INFO:\t" + "Found matching fragment: " + str(f_formula_dict))
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

        self.molecular_ion_prediction = class_Prediction.Prediction(self.ms_file, self.mass, best_molecular_ion_spec, spec_before=self.mi_spec_before, spec_after=self.mi_spec_after,
                                                   absolute_pred_folder=os.path.join(self.kwargs["one_analysis_folder"], "predictions"),
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
            spec_before = class_Spec.Spec(self.ms_file, 
                               self.peak_index - len(self.ms_file.available_modes), 
                               spec_requested_filter_mode=best_molecular_ion_spec.filter_mode, 
                               absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"),
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
            spec_after = class_Spec.Spec(self.ms_file, 
                              self.peak_index + len(self.ms_file.available_modes), 
                              spec_requested_filter_mode=best_molecular_ion_spec.filter_mode, 
                              absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"),
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
                class_Prediction.Prediction(self.ms_file, self.mass, curr_spec,
                                    absolute_pred_folder=os.path.join(self.kwargs["one_analysis_folder"], "predictions"),
                                    pred_formula_cache_folder_path=self.kwargs["pred_formula_cache_folder_path"],
                                    pred_save_matplotlib_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_matplotlib_plot_of_isotopologues"],
                                    pred_save_go_plot_of_isotopologues=self.kwargs["oa_molecular_ion_pred_save_go_plot_of_isotopologues"],
                                    pred_save_xic_plot=self.kwargs["oa_molecular_ion_pred_save_xic_plot"],
                                    pred_save_detailed_log=self.kwargs["oa_molecular_ion_pred_save_detailed_log"],
                                    pred_return_if_no_peak_is_found=self.kwargs["oa_molecular_ion_only_calc_prediction_if_molecular_ion_peak_is_found"],
                                    **additional_kwargs)
            self.make_oa_log_entry("INFO:\t" + "Old prediction did not match with the prediction of multiple spectra!")
            self.make_oa_log_entry("INFO:\t" + "Plots are created for every spectrum next to the original index. Plots will be available.")
            with open(os.path.join(self.kwargs["one_analysis_folder"], "predictions", "BEST_FORMULA_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0]).items()])) + ".txt"), "a") as txt_file:
                txt_file.write("Best prediction according to three spectra which are located next to the original peak: " + str(ast.literal_eval(list(self.summarized_molecular_ion_formula_score_dict.keys())[0])))
                txt_file.write("\n")
                log_summarized_molecular_ion_formula_score_dict = {float(dictkey) if isinstance(dictkey, np.float64) else dictkey : float(dictvalue) if isinstance(dictvalue, np.float64) else dictvalue for dictkey, dictvalue in self.summarized_molecular_ion_formula_score_dict.items()}
                txt_file.write("Summarized formula score dict: " + str(log_summarized_molecular_ion_formula_score_dict))
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

        self.summarized_molecular_ion_formula_score_dict = {f: s for f, s in self.summarized_molecular_ion_formula_score_dict.items() if s > (self.kwargs["oa_reject_formula_if_score_lower_than"]/5)}
        self.summarized_molecular_ion_formula_score_dict = dict(sorted(self.summarized_molecular_ion_formula_score_dict.items(), key=lambda x: x[1], reverse=True))
        log_summarized_molecular_ion_formula_score_dict = {float(dictkey) if isinstance(dictkey, np.float64) else dictkey: float(dictvalue) if isinstance(dictvalue, np.float64) else dictvalue for dictkey, dictvalue in self.summarized_molecular_ion_formula_score_dict.items()}
        self.make_oa_log_entry("INFO:\t" + "Summarized formula score dict: " + str(log_summarized_molecular_ion_formula_score_dict))
        print("Summarized formula score dict: " + str(log_summarized_molecular_ion_formula_score_dict))

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
            if e == "list index out of range":
                print("No formula found! So no best formula prediction could be chosen!")
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Formula predicted with " + str(len(available_specs)) + " spectra.")
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "All formula score dicts: ")
        for i, f_score_dict in enumerate(all_formula_score_dicts):
            log_f_score_dict = {chemical_formula_parsing.get_formula_string_from_dict(dictkey) if isinstance(dictkey, dict) else dictkey: float(dictvalue) if isinstance(dictvalue, np.float64) else dictvalue for dictkey, dictvalue in f_score_dict.items()}
            self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Formula score dict " + str(i) + ": " + str(log_f_score_dict))
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Best formula prediction: " + str(self.best_molecular_ion_prediction) + " Score: " + str(self.score_of_best_molecular_ion_prediction))
        log_summarized_molecular_ion_formula_score_dict = {chemical_formula_parsing.get_formula_string_from_dict(dictkey) if isinstance(dictkey, dict) else dictkey: float(dictvalue) if isinstance(dictvalue, np.float64) else dictvalue for dictkey, dictvalue in self.summarized_molecular_ion_formula_score_dict.items()}
        self.molecular_ion_prediction.make_op_log_entry("INFO:\t" + "Summarized formula score dict: " + str(log_summarized_molecular_ion_formula_score_dict))
        self.molecular_ion_prediction.make_op_log_entry("=============================================================")
        self.make_oa_log_entry("INFO:\t" + "Molecular Ion Prediction finished. \nBest formula prediction: " + str(self.best_molecular_ion_prediction) + " single score: " + str(self.score_of_best_molecular_ion_prediction))
        return self.molecular_ion_prediction
    
    def get_formula_score_dict_with_multiple_specs(self, specs, mass, return_if_no_peak_found=True, absolute_pred_subfolder_path="predictions"):
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
                pred = class_Prediction.Prediction(self.ms_file, mass, spec,
                                  absolute_pred_folder=os.path.join(self.kwargs["one_analysis_folder"], absolute_pred_subfolder_path),
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
                spec_before = class_Spec.Spec(self.ms_file, 
                                   self.peak_index - len(self.ms_file.available_modes), 
                                   spec_requested_filter_mode=self.best_frag_spec.filter_mode, 
                                   absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"), 
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
                spec_after = class_Spec.Spec(self.ms_file, 
                                  self.peak_index + len(self.ms_file.available_modes), 
                                  spec_requested_filter_mode=self.best_frag_spec.filter_mode, 
                                  absolute_spec_folder=os.path.join(self.kwargs["one_analysis_folder"], "spectra"), 
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

        self.possible_fragment_masses = [m for m in self.best_frag_spec.summarized_masses if m < self.mass-0.1 and self.best_frag_spec.summarized_intensities[self.best_frag_spec.summarized_masses.index(m)] > (self.intensity_of_molecular_ion * self.kwargs["oa_include_frag_intensity_noise_multiplier"])]
        self.possible_fragment_masses = [float(m) for m in self.possible_fragment_masses if isinstance(m, np.float64)]
        print(self.possible_fragment_masses)
        self.possible_fragment_masses = sorted(self.possible_fragment_masses, key=lambda m: self.best_frag_spec.summarized_intensities[self.best_frag_spec.summarized_masses.index(m)], reverse=True)
        print(self.possible_fragment_masses)
        self.possible_fragment_masses = [self.possible_fragment_masses[m] for m in range(len(self.possible_fragment_masses)) if m <= self.kwargs["oa_fragments_absolute_max_number_of_fragment_masses"]]
        log_possible_fragment_masses = [float(m) for m in self.possible_fragment_masses if isinstance(m, np.float64)]
        print("Possible fragment masses: " + str(log_possible_fragment_masses))
        self.make_oa_log_entry("INFO:\t" + "Possible fragment masses: " + str(log_possible_fragment_masses))
        self.fragment_predictions = {}
        self.fragment_predictions_formula_score_dicts = {}

        while len(self.possible_fragment_masses) >= 1:
            frag_mass = self.possible_fragment_masses[0]
            self.possible_fragment_masses.pop(self.possible_fragment_masses.index(frag_mass))

            try:
                if self.kwargs["oa_fragments_do_good_peak_comparison_with_area_between_curves"]:
                    os.makedirs(os.path.join(self.kwargs["one_analysis_folder"], "predictions", "fragments", "peak_matching"), exist_ok=True)
                    fragment_xic = MS_functions.get_xic(self.ms_file.rawdata, frag_mass, round(((self.kwargs["mass_deviation"]*self.mass) / 1000000), 4), requested_filter_mode=self.best_frag_spec.filter_mode)

                    print("Fragment mass:   " + str(frag_mass))

                    frag_type, ratio_frag_mi_int, mi_spec_frag_int_normalized, frag_spec_frag_int_normalized = self.get_likely_type_of_frag_mass_with_respect_to_mi_mass(self.mass, frag_mass, self.best_molecular_ion_spec, self.best_frag_spec)
                    if not frag_type == "fragment_ion":
                        print("Normalized intensity of the fragment is too low in relation to the molecular ion. The current mass cannot be a fragment mass. Frag mass: " + str(frag_mass))
                        print("Ratio Fragment/MI normalized intensity: " + str(ratio_frag_mi_int))
                        self.make_oa_log_entry("INFO:\t" + "Normalized intensity of the fragment is too low in relation to the molecular ion. The current mass cannot be a fragment mass. Frag mass: " + str(frag_mass))
                        #continue
                    else:
                        print("Current fragment mass can be a fragment of the selected molecular ion. Normalized intensity ratio looks good: " + str(ratio_frag_mi_int))

                    ratio_in_spec, int_mi_spec, int_frag_spec, type_of_ion = self.get_most_likely_type_of_ion(frag_mass, self.best_molecular_ion_spec, self.best_frag_spec)
                    if type_of_ion == "fragment_ion": #mass is a fragment of some sort
                        print("Frag mass can be a fragment: Intensity ratio of the potential fragment mass in mi_spec and frag_spec: " + str(ratio_in_spec))
                    else:
                        print("Frag mass cannot be a fragment. Intensity ratio does not match: " + str(ratio_in_spec))
                        #continue

                    area_between_curves, peak1_rt, peakintensity1, peak2_rt, peakintensity2 = MS_functions.compare_peak_shape_similarity(xic1=self.xic, xic2=fragment_xic, peak_rt=self.rt, debug_output=True, peakwidth=10)
                    if area_between_curves < 0:
                        print("Error in calculating area between the two curves.")
                        print("Peakintensity1: " + str(peakintensity1))
                        print("Peakintensity2: " + str(peakintensity2))
                        print("Peak RT 1: " + str(peak1_rt))
                        print("Peak RT 2: " + str(peak2_rt))
                        continue

                    if area_between_curves > self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"]:
                        print("No good fragment peak shape was detected. Continuing with the next fragment...")
                        self.make_oa_log_entry("INFO:\t" + "Fragment peak with mass: " + str(frag_mass) + "  -> Does not have a good fragment peak shape. Area between curves too high: " + str(area_between_curves) + " Continuing....")
                        title = "Peak matching evaluation" + str(round(self.mass, 4)) + " / " + str(round(frag_mass, 4))
                        save_filepath = os.path.join(self.kwargs["one_analysis_folder"], "predictions", "fragments", "peak_matching", "FALSE_Fragmass_" + str(round(frag_mass, 4)) + ".png")
                        plotting.create_xic_matching_plot(peak1_rt, peakintensity1, peak2_rt, peakintensity2, title, area_between_curves, save_filepath)
                        continue
                    else:
                        title = "Peak matching evaluation" + str(round(self.mass, 4)) + " / " + str(round(frag_mass, 4))
                        save_filepath = os.path.join(self.kwargs["one_analysis_folder"], "predictions", "fragments", "peak_matching", "TRUE_Fragmass_" + str(round(frag_mass, 4)) + ".png")
                        plotting.create_xic_matching_plot(peak1_rt, peakintensity1, peak2_rt, peakintensity2, title, area_between_curves, save_filepath)
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

                curr_prediction = class_Prediction.Prediction(self.ms_file, frag_mass, self.best_frag_spec, spec_before=self.frag_spec_before, spec_after=self.frag_spec_after,
                                                                  absolute_pred_folder=os.path.join(self.kwargs["one_analysis_folder"], "predictions", "fragments"),
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


    def get_likely_type_of_frag_mass_with_respect_to_mi_mass(self, mi_mass, frag_mass, mi_spec, frag_spec):
        closest_mi_mass_mi_spec = min(list(mi_spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(((mi_mass - x) / mi_mass) * 1000000))
        mi_intensity_mi_spec = mi_spec.summarized_mass_intensity_dict[closest_mi_mass_mi_spec]

        closest_mi_mass_frag_spec = min(list(frag_spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(((mi_mass - x) / mi_mass) * 1000000))
        mi_intensity_frag_spec = frag_spec.summarized_mass_intensity_dict[closest_mi_mass_frag_spec]

        closest_frag_mass_mi_spec = min(list(mi_spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(((frag_mass - x) / frag_mass) * 1000000))
        frag_intensity_mi_spec = mi_spec.summarized_mass_intensity_dict[closest_frag_mass_mi_spec]

        closest_frag_mass_frag_spec = min(list(frag_spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(((frag_mass - x) / frag_mass) * 1000000))
        frag_intensity_frag_spec = frag_spec.summarized_mass_intensity_dict[closest_frag_mass_frag_spec]
        print("MI Spec: MI_int: " + str(mi_intensity_mi_spec) + "   Frag_int: " + str(frag_intensity_mi_spec))
        mi_spec_frag_int_normalized = (frag_intensity_mi_spec / mi_intensity_mi_spec)
        frag_spec_frag_int_normalized = (frag_intensity_frag_spec / mi_intensity_frag_spec)
        ratio = (frag_spec_frag_int_normalized) / (mi_spec_frag_int_normalized) # if lower than 1: highly unlikely to be a true fragment.; If higher than 1 likely to be a fragment.
        print("MI Spec Frag Intensity Normalized: " + str(mi_spec_frag_int_normalized))
        print("Frag Spec: MI_int: " + str(mi_intensity_frag_spec) + "   Frag_int: " + str(frag_intensity_frag_spec))
        print("Frag Spec Frag Intensity Normalized: " + str(frag_spec_frag_int_normalized))
        print("Intensity Ratio: " + str(ratio))
        if ratio >= 1:
            frag_type = "fragment_ion"
        else:
            frag_type = "non_fragment"
        return frag_type, ratio, mi_spec_frag_int_normalized, frag_spec_frag_int_normalized

    def get_most_likely_type_of_ion(self, mass_of_ion, mi_spec, frag_spec):
        closest_mass_in_mi_spec = min(list(mi_spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(((mass_of_ion - x) / mass_of_ion) * 1000000))
        closest_mass_in_frag_spec = min(list(frag_spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(((mass_of_ion - x) / mass_of_ion) * 1000000))
        mi_spec_intensity = mi_spec.summarized_mass_intensity_dict[closest_mass_in_mi_spec]
        frag_spec_intensity = frag_spec.summarized_mass_intensity_dict[closest_mass_in_frag_spec]
        if ((abs(closest_mass_in_mi_spec-mass_of_ion)/mass_of_ion)*1000000) >= self.kwargs["mass_deviation"]:
            print("Mass found in MI Spec does not match with the given mass.")
            mi_spec_intensity = 0.0001
        if ((abs(closest_mass_in_frag_spec-mass_of_ion)/mass_of_ion)*1000000) >= self.kwargs["mass_deviation"]:
            print("Mass found in Frag Spec does not match with the given mass.")
            frag_spec_intensity = 0.0001
        return_ratio = mi_spec_intensity / frag_spec_intensity # if >1 (more intensity of the mass in the MI spec) --> mass is molecular ion; if <1: mass is a fragment.
        if return_ratio >= 1:
            type_of_ion = "molecular_ion"
        else:
            type_of_ion = "fragment_ion"
        return return_ratio, mi_spec_intensity, frag_spec_intensity, type_of_ion

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
    print("This is the UVenture module. It is not meant to be run directly.")
