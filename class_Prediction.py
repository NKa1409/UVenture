import ast
import copy
import datetime
import math
import os
import traceback


import UVenture.formula_cache as formula_cache
import UVenture.chemical_formula_parsing as chemical_formula_parsing
import UVenture.MS_functions as MS_functions
import UVenture.formula_calculations as formula_calculations
import UVenture.plotting as plotting



class Prediction:
    def __init__(self, ms_file, mass, spec, spec_before=None, spec_after=None, **kwargs):
        self.debug_output = True

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
        print("Starting prediction class")
        
        if "prediction_subfolder" in kwargs:
            pred_folder = os.path.join(self.ms_file.parentfolder, "predictions", kwargs["prediction_subfolder"])
        else:
            pred_folder = os.path.join(self.ms_file.parentfolder, "predictions", str(self.spec.index) + "_" + str(round(self.mass, 4)))
        if "absolute_pred_folder" in kwargs:
            pred_folder = kwargs["absolute_pred_folder"]
        if "prediction_subfolder" in kwargs and "absolute_pred_folder" in kwargs:
            print("ERROR: You can only use either 'prediction_subfolder=' or 'absolute_pred_folder=', not both!")
            print("Using absolute_pred_folder now...")
            print(kwargs["absolute_pred_folder"])
        
        default_kwargs = {"pred_folder":pred_folder,  
                          "pred_log_filepath": os.path.join(pred_folder, "prediction_log_for_mass_" + str(round(self.mass, 4)) + ".txt"),
                          "pred_save_matplotlib_plot_of_isotopologues":False,
                          "pred_save_go_plot_of_isotopologues":False,
                          "pred_save_xic_plot":False,
                          "pred_save_detailed_log":True,
                          "pred_formula_cache_folder_path":"U://MyFolder//MONOTONS//Filtermessungen//Formula_Predictions//",

                          "pred_atoms_to_keep_in_prediction": ['C', 'H', 'N', 'O', 'S', 'Cl', 'Br'],

                          "mass_deviation":11,
                          "charge_of_measured_mass":-1,

                          "pred_minimum_assumed_noise":10000,
                          "pred_noise_divisor_for_isotopologue_calculation":5,
                          "pred_add_value_to_score_if_isotopo_was_found":250,
                          "pred_isotopologue_score_e_function_exponent":0.4,
                          "pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score":30,
                          "pred_ratio_influence_ppm_deviation_vs_intensity": 0.5,
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


        self.intensity_of_ion = self.spec.summarized_mass_intensity_dict[self.mass]
        #self.intensity_of_ion = sum([self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if abs(((self.spec.summarized_masses[i] - self.mass)/self.mass)*1000000) <= self.kwargs["mass_deviation"]])
        self.make_op_log_entry("")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("=============================================================")   
        self.make_op_log_entry("INFO:\t" + "Creating prediction object for mass: " + str(self.mass) + " at retention time: " + str(self.rt) + " at index: " + str(self.spec.index))    
        self.make_op_log_entry("INFO:\t" + "Intensity of the molecular ion: " + str(self.intensity_of_ion))
        self.make_op_log_entry("INFO:\t" + "Assuming a mass deviation of: " + str(self.kwargs["mass_deviation"]))
        print("Intensity of ion for prediction: " + str(self.intensity_of_ion))
        self.xic = self.ms_file.get_xic(self.mass, ((self.kwargs["mass_deviation"]*self.mass)/1000000), requested_filter_mode=self.spec.filter_mode)
        self.make_op_log_entry("INFO:\t" + "XIC calculated.")
        self.identified_peaks_for_mass = MS_functions.get_peaks_in_xy_series(self.xic[0], self.xic[1], sg_window=10, sg_order=3)
        self.make_op_log_entry("INFO:\t" + "Number of peaks found in XIC: " + str(len(self.identified_peaks_for_mass)))
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
            print("Peak is detected where prediction should be made.")
            self.peak_found = True
            self.make_op_log_entry("INFO:\t" + "Peak is detected where a prediction should be made!")


        neutral_mass = self.mass + (0.000548 * self.kwargs["charge_of_measured_mass"])
        _, possible_formulas = formula_cache.get_formula_from_cache(self.kwargs["pred_formula_cache_folder_path"], neutral_mass, self.kwargs["mass_deviation"], debug_output=self.debug_output)
        self.make_op_log_entry("INFO:\t" + "Retrieved formulas from cache.")
        if self.kwargs["charge_of_measured_mass"] < 0:
            self.make_op_log_entry("INFO:\t" + "Charge of the measured mass is negative. Removing formulas with Li, Na and K...")
            possible_formulas = {key: value for key, value in possible_formulas.items() if "Na" not in key and "K" not in key and "Li" not in key and "Mg" not in key and "Ca" not in key and "Be" not in key}

        atoms_to_keep_in_prediction = self.kwargs["pred_atoms_to_keep_in_prediction"]
        possible_formulas_dict_list = [chemical_formula_parsing.get_formula_to_dict(formula_string=f) for f, dev in possible_formulas.items()]
        possible_deviations_list = [dev for f, dev in possible_formulas.items()]
        new_possible_formulas = {}
        for i in range(len(possible_deviations_list)):
            if set(list(possible_formulas_dict_list[i].keys())) - set(atoms_to_keep_in_prediction):
                continue
            else:
                new_possible_formulas[chemical_formula_parsing.get_formula_string_from_dict(possible_formulas_dict_list[i])] = possible_deviations_list[i]
        possible_formulas = new_possible_formulas

        print(possible_formulas)

        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("OUTPUT OF SIMPLE FORMULA PREDICTION (ONLY MASS DEVIAITON):")
        for formula, deviation in possible_formulas.items():
            self.make_op_log_entry("{:>20} \t {:>20}".format(str(formula), str(round(deviation, 2))))
        
        molecule_formulas = [[chemical_formula_parsing.get_formula_to_dict(formula), deviation] for formula, deviation in possible_formulas.items()]
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
        print("Formula score dict: " + str(self.formula_score_dict))
        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICTIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>20} \t {:>10}".format(str(chemical_formula_parsing.get_formula_string_from_dict(item[0])), str(round(item[1], 2))))
        
        if self.kwargs["pred_include_likelyhood_of_formula"] == True:
            self.formula_score_dict = self.add_likelyhood_of_formula_to_score(self.formula_score_dict)


        self.formula_score_dict = {f: s for f, s in self.formula_score_dict.items() if s > self.kwargs["pred_reject_formula_if_score_lower_than"]}
        self.formula_score_dict = dict(sorted(self.formula_score_dict.items(), key=lambda x: x[1], reverse=True))

        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICTIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>20} \t {:>10}".format(str(chemical_formula_parsing.get_formula_string_from_dict(item[0])), str(round(item[1], 2))))

        if self.kwargs["pred_include_peak_matching_of_isotopologues"]:
            self.formula_score_dict = self.get_peak_matching_of_isotopologues(self.formula_score_dict)
            self.formula_score_dict = {f: s for f, s in self.formula_score_dict.items() if s > self.kwargs["pred_reject_formula_if_score_lower_than"]}
            self.formula_score_dict = dict(sorted(self.formula_score_dict.items(), key=lambda x: x[1], reverse=True))

        self.make_op_log_entry("=============================================================")
        self.make_op_log_entry("FORMULA SCORE DICTIONARY:")
        for item in list(self.formula_score_dict.items()):
            self.make_op_log_entry("{:>20} \t {:>10}".format(str(chemical_formula_parsing.get_formula_string_from_dict(item[0])), str(round(item[1], 2))))

        # If the user wants to save the XIC plot, create it and save it
        if self.kwargs["pred_save_xic_plot"] == True:
            self.make_op_log_entry("INFO:\t" + "Start saving XIC plot...")
            self.xic_plot_filepath = os.path.join(self.kwargs["pred_folder"], "xic_" + str(round(self.mass, 4)) + "+-" + str(self.kwargs["mass_deviation"]) + "_" + str(self.spec.filter_mode) + ".png")
            title = "XIC for mass: " + str(round(self.mass, 4)) + " at RT: " + str(round(self.spec.rt, 2)) + " seconds\nFilter mode: " + str(self.spec.filter_mode)
            plotting.create_xic(self.xic[0], self.xic[1], title, self.xic_plot_filepath, retention_time=self.spec.rt)
            self.make_op_log_entry("INFO:\t" + "Finished saving XIC plot...")
        
        # If the user wants to save the matplotlib plot of isotopologues, create it and save it
        if self.kwargs["pred_save_matplotlib_plot_of_isotopologues"] == True and len(list(self.formula_score_dict.keys())) > 0:
            self.make_op_log_entry("INFO:\t" + "Start saving matplotlib plot of isotopologues...")
            self.matplotlib_plot_filepath = os.path.join(self.kwargs["pred_folder"], "isotopo_matplotlib_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.formula_score_dict.keys())[0]).items()])) + ".png")
            plotting.create_isotopo_plot(self.spec.summarized_masses, self.spec.summarized_intensities, ast.literal_eval(list(self.formula_score_dict.keys())[0]), self.matplotlib_plot_filepath)
            self.make_op_log_entry("INFO:\t" + "Finished saving matplotlib plot of isotopologues...")
        
        # If the user wants to save the go plot of isotopologues, create it and save it
        if self.kwargs["pred_save_go_plot_of_isotopologues"] == True and len(list(self.formula_score_dict.keys())) > 0:
            self.make_op_log_entry("INFO:\t" + "Start saving go plot of isotopologues...")
            self.go_plot_filepath = os.path.join(self.kwargs["pred_folder"], "isotopo_go_plot_" + str(round(self.mass, 4)) + "_" + str("".join([str(a) + str(n) for a, n in ast.literal_eval(list(self.formula_score_dict.keys())[0]).items()])) + ".html")
            plotting.create_isotopo_plot_with_go(self.spec.summarized_masses, self.spec.summarized_intensities, ast.literal_eval(list(self.formula_score_dict.keys())[0]), self.go_plot_filepath, mass_deviation=self.kwargs["mass_deviation"])
            self.make_op_log_entry("INFO:\t" + "Finished saving go plot of isotopologues...")
        
        try:
            if len(list(self.formula_score_dict.keys())) >= 1:
                self.make_op_log_entry("INFO:\t" + "Setting best formula prediction...")
                self.best_formula_prediction = ast.literal_eval(list(self.formula_score_dict.keys())[0])
                self.score_of_best_formula = self.formula_score_dict[(list(self.formula_score_dict.keys())[0])]
                self.simulated_isotopologue_pattern_for_best_formula = formula_calculations.simulate_isotope_pattern_of_formula(self.best_formula_prediction, mass_resolution_ppm=self.kwargs["mass_deviation"], debug_output=self.debug_output)
                #self.simulated_isotopologue_pattern_for_best_formula = {mass: abundance, ...}
                self.make_op_log_entry("INFO:\t" + "Best formula prediction set. " + str(chemical_formula_parsing.get_formula_string_from_dict(self.best_formula_prediction)) + "   " + str(self.score_of_best_formula))
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
                formula_likelyness = formula_calculations.calculate_likelyhood_of_formula_dict(f_dict,
                                                                            charge_of_measured_mass=self.kwargs["charge_of_measured_mass"],
                                                                            score_subst_rdbe_non_integer=self.kwargs["pred_formula_likelihood_substract_score_rdbe_non_integer"],
                                                                            score_subst_senior_rule=self.kwargs["pred_formula_likelihood_substract_score_rdbe_higher_than_senior_rule"],
                                                                            debug_output=self.debug_output)
                formula_score_dict[formula_dict_str] = float(formula_score_dict[formula_dict_str]) + formula_likelyness
                self.make_op_log_entry("INFO:\t" + "General likelyness of formula: " + str(chemical_formula_parsing.get_formula_string_from_dict(f_dict)) + " calculated to be: " + str(round(formula_likelyness, 2)))
            except Exception as e:
                print("ERROR in formula score likelyhood: " + str(e))    
        self.make_op_log_entry("Finished adding likelyhood of formula to the score...")
        return formula_score_dict

    def check_for_isotope_pattern(self, formula_dict):
        influence_ppm_deviation = self.kwargs["pred_ratio_influence_ppm_deviation_vs_intensity"]
        influence_intensity_changes = 1-self.kwargs["pred_ratio_influence_ppm_deviation_vs_intensity"]

        isotope_pattern_dict = formula_calculations.simulate_isotope_pattern_of_formula(formula_dict, debug_output=self.debug_output)
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
            closest_isotopo_mass = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs(isotopo_mass - x))
            closest_isotopo_mass_before = min(list(self.spec_before.summarized_mass_intensity_dict.keys()), key=lambda x: abs(isotopo_mass - x))
            closest_isotopo_mass_after = min(list(self.spec_after.summarized_mass_intensity_dict.keys()), key=lambda x: abs(isotopo_mass - x))

            self.intensity_of_ion = 0
            if ( ( abs(closest_isotopo_mass - isotopo_mass) / isotopo_mass ) * 1000000 ) <= self.kwargs["mass_deviation"]:
                self.intensity_of_ion = self.spec.summarized_mass_intensity_dict[closest_isotopo_mass]
            if ((abs(closest_isotopo_mass_before - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                self.intensity_of_ion += self.spec_before.summarized_mass_intensity_dict[closest_isotopo_mass_before]
            if ((abs(closest_isotopo_mass_after - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                self.intensity_of_ion += self.spec_after.summarized_mass_intensity_dict[closest_isotopo_mass_after]
        else:
            self.intensity_of_ion = 0
            if ((abs(self.mass - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                self.intensity_of_ion = self.spec.summarized_mass_intensity_dict[self.mass]
            #self.intensity_of_ion = sum([self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if abs(((self.spec.summarized_masses[i] - self.mass) / self.mass) * 1000000) <= self.kwargs["mass_deviation"]])
        print("Intensity of the ion to check: " + str(self.intensity_of_ion))

        for isotopo in list(isotope_pattern_dict.items()):
            break_the_isotopo_prediction = False
            isotopo_mass = isotopo[0]
            measured_mass_with_minimal_deviation = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
            deviation = ((isotopo_mass - measured_mass_with_minimal_deviation) / isotopo_mass) * 1000000
            theoretical_intensity = (self.intensity_of_ion / (list(isotope_pattern_dict.items())[0][1])) * isotopo[1]
            if theoretical_intensity <= 0.01:
                theoretical_intensity = 0.01
            noise = sorted(list(self.spec.summarized_mass_intensity_dict.values()))[int(len(list(self.spec.summarized_mass_intensity_dict.items())) / self.kwargs["pred_noise_divisor_for_isotopologue_calculation"])] + self.kwargs["pred_minimum_assumed_noise"]

            if (not self.spec_before is None) and (not self.spec_after is None):
                print("Spec before and spec after is used.")
                closest_isotopo_mass = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
                closest_isotopo_mass_before = min(list(self.spec_before.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
                closest_isotopo_mass_after = min(list(self.spec_after.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
                measured_intensity = 0
                if ((abs(closest_isotopo_mass - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                    measured_intensity = self.spec.summarized_mass_intensity_dict[closest_isotopo_mass]
                if ((abs(closest_isotopo_mass_before - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                    measured_intensity += self.spec_before.summarized_mass_intensity_dict[closest_isotopo_mass_before]
                if ((abs(closest_isotopo_mass_after - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                    measured_intensity += self.spec_after.summarized_mass_intensity_dict[closest_isotopo_mass_after]
            else:
                print("Only one spec is used.")
                measured_intensity = 0
                closest_isotopo_mass = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
                if ((abs(closest_isotopo_mass - isotopo_mass) / isotopo_mass) * 1000000) <= self.kwargs["mass_deviation"]:
                    closest_isotopo_mass = min(list(self.spec.summarized_mass_intensity_dict.keys()), key=lambda x: abs((((isotopo_mass - x) / isotopo_mass) * 1000000) - initial_deviation))
                    measured_intensity = self.spec.summarized_mass_intensity_dict[closest_isotopo_mass]

            print("Measured intensity for ion: " + str(measured_intensity))

            if (deviation < self.kwargs["mass_deviation"]) and (abs(deviation - initial_deviation) <= (self.kwargs["mass_deviation"] / 2)) and (abs(measured_mass_with_minimal_deviation-isotopo_mass) <= (self.kwargs["mass_deviation"]/1000000)*150 ):
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
                    score = score ** (0.9/(2-self.kwargs["pred_isotopologue_score_e_function_exponent"]))

                score = score.real
                score = score - 0.5
                score = ((score * 100) + 0.001)
                print("Score1: " + str(score))

                score = score * (influence_intensity_changes / 0.5)

                score_multiplier_by_intensity_and_noise = theoretical_intensity / (theoretical_intensity + noise)
                if isotopo == list(isotope_pattern_dict.items())[0]:
                    print("First isotopo, continuing with next isotopo")
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
                    ppm_influence = ppm_influence * (influence_ppm_deviation / 0.5)
                    score = (score) - ppm_influence
                    print("Score4: " + str(score))
                    isotopes_found[isotopo[0]].append(score)
                    isotopes_found[isotopo[0]].append(measured_mass_with_minimal_deviation)
                    if (score - abs(self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]) >= previous_score):
                        isotopes_found[isotopo[0]][-2] = isotopes_found[isotopo[0]][-2] * 0.1
                    break
                else:
                    ppm_influence = ((2 ** ((self.kwargs["pred_isotopologue_score_e_function_exponent"] * 2) * (abs(deviation) - (self.kwargs["mass_deviation"] / 3)))) * 2)
                    ppm_influence = abs(ppm_influence * self.kwargs["pred_multiplier_isotopologue_influence_of_ppm_deviation_on_score"] * score_multiplier_by_intensity_and_noise)
                    if ppm_influence > 500:
                        ppm_influence = 500
                    ppm_influence = ppm_influence * (influence_ppm_deviation/0.5)
                    score = (score) - ppm_influence
                    print("Score4: " + str(score))
                    isotopes_found[isotopo[0]].append(score)
                    isotopes_found[isotopo[0]].append(measured_mass_with_minimal_deviation)
                
                if (score < self.kwargs["pred_stop_isotopologue_search_if_score_lower_than"]):
                    break

                if (score - ( abs(previous_score) / 10 ) >= previous_score):
                    isotopes_found[isotopo[0]][-2] = isotopes_found[isotopo[0]][-2] * 0.1
                    break
                previous_score = score
                if break_the_isotopo_prediction:
                    break
            else:
                isotopes_found[isotopo[0]] = [False, theoretical_intensity, measured_intensity]
                try:
                    multiplier_neg_score = abs(previous_score)
                    if multiplier_neg_score >= 50:
                        multiplier_neg_score = 50
                except UnboundLocalError:
                    multiplier_neg_score = 50
                negative_score = multiplier_neg_score * (theoretical_intensity / (theoretical_intensity + noise))
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
        self.intensity_of_ion = self.spec.summarized_mass_intensity_dict[self.mass]
        #self.intensity_of_ion = sum([self.spec.summarized_intensities[i] for i in range(len(self.spec.summarized_masses)) if abs(((self.spec.summarized_masses[i] - self.mass) / self.mass) * 1000000) <= self.kwargs["mass_deviation"]])
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
            simulated_isotopo_abundances_dict = formula_calculations.simulate_isotope_pattern_of_formula(formula_dict, debug_output=self.debug_output)

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
                isotopo_xic = self.ms_file.get_xic(isotopo_masse, mass_deviation_for_xic, requested_filter_mode=self.spec.filter_mode)
                area, peak1_rt, peakintensity1, peak2_rt, peakintensity2 = MS_functions.compare_peak_shape_similarity(mi_xic, isotopo_xic, self.rt, debug_output=self.debug_output)
                break_formula_evaluation = False
                if area < 0:
                    print("Area was not calculated correct. Setting it to: self.kwargs['oa_fragments_do_peak_computation_if_area_higher_than'] * 30")
                    area = self.kwargs["oa_fragments_do_peak_computation_if_area_higher_than"] * 30

                print("Area for isotopologue with mass: " + str(isotopo_masse) + "   ; area = " + str(area))
                theoretical_intensity_of_isotopo_peak = (intensity_of_mi / (list(simulated_isotopo_abundances_dict.items())[0][1])) * isotopo_abundance
                within_deviation_mass_list = [masse for masse in list(self.spec.summarized_mass_intensity_dict.keys()) if abs((((isotopo_masse - masse) / isotopo_masse) * 1000000) - initial_deviation) < (self.kwargs["mass_deviation"] / 2)]
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
                    delta_score = delta_score * 2   #previously * 1.3
                    break_formula_evaluation = True
                if break_formula_evaluation:
                    rel_isotopo_abundance = (sum(list(simulated_isotopo_abundances_dict.values())[list(simulated_isotopo_abundances_dict.values()).index(isotopo_abundance):]) / (list(simulated_isotopo_abundances_dict.items())[0][1] ))
                isotopo_xic_matching_multiplier = 1
                change_score = change_score + (delta_score * (score*rel_isotopo_abundance) * isotopo_xic_matching_multiplier)
                if rel_isotopo_abundance <= 0.001 or break_formula_evaluation:
                    break
            print("Change score for formula: " + str(chemical_formula_parsing.get_formula_string_from_dict(formula)) + "    ; Change score: " + str(change_score) + "    ; New score = " + str(score + change_score))
            new_score = score + change_score
            new_formula_score_dict[formula] = new_score
            self.make_op_log_entry("Change score for formula: " + str(chemical_formula_parsing.get_formula_string_from_dict(formula)) + "    ; Change score: " + str(change_score) + "    ; New score = " + str(score + change_score))
        return new_formula_score_dict