# formula_calculations.py
from UVenture.chemical_formula_parsing import get_formula_to_dict, get_formula_string_from_dict
from UVenture.MS_functions import summarize_mass_intensity_dict_for_isotopo_simulation
import traceback
import pyteomics.mass

def calc_dbe(formula):
    try:
        if isinstance(formula, str):
            formula = get_formula_to_dict(formula)
        if isinstance(formula, dict):
            pass
        dbe = ((formula.get("Si", 0) * 2) + (formula.get("C", 0)*2) + 2 + formula.get("N", 0) + formula.get("P", 0) - formula.get("H", 0) - formula.get("F", 0) - formula.get("Cl", 0) - formula.get("Br", 0) - formula.get("I", 0) ) / 2
        return dbe
    except Exception as e:
        print("Error calculating DBE in formula_calculations: " + str(e))
        print(traceback.format_exc())
        return 0
    

def get_mass_of_most_abundant_isotopologue_formula(formula):
    if len(formula) == 0:
        return None
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


def simulate_isotope_pattern_of_formula(formula, mass_resolution_ppm=10, debug_output=False):
    if len(formula) == 0:
        return None
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


def calculate_likelyhood_of_formula_dict(f_dict, charge_of_measured_mass=None, score_subst_rdbe_non_integer=10, score_subst_senior_rule=10, debug_output=False):
    formula_likelyness = 0
    try:
        # X/C ratio check from #https://pmc.ncbi.nlm.nih.gov/articles/PMC1851972/ The numbers are representing 99.7% of all available molecular formulas.
        # Not true anymore. It has been made less strict
        if "C" in list(f_dict.keys()):
            if (f_dict.get("H", 0) / f_dict.get("C", 0) < 0.2) or (f_dict.get("H", 0) / f_dict.get("C", 0) > 4):
                formula_likelyness = formula_likelyness -40
            if (f_dict.get("F", 0) / f_dict.get("C", 0) < 0) or (f_dict.get("F", 0) / f_dict.get("C", 0) > 4):
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
        if heteroatom_count > 5 and f_dict.get("C", 0) <= int(heteroatom_count / 4): # ... and more than 4x more heteroatoms than carbon atoms
            formula_likelyness = formula_likelyness - float((heteroatom_count / 4) - f_dict.get("C", 0)) * 50

        if ("C" in list(f_dict.keys())) and ("N" in list(f_dict.keys())):
            if charge_of_measured_mass < 0:
                if int(f_dict["N"]) > 2 and (int(f_dict["C"]) / int(f_dict["N"]) <= 2):
                    formula_likelyness = formula_likelyness - ((float(f_dict["N"]) - 2) ** 2) * 50
            else:
                if int(f_dict["N"]) > 2 and (int(f_dict["C"]) / int(f_dict["N"]) <= 0.74):
                    formula_likelyness = formula_likelyness - ((float(f_dict["N"]) - 2) ** 2) * 50

        
        if ("C" in list(f_dict.keys())) and ("H" in list(f_dict.keys())):
            if (float(f_dict["H"]) / float(f_dict["C"]) >= 2):
                formula_likelyness = formula_likelyness - (((float(f_dict["H"]) / float(f_dict["C"])) - 2) ** 3) * 50

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
