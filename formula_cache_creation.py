import copy
import itertools
import os
from pyteomics import mass

def get_most_abundant_mass(formula, mass_resolution_ppm=15):
    if len(formula) == 0:
        return None
    # Get the most abundant isotope mass of the given formula
    import pyteomics.mass
    isotopologues_list = list(pyteomics.mass.mass.isotopologues(formula))

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
        except Exception as e:
            print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
            print(e)
            break
    sorted_masses_according_to_abundance = [m for _,m in sorted(zip(new_abundances_list, new_masses_list), reverse=True)]
    sorted_abundance_list = sorted(new_abundances_list, reverse=True)
    abundance_dict = dict(zip(sorted_masses_according_to_abundance, sorted_abundance_list))
    return list(abundance_dict.keys())[0]


def get_formula_to_dict(formula_string):
    if len(formula_string) == 0:
        return {}
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

# Function to convert an index to a combination
def index_to_combination(index, max_number, num_atoms):
    combination = []
    for _ in range(num_atoms):
        combination.append(index % (max_number + 1))
        index //= (max_number + 1)
    return combination[::-1]


def combination_to_index(combination, max_number):
    index = 0
    base = max_number + 1
    for num in combination:
        index = index * base + num
    return index

# Define the set of 5 different atoms
atoms = ['H', 'C', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', "Li", "Na", "K", 'Cr', 'B', 'Fe', 'Si', 'Al', 'Be', 'Mg', 'Ca', 'Ti', 'Mn', 'Sn', 'Zn', 'As', 'Se', 'Pb', 'Hg']
atoms.reverse()
atoms_starting_from_c = copy.deepcopy(atoms)
atoms_starting_from_c.reverse()


formula_cache_basefilepath = "Formula_Predictions//formulas_for_mass_"
if not os.path.exists("Formula_Predictions//"):
    os.makedirs("Formula_Predictions")
index_tracking_filepath = "Formula_Predictions//index_tracking.txt"


# Define the maximum number of each atom
max_number = 120

# Example last processed combination
start_calculation_at_index = 0
start_combination = index_to_combination(start_calculation_at_index, max_number, len(atoms))

#absolute maximum atom numbers. Ordered from Hg to H. So Hydrogen is the last element in the list.
absolute_maximum_combination = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 3, 2, 1, 1, 1, 3, 5, 10, 30, 3, 8, 12, 30, 20, 70]



while True:
    def resume_from2(combination, max_number, num_atoms):
        start_index = sum(c * ((max_number + 1) ** i) for i, c in enumerate(reversed(combination)))
        for index in range(start_index, (max_number + 1) ** num_atoms):
            combo = index_to_combination(index, max_number, num_atoms)
            yield combo

    combinations = resume_from2(start_combination, max_number, len(atoms))

    for combo in combinations:
        reversed_combo = list(combo)
        reversed_combo.reverse()
        reversed_combo = tuple(reversed_combo)

        formula = ''.join(f"{atom}{num}" for atom, num in zip(atoms_starting_from_c, reversed_combo) if num > 0)
        formula_dict = get_formula_to_dict(formula)

        if (formula_dict.get("H", 0) > formula_dict.get("C", 1) * 6) and (formula_dict.get("C", 0) > 2):
            curr_index = combination_to_index(combo, max_number)
            missing_to_wrap_around = max_number + 1 - combo[-1]
            start_combination = index_to_combination(curr_index+missing_to_wrap_around, max_number, len(atoms))
            break
        if (formula_dict.get("H", 1)*5 < formula_dict.get("C", 0)) and (formula_dict.get("C", 0) > 3):
            continue
        masse = get_most_abundant_mass(formula)
        if masse is None:
            continue

        with open(formula_cache_basefilepath + str(int(masse)) + ".txt", "a") as file:
            file.write(str(formula) + "," + str(round(masse, 5)) + "\n")
        if formula_dict.get("H", 0) == 1:
            with open(index_tracking_filepath, "w") as file:
                file.write("Curr index: " + str(combination_to_index(combo, max_number)) + "\n")
                file.write("Max number: " + str(max_number) + "\n")
                file.write("Num atoms: " + str(len(atoms)) + "\n")
                file.write("Curr combination: " + str(combo) + "\n")
                file.write("Atoms: " + str(atoms) + "\n")
        
        if formula_dict.get("C", 0) == 1 and formula_dict.get("H", 0) == 1:
            print(formula)
            print("Curr index: " + str(combination_to_index(combo, max_number)) + " of " + str((max_number + 1) ** len(atoms)))

        start_new_iteration = False
        for i in range(len(combo)):
            if combo[i] > absolute_maximum_combination[i]:
                curr_index = combination_to_index(combo, max_number)
                missing_to_wrap_around = max_number + 1 - combo[-1]
                new_combination = copy.deepcopy(combo)
                new_combination[i-1] = new_combination[i-1] + 1
                new_combination[i] = 0
                start_combination = copy.deepcopy(new_combination)
                start_new_iteration = True
                break
        if start_new_iteration:
            break

