import argparse
import os
import sys
import time
import pyteomics.mass
import copy
import multiprocessing
import sys
import time
import psutil
import inspect

class MyMultiprocessing:
    def __init__(self, only_multithread=False):
        self.only_multithread = only_multithread
        self.functions = []
        self.args = []
        self.pool = None
        self.pool_list = []
        self.pool_list_starmap = []
        self.ret_values = []
        self.max_cores = psutil.cpu_count() - 1
        self.timeout = 0.02
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        self.caller = mod.__name__
        self.caller_file_location = mod.__file__
        if self.caller == "__main__":
            pass
        else:
            print("You need to put the class into the following if-statement: \n"
                  "if __name__ == '__main__':\n"
                  "    mp = MultiprocessingClass.MyMultiprocessing()\n"
                  "    mp.add_function(fn, arglist)\n"
                  "    mp.run_pool()\n"
                  "    ret = mp.get_return_values()")
            time.sleep(1)
            sys.exit()
        if self.only_multithread == False:
            self.pool = multiprocessing.Pool(processes=self.max_cores)
        else:
            from multiprocessing.pool import ThreadPool
            self.pool = ThreadPool(processes=self.max_cores)
            print("onlyMultithread")

    def add_function(self, function, arg_list):
        if hasattr(function, "__self__"):
            print("Function is bound method!\n"
                  "The function you want to add is not static. If the function is inside a class, consider adding the "
                  "'@staticmethod' decorator.\n"
                  "Otherwise use a function that is not inside the class "
                  "(that is not relying on the 'self' argument). \n"
                  "Only functions that do not have the 'self' argument are valid functions to use with this "
                  "multiprocessing module.\n"
                  "The given function was therefore NOT added!")
            return None
        self.functions.append(function)

        if arg_list is None:
            print("arg_list must not be empty!\n"
                  "At least one argument is needed in the function for the class to work properly!!!")
        if arg_list == []:
            print("arg_list must not be empty!\n"
                  "At least one argument is needed in the function for the class to work properly!!!")
        if arg_list == ():
            print("arg_list must not be empty!\n"
                  "At least one argument is needed in the function for the class to work properly!!!")
        if isinstance(arg_list, list):
            pass
        elif isinstance(arg_list, int):
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        elif isinstance(arg_list, str):
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        elif isinstance(arg_list, float):
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        else:
            temp_arg_list = []
            temp_arg_list.append(arg_list)
            arg_list = temp_arg_list
        self.args.append(arg_list)


    def get_return_values(self):
        start_time = time.time()
        while len(self.ret_values) < len(self.functions):
            for pool in list(self.pool_list):
                try:
                    self.ret_values.append(pool.get(timeout=self.timeout))
                    self.pool_list.remove(pool)
                except:
                    pass
            for pool in list(self.pool_list_starmap):
                try:
                    self.ret_values.append(pool.get(timeout=self.timeout))
                    self.pool_list_starmap.remove(pool)
                except:
                    pass
            if start_time < (time.time() - self.timeout):
                break
        return self.ret_values

    def run_pool(self, wait_for_finish=True):
        for function in self.functions:
            self.pool_list.append(self.pool.apply_async(func=function, args=tuple(self.args[len(self.pool_list)])))
            self.pool_list_starmap.append(self.pool.starmap_async(function, tuple(self.args[len(self.pool_list_starmap)])))
        if wait_for_finish:
            self.pool.close()
            self.pool.join()
            start_time = time.time()
            while len(self.ret_values) < len(self.functions):
                for pool in list(self.pool_list):
                    try:
                        self.ret_values.append(pool.get(timeout=self.timeout))
                        self.pool_list.remove(pool)
                    except:
                        pass
                for pool in list(self.pool_list_starmap):
                    try:
                        self.ret_values.append(pool.get(timeout=self.timeout))
                        self.pool_list_starmap.remove(pool)
                    except:
                        pass
                if start_time < (time.time()-self.timeout):
                    break
        else:
            pass





def get_formula_to_dict(formula_string):
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


def get_formula_string_from_dict(formula_dict):
    return "".join([str(a) + str(n) for a, n in formula_dict.items()])


def simulate_isotope_pattern_of_formula(formula, mass_resolution_ppm=10):
    import pyteomics.mass
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

            # indices_list = [masses_list.index(m) for m in curr_mass_plus_deviation_list]
            #
            # for i in range(len(masses_list)-1, -1, -1):
            #     if i in indices_list:
            #         del masses_list[i]
            #         del abundances_list[i]
            #     print(masses_list)
            # print(curr_mass_plus_deviation_list)
        except Exception as e:
            print("EXCEPTION IN simulate_isotope_pattern_of_formula()!!!")
            print(e)
            break
    sorted_masses_according_to_abundance = [m for _,m in sorted(zip(new_abundances_list, new_masses_list), reverse=True)]
    sorted_abundance_list = sorted(new_abundances_list, reverse=True)
    abundance_dict = dict(zip(sorted_masses_according_to_abundance, sorted_abundance_list))
    return abundance_dict


def add_1_and_handle_carry(combined_string):
    # Split the combined string into individual two-digit components
    digits = [int(combined_string[i:i+2]) for i in range(0, len(combined_string), 2)]
    bases = [27, 71, 6, 16, 4, 4, 21, 6, 4, 3, 6, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    bases.reverse()
    #C, H, N, O, S, P, F, Cl, Br, I, Si, B, Li, Na, K, Mg, Ca, Mn, Fe, Cu, Zn, Se, Mo, As, Co, Ni, Cr

    # Add 1 to the last digit (base 3)
    digits[-1] += 1

    # Handle carry-over using a loop
    for i in range(len(digits) - 1, -1, -1):
        if digits[i] >= bases[i]:
            digits[i] = 0
            if i > 0:
                digits[i - 1] += 1
            else:
                raise ValueError("Number1 exceeded two-digit limit")

    # Combine the numbers back into a string with two digits each
    return ''.join(f"{d:02}" for d in digits)


def generate_nth_combined_string(n):
    digits = []
    bases = [27, 71, 6, 16, 4, 4, 21, 6, 4, 3, 6, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    bases.reverse()
    max_number = 1
    for base in range(len(bases)):
        if base == 0:
            max_number = bases[base]
        else:
            max_number = max_number * bases[base]
    # Iterate through each base in reverse order to calculate the digit values
    for base in reversed(bases):
        digit = n % base
        digits.append(digit)
        n //= base
    
    # Reverse the digits list to match the original order
    digits.reverse()
    
    # Combine the numbers back into a string with two digits each
    combined_string = ''.join(f"{digit:02}" for digit in digits)
    return combined_string


def combined_string_to_formula(combined_string):
    position_element_dict = {26: "C", 25: "H", 24: "N", 23: "O", 22: "S", 21: "P", 20: "F", 19: "Cl", 18: "Br", 17: "I", 16: "Si", 15: "B", 14: "Li", 13: "Na", 12: "K", 11: "Mg", 10: "Ca", 9: "Mn", 8: "Fe", 7: "Cu", 6: "Zn", 5: "Cr", 4: "Mo", 3: "As", 2: "Co", 1: "Ni", 0: "Se"}
    # Split the combined string into individual two-digit components
    atom_counts = [int(combined_string[i:i+2]) for i in range(0, len(combined_string), 2)]
    formula_dict = {}
    for i in range(len(atom_counts)):
        if atom_counts[i] > 0:
            formula_dict[position_element_dict[i]] = atom_counts[i]
    formula_dict = {k: v for k, v in reversed(formula_dict.items())}
    formula_string = "".join([str(a) + str(n) for a, n in formula_dict.items()])
    return formula_dict


def formula_creation(startstring="000000000000000000000000000000000000000000000000000000", iterations=1000):
    combined_string = startstring
    formulas_to_write = []
    for i in range(iterations):
        write_str = ""
        combined_string = add_1_and_handle_carry(combined_string)
        f_dict = combined_string_to_formula(combined_string)
        if f_dict.get("H", 0) / f_dict.get("C", 999) > 8:
            continue
        if f_dict.get("O", 0) / f_dict.get("C", 999) > 4:
            continue
        if f_dict.get("N", 0) / f_dict.get("C", 999) > 4:
            continue
        if f_dict.get("S", 0) / f_dict.get("C", 999) > 4:
            continue
        if f_dict.get("P", 0) / f_dict.get("C", 999) > 2:
            continue
        if f_dict.get("F", 0) / f_dict.get("C", 999) > 8:
            continue
        dbe = f_dict.get("C", 0) - (f_dict.get("H", 0) / 2) - (f_dict.get("F", 0)/2) - (f_dict.get("Cl", 0)/2) - (f_dict.get("Br", 0)/2) - (f_dict.get("I", 0)/2) + (f_dict.get("N", 0) / 2) + 1
        if dbe < -5 and f_dict.get("C", 0) >= 1:
            continue

        write_str = get_formula_string_from_dict(f_dict)
        write_str = write_str + ","
        if (f_dict.get("Cl", 0) >= 2) or \
                (f_dict.get("Br", 0) >= 1) or \
                    (f_dict.get("C", 0) >= 20) or \
                        (f_dict.get("S", 0) >= 4) or \
                            (f_dict.get("B", 0) >= 1) or \
                                (f_dict.get("Si", 0) >= 3) or \
                                    (f_dict.get("Mg", 0) >= 1) or \
                                        (f_dict.get("Cu", 0) >= 1) or \
                                            (f_dict.get("Zn", 0) >= 1) or \
                                                (f_dict.get("Se", 0) >= 1) or \
                                                    (f_dict.get("Mo", 0) >= 1) or \
                                                        (f_dict.get("Cr", 0) >= 1):
            write_str = write_str + str(round((list(simulate_isotope_pattern_of_formula(f_dict).keys())[0]), 6))
        else:
            write_str = write_str + str(round(calculate_mass_most_abundant(f_dict), 6))
        formulas_to_write.append(write_str)
    return formulas_to_write
    

def formula_creations_save(parentfolder_name="formula_predictions_process_x//", startstring="", iterations=1000):
    os.makedirs(parentfolder_name, exist_ok=True)
    formulas_to_write = formula_creation(startstring=startstring, iterations=iterations)
    for line in formulas_to_write:
        with open(parentfolder_name + "formulas_for_mass_" + str(int(float(line.strip().split(",")[1]))) + ".txt", "a") as f:
            f.write(line + "\n")
    return None


def do_work_singleproc(startnumber, end_number, iterations_per_save=1000, parentfolder_name="formula_predictions_process_x//"):
    for i in range(startnumber, end_number+1, iterations_per_save):
        formula_creations_save(parentfolder_name=parentfolder_name,
                                startstring=generate_nth_combined_string(i),
                                  iterations=iterations_per_save)
        print("Progress in percent: " + str(round(((i-startnumber)/(end_number-startnumber))*100, 2)))
    return None


def calculate_mass_most_abundant(formula_dict):
    if isinstance(formula_dict, str):
        formula_dict = get_formula_to_dict(formula_dict)
        print("Formula was given as a string. Try to give it as a dictionary next time!")
    #starttime = time.time()
    atom_weight_dict = {"H": 1.007825, 
                        "C": 12.00000, 
                        "N": 14.003074, 
                        "O": 15.994915, 
                        "S": 31.972072, 
                        "P": 30.973763, 
                        "F": 18.998403, 
                        "Cl": 34.968853, 
                        "Br": 78.918336, 
                        "I": 126.904477, 
                        "Si": 27.976928, 
                        "B": 11.009305, 
                        "Li": 7.016005, 
                        "Na": 22.989770, 
                        "K": 38.963708, 
                        "Mg": 23.985045, 
                        "Ca": 39.962591, 
                        "Mn": 54.938046, 
                        "Fe": 55.934939, 
                        "Cu": 62.929599, 
                        "Zn": 63.929145, 
                        "Se": 79.916521, 
                        "Mo": 97.905405, 
                        "As": 74.921596, 
                        "Co": 58.933198, 
                        "Ni": 57.935347, 
                        "Cr": 51.940510}
    mass = 0
    for atom, count in formula_dict.items():
        mass += atom_weight_dict[atom] * count
    #print("Time elapsed: " + str(time.time()-starttime))
    #starttime = time.time()
    #mass2 = round((pyteomics.mass.calculate_mass(formula_dict)), 6)
    #print("Time elapsed: " + str(time.time()-starttime))
    return mass


def do_work_multiprocessing(startnumber, end_number, parentfolder_name="formula_predictions_multiproc//"):
    mp = MyMultiprocessing()
    iterations_per_multiproc = int(end_number / mp.max_cores-1)+1
    for i in range(startnumber, end_number+1, iterations_per_multiproc):
        mp.add_function(do_work_singleproc, [i, i+iterations_per_multiproc, 1000, parentfolder_name + "process_" + str(i) + "//"])
    mp.run_pool()
    return None


def set_as_main_func_to_be_executable_via_terminal_argparse():
    parser = argparse.ArgumentParser(description="Formula creation for mass prediction")
    parser.add_argument("-s", "--startnumber", type=int, help="Start number for the formula creation")
    parser.add_argument("-e", "--endnumber", type=int, help="End number for the formula creation")
    parser.add_argument("-i", "--iterations", type=int, help="Iterations per save")
    parser.add_argument("-p", "--parentfolder", type=str, help="Parent folder name")
    args = parser.parse_args()
    if args.startnumber is None:
        print("Start number is missing!")
        sys.exit()
    if args.endnumber is None:
        print("End number is missing!")
        sys.exit()
    if args.iterations is None:
        args.iterations = 1000
    if args.parentfolder is None:
        args.parentfolder = "formula_predictions_task_x//"
    do_work_singleproc(args.startnumber, args.endnumber, args.iterations, args.parentfolder)


def _use_this_func_to_start_multiprocessing_via_terminal():
    import subprocess
    name_of_the_python_script = "FormulaCreation.py"
    # List of args to run in parallel
    tasks = [(0, 10000, 1000, "task1//"), (10000, 20000, 1000, "task2//"), (20000, 30000, 1000, "task3//")]
    # List to hold subprocess references
    processes = []
    for start, end, iterations, foldername in tasks:
        # Call the script with the specific task ID and duration
        process = subprocess.Popen(["python", name_of_the_python_script, str(start), str(end), str(iterations), str(foldername)])
        processes.append(process)

    # Wait for all processes to complete
    for process in processes:
        process.wait()


def combined_string_to_number(combined_string, bases):
# Ensure that the combined string length matches the total digits required
    if len(combined_string) != sum([len(f"{base-1:02}") for base in bases]):
        raise ValueError("Combined string length does not match expected length based on bases.")
    
    # Parse the combined string into individual components
    digits = []
    start_index = 0
    for base in bases:
        digit_length = len(f"{base-1:02}")  # length of each digit based on its base
        digits.append(int(combined_string[start_index:start_index + digit_length]))
        start_index += digit_length
    
    # Convert the individual components into a single number
    number = 0
    multiplier = 1
    for digit, base in zip(reversed(digits), reversed(bases)):
        number += digit * multiplier
        multiplier *= base
    
    return number


if __name__ == '__main__':
    #starttime = time.time()
    #do_work_singleproc(0, 16000, 1000, parentfolder_name="test1//") # for 16000 iterations 378 seconds
    #print("Time elapsed: " + str(time.time()-starttime))
    #print(generate_nth_combined_string(16000))
    #        C, H, N, O, S, P, F  Cl Br I  Si B  Li Na K  Mg Ca Mn Fe Cu Zn Se Mo As Co Ni Cr
    mystr = "00 00 00 00 00 00 00 00 00 00 00 00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00"
    mystr = mystr.replace(" ", "")
    mystr = mystr[::-1]
    print(mystr)
    bases = [27, 71, 6, 16, 4, 4, 21, 6, 4, 3, 6, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    bases.reverse()
    print(combined_string_to_number("000000000000000000000000000000000000000000000000000003", bases=bases))
    print(combined_string_to_number("000000000000000000000000000000000000000000000000000100", bases=bases))
    print(combined_string_to_formula(mystr))
    mystr = add_1_and_handle_carry(mystr)
    print(combined_string_to_number(mystr, bases=bases))
    print(combined_string_to_formula(mystr))
    mystr = add_1_and_handle_carry(mystr)
    print(combined_string_to_number(mystr, bases=bases))
    print(combined_string_to_formula(mystr))
    mystr = add_1_and_handle_carry(mystr)
    print(combined_string_to_number(mystr, bases=bases))
    print(combined_string_to_formula(mystr))
    mystr = add_1_and_handle_carry(mystr)
    print(combined_string_to_number(mystr, bases=bases))
    print(combined_string_to_formula(mystr))
    #starttime = time.time()
    #do_work_multiprocessing(0, 16000, parentfolder_name="test3//") # for 16000 iterations 142 seconds (with only multithread on default)
    #print("Time elapsed: " + str(time.time()-starttime))
    #starttime = time.time()
    #do_work_multiprocessing(0, 6000, parentfolder_name="test4//") # for 6000 iterations ~300 seconds (with only multithread on default)







