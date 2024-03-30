import ast
import pandas as pd
import matplotlib.pyplot as plt
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




def prepare_subsets(x, y, area_list, label_list):
    subset_color_dict = {"CH": "blue",
                         "CHO": "red",
                         "CHN": "lime",
                         "CHNO": "orange",
                         "CHS": "black",
                         "CHOS": "peru",
                         "CHNS": "black",
                         "CHNOS": "deeppink",
                         "else": "black"}
    subsets = []
    for label in list(set(label_list)):
        x_subset = []
        y_subset = []
        size_subset = []
        try:
            color = subset_color_dict[label]
        except:
            color = "black"
        for entry in range(len(label_list)):
            if label_list[entry] == label:
                x_subset.append(x[entry])
                y_subset.append(y[entry])
                size_subset.append(area_list[entry])
            else:
                pass
        subsets.append([x_subset, y_subset, size_subset, color, label])
    return subsets



def blank_substraction(df_werte, df_blank):
    identified_peak_counter = 0
    inaccuracy = 0.001  # ((201.001-201.000)/201) * 1 000 000 = 4.97ppm
    max_rt_deviation = 2
    max_peak_width = 4
    max_peak_width_deviation = 1
    # Calculate blank substracted area where blank m/z available
    df_werte["blank_subst_area"] = df_werte["area"]
    for mz_wert in range(len(df_werte["mz"].tolist())):
        blank_peak_counter = 0
        for mz_blank in range(len(df_blank["mz"].tolist())):
            if (df_werte["mz"][mz_wert] - inaccuracy < df_blank["mz"][mz_blank] < df_werte["mz"][
                mz_wert] + inaccuracy) and \
                    (df_werte["rt"][mz_wert] - max_rt_deviation < df_blank["rt"][mz_blank] < df_werte["rt"][
                        mz_wert] + max_rt_deviation) and \
                    ((df_werte["rt_range:max"][mz_wert] - df_werte["rt_range:min"][mz_wert]) <= max_peak_width) and \
                    ((df_blank["rt_range:max"][mz_blank] - df_blank["rt_range:min"][mz_blank]) <= max_peak_width) and \
                    (abs(((df_werte["rt_range:max"][mz_wert] - df_werte["rt_range:min"][mz_wert]) - (
                            df_blank["rt_range:max"][mz_blank] - df_blank["rt_range:min"][
                        mz_blank]))) <= max_peak_width_deviation):
                if blank_peak_counter >= 1:
                    print("blank_peak_counter > 1. There are more than one entries with the mass of " + str(
                        df_werte["mz"][mz_wert]) + " in the blank spectrum")
                blank_peak_counter = blank_peak_counter + 1
                identified_peak_counter = identified_peak_counter + 1
                blank_subst_value = df_werte["area"][mz_wert] - df_blank["area"][mz_blank]
                df_werte.at[mz_wert, "blank_subst_area"] = blank_subst_value
                print("m/z = " + str(df_werte["mz"][mz_wert]) + ": " + str(blank_subst_value))
    return df_werte


def get_approx_measurement_error(df_werte):
    # Calculate approximated error of measurements
    get_percentile = 0.95
    under_zero_list = []
    for element in range(len(df_werte["blank_subst_area"].tolist())):
        if df_werte["blank_subst_area"][element] < 0:
            under_zero_list.append(abs(df_werte["blank_subst_area"][element]))
    mean = sum(under_zero_list) / len(under_zero_list)
    variance = sum((x - mean) ** 2 for x in under_zero_list) / (len(under_zero_list) - 1)
    std_deviation = variance ** 0.5
    print("Standard deviation of under zeros: " + str(std_deviation))
    under_zero_list.sort()
    rank = get_percentile * (len(under_zero_list) - 1)
    result = under_zero_list[int(rank)]
    print(str(get_percentile * 100) + "´s percentile: " + str(result))

    # Initiate area_error with the x percentile:
    df_werte["area_error"] = result
    return df_werte


def detect_unlogical_masses(df_werte):
    # mark rows with unlogical masses (x,35 - x,75)
    for entry in range(len(df_werte["mz"].tolist())):
        exact_mass = df_werte["mz"][entry]
        if 0.35 <= (exact_mass - int(exact_mass)) <= 0.75:
            df_werte.at[entry, "logical_mass"] = "False"
        else:
            df_werte.at[entry, "logical_mass"] = "True"
    return df_werte


def create_df_with_valid_entries(df_werte):
    # create new df only with valid values
    reject_at_average_height = 5000
    reject_at_height_area_ratio = 8
    df_valid = df_werte.copy(deep=True)
    delete_list = []
    # Eliminate unlogical mass ranges
    # eliminate rows where area < 3*error
    # eliminate rows where blank_subst_area < 0
    # eliminate rows where the average peak height is under a certain limit.
    # eliminate rows where the ratio of height/area is under a certain limit
    for element in range(len(df_werte["mz"].tolist())):
        if df_werte["logical_mass"][element] == "False" or \
                df_werte["area"][element] <= 3 * df_werte["area_error"][element] or \
                df_werte["blank_subst_area"][element] < 0 or \
                df_werte["area"][element] / (
                df_werte["rt_range:max"][element] - df_werte["rt_range:min"][element]) < reject_at_average_height or \
                df_werte["height"][element] / df_werte["area"][element] < reject_at_height_area_ratio:
            delete_list.append(element)
    df_valid.drop(delete_list, axis=0, inplace=True)
    df_valid.index = list(range(0, len(df_valid), 1))
    print(f"Dataframe was reduced from {str(len(df_werte))} to {str(len(df_valid))} entries.")
    return df_valid


def delete_rows_without_formula_approximation(df_valid):
    # Delete rows where no formula could be found
    delete_list = []
    for element in range(len(df_valid["mz"].tolist())):
        if df_valid["formula_approximations"][element] == "{}":
            delete_list.append(element)
    df_valid.drop(delete_list, axis=0, inplace=True)
    df_valid.index = list(range(0, len(df_valid), 1))
    return df_valid


def get_best_formula_approximation_into_dict_format(df_valid):
    #Get best formula approximation into dict format
    x_val_list = []
    y_val_list = []
    area = []
    for index, row in df_valid.iterrows():
        print(row["formula_approximations"])
        print(row["mz"])
        formula_dict = ast.literal_eval(row["formula_approximations"])
        best_approx = max(formula_dict, key=formula_dict.get)
        score = max(formula_dict.values())
        atom_dict = {}
        for entry in range(len(best_approx)):
            atom_name = ""
            atom_count = ""
            iterator = 0
            if not best_approx[entry - 1].isnumeric():
                continue
            if best_approx[entry].isnumeric():
                continue
            if (not best_approx[entry].isnumeric()) and (best_approx[entry + 1].isnumeric()):
                atom_name = best_approx[entry]
                try:
                    while best_approx[entry + 1 + iterator].isnumeric():
                        atom_count = atom_count + best_approx[entry + 1 + iterator]
                        iterator = iterator + 1
                except:
                    pass
            elif (not best_approx[entry].isnumeric()) and (not best_approx[entry + 1].isnumeric()):
                atom_name = best_approx[entry] + best_approx[entry + 1]
                try:
                    while best_approx[entry + 2 + iterator].isnumeric():
                        atom_count = atom_count + best_approx[entry + 2 + iterator]
                        iterator = iterator + 1
                except:
                    pass
            atom_dict[atom_name] = int(atom_count)
        df_valid.at[index, "best_formula_approximation"] = str(atom_dict)
        df_valid.at[index, "best_formula_approximation_score"] = score
        df_valid.at[index, "best_formula_approximation_score_best"] = score
        df_valid.at[index, "best_formula_approximation_score_worst"] = score
    return df_valid


def combine_entries_with_same_formula(df_valid):
    # Combine entries with the same best formula
    agg_functions_2 = {}
    for entry in df_valid.columns:
        if entry == "blank_subst_area":
            agg_functions_2[entry] = "sum"
            continue
        elif entry == "area":
            agg_functions_2[entry] = "sum"
            continue
        elif entry == "area_error":
            agg_functions_2[entry] = "mean"
            continue
        elif entry == "rt":
            agg_functions_2[entry] = "mean"
            continue
        elif entry == "best_formula_approximation_score_best":
            agg_functions_2[entry] = "max"
        elif entry == "best_formula_approximation_score_worst":
            agg_functions_2[entry] = "min"
        else:
            agg_functions_2[entry] = "max"
    old_length = len(df_valid)
    df_valid = df_valid.groupby(df_valid["best_formula_approximation"]).aggregate(agg_functions_2)
    df_valid.index = list(range(len(df_valid)))
    print(f"A total of {str(round(old_length - len(df_valid)))} entries were combined.")
    return df_valid





if __name__ == "__main__":
    df_blank = pd.read_csv("Blank.txt")
    df_werte = pd.read_csv("Nitrophenole.txt")
    print(df_blank.columns)


    df_werte = blank_substraction(df_werte, df_blank)

    df_werte = get_approx_measurement_error(df_werte)


    #Round mz values
    max_mz_deviation = 0.05
    print(len(df_werte["mz"].unique().tolist()))
    for entry in range(len(df_werte["mz"].tolist())):
        df_werte.at[entry, "mz_rounded"] = max_mz_deviation * round(df_werte["mz"][entry] / max_mz_deviation)
    print(len(df_werte["mz_rounded"].unique().tolist()))


    #Round rt values
    max_rt_deviation = 0.1
    print(len(df_werte["rt"].unique().tolist()))
    for entry in range(len(df_werte["rt"].tolist())):
        df_werte.at[entry, "rt_rounded"] = max_rt_deviation * round(df_werte["rt"][entry] / max_rt_deviation)
    print(len(df_werte["rt_rounded"].unique().tolist()))


    df_werte = detect_unlogical_masses(df_werte)

    df_valid = create_df_with_valid_entries(df_werte)



    #Calculate formula approximations for every row
    exact_masses = []
    for entry in range(len(df_valid["mz"].tolist())):
        exact_masses.append(df_valid["mz"][entry] + 1.007825)

    mp = MyMultiprocessing(only_multithread=False)
    for exact_mass in exact_masses:
        mp.add_function(get_one_formula, [exact_mass])
        print("Function added")
    print("Starting...")
    mp.run_pool(wait_for_finish=True)
    mp.get_return_values()
    time.sleep(5)
    formula_dicts = mp.get_return_values()
    print("Ending.")
    for formula_dict in formula_dicts:
        for row_index in range(len(df_valid["mz"].tolist())):
            if (round(df_valid["mz"][row_index]+1.007825, 6)) == (round(formula_dict[0], 6)):
                df_valid.at[row_index, "formula_approximations"] = str(formula_dict[1])
            else:
                pass



    df_valid = delete_rows_without_formula_approximation(df_valid)

    df_valid = get_best_formula_approximation_into_dict_format(df_valid)


    df_valid = combine_entries_with_same_formula(df_valid)



    df_valid.to_csv("NenaTailfingen_Processed.csv")


    #Create bar chart of absolute mz intensities
    # fig, ax = plt.subplots(1, 2, layout="tight")
    # avg_mass_list = []
    # abs_area_of_mass_list = []
    # for avg_mass in df_werte["mz"].unique().tolist():
    #     avg_mass_list.append(avg_mass)
    #     for index, row in df_werte.iterrows():
    #         if row["mz"] == avg_mass:
    #             abs_area_of_mass_list.append(row["area"])
    #             break
    #         else:
    #             pass
    # ax[0].bar(avg_mass_list, abs_area_of_mass_list)
    #
    # avg_mass_list = []
    # abs_area_of_mass_list = []
    # for avg_mass in df_valid["mz"].unique().tolist():
    #     avg_mass_list.append(avg_mass)
    #     for index, row in df_valid.iterrows():
    #         if row["mz"] == avg_mass:
    #             abs_area_of_mass_list.append(row["area"])
    #             break
    #         else:
    #             pass
    # ax[1].bar(avg_mass_list, abs_area_of_mass_list)
    # plt.show()


    #Create bar chart with compared areas of df_valid and df_werte
    fig, ax = plt.subplots(2, 1, layout="constrained")
    total_area_df_werte = df_werte["blank_subst_area"].sum()
    summed_error_df_werte = df_werte["area_error"].sum()
    total_area_df_valid = df_valid["blank_subst_area"].sum()
    summed_error_df_valid = df_valid["area_error"].sum()
    ax[0].title.set_text("Total area")
    bar0 = ax[0].bar(["df_werte", "df_valid"], [total_area_df_werte, total_area_df_valid], yerr=[summed_error_df_werte, summed_error_df_valid], alpha=0.5, ecolor="black", width=0.7)
    ax[0].bar_label(bar0)
    ax[0].set_ylim(0, 1.5*total_area_df_werte)
    ax[0].grid(axis="y")
    ax[1].title.set_text("Entries")
    bar1 = ax[1].bar(["df_werte", "df_valid"], [len(df_werte), len(df_valid)], alpha=0.5, ecolor="black", width=0.7)
    ax[1].bar_label(bar1)
    ax[1].set_ylim(0, 1.5*len(df_werte))
    ax[1].grid(axis="y")
    plt.show()



    #van krevelen plot
    cos_ratio = []
    hc_ratio = []
    oc_ratio = []
    n_list = []
    o_list = []
    c_list = []
    nc_ratio = []
    mcr_list = []
    mz_list = []
    score_list = []
    blank_subst_area_list = []
    area_list = []
    label_list = []
    rt_list = []
    ppm_deviation_list = []

    for index, row in df_valid.iterrows():
        formula_dict = ast.literal_eval(row["best_formula_approximation"])
        try:
            n_atoms = formula_dict["N"]
        except:
            n_atoms = 0
        try:
            o_atoms = formula_dict["O"]
        except:
            o_atoms = 0
        try:
            c_atoms = formula_dict["C"]
        except:
            c_atoms = 0
        try:
            h_atoms = formula_dict["H"]
        except:
            h_atoms = 0
        try:
            s_atoms = formula_dict["S"]
        except:
            s_atoms = 0
        try:
            p_atoms = formula_dict["P"]
        except:
            p_atoms = 0
        try:
            cl_atoms = formula_dict["Cl"]
        except:
            cl_atoms = 0

        cos_ratio.append(2 * (o_atoms / c_atoms) - (h_atoms / c_atoms))
        try:
            mcr_list.append((1 + c_atoms - 0.5*h_atoms + 0.5*n_atoms) / o_atoms)
        except ZeroDivisionError:
            mcr_list.append(1)
        hc_ratio.append(h_atoms / c_atoms)
        oc_ratio.append(o_atoms / c_atoms)
        nc_ratio.append(n_atoms / c_atoms)
        n_list.append(n_atoms)
        o_list.append(o_atoms)
        c_list.append(c_atoms)
        mz_list.append(row["mz"])
        area_list.append(row["blank_subst_area"])
        blank_subst_area_list.append(row["blank_subst_area"])
        score_list.append(row["best_formula_approximation_score_best"])
        ppm_deviation_list.append(1/row["best_formula_approximation_score_best"])
        rt_list.append(row["rt"])

        if (n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CH")
        elif (n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHO")
        elif (not n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHN")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNO")
###
        elif (n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHS")
        elif (n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHOS")
        elif (not n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNS")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNOS")
###
        elif (n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHP")
        elif (n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHOP")
        elif (not n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNP")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNOP")

        elif (n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHSP")
        elif (n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHOSP")
        elif (not n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNSP")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (cl_atoms == 0):
            label_list.append("CHNOSP")
###
        elif (n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHCl")
        elif (n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHOCl")
        elif (not n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNCl")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNOCl")

        elif (n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHSCl")
        elif (n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHOSCl")
        elif (not n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNSCl")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNOSCl")

        elif (n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHPCl")
        elif (n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHOPCl")
        elif (not n_atoms == 0) and (o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNPCl")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNOPCl")

        elif (n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHSPCl")
        elif (n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHOSPCl")
        elif (not n_atoms == 0) and (o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNSPCl")
        elif (not n_atoms == 0) and (not o_atoms == 0) and (not s_atoms == 0) and (not p_atoms == 0) and (not cl_atoms == 0):
            label_list.append("CHNOSPCl")

        else:
            label_list.append("else")

    max_area_list = max(area_list)
    for element in range(len(area_list)):
        area_list[element] = ((area_list[element] / max_area_list)**1.5) * 100




    #######################################################################
    #######################################################################
    rows = 2
    cols = 4
    fig, ax = plt.subplots(rows, cols, figsize=[cols*5, rows*5])
    #plt.subplot_tool()
    plt.subplots_adjust(left=0.045, right=0.97, bottom=0.125, top=0.90, wspace=0.25, hspace=0.25)
    sc = []
    for axrow in range(rows):
        sc.append([])
        for axentry in range(cols):
            sc[axrow].append(0)
    legends = []

    subsets = prepare_subsets(nc_ratio, hc_ratio, area_list, label_list)
    sc[0][0] = (ax[0][0].scatter(nc_ratio, hc_ratio, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[0][0].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[0][0].set_xlabel("N/C ratio")
    ax[0][0].set_ylabel("H/C ratio")
    ax[0][0].set_title("H/C vs. N/C Plot")
    legends.append(ax[0][0].legend())


    subsets = prepare_subsets(c_list, cos_ratio, area_list, label_list)
    sc[0][1] = (ax[0][1].scatter(c_list, cos_ratio, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[0][1].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[0][1].arrow(15, -2, -10, 3, alpha=0.2, width=0.5, color="lightblue")
    ax[0][1].text(15, -2, "way of organics in atmosphere", color="lightblue")
    ax[0][1].invert_xaxis()
    ax[0][1].set_xlabel("C atom number")
    ax[0][1].set_ylabel("COS")
    ax[0][1].set_title("COS vs. C atom number Plot")
    legends.append(ax[0][1].legend())


    subsets = prepare_subsets(oc_ratio, hc_ratio, area_list, label_list)
    sc[0][2] = (ax[0][2].scatter(oc_ratio, hc_ratio, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[0][2].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.75))

    box_props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
    text = ("very highly oxidized \n" +
            "highly oxidized \n" +
            "intermediately oxidized \n" +
            "oxidized unsaturated \n" +
            "highly unsaturated")
    ax[0][2].text(0.98, 0.98, text, transform=ax[0][2].transAxes, fontsize=9,
                  verticalalignment="top", bbox=box_props, ha="right", va="top")

    ax[0][2].plot([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 0.35*(max(oc_ratio)/1.2)], color="black")
    ax[0][2].fill_between([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 0.35*(max(oc_ratio)/1.2)], [max(hc_ratio), max(hc_ratio)], alpha=0.2)

    ax[0][2].plot([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 1.1*(max(oc_ratio)/1.2)], color="black")
    ax[0][2].fill_between([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 1.1*(max(oc_ratio)/1.2)], [2.2, 2.2 - 0.35*(max(oc_ratio)/1.2)], alpha=0.2)

    ax[0][2].plot([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 2.2*(max(oc_ratio)/1.2)], color="black")
    ax[0][2].fill_between([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 2.2*(max(oc_ratio)/1.2)], [2.2, 2.2 - 1.1*(max(oc_ratio)/1.2)], alpha=0.2)

    ax[0][2].plot([min(oc_ratio), 0.8], [1.6, 0], color="black")
    ax[0][2].fill_between([min(oc_ratio), max(oc_ratio)], [2.2, 2.2 - 2.2*(max(oc_ratio)/1.2)], [1.6, 1.6-1.6*(max(oc_ratio)/0.8)], alpha=0.2)

    ax[0][2].fill_between([min(oc_ratio), max(oc_ratio)], [1.6, 1.6 - 1.6 * (max(oc_ratio) / 0.8)], [1.6, -99999999], alpha=0.2)

    ax[0][2].set_ylim([0, 3])
    ax[0][2].set_xlabel("O/C ratio")
    ax[0][2].set_ylabel("H/C ratio")
    ax[0][2].set_title("H/C vs. O/C Plot")
    legends.append(ax[0][2].legend(loc="upper left"))



    subsets = prepare_subsets(oc_ratio, hc_ratio, area_list, label_list)
    sc[0][3] = (ax[0][3].scatter(oc_ratio, hc_ratio, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[0][3].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[0][3].set_xlabel("O/C ratio")
    ax[0][3].set_ylabel("H/C ratio")
    ax[0][3].set_title("H/C vs. O/C Plot")
    legends.append(ax[0][3].legend())


    subsets = prepare_subsets(c_list, mcr_list, area_list, label_list)
    sc[1][0] = (ax[1][0].scatter(c_list, mcr_list, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[1][0].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[1][0].set_xlabel("C atom number")
    ax[1][0].set_ylabel("MCR maximum carbonyl ratio")
    ax[1][0].set_title("MCR vs. C atom number Plot")
    legends.append(ax[1][0].legend())


    subsets = prepare_subsets(mz_list, ppm_deviation_list, area_list, label_list)
    sc[1][1] = (ax[1][1].scatter(mz_list, ppm_deviation_list, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[1][1].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[1][1].set_xlabel("mz / u")
    ax[1][1].set_ylabel("delta mass / ppm")
    ax[1][1].set_title("Mass accuracy vs. absolute mass")
    ax[1][1].set_yscale('log')
    legends.append(ax[1][1].legend())


    subsets = prepare_subsets(mz_list, rt_list, area_list, label_list)
    sc[1][2] = (ax[1][2].scatter(mz_list, rt_list, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[1][2].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[1][2].set_xlabel("mz / u")
    ax[1][2].set_ylabel("rt / min")
    ax[1][2].set_title("rt vs. absolute mass Plot")
    legends.append(ax[1][2].legend())


    subsets = prepare_subsets(blank_subst_area_list, c_list, area_list, label_list)
    piechart_dict = {}
    print(subsets)
    for entry in subsets:
        piechart_dict[entry[4]] = sum(entry[0])
    total_area = sum(piechart_dict.values())
    dont_show_if_under_percent = 0.05
    delete_list = []
    others_area = 0
    print(piechart_dict)
    for entry in range(len(piechart_dict)):
        if piechart_dict[list(piechart_dict.keys())[entry]] <= (total_area * dont_show_if_under_percent):
            print(piechart_dict[list(piechart_dict.keys())[entry]])
            print(total_area * dont_show_if_under_percent)
            print(list(piechart_dict.keys())[entry])
            print(others_area)
            print()
            others_area = others_area + piechart_dict[list(piechart_dict.keys())[entry]]
            delete_list.append(list(piechart_dict.keys())[entry])
    for entry in delete_list:
        print(entry)
        del piechart_dict[entry]
    piechart_dict["others"] = others_area
    print(piechart_dict)
    ax[1][3].pie(piechart_dict.values(), labels=piechart_dict.keys(), autopct="%.1f%%",
                 wedgeprops={"linewidth": 3.0, "edgecolor": "white"})
    ax[1][3].set_title("Subset areas")
    legends.append(ax[1][3].legend())
    ax[1][3].get_legend().remove()



    for curr_legend in legends:
        for legend_entry in curr_legend.legend_handles:
            legend_entry._sizes = [30]

    annot = []
    for axrow in range(len(ax)):
        annot.append([])
        for axentry in range(len(ax[axrow])):
            annot[axrow].append(0)
    for ax_row in range(len(ax)):
        for ax_entry in range(len(ax[ax_row])):
            annot[ax_row][ax_entry] = (ax[ax_row][ax_entry].annotate("", xy=(0, 0), xytext=(10, 10), xycoords="data", textcoords="axes pixels",
                                        bbox=dict(boxstyle="round, pad=0.2", fc="w"), arrowprops=dict(arrowstyle="->")))
    for annot_row in range(len(annot)):
        for annot_entry in range(len(annot[annot_row])):
            annot[annot_row][annot_entry].set_visible(False)

    names = []
    for index, row in df_valid.iterrows():
        formula_dict = ast.literal_eval(row["best_formula_approximation"])
        atom_names = list(formula_dict.keys())
        atom_numbers = list(formula_dict.values())
        molecule_string = ""
        for element in range(len(atom_names)):
            molecule_string = molecule_string + str(atom_names[element]) + str(atom_numbers[element])
        print(molecule_string)
        text = "Formula: " + str(molecule_string) + "\n"
        placing = sorted(df_valid["blank_subst_area"].unique().tolist(), reverse=True).index(row["blank_subst_area"]) + 1
        text = text + "Area=" + str(round(row["blank_subst_area"])) + "    " + str(placing) + "/" + str(len(df_valid)) + "\n"
        text = text + "m/z=" + str(round(row["mz"], 4)) + "    " + str(round(row["best_formula_approximation_score_best"], 3)) + " > " + str(round(row["best_formula_approximation_score_worst"], 3)) + "\n"
        text = text + "rt=" + str(round(row["rt"], 2)) + "    Height/Area:" + str(round((row["height"] / row["area"]), 1)) + "\n"
        names.append(text)




    def update_annot(ind):
        for row in range(len(annot)):
            for entry in range(len(annot[row])):
                try:
                    pos = sc[row][entry].get_offsets()[ind["ind"][0]]
                except:
                    return
                annot[row][entry].xy = pos
                #text = "{}: {}".format(" ".join(list(map(str, ind["ind"])))," ".join([names[n] for n in ind["ind"]]))

                text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
                text = text.strip()
                annot[row][entry].set_text(text)
                #annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
                annot[row][entry].get_bbox_patch().set_alpha(0.4)

    def hover(event):
        for row in range(len(annot)):
            for entry in range(len(annot[row])):

                vis = annot[row][entry].get_visible()
                if event.inaxes == ax[row][entry]:
                    try:
                        cont, ind = sc[row][entry].contains(event)
                    except:
                        return
                    if cont:
                        update_annot(ind)
                        annot[row][entry].set_visible(True)
                        fig.canvas.draw_idle()
                    else:
                        if vis:
                            annot[row][entry].set_visible(False)
                            fig.canvas.draw_idle()


    fig.canvas.mpl_connect("motion_notify_event", hover)

    plt.show()




