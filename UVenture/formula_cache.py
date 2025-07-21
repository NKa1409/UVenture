import gzip
import traceback
from UVenture.chemical_formula_parsing import get_formula_to_dict, get_formula_string_from_dict
from UVenture.formula_calculations import calculate_likelyhood_of_formula_dict


def get_formula_from_cache(cache_folder, mass, max_ppm_dev, max_unlikeliness_of_formula=10, min_unlikelyness_of_formula=0, file_basename="formulas_for_mass_", filetype=".txt.gz", max_mass=1200, remove_unlogical_formulas=False, debug_output=True):
    if debug_output == True:
        print("Starting formula generation from cache...")
    if mass >= max_mass:
        #on the fly prediction
        result = get_one_formula(mass, max_dev=0.005, charge_of_measured_mass=0, dbe_range=(-10, 150), c_oxidation_state_range=(-4, 4),
                    heteroatom_ratio_range=(0, 99999), aromaticity_eq_range=(-20, 100), h_c_ratio=(0, 4),
                    max_h=100, max_c=30, max_n=30, max_o=20, max_s=3, max_p=2, max_cl=4, max_br=2, max_i=0, print_progress=True, return_cache_dict=False)
        pred_dict = result[1]

        #lookup in cache files
        files_to_search = [cache_folder + file_basename + str(int(mass)) + filetype,
                           cache_folder + file_basename + str(int(mass) - 1) + filetype,
                           cache_folder + file_basename + str(int(mass) + 1) + filetype]
        all_dict = {}
        for file in files_to_search:
            try:
                with open(file, "r") as f:
                    lines = f.readlines()
                    formulas = []
                    deviations = []
                    for line_index in range(len(lines)):
                        lines[line_index] = lines[line_index].strip()
                        if lines[line_index] == "":
                            continue
                        formulas.append(lines[line_index].split(",")[0])
                        deviations.append((((float(lines[line_index].split(",")[-1])) - mass) / mass) * 1000000)
                    mydict = dict(zip(formulas, deviations))
                all_dict.update(mydict)
            except Exception as e:
                print("Error with formula determination! " + str(e))
                print(traceback.format_exc())
                if debug_output == True:
                    print(all_dict)
                continue
        all_dict = {f: dev for f, dev in all_dict.items() if abs(dev) <= max_ppm_dev}
        all_dict = {k: v for k, v in sorted(all_dict.items(), key=lambda item: abs(item[1]), reverse=False)}
        all_dict.update(pred_dict)
    else:
        files_to_search = [cache_folder + file_basename + str(int(mass)) + filetype,
                           cache_folder + file_basename + str(int(mass)-1) + filetype,
                           cache_folder + file_basename + str(int(mass)+1) + filetype]
        all_dict = {}
        for file in files_to_search:
            try:
                with gzip.open(file, "rt") as f:
                    lines = f.readlines()
                    formulas = []
                    deviations = []
                    for line_index in range(len(lines)):
                        lines[line_index] = lines[line_index].strip()
                        if lines[line_index] == "":
                            continue
                        formulas.append(lines[line_index].split(",")[0])
                        deviations.append( (((float(lines[line_index].split(",")[-1])) - mass) / mass) * 1000000 )
                    mydict = dict(zip(formulas, deviations))
                all_dict.update(mydict)
            except Exception as e:
                print("Error with formula determination! " + str(e))
                print(traceback.format_exc())
                if debug_output == True:
                    print(all_dict)
                continue
        all_dict = {f: dev for f, dev in all_dict.items() if abs(dev) <= max_ppm_dev}
        all_dict = {k: v for k, v in sorted(all_dict.items(), key=lambda item: abs(item[1]), reverse=False)}
    formula_dict_list = [get_formula_to_dict(f) for f, dev in all_dict.items()]
    formula_dict_list = [dict(t) for t in {tuple(d.items()) for d in formula_dict_list}]
    formula_dict_list = [get_formula_string_from_dict(f) for f in formula_dict_list]
    all_dict = {f: dev for f, dev in all_dict.items() if f in formula_dict_list}

    if remove_unlogical_formulas == True:
        for formula in list(all_dict.keys()):
            curr_formula_dict = get_formula_to_dict(formula)
            f_likelyness = calculate_likelyhood_of_formula_dict(curr_formula_dict)
            f_unlikelyness = -1 * f_likelyness
            if f_unlikelyness >= max_unlikeliness_of_formula or f_unlikelyness < min_unlikelyness_of_formula:
                if debug_output == True:
                    print("Unlikelyness of current formula: " + str(f_unlikelyness) + "Formula: " + str(formula))
                del all_dict[formula]

    return [mass, all_dict]


def get_one_formula(measured_mass, max_dev=0.005, charge_of_measured_mass=-1, dbe_range=(-1, 15), c_oxidation_state_range=(-4, 4),
                    heteroatom_ratio_range=(0, 0.75), aromaticity_eq_range=(-20, 100), h_c_ratio=(0.3, 3),
                    max_h=100, max_c=30, max_n=30, max_o=20, max_s=5, max_p=3, max_cl=5, max_br=1, max_i=1, print_progress=True, return_cache_dict=False):
    mass_dict = {"H": 1.007825,
                 "B": 11.009305,
                 "C": 12.00000,
                 "N": 14.003074,
                 "O": 15.994915,
                 "F": 18.998403,
                 "Cl": 34.968853,
                 "P": 30.973763,
                 "S": 31.972072,
                 "Br": 80.916290,
                 "I": 126.904477,
                 "Fe": 55.934939,
                 "Si": 27.976928,
                 "As": 74.921596,
                 "Se": 79.916521}
    cache_dict = {}

    if not charge_of_measured_mass == 0:
        mass = abs(measured_mass * charge_of_measured_mass)
    else:
        mass = measured_mass

    mass = mass + (charge_of_measured_mass * 0.0005485)

    h_atoms = [0, max_h]
    c_atoms = [0, max_c]
    n_atoms = [0, max_n]
    o_atoms = [0, max_o]
    s_atoms = [0, max_s]
    p_atoms = [0, max_p]
    cl_atoms = [0, max_cl]
    br_atoms = [0, max_br]
    i_atoms = [0, max_i]
    dbe_range = dbe_range
    max_deviation_mz = max_dev
    formulas = {}
    if c_atoms[-1] > int(mass / mass_dict["C"]) + 1:
        c_atoms[-1] = int(mass / mass_dict["C"])
    if n_atoms[-1] > int(mass / mass_dict["N"]) + 1:
        n_atoms[-1] = int(mass / mass_dict["N"])
    if o_atoms[-1] > int(mass / mass_dict["O"]) + 1:
        o_atoms[-1] = int(mass / mass_dict["O"])
    if s_atoms[-1] > int(mass / mass_dict["S"]) + 1:
        s_atoms[-1] = int(mass / mass_dict["S"])
    if p_atoms[-1] > int(mass / mass_dict["P"]) + 1:
        p_atoms[-1] = int(mass / mass_dict["P"])
    if cl_atoms[-1] > int(mass / mass_dict["Cl"]) + 1:
        cl_atoms[-1] = int(mass / mass_dict["Cl"])
    if br_atoms[-1] > int(mass / mass_dict["Br"]) + 1:
        br_atoms[-1] = round(mass / mass_dict["Br"])
    if i_atoms[-1] > int(mass / mass_dict["I"]) + 1:
        i_atoms[-1] = round(mass / mass_dict["I"])
    if h_atoms[-1] > int(mass / mass_dict["H"]) + 1:
        h_atoms[-1] = int(mass / mass_dict["H"])
    if 0.25 <= (mass - int(mass)) <= 0.75:
        return [mass, formulas]
    for c_count in range(c_atoms[0], c_atoms[-1] + 1, 1):
        if print_progress:
            print(str(c_count) + " / " + str(c_atoms[-1]))
        for n_count in range(n_atoms[0], n_atoms[-1] + 1, 1):
            for br_count in range(br_atoms[0], br_atoms[1]+1, 1):
                for i_count in range(i_atoms[0], i_atoms[1]+1, 1):
                    for o_count in range(o_atoms[0], o_atoms[-1] + 1, 1):
                        for s_count in range(s_atoms[0], s_atoms[-1] + 1, 1):
                            for p_count in range(p_atoms[0], p_atoms[-1] + 1, 1):
                                for cl_count in range(cl_atoms[0], cl_atoms[-1] + 1, 1):
                                    for h_count in range(h_atoms[0], 1 + int(2 * (c_count + (2 / 3 * n_count) + (1 / 2 * s_count) + (1 / 2 + o_count)) + 2) + 1, 1):
                                        halogen_count = 0
                                        heavy_atom_count = 0
                                        if c_count == 0:
                                            current_heteroatom_ratio = 999
                                            current_cos = 0
                                            current_h_c_ratio = 999
                                        else:
                                            current_heteroatom_ratio = (n_count + o_count + s_count + p_count + cl_count + br_count + i_count) / (c_count + (n_count + o_count + s_count + p_count + cl_count + br_count + i_count))
                                            current_cos = 2 * (o_count / c_count) - (h_count / c_count)  # see https://doi.org/10.1016/j.atmosenv.2018.06.036
                                            current_h_c_ratio = h_count / c_count

                                        dbe = (c_count + 1) - ((h_count - n_count + halogen_count) / 2)
                                        try:
                                            aromaticity_equivalent = ((c_count - (h_count - c_count)) / (
                                                dbe)) + 1  # see https://doi.org/10.1002/rcm.7038
                                        except ZeroDivisionError:
                                            aromaticity_equivalent = 0
                                        if (dbe_range[0] <= dbe <= dbe_range[-1] and
                                                (c_oxidation_state_range[0] <= current_cos <= c_oxidation_state_range[-1]) and
                                                (heteroatom_ratio_range[0] <= current_heteroatom_ratio <= heteroatom_ratio_range[-1]) and
                                                (aromaticity_eq_range[0] <= aromaticity_equivalent <= aromaticity_eq_range[-1]) and
                                                (h_c_ratio[0] <= current_h_c_ratio <= h_c_ratio[-1])):
                                            calc_mass = h_count * mass_dict["H"] + \
                                                        c_count * mass_dict["C"] + \
                                                        n_count * mass_dict["N"] + \
                                                        o_count * mass_dict["O"] + \
                                                        s_count * mass_dict["S"] + \
                                                        p_count * mass_dict["P"] + \
                                                        cl_count * mass_dict["Cl"] + \
                                                        br_count * mass_dict["Br"] + \
                                                        i_count * mass_dict["I"]
                                            if calc_mass >= mass+2:
                                                continue
                                            if return_cache_dict:
                                                print_str = ""
                                                if c_count > 0:
                                                    print_str = print_str + "C" + str(c_count)
                                                if h_count > 0:
                                                    print_str = print_str + "H" + str(h_count)
                                                if n_count > 0:
                                                    print_str = print_str + "N" + str(n_count)
                                                if o_count > 0:
                                                    print_str = print_str + "O" + str(o_count)
                                                if s_count > 0:
                                                    print_str = print_str + "S" + str(s_count)
                                                if p_count > 0:
                                                    print_str = print_str + "P" + str(p_count)
                                                if cl_count > 0:
                                                    print_str = print_str + "Cl" + str(cl_count)
                                                if br_count > 0:
                                                    print_str = print_str + "Br" + str(br_count)
                                                if i_count > 0:
                                                    print_str = print_str + "I" + str(i_count)
                                                cache_dict[print_str] = calc_mass
                                            if calc_mass - max_deviation_mz <= mass <= calc_mass + max_deviation_mz:
                                                # calculate score
                                                score = 0
                                                ppm_deviation = ((calc_mass - mass) / mass) * 1000000
                                                if ppm_deviation == 0:
                                                    score = 99999
                                                else:
                                                    try:
                                                        score = (1 / ppm_deviation)
                                                    except ZeroDivisionError:
                                                        score = 999999
                                                print_str = ""
                                                if c_count > 0:
                                                    print_str = print_str + "C" + str(c_count)
                                                if h_count > 0:
                                                    print_str = print_str + "H" + str(h_count)
                                                if n_count > 0:
                                                    print_str = print_str + "N" + str(n_count)
                                                if o_count > 0:
                                                    print_str = print_str + "O" + str(o_count)
                                                if s_count > 0:
                                                    print_str = print_str + "S" + str(s_count)
                                                if p_count > 0:
                                                    print_str = print_str + "P" + str(p_count)
                                                if cl_count > 0:
                                                    print_str = print_str + "Cl" + str(cl_count)
                                                if br_count > 0:
                                                    print_str = print_str + "Br" + str(br_count)
                                                if i_count > 0:
                                                    print_str = print_str + "I" + str(i_count)
                                                formulas[print_str] = ppm_deviation
                                        else:
                                            pass

    formulas = {k: v for k, v in sorted(formulas.items(), key=lambda item: abs(item[1]), reverse=False)}
    if len(formulas) == 0:
        formulas = {"": 99999999}
    if return_cache_dict:
        cache_dict = {k: v for k, v in sorted(cache_dict.items(), key=lambda item: abs(item[1]), reverse=False)}
        return cache_dict
    return [mass, formulas]