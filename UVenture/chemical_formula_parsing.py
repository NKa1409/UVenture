import ast


def get_formula_to_dict(formula_string):
    if formula_string == "":
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


def get_formula_string_from_dict(formula_dict):
    if isinstance(formula_dict, str):
        try:
            eventual_formula_dict = ast.literal_eval(formula_dict)
            if isinstance(eventual_formula_dict, dict):
                formula_dict = eventual_formula_dict
            else:
                pass
        except:
            pass
    if isinstance(formula_dict, str):
        try:
            formula_dict = get_formula_to_dict(formula_dict)
        except:
            pass
    if len(formula_dict) == 0:
        return ""
    return "".join([str(a) + str(n) for a, n in formula_dict.items()])