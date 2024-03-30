import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MS_functions



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


def get_formula_into_dict_format(df_werte, formula_approximation_colname):
    print(df_werte[formula_approximation_colname])
    for index, row in df_werte.iterrows():
        formula = row[formula_approximation_colname]
        if str(formula) == "nan":
            print("Formula was nan")
            continue
        formula_new = ""
        if not formula[-1].isdigit():
            formula = formula + "1"
        for entry in range(len(formula)):
            if entry == len(formula):
                break
            if not formula[entry].isdigit() and formula[entry+1].isupper():
                formula_new = formula_new + formula[entry] + str(1)
            else:
                formula_new = formula_new + formula[entry]
        df_werte.at[index, formula_approximation_colname] = formula_new
    print(df_werte[formula_approximation_colname])

    df_werte["best_formula_approximation"] = pd.Series()
    print(df_werte)
    for index, row in df_werte.iterrows():
        formula = row[formula_approximation_colname]
        if str(formula).lower() == "nan":
            df_werte.at[index, "best_formula_approximation"] = {}
            continue
        atom_dict = {}
        for entry in range(len(formula)):
            atom_name = ""
            atom_count = ""
            iterator = 0
            if not formula[entry - 1].isnumeric():
                continue
            if formula[entry].isnumeric():
                continue
            if (not formula[entry].isnumeric()) and (formula[entry + 1].isnumeric()):
                atom_name = formula[entry]
                try:
                    while formula[entry + 1 + iterator].isnumeric():
                        atom_count = atom_count + formula[entry + 1 + iterator]
                        iterator = iterator + 1
                except:
                    pass
            elif (not formula[entry].isnumeric()) and (not formula[entry + 1].isnumeric()):
                atom_name = formula[entry] + formula[entry + 1]
                try:
                    while formula[entry + 2 + iterator].isnumeric():
                        atom_count = atom_count + formula[entry + 2 + iterator]
                        iterator = iterator + 1
                except:
                    pass
            if atom_count == "":
                atom_count = 1
            else:
                pass
            atom_dict[atom_name] = int(atom_count)
        print(atom_dict)
        df_werte.at[index, "best_formula_approximation"] = str(atom_dict)
    return df_werte


def calc_blank_subst_area(df_werte, blank_area_colname, sample_area_colname, blank_subst_area_colname):
    for index, row in df_werte.iterrows():
        try:
            df_werte.at[index, blank_subst_area_colname] = row[sample_area_colname] - row[blank_area_colname]
        except:
            df_werte.at[index, blank_subst_area_colname] = 0
    return df_werte


def filter_where_blank_subst_area_under_zero(df_werte, blank_subst_area_colname, sample_area_colname, blank_area_colname):
    delete_row_list = []
    for index, row in df_werte.iterrows():
        if row[blank_subst_area_colname] <= 0:
            delete_row_list.append(index)
        else:
            pass
    df_werte = df_werte.drop(delete_row_list)
    df_werte.index = list(range(len(df_werte)))
    return df_werte


if __name__ == "__main__":
    summary_filepath="U://MyFolder//MONOTONS//Filtermessungen//ACROSS Filter//Neue Messungen//F7//F7_HRAIF_1_outputfolder//SUMMARY.txt"
    df = MS_functions.read_summary_to_df(summary_filepath, sep="\t", colnames=("MASS", "FORMULA", "INDEX", "RT", "INTENSITY", "SCORE", "FRAGS", "F_SCORES", "F_INTENSITIES", "NLS", "NLS_DEV"))



    df_werte = pd.read_csv("CRDS_MZmineAusertung.txt")

    formula_approximation_colname = "row identity (main ID)"
    row_mz_colname = "row m/z"
    blank_area_colname = "FilterBlank_keinAIF_NEG_APCI_PFPSaeule.raw Peak area"
    sample_area_colname = "FilterCRDS_keinAIF_NEG_APCI_PFPSaeule.raw Peak area"
    blank_subst_area_colname = "blank_subst_area"
    row_rt_colname = "row retention time"
    formula_approx_dict_colname = "best_formula_approximation"





    print("START")
    df_werte = get_formula_into_dict_format(df_werte, formula_approximation_colname)
    print("111111111111111111111")
    df_werte = calc_blank_subst_area(df_werte, blank_area_colname, sample_area_colname, blank_subst_area_colname)
    print("2222222222222222222222")
    df_werte = filter_where_blank_subst_area_under_zero(df_werte, blank_subst_area_colname, sample_area_colname, blank_area_colname)
    print("333333333333333")

    #Create bar chart of absolute mz intensities
    fig, ax = plt.subplots(1, 1, layout="tight")
    avg_mass_list = []
    abs_area_of_mass_list = []
    for avg_mass in df_werte[row_mz_colname].unique().tolist():
        avg_mass_list.append(avg_mass)
        for index, row in df_werte.iterrows():
            if row[row_mz_colname] == avg_mass:
                abs_area_of_mass_list.append(row[blank_subst_area_colname])
                break
            else:
                pass
    ax.bar(avg_mass_list, abs_area_of_mass_list)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    plt.show()




    #Plots
    cos_ratio = []
    hc_ratio = []
    oc_ratio = []
    n_list = []
    o_list = []
    c_list = []
    nc_ratio = []
    mcr_list = []
    mz_list = []
    blank_subst_area_list = []
    area_list = []
    label_list = []
    rt_list = []

    for index, row in df_werte.iterrows():
        print(row[formula_approx_dict_colname])
        try:
            formula_dict = ast.literal_eval(row[formula_approx_dict_colname])
        except:
            formula_dict = {}
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

        try:
            cos_ratio.append(2 * (o_atoms / c_atoms) - (h_atoms / c_atoms))
        except:
            cos_ratio.append(None)
        try:
            mcr_list.append((1 + c_atoms - 0.5*h_atoms + 0.5*n_atoms) / o_atoms)
        except ZeroDivisionError:
            mcr_list.append(None)
        try:
            hc_ratio.append(h_atoms / c_atoms)
        except:
            hc_ratio.append(None)
        try:
            oc_ratio.append(o_atoms / c_atoms)
        except:
            oc_ratio.append(None)
        try:
            nc_ratio.append(n_atoms / c_atoms)
        except:
            nc_ratio.append(None)
        n_list.append(n_atoms)
        o_list.append(o_atoms)
        c_list.append(c_atoms)
        mz_list.append(row[row_mz_colname])
        area_list.append(row[blank_subst_area_colname])
        blank_subst_area_list.append(row[blank_subst_area_colname])
        rt_list.append(row[row_rt_colname])

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
        area_list[element] = ((area_list[element] / max_area_list)) * 150


    print(nc_ratio)
    print(hc_ratio)
    print(area_list)
    for entry in range(len(area_list)):
        if isinstance(area_list[entry], complex):
            area_list[entry] = area_list[entry].real
    print(area_list)

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

    subsets = prepare_subsets(oc_ratio, hc_ratio, area_list, label_list)
    sc[0][0] = (ax[0][0].scatter(oc_ratio, hc_ratio, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[0][0].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[0][0].set_xlabel("O/C ratio")
    ax[0][0].set_ylabel("H/C ratio")
    ax[0][0].set_title("H/C vs. O/C Plot")
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

    try:
        min_oc_ratio = min(oc_ratio)
    except TypeError:
        min_oc_ratio = 0
    try:
        max_oc_ratio = max(oc_ratio)
    except TypeError:
        max_oc_ratio = 2
    try:
        max_hc_ratio = max(hc_ratio)
    except TypeError:
        max_hc_ratio = 2.5
    try:
        min_hc_ratio = min(hc_ratio)
    except TypeError:
        min_hc_ratio = 0
    ax[0][2].plot([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 0.35*(max_oc_ratio/1.2)], color="black")
    ax[0][2].fill_between([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 0.35*(max_oc_ratio/1.2)], [max_hc_ratio, max_hc_ratio], alpha=0.2)

    ax[0][2].plot([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 1.1*(max_oc_ratio/1.2)], color="black")
    ax[0][2].fill_between([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 1.1*(max_oc_ratio/1.2)], [2.2, 2.2 - 0.35*(max_oc_ratio/1.2)], alpha=0.2)

    ax[0][2].plot([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 2.2*(max_oc_ratio/1.2)], color="black")
    ax[0][2].fill_between([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 2.2*(max_oc_ratio/1.2)], [2.2, 2.2 - 1.1*(max_oc_ratio/1.2)], alpha=0.2)

    ax[0][2].plot([min_oc_ratio, 0.8], [1.6, 0], color="black")
    ax[0][2].fill_between([min_oc_ratio, max_oc_ratio], [2.2, 2.2 - 2.2*(max_oc_ratio/1.2)], [1.6, 1.6-1.6*(max_oc_ratio/0.8)], alpha=0.2)

    ax[0][2].fill_between([min_oc_ratio, max_oc_ratio], [1.6, 1.6 - 1.6 * (max_oc_ratio / 0.8)], [1.6, -99999999], alpha=0.2)

    ax[0][2].set_ylim([0, max_hc_ratio+0.1])
    ax[0][2].set_xlim([min_oc_ratio-0.1, max_oc_ratio+0.1])
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

    subsets = prepare_subsets(c_list, mcr_list, area_list, label_list)
    sc[1][1] = (ax[1][1].scatter(c_list, mcr_list, s=area_list, alpha=0.5))
    for entry in subsets:
        (ax[1][1].scatter(entry[0], entry[1], s=entry[2], c=entry[3], label=entry[4], alpha=0.5))
    ax[1][1].set_xlabel("C atom number")
    ax[1][1].set_ylabel("MCR maximum carbonyl ratio")
    ax[1][1].set_title("MCR vs. C atom number Plot")
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
    for index, row in df_werte.iterrows():
        try:
            formula_dict = ast.literal_eval(row[formula_approx_dict_colname])
        except:
            formula_dict = {}
        atom_names = list(formula_dict.keys())
        atom_numbers = list(formula_dict.values())
        molecule_string = ""
        for element in range(len(atom_names)):
            molecule_string = molecule_string + str(atom_names[element]) + str(atom_numbers[element])
        print(molecule_string)
        text = "Formula: " + str(molecule_string) + "\n"
        placing = sorted(df_werte[blank_subst_area_colname].unique().tolist(), reverse=True).index(row[blank_subst_area_colname]) + 1
        text = text + "Area=" + str(round(row[blank_subst_area_colname])) + "    " + str(placing) + "/" + str(len(df_werte)) + "\n"
        text = text + "m/z=" + str(round(row[row_mz_colname], 4)) + "\n"
        text = text + "rt=" + str(round(row[row_rt_colname], 2)) + "\n"
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




