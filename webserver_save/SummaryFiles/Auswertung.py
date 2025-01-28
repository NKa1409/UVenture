import pandas as pd
import matplotlib.pyplot as plt
import ast
from scipy.stats import gaussian_kde
import numpy as np
import traceback
import pyteomics.mass
import mplcursors
import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from matplotlib.cm import get_cmap
import copy
import os

def calculate_mass(formula):
    masse = pyteomics.mass.mass.calculate_mass(formula)
    return masse

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


def create_shaded_areas(ax, x_vals, y_vals, sizes_vals, colors_vals, subset_vals, proportion=0.5, number_of_groups=3):
    unique_labels = np.unique(np.array(subset_vals))
    subset_indices = {label: np.where(np.array(subset_vals) == label)[0] for label in unique_labels}

    lenghts_subsets = {k: len(v) for k, v in subset_indices.items()}
    lenghts_subsets = {k: v for k, v in sorted(lenghts_subsets.items(), key=lambda item: item[1], reverse=True)}
    top_subsets = {k: lenghts_subsets[k] for k in list(lenghts_subsets.keys())[:number_of_groups]}
    subset_indices = {k: v for k, v in subset_indices.items() if k in list(top_subsets.keys())}
    # subset_indices = {k: v for k, v in subset_indices.items() if k in ["CHO", "CHON", "CHOS", "CHOP", "COF"]}
    for label, indices in subset_indices.items():
        try:
            x_subset = np.array(x_vals)[indices]
            y_subset = np.array(y_vals)[indices]
            sizes_subset = np.array(sizes_vals)[indices]
            color_subset = np.array(colors_vals)[indices]
            color = color_subset[0]
            plot_kde_subset(ax, x_subset, y_subset, sizes_subset, color=color, label=label, proportion=proportion)
        except Exception as e:
            print(traceback.format_exc())

# Function to calculate KDE and plot the shaded area
def plot_kde_subset(ax, x_subset, y_subset, weights, label, color="red", proportion=0.5):
    data_subset = np.vstack([x_subset, y_subset])
    kde = gaussian_kde(data_subset, weights=weights)
    xmin, xmax = x_subset.min(), x_subset.max()
    ymin, ymax = y_subset.min(), y_subset.max()

    # Create a grid of points to evaluate the KDE on
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    density = kde(positions).reshape(xx.shape)

    # Find the value of the density that corresponds to the 50% level
    sorted_density = np.sort(density.ravel())
    cumulative_density = np.cumsum(sorted_density)
    cumulative_density /= cumulative_density[-1]
    density_threshold = sorted_density[np.searchsorted(cumulative_density, proportion)]

    # Overlay the shaded area for the 50% contour
    ax.contourf(xx, yy, density, levels=[density_threshold, density.max()], alpha=0.3, colors=color)

def get_cleaned_lists(l1, l2, size_list, subset_list, point_info, defined_subset_colors=True):
    size_multiplicator = 500
    if not len(l1) == len(l2):
        print("Length of the two provided lists is not the same!")
        print(l1)
        print(l2)
        print("Returning...")
        return
    rm_index_lst = []
    for i in range(len(l1)):
        if (not isinstance(l1[i], float)) and (not isinstance(l1[i], int)):
            rm_index_lst.append(i)
        if (not isinstance(l2[i], float)) and (not isinstance(l2[i], int)):
            rm_index_lst.append(i)
    rm_index_lst = list(set(rm_index_lst))
    for index in sorted(rm_index_lst, reverse=True):
        del l1[index]
        del l2[index]
        del subset_list[index]
        del size_list[index]
        del point_info[index]
    color_list = []
    if defined_subset_colors==True:
        for entry in subset_list:
            if entry == "CHO":
                color_list.append("blue")
            elif entry == "CHON":
                color_list.append("green")
            elif entry == "CHOS":
                color_list.append("red")
            elif "F" in entry or "Cl" in entry or "Br" in entry:
                color_list.append("yellow")
            elif "B_" in entry:
                color_list.append("gray")
            else:
                color_list.append("black")
    else:
        # Create a color map for the subsets
        #subset_colors = {subset: color for subset, color in zip(list(set(subset_list)), plt.cm.tab10.colors)}
        num_colors = len(list(set(subset_list)))
        cmap = get_cmap('tab20', num_colors)  # You can choose another colormap if you want
        # Create a color map for the subsets
        subset_colors = {subset: cmap(i) for i, subset in enumerate(list(set(subset_list)))}
        for entry in subset_list:
            color_list.append(subset_colors[entry])
    max_size_list = max(size_list)
    for i in range(len(size_list)):
        size_list[i] = (size_list[i] / max_size_list) * size_multiplicator
    point_info_list = point_info
    return l1, l2, size_list, color_list, subset_list, point_info_list

def get_cleaned_lists_3d(l1, l2, l3, size_list, subset_list, point_info, defined_subset_colors=True):
    size_multiplicator = 500
    if not len(l1) == len(l2):
        print("Length of the two provided lists is not the same!")
        print(l1)
        print(l2)
        print("Returning...")
        return
    rm_index_lst = []
    for i in range(len(l1)):
        if (not isinstance(l1[i], float)) and (not isinstance(l1[i], int)):
            rm_index_lst.append(i)
        if (not isinstance(l2[i], float)) and (not isinstance(l2[i], int)):
            rm_index_lst.append(i)
    rm_index_lst = list(set(rm_index_lst))
    for index in sorted(rm_index_lst, reverse=True):
        del l1[index]
        del l2[index]
        del l3[index]
        del subset_list[index]
        del size_list[index]
        del point_info[index]
    color_list = []
    if defined_subset_colors==True:
        for entry in subset_list:
            if entry == "CHO":
                color_list.append("blue")
            elif entry == "CHON":
                color_list.append("green")
            elif entry == "CHOS":
                color_list.append("red")
            elif "F" in entry or "Cl" in entry or "Br" in entry:
                color_list.append("yellow")
            elif "B_" in entry:
                color_list.append("gray")
            else:
                color_list.append("black")
    else:
        # Create a color map for the subsets
        #subset_colors = {subset: color for subset, color in zip(list(set(subset_list)), plt.cm.tab10.colors)}
        num_colors = len(list(set(subset_list)))
        cmap = get_cmap('tab20', num_colors)  # You can choose another colormap if you want
        # Create a color map for the subsets
        subset_colors = {subset: cmap(i) for i, subset in enumerate(list(set(subset_list)))}
        for entry in subset_list:
            color_list.append(subset_colors[entry])
    max_size_list = max(size_list)
    for i in range(len(size_list)):
        size_list[i] = (size_list[i] / max_size_list) * size_multiplicator
    point_info_list = point_info
    return l1, l2, l3, size_list, color_list, subset_list, point_info_list

def insert_interactive_scatter(scatter, info):
    cursor = mplcursors.cursor(scatter, hover=True)
    @cursor.connect("add")
    def on_add(sel):
        index = sel.index
        sel.annotation.set(text=info[index], position=(0, 20), anncoords="offset points")
        sel.annotation.get_bbox_patch().set(fc="white", alpha=0.6)

def get_processed_df(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in range(len(lines)):
            lines[line] = lines[line].strip()
    entries = []
    for line in lines:
        row_entry = line.split("\t")
        new_row_entry = []
        new_row_entry.append(float(row_entry[0]))
        row_entry.pop(0)
        for entry in row_entry:
            new_list = ast.literal_eval(entry)
            new_row_entry.append(new_list)
        entries.append(new_row_entry)
    print(entries[0])
    df = pd.DataFrame(columns=["rt", "mz", "approximation", "intensity", "score", "fragments", "neutrallosses"])
    for entry in entries:
        new_entry = {"rt": entry[0]}
        new_entry["mz"] = entry[1][0]
        new_entry["kendrick_mass_ch2"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((14 / 14.01565)), 6)
        new_entry["kendrick_mass_ch"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((13 / 13.007825)), 6)
        new_entry["kendrick_mass_no3"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((62 / 61.987819)), 6)
        new_entry["kendrick_mass_no2"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((46 / 45.992904)), 6)
        new_entry["kendrick_mass_so4"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((96 / 95.951732)), 6)
        new_entry["kendrick_mass_co2"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((44 / 43.989830)), 6)
        new_entry["kendrick_mass_no"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((30 / 29.997989)), 6)
        new_entry["kendrick_mass_h2o"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((18 / 18.010565)), 6)
        new_entry["kendrick_mass_oh"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((17 / 17.002740)), 6)
        new_entry["kendrick_mass_o"] = round(int(calculate_mass(entry[1][3])) - calculate_mass(entry[1][3]) * ((16 / 15.994915)), 6)
        new_entry["approximation"] = entry[1][3]
        new_entry["intensity"] = entry[1][2]
        new_entry["score"] = entry[1][1]
        new_entry["mi_ppm_deviation"] = ((calculate_mass(entry[1][3]) - entry[1][0]) / entry[1][0]) * 1000000
        available_frags = []
        frags_score = []
        frags_weight = []
        frags_intensity = []
        available_nls = []
        nls_deviation = []
        nls_weight = []
        for frag in range(len(entry)):
            if frag == 0:
                continue
            available_frags.append(entry[frag][7])
            if str(entry[frag][7]) != "None":
                available_nls.append(entry[frag][8])
                nls_deviation.append(entry[frag][9])
                frags_score.append(entry[frag][5])
                frags_weight.append(entry[frag][4])
                nls_weight.append(entry[frag][0] - entry[frag][4])
                frags_intensity.append(entry[frag][2])
            else:
                available_nls.append(None)
                nls_deviation.append(None)
                frags_score.append(None)
                frags_weight.append(None)
                nls_weight.append(None)
                frags_intensity.append(None)
        new_entry["fragments"] = [available_frags]
        new_entry["neutrallosses"] = [available_nls]
        new_entry["nls_deviation"] = [nls_deviation]
        new_entry["nls_weight"] = [nls_weight]
        new_entry["frags_score"] = [frags_score]
        new_entry["frags_weight"] = [frags_weight]
        new_entry["frags_intensity"] = [frags_intensity]
        new_df = pd.DataFrame(new_entry)
        df = pd.concat([df, new_df], ignore_index=True)
    print(len(df))
    return df


def process_fragments(df):
    # Function to check cluster widths
    def max_cluster_width(df, labels, group_col_name):
        clusters = [df[df['cluster'] == label][group_col_name] for label in np.unique(labels)]
        widths = [cluster.max() - cluster.min() for cluster in clusters]
        return max(widths)

    # Define a function to select the row to keep based on the conditions
    def select_row(cluster_df):
        all_available_fragments = cluster_df["fragments"].to_list()
        all_available_fragments_strings = []
        for i in range(len(all_available_fragments)):
            all_available_fragments_strings.extend(all_available_fragments[i])
        all_available_fragments = list(set(all_available_fragments_strings))
        keep_rows = [cluster_df["mz"].idxmax()]

        # Check for the second condition
        for i, row in cluster_df.iterrows():
            if not row["approximation"] in all_available_fragments:
                keep_rows.append(i)
            else:
                print(str(row["approximation"]) + "  /  " + str(all_available_fragments))
        keep_rows = list(set(keep_rows))
        return keep_rows


    # Initialize the number of clusters
    num_bins = 2
    max_width = 5

    # Perform clustering and check cluster widths
    while True:
        kmeans = KMeans(n_clusters=num_bins, random_state=0).fit(df[['rt']])
        df['cluster'] = kmeans.labels_
        if max_cluster_width(df, df['cluster'], "rt") <= max_width:
            break
        num_bins += 1


    # Apply the function to each cluster
    result_indices = df.groupby('cluster').apply(lambda cluster_df: select_row(cluster_df)).explode().astype(int).tolist()


    # Select the rows with these indices
    result_df = df.loc[result_indices].drop(columns='cluster')
    result_df.sort_index(inplace=True)

    df = result_df
    return df


def drop_unlogical_entries(df):
    df = df.drop(df[df.score < 5].index)
    # df = df.drop(df[df.approximation.str.contains("F")].index)
    # df = df.drop(df[df.approximation.str.contains("P")].index)
    # df = df.drop(df[df.approximation.str.contains("B")].index)
    # df = df.drop(df[df.approximation.str.contains("Cr")].index)
    # df = df.drop(df[df.approximation.str.contains("S")].index)
    mask = df["approximation"].str.contains("O", case=False, na=False)
    df = df[mask]
    mask = df["approximation"].str.contains("C", case=False, na=False)
    df = df[mask]
    df = df.drop(df[df.mi_ppm_deviation > 7].index)
    df = df.drop(df[df.mi_ppm_deviation < -7].index)
    print(df.shape)
    return df


def get_atom_ratios(df, approx_row_colname, atom_zaehler, atom_nenner):
    ratiolist = []
    for index, row in df.iterrows():
        try:
            curr_approx_dict = row[approx_row_colname]
            zaehler_atom_number = float(curr_approx_dict.get(atom_zaehler, 0))
            nenner_atom_number = float(curr_approx_dict.get(atom_nenner, 0))
            ratio = float(zaehler_atom_number) / float(nenner_atom_number)
            ratiolist.append(ratio)
        except Exception as e:
            print(e)
            ratiolist.append(None)
    return ratiolist

def get_atom_numbers(df, approx_row_colname, atomstring):
    atomnumberlist = []
    for index, row in df.iterrows():
        try:
            atomnumberlist.append(row[approx_row_colname].get(atomstring, 0))
        except Exception as e:
            print(e)
            atomnumberlist.append(0)
    return atomnumberlist
def get_subset_list(df):
    subset_list = []
    for index, row in df.iterrows():
        try:
            curr_subset = ""
            if "C" in row["approx_dict"]:
                curr_subset += "C"
            if "H" in row["approx_dict"]:
                curr_subset += "H"
            if "O" in row["approx_dict"]:
                curr_subset += "O"
            if "N" in row["approx_dict"]:
                curr_subset += "N"
            if "S" in row["approx_dict"]:
                curr_subset += "S"
            if "P" in row["approx_dict"]:
                curr_subset += "P"
            if "F" in row["approx_dict"]:
                curr_subset += "F"
            if "B" in row["approx_dict"]:
                curr_subset += "B_"
        except:
            curr_subset = "NAN"
        subset_list.append(curr_subset)
    return subset_list


def filter_dataframe(df, filter):
    fragments_filtered_df = pd.DataFrame()
    nl_filtered_df = pd.DataFrame()
    must_contain_filtered_df = pd.DataFrame()
    if "fragment=(" in filter:
        currfragment = filter.split("fragment=(")[1]
        currfragment = currfragment.split(")")[0]
        currfragment = currfragment.split(",")
        if "" in currfragment:
            print("Remove the unneccessary ',' from the filter!!!")
            currfragment = [e for e in currfragment if not e == ""]
        allowed_fragments = set(currfragment)
        fragments_filtered_df = df[df["fragments"].apply(lambda x: any(item in allowed_fragments for item in x))]
        print(fragments_filtered_df)
    if "neutralloss=(" in filter:
        currnl = filter.split("neutralloss=(")[1]
        currnl = currnl.split(")")[0]
        currnl = currnl.split(",")
        if "" in currnl:
            print("Remove the unneccessary ',' from the filter!!!")
            currnl = [e for e in currnl if not e == ""]
        allowed_nls = set(currnl)
        nl_filtered_df = df[df["neutrallosses"].apply(lambda x: any(item in allowed_nls for item in x))]
        print(nl_filtered_df)
    if "mustcontain_or=(" in filter:
        curr_must_contain = filter.split("mustcontain_or=(")[1]
        curr_must_contain = curr_must_contain.split(")")[0]
        curr_must_contain = curr_must_contain.split(",")
        if "" in curr_must_contain:
            print("Remove the unneccessary ',' from the filter!!!")
            curr_must_contain = [e for e in curr_must_contain if not e == ""]
        allowed_must_contains = set(curr_must_contain)
        must_contain_filtered_df = df[df["approximation"].apply(lambda x: any(item in allowed_must_contains for item in x))]
        print(must_contain_filtered_df)

    df = pd.concat([fragments_filtered_df, nl_filtered_df, must_contain_filtered_df])
    return df



def plot_vk_plot(df, x_vals=("O", "C"), y_vals=("H", "C"), add_filter="", approx_row_colname="approx_dict", title=""):
    #syntax >>>   "fragment=(NO3,SO4)neutralloss=(NO3)"
    if add_filter != "":
        df = filter_dataframe(df, add_filter)


    print(df)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 3], layout="tight", dpi=150)
    fig.suptitle("Van Krevelen")
    if len(x_vals) == 2:
        x_ratio = get_atom_ratios(df, approx_row_colname, x_vals[0], x_vals[1])
    elif len(x_vals) == 1:
        x_ratio = get_atom_numbers(df, approx_row_colname, x_vals[0])

    if len(y_vals) == 2:
        y_ratio = get_atom_ratios(df, approx_row_colname, y_vals[0], y_vals[1])
    elif len(y_vals) == 1:
        y_ratio = get_atom_numbers(df, approx_row_colname, y_vals[0])

    subset_list = get_subset_list(df)
    size_list = list(df["intensity"])
    point_infos = []
    for index, row in df.iterrows():
        point_infos.append(str(row["approximation"]) + "\n" + str(row["neutrallosses"]))


    y_ratio, x_ratio, size_list, colorlist, subset_list, point_info = get_cleaned_lists(y_ratio, x_ratio, size_list, subset_list, point_infos)
    scatterpoints = ax.scatter(x_ratio, y_ratio, size_list, colorlist, alpha=0.5)
    insert_interactive_scatter(scatterpoints, point_info)
    #create_shaded_areas(ax, x_ratio, y_ratio, size_list, colorlist, subset_list, proportion=0.5, number_of_groups=3)
    fig.suptitle(title)
    ax.set_xlabel("O/C")
    ax.set_ylabel("H/C")
    ax.set_xlim([0, 1.6])
    ax.set_ylim([0, 2.7])
    ax.legend()
    return fig, ax


def make_summary_of_summary(df, approx_row_colname="approx_dict"):
    total_compounds_identified = len(df)

    oc_ratio = get_atom_ratios(df, approx_row_colname, "O", "C")
    oc_ratio = [e for e in oc_ratio if not e == 0]
    entries_with_o = len(oc_ratio)
    avg_oc_ratio = sum(oc_ratio) / len(oc_ratio)

    hc_ratio = get_atom_ratios(df, approx_row_colname, "H", "C")
    hc_ratio = [e for e in hc_ratio if not e == 0]
    entries_with_h = len([e for e in hc_ratio if not e == 0])
    avg_hc_ratio = sum(hc_ratio) / len(hc_ratio)

    sc_ratio = get_atom_ratios(df, approx_row_colname, "S", "C")
    sc_ratio = [e for e in sc_ratio if not e == 0]
    entries_with_s = len([e for e in sc_ratio if not e == 0])
    avg_sc_ratio = sum(sc_ratio) / len(sc_ratio)

    nc_ratio = get_atom_ratios(df, approx_row_colname, "N", "C")
    nc_ratio = [e for e in nc_ratio if not e == 0]
    entries_with_n = len([e for e in nc_ratio if not e == 0])
    avg_nc_ratio = sum(nc_ratio) / len(nc_ratio)

    oc_nc_slope = np.polyfit(oc_ratio, nc_ratio, 1)
    oc_sc_slope = np.polyfit(oc_ratio, sc_ratio, 1)
    oc_hc_slope = np.polyfit(oc_ratio, hc_ratio, 1)

    hc_oc_slope = np.polyfit(hc_ratio, oc_ratio, 1)
    hc_sc_slope = np.polyfit(hc_ratio, sc_ratio, 1)
    hc_nc_slope = np.polyfit(hc_ratio, nc_ratio, 1)


    subset_list = get_subset_list(df)
    cho_subset_number = len([e for e in subset_list if e == "CHO"])
    chno_subset_number = len([e for e in subset_list if e == "CHON"])
    chnos_subset_number = len([e for e in subset_list if e == "CHONS"])
    chos_subset_number = len([e for e in subset_list if e == "CHOS"])
    cho_percent = cho_subset_number / total_compounds_identified
    chno_percent = chno_subset_number / total_compounds_identified
    chnos_percent = chnos_subset_number / total_compounds_identified
    chos_percent = chos_subset_number / total_compounds_identified


    ## N dataframe
    n_df = filter_dataframe(df, "mustcontain_or=(N)")
    length_n_df = len(n_df)
    percent_with_n = length_n_df / total_compounds_identified

    ## S dataframe
    s_df = filter_dataframe(df, "mustcontain_or=(S)")
    length_s_df = len(s_df)
    percent_with_s = length_s_df / total_compounds_identified

    ## O dataframe
    o_df = filter_dataframe(df, "mustcontain_or=(O)")
    length_o_df = len(o_df)
    percent_with_o = length_o_df / total_compounds_identified


    ## Nitrate compounds
    nitrate_df = filter_dataframe(df, "fragment=(N1O3)neutralloss=(N1O3)mustcontain_or=()")
    nitrate_compound_number = len(nitrate_df)

    ## Nitroaromatics
    nitroaromatics = filter_dataframe(df, "fragment=(N1O1,N1O2)neutralloss=(N1O1,N1O2)")
    nitroaromatics_compound_number = len(nitroaromatics)

    ## Sulfate/Sulfon compounds
    sulfate_df = filter_dataframe(df, "fragment=(O3S1,O4S1)neutralloss=(O3S1,O4S1)mustcontain_or=()")
    sulfate_compound_number = len(sulfate_df)

    ## Carboxy compounds
    carboxy_df = filter_dataframe(df, "fragment=(C1O2,C1H2O3,C2O4,H2O1)neutralloss=(C1O2,C1H2O3,C2O4,H2O1)mustcontain_or=()")
    carboxy_compound_number = len(carboxy_df)




if __name__ == "__main__":
    mypath = "./"
    onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    onlyfiles = [f for f in onlyfiles if f.endswith(".txt")]
    print(onlyfiles)
    dfs = []
    for file in onlyfiles:
        dfs.append(get_processed_df(file))

    for i in range(len(dfs)):
        dfs[i] = process_fragments(dfs[i])


    #molecs = [(row["rt"], row["approximation"], row["fragments"]) for index, row in dfs[0].iterrows()]
    #for i in molecs:
    #    print(i)
    #print(len([i for i in molecs if not i[2][0] == "None"]))

    for i in range(len(dfs)):
        dfs[i] = drop_unlogical_entries(dfs[i])

    for i in range(len(dfs)):
        f_dict_list = []
        for index, row in dfs[i].iterrows():
            f_dict_list.append(get_formula_to_dict(row["approximation"]))
        dfs[i]["approx_dict"] = f_dict_list

    print(dfs[0].columns)


    #reihenfolge: CHNOS
    #fig, ax = plot_vk_plot(dfs[0], add_filter="fragment=()neutralloss=(C1O2,C1H2O3,C2O4,H2O1)mustcontain_or=()", title=str(onlyfiles[0].split("/")[-1]))
    #fig.show(block=False)
    #fig2, ax2 = plot_vk_plot(dfs[1], add_filter="fragment=()neutralloss=(C1O2,C1H2O3,C2O4,H2O1)mustcontain_or=()", title=str(onlyfiles[1].split("/")[-1]))
    fig3, ax3 = plot_vk_plot(dfs[1], add_filter="fragment=(N1O3)neutralloss=(N1O3)mustcontain_or=()", title=str(onlyfiles[1].split("/")[-1]))
    fig4, ax4 = plot_vk_plot(dfs[1], add_filter="fragment=(N1O3)neutralloss=(N1O3)mustcontain_or=()", title=str(onlyfiles[1].split("/")[-1]))
    plt.show()

    for file in onlyfiles:
        foldername = str(file.split("/")[-1].split(".")[0] + "/")
        os.makedirs(foldername, exist_ok=True)



