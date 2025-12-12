# plotting.py
import math
import os
import PIL
import matplotlib
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import traceback
import matplotlib.figure

import UVenture.formula_calculations as formula_calculations
import UVenture.MS_functions as MS_functions



DPI = 300
matplotlib.use("Agg")  # Use a non-interactive backend for matplotlib
print(matplotlib.__file__)   # should point into site-packages
print(matplotlib.__version__)



def create_2d_massspec_plot(ms_file_object, filter_mode="Full scan", save="", max_dim=5000):
    """
    Creates a 2D mass spectrum plot using matplotlib.
    """
    Y = [ms_file_object.rawdata[i]["m/z array"] for i in range(len(ms_file_object.rawdata)) if (ms_file_object.all_modes[i] == filter_mode)]
    intensity = [ms_file_object.rawdata[i]["intensity array"] for i in range(len(ms_file_object.rawdata)) if (ms_file_object.all_modes[i] == filter_mode)]
    # Flatten data
    X_vals, Y_vals, I_vals = [], [], []
    for x, (y_row, i_row) in enumerate(zip(Y, intensity)):
        for y, i in zip(y_row, i_row):
            X_vals.append(x)
            Y_vals.append(y)
            I_vals.append(i)
    X_vals = np.array(X_vals, dtype=float)
    Y_vals = np.array(Y_vals, dtype=float)
    I_vals = np.array(I_vals, dtype=float)
    # Determine data bounds
    x_min, x_max = X_vals.min(), X_vals.max()
    y_min, y_max = Y_vals.min(), Y_vals.max()

    # Target image resolution
    max_dim = max_dim
    # Compute scaling factors
    x_range = x_max - x_min
    y_range = y_max - y_min
    # Determine scale to fit date within 5000x5000 while preserving aspect ratio
    if len(Y) < max_dim:
        max_dim_x = len(Y)-1
    else:
        max_dim_x = max_dim
    scalex = max(x_range / max_dim_x, 1e-9)
    scaley = max(y_range / max_dim, 1e-9)
    width = int(np.ceil(x_range / scalex)) + 1
    height = int(np.ceil(y_range / scaley)) + 1
    # Shift and scale coordinates into image space
    X_img = ((X_vals - x_min) / scalex).astype(int)
    Y_img = ((y_max - Y_vals) / scaley).astype(int)  # Invert Y

    # Step 1: Create a 2D array (single channel) for grayscale image
    gray_image = np.zeros((height, width), dtype=np.float32)  # or np.uint16 if preferred
    # Sum intensities per pixel
    #from collections import defaultdict
    #pixel_intensity = defaultdict(float)
    pixel_intensity = {}
    for x, y, i in zip(X_img, Y_img, I_vals):
        try:
            pixel_intensity[(y, x)] += i
        except:
            pixel_intensity[(y, x)] = i
    # Map intensity to colormap (viridis)
    for (y, x), val in pixel_intensity.items():
        gray_image[y, x] = val
    if not save == "":
        foldername = os.path.normpath(save)
        foldername = os.path.dirname(foldername)
        os.makedirs(foldername, exist_ok=True)
        myimage = PIL.Image.fromarray(gray_image, mode='I;16')
        myimage.save(save)
    
    return True


def create_barchart_massspec_with_go(masses, intensities, plot_filepath):
    # Sort the data by intensities
    data = sorted(zip(masses, intensities), key=lambda x: x[1], reverse=True)
    # Select the top 10
    top_10_data = data[:10]
    # Create the figure
    go_fig = go.Figure(data=go.Bar(x=masses, y=intensities, marker=dict(color='black', opacity=1), width=0.1))
    # Add data labels for the top 10 values
    for mass, intensity in top_10_data:
        go_fig.add_annotation(x=mass, y=intensity, text=str(mass), showarrow=False, font=dict(size=12, color="Black"), bgcolor="White", opacity=0.8, textangle=-90)
    go_fig.update_layout(plot_bgcolor='white')
    plot_filepath = os.path.normpath(plot_filepath)
    os.makedirs(os.path.dirname(plot_filepath), exist_ok=True)
    go_fig.write_html(plot_filepath)


def create_barchart_massspec(masses, intensities, title, image_filepath):
    fig = matplotlib.figure.Figure(layout="tight")
    ax = fig.subplots()
    ax.bar(masses, intensities, label="measured ions", width=0.2, color="blue")
    # Sort the data by intensities
    data = sorted(zip(masses, intensities), key=lambda x: x[1], reverse=True)
    # Select the top 4
    top_10_data = data[:4]
    # Add labels for the top 4 values
    for mass, intensity in top_10_data:
        ax.text(mass, (intensity), str(round(mass, 4)), ha='center', va='bottom', rotation=0)
    try:
        ax.get_legend().remove()
    except:
        pass
    ax.set_xlabel("masses / Da")
    ax.set_ylabel("intensity / a.u.")

    ax.set_title(title)
    matplotlib.rcParams.update({'figure.autolayout': True})
    image_filepath = os.path.normpath(image_filepath)
    os.makedirs(os.path.dirname(image_filepath), exist_ok=True)
    fig.savefig(image_filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()


def create_isotopo_plot(spec_masses, spec_intensities, formula_to_simulate, filepath, mass_deviation=10):
    print("Making plot of isotopologues for formula: " + str(formula_to_simulate))
    isotope_simulation = formula_calculations.simulate_isotope_pattern_of_formula(formula_to_simulate, debug_output=False)
    print("Isotope simulation: " + str(isotope_simulation))
    
    intensities_to_include = []
    for entry in list(isotope_simulation.keys()):
        mass_with_min_deviation = min(spec_masses, key=lambda x: abs(entry - x))
        if abs((mass_with_min_deviation - entry) / entry) * 1000000 < mass_deviation:
            summed_intensity = spec_intensities[spec_masses.index(mass_with_min_deviation)]
            intensities_to_include.append(summed_intensity)
        else:
            break
    if len(intensities_to_include) == 0:
        print("ERROR: NO intensities to include in isotopo plot! see plotting.py")
        intensities_to_include = [1]
    sum_intensities = sum(intensities_to_include)
    intensities_to_include = [(i/sum_intensities) for i in intensities_to_include] # normalize intensities
    
    fig = matplotlib.figure.Figure()
    gs = fig.add_gridspec(1, len(intensities_to_include), hspace=0, wspace=0)
    ax = gs.subplots(sharex="col", sharey="row")

    for entry in range(len(intensities_to_include)):
        cmap = matplotlib.colormaps.get_cmap('Greens')
        colors = cmap((spec_intensities - min(spec_intensities)) / (max(spec_intensities) - min(spec_intensities)) * 0.3 + 0.7)                
        ax[entry].bar(spec_masses, spec_intensities, color=colors, label="all ions", width=0.0005, alpha=1)
        ax[entry].bar(list(isotope_simulation.keys())[entry], intensities_to_include[entry], label="measured ions", color="blue", width=0.005, alpha=0.5)
        ax[entry].bar(list(isotope_simulation.keys())[entry], -1 * list(isotope_simulation.values())[entry], label="simulated intensity", color="red", width=0.005, alpha=0.5)
        ax[entry].axhline(0, color='black', linewidth=1)
        ax[entry].set_xlim(list(isotope_simulation.keys())[entry] - 0.1, list(isotope_simulation.keys())[entry] + 0.1)
        ax[entry].set_ylim(-1, 1)
        if entry == 0:  # only set the y-label for the first (leftmost) subplot
            ax[entry].set_ylabel("intensity / a.u.")
        if entry == len(intensities_to_include) // 2:  # only set the x-label for the middle subplot
            ax[entry].set_xlabel("masses / Da")
        ax[entry].set_title(str(round(list(isotope_simulation.keys())[entry], 4)), rotation='vertical')
    for a in fig.get_axes():
        a.label_outer()
    #set the title of the whole figure so that it will be displayed above the subplots
    formula_string = formula_calculations.get_formula_string_from_dict(formula_to_simulate)
    fig.suptitle("Isotopologues plot for formula: " + str(formula_string))
    for a in ax:
        a.legend()
    matplotlib.rcParams.update({'figure.autolayout': True})
    filepath = os.path.normpath(filepath)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    fig.savefig(filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()
    return True


def create_isotopo_plot_with_go(spec_masses, spec_intensities, formula_to_simulate, filepath, mass_deviation=10):
    isotope_simulation = formula_calculations.simulate_isotope_pattern_of_formula(formula_to_simulate, debug_output=False)
    fig = make_subplots(rows=1, cols=len(isotope_simulation), shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0, vertical_spacing=0)
    intensities_to_include = []
    for entry in list(isotope_simulation.keys()):
        summed_intensity = sum([spec_intensities[i] for i in range(len(spec_masses)) if abs(((spec_masses[i] - entry) / entry)*1000000) <= (mass_deviation)])
        intensities_to_include.append(summed_intensity)
    intensities_to_include = [(i/sum(intensities_to_include)) for i in intensities_to_include] # normalize intensities

    for entry in range(len(intensities_to_include)):
        curr_x_values = [spec_masses[i] for i in range(len(spec_masses)) if list(isotope_simulation.keys())[entry]-0.1 <= spec_masses[i] <= list(isotope_simulation.keys())[entry]+0.1] 
        curr_y_values = [spec_intensities[i] for i in range(len(spec_masses)) if list(isotope_simulation.keys())[entry]-0.1 <= spec_masses[i] <= list(isotope_simulation.keys())[entry]+0.1]  
        curr_y_values = [(i/sum(curr_y_values)) for i in curr_y_values]
        fig.add_trace(go.Bar(x=curr_x_values, y=curr_y_values, name="all ions", marker=dict(color='green', opacity=1), width=0.001), row=1, col=entry+1)
        fig.add_trace(go.Bar(x=[list(isotope_simulation.keys())[entry]], y=[intensities_to_include[entry]], name="measured ions", marker=dict(color='blue', opacity=0.65), width=0.01), row=1, col=entry+1)
        fig.add_trace(go.Bar(x=[list(isotope_simulation.keys())[entry]], y=[-1 * list(isotope_simulation.values())[entry]], name="simulated intensity", marker=dict(color='red', opacity=0.65), width=0.01), row=1, col=entry+1)
        fig.update_xaxes(range=[(list(isotope_simulation.keys())[entry]) - 0.1, (list(isotope_simulation.keys())[entry]) + 0.1], row=1, col=entry+1)
        fig.add_shape(type="line", x0=0, x1=400, y0=-0.01, y1=0.01, line=dict(color="black", width=10), row=1, col=entry+1)
    fig.update_yaxes(range=[-1, 1], showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor="Gray", )
    fig.update_layout(shapes=[dict(type="rect", xref="paper", yref="paper", x0=0, y0=0, x1=1, y1=1, line=dict(color="Black", width=4))],
                        barmode='overlay', title_text="Isotopologues Plot", xaxis_title="masses / Da", yaxis_title="intensity / a.u.", plot_bgcolor='white')
    filepath = os.path.normpath(filepath)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    fig.write_html(filepath)
    return True


def create_xic(rt_list, intensities, title, filepath, retention_time=0, width_observed_time_seconds=-1):
    """
    Creates an Extracted Ion Chromatogram (XIC) plot using matplotlib.
    width_observed_time_seconds can be either of the three cases:
        -1 for automatic detection of the red bar in the xic, 
        0 for no red bar in the XIC,
        arb. pos. number for a red bar in the xic of the specified length.

    """
    fig = matplotlib.figure.Figure()
    ax = fig.subplots()
    identified_peaks = MS_functions.get_peaks_in_xy_series(rt_list, intensities, sg_window=10, sg_order=3)
    identified_peak_times = [entry[0] for entry in identified_peaks]
    plot_heights = [max(intensities)/4 for i in range(len(identified_peak_times))]
    ax.scatter(identified_peak_times, plot_heights, color="green", label="identified peak", s=20, alpha=0.5)
    if not retention_time == 0:
        if width_observed_time_seconds == -1:
            width_observed_time_seconds = 10
            width_observed_time_seconds = (max(rt_list) - min(rt_list)) / 30
        if not width_observed_time_seconds == 0:
            ax.bar(retention_time, max(intensities), color="red", label="observed time", alpha=0.5, width=width_observed_time_seconds)
        xic_peak_index = rt_list.index(min(rt_list, key=lambda x: abs(retention_time - x)))
        ax.scatter(retention_time, intensities[xic_peak_index], color="red", marker="x", s=100)
    ax.plot(rt_list, intensities, color="blue", label="XIC")
    ax.set_xlabel("retention time / s")
    ax.set_ylabel("intensity / a.u.")
    ax.set_title(title, wrap=True)
    ax.legend()
    filepath = os.path.normpath(filepath)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    fig.savefig(filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()
    return True


def create_oa_summary_plot(oa_summary_object, true_fragment_list, best_frag_spec, fragment_predictions, save_filepath):
    """
    Creates a summary plot for the OA (one analysis) object.
    """
    ncols = 3
    nrows = 1
    frag_masses = [frag[4] for frag in true_fragment_list if frag[4] >= 1]
    nrows = nrows + int( (len(frag_masses)+1)/3 + 0.667 )
    fig = matplotlib.figure.Figure(constrained_layout=True, figsize=(5*ncols, 5*nrows), dpi=DPI)
    gs = matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    ax = [[None, None, None]]
    ax[0][0] = fig.add_subplot(gs[0, 0])
    ax[0][1] = fig.add_subplot(gs[0, 1])
    ax[0][2] = fig.add_subplot(gs[0, 2])
    tic_ax = fig.add_subplot(gs[1,0])
    frag_axs = []
    for i in range(1, len(frag_masses)+1, 1):
        row_index = int(i/3 + 1)
        col_index = int(i%3)
        frag_axs.append(fig.add_subplot(gs[row_index, col_index]))

    fig.suptitle("Summary Plot for Mass: " + str(round(oa_summary_object.mass, 4)) + " at RT: " + str(round(oa_summary_object.rt, 2)) + "\nPredicted Formula: " + str(oa_summary_object.best_molecular_ion_prediction), wrap=True, fontsize=16, fontweight="bold")

    ax[0][0].bar(best_frag_spec.summarized_masses, best_frag_spec.summarized_intensities, color="gray", label="Spectrum", alpha=0.2)
    intensity_of_molecular_ion_in_frag_spec = sum([best_frag_spec.summarized_intensities[i] for i in range(len(best_frag_spec.summarized_intensities)) if abs(((best_frag_spec.summarized_masses[i] - oa_summary_object.mass)/oa_summary_object.mass)*1000000) <= oa_summary_object.kwargs["mass_deviation"]])
    ax[0][0].bar(oa_summary_object.mass, intensity_of_molecular_ion_in_frag_spec, color="red", width=1.5, label="Molecule")
    ax[0][0].text(oa_summary_object.mass, (intensity_of_molecular_ion_in_frag_spec), str(true_fragment_list[0][3]), ha='center', va='bottom', rotation=0)
    ax[0][0].text(0.01, 0.99, "Max int in spec.: " + str(max(best_frag_spec.summarized_intensities)), transform=ax[0][0].transAxes, ha="left", va="top", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    for fragment in true_fragment_list:
        ax[0][0].bar(fragment[4], fragment[6], color="blue", label=str(fragment[7]))
        ax[0][0].text(fragment[4], (fragment[6]), str(fragment[7]), ha='center', va='bottom', rotation=0)
    ax[0][0].set_title("Fragment spectrum\n" + " RT: " + str(best_frag_spec.rt), wrap=True)
    ax[0][0].set_xlabel("ion mass / u")
    ax[0][0].set_ylabel("intensity / a.u.")
    max_y_scale = max([frag[6] for frag in true_fragment_list])
    max_y_scale = max([max_y_scale, intensity_of_molecular_ion_in_frag_spec])
    ax[0][0].set_ylim([0, max_y_scale*1.3])
    ax[0][0].legend()

    ax[0][1].bar(oa_summary_object.best_molecular_ion_spec.summarized_masses, oa_summary_object.best_molecular_ion_spec.summarized_intensities, color="gray", label="Spectrum", alpha=0.2)
    intensity_of_molecular_ion_in_mi_spec = sum([oa_summary_object.best_molecular_ion_spec.summarized_intensities[i] for i in range(len(oa_summary_object.best_molecular_ion_spec.summarized_intensities)) if abs(((oa_summary_object.best_molecular_ion_spec.summarized_masses[i] - oa_summary_object.mass)/oa_summary_object.mass)*1000000) <= oa_summary_object.kwargs["mass_deviation"]])
    ax[0][1].bar(oa_summary_object.mass, intensity_of_molecular_ion_in_mi_spec, color="red", width=1.5, label="Molecule")
    ax[0][1].text(oa_summary_object.mass, (intensity_of_molecular_ion_in_mi_spec), str(true_fragment_list[0][3]), ha='center', va='bottom', rotation=0)
    ax[0][1].text(0.01, 0.99, "Max int in spec.: " + str(max(oa_summary_object.best_molecular_ion_spec.summarized_intensities)), transform=ax[0][1].transAxes, ha="left", va="top", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    frags_in_mi_spec = {}
    for fragment in true_fragment_list:
        frag_mass = fragment[4]
        frag_formula = fragment[7]
        within_mi_spec = any(math.isclose(frag_mass, v, rel_tol=0.0, abs_tol=( (oa_summary_object.kwargs["mass_deviation"] * frag_mass) / 1000000 )) for v in oa_summary_object.best_molecular_ion_spec.summarized_masses)
        frag_mass_in_mi_spec = min(oa_summary_object.best_molecular_ion_spec.summarized_masses, key=lambda v: abs(v - frag_mass))
        frag_int_in_mi_spec = oa_summary_object.best_molecular_ion_spec.summarized_mass_intensity_dict[frag_mass_in_mi_spec]
        if within_mi_spec:
            frags_in_mi_spec[frag_mass_in_mi_spec] = frag_int_in_mi_spec
            ax[0][1].bar(frag_mass_in_mi_spec, frag_int_in_mi_spec, color="blue", label=frag_formula)
            ax[0][1].text(frag_mass_in_mi_spec, (frag_int_in_mi_spec), str(fragment[7]), ha='center', va='bottom', rotation=0)
    ax[0][1].set_title("Molecular ion spectrum\n" + " RT: " + str(oa_summary_object.best_molecular_ion_spec.rt), wrap=True)
    ax[0][1].set_xlabel("ion mass / u")
    ax[0][1].set_ylabel("intensity / a.u.")
    try:
        max_y_scale = max([v for k, v in frags_in_mi_spec.items()])
    except ValueError:
        max_y_scale = 0
    max_y_scale = max([max_y_scale, intensity_of_molecular_ion_in_mi_spec])
    ax[0][1].set_ylim([0, max_y_scale*1.3])
    highest_signals = top5_items = sorted(oa_summary_object.best_molecular_ion_spec.summarized_mass_intensity_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[:5]
    mi_text = "Highest signals:\n"
    for ma, i in highest_signals:
        mi_text = mi_text + "-   " + str(round(ma, 4)) + " \n"
    ax[0][1].text(0.01, 0.93, mi_text, transform=ax[0][1].transAxes, ha="left", va="top", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    ax[0][1].legend()

    xic = oa_summary_object.xic
    peak_index = oa_summary_object.peak_index    
    identified_peaks = MS_functions.get_peaks_in_xy_series(xic[0], xic[1], sg_window=10, sg_order=3)
    identified_peak_times = [entry[0] for entry in identified_peaks]
    plot_heights = [max(xic[1])/4 for i in range(len(identified_peak_times))]
    ax[0][2].scatter(identified_peak_times, plot_heights, color="green", label="identified peak", s=20, alpha=0.5)
    width_observed_time_seconds = -1
    if not oa_summary_object.rt == 0:
        if width_observed_time_seconds == -1:
            width_observed_time_seconds = 10
            width_observed_time_seconds = (max(xic[0]) - min(xic[0])) / 30
        if not width_observed_time_seconds == 0:
            ax[0][2].bar(oa_summary_object.rt, max(xic[1]), color="red", label="observed time", alpha=0.5, width=width_observed_time_seconds)
        xic_peak_index = xic[0].index(min(xic[0], key=lambda x: abs(oa_summary_object.rt - x)))
        ax[0][2].scatter(oa_summary_object.rt, xic[1][xic_peak_index], color="red", marker="x", s=100)
    ax[0][2].plot(xic[0], xic[1], color="blue", label="XIC")
    ax[0][2].set_xlabel("retention time / s")
    ax[0][2].set_ylabel("intensity / a.u.")
    title = "XIC m=" + str(round(oa_summary_object.mass, 4)) + " +- " + str( round(((oa_summary_object.kwargs["mass_deviation"]*oa_summary_object.mass) / 1000000), 4) ) + " [" + str( oa_summary_object.kwargs["oa_xic_requested_filter_mode"] ) + "]"
    ax[0][2].set_title(title, wrap=True)
    ax[0][2].legend()


    xic = oa_summary_object.xic
    tic = oa_summary_object.ms_file.tic
    rts = oa_summary_object.ms_file.rt_list
    width_observed_time_seconds = -1
    if not oa_summary_object.rt == 0:
        if width_observed_time_seconds == -1:
            width_observed_time_seconds = 10
            width_observed_time_seconds = (max(xic[0]) - min(xic[0])) / 30
        if not width_observed_time_seconds == 0:
            tic_ax.bar(oa_summary_object.rt, max(xic[1]), color="red", label="observed time", alpha=0.5, width=width_observed_time_seconds)
        xic_peak_index = xic[0].index(min(xic[0], key=lambda x: abs(oa_summary_object.rt - x)))
        tic_ax.scatter(oa_summary_object.rt, xic[1][xic_peak_index], color="red", marker="x", s=100)
    tic_ax.plot(rts, tic, color="black", label="TIC")
    tic_ax.plot(xic[0], xic[1], color="blue", label="Molecular ion XIC", alpha=0.7)
    tic_ax.set_xlabel("retention time / s")
    tic_ax.set_ylabel("intensity / a.u.")
    title = "TIC and XIC of molecular ion (" + str(round(oa_summary_object.mass, 4)) + ")"
    tic_ax.set_title(title, wrap=True)
    tic_ax.legend()


    for i, ax in enumerate(frag_axs):
        frag_score = true_fragment_list[i][5]
        nl = true_fragment_list[i][8]
        fragformula = true_fragment_list[i][7]
        frag_int = true_fragment_list[i][6]
        curr_xic = oa_summary_object.ms_file.get_xic(frag_masses[i], round(((oa_summary_object.kwargs["mass_deviation"]*oa_summary_object.mass) / 1000000), 4), requested_filter_mode=oa_summary_object.best_frag_spec.filter_mode)
        width_observed_time_seconds = (max(curr_xic[0]) - min(curr_xic[0])) / 30
        ax.bar(oa_summary_object.ms_file.rt_list[peak_index], max(curr_xic[1]), color="red", label="observed time", alpha=0.5, width=width_observed_time_seconds)
        ax.plot(curr_xic[0], curr_xic[1], color="blue", label="XIC " + str(round(frag_masses[i], 4)))
        ax.set_xlabel("retention time / s")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("XIC for fragment: " + str(round(frag_masses[i], 4)) + " at RT: " + str(round(oa_summary_object.rt, 2)) + " \nFormula: " + str(fragformula), wrap=True)
        frag_text = str(fragformula) + "\n" + \
                    "Score: " + str(round(frag_score, 1)) + "\n" + \
                    "Intensity: " + str(round(frag_int, 1)) + "\n" + \
                    "NL: " + str(nl)
        ax.text(0.01, 0.99, frag_text, transform=ax.transAxes, ha="left", va="top", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
        ax.legend()

    matplotlib.rcParams.update({'figure.autolayout': True})
    save_filepath = os.path.normpath(save_filepath)
    os.makedirs(os.path.dirname(save_filepath), exist_ok=True)
    fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()


def create_xic_matching_plot(peak1_rt, peakintensity1, peak2_rt, peakintensity2, title, area_between_curves, save_filepath):
    fig = matplotlib.figure.Figure()
    ax = fig.subplots()
    ax.plot(peak1_rt, peakintensity1, color="blue", label="Normalized molecular ion peak")
    ax.plot(peak2_rt, peakintensity2, color="red", label="Normalized fragment peak")
    ax.text(min(peak1_rt), max(peakintensity1), "Area between curves: " + str(round(area_between_curves, 2)))
    ax.set_title(title, wrap=True)
    ax.set_xlabel("retention time / seconds")
    ax.set_ylabel("normalized intensity / a.u.")
    ax.legend()
    matplotlib.rcParams.update({'figure.autolayout': True})
    save_filepath = os.path.normpath(save_filepath)
    os.makedirs(os.path.dirname(save_filepath), exist_ok=True)
    fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()





