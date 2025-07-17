# plotting.py
import os
import PIL
import matplotlib
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import traceback

import UVenture.formula_calculations as formula_calculations
import UVenture.MS_functions as MS_functions



DPI = 300
matplotlib.use("Agg")  # Use a non-interactive backend for matplotlib




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
    fig.write_html(filepath)
    return True


def create_xic(rt_list, intensities, title, filepath, retention_time=0):
    """
    Creates an Extracted Ion Chromatogram (XIC) plot using matplotlib.
    """
    fig = matplotlib.figure.Figure()
    ax = fig.subplots()
    identified_peaks = MS_functions.get_peaks_in_xy_series(rt_list, intensities, sg_window=10, sg_order=3)
    identified_peak_times = [entry[0] for entry in identified_peaks]
    plot_heights = [max(intensities)/4 for i in range(len(identified_peak_times))]
    ax.scatter(identified_peak_times, plot_heights, color="green", label="identified peak", s=20, alpha=0.5)
    if not retention_time == 0:
        width_observed_time = (max(rt_list) / len(rt_list)) * 40
        ax.bar(retention_time, max(intensities), color="red", label="observed time", alpha=0.5, width=width_observed_time)
        xic_peak_index = rt_list.index(min(rt_list, key=lambda x: abs(retention_time - x)))
        ax.scatter(retention_time, intensities[xic_peak_index], color="red", marker="x", s=100)
    ax.plot(rt_list, intensities, color="blue", label="XIC")
    ax.set_xlabel("retention time / s")
    ax.set_ylabel("intensity / a.u.")
    ax.set_title(title, wrap=True)
    ax.legend()
    fig.savefig(filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()
    return True


def create_oa_summary_plot(oa_summary_object, true_fragment_list, best_frag_spec, fragment_predictions, save_filepath):
    """
    Creates a summary plot for the OA (one analysis) object.
    """
    ncols = 2
    nrows = 2
    frag_masses = [frag[4] for frag in true_fragment_list if frag[4] >= 1]
    nrows = nrows + int( len(frag_masses)/2 + 0.5 )
    fig = matplotlib.figure.Figure(constrained_layout=True, figsize=(5*ncols, 5*nrows), dpi=DPI)
    gs = matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    ax = [[None, None],
            [None      ]]
    ax[0][0] = fig.add_subplot(gs[0, 0])
    ax[0][1] = fig.add_subplot(gs[0, 1])
    ax[1][0] = fig.add_subplot(gs[1, :])
    ax[1][0].axis("off")
    frag_axs = []
    for i in range(0, len(frag_masses), 1):
        row_index = int(i/2 + 2)
        col_index = int(i%2)
        frag_axs.append(fig.add_subplot(gs[row_index, col_index]))

    fig.suptitle("Summary Plot for Mass: " + str(round(oa_summary_object.mass, 4)) + " at RT: " + str(round(oa_summary_object.rt, 2)), wrap=True, fontsize=16, fontweight="bold")

    ax[0][0].bar(best_frag_spec.summarized_masses, best_frag_spec.summarized_intensities, color="gray", label="Spectrum", alpha=0.2)
    intensity_of_molecular_ion_in_frag_spec = sum([best_frag_spec.summarized_intensities[i] for i in range(len(best_frag_spec.summarized_intensities)) if abs(((best_frag_spec.summarized_masses[i] - oa_summary_object.mass)/oa_summary_object.mass)*1000000) <= oa_summary_object.kwargs["mass_deviation"]])
    ax[0][0].bar(oa_summary_object.mass, intensity_of_molecular_ion_in_frag_spec, color="red", width=1.5, label="Molecule")
    ax[0][0].text(oa_summary_object.mass, (intensity_of_molecular_ion_in_frag_spec), str(true_fragment_list[0][3]), ha='center', va='bottom', rotation=0)
    for fragment in true_fragment_list:
        ax[0][0].bar(fragment[4], fragment[6], color="blue", label="Fragment")
        ax[0][0].text(fragment[4], (fragment[6]), str(fragment[7]), ha='center', va='bottom', rotation=0)
    ax[0][0].set_title("Fragment spectrum\n" + str(best_frag_spec.filter) + "\nIndex: " + str(best_frag_spec.index) + " RT: " + str(best_frag_spec.rt), wrap=True)
    ax[0][0].set_xlabel("ion mass / u")
    ax[0][0].set_ylabel("intensity / a.u.")
    ax[0][0].legend()

    text_str = ""
    text_str = text_str + "{:<12}|{:<15}|{:<12}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format("Comb. score", "Frag mass", "Intensity", "Frag formula", "NL formula", "Molec formula", "score F", "dev. NL", "score M")
    text_str = text_str + "{:_<12}|{:_<15}|{:_<12}|{:_<15}_{:_<15}___{:_<15}|{:_>10}|{:_>10}|{:_>10}\n".format("", "", "", "", "", "", "", "", "")
    for fragment in oa_summary_object.matching_fragments_list:
        comb_mi_score = sum([frag[1] for frag in oa_summary_object.matching_fragments_list if frag[3] == fragment[3]])
        text_str = text_str + "{:<12}|{:<15}|{:<12}|{:<15} {:<15}=> {:<15}|{:>10}|{:>10}|{:>10}\n".format(str(round(comb_mi_score, 1)), 
                                                                                                            str(round(fragment[4], 5)),
                                                                                                        str(round(fragment[6], 1)),
                                                                                                        str(fragment[7]),
                                                                                                        str(fragment[8]),
                                                                                                        str(fragment[3]),
                                                                                                        str(round(fragment[5], 1)),
                                                                                                        str(round(fragment[9], 3)),
                                                                                                        str(round(fragment[1], 1)))
    text_str = text_str + "================================================================================\n" + "All fragment predictions for spec: \n"

    text_str_all = "{:<10}|{:<15}|{:<12}\n".format("Mass", "Formula", "Score")
    text_str_all = text_str_all + "{:_<10}|{:_<15}|{:_<12}\n".format("", "", "")
    for pred_mass, pred in list(fragment_predictions.items()):
        try:
            try:
                pred_formula = "".join([str(a) + str(n) for a, n in pred.best_formula_prediction.items()])
            except AttributeError:
                pred_formula = "None"
            pred_score = pred.score_of_best_formula
            if pred_score is None:
                pred_score = 0
            text_str_all = text_str_all + "{:<10}|{:<15}|{:<12}\n".format(str(round(pred_mass, 4)), str(pred_formula), str(round(pred_score, 1)))
        except Exception as e:
            print("Error creating text_str_all for summary plot.")
            print(e)
            print(traceback.format_exc())
    
    text_str = text_str + text_str_all

    props = dict(boxstyle='round', facecolor='grey', alpha=0.05)  # bbox features
    text = ax[1][0].text(0.02, 0.98, text_str, fontfamily="monospace", transform=ax[1][0].transAxes, fontsize=8,
                verticalalignment="top", bbox=props)

    xic = oa_summary_object.xic
    peak_index = oa_summary_object.peak_index
    xic_peak_index = xic[2].index(min(xic[2], key=lambda x: abs(peak_index - x)))
    ax[0][1].scatter(oa_summary_object.ms_file.rt_list[peak_index], xic[1][xic_peak_index], color="red", marker="x", s=100)
    ax[0][1].plot(xic[0], xic[1], label="XIC measured")
    ax[0][1].set_title("xic_" + str(round(oa_summary_object.mass, 4)) + "+-" + str( round(((oa_summary_object.kwargs["mass_deviation"]*oa_summary_object.mass) / 1000000), 4) ) + "_" + str(oa_summary_object.kwargs["oa_xic_requested_filter_mode"]), wrap=True)
    ax[0][1].set_xlabel("retention time / seconds")
    ax[0][1].set_ylabel("intensity / a.u.")

    for i, ax in enumerate(frag_axs):
        curr_xic = MS_functions.get_xic(oa_summary_object.ms_file.rawdata, frag_masses[i], round(((oa_summary_object.kwargs["mass_deviation"]*oa_summary_object.mass) / 1000000), 4), requested_filter_mode=oa_summary_object.best_frag_spec.filter_mode)
        observed_window_width = (max(oa_summary_object.ms_file.rt_list) / len(oa_summary_object.ms_file.rt_list)) * 20
        ax.bar(oa_summary_object.ms_file.rt_list[peak_index], max(curr_xic[1]), color="red", label="observed time", alpha=0.5, width=observed_window_width)
        ax.plot(curr_xic[0], curr_xic[1], color="blue", label="XIC " + str(round(frag_masses[i], 4)))
        ax.set_xlabel("retention time / s")
        ax.set_ylabel("intensity / a.u.")
        ax.set_title("Extracted Ion Chromatogram for fragment: " + str(round(frag_masses[i], 4)) + " u at RT: " + str(round(oa_summary_object.rt, 2)) + " s with mode: " + str(oa_summary_object.best_frag_spec.filter_mode), wrap=True)
        ax.legend()

    matplotlib.rcParams.update({'figure.autolayout': True})
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
    fig.savefig(save_filepath, bbox_inches='tight', dpi=DPI)
    fig.clf()
    fig.clear()





