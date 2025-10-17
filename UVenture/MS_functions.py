# MS_functions.py
import copy
import math
import sys
import traceback
import similaritymeasures
import scipy
import numpy as np
from pyteomics import mzml
import pandas as pd
import UVenture.GaussianCompressor as GaussianCompressor






def _rt_window_indices_scan(rts, center, half_width=50):
    # Returns a list of indices that fall within this window.
    # E.g.: rts = [2,4,5,6,7,9,10,11,12,13], center = 7, half_width=2
    # OUTPUT: [2,3,4,5] <-- those are the indices of the retention times that fall within the specified window.
    lo_v, hi_v = center - half_width, center + half_width
    out = []
    for i, v in enumerate(rts):
        if v < lo_v:
            continue
        if v > hi_v:
            break
        out.append(i)
    return out



def read_mzml_file(mzml_filename):
    f = mzml.read(mzml_filename)
    f = list(f)
    return f


def get_xic(f, mass, mass_deviation, requested_filter_mode="Full scan"):
    # The mass deviation is defined as the requested mass +1x the mass deviation and -1x the mass deviation.
    # If the requested mass is 1000 and the mass deviation is 0.005, the range is 999.995 to 1000.005.
    if (not requested_filter_mode == "Full scan") and (not requested_filter_mode == "AIF") and (not requested_filter_mode == "MS/MS"):
        requested_filter_mode = "Full scan"
    rt_list = []
    for element in f:
        rt_list.append(element["scanList"]["scan"][0]["scan time"])
    intensity_list = []
    index_with_wrong_ms_level_list = []
    for entry in f:
        filter = entry["scanList"]["scan"][0]["filter string"]
        filter_mode = ""
        if (" d " in filter) and ("@hcd" in filter):
            filter_mode = "MS/MS"
        if (not " d " in filter) and (not "hcd" in filter):
            filter_mode = "Full scan"
        if (not " d " in filter) and ("hcd" in filter):
            filter_mode = "AIF"
        if not filter_mode == requested_filter_mode:
            index_with_wrong_ms_level_list.append(f.index(entry))
        mass_indices = (i for i in range(len((entry["m/z array"]))) if mass-mass_deviation <= entry["m/z array"][i] <= mass+mass_deviation)
        try:
            curr_sum_int = 0
            for mass_index in mass_indices:
                curr_sum_int = curr_sum_int + entry["intensity array"][mass_index]
            intensity_list.append(curr_sum_int)
        except:
            intensity_list.append(0)
    original_index_list = list(range(len(f)))
    if not requested_filter_mode == "all":
        for index in sorted(index_with_wrong_ms_level_list, reverse=True):
            del rt_list[index]
            del intensity_list[index]
            del original_index_list[index]
    return [rt_list, intensity_list, original_index_list]


def _odd_below(n):
    return n-1 if n % 2 == 0 else n

def _gaussian_overlap(mu1, s1, mu2, s2):
    # Closed-form overlap of two normalized Gaussians (0..1)
    num = math.sqrt(2*s1*s2/(s1*s1 + s2*s2))
    den = 4*(s1*s1 + s2*s2)
    return num * math.exp(-((mu1 - mu2)**2)/den)
   
def compare_peak_shape_similarity_Gauss(xic1_original, xic2_original, peak_rt, peakwidth=10, debug_output=False):
    import copy
    import numpy as np
    from scipy.signal import savgol_filter, find_peaks, peak_widths
    from scipy.optimize import curve_fit

    # ----------------- helpers -----------------
    def _odd(n): 
        n = int(n)
        return n if n % 2 == 1 else max(1, n-1)

    def _smooth_nonneg(y, win, order=2):
        y = np.asarray(y, float)
        if len(y) >= 5 and win >= 3:
            win = min(_odd(len(y)-1), max(_odd(order+3), _odd(win)))
            ys = savgol_filter(y, win, order, mode="nearest")
        else:
            ys = y.copy()
        ys[ys < 0] = 0.0
        return ys

    def _idx_near(t, target):
        return int(np.argmin(np.abs(t - target)))

    def _cut_by_index(t, y, idx, halfw):
        L = max(0, idx - halfw)
        R = min(len(t) - 1, idx + halfw)
        return t[L:R+1], y[L:R+1]

    # Gaussian and mixture
    def _gauss(t, A, mu, sigma):
        s = max(float(sigma), 1e-12)
        return float(A) * np.exp(-0.5 * ((t - float(mu)) / s) ** 2)

    def _gmix(t, *p):
        b = p[0]
        k = (len(p) - 1) // 3
        y = np.full_like(t, b, dtype=float)
        for i in range(k):
            A, mu, s = p[1 + 3*i : 1 + 3*(i+1)]
            y += _gauss(t, A, mu, s)
        return y

    def _init_gmm(tw, yw, kmax=3):
        b0 = float(np.percentile(yw, 5)) if len(yw) else 0.0
        ywb = yw - b0
        ywb[ywb < 0] = 0.0

        prom = 0.05 * (np.max(ywb) if len(ywb) and np.max(ywb) > 0 else 1.0)
        if len(ywb) >= 3:
            peaks, _ = find_peaks(ywb, prominence=prom)
        else:
            peaks = np.array([], dtype=int)
        if len(peaks) == 0:
            peaks = np.array([int(np.argmax(ywb))])

        order = np.argsort(ywb[peaks])[::-1][:kmax]
        peaks = peaks[order]

        widths_pts = peak_widths(ywb, peaks, rel_height=0.5)[0] if len(peaks) else np.array([max(3, len(ywb)//10)])
        dt = np.median(np.diff(tw)) if len(tw) > 1 else 1.0

        p0 = [b0]; lo = [0.0]; hi = [np.inf]
        for i, pk in enumerate(peaks):
            A = float(max(ywb[pk], 1e-6))
            mu = float(tw[pk])
            w = widths_pts[i] if i < len(widths_pts) else max(3, len(ywb)//10)
            sigma = float((w * dt) / (2.0 * np.sqrt(2.0 * np.log(2.0))))
            sigma = max(sigma, 0.5 * dt)
            p0 += [A, mu, sigma]
            lo += [0.0, tw[0], 0.25 * dt]
            hi += [np.inf, tw[-1], max(1e3*dt, tw[-1]-tw[0])]
        return np.array(p0, float), (np.array(lo, float), np.array(hi, float))

    def _fit_gmm(tw, yw, kmax=3):
        if len(tw) < 3:
            return np.array([0.0, np.max(yw) if len(yw) else 0.0, float(tw[len(tw)//2]) if len(tw) else 0.0, 1.0]), 1
        p0, bounds = _init_gmm(tw, yw, kmax=kmax)
        try:
            popt, _ = curve_fit(_gmix, tw, yw, p0=p0, bounds=bounds, maxfev=20000)
        except Exception:
            p0, bounds = _init_gmm(tw, yw, kmax=1)
            try:
                popt, _ = curve_fit(_gmix, tw, yw, p0=p0, bounds=bounds, maxfev=20000)
            except Exception:
                popt = p0
        k = (len(popt) - 1) // 3
        return popt, k

    def _select_component(popt, peak_rt):
        k = (len(popt) - 1) // 3
        mus = np.array([popt[1 + 3*i + 1] for i in range(k)], float)
        idx = int(np.argmin(np.abs(mus - peak_rt)))
        A = float(popt[1 + 3*idx + 0])
        mu = float(popt[1 + 3*idx + 1])
        s  = float(max(popt[1 + 3*idx + 2], 1e-12))
        return A, mu, s

    def _shapes_on_common_axis(mu1, s1, mu2, s2, span=4.0, N=201, shift_grid=0.5):
        """
        Common normalized axis x in [-1,1] mapped to RT via mu1 and max(s1,s2).
        Returns the area-minimizing pair by allowing a small mu2 shift.
        """
        x = np.linspace(-1.0, 1.0, N)
        scale = span * max(s1, s2)
        rt = mu1 + x * scale

        # try small shifts on mu2 to absorb tiny RT drift
        deltas = np.linspace(-shift_grid*s1, shift_grid*s1, 9)
        best = (np.inf, rt, None, None)
        for d in deltas:
            y1 = np.exp(-0.5 * ((rt - mu1) / s1) ** 2)
            y2 = np.exp(-0.5 * ((rt - (mu2 + d)) / s2) ** 2)
            y1 /= max(np.max(y1), 1e-12)
            y2 /= max(np.max(y2), 1e-12)
            area = float(np.trapz(np.abs(y1 - y2), x))  # dimensionless
            if area < best[0]:
                best = (area, rt, y1, y2)
        return best  # area, rt_grid, y1, y2

    # ----------------- main -----------------
    xic1 = copy.deepcopy(xic1_original)
    xic2 = copy.deepcopy(xic2_original)
    t1, y1 = np.asarray(xic1[0], float), np.asarray(xic1[1], float)
    t2, y2 = np.asarray(xic2[0], float), np.asarray(xic2[1], float)

    idx1 = _idx_near(t1, peak_rt)
    idx2 = _idx_near(t2, peak_rt)

    sg_win = _odd(2*int(peakwidth)+1)
    y1s = _smooth_nonneg(y1, sg_win, order=2)
    y2s = _smooth_nonneg(y2, sg_win, order=2)

    t1w, y1w = _cut_by_index(t1, y1s, idx1, int(peakwidth))
    t2w, y2w = _cut_by_index(t2, y2s, idx2, int(peakwidth))

    # baseline to >=0 for fitting
    if len(y1w): y1w = y1w - np.min(y1w); y1w[y1w < 0] = 0.0
    if len(y2w): y2w = y2w - np.min(y2w); y2w[y2w < 0] = 0.0

    popt1, _ = _fit_gmm(t1w, y1w, kmax=3)
    popt2, _ = _fit_gmm(t2w, y2w, kmax=3)

    A1, mu1, s1 = _select_component(popt1, peak_rt)
    A2, mu2, s2 = _select_component(popt2, peak_rt)

    # evaluate both components on a common axis tied to the reference component
    area, rt_grid, y1c, y2c = _shapes_on_common_axis(mu1, s1, mu2, s2, span=4.0, N=201, shift_grid=0.5)

    if not np.isfinite(area):
        area = -1.0

    # return arrays used for the scoring (equal length, common RT grid)
    return float(area), list(rt_grid), list(y1c), list(rt_grid), list(y2c)

def compare_peak_shape_similarity(xic1_original, xic2_original, peak_rt, peakwidth=10, debug_output=False):
    # xic = (t, y, meta?)  ; uses t and y only
    xic1 = copy.deepcopy(xic1_original)
    xic2 = copy.deepcopy(xic2_original)
    t1, y1 = xic1[0], xic1[1]
    t2, y2 = xic2[0], xic2[1]

    # index near target
    index = xic1[0].index(min(xic1[0], key=lambda x: abs(peak_rt - x)))

    # Savitzky-Golay (ensure odd window and within bounds)
    sg_order = 2
    sg_window = max(5, 2*max(1, peakwidth//2)+1)
    sg_window = min(_odd_below(len(y1)), sg_window)
    sg_window = max(sg_order+2 | 1, sg_window)  # ensure >= order+2 and odd
    y1s = scipy.signal.savgol_filter(y1, sg_window, sg_order, mode="nearest")
    y2s = scipy.signal.savgol_filter(y2, sg_window, sg_order, mode="nearest")
    y1s = [y if y>0 else 0 for y in y1s]
    y2s = [y if y>0 else 0 for y in y2s]
    xic1[1] = y1s
    xic2[1] = y2s

    # cutouts (index-based window)
    front = max(1, index - peakwidth)
    back  = min(len(t1)-1, index + peakwidth)
    xic1_cutout = (t1[front:back], y1s[front:back])
    xic2_cutout = (t2[front:back], y2s[front:back])
    

    peak1_rt = xic1_cutout[0]
    peak2_rt = xic2_cutout[0]

    intensity1_at_peak_rt = y1s[index] if y1s[index] >= 0 else 0.000001
    intensity2_at_peak_rt = y2s[index] if y2s[index] >= 0 else 0.000001

    try:
        min_in_window1 = min(xic1_cutout[1])
        try:
            peakintensity1 = [((i - min_in_window1) / (y1s[index] - min_in_window1)) for i in xic1_cutout[1]]
        except:
            max_peakint1 = max(xic1_cutout[1]) if not max(xic1_cutout[1]) == 0 else 0.0001
            peakintensity1 = [((i - min_in_window1) / (max_peakint1)) for i in xic1_cutout[1]]

        min_in_window2 = min(xic2_cutout[1])
        try:
            peakintensity2 = [((i - min_in_window2) / (y2s[index] - min_in_window2)) for i in xic2_cutout[1]]
        except:
            max_peakint2 = max(xic2_cutout[1]) if not max(xic2_cutout[1]) == 0 else 0.0001
            peakintensity2 = [[((i - min_in_window2) / (max_peakint2)) for i in xic2_cutout[1]]]
    except Exception as e:
        print("Error: Going into except statement in compare_peak_shape_similarity in MS_functions.py because of: " + str(e))
        if index - peakwidth <= 2:
            peakwidth = index - 2
        if index + peakwidth >= len(y1s):
            peakwidth = (len(y1s) - index - 2)
        peakintensity1 = [i / intensity1_at_peak_rt for i in peakintensity1]
        peak1_rt = xic1[0][int(index - peakwidth):int(index + peakwidth)]
        peakintensity2 = [i / intensity2_at_peak_rt for i in peakintensity2]
        peak2_rt = xic2[0][int(index - peakwidth):int(index + peakwidth)]
    try:
        P = np.array([peak1_rt, peakintensity1]).T
        Q = np.array([peak2_rt, peakintensity2]).T
        area = similaritymeasures.area_between_two_curves(P, Q)
    except Exception as e:
        print("Error: Problem with similaritymeasures in MS_functions: " + str(e))
        area = -1
    if not (isinstance(area, float) or isinstance(area, int)):
        try:
            area = float(area)
        except Exception as e2:
            print("Error: Problem with similaritymeasures in MS_functions.py. Not instance of... " + str(e2))
            area = -1
    if area == np.nan or (str(area).lower() == "nan"):
        area = -1
        if debug_output == True:
            print("area was nan. Chaning area to: " + str(area))

    return area, peak1_rt, peakintensity1, peak2_rt, peakintensity2

def get_mode_of_spec(filter_string):
    if " d " in filter_string and "@hcd" in filter_string:
        ms_ms_masses = []
        filter_parsed = filter_string.split(" ")
        filter_parsed = [x for x in filter_parsed if "hcd" in x]
        for element in filter_parsed:
            ms_ms_masses.append(round(float(element.split("@")[0]), 2))
        return "MS/MS"
    elif " d " not in filter_string and "hcd" not in filter_string:
        return "Full scan"
    elif " d " not in filter_string and "hcd" in filter_string:
        return "AIF"


def min_deviation_between_list_elements(input_list):
    # Sort the list in ascending order
    sorted_list = sorted(input_list)
    # Initialize the minimum deviation with a large value
    min_deviation = float('inf')
    # Initialize the pair of values where the minimum deviation occurs
    min_deviation_values = None
    # Iterate over the sorted list and compare adjacent elements
    for i in range(len(sorted_list) - 1):
        deviation = abs(sorted_list[i + 1] - sorted_list[i])
        if deviation < min_deviation:
            min_deviation = deviation
            min_deviation_values = (sorted_list[i], sorted_list[i + 1])
    return [min_deviation, min_deviation_values]


def get_best_approx_for_ppm_spacing_within_peak(mass_list, worst_expected_ppm_deviation=20):
    mass_list = sorted(mass_list)
    ppm_spacing_list = [(((mass_list[i+1] - mass_list[i]) / mass_list[i]) * 1000000) for i in range(len(mass_list)-1)]
    ppm_spacing_list = [ppm for ppm in ppm_spacing_list if ppm < worst_expected_ppm_deviation]
    avg_ppm = sum(ppm_spacing_list) / len(ppm_spacing_list)
    return avg_ppm

def get_peaks_in_xy_series(x, y, sg_window=10, sg_order=3):
        if len(x) != len(y):
            print("Error in MS_functions.get_peaks_in_xy_series(): x and y have different lengths!")
            return None
        y = scipy.signal.savgol_filter(y, sg_window, sg_order, mode="nearest")
        peak_properties = scipy.signal.find_peaks(y, height=max(y)/100, distance=2, prominence=max(y)/100, width=(2, len(y)/10))
        identified_peaks = []
        for element in range(len(peak_properties[1]["peak_heights"])):
            one_peak = []
            one_peak.append(x[int(peak_properties[1]["left_ips"][element] + (peak_properties[1]["widths"][element]/2))])
            one_peak.append(peak_properties[1]["peak_heights"][element])
            one_peak.append(x[int(peak_properties[1]["left_ips"][element])])
            one_peak.append(x[int(peak_properties[1]["right_ips"][element])])
            one_peak.append(peak_properties[1]["prominences"][element])
            identified_peaks.append(one_peak)
        #[[time, height, lefttime, righttime, prominence], [time, height, left, right, prominence], ...]
        return identified_peaks


def summarize_mass_intensity_dict(dictio, deviation=11, debug_output=True):
        dictio = {k: v for k, v in dictio.items() if v > 0}
        if debug_output == True:
            print("summarizing dict according to new method. old length of start dictio:" + str(len(dictio)))
        #sort the dictio by its keys
        dictio = dict(sorted(dictio.items(), key=lambda item: item[0]))
        old_masses_list = list(dictio.keys())
        old_abundances_list = list(dictio.values())
        new_masses_list = []
        new_abundances_list = []

        best_approx_ppm_spacing_within_peak = get_best_approx_for_ppm_spacing_within_peak(old_masses_list)
        if debug_output == True:
            print("best approx for ppm spacing within peak: " + str(best_approx_ppm_spacing_within_peak))

        while len(old_masses_list) > 0:
            try:
                remove_all_lower = False
                remove_all_upper = False
                index_of_highest_abundance = old_abundances_list.index(max(old_abundances_list))
                curr_mass = old_masses_list[index_of_highest_abundance]
                mass_lower_border = curr_mass - 4*((deviation*curr_mass)/1000000)
                mass_upper_border = curr_mass + 4*((deviation*curr_mass)/1000000)
                iteration_step_lower = 0
                integration_step_upper = 0

                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                expected_min_spacing_between_measurement_points = (best_approx_ppm_spacing_within_peak * curr_mass) / 1000000
                try:
                    while (mass_lower_border < old_masses_list[index_of_highest_abundance-(iteration_step_lower+1)]) and \
                            (last_existing_mass_within_border - old_masses_list[index_of_highest_abundance-(iteration_step_lower+1)] <= 3.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance-(iteration_step_lower+1)] <= (1.1 * last_abundance)):

                        iteration_step_lower += 1
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance - iteration_step_lower]
                        last_abundance = old_abundances_list[index_of_highest_abundance - iteration_step_lower]
                except IndexError:
                    remove_all_lower = True
                    iteration_step_lower = 0

                
                last_existing_mass_within_border = curr_mass
                last_abundance = old_abundances_list[index_of_highest_abundance]
                try:
                    while (mass_upper_border > old_masses_list[index_of_highest_abundance+integration_step_upper+1]) and \
                            (old_masses_list[index_of_highest_abundance+integration_step_upper+1] - last_existing_mass_within_border <= 3.2*expected_min_spacing_between_measurement_points) and \
                                (old_abundances_list[index_of_highest_abundance+integration_step_upper+1] <= (1.1 * last_abundance)):
                        integration_step_upper += 1
                        last_existing_mass_within_border = old_masses_list[index_of_highest_abundance + integration_step_upper]
                        last_abundance = old_abundances_list[index_of_highest_abundance + integration_step_upper]
                except IndexError:
                    remove_all_upper = True
                    integration_step_upper = 0
                
                if remove_all_lower == True and remove_all_upper == True:
                    break
                if remove_all_upper:
                    summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]) / summed_intensity
                    if summed_intensity >= 1:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_masses_list), 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, len(old_abundances_list), 1)]
                    continue
                if remove_all_lower:
                    summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper+1, 1)])
                    weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]) / summed_intensity
                    if summed_intensity >= 1:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance+integration_step_upper+1, 1)]
                    continue

                if iteration_step_lower == 0 and integration_step_upper == 0:
                    summed_intensity = old_abundances_list[index_of_highest_abundance]
                    weighted_mass_average = old_masses_list[index_of_highest_abundance]
                    if summed_intensity >= 1:
                        new_masses_list.append(weighted_mass_average)
                        new_abundances_list.append(summed_intensity)
                    old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance ]
                    old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance ]
                    continue

                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]) / summed_intensity
                if summed_intensity >= 1:
                    new_masses_list.append(weighted_mass_average)
                    new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance-iteration_step_lower, index_of_highest_abundance+integration_step_upper+1, 1)]
            except Exception as e:
                print("EXCEPTION IN summarize_mass_intensity_dict()!!!")
                print(traceback.format_exc())
                break
        outdict = dict(zip(new_masses_list, new_abundances_list))
        if debug_output == True:
            print("New length of summarized dictio: " + str(len(outdict)))
        return outdict



def summarize_mass_intensity_dict_for_isotopo_simulation(dictio, deviation=30, debug_output=True):
    dictio = dict(sorted(dictio.items(), key=lambda item: item[0]))
    old_masses_list = list(dictio.keys())
    old_abundances_list = list(dictio.values())
    new_masses_list = []
    new_abundances_list = []

    while len(old_masses_list) > 0:
        try:
            remove_all_lower = False
            remove_all_upper = False
            index_of_highest_abundance = old_abundances_list.index(max(old_abundances_list))
            curr_mass = old_masses_list[index_of_highest_abundance]
            mass_lower_border = curr_mass - ((deviation * curr_mass) / 1000000)
            mass_upper_border = curr_mass + ((deviation * curr_mass) / 1000000)

            iteration_step_lower = 0
            iteration_step_upper = 0

            try:
                while (mass_lower_border < old_masses_list[index_of_highest_abundance - (iteration_step_lower + 1)]):
                    iteration_step_lower += 1
            except IndexError:
                remove_all_lower = True
                iteration_step_lower = 0

            try:
                while (mass_upper_border > old_masses_list[index_of_highest_abundance + iteration_step_upper + 1]):
                    iteration_step_upper += 1
            except IndexError:
                remove_all_upper = True
                iteration_step_upper = 0

            if remove_all_lower == True and remove_all_upper == True:
                break
            if remove_all_upper:
                summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, len(old_abundances_list), 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, len(old_abundances_list), 1)]) / summed_intensity
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, len(old_masses_list), 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, len(old_abundances_list), 1)]
                continue
            if remove_all_lower:
                summed_intensity = sum([old_abundances_list[i] for i in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)])
                weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)]) / summed_intensity
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(0, index_of_highest_abundance + iteration_step_upper + 1, 1)]
                continue

            if iteration_step_lower == 0 and iteration_step_upper == 0:
                summed_intensity = old_abundances_list[index_of_highest_abundance]
                weighted_mass_average = old_masses_list[index_of_highest_abundance]
                new_masses_list.append(weighted_mass_average)
                new_abundances_list.append(summed_intensity)
                old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if not i == index_of_highest_abundance]
                old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if not i == index_of_highest_abundance]
                continue

            summed_intensity = sum([old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)])
            weighted_mass_average = sum([old_masses_list[i] * old_abundances_list[i] for i in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)]) / summed_intensity
            new_masses_list.append(weighted_mass_average)
            new_abundances_list.append(summed_intensity)
            old_masses_list = [old_masses_list[i] for i in range(len(old_masses_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)]
            old_abundances_list = [old_abundances_list[i] for i in range(len(old_abundances_list)) if i not in range(index_of_highest_abundance - iteration_step_lower, index_of_highest_abundance + iteration_step_upper + 1, 1)]
        except Exception as e:
            print("EXCEPTION IN summarize_mass_intensity_dict_for_isotopo_simulation()!!!")
            print(traceback.format_exc())
            break
    outdict = dict(zip(new_masses_list, new_abundances_list))
    return outdict


def read_summary_to_df(filepath, sep="\t", colnames=("MASS", "FORMULA", "INDEX", "RT", "INTENSITY", "SCORE", "FRAGS", "F_SCORES", "F_INTENSITIES", "NLS", "NLS_DEV"), header=None):
    df = pd.read_csv(filepath, sep=sep, header=header)
    new_colnames = []
    for entry in range(len(df.columns)):
        try:
            new_colnames.append(colnames[entry])
        except Exception as e:
            print("Error setting columnames in read_summary_to_df(): " + str(e))
            print(traceback.format_exc())
            new_colnames.append(entry)
    df.columns = new_colnames
    return df






if __name__ == "__main__":
    print("This module is not meant to be run directly. Please import it in your script.")
    sys.exit()

