import copy
import os
import sys
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import UVenture.class_Spec as class_Spec
import UVenture.peakdetection_funcs as peakdetection_funcs
import UVenture.MS_functions as MS_functions
import UVenture.GaussianCompressor as GaussianCompressor


def remove_isotopo_signals(peak_df, mass_deviation_ppm=12, height_deviation_isotopo=0.5):
    # Adjust the mass deviaton. 
    # Usually the mass deviation is x but because for isotopologue measurements, the mass is more or less the same, 
    # and it is in one spectrum, the mass deviation between the isotope and isotopologue is smaller.
    mass_deviation_ppm = mass_deviation_ppm / 4
    
    # Remove isotopo signals
    isotopo_signal_rows = []
    for i in peak_df["peak_bins"].unique():
        subgroup_df = peak_df[peak_df["peak_bins"] == i]
        subgroup_df = subgroup_df.sort_values(by=["height"], ascending=[True])
        for index, row in subgroup_df.iterrows():
            curr_row = row
            curr_mass = row["mass"]
            curr_height = row["height"]
            #Start checking if an isotopo signal exists for curr_row
            for index2, row2 in subgroup_df.iterrows():
                if index == index2:
                    continue
                mass_deviation_isotopo = (mass_deviation_ppm / 1000000) * curr_mass
                if (curr_mass+1.00336-mass_deviation_isotopo) <= row2["mass"] <= (curr_mass+1.00336+mass_deviation_isotopo):
                    height_ratio = row2["height"] / curr_height
                    theo_height_ratio = (((curr_mass*0.85)*0.7 ) / 12) * 0.01082 # first multiplicator is for hydrogen atoms and second is for heteroatoms
                    if theo_height_ratio * (1 - height_deviation_isotopo) <= height_ratio <= theo_height_ratio * (1 + height_deviation_isotopo):
                        print("Found isotopo signal: " + str(row) + " " + str(row2))
                        isotopo_signal_rows.append(row2)
    print("Isotopo signal rows: " + str(isotopo_signal_rows))
    indices_to_delete = peak_df.index[peak_df.apply(tuple, axis=1).isin([tuple(row) for row in isotopo_signal_rows])]
    indices_to_delete = list(set(indices_to_delete))
    print("Indices to delete: " + str(indices_to_delete))
    peak_df.drop(indices_to_delete, inplace=True)
    return peak_df

def remove_duplicates_from_peak_df(peak_df, rt_bins=100, mass_deviation_ppm=15):
    rt_threshold = 2
    def calculate_mass_threshold(mass):
        return mass * (mass_deviation_ppm / 1000000)
    peak_df = peak_df.sort_values(by=["mass", "rt"]).reset_index(drop=True)
    keep_mask = [True] * len(peak_df)

    for i, row_i in peak_df.iterrows():
        if not keep_mask[i]:
            continue
        for j in range(i + 1, len(peak_df)):  # Compare only subsequent rows
            if not keep_mask[j]:
                continue  # Skip rows already marked as duplicates
            row_j = peak_df.iloc[j]
            mass_diff = abs(row_i["mass"] - row_j["mass"])
            rt_diff = abs(row_i["rt"] - row_j["rt"])
            if rt_diff <= rt_threshold and mass_diff <= calculate_mass_threshold(row_i["mass"]):
                keep_mask[j] = False
    print("Keep mask: " + str(keep_mask))
    peak_df = peak_df[keep_mask].reset_index(drop=True)
    
    peak_df["peak_bins"] = pd.cut(peak_df["rt"], bins=rt_bins, labels=False)
    peak_df = peak_df.sort_values(by=["rt"], ascending=[True])
    return peak_df

def dedupe_peaks(peaks, dm=0.001, drt=3.0, keep="max_area"):
    """
    peaks: list of (mass, rt, area, height)
    dm: mass tolerance (±)
    drt: rt tolerance (±)
    keep: "first" | "max_area" | "max_height"
    """
    if not peaks:
        return []

    peaks = np.asarray(peaks, dtype=float)
    m = peaks[:, 0]
    rt = peaks[:, 1]
    area = peaks[:, 2]
    height = peaks[:, 3]
    n = len(peaks)

    # Disjoint set
    parent = np.arange(n)
    rank = np.zeros(n, dtype=int)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    # Sweep by mass to avoid O(n^2)
    idx = np.argsort(m)
    left = 0
    for r in range(n):
        i = idx[r]
        # shrink window so mass difference <= dm
        while m[i] - m[idx[left]] > dm:
            left += 1
        # compare to candidates in window
        for k in range(left, r):
            j = idx[k]
            if abs(rt[i] - rt[j]) <= drt:
                union(i, j)

    # Choose representative per cluster
    groups = {}
    for i in range(n):
        root = find(i)
        groups.setdefault(root, []).append(i)

    reps = []
    for root, members in groups.items():
        if keep == "first":
            rep = min(members)
        elif keep == "max_area":
            rep = max(members, key=lambda i: (area[i], height[i]))
        elif keep == "max_height":
            rep = max(members, key=lambda i: (height[i], area[i]))
        else:
            raise ValueError("keep must be one of: first, max_area, max_height")
        reps.append(rep)

    reps.sort()  # preserve original order
    return [tuple(peaks[i]) for i in reps]

def remove_isotopologues(peak_df, 
                         mass_col="mass",
                         rt_col="rt",
                         height_col="height",
                         delta=1.0036,
                         mass_deviation_ppm=12,
                         rt_tol=2):
    """
    Drop isotopologue peaks:
      - candidate mass in [mass + delta - mass_tol, mass + delta + mass_tol]
      - same retention time (exact if rt_tol == 0, else |Δrt| <= rt_tol)
      - height ratio h_iso / h_base within [0.5*theo, 1.5*theo],
        theo = (((mass_base*0.85)*0.7)/12) * 0.01082

    Returns a copy of df with isotopologue rows removed.
    """
    mass_tol = (mass_deviation_ppm / 1000000) * peak_df[mass_col].to_numpy(float)
    if peak_df.empty:
        return peak_df.copy()

    d = peak_df.reset_index(drop=False).rename(columns={"index": "_orig"})
    m = d[mass_col].to_numpy(float)
    rt = d[rt_col].to_numpy(float)
    h  = d[height_col].to_numpy(float)

    n = len(d)
    keep = np.ones(n, dtype=bool)

    # sort by mass for windowed search
    sidx = np.argsort(m)
    m_s  = m[sidx]
    rt_s = rt[sidx]
    h_s  = h[sidx]
    idx_s = sidx  # map sorted position -> original row

    # vectorized candidate windows for +delta
    low  = m_s + (delta - mass_tol)
    high = m_s + (delta + mass_tol)
    L = np.searchsorted(m_s, low, side="left")
    R = np.searchsorted(m_s, high, side="right")

    # iterate bases in increasing mass
    for ii in range(n):
        if not keep[idx_s[ii]]:
            continue  # base already removed earlier (rare)
        mass_base = m_s[ii]
        rt_base   = rt_s[ii]
        h_base    = h_s[ii]
        if not np.isfinite(h_base) or h_base <= 0:
            continue

        # theoretical height ratio for this base
        theo = (((mass_base * 0.85) * 0.7) / 12.0) * 0.01082
        lo_ratio = 0.5 * theo
        hi_ratio = 1.5 * theo

        for jj in range(L[ii], R[ii]):
            j_orig = idx_s[jj]
            if not keep[j_orig]:
                continue
            # retention time check
            if rt_tol == 0.0:
                same_rt = (rt_s[jj] == rt_base)
            else:
                same_rt = (abs(rt_s[jj] - rt_base) <= rt_tol)
            if not same_rt:
                continue

            # height ratio check
            ratio = h_s[jj] / h_base if h_base != 0 else np.inf
            if lo_ratio <= ratio <= hi_ratio:
                # jj is isotopologue of ii -> remove heavier peak
                keep[j_orig] = False

    out = d[keep].drop(columns=["_orig"]).reset_index(drop=True)
    return out

def _gauss(t, A, mu, sigma, c):
    return A * np.exp(-0.5 * ((t - mu) / sigma) ** 2) + c

class FindPeaks:
    def __init__(self, ms_file, settings_dict, mass_range=1, 
                    threshold_area=400000, threshold_intensity=100000,
                    min_peak_width=4, max_peak_width=40,
                    peaklist_filename="", mzrt_filename="", **kwargs) -> None:
        default_kwargs = {
            "mass_deviation": 15,
            "deduplication_retentiontime": 5,
            "height_deviation_isotopo": 0.5,
            "rt_deviation_isotopo": 2,
            "max_gaussian_fits": 10,
            "target_rel_area": 0.0002,
            "min_neighbour_spectra_above_noise": 1
        }
        self.kwargs = {**default_kwargs, **kwargs}
        self.ms_file = ms_file
        if not isinstance(mass_range, int):
            print("mass_range must be an integer. Setting mass_range to 1")
            self.mass_range = 1
        else:
            self.mass_range = mass_range
        self.threshold_area = threshold_area
        self.threshold_intensity = threshold_intensity
        self.max_peak_width = max_peak_width
        self.min_peak_width = min_peak_width


        min_mz, max_mz = ms_file.mz_range
        self.min_mz = int(min_mz)
        self.max_mz = int(max_mz)

        self.peak_df = pd.DataFrame(columns=["mass", "rt", "height", "area"])

        for mass in range(self.max_mz, self.min_mz-1, -self.mass_range):
            print("Processing mass: " + str(mass))
            xic = ms_file.get_xic(mass=mass, mass_deviation=mass_range, requested_filter_mode="Full scan")
            times = np.array(xic[0])
            intensities = np.array(xic[1])
            true_indices_of_entries = np.array(xic[2])

            peak_areas_gaussian_compressor = self.get_peaks_gaussian_compressor(mass, times, intensities)

            peak_areas_scipy_method = self.get_peaks_scipy_method(mass, times, intensities, true_indices_of_entries)

            peak_areas = peak_areas_scipy_method + peak_areas_gaussian_compressor

            print("Combined peak areas: " + str(peak_areas))
            peak_areas = self.remove_duplicates_in_peak_areas(mass, peak_areas)
            print("Combined peak areas after first deduplication: " + str(peak_areas))

            # Add the peak areas to the DataFrame
            print("Adding peak areas to DataFrame.")
            if len(peak_areas) > 0:
                new_peak_df = pd.DataFrame(peak_areas, columns=["mass", "rt", "height", "area"])
                if len(self.peak_df) >= 1:
                    self.peak_df = pd.concat([self.peak_df, new_peak_df], ignore_index=True)
                elif len(self.peak_df) == 0:
                    self.peak_df = new_peak_df

            if len(self.peak_df) <= 2:
                continue

            # Remove duplicates
            #print("Peak df before removing duplicates: " + str(peak_df))
            print("Length of peak df before removing duplicates: " + str(len(self.peak_df)))
            peak_areas_list = self.peak_df[["mass", "rt", "height", "area"]].values.tolist()
            peak_areas_list = self.remove_duplicates_in_peak_areas(mass, peak_areas_list)
            self.peak_df = pd.DataFrame(peak_areas_list, columns=["mass", "rt", "height", "area"])
            print("Length of peak df after removing duplicates: " + str(len(self.peak_df)))


            # Remove isotopo signals
            print("Peak df before removing isotopo signals: " + str(self.peak_df))
            self.peak_df = remove_isotopologues(self.peak_df, mass_col="mass", rt_col="rt", height_col="height", delta=1.0036, mass_deviation_ppm=self.kwargs["mass_deviation"], rt_tol=self.kwargs["rt_deviation_isotopo"])
            print("Peak df after removing duplicates and isotopo signals: " + str(self.peak_df))

            
            # Remove peaks that do not have at least 3 spectra above the intensity threshold
            self.peak_df = self.remove_peaks_with_low_neighbors(self.peak_df, self.threshold_intensity, min_neighbors=self.kwargs["min_neighbour_spectra_above_noise"])

            new_peak_df = self.peak_df.copy()
            for index, row in self.peak_df.iterrows():
                if not row["mass"] >= mass + 5:
                    print("Skipping mass: " + str(row["mass"]) + " because it is not >= " + str(mass +5))
                    continue

                #Check if a peak exists at row["rt"]
                peak_exists, properties = self.check_if_peak_exists(self.ms_file, row["mass"], row["rt"])
                peakheight = properties[2] if peak_exists else 0
                peakarea = properties[3] if peak_exists else 0
                new_peak_df = new_peak_df.drop(index)


                if peak_exists == False:
                    print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    continue
                elif peakheight < self.threshold_intensity:
                    print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " because peak height is below threshold.")
                    continue
                else:
                    print("Peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    if not peaklist_filename == "":
                        # Save the peak to the peaklist file
                        with open(peaklist_filename, "a") as f:
                            f.write(f"{row['mass']}\t{row['rt']}\t{round(peakarea, 2)}\t{round(peakheight, 2)}\n")
                        print("Peak saved to peaklist: " + str(peaklist_filename) + ".")
                    if not mzrt_filename == "":
                        # Save the peak to the mzrt file to be processed directly
                        with open(mzrt_filename, "a") as f:
                            file = os.path.normpath(ms_file.filename)
                            file = file.split(os.sep)[-1]
                            print("Taskstorage file: " + str(file))
                            f.write(f"{file}\t{row['mass']}\t{row['rt']}\t{settings_dict}\n")
                            print("Peak saved to mzrt file: " + str(mzrt_filename) + ".")
                        continue
            self.peak_df = new_peak_df.copy()
            print("Length of peak df at the end of mass " + str(mass) + ": " + str(len(self.peak_df)))


    def assess_gaussian_peak(self, t, y, rt, window=None, min_points=4):
        """
        Assess whether the signal near `rt` is Gaussian-like.

        Parameters
        ----------
        t : array-like
            Timestamps (1D, monotonic).
        y : array-like
            Intensities (1D, same length as t).
        rt : float
            Retention time to assess.
        window : float, optional
            Half-width of the fitting region in same units as t.
            If None, uses 0.05 * (t.max() - t.min()).
        min_points : int
            Minimum points required in the window.

        Returns
        -------
        result : dict
            {
            'score': float in [0,1], higher is more Gaussian-like,
            'r2': coefficient of determination of the fit,
            'fwhm': estimated full width at half maximum,
            'sigma': sigma of Gaussian,
            'amplitude': A,
            'baseline': c,
            'mu': fitted center,
            'n_points': number of samples used
            }
            On failure, returns NaNs with score=0 and r2=-inf.
        """
        t = np.asarray(t, float)
        y = np.asarray(y, float)
        if t.ndim != 1 or y.ndim != 1 or t.size != y.size:
            raise ValueError("t and y must be 1D arrays of the same length.")

        if window is None:
            window = 0.05 * (t.max() - t.min())
            if window <= 0:
                return {'score': 0.0, 'r2': float('-inf'), 'fwhm': np.nan, 'sigma': np.nan,
                        'amplitude': np.nan, 'baseline': np.nan, 'mu': np.nan, 'n_points': 0}

        mask = np.abs(t - rt) <= window
        if mask.sum() < min_points:
            return {'score': 0.0, 'r2': float('-inf'), 'fwhm': np.nan, 'sigma': np.nan,
                    'amplitude': np.nan, 'baseline': np.nan, 'mu': np.nan, 'n_points': int(mask.sum())}

        tx = t[mask]
        yx = y[mask]

        # Initial guesses
        i_peak = np.argmax(yx)
        A0 = max(yx[i_peak] - np.median(yx), np.ptp(yx) * 0.5)
        mu0 = float(tx[i_peak])
        # crude sigma guess: span where signal above half dynamic range
        half = np.median(yx) + 0.5 * (np.max(yx) - np.median(yx))
        above = tx[yx >= half]
        if above.size >= 2:
            sigma0 = (above[-1] - above[0]) / 2.355 if (above[-1] > above[0]) else max(window * 0.2, 1e-6)
        else:
            sigma0 = max(window * 0.2, 1e-6)
        c0 = float(np.median(yx))

        p0 = [A0, mu0, sigma0, c0]
        # Bounds: positive amplitude, sigma>0; allow small center drift within window
        bounds = ([0.0, rt - window, 1e-6, -np.inf],
                [np.inf, rt + window, np.inf, np.inf])

        try:
            popt, _ = curve_fit(_gauss, tx, yx, p0=p0, bounds=bounds, maxfev=20000)
            A, mu, sigma, c = map(float, popt)
            yfit = _gauss(tx, A, mu, sigma, c)

            # R^2
            ss_res = float(np.sum((yx - yfit) ** 2))
            ss_tot = float(np.sum((yx - np.mean(yx)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else -np.inf

            # Normalized RMSE score in [0,1]
            denom = max(np.max(yx) - np.min(yx), 1e-12)
            nrmse = np.sqrt(np.mean((yx - yfit) ** 2)) / denom
            score = float(np.clip(1.0 - nrmse, 0.0, 1.0))

            fwhm = 2.354820045 * sigma
            return {
                'score': score,
                'r2': r2,
                'fwhm': fwhm,
                'sigma': sigma,
                'amplitude': A,
                'baseline': c,
                'mu': mu,
                'n_points': int(tx.size),
            }
        except Exception:
            return {'score': 0.0, 'r2': float('-inf'), 'fwhm': np.nan, 'sigma': np.nan,
                    'amplitude': np.nan, 'baseline': np.nan, 'mu': np.nan, 'n_points': int(tx.size)}

    def check_if_peak_exists(self, ms_file, mass, rt, gaussian_threshold=0.93):
        mass_dev_xic = (self.kwargs["mass_deviation"] * mass) / 1000000
        xic = ms_file.get_xic(mass=mass, mass_deviation=mass_dev_xic, requested_filter_mode="Full scan")
        times = np.array(xic[0])
        intensities = np.array(xic[1])
        true_indices_of_entries = np.array(xic[2])

        gaussian_fit = self.assess_gaussian_peak(times, intensities, rt)
        print("Gaussian fit results for mass: " + str(mass) + " RT: " + str(rt))
        print(gaussian_fit)
        if gaussian_fit["score"] >= gaussian_threshold and \
            gaussian_fit["fwhm"]*1.5 >= self.min_peak_width and \
                gaussian_fit["fwhm"] <= self.max_peak_width and \
                    gaussian_fit["amplitude"] >= self.threshold_intensity:
            peakarea = gaussian_fit["amplitude"] * gaussian_fit["fwhm"] * (np.sqrt(2.0 * np.pi) / 2.0)
            peakheight = gaussian_fit["amplitude"]
            return True, [mass, rt, peakheight, peakarea]
        else:
            return False, [mass, rt, 0, 0]

    def remove_peaks_with_low_neighbors(self, peak_df, threshold_intensity, min_neighbors=1):
        if not isinstance(min_neighbors, int) or min_neighbors < 1:
            print("min_neighbors must be a positive integer. Setting min_neighbors to 1.")
            min_neighbors = 1
        new_peak_df = peak_df.copy()
        for index, row in peak_df.iterrows():
            mass_deviation_amu = (self.kwargs["mass_deviation"] * row["mass"]) / 1000000
            xic = self.ms_file.get_xic(mass=row["mass"], mass_deviation=mass_deviation_amu, requested_filter_mode="Full scan")
            times = np.array(xic[0])
            intensities = np.array(xic[1])
            rt_index_in_times = np.argmin(np.abs(times - row["rt"]))
            neighboring_intensities = []
            for i in range(-min_neighbors, min_neighbors+1, 1):
                neighboring_intensities.append(intensities[rt_index_in_times + i])
            
            if not all(i >= threshold_intensity for i in neighboring_intensities):
                print("Removing peak at mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " because neighboring intensities are below threshold.")
                new_peak_df = new_peak_df.drop(index)
        return new_peak_df

    def get_peaks_gaussian_compressor(self, mass, times, intensities):
        gaussian_compressor_results = GaussianCompressor.get_peaks_in_chromatogram(times, intensities,
                                                                                   max_peaks=self.kwargs["max_gaussian_fits"], target_rel_area=self.kwargs["target_rel_area"])
        # A, mu, sigma
        # Remove gaussian peaks that are too small or too wide
        remaining_params = []
        peak_areas_gaussian_compressor = []
        for p in gaussian_compressor_results:
            height_of_peak = p[0] / (p[2] * np.sqrt(2.0 * np.pi))
            fwhm_of_peak = 2.355 * p[2]
            area = p[0]
            rt = p[1]
            if height_of_peak >= self.threshold_intensity and fwhm_of_peak <= (self.max_peak_width/2) and fwhm_of_peak >= self.min_peak_width:
                index_of_rt = class_Spec.get_index_of_spectrum_closest_to_rt(self.ms_file, rt)
                print("Index of RT: " + str(index_of_rt) + " RT: " + str(rt) + " Area: " + str(area))
                spec = class_Spec.Spec(self.ms_file, index_of_rt, spec_requested_filter_mode="Full scan", mass_deviation=self.kwargs["mass_deviation"])
                interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(self.mass_range*1.1) and k < mass+(self.mass_range*1.1)}
                interesting_range = {k: v for k, v in interesting_range.items() if v > self.threshold_intensity}
                for k, v in interesting_range.items():
                    print("Mass: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                    peak_areas_gaussian_compressor.append((k, rt, v, area))
                remaining_params.append(p)
        print("Peak areas from Gaussian Compressor: " + str(peak_areas_gaussian_compressor))
        return peak_areas_gaussian_compressor

    def get_peaks_scipy_method(self, mass, times, intensities, true_indices_of_entries):
        peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=self.min_peak_width, max_width_seconds=self.max_peak_width)
        print("Found peaks: " + str(peaks) + ".  With the old method")
        #Calculate peak area
        peak_areas = []
        for peak in range(len(peaks)):
            left = int(properties["left_ips"][peak])
            right = int(properties["right_ips"][peak])
            area = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
            if area <= self.threshold_area:
                continue
            spec = class_Spec.Spec(self.ms_file, true_indices_of_entries[peaks[peak]], spec_requested_filter_mode="Full scan", mass_deviation=self.kwargs["mass_deviation"])
            interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(self.mass_range*1.1) and k < mass+(self.mass_range*1.1)}
            interesting_range = {k: v for k, v in interesting_range.items() if v > self.threshold_intensity}
            for k, v in interesting_range.items():
                print("Mass: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                peak_areas.append((k, times[peaks[peak]], v, area))
        print("Peak areas from peak detection: " + str(peak_areas))
        return peak_areas

    def remove_duplicates_in_peak_areas(self, mass, peak_areas):
        # Remove duplicates in peak_areas based on mass and rt
        scaling_factor_mass_resolution = ((mass - 100)/100) ** 2
        if scaling_factor_mass_resolution < 1:
            scaling_factor_mass_resolution = 1
        peak_areas = dedupe_peaks(peak_areas, dm=((self.kwargs["mass_deviation"] * (mass*scaling_factor_mass_resolution)) / 1000000), drt=self.kwargs["deduplication_retentiontime"], keep="max_area")
        return peak_areas



def __old__get_all_possible_peaks(ms_file, settings_dict, mass_range=1, 
                             threshold_area=400000, threshold_intensity=100000,
                             rt_bins=400, mass_deviation_isotopo=2, height_deviation_isotopo=0.5,
                             min_peak_width=4, max_peak_width=40,
                             peaklist_filename="", mzrt_filename="", **kwargs):
    
    default_kwargs = {
        "mass_deviation": 15,
        "deduplication_retentiontime": 5,
    }
    kwargs = {**default_kwargs, **kwargs}
    
    if not isinstance(mass_range, int):
        print("mass_range must be an integer. Setting mass_range to 1")
        mass_range = 1
    min_mz, max_mz = ms_file.mz_range
    min_mz = int(min_mz)
    max_mz = int(max_mz)

    peak_df = pd.DataFrame(columns=["mass", "rt", "height", "area"])
    for mass in range(max_mz, min_mz-1, -mass_range):
        print("Processing mass: " + str(mass))
        xic = ms_file.get_xic(mass=mass, mass_deviation=mass_range, requested_filter_mode="Full scan")
        times = np.array(xic[0])
        intensities = np.array(xic[1])
        true_indices_of_entries = np.array(xic[2])
        

        gaussian_compressor_results = GaussianCompressor.get_peaks_in_chromatogram(times, intensities,
                                                                                   max_peaks=10, target_rel_area=0.0002)
        # A, mu, sigma
        # Remove gaussian peaks that are too small or too wide
        remaining_params = []
        peak_areas_gaussian_compressor = []
        for p in gaussian_compressor_results:
            height_of_peak = p[0] / (p[2] * np.sqrt(2.0 * np.pi))
            fwhm_of_peak = 2.355 * p[2]
            area = p[0]
            rt = p[1]
            if height_of_peak >= threshold_intensity and fwhm_of_peak <= (max_peak_width/2) and fwhm_of_peak >= min_peak_width:
                index_of_rt = class_Spec.get_index_of_spectrum_closest_to_rt(ms_file, rt)
                print("Index of RT: " + str(index_of_rt) + " RT: " + str(rt) + " Area: " + str(area))
                spec = class_Spec.Spec(ms_file, index_of_rt, spec_requested_filter_mode="Full scan", mass_deviation=kwargs["mass_deviation"])
                interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(mass_range*1.1) and k < mass+(mass_range*1.1)}
                interesting_range = {k: v for k, v in interesting_range.items() if v > threshold_intensity}
                for k, v in interesting_range.items():
                    print("Mass: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                    peak_areas_gaussian_compressor.append((k, rt, v, area))
                remaining_params.append(p)
        print("Found peaks: " + str(len(gaussian_compressor_results)) + ". With the Gaussian Compressor method")
        print("Peak areas from Gaussian Compressor: " + str(peak_areas_gaussian_compressor))

        
        peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=min_peak_width, max_width_seconds=max_peak_width)
        print("Found peaks: " + str(peaks) + ".  With the old method")
        #Calculate peak area
        peak_areas = []
        for peak in range(len(peaks)):
            left = int(properties["left_ips"][peak])
            right = int(properties["right_ips"][peak])
            area = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
            if area <= threshold_area:
                continue
            spec = class_Spec.Spec(ms_file, true_indices_of_entries[peaks[peak]], spec_requested_filter_mode="Full scan", mass_deviation=kwargs["mass_deviation"])
            interesting_range = {k: v for k, v in spec.summarized_mass_intensity_dict.items() if k > mass-(mass_range*1.1) and k < mass+(mass_range*1.1)}
            interesting_range = {k: v for k, v in interesting_range.items() if v > threshold_intensity}
            for k, v in interesting_range.items():
                print("Mass: " + str(k) + " Height: " + str(v) + " Area: " + str(area) + " RT: " + str(spec.rt))
                peak_areas.append((k, times[peaks[peak]], v, area))
        print("Peak areas from peak detection: " + str(peak_areas))

        peak_areas = peak_areas + peak_areas_gaussian_compressor

        # Remove duplicates in peak_areas based on mass and rt
        print("Length of peak areas before deduplication: " + str(len(peak_areas)))
        scaling_factor_mass_resolution = ((mass - 100)/100) ** 2
        if scaling_factor_mass_resolution < 1:
            scaling_factor_mass_resolution = 1
        peak_areas = dedupe_peaks(peak_areas, dm=((kwargs["mass_deviation"] * (mass*scaling_factor_mass_resolution)) / 1000000), drt=kwargs["deduplication_retentiontime"], keep="max_area")
        print("Length of peak areas after deduplication: " + str(len(peak_areas)))
        
        # Add the peak areas to the DataFrame
        if len(peak_areas) > 0:
            new_peak_df = pd.DataFrame(peak_areas, columns=["mass", "rt", "height", "area"])
            if len(peak_df) >= 1:
                peak_df = pd.concat([peak_df, new_peak_df], ignore_index=True)
            elif len(peak_df) == 0:
                peak_df = new_peak_df

        if len(peak_df) <= 2:
            continue

        # Remove duplicates
        #print("Peak df before removing duplicates: " + str(peak_df))
        peak_df = remove_duplicates_from_peak_df(peak_df, rt_bins=rt_bins, mass_deviation_ppm=7)
        # Remove isotopo signals
        peak_df = remove_isotopo_signals(peak_df, mass_deviation_ppm=kwargs["mass_deviation"], height_deviation_isotopo=height_deviation_isotopo)
        print("Peak df after removing duplicates and isotopo signals: " + str(peak_df))

        
        new_peak_df = peak_df.copy()
        for index, row in peak_df.iterrows():
            if not row["mass"] >= mass + 5:
                continue
            new_peak_df = new_peak_df.drop(index)
            print("Mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " Height: " + str(row["height"]))
            xic_mass_deviation = (mass_deviation_isotopo * row["mass"]) / 1000000
            xic_mass_deviation = xic_mass_deviation * 3
            xic = ms_file.get_xic(mass=row["mass"], mass_deviation=xic_mass_deviation, requested_filter_mode="Full scan")
            times = np.array(xic[0])
            intensities = np.array(xic[1])
            true_indices_of_entries = np.array(xic[2])

            # Check if the neighboring scans are also above the intensity threshold
            rt_index_in_times = np.argmin(np.abs(times - row["rt"]))
            neighboring_intensities = [intensities[rt_index_in_times-1], intensities[rt_index_in_times], intensities[rt_index_in_times+1]]
            if not all(i >= threshold_intensity for i in neighboring_intensities):
                print("Skipping peak at mass: " + str(row["mass"]) + " RT: " + str(row["rt"]) + " because neighboring intensities are below threshold.")
                continue

            peaks, properties, smoothed_intensity = peakdetection_funcs.get_peaks_with_smooth_and_bgsubst(times, intensities, min_width_seconds=min_peak_width, max_width_seconds=max_peak_width)
            
            #Check if a peak exists at row["rt"]
            if len(peaks) == 0:
                print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                continue
            else:
                peak_found = False
                for peak in range(len(peaks)):
                    left = int(properties["left_ips"][peak])
                    right = int(properties["right_ips"][peak])
                    if times[left] <= row["rt"] <= times[right]:
                        # A peak was found at the specified RT
                        # Peak will be added to the peaklist
                        print("Peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                        new_row = copy.deepcopy(row)
                        new_row["area"] = np.trapz(smoothed_intensity[left:right], dx=(times[1] - times[0]))
                        peak_found = True
                        if not peaklist_filename == "":
                            # Save the peak to the peaklist file
                            with open(peaklist_filename, "a") as f:
                                f.write(f"{new_row['mass']}\t{new_row['rt']}\t{round(new_row['area'], 2)}\t{round(new_row['height'], 2)}\t{times[left]}\t{times[right]}\n")
                            print("Peak saved to peaklist: " + str(peaklist_filename) + ".")
                        if not mzrt_filename == "":
                            # Save the peak to the mzrt file to be processed directly
                            with open(mzrt_filename, "a") as f:
                                file = os.path.normpath(ms_file.filename)
                                file = file.split(os.sep)[-1]
                                print("Taskstorage file: " + str(file))
                                f.write(f"{file}\t{new_row['mass']}\t{new_row['rt']}\t{settings_dict}\n")
                            print("Peak saved to mzrt file: " + str(mzrt_filename) + ".")
                        break
                if not peak_found:
                    print("No peak found for mass: " + str(row["mass"]) + " RT: " + str(row["rt"]))
                    continue
        peak_df = new_peak_df.copy()
    return True


if __name__ == "__main__":

    print("This module is not meant to be run directly. Please use it as part of the UVenture package.")
    sys.exit()
