import copy
import numpy as np
import UVenture.class_MS_file as class_MS_file
import UVenture.MS_functions as MS_functions


def sum_multiple_mzmlfiles(mzml_filename_list, rt_tolerance=None, mass_tolerance=1, log_level="-"):
    """
    Combine N mzML files into one MS_File.
    - Sum spectra when retention times match within rt_tolerance and the mode matches.
    - Sum total ion current.
    - Append unmatched spectra only if they occur after the end of the first file.
    - Recompute rt_list, tic, filters, modes, method_duration, rt_range, mz_range.
    mass_tolerance should be given in ppm
    """
    spectra_weight = len(mzml_filename_list)

    if not mzml_filename_list:
        raise ValueError("mzml_filename_list is empty")

    base = class_MS_file.MS_File(mzml_filename_list[0], log_level=log_level)
    for spec in range(len(base.rawdata)):
        base.rawdata[spec]["intensity array"] = base.rawdata[spec]["intensity array"] / len(mzml_filename_list)
        base.rawdata[spec]["total ion current"] = base.rawdata[spec]["total ion current"] / len(mzml_filename_list)

    base_mz_array_single_weighed_start = []
    base_intensity_array_single_weighed_start = []
    base_tic_single_weighed_start = []
    for spec in range(len(base.rawdata)):
        base_mz_array_single_weighed_start.append(base.rawdata[spec]["m/z array"])
        base_intensity_array_single_weighed_start.append(base.rawdata[spec]["intensity array"])
        base_tic_single_weighed_start.append(base.rawdata[spec]["total ion current"])

    def _nearest_index(sorted_array, value):
        arr = np.asarray(sorted_array)
        j = int(np.searchsorted(arr, value))
        cand = []
        if j < arr.size:
            cand.append(j)
        if j > 0:
            cand.append(j - 1)
        if not cand:
            return None
        return min(cand, key=lambda idx: abs(arr[idx] - value))

    def _merge_peaks(mz_a, int_a, mz_b, int_b, tol):
        # tol should be given in ppm. 
        if len(mz_a) == 0:
            return list(mz_b), list(int_b)
        if len(mz_b) == 0:
            return list(mz_a), list(int_a)
        mz = np.concatenate([np.asarray(mz_a), np.asarray(mz_b)])
        inten = np.concatenate([np.asarray(int_a), np.asarray(int_b)])
        order = np.argsort(mz)
        mz = mz[order]
        inten = inten[order]
        out_mz, out_int = [], []
        cur_mz = float(mz[0])
        cur_int = float(inten[0])
        for k in range(1, mz.size):
            if abs(mz[k] - cur_mz) <= ((tol * cur_mz)/1000000):
                total_int = cur_int + float(inten[k])
                if total_int > 0:
                    cur_mz = (cur_mz * cur_int + float(mz[k]) * float(inten[k])) / total_int
                cur_int = total_int
            else:
                out_mz.append(cur_mz)
                out_int.append(cur_int)
                cur_mz = float(mz[k])
                cur_int = float(inten[k])
        out_mz.append(cur_mz)
        out_int.append(cur_int)
        return out_mz, out_int

    for fname in mzml_filename_list[1:]:
        new = class_MS_file.MS_File(fname, log_level=log_level)

        base_rt = base.rt_list
        new_rt = new.rt_list

        if rt_tolerance is None:
            d1 = np.median(np.diff(base_rt)) if len(base_rt) > 1 else 0.0
            d2 = np.median(np.diff(new_rt)) if len(new_rt) > 1 else 0.0
            curr_rt_tol = 2 * max(d1, d2, 1e-9)
        else:
            curr_rt_tol = float(rt_tolerance)

        used_new = np.zeros(len(new_rt), dtype=bool)

        # sum overlapping spectra by nearest RT and same mode
        for i, curr_rt_base in enumerate(base_rt):
            base_filter = base.rawdata[i]["scanList"]["scan"][0]["filter string"]

            # Get the best matching spectrum of the new MS File for the current spectrum of the base MS File.
            j = _nearest_index(new_rt, curr_rt_base)
            best_idx = 0
            for new_spectra_idx in range(0, len(new.rawdata)-1, 1):
                if new_spectra_idx >= len(base_rt) or new_spectra_idx >= len(new.rawdata) or new_spectra_idx <= 0:
                    continue
                new_filter = new.rawdata[new_spectra_idx]["scanList"]["scan"][0]["filter string"]
                if not new_filter == base_filter:
                    continue
                else:
                    if abs(new_spectra_idx - j) < abs(best_idx - j):
                        best_idx = new_spectra_idx
            j = best_idx
            if abs(new_rt[j] - curr_rt_base) > curr_rt_tol or j == 0: # If RT deviation is too large, go to next iteration of base file
                mz_m, int_m = _merge_peaks(base.rawdata[i]["m/z array"], base.rawdata[i]["intensity array"], base_mz_array_single_weighed_start[i], base_intensity_array_single_weighed_start[i], mass_tolerance)
                base.rawdata[i]["intensity array"] = int_m
                base.rawdata[i]["m/z array"] = mz_m
                base.rawdata[i]["total ion current"] = base.rawdata[i]["total ion current"] + base_tic_single_weighed_start[i]
                continue

            # Now that the best matching new spectrum has been found, calculate the average of the two given spectra. 
            # Take the weight given in the arguments to account for multiple spectras being added up.
            # If the weight is 4, it means that the output mzml file consists of (3/4*base + 1/4*new). This is performed for every calculation

            mz_a = base.rawdata[i]["m/z array"]
            int_a = base.rawdata[i]["intensity array"]
            mz_b = new.rawdata[j]["m/z array"]
            int_b = new.rawdata[j]["intensity array"] / len(mzml_filename_list)

            mz_m, int_m = _merge_peaks(mz_a, int_a, mz_b, int_b, mass_tolerance)
            base.rawdata[i]["m/z array"] = mz_m
            base.rawdata[i]["intensity array"] = int_m
            base.rawdata[i]["total ion current"] = float(base.rawdata[i]["total ion current"]) + (float(new.rawdata[j]["total ion current"]) / spectra_weight)
            used_new[j] = True

        # append trailing spectra if new file is longer
        base_end = base_rt[-1] if base_rt else -np.inf
        for j, rt in enumerate(new_rt):
            if used_new[j]:
                continue
            if rt > base_end + curr_rt_tol:
                base.rawdata.append(copy.deepcopy(new.rawdata[j]))

        # recompute derived properties
        base.rt_list = [e["scanList"]["scan"][0]["scan time"] for e in base.rawdata]
        base.tic = [e["total ion current"] for e in base.rawdata]
        base.all_filters = [e["scanList"]["scan"][0]["filter string"] for e in base.rawdata]

        base.all_modes = []
        for f in base.all_filters:
            curr_mode = MS_functions.get_mode_of_spec(f)
            base.all_modes.append(curr_mode)
        base.available_modes = list(sorted(set(m for m in base.all_modes if m != "Unknown")))

        base.rt_range = [min(base.rt_list), max(base.rt_list)]
        base.method_duration = base.rt_range[1] - base.rt_range[0]

        # mz range over all scans
        mz_min, mz_max = None, None
        for e in base.rawdata:
            arr = e.get("m/z array", [])
            if len(arr) == 0:
                continue
            mn, mx = min(arr), max(arr)
            mz_min = mn if mz_min is None or mn < mz_min else mz_min
            mz_max = mx if mz_max is None or mx > mz_max else mz_max
        if mz_min is not None and mz_max is not None:
            base.mz_range = [float(mz_min), float(mz_max)]

        base.aif_background_spectrum = None
        base.ms1_background_spectrum = None

        try:
            base.save_ms_file_log_entry("vINFO:\tMerged file: {}".format(fname))
            base.save_ms_file_log_entry("vINFO:\tScans: {} Duration(s): {}".format(len(base.rawdata), base.method_duration))
        except Exception:
            pass

    return base









