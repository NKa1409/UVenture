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


def calculate_averaged_spectrum(mzml_file, indices_of_spectra, mass_deviation=15):
    def _ppm_tol(m):
        return (mass_deviation * float(m)) / 1e6
    # akzeptiere entweder eine MS_File Instanz oder einen Dateinamen
    if isinstance(mzml_file, class_MS_file.MS_File):
        ms = mzml_file
    else:
        ms = class_MS_file.MS_File(mzml_file)
    indices = indices_of_spectra
    all_mz = []
    all_int = []
    for idx in indices:
        mz_arr = np.asarray(ms.rawdata[idx].get("m/z array", []), dtype=float)
        int_arr = np.asarray(ms.rawdata[idx].get("intensity array", []), dtype=float)
        if mz_arr.size == 0:
            continue
        all_mz.append(mz_arr)
        all_int.append(int_arr)
    if not all_mz:
        return [], []
    all_mz = np.concatenate(all_mz)
    all_int = np.concatenate(all_int)
    order = np.argsort(all_mz)
    all_mz = all_mz[order]
    all_int = all_int[order]

    merged_mz = []
    merged_int = []

    cur_mz = float(all_mz[0])
    cur_int = float(all_int[0])
    for k in range(1, all_mz.size):
        mz_k = float(all_mz[k])
        int_k = float(all_int[k])
        if abs(mz_k - cur_mz) <= _ppm_tol(cur_mz):
            total_int = cur_int + int_k
            if total_int > 0:
                cur_mz = (cur_mz * cur_int + mz_k * int_k) / total_int
            cur_int = total_int
        else:
            merged_mz.append(cur_mz)
            merged_int.append(cur_int)
            cur_mz = mz_k
            cur_int = int_k
    merged_mz.append(cur_mz)
    merged_int.append(cur_int)
    n_bg = float(len(indices))
    merged_int = [i / n_bg for i in merged_int]
    summarized_dict = MS_functions.summarize_mass_intensity_dict(dictio=dict(zip(merged_mz, merged_int)), deviation=11, debug_output=True)
    return list(summarized_dict.keys()), list(summarized_dict.values())



def remove_background(
    mzml_file,
    background_signals=None,
    mass_deviation=5,
    remove_completely=True,
    death_time=10,
    save_mzml_filepath="",
    max_signals=10
):
    # mass_deviation given in ppm
    # background_signals:
    #   - list: [m1, m2, ...]  -> gleiche Liste für alle Modi
    #   - dict: {mode: [m1, m2, ...]} -> pro Modus (z.B. "Full scan", "MS/MS", "AIF")
    #   - None: Hintergrund wird aus frühen Scans (rt <= death_time) pro Modus abgeleitet
    # Wenn background_signals None ist, wird remove_completely intern auf False gesetzt.
    # remove_completely:
    #   - True: Peaks an den Hintergrundmassen werden komplett entfernt.
    #   - False: gemittelter Hintergrund der ersten (max. 10) Scans pro Modus wird subtrahiert.
    # death_time in Sekunden.
    # save_mzml_filepath: wenn nicht leer, wird die modifizierte Datei gespeichert.

    # akzeptiere entweder eine MS_File Instanz oder einen Dateinamen
    if isinstance(mzml_file, class_MS_file.MS_File):
        ms = mzml_file
    else:
        ms = class_MS_file.MS_File(mzml_file)

    def _ppm_tol(m):
        return (mass_deviation * float(m)) / 1e6

    # RT Liste
    try:
        rt_list = list(ms.rt_list)
    except AttributeError:
        rt_list = [e["scanList"]["scan"][0]["scan time"] for e in ms.rawdata]

    # Mode-Liste (bevorzugt) oder Fallback auf Filterstring
    try:
        mode_list = list(ms.all_modes)
        if len(mode_list) != len(ms.rawdata):
            raise ValueError
    except Exception:
        mode_list = [e["scanList"]["scan"][0]["filter string"] for e in ms.rawdata]

    # frühe Scans als Hintergrund
    bg_scan_indices = [i for i, rt in enumerate(rt_list) if rt <= death_time]
    if not bg_scan_indices:
        if save_mzml_filepath:
            if hasattr(ms, "save_to_mzml_file"):
                ms.save_to_mzml_file(save_mzml_filepath)
            else:
                raise AttributeError(
                    "MS_File object does not implement 'save_to_mzml_file'. Adapt the saving logic in remove_background()."
                )
        return ms

    # max. 10 Hintergrundscans
    bg_scan_indices = bg_scan_indices[:10]

    # Hintergrundscans nach Modus gruppieren
    from collections import defaultdict
    mode_to_bg_indices = defaultdict(list)
    for idx in bg_scan_indices:
        m = mode_list[idx]
        mode_to_bg_indices[m].append(idx)

    # pro Modus Hintergrundmassen/-intensitäten
    bg_mass_int_dict_by_mode = {}

    for m, idxs in mode_to_bg_indices.items():
        bg_mz, bg_int = calculate_averaged_spectrum(ms, idxs, mass_deviation=mass_deviation)
        bg_mass_int_dict_by_mode[m] = dict(zip(list(bg_mz), list(bg_int)))
        bg_mass_int_dict_by_mode[m] = {m:i for m, i in bg_mass_int_dict_by_mode[m].items() if i >= 200}
        try:
            bg_mass_int_dict_by_mode[m] = dict(sorted(bg_mass_int_dict_by_mode[m].items(), key=lambda x: x[1], reverse=True)[:max_signals])
        except:
            bg_mass_int_dict_by_mode[m] = bg_mass_int_dict_by_mode[m]
    
    if background_signals is not None:
        for m in list(bg_mass_int_dict_by_mode.keys()):
            bg_mass_int_dict_by_mode[m] = {ma:i for ma,i in bg_mass_int_dict_by_mode[m].items() if any(abs(ma - target) <= target * mass_deviation * 1e-6 for target in background_signals)}
            for masse in background_signals:
                try:
                    bg_mass_int_dict_by_mode[m][masse] += 1
                except:
                    bg_mass_int_dict_by_mode[m][masse] = 1
    print(bg_mass_int_dict_by_mode)


    bg_mz_by_mode = {}
    bg_int_by_mode = {}
    for mode, bg_dict in bg_mass_int_dict_by_mode.items():
        bg_mz = np.fromiter(bg_dict.keys(), dtype=float)
        bg_int = np.fromiter(bg_dict.values(), dtype=float)
        order = np.argsort(bg_mz)
        bg_mz_by_mode[mode] = bg_mz[order]
        bg_int_by_mode[mode] = bg_int[order]

    # Hintergrundentfernung/Subtraktion pro Spektrum, nach Modus getrennt
    for i in range(len(ms.rawdata)):
        spec = ms.rawdata[i]
        m = mode_list[i]
        progress = i / len(ms.rawdata)
        filled = int(20 * progress)
        bar = "#" * filled + "-" * (20 - filled)
        print(f"\rProcessing: [{bar}] {i}/{len(ms.rawdata)} ({progress * 100:.1f}%)", end="", flush=True)
        mz_arr = np.asarray(spec.get("m/z array", []), dtype=float)
        int_arr = np.asarray(spec.get("intensity array", []), dtype=float)
        if mz_arr.size == 0:
            continue

        tic = float(spec.get("total ion current", 0.0))


        if not m in bg_mz_by_mode:
            spec["m/z array"] = mz_arr.tolist()
            spec["intensity array"] = int_arr.tolist()
            spec["total ion current"] = tic
            continue

        if remove_completely:
            keep_mask = np.ones(mz_arr.shape, dtype=bool)
            if background_signals is not None:
                masses = background_signals
            else:
                masses = bg_mass_int_dict_by_mode[m]
            for m_bg in masses:
                tol = _ppm_tol(m_bg)
                mask = np.abs(mz_arr - m_bg) <= tol
                if not np.any(mask):
                    continue
                tic -= float(np.sum(int_arr[mask]))
                keep_mask &= ~mask
            mz_arr = mz_arr[keep_mask]
            int_arr = int_arr[keep_mask]


        else:
            tol = mz_arr * mass_deviation * 1e-6
            # index range in bg_mz for each mz_arr[i]
            left = np.searchsorted(bg_mz, mz_arr - tol, side="left")
            right = np.searchsorted(bg_mz, mz_arr + tol, side="right")

            # subtract background
            bg_sub = np.zeros_like(int_arr)

            for i in range(mz_arr.size):
                lo = left[i]
                hi = right[i]
                if lo < hi:
                    # choose the closest background mass in the window
                    window = bg_mz[lo:hi]
                    j = lo + np.argmin(np.abs(window - mz_arr[i]))
                    bg_sub[i] = bg_int[j]

            int_corr = int_arr - bg_sub

            spec["m/z array"] = mz_arr.tolist()
            spec["intensity array"] = int_corr.tolist()
            # if you want TIC after subtraction:
            spec["total ion current"] = float(int_corr.sum())
            spec_massses = spec["m/z array"]
            spec_intensities = spec["intensity array"]
            spec_massint_dict = dict(zip(list(spec_massses), list(spec_intensities)))
            corrected = {
                ma: I - next(
                    (bgI for mbg, bgI in bg_mass_int_dict_by_mode[m].items()
                    if abs(ma - mbg) <= ma * mass_deviation * 1e-6),
                    0
                )
                for ma, I in spec_massint_dict.items()
            }
            mz_arr = list(corrected.keys())
            int_arr = list(corrected.values())

        spec["m/z array"] = mz_arr
        spec["intensity array"] = int_arr
        spec["total ion current"] = tic

    # TIC und mz_range neu berechnen
    try:
        ms.tic = [e["total ion current"] for e in ms.rawdata]
    except Exception:
        pass

    mz_min = None
    mz_max = None
    for e in ms.rawdata:
        arr = np.asarray(e.get("m/z array", []), dtype=float)
        if arr.size == 0:
            continue
        mn = float(arr.min())
        mx = float(arr.max())
        mz_min = mn if mz_min is None or mn < mz_min else mz_min
        mz_max = mx if mz_max is None or mx > mz_max else mz_max
    if mz_min is not None and mz_max is not None:
        ms.mz_range = [float(mz_min), float(mz_max)]

    # evtl. gecachte Hintergrundspektren invalidieren
    if hasattr(ms, "aif_background_spectrum"):
        ms.aif_background_spectrum = None
    if hasattr(ms, "ms1_background_spectrum"):
        ms.ms1_background_spectrum = None

    # optionales Speichern
    if save_mzml_filepath:
        if hasattr(ms, "save_to_mzml_file"):
            ms.save_to_mzml_file(save_mzml_filepath)
        else:
            raise AttributeError(
                "MS_File object does not implement 'save_to_mzml_file'. Adapt the saving logic in remove_background()."
            )

    return ms
