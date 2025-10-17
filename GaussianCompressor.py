#!/usr/bin/env python3

"""
HPLC-MS chromatogram compressor using up to 10 Gaussian peaks.

Goal
-----
Approximate a 1D time series y(t) with a sum of K Gaussians (K ≤ 10),
selecting the smallest K that achieves a target normalized L1 area error.
Only the Gaussian parameters are saved to disk to reduce storage.

Quality metric
--------------
Area error = ∫ |y(t) - \hat{y}(t)| dt
Normalized area error = Area error / ∫ |y(t)| dt

CLI
---
python hplc_ms_gaussian_compressor.py \
  --input chromatogram.csv \
  --max-peaks 10 \
  --target-rel-area 0.02 \
  --output-json params.json \
  --save-reconstruction recon.csv \
  --plot

Input format
------------
CSV or TSV with two columns: time, intensity. Header optional.

Notes
-----
- Uses robust nonlinear least squares (soft L1 loss) with bounds.
- Peak seeds come from scipy.signal.find_peaks on a lightly smoothed trace.
- Amplitudes constrained to be nonnegative, sigmas positive within sensible bounds.
- No baseline term is fitted by default to keep parameter count minimal.
  If your data have strong baseline drift, pre-correct it upstream.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
import time
from typing import List, Tuple

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks, savgol_filter, peak_widths


import numpy as np




@dataclass
class PeakParams:
    A: float
    mu: float
    sigma: float

    def as_list(self) -> List[float]:
        return [float(self.A), float(self.mu), float(self.sigma)]


def _smooth(y):
    n = y.size
    win = max(5, int(0.005*n) | 1)
    poly = 2 if win >= 7 else 1
    return savgol_filter(y, window_length=win, polyorder=poly, mode="interp")

def detect_next_seed_from_residual(t, resid):
    r = _smooth(np.maximum(resid, 0.0))
    peaks, props = find_peaks(r, prominence=max(float(r.max())*0.05, 1e-12))
    if peaks.size == 0:
        return None
    order = np.argsort(props.get("prominences", np.ones_like(peaks)))[::-1]
    idx = int(peaks[order][0])
    w = peak_widths(r, np.array([idx]), rel_height=0.5)[0][0] if peaks.size else 0.0
    dt = float(np.median(np.diff(t)))
    fwhm_to_sigma = 1.0/(2.0*np.sqrt(2.0*np.log(2.0)))
    sigma0 = max((w*dt)*fwhm_to_sigma if w>0 else (t[-1]-t[0])/200.0, dt)
    return PeakParams(A=float(r[idx]), mu=float(t[idx]), sigma=float(sigma0))


def gauss(t: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    # Guard against sigma underflow during intermediate evaluations
    sigma = max(sigma, np.finfo(float).eps)
    z = (t - mu) / sigma
    return A * np.exp(-0.5 * z * z)


def sum_of_gaussians(t: np.ndarray, params: List[PeakParams]) -> np.ndarray:
    y = np.zeros_like(t, dtype=float)
    for p in params:
        y += gauss(t, p.A, p.mu, p.sigma)
    return y


def flatten_params(params: List[PeakParams]) -> np.ndarray:
    vec = []
    for p in params:
        vec.extend([p.A, p.mu, p.sigma])
    return np.asarray(vec, dtype=float)


def unflatten_params(vec: np.ndarray) -> List[PeakParams]:
    assert vec.size % 3 == 0
    out = []
    for i in range(0, vec.size, 3):
        out.append(PeakParams(A=vec[i], mu=vec[i + 1], sigma=vec[i + 2]))
    return out


def area_error(t: np.ndarray, y_true: np.ndarray, y_fit: np.ndarray) -> Tuple[float, float]:
    # Trapezoidal integral of absolute difference, normalized by area of |y_true|
    abs_diff = np.abs(y_true - y_fit)
    total = np.trapz(abs_diff, t)
    denom = np.trapz(np.abs(y_true), t)
    rel = total / denom if denom > 0 else math.inf
    return total, rel


def estimate_noise_sigma(y: np.ndarray) -> float:
    # Median absolute deviation as a robust noise scale estimate
    med = np.median(y)
    mad = np.median(np.abs(y - med))
    # Convert MAD to sigma for Gaussian noise
    return 1.4826 * mad if mad > 0 else max(1e-12, np.std(y) * 0.25)


def seed_peaks(t, y, max_peaks):
    y_s = _smooth(y)
    peaks, props = find_peaks(y_s, prominence=max(float(y_s.max())*0.02, 1e-12))
    if peaks.size == 0:
        idx = int(np.argmax(y_s))
        est_sigma = max((t[-1]-t[0])/100.0, np.finfo(float).eps)
        return [PeakParams(A=float(max(y[idx],0.0)), mu=float(t[idx]), sigma=float(est_sigma))]
    prom = props.get("prominences", np.ones_like(peaks, float))
    order = np.argsort(prom)[::-1]
    peaks = peaks[order][:max_peaks]
    widths = peak_widths(y_s, peaks, rel_height=0.5)[0]
    dt = float(np.median(np.diff(t)))
    f2s = 1.0/(2.0*np.sqrt(2.0*np.log(2.0)))
    seeds = []
    for i, idx in enumerate(peaks):
        A0 = max(float(y[idx]), 0.0)
        mu0 = float(t[idx])
        sigma0 = max((float(widths[i])*dt)*f2s if i < len(widths) and widths[i]>0 else (t[-1]-t[0])/200.0, dt)
        seeds.append(PeakParams(A=A0, mu=mu0, sigma=sigma0))
    return seeds


def _resid_jac_vec(t: np.ndarray, y: np.ndarray, v: np.ndarray):
    k = v.size // 3
    A = v[0::3]
    mu = v[1::3]
    sigma = np.maximum(v[2::3], np.finfo(float).eps)

    T = t[None, :]
    MU = mu[:, None]
    SG = sigma[:, None]

    Z = (T - MU) / SG
    E = np.exp(-0.5 * Z * Z)            # (k,n)
    G = A[:, None] * E                  # (k,n)

    y_hat = G.sum(axis=0)               # (n,)
    r = y_hat - y                       # (n,)

    J = np.empty((t.size, v.size), dtype=float)
    J[:, 0::3] = E.T                    # d/dA
    J[:, 1::3] = (G * (Z / SG)).T       # d/dmu
    J[:, 2::3] = (G * ((Z * Z) / SG)).T # d/dsigma
    return r, J

def fit_k_gaussians(
    t: np.ndarray,
    y: np.ndarray,
    seeds: List[PeakParams],
    k: int,
    min_sigma: float,
    max_sigma: float,
    loss: str = "soft_l1",
    f_scale: float | None = None,
) -> Tuple[List[PeakParams], np.ndarray]:
    seeds_k = sorted(seeds, key=lambda p: p.A, reverse=True)[:k]
    seeds_k = sorted(seeds_k, key=lambda p: p.mu)

    x0 = flatten_params(seeds_k)

    # bounds
    t_min, t_max = float(np.min(t)), float(np.max(t))
    lb = np.tile([0.0, t_min, min_sigma], len(seeds_k)).astype(float)
    ub = np.tile([np.inf, t_max, max_sigma], len(seeds_k)).astype(float)

    if f_scale is None:
        f_scale = estimate_noise_sigma(y)
        if f_scale <= 0:
            f_scale = max(1e-12, float(np.std(y) * 0.1))

    cache = {"x": None, "r": None, "J": None}

    def _fun(v):
        if cache["x"] is not None and np.array_equal(v, cache["x"]):
            return cache["r"]
        r, J = _resid_jac_vec(t, y, v)
        cache["x"], cache["r"], cache["J"] = v.copy(), r, J
        return r

    def _jac(v):
        if cache["x"] is not None and np.array_equal(v, cache["x"]):
            return cache["J"]
        r, J = _resid_jac_vec(t, y, v)
        cache["x"], cache["r"], cache["J"] = v.copy(), r, J
        return J

    res = least_squares(
        _fun, x0, jac=_jac, bounds=(lb, ub),
        loss=loss, f_scale=f_scale, method="trf",
        max_nfev=2000, xtol=1e-8, ftol=1e-8, gtol=1e-8, verbose=0,
    )

    params = unflatten_params(res.x)
    for p in params:
        if p.A < 0:
            p.A = 0.0
    params.sort(key=lambda p: p.mu)
    y_fit = sum_of_gaussians_fast(t, params)
    return params, y_fit


def _peak_curve(t, p):
    inv2s2 = 0.5 / (p.sigma * p.sigma)
    dt = t - p.mu
    return p.A * np.exp(-dt * dt * inv2s2)

def sum_of_gaussians_fast(t, params):
    if not params:
        return np.zeros_like(t, dtype=float)
    y = np.zeros_like(t, dtype=float)
    for p in params:
        y += _peak_curve(t, p)
    return y

def choose_model(t, y, max_peaks, target_rel_area, fit_stride=None):
    start_time = time.time()

    if fit_stride is None:
        # aim for about 2000 points during solves
        fit_stride = max(1, int(len(t) // 2000))
    sl = slice(None, None, fit_stride)
    t_fit = t[sl]
    y_fit_in = y[sl]

    seeds0 = seed_peaks(t_fit, y_fit_in, max_peaks=max_peaks)

    dt = float(np.median(np.diff(t_fit)))
    span = float(t_fit[-1] - t_fit[0]) if t_fit.size > 1 else 1.0
    min_sigma = max(dt, span/2000.0)
    max_sigma = max(span/3.0, min_sigma*2.0)
    noise_scale = estimate_noise_sigma(y_fit_in)

    rng = np.random.default_rng(42)
    params_running, history = [], []
    best_params, best_fit, best_rel = None, None, np.inf

    if seeds0:
        params_running.append(sorted(seeds0, key=lambda p: p.A, reverse=True)[0])

    y_fit = sum_of_gaussians_fast(t_fit, params_running)
    resid = y_fit_in - y_fit

    for k in range(1, max_peaks + 1):
        # detect on thinned grid
        if k > len(params_running):
            new_seed = detect_next_seed_from_residual(t_fit, resid)
            if new_seed is None:
                break
            params_running.append(new_seed)
            y_fit += _peak_curve(t_fit, new_seed)
            resid = y_fit_in - y_fit

        n_starts = 2 if k <= 3 else 1
        base = params_running

        def jitter_once():
            out = []
            for p in base:
                mu_j = float(np.clip(p.mu + rng.normal(0, 0.25*max(dt, p.sigma)), t_fit[0], t_fit[-1]))
                sigma_j = float(np.clip(p.sigma * np.exp(rng.normal(0, 0.25)), min_sigma, max_sigma))
                A_j = max(0.0, p.A * np.exp(rng.normal(0, 0.25)))
                out.append(PeakParams(A=A_j, mu=mu_j, sigma=sigma_j))
            return out

        starts = [base] + [jitter_once() for _ in range(n_starts)]

        best_k, best_fit_k_fitgrid, best_rel_k = None, None, np.inf
        for idx, st in enumerate(starts):
            params_k, y_fit_k_fitgrid = fit_k_gaussians(
                t_fit, y_fit_in, st, k=len(st),
                min_sigma=min_sigma, max_sigma=max_sigma,
                loss="soft_l1", f_scale=noise_scale
            )
            # compute rel area on FULL grid once per candidate
            y_fit_full = sum_of_gaussians_fast(t, params_k)
            _, rel_full = area_error(t, y, y_fit_full)

            if rel_full < best_rel_k:
                best_rel_k, best_k = rel_full, params_k
                best_fit_k_fitgrid = y_fit_k_fitgrid

            if target_rel_area is not None and rel_full <= target_rel_area:
                break
            if idx == 0 and n_starts > 0 and best_rel_k <= best_rel * 1.01:
                break

        params_running = best_k
        y_fit = best_fit_k_fitgrid
        resid = y_fit_in - y_fit

        area_abs, area_rel = area_error(t, y, sum_of_gaussians_fast(t, params_running))
        history.append({"k": len(params_running), "area_abs": area_abs, "area_rel": area_rel})
        if len(history) >= 2 and history[-2]["area_rel"] - area_rel < 1e-4:
            if target_rel_area is None or area_rel <= max(2*estimate_noise_sigma(y)/np.trapz(np.abs(y), t), target_rel_area):
                break
        

        if area_rel < best_rel:
            best_params, best_fit, best_rel = params_running, sum_of_gaussians_fast(t, params_running), area_rel
        if target_rel_area is not None and area_rel <= target_rel_area:
            break

        if len(params_running) < max_peaks:
            new_seed = detect_next_seed_from_residual(t_fit, resid)
            if new_seed is not None:
                params_running = params_running + [new_seed]
                y_fit += _peak_curve(t_fit, new_seed)
                resid = y_fit_in - y_fit

    if best_params is None:
        best_params, best_fit = fit_k_gaussians(t_fit, y_fit_in, seeds0[:1], 1, min_sigma, max_sigma)
        best_rel = area_error(t, y, sum_of_gaussians_fast(t, best_params))[1]

    info = {"area_rel_best": best_rel, "history": history,
            "min_sigma": min_sigma, "max_sigma": max_sigma,
            "target_rel_area": target_rel_area}
    return best_params, sum_of_gaussians_fast(t, best_params), info



def load_xy(path: str) -> Tuple[np.ndarray, np.ndarray]:
    # Try CSV then fallback to generic whitespace
    try:
        data = np.genfromtxt(path, delimiter=",", comments="#", dtype=float)
        if data.ndim == 1:
            data = data.reshape(-1, 2)
        if data.shape[1] < 2 or np.any(np.isnan(data)):
            raise ValueError("invalid CSV parse")
    except Exception:
        data = np.genfromtxt(path, delimiter=None, comments="#", dtype=float)
        if data.ndim == 1:
            data = data.reshape(-1, 2)
    if data.shape[1] < 2:
        raise ValueError("input must have at least two columns: time, intensity")
    t = np.asarray(data[:, 0], dtype=float)
    y = np.asarray(data[:, 1], dtype=float)
    # Ensure strictly increasing time
    order = np.argsort(t)
    t = t[order]
    y = y[order]
    return t, y


def save_params_json(path: str, params: List[PeakParams]) -> None:
    payload = {
        "model": "sum_of_gaussians",
        "version": 1,
        "params": [
            {"A": float(p.A), "mu": float(p.mu), "sigma": float(p.sigma)} for p in params
        ],
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def save_reconstruction_csv(path: str, t: np.ndarray, y_true: np.ndarray, y_fit: np.ndarray) -> None:
    arr = np.column_stack([t, y_true, y_fit])
    header = "time,intensity,reconstruction"
    np.savetxt(path, arr, delimiter=",", header=header, comments="", fmt="%.10g")


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Compress HPLC-MS chromatograms with up to 10 Gaussian peaks.")
    p.add_argument("--input", required=True, help="Path to CSV/TSV with columns: time,intensity")
    p.add_argument("--max-peaks", type=int, default=10, help="Maximum number of Gaussians")
    p.add_argument(
        "--target-rel-area",
        type=float,
        default=0.02,
        help="Stop at smallest K with normalized area error ≤ this value. Use 0 to always use max K.",
    )
    p.add_argument("--output-json", required=True, help="Where to write Gaussian parameters as JSON")
    p.add_argument("--save-reconstruction", default=None, help="Optional CSV with time,intensity,reconstruction")
    p.add_argument("--plot", action="store_true", help="Show a quick plot for visual check")

    args = p.parse_args(argv)

    t, y = load_xy(args.input)
    # Clip tiny negatives which may occur with noisy MS baselines
    y = np.maximum(y, 0.0)

    target = None if args.target_rel_area and args.target_rel_area <= 0 else args.target_rel_area

    params, y_fit, info = choose_model(t, y, max_peaks=args.max_peaks, target_rel_area=target)

    save_params_json(args.output_json, params)

    if args.save_reconstruction:
        save_reconstruction_csv(args.save_reconstruction, t, y, y_fit)

    # Console report
    print("Chosen K:", len(params))
    print("Normalized area error:", f"{info['area_rel_best']:.6g}")
    print("K vs rel area:")
    for row in info["history"]:
        print(f"  K={row['k']:2d}  rel_area={row['area_rel']:.6g}")
    print("Parameter list (A, mu, sigma):")
    for i, p_i in enumerate(params, 1):
        print(f"  {i:2d}: {p_i.A:.6g}, {p_i.mu:.6g}, {p_i.sigma:.6g}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt

            plt.figure()
            plt.plot(t, y, label="data")
            plt.plot(t, y_fit, label="fit", linewidth=2)
            # Also draw individual components for inspection
            for i, p_i in enumerate(params, 1):
                plt.plot(t, gauss(t, p_i.A, p_i.mu, p_i.sigma), linestyle=":", label=f"g{i}")
            plt.xlabel("time")
            plt.ylabel("intensity")
            plt.legend()
            plt.tight_layout()
            plt.show()
        except Exception as e:
            print("Plotting failed:", e, file=sys.stderr)

    return 0


def merge_peaks_by_resolution(params: List[PeakParams], rs_thresh: float = 1.0) -> List[PeakParams]:
    """
    Merge adjacent Gaussians whose chromatographic resolution R_s < rs_thresh.
    Preserves total area and matches first two moments to yield one Gaussian per cluster.
    Returns a new list of PeakParams sorted by mu.
    """
    FWHM = 2.354820045
    if not params:
        return []

    ps = sorted(params, key=lambda p: p.mu)

    def _rs(p1: PeakParams, p2: PeakParams) -> float:
        w1 = FWHM * max(p1.sigma, np.finfo(float).eps)
        w2 = FWHM * max(p2.sigma, np.finfo(float).eps)
        return 2.0 * abs(p2.mu - p1.mu) / (w1 + w2)

    # build clusters by adjacency
    clusters, cur = [], [ps[0]]
    for p in ps[1:]:
        if _rs(cur[-1], p) < rs_thresh:
            cur.append(p)
        else:
            clusters.append(cur)
            cur = [p]
    clusters.append(cur)

    out: List[PeakParams] = []
    for cl in clusters:
        if len(cl) == 1:
            out.append(cl[0])
            continue

        # area weights
        areas = np.array([p.A * max(p.sigma, np.finfo(float).eps) * np.sqrt(2*np.pi) for p in cl], dtype=float)
        Atot = float(areas.sum())
        if Atot <= 0:
            out.append(max(cl, key=lambda p: p.A))
            continue

        mus = np.array([p.mu for p in cl], dtype=float)
        sig2 = np.array([p.sigma**2 for p in cl], dtype=float)
        w = areas / Atot

        mu_bar = float(np.sum(w * mus))
        var_bar = float(np.sum(w * (sig2 + (mus - mu_bar)**2)))
        sigma_bar = float(np.sqrt(max(var_bar, np.finfo(float).eps)))
        A_eff = Atot / (sigma_bar * np.sqrt(2*np.pi))  # area-preserving

        out.append(PeakParams(A=float(A_eff), mu=mu_bar, sigma=sigma_bar))

    return out


def get_peaks_in_chromatogram(timevals, intensityvals, max_peaks=50, target_rel_area=2e-7):
    t = np.asarray(timevals, dtype=float)
    y = np.asarray(intensityvals, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or t.size != y.size or t.size < 5:
        raise ValueError("input time and intensity must be 1D arrays of same length >=5")
    

    # Ensure strictly increasing time
    order = np.argsort(t)
    t = t[order]
    y = y[order]
    # Clip tiny negatives which may occur with noisy MS baselines
    y = np.maximum(y, 0.0)

    target = None if target_rel_area and target_rel_area <= 0 else target_rel_area
    max_peaks = max(1, min(300, int(max_peaks)))

    print("Fitting up to", max_peaks, "peaks with target rel area", target)
    params, y_fit, info = choose_model(t, y, max_peaks=max_peaks, target_rel_area=target)
    params = merge_peaks_by_resolution(params)

    params_list = []
    for p in params:
        params_list.append((p.A, p.mu, p.sigma))
    return params_list # params contains list of PeakParams objects with A, mu, sigma attributes



if __name__ == "__main__":
    # Direct usage without CLI args
    # Set your input path
    path = "test_chromatogram.txt"
    t, y = load_xy(path)

    # Fit up to 10 peaks and stop when target normalized area error is reached
    params, y_fit, info = choose_model(t, np.maximum(y, 0.0), max_peaks=10, target_rel_area=0.00000002)

    print("Fitted parameters:")
    for i, p in enumerate(params, 1):
        print(f"  Peak {i}: A={p.A:.6g}, mu={p.mu:.6g}, sigma={p.sigma:.6g}")
    print()

    # Combine individual params if there are multiple gauss peaks for a single retention time:
    params = merge_peaks_by_resolution(params)
    print("Fitted parameters:")
    for i, p in enumerate(params, 1):
        print(f"  Peak {i}: A={p.A:.6g}, mu={p.mu:.6g}, sigma={p.sigma:.6g}")
    print()

    # Remove gaussian peaks that are too small or too wide
    print(len(params), "peaks before filtering")
    abs_duration_of_timeseries = t[-1] - t[0]
    #params = [p for p in params if p.A > 5e6 and p.sigma < abs_duration_of_timeseries/4.0]
    print(len(params), "peaks after filtering")
    
    # Save compact outputs
    save_params_json("params.json", params)
    y_fit = sum_of_gaussians(t, params)

    # Report
    print("Chosen K:", len(params))
    print("Normalized area error:", f"{info['area_rel_best']:.6g}")
    for row in info["history"]:
        print(f"  K={row['k']:2d}  rel_area={row['area_rel']:.6g}")

    # Optional plot
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(t, y, label="data")
        plt.plot(t, y_fit, label="fit", linewidth=2)
        for i, p_i in enumerate(params, 1):
            plt.plot(t, gauss(t, p_i.A, p_i.mu, p_i.sigma), linestyle=":", label=f"g{i}")
        plt.xlabel("time")
        plt.ylabel("intensity")
        #plt.legend()
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print("Plotting failed:", e)
