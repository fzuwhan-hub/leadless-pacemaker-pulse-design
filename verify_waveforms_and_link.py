import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2

B_REF = 10.0e6
DT_SYNTH = 0.05e-9
T_MAX = 16.0e-6
N_FFT = 2**20
DT_BER = 25.0e-9
T_GUARD = 2.925e-6


def causal_window(t_arr, tw, delta):
    w = np.ones_like(t_arr)
    mask = (t_arr > tw - delta) & (t_arr <= tw)
    w[mask] = 0.5 * (1.0 + np.cos(np.pi * (t_arr[mask] - tw + delta) / delta) ** 3)
    w[t_arr > tw] = 0.0
    return w


def synthesize_candidate(family, shape_val, t_arr, freqs):
    te = 35.0e-9 / np.log(9.0)
    if family == "DSBP":
        tp = 1.0e-6
        g = np.where(
            t_arr < shape_val * tp,
            1.0 - np.exp(-t_arr / te),
            (1.0 - np.exp(-shape_val * tp / te)) * np.exp(-(t_arr - shape_val * tp) / te),
        )
        eg = np.sum(g**2) * DT_SYNTH
        tc = float(t_arr[np.searchsorted(np.cumsum(g**2) * DT_SYNTH / eg, 1.0 - 1.0e-6)])
        tw_active = tc
    elif family == "CR-RC":
        tau1 = 5.0e-9 + 13.0 * 25.0e-9 / 79.0
        tau2 = shape_val * tau1
        g0 = np.exp(-t_arr / tau1) - np.exp(-t_arr / tau2)
        e0 = np.sum(g0**2) * DT_SYNTH
        tc = float(t_arr[np.searchsorted(np.cumsum(g0**2) * DT_SYNTH / e0, 1.0 - 1.0e-6)])
        tw_active = 1.1 * tc
        g = g0 * causal_window(t_arr, tw_active, 0.1 * tc)
        eg = np.sum(g**2) * DT_SYNTH
    else:
        f0 = 9.0e6
        g0 = np.exp(-np.pi * f0 * t_arr / shape_val) * np.sin(2.0 * np.pi * f0 * t_arr)
        e0 = np.sum(g0**2) * DT_SYNTH
        tc = float(t_arr[np.searchsorted(np.cumsum(g0**2) * DT_SYNTH / e0, 1.0 - 1.0e-6)])
        tw_active = 1.1 * tc
        g = g0 * causal_window(t_arr, tw_active, 0.1 * tc)
        eg = np.sum(g**2) * DT_SYNTH

    t99 = float(t_arr[np.searchsorted(np.cumsum(g**2) * DT_SYNTH / eg, 0.99)])
    t_sep = 1.2 * t99
    shift = int(round(t_sep / DT_SYNTH))
    q = g.copy()
    q[shift:] -= g[:-shift]

    qf = np.fft.rfft(q, n=N_FFT) * DT_SYNTH
    psd = np.abs(qf) ** 2
    b99_unscaled = float(np.interp(0.99, np.cumsum(psd) / np.sum(psd), freqs))
    alpha = B_REF / b99_unscaled

    mean_f2 = float(np.sum((freqs**2) * psd) / np.sum(psd))
    cei = float(np.sqrt((alpha**2 * mean_f2) / (B_REF**2)))

    eq = float(np.sum(q**2) * DT_SYNTH)
    theta = float((eq / (np.max(np.abs(q)) ** 2)) / alpha)

    t_supp = (tw_active + t_sep) / alpha
    t_rec = float(np.ceil(t_supp / DT_BER) * DT_BER)

    return {
        "alpha": alpha,
        "cei": cei,
        "theta_us": theta * 1e6,
        "Trec_us": t_rec * 1e6,
        "f_osc_MHz": (alpha * 9.0) if family == "RLC" else np.nan,
    }


def ber_2ppm(rho_db, n_samples):
    rho = 10.0 ** (rho_db / 10.0)
    nc = 2.0 * rho
    sd = np.sqrt(2.0 * n_samples)
    lo = max(0.0, n_samples - 10.0 * sd)
    hi = n_samples + 14.0 * sd
    val, _ = quad(lambda y: ncx2.cdf(y, df=n_samples, nc=nc) * chi2.pdf(y, df=n_samples), lo, hi, limit=200)
    return float(val)


if __name__ == "__main__":
    opt = pd.read_csv("optimization_data.csv")
    t_arr = np.arange(0.0, T_MAX, DT_SYNTH)
    freqs = np.fft.rfftfreq(N_FFT, d=DT_SYNTH)

    synth = [synthesize_candidate(r["family"], r["shape_value"], t_arr, freqs) for _, r in opt.iterrows()]
    ref_idx = opt.index[opt["candidate_id"] == "P2"][0]
    cei_rel = np.array([s["cei"] for s in synth]) / synth[ref_idx]["cei"]

    rel_err = float(np.max(np.abs(cei_rel - opt["C_relative_to_DS_ref"].values) / opt["C_relative_to_DS_ref"].values))
    assert rel_err < 0.005
    print(f"Max relative C_EI difference across 131 records: {rel_err:.2%}\n")

    reps = [
        ("P2", "DS-ref."),
        ("P6", "CR-RC Rep."),
        ("D50", "Knee Point (D50)"),
        ("P16", "RLC-mid"),
        ("P13", "RLC-end"),
    ]
    print("Representative waveforms:")
    for cid, label in reps:
        idx = opt.index[opt["candidate_id"] == cid][0]
        r = opt.iloc[idx]
        s = synth[idx]
        extra = f", f_osc={s['f_osc_MHz']:.2f} MHz" if not np.isnan(s["f_osc_MHz"]) else ""
        print(
            f"  {label:16s} ({cid:3s}): {r['shape_parameter']}={r['shape_value']:.4f}, "
            f"alpha={s['alpha']:.4f}, J={r['J_relative_to_DS_ref']:.2f}, "
            f"C={r['C_relative_to_DS_ref']:.2f}, R_tail={r['R_tail_percent']:.3f}%, "
            f"T_rec={s['Trec_us']:.3f} us{extra}"
        )

    print("\nAnalytical received Eb/N0 thresholds (N = 583):")
    for ber_target in [1e-2, 1e-3, 1e-4]:
        thr = brentq(lambda db: ber_2ppm(db, 583) - ber_target, 5.0, 25.0)
        print(f"  BER = {ber_target:.0e}: {thr:.2f} dB")

    gamma_583 = brentq(lambda db: ber_2ppm(db, 583) - 1e-4, 5.0, 25.0)
    common_ber = [("P2", "DS-ref.", 17.19), ("P6", "CR-RC Rep.", 16.89),
                  ("P16", "RLC-mid", 15.61), ("P13", "RLC-end", 15.30)]

    print("\nWaveform-specific timing and rate-energy benchmarks (BER = 1e-4):")
    for cid, label, g_comm in common_ber:
        idx = opt.index[opt["candidate_id"] == cid][0]
        t_rec_us = synth[idx]["Trec_us"]
        t_slot_us = t_rec_us + T_GUARD * 1e6
        n_i = int(round(t_slot_us / (DT_BER * 1e6)))
        t_s_us = 2.0 * t_slot_us
        rate = 1000.0 / t_s_us
        g_ni = brentq(lambda db: ber_2ppm(db, n_i) - 1e-4, 5.0, 25.0)
        g_own = g_comm + g_ni - gamma_583
        print(
            f"  {label:11s}: T_rec={t_rec_us:.3f} us, N_i={n_i}, T_s={t_s_us:5.2f} us, "
            f"Rate={rate:5.1f} kbit/s, Common={g_comm:.2f} dB, Specific={g_own:.2f} dB"
        )

    idx_ds = opt.index[opt["candidate_id"] == "P2"][0]
    idx_rlc = opt.index[opt["candidate_id"] == "P13"][0]
    ratio_th = synth[idx_rlc]["theta_us"] / synth[idx_ds]["theta_us"]
    vpk_db = (15.30 - 17.19) - 10.0 * np.log10(ratio_th)
    print(
        f"\nEqual peak-amplitude comparison (BER = 1e-4):\n"
        f"  Theta_DS-ref={synth[idx_ds]['theta_us']:.4f} us, "
        f"Theta_RLC-end={synth[idx_rlc]['theta_us']:.4f} us, "
        f"ratio={ratio_th:.4f}, Vpk_RLC/Vpk_DS={vpk_db:+.2f} dB"
    )
