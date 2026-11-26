import numpy as np
import pandas as pd


def check_grid_grouping(opt_df):
    d_grid = np.linspace(0.10, 0.50, 41)
    tau1_grid = np.linspace(5.0, 30.0, 80)
    r_grid = np.linspace(0.25, 4.00, 80)
    q_grid = np.linspace(1.00, 5.00, 10)
    n_nominal = len(d_grid) + len(tau1_grid) * len(r_grid) + len(q_grid) * 360

    ds_rows = opt_df[opt_df["family"] == "DSBP"]
    cr_rows = opt_df[opt_df["family"] == "CR-RC"]
    rl_rows = opt_df[opt_df["family"] == "RLC"]

    assert n_nominal == 10041
    assert not np.any(np.isclose(r_grid, 1.0))
    assert len(ds_rows) == 41 and np.allclose(ds_rows["shape_value"].values, d_grid, atol=1e-5)
    assert len(cr_rows) == 80 and np.allclose(cr_rows["shape_value"].values, r_grid, atol=1e-5)
    assert len(rl_rows) == 10 and np.allclose(rl_rows["shape_value"].values, q_grid, atol=1e-5)
    print(f"Nominal combinations: {n_nominal} -> {len(opt_df)} scale-grouped records (41 DSBP, 80 CR-RC, 10 RLC)")


def check_pareto_and_knee(opt_df):
    x = np.column_stack([
        -opt_df["J_relative_to_DS_ref"].values,
        opt_df["C_relative_to_DS_ref"].values,
        opt_df["R_tail_percent"].values,
    ])

    unique_coords = np.unique(np.round(x, 6), axis=0)
    c01 = x[opt_df["candidate_id"] == "C01"][0]
    c80 = x[opt_df["candidate_id"] == "C80"][0]
    assert np.allclose(c01, c80)

    n = len(opt_df)
    nondom = np.ones(n, dtype=bool)
    for i in range(n):
        le_all = np.all(x <= x[i] + 1e-12, axis=1)
        lt_any = np.any(x < x[i] - 1e-12, axis=1)
        if np.any(le_all & lt_any):
            nondom[i] = False

    assert np.array_equal(nondom, (opt_df["pareto_status"] == "YES").values)

    z = (x - x.min(axis=0)) / (x.max(axis=0) - x.min(axis=0))
    dist = np.sqrt(np.sum(z**2, axis=1))
    p_idx = np.where(nondom)[0]
    knee_idx = p_idx[np.argmin(dist[p_idx])]
    knee = opt_df.iloc[knee_idx]
    assert opt_df[opt_df["nominal_knee"] == "YES"]["candidate_id"].tolist() == [knee["candidate_id"]]

    counts = opt_df[nondom]["family"].value_counts().to_dict()
    print(f"Distinct objective vectors: {len(unique_coords)}")
    print(f"Pareto-optimal set: { int(nondom.sum()) } records {counts}")
    print(
        f"Nominal Knee Point: {knee['candidate_id']} ({knee['family']}, "
        f"{knee['shape_parameter']}={knee['shape_value']:.2f}, "
        f"J={knee['J_relative_to_DS_ref']:.2f}, C={knee['C_relative_to_DS_ref']:.2f}, "
        f"R_tail={knee['R_tail_percent']:.3f}%, D={dist[knee_idx]:.4f})"
    )


def check_fig3_tables():
    amp = pd.read_csv("amplitude_response.csv")
    res = pd.read_csv("reconstruction_residuals.csv")
    t0 = res[np.isclose(res["time_s"], 0.0)].sort_values("frequency_Hz")
    t54 = res[np.isclose(res["time_s"], 0.54)].sort_values("frequency_Hz")
    assert np.allclose(amp["measured_diastole_dB"].values, t0["measured_dB"].values)
    assert np.allclose(amp["measured_systole_dB"].values, t54["measured_dB"].values)
    rms_db = float(np.sqrt(np.mean(res["absolute_residual_dB"].values ** 2)))
    print(f"Fig. 3 reconstruction RMS residual: {rms_db:.3f} dB")


if __name__ == "__main__":
    df = pd.read_csv("optimization_data.csv")
    check_grid_grouping(df)
    check_pareto_and_knee(df)
    check_fig3_tables()
