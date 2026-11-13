# Waveform Optimization Data

| File | Contents |
| --- | --- |
| [candidate_objectives.csv](candidate_objectives.csv) | Three objective values for each of the 131 candidate records. `shown_in_Fig5` and `shown_in_Fig6` identify the plotted records. |
| [parameter_shape_mapping.csv](parameter_shape_mapping.csv) | The 10,041 nominal parameter combinations. `shape_id` joins to `candidate_id` in the objective table. |
| [verify_pareto.py](verify_pareto.py) | Recomputes nondominance and the nominal Knee Point from the objective values and checks the parameter-to-shape mapping. |

The parameter mapping groups equivalent scales under the windowing and bandwidth normalization rules described in the manuscript. It contains 41 DSBP records indexed by `d`, 80 CR–RC records indexed by `r`, and 10 RLC records indexed by `Q`. The final CR–RC pulses at `r = 0.25` (C01) and `r = 4` (C80) differ only in overall polarity. They share the objective vector from C80, so the table contains 130 distinct objective vectors. For the RLC combinations, `original_frequency_grid_index` identifies one of 360 nominal frequency positions. `f0_MHz` is blank for these entries.

Figure 5 shows the family trajectories at fixed unscaled parameters. DSBP uses `Tp_us = 1`, CR–RC uses `tau1_ns = 9.113924`, and RLC uses `f0_MHz = 9.0`. C01 and C80 coincide in Fig. 5(b). Figure 6 shows eight nondominated DSBP representatives, six nondominated RLC representatives, and eight dominated CR–RC controls.

Pareto screening maximizes `J_relative_to_DS_end` and minimizes `C_relative_to_DS_end` and `R_tail_percent`. The spectral matching and edge intensity columns use the common DS-end reference. The Knee Point minimizes the Euclidean distance to the ideal point among nondominated candidates after equal-weight min–max normalization over all 131 records.

Run the verification with Python 3 from the repository root:

```bash
python waveform-optimization/verify_pareto.py
```

The script uses the Python standard library. It reports 43 nondominated shapes (34 DSBP and 9 RLC) and identifies P16 as the nominal Knee Point.
