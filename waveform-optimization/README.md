# Waveform Optimization Data

| File | Contents |
| --- | --- |
| [candidate_objectives.csv](candidate_objectives.csv) | Three objective values for each of the 131 shapes. `shown_in_Fig5` and `shown_in_Fig6` identify the plotted records. |
| [parameter_shape_mapping.csv](parameter_shape_mapping.csv) | The 10,041 nominal parameter combinations. `shape_id` joins to `candidate_id` in the objective table. |
| [verify_pareto.py](verify_pareto.py) | Recomputes nondominance and the nominal Knee Point from the objective values and checks the parameter-to-shape mapping. |

The parameter mapping groups equivalent scales under the windowing and bandwidth normalization rules described in the manuscript. It contains 41 DSBP shapes indexed by `d`, 80 CR–RC shapes indexed by `r`, and 10 RLC shapes indexed by `Q`. For the RLC combinations, `original_frequency_grid_index` identifies one of 360 nominal frequency positions. `f0_MHz` is blank for these entries.

Figure 5 shows the family trajectories at fixed unscaled parameters: `Tp_us = 1` for DSBP, `tau1_ns = 9.113924` for CR–RC, and `f0_MHz = 9.0` for RLC. Figure 6 shows eight nondominated DSBP representatives, six nondominated RLC representatives, and eight dominated CR–RC controls.

Pareto screening maximizes `J_relative_to_DS_end` and minimizes `C_relative_to_DS_end` and `R_tail_percent`. The spectral matching and edge intensity columns use the common DS-end reference. The Knee Point minimizes the Euclidean distance to the ideal point among nondominated shapes after equal-weight min–max normalization over all 131 shapes.

Run the verification with Python 3 from the repository root:

```bash
python waveform-optimization/verify_pareto.py
```

The script uses the Python standard library. It reports 43 nondominated shapes (34 DSBP and 9 RLC) and identifies P16 as the nominal Knee Point.
