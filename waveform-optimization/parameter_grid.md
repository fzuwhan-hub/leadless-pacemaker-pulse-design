# Waveform parameter grid

[parameter_shape_mapping.csv](parameter_shape_mapping.csv) links each combination in the nominal scan (`raw_id`) to a normalized waveform candidate (`shape_id`).

| Family | Shape grid | Template scale | Scan rows | Candidates |
| --- | --- | --- | ---: | ---: |
| DSBP | 41 values of $d$ from 0.10 to 0.50 | $T_p=1$ µs | 41 | 41 |
| CR–RC | 80 values of $r$ from 0.25 to 4 | 80 values of $\tau_1$ from 5 to 30 ns | 6,400 | 80 |
| RLC | 10 values of $Q$ from 1 to 5 | 360 frequency positions | 3,600 | 10 |
| Total | | | 10,041 | 131 |

The grids for $d$, $r$, $\tau_1$, and $Q$ are evenly spaced. The CSV gives the DSBP template period in microseconds (`Tp_us`) and the CR–RC time constant in nanoseconds (`tau1_ns`). The scale index (`raw_scale_index`) starts at 1 for each shape value.

For RLC, `original_frequency_grid_index` gives the position in the 3 dB passband grid. The `f0_MHz` field is blank for these rows. Converting the index to MHz requires the numerical frequency grid.

The final pulses use a one-sided 99% energy bandwidth of 10 MHz with a tolerance of 0.5%, followed by unit-energy normalization. The normalized CR–RC pulses at $r=0.25$ and $r=4$ differ by an overall sign. The CSV gives each of them its own candidate ID.

## RLC frequency scaling

The final oscillation frequency is $f_{\mathrm{osc}}=\alpha f_0$. The values in the last column use the rounded scale factors shown here.

| Waveform | $Q$ | $f_0$ (MHz) | $\alpha$ | $f_{\mathrm{osc}}$ (MHz) |
| --- | ---: | ---: | ---: | ---: |
| RLC representative | $25/9$ | 9.40 | 0.493 | 4.63 |
| RLC-end | 5 | 9.80 | 0.548 | 5.37 |

The waveform construction and grouping rules are given in manuscript Section 3.1 and Supplementary Note S2.1.
