# Reproducibility map

Every numerical statement of the simulation study (Section 7 of the paper), with the script that produces it and
the shipped reference data. "quick" = checked by the default quick mode; "full" = only in full mode.
Paths: scripts in `code/`, data in `results/`.

## Common setup

| Item | Value | Where |
|---|---|---|
| Academic plant, θ0 = (0.5, 1), R = 0.25, c = (1.2, 0.8), u_max = 2.5, ε = 0.02 | | `params_lbg.m` + overrides at the top of each script |
| LBG gains k1 = k2 = 1, Γθ = 20 I, κ = 20, D_max = 1, η = 0.01 | | same; `Gt.Dmax` |
| w = -(z1 + k1 z2) | | `prestab.m` (output `dVdv`) |
| Manipulator J = 0.1, b = 10, θ0 = (0.5, 1), R = 0.3, c = (1.4, 1.2), u_max = 1.5 | | `params_arm.m` |
| RK4, dt = 2e-3 s; settling band 0.02 | | `simulate.m`, `simulate_blf.m` |
| Envelope: 5 x 5 grid of x(0) satisfying Assumption 2 (12 points), θ on 8 boundary points + centre, r ∈ {-0.6, -0.3, 0.3, 0.5, 0.7} → 540 runs | | `run_baselines2.m`; x(0) set in `reference_results/envelope_cert_settle.mat` |
| Transfer setpoints {-0.9, -0.8, -0.45, 0.4, 0.6, 0.75, 0.8, 0.85} | | `run_baselines2.m` (`Rtr`) |

## Certificate

| Statement | Source | Mode |
|---|---|---|
| 175 cells of width 0.01 on [-0.9, 0.85], about 8 h single thread | `build_cert_table.m`; `logs/build_cert_table.log` (28 616 s) | full |
| Width 0.02, 4 workers: 90-94 s per table | `run_stress_scen.m` / `build_all_tables.m`; `logs/stress.log` | full |
| Manipulator table on [-1.2, 1.2]: 79 s | `run_arm_p1.m`; `logs/arm_p1.log` | full |
| All 591 cells of the five tables pass the outward-rounded test unchanged | `verify_tables_rig.m`; `logs/verify_rig.log`; `*_rig.mat` | quick |
| Violations found by multistart optimisation at 1.02 x the certified level at v ∈ {-0.9, 0, 0.3} | `test_certify4.m`; `logs/test_certify4.log` (worst/limit 1.0016, 1.0117, 1.0103) | quick |
| Γ_c drops below F near v = -0.9 and v = 0.85 (Fig. 2) | `Gamma_cert.mat`; `logs/build_cert_table.log` | quick (figure) |
| det Ω0 = v² sin v | analytic (paper) | - |
| Lemma 1 for both examples | `sym_check.m`, `sym_check_arm.m` (residual 0) | quick |
| Proposition 2 (counterexample) | `check_cex1d.m` | quick |

## Academic example (Table 1, Figs. 3-5)

| Statement | Source | Mode |
|---|---|---|
| PF+ tuned on a grid of 432 designs | `run_envelope.m` (3 x 3 x 6 x 8 grid); design in `envelope_cert_settle.mat` (`bp`); re-run reproduces it exactly, `logs/run_envelope.log` | full |
| BLF tuned on 288 designs (162 in S1/S2 and for the manipulator) | `run_baselines2.m`; `run_stress_scen.m`; `run_arm_p1.m` | full |
| Table 1, all rows (median/worst ts, violating runs, transfer violations) | `run_baselines2.m`; `logs/baselines2.log`; `baselines2.mat` | quick: `eval_academic_deployed.m` |
| ERG-N: 204 of 540 runs never settle (38 %) | same (design k = 3, κ = 500) | quick |
| CE ablation: 7.4 % violating runs, up to 32.4 % at transfer setpoints | same | quick |
| Fig. 3: max abs(u) = 0.87, no saturation; s - V constant | `make_figs.m`; `logs/make_figs.log` | quick |
| Fig. 4: run r = 0.7, θ = (0.75, 1), x(0) = (0.5, -0.3) | `make_figs.m` (first violating CE run); `logs/make_figs.log` | quick |
| Fig. 5(a): PF+ up to 20.4 % at transfer setpoints, LBG 0 % | `make_figs.m`, `run_baselines2.m` | quick |
| u_max = 1.8: table in 94 s; BLF re-tuned on 288 runs; BLF 12.5-37.5 % at 7 of 9 transfer setpoints, LBG 0 % | `run_stress_scen('S1', 1.8, 0.25)`; `logs/stress.log`; `stress_S1.mat` | quick (evaluation), full (tuning) |
| R = 0.4: LBG no violation on envelope and transfer | `run_stress_scen('S2', 2.5, 0.40)`; `logs/stress.log`; `stress_S2.mat` | quick |
| κ = 500 violations at dt = 2e-3 vanish at dt = 2e-4 and 5e-5 (integration-step caveat) | `check_k500.m`; `logs/check_k500.log` | full |

## Manipulator (Fig. 6)

| Statement | Source | Mode |
|---|---|---|
| 13 admissible x(0), 5 payloads, r ∈ {±0.3, ±0.6, ±1} → 390 runs; transfer {±0.45, 0.8, -0.9, 1.1, ±1.2} | `run_arm_p1.m`; `logs/arm_p1.log` | quick: `eval_arm_deployed.m` |
| LBG: no violation on envelope and transfer; every envelope run settled; median ts 3.6-11.7 s up to abs(r) = 1.1 | same | quick |
| r = ±1.2: Γ_c < F, 39 of the 65 runs at each of the two setpoints r = -1.2 and r = +1.2 rad stop safely short of the target | same (`eval_arm_deployed.m` prints the unsettled count) | quick |
| ERG-N: 98 of 390 violating, 312 not settled; BLF 7.7 % at r = ±1.2 | same | quick |
| Fig. 6: torque below 1.1 N m, speed below 0.4 rad/s | `make_figs.m`; `logs/make_figs.log` | quick |
| p = 2: all 702 runs safe, 25 % (174) stop short with s∞ = Γ_c(v∞) | `run_arm.m`; `diag_arm_stall.m`; `logs/arm.log` | quick (counts), full (diagnosis) |
| Changing the adaptation gain does not remove the stalls | `arm_gain_scan.m`; `logs/arm_gain_scan.log`; `Gamma_cert_arm_k*_g*.mat` | full |

## Totals

"More than 4000 LBG runs, two plants, two actuator limits, two parameter sets, 26 setpoints": LBG runs in the
reference results are 540 + 8 x 108 (academic), 4 x 72 + 9 x 72 (u_max = 1.8), 4 x 99 + 9 x 99 (R = 0.4),
390 + 7 x 65 (manipulator, p = 1) and 702 + 7 x 117 (p = 2), i.e. 5993 runs.
