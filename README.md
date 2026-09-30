# Lyapunov-Budget Reference Governor

This repository contains the MATLAB code and certified numerical results accompanying the paper

> "Lyapunov-budget reference governor: certified state and input constraint satisfaction for adaptive
> strict-feedback systems", Bingcheng Liu (manuscript submitted to *Transactions of the Institute of
> Measurement and Control*).

It contains the governor and its adaptive backstepping pre-stabilizer, the interval branch-and-bound certificate
of the admissible level Γ_c(v), the certified threshold tables used in the paper, the baselines (PF+, BLF, ERG-N,
CE ablation), and scripts that reproduce every figure and table of the simulation study.

## 1. Repository contents

| Folder | Content |
|---|---|
| `code/core/` | `prestab.m` (tuning-function adaptive backstepping with frozen setpoint, both plants), `simulate.m` (vectorized RK4 closed loop: LBG, PF+, ERG-N/CE), `params_lbg.m`, `params_arm.m` |
| `code/baselines/` | `simulate_blf.m` (BLF adaptive backstepping), `compute_Gamma.m`, `make_table.m` (sampled, non-certified threshold; used for the ERG-N baseline of the academic example) |
| `code/certificate/` | `certify_level.m` (B&B interval certificate of one cell), `certify_level_rig.m` (outward-rounded version), `build_cert_table.m`, `build_cert_generic.m`, `build_all_tables.m`, `verify_tables_rig.m`, `test_certify4.m` |
| `code/academic_example/` | `run_envelope.m` (PF+ tuning), `run_baselines2.m` (Table 1), `run_stress_scen.m` (u_max = 1.8, R = 0.4), `check_k500.m`, `eval_academic_deployed.m` (quick check) |
| `code/manipulator/` | `run_arm_p1.m` (payload unknown), `run_arm.m` (payload and friction unknown), `diag_arm_stall.m`, `arm_gain_scan.m`, `eval_arm_deployed.m` (quick check) |
| `code/theory_checks/` | `sym_check.m`, `sym_check_arm.m` (Lemma 1, symbolic), `check_cex1d.m` (counterexample of Proposition 2) |
| `code/figures/` | `make_figs.m` (Figs. 2-6) |
| `code/utils/` | path, pool, working-folder and logging helpers |
| `scripts/` | entry points `reproduce_certificate.m`, `reproduce_academic.m`, `reproduce_manipulator.m`, `reproduce_figures.m`, `reproduce_all.m` |
| `results/certified_tables/` | certified threshold tables (`Gamma_cert*.mat`) and their outward-rounded versions (`*_rig.mat`) |
| `results/reference_results/` | frozen outputs of the tuning and evaluation runs reported in the paper |
| `results/logs/` | console output of the runs reported in the paper |
| `docs/REPRODUCIBILITY.md` | claim-by-claim map from the paper to scripts, data and logs |

## 2. Requirements

- MATLAB. All results were produced with MATLAB R2024b (24.2) on Windows 11 with 4 CPU cores. Older releases
  were not tested; the code uses `exportgraphics` (R2020a+) and `yline` (R2018b+).
- Symbolic Math Toolbox: only for `sym_check.m` and `sym_check_arm.m`.
- Optimization Toolbox: only for `test_certify4.m` (`fmincon`, tightness test of the certificate, 2 min).
- Parallel Computing Toolbox: optional. `parfor` loops (table construction and verification, BLF/PF+ grids of
  the scenario and manipulator scripts) run serially without it, with identical results but longer run time.
- No third-party packages. The interval arithmetic, including the outward rounding, is implemented in
  `certify_level_rig.m`.

## 3. Quick start

```matlab
cd('path/to/LBG_TIMC_public/scripts')
reproduce_certificate      % ~9 min: Lemma 1, Proposition 2, outward-rounded check of all 591 cells, tightness
reproduce_academic         % ~10 min: Table 1, transfer setpoints, u_max = 1.8 and R = 0.4
reproduce_manipulator      % ~4 min: manipulator, p = 1 and p = 2
reproduce_figures          % ~3 min: Figs. 2-6
% or all of them:
reproduce_all
```

Every script works in its own folder `output/<step>/` (created on first use, excluded from git), copies the
reference files it needs into `output/<step>/results/`, and writes logs to `output/<step>/logs/`. The shipped
files in `results/` are never modified. Each quick check prints `MATCH`/`DIFFERS` against the reference results.

**Quick vs. full mode.** The default `'quick'` mode loads the certified tables and the deployed designs used in
the paper, re-verifies the certificate, re-runs every reported simulation of the deployed designs and draws the
figures. It skips the offline tuning searches of the baselines and the construction of the tables from scratch.
`reproduce_all('full')` (or `'full'` for a single step) also repeats these: the 175-cell table
`Gamma_cert.mat` (~8 h single thread), the other tables (~1.5 min each with 4 workers), the PF+ grid (432
designs), the BLF grids (288 and 162 designs), the manipulator grids and the gain scan.

## 4. Reproducing paper results

| Paper item | Produced by | Reference data / log |
|---|---|---|
| Fig. 1 (block diagram) | drawn in the LaTeX source, no code | - |
| Fig. 2 `fig_gamma.pdf` | `make_figs.m` | `results/logs/make_figs.log` |
| Fig. 3 `fig_response.pdf` | `make_figs.m` | same |
| Fig. 4 `fig_ablation.pdf` | `make_figs.m` | same |
| Fig. 5 `fig_transfer.pdf` | `make_figs.m` | same |
| Fig. 6 `fig_arm.pdf` | `make_figs.m` | same |
| Table 1 | full: `run_envelope.m` (PF+ design) + `run_baselines2.m`; quick: `eval_academic_deployed.m` | `baselines2.mat`, `envelope_cert_settle.mat`, `logs/baselines2.log` |
| Tighter actuator / larger uncertainty | `run_stress_scen('S1',1.8,0.25)`, `run_stress_scen('S2',2.5,0.40)` | `stress_S1.mat`, `stress_S2.mat`, `logs/stress.log` |
| Manipulator (p = 1) | `run_arm_p1.m` | `arm_p1.mat`, `logs/arm_p1.log` |
| Theorem 2(c) illustration (p = 2) | `run_arm.m`, `diag_arm_stall.m`, `arm_gain_scan.m` | `arm.mat`, `logs/arm.log`, `logs/arm_gain_scan.log` |
| Certificate, 8 h / 175 cells | `build_cert_table.m` | `logs/build_cert_table.log` |
| Outward rounding, 591 cells | `verify_tables_rig.m` | `*_rig.mat`, `logs/verify_rig.log` |
| Tightness (violation at 1.02 x the level) | `test_certify4.m` | `logs/test_certify4.log` |
| Lemma 1 (exact identity) | `sym_check.m`, `sym_check_arm.m` | residual 0 |
| Proposition 2 (counterexample) | `check_cex1d.m` | printed stall at v = -2.950 |

`docs/REPRODUCIBILITY.md` lists every numerical statement of Section 7 with its source.

## 5. Certified threshold tables

All tables are piecewise-linear functions Γ_c(v) stored as a struct `Gt` with nodes `Gt.v`, node values
`Gt.G`, certified cell levels `Gt.b`, `Gt.Dmax = 1` and the Lipschitz constant `Gt.LG`
(`Gt.G(j) = min(b(j-1), b(j))`, Lemma 3 of the paper).

| File | Scenario | Cells | Used for |
|---|---|---|---|
| `Gamma_cert.mat` | academic, u_max = 2.5, R = 0.25, width 0.01 on [-0.90, 0.85] | 175 | Table 1, Figs. 2-5 |
| `Gamma_cert_S1.mat` | academic, u_max = 1.8, width 0.02 | 88 | Fig. 2, Fig. 5(b) |
| `Gamma_cert_S2.mat` | academic, R = 0.40, width 0.02 | 88 | larger-uncertainty scenario |
| `Gamma_cert_arm_p1.mat` | manipulator, payload unknown, width 0.02 on [-1.2, 1.2] | 120 | Fig. 2(b), Fig. 6, manipulator results |
| `Gamma_cert_arm.mat` | manipulator, payload and friction unknown | 120 | Theorem 2(c) illustration |
| `*_rig.mat` (5 files) | the five tables above re-verified by `certify_level_rig.m` | 591 | proof that no cell had to be reduced |
| `Gamma_cert_arm_p1_nom.mat`, `Gamma_cert_arm_nom.mat` | nominal model (R = 0, no adaptation) | 120 | ERG-N baseline of the manipulator |
| `Gamma_cert_arm_k*_g*.mat` | manipulator (p = 2) for other gains (k, g), width 0.04 | 60 | gain scan (`arm_gain_scan.m`) |

The tables were built with `certify_level.m` (interval bounds inflated by a relative 1e-12). All 591 cells of
the five tables used in the paper then passed `certify_level_rig.m`, in which every floating-point operation is
rounded outward by one ulp, the decimal data are rounded conservatively and `sin` is widened by 2 ulp (assuming
a library sine accurate to 1 ulp). The `*_rig.mat` tables are therefore identical to the originals; the
simulations load the originals.

## 6. Runtime

Times reported in the paper and in the logs (MATLAB R2024b, 4 cores): 175-cell table `Gamma_cert.mat` 28 616 s
single thread; S1/S2 tables 94 s and 90 s, manipulator table 79 s, each with 4 workers; outward-rounded
verification of all 591 cells 393 s; BLF grid of `run_baselines2.m` 11 615 s; PF+ grid of `run_envelope.m`
39 689 s (single thread, re-run for this release, `logs/run_envelope.log`). The quick steps took about
9, 9, 4 and 3 min (certificate, academic, manipulator, figures) on the same machine.

## 7. Reproducibility notes

- Integration: fixed-step RK4 with dt = 2e-3 s (`P.dt`), horizon 30 s (academic) and 20 s (manipulator), unless
  a script sets otherwise (`check_k500.m` uses 2e-4 and 5e-5 s). The theorems are stated in continuous time; the
  sampled-data effect of the integration step is the only approximation in the simulated safety results.
- Settling time: last exit from the band |x1 - r| <= 0.02, `Inf` if the run has not settled at the horizon.
  A run is counted as violating if a state bound is exceeded by more than 1e-6.
- Randomness: the simulations are deterministic (grids of initial states, parameters and setpoints). Random
  numbers appear only in `compute_Gamma.m` (sampled nominal threshold, `rng(1)`) and in `test_certify4.m`
  (multistart points, `rng(2)`), both seeded.
- `parfor` results do not depend on the number of workers. Results should be bit-identical on the same MATLAB
  release and platform; other releases or CPUs may change the last digits of floating-point results, and hence
  in rare cases a settling time by one recording interval (0.02 s).
- `results/reference_results/tight_lbgtune.mat` is a frozen input from the development history: a sampled
  (non-certified) threshold used only as the initial bisection bracket in `build_cert_table.m` and to define the
  grid of initial states in `run_envelope.m`. The resulting 12 initial states are exactly the grid points that
  satisfy Assumption 2 with the certified `Gamma_cert.mat`; `eval_academic_deployed.m` checks this.

## 8. Citation

If you use this code, please cite:

> Bingcheng Liu. Lyapunov-budget reference governor: certified state and input constraint satisfaction for
> adaptive strict-feedback systems. Manuscript submitted to *Transactions of the Institute of Measurement and
> Control*.

The journal reference will be added after publication.

## 9. License

The code and data in this repository are released under the MIT License, see [LICENSE](LICENSE).
