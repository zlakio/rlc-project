# Eigenvalue-Anchored PCA for Pre-Failure Fault Detection in N-Stage RLC Networks

**Full title:** *Eigenvalue-Anchored PCA for Pre-Failure Fault Detection Under Aleatory and Epistemic Uncertainty in N-Stage RLC Networks*
**Context:** Control Systems course project, VIT Vellore
**SDG alignment:** SDG 9 (primary), SDG 12 (secondary), SDG 7 (supporting)

A simulation study that detects slow component degradation in an RLC ladder network *before* it becomes a hard failure, and shows how much of that detection survives measurement noise and manufacturing tolerance.

---

## 1. What this project does

1. Model an N-stage RLC ladder as a linear state-space system.
2. Slowly degrade one capacitor (a "pre-failure" fault: 100% down to 60% of nominal capacitance).
3. Simulate noisy step responses at every severity.
4. Compress the responses with **PCA** to get a low-dimensional "latent signature trajectory".
5. **Anchor** that trajectory to physics: track how the system's eigenvalues migrate under the same fault and show the PCA movement follows it.
6. Turn the latent space into a **detector** (Mahalanobis distance / T², plus the Q statistic) and measure it honestly: false alarm rate, detection rate vs. severity, ROC/AUC, minimum detectable fault.
7. Repeat under two kinds of uncertainty:

| Uncertainty | Meaning | How it is modelled |
|---|---|---|
| **Aleatory** | Random measurement noise | Additive white Gaussian noise (15 dB SNR) |
| **Epistemic** | Not knowing exact component values | Monte Carlo draws of R, L, C within manufacturing tolerance (R ±5%, L ±10%, C ±10%) |

---

## 2. Pipeline

```
 component values ──► build_rlc_ladder ──► A, B, C, D (state space)
        │                                     │
        │ (fault: capacitor of middle         ├──► eig(A) ──► eigenvalue migration map
        │  stage scaled 1.0 → 0.6)            │
        │                                     └──► step response y(t)
        │                                              │
        └─ (optional: random tolerance)                ▼
                                              + AWGN (15 dB SNR)
                                                       │
                                                       ▼
                                  data matrix: one noisy response per row
                                                       │
                                                       ▼
                                           PCA ──► scores (PC1, PC2, ...)
                                                       │
                          eigenvalue drift  ◄─ correlate ─► PC1 score   ("anchoring")
                                                       │
                                                       ▼
                       healthy baseline ──► T² (Mahalanobis) and Q (residual) ──► alarm
```

## 3. The model

Each stage *k* has a series resistor `R_k` and inductor `L_k` and a shunt capacitor `C_k` to ground. A Thevenin source (internal resistance `Rs`) drives stage 1; a load resistor `RL` terminates stage N.

```
L_k · di_k/dt = v_(k-1) − v_k − R_k · i_k
C_k · dv_k/dt = i_k − i_(k+1)            (last stage also drains through RL)
```

State vector `x = [iL1; vC1; ...; iLN; vCN]` (2N states), input `u` = source voltage, output `y` = voltage across the last capacitor. The eigenvalues of `A` are the circuit's natural modes; a drifting component moves them. Default: N = 5, R = 10 Ω, L = 10 mH, C = 100 µF, RL = 50 Ω, Rs = 5 Ω, fault in stage 3.

## 4. Design decisions and what we learned along the way

This section records the thinking, including the things that were wrong at first.

**Why PCA, and why "anchor" it?** A step response is a 2001-sample vector, but the fault only moves it along a one-parameter curve. PCA finds that low-dimensional structure. PCA alone is purely statistical and has no idea about circuits, so we tie it to physics: if the PC1 score tracks the dominant eigenvalue's drift, the latent trajectory is interpretable as "the circuit's natural mode is moving".

**Bug found: eigenvalue ordering.** `eig()` does not guarantee a consistent order, and complex-conjugate pairs can swap. In the first run the "dominant" eigenvalue jumped from `-562 + 1188i` to `-560 - 1287i` (a jump of 2457 between consecutive steps), so the correlation of 0.43 was meaningless. Fix: track each mode by matching it to its nearest neighbour from the previous severity step. The largest jump dropped to 7.6 and the correlation became **r = 0.99**. Caveat: two smoothly trending quantities correlate easily, so this is supporting evidence. The analytical sensitivity check (Section 5, `eigen_sensitivity.m`) is the stronger test and agrees with the simulation.

**The first detector was optimistic.** The initial approach fit the healthy baseline on 20 healthy trials, scored those same trials, and flagged a severity when the *mean* distance crossed the threshold. That reported detection at 94.6% capacitance. Problems: the baseline was evaluated on its own training data, and a mean cannot give a false alarm rate, detection rate, or ROC. The metrics script (README item 2.1) fixes this: fit on separate healthy training trials, score fresh trials, and alarm per trial.

**Baseline size matters.** The first held-out evaluation fit the healthy baseline on only 200 trials and gave a false alarm rate of 2.9% and a minimum detectable fault of 89.2% capacitance. Refitting on 1000 trials dropped the false alarm rate to 0.3% and moved the minimum detectable fault to 87.1%. The earlier figure was flattered by a loose, noisy baseline. Detection numbers are only trustworthy once the baseline is large enough. (`detection_metrics.m` now uses 1000 baseline trials as well, so its outputs match `detection_metrics_q.m`.)

**Closed-form threshold.** For 2 PCs the Mahalanobis distance squared follows a chi-square distribution with 2 degrees of freedom, whose quantile is `-2·ln(1 − p)`. So the 99% threshold is `sqrt(-2·ln(0.01)) ≈ 3.03` with no toolbox needed.

**PCA does not beat a plain energy threshold.** The baseline comparison (Section 5) shows that with noise only, a one-line raw-energy threshold detects as well as the swept-PCA T² detector, and slightly better at the smallest faults. PCA fit on healthy data only is blind with 2 components (its components are just noise directions). The value of the PCA trajectory is that it can be anchored to eigenvalue migration, not that it detects smaller faults than simple alternatives.

**Tolerance, not noise, is the bottleneck.** Eigenvalue sensitivities show that a 10% capacitance loss moves the dominant mode by about the same amount as the spread caused by component tolerance alone (42 vs 44 rad/s). The SNR x tolerance map confirms it: once tolerance is present, a higher SNR stops helping beyond roughly 15-20 dB.

**Reproducibility.** Every shared parameter lives in `get_params.m`, every script seeds the RNG (`rng(42)`), and results and 300-dpi figures are saved to `results/`. `test_build_rlc_ladder.m` checks the model against closed-form eigenvalues for N = 1.

## 5. Results so far

All numbers: seed 42, 15 dB SNR, 99% confidence, PCA with 2 retained components. Detector baselines are fit on 1000 healthy training trials; every rate below is measured on fresh trials (1000 healthy, 100 per severity).

**Eigenvalue anchoring (aleatory-only sweep):** the dominant mode moves from `-508.4 - 1849.7i` to `-502.3 - 2114.1i` as capacitance falls from 100% to 60% (mostly a frequency shift; damping barely changes). Correlation with mean PC1 score: **r = 0.99**.

**Detection performance:**

| | Aleatory only | | | Aleatory + Epistemic | | |
|---|---|---|---|---|---|---|
| | T² | Q | T² OR Q | T² | Q | T² OR Q |
| False alarm rate | 0.3% | 1.4% | 1.7% | 1.4% | 1.2% | 2.3% |
| Minimum detectable fault (>=95% detection) | **87.1% cap. (12.9% loss)** | 60.0% (40% loss) | 87.1% (12.9% loss) | not reached | not reached | not reached |
| AUC at 95% capacitance | 0.805 | 0.670 | 0.722 | 0.494 | 0.608 | 0.594 |
| AUC at 90% capacitance | 0.995 | 0.765 | 0.954 | 0.529 | 0.614 | 0.605 |
| AUC at 80% capacitance | 1.000 | 0.905 | 1.000 | 0.665 | 0.781 | 0.776 |
| AUC at 60% capacitance | 1.000 | 0.992 | 1.000 | 0.917 | 0.944 | 0.942 |

**Detection rate (%) at selected severities:**

| Capacitance remaining | Aleat. T² | Aleat. Q | Aleat. OR | Comb. T² | Comb. Q | Comb. OR |
|---|---|---|---|---|---|---|
| 98% | 1 | 2 | 3 | 1 | 2 | 3 |
| 95% | 12 | 9 | 21 | 1 | 3 | 3 |
| 90% | 75 | 10 | 78 | 2 | 3 | 5 |
| 85% | 100 | 28 | 100 | 10 | 9 | 15 |
| 80% | 100 | 39 | 100 | 8 | 10 | 16 |
| 70% | 100 | 74 | 100 | 21 | 22 | 29 |
| 60% | 100 | 96 | 100 | 52 | 42 | 61 |

Reading:
- **Noise alone:** T² reliably catches faults of roughly 13% capacitance loss. Q on its own is much weaker; adding it to T² changes the minimum detectable fault very little but raises detection slightly just below the threshold (e.g. 12% -> 21% at 95% capacitance).
- **Noise plus component tolerance:** no method reaches 95% detection within a 40% capacitance loss. A +/-10% capacitor tolerance is as large as the fault being sought, so small faults look like normal manufacturing spread.
- **Q helps under tolerance:** it has higher AUC than T² at every severity tested (e.g. 0.608 vs 0.494 at 95%), consistent with part of the fault signature lying outside the first two PCs. T² OR Q roughly doubles detection at mid severities (80%: 16% vs 8%) but is still far from the 95% target.
- Detection rates use 100 trials per severity, so individual points carry roughly +/-5 percentage points of sampling noise.

### Baseline comparison (`baseline_comparison.m`)

Same protocol as above (seed 42, 15 dB, 99% confidence, 2 PCs): detectors fit on 1000 healthy training trials, thresholds for Q and the energy statistics from a separate 1000-trial healthy calibration set, every rate measured on fresh trials. Methods: **swept PCA T2** (the project's detector), **healthy-only PCA** (T2, Q, T2 OR Q; basis fit on healthy trials only), **raw energy** (`mean(y^2)`, two-sided), **deviation energy** (`mean((y - healthy mean response)^2)`).

| Method | Noise only: FAR | Min. detectable fault | AUC at 90% cap. | Noise + tolerance: FAR | Min. detectable fault | AUC at 90% cap. |
|---|---|---|---|---|---|---|
| Swept PCA: T2 | 0.2% | 87.1% cap. (12.9% loss) | 0.999 | 1.3% | not reached | 0.492 |
| Healthy PCA: T2 | 0.0% | not reached | 0.521 | 0.2% | not reached | 0.524 |
| Healthy PCA: Q | 0.3% | 74.2% cap. (25.8% loss) | 0.767 | 1.1% | not reached | 0.629 |
| Healthy PCA: T2 OR Q | 0.3% | 74.2% cap. (25.8% loss) | 0.767 | 1.2% | not reached | 0.628 |
| Raw energy | 0.8% | 89.2% cap. (10.8% loss) | 0.999 | 0.6% | not reached | 0.514 |
| Deviation energy | 0.3% | 74.2% cap. (25.8% loss) | 0.768 | 1.2% | not reached | 0.623 |

Detection rate (%) at selected severities:

| Capacitance remaining | Noise only: swept T2 | Raw energy | Healthy Q | Noise + tol.: swept T2 | Raw energy | Healthy Q | Deviation energy |
|---|---|---|---|---|---|---|---|
| 95% | 12 | 28 | 5 | 2 | 1 | 1 | 2 |
| 90% | 86 | 92 | 10 | 4 | 6 | 7 | 11 |
| 80% | 100 | 100 | 61 | 10 | 16 | 12 | 23 |
| 70% | 100 | 100 | 100 | 31 | 38 | 38 | 47 |
| 60% | 100 | 100 | 100 | 61 | 64 | 58 | 77 |

Reading:
- **Noise only:** a raw-energy threshold matches the swept-PCA T2 detector (AUC 0.999 for both, 89.2% vs 87.1% minimum detectable capacitance) and is better at the smallest faults (28% vs 12% detection at 95% capacitance). PCA adds no detection power here.
- **Healthy-only PCA with 2 components never alarms** (0% at every severity, AUC 0.52): with noise-only healthy data the retained components are noise directions that do not contain the fault. Its Q statistic and deviation energy do detect, but only larger faults (25.8% loss). They behave almost identically, consistent with the healthy components removing almost nothing. Both respond to the square of the fault size, which is weak for small faults, whereas raw energy and the swept PC1 respond linearly (an interpretation, not separately tested).
- **Noise plus tolerance:** no method reaches 95% detection. T2 (swept or healthy) and raw energy fall to chance at 90% capacitance. Residual-type statistics do best (Q and deviation energy, AUC about 0.62-0.63); deviation energy detects the most at 60% capacitance (77% vs 61% for swept T2). Raw energy, the best method without tolerance, becomes useless because tolerance changes total energy as much as the fault does.
- Rates use 100 trials per severity (about +/-5 percentage points) and one seed. They differ by a few points from the Q-statistic table above because each script draws its own fresh trials.

### Analytical eigenvalue sensitivity (`eigen_sensitivity.m`)

Capacitor `C_k` enters `A` only through row `2k`, scaling the whole row by `1/C_k`. So for `C_k -> s*C_k`, `A(s) = A0 + (1/s - 1)*E_k` exactly (`E_k` = `A0` with all other rows zeroed), and each simple eigenvalue moves at `d(lambda_i)/ds = w_i' (dA/ds) v_i / (w_i' v_i)` with left/right eigenvectors `w_i`, `v_i`. The script checks the identity numerically, then predicts the migration two ways and compares both with the simulated eigenvalue tracks from `main_simulation.m`.

| Capacitance | Simulated shift of dominant mode | Error, first-order (linearised once at healthy point) | Error, integrated (re-linearised along the path) |
|---|---|---|---|
| 95% | 18.3 rad/s | 1.1 rad/s | 0.0006 rad/s |
| 90% | 42.2 rad/s | 5.2 rad/s (12%) | 0.001 rad/s |
| 80% | 93.3 rad/s | 21.8 rad/s (23%) | 0.003 rad/s |
| 60% | 264.5 rad/s | 119 rad/s (45%) | 0.011 rad/s |

- The integrated sensitivity reproduces all ten simulated eigenvalue paths to within 0.011 rad/s over the whole sweep. The one-shot first-order prediction is only good for small faults.
- Correlation of dominant-mode drift with mean PC1: **r = 0.9899** for the simulated eigenvalues and **r = 0.9899** for the analytic (integrated) prediction, so PC1 is explained equally well without any eigenvalue simulation. The first-order drift gives r = 0.9988 because it is exactly linear in capacitance; that mainly shows PC1 is nearly linear in capacitance.
- **Which capacitor the dominant mode feels** (shift per 1% capacitance change, stages 1 to 5): 0.91, 2.92, **3.64**, 2.19, 0.33 rad/s. The simulated fault stage (3) is the most observable one; a fault in stage 5 would be about 10x harder to see.
- **Fault vs tolerance** (shift at 90% capacitance vs 1-sigma spread caused by component tolerance alone, linearised / Monte Carlo):

| Mode | Fault shift | Tolerance spread | Shift / spread | Share of spread from the faulted capacitor's own tolerance |
|---|---|---|---|---|
| -508 +/- 1850i (dominant) | 42.2 rad/s | 43.6 / 43.8 rad/s | 0.96 | 23% |
| -563 +/- 1188i | 24.3 rad/s | 32.2 / 34.1 rad/s | 0.71 | 18% |
| -135 (real) | 3.7 rad/s | 4.3 / 4.4 rad/s | 0.83 | 23% |
| remaining modes | 1.5 to 2.6 rad/s | 26 to 59 rad/s | <= 0.07 | <= 1% |

Even the best mode shifts by only about one tolerance standard deviation at 10% capacitance loss, which explains analytically why detection collapses under epistemic uncertainty. The tolerance of the faulted capacitor itself moves modes along the same direction as the fault, so that part cannot be separated from it.

### SNR x tolerance map (`snr_tolerance_heatmap.m`)

Every cell is a full experiment (swept-PCA basis, 500 healthy training trials, 500 fresh healthy and 100 fresh trials per severity on a 2% severity grid). Tolerance multiplier 1 is the R +/-5%, L +/-10%, C +/-10% bands above; 0 means noise only.

Minimum detectable capacitance loss, T2 detector (n.r. = 95% detection not reached within 40% loss):

| SNR | no tol. | 0.25x | 0.5x | 1x | 1.5x | 2x |
|---|---|---|---|---|---|---|
| 5 dB | n.r. | n.r. | n.r. | n.r. | n.r. | n.r. |
| 10 dB | 32% | 34% | n.r. | n.r. | n.r. | n.r. |
| 15 dB | 14% | 22% | 34% | n.r. | n.r. | n.r. |
| 20 dB | 8% | 18% | 32% | n.r. | n.r. | n.r. |
| 25 dB | 4% | 16% | 32% | n.r. | n.r. | n.r. |
| 30 dB | 2% | 18% | 32% | n.r. | n.r. | n.r. |

Detection rate (%) at 90% capacitance, T2 detector:

| SNR | no tol. | 0.25x | 0.5x | 1x | 1.5x | 2x |
|---|---|---|---|---|---|---|
| 5 dB | 3 | 1 | 1 | 1 | 1 | 4 |
| 10 dB | 14 | 11 | 6 | 2 | 4 | 2 |
| 15 dB | 70 | 20 | 5 | 0 | 2 | 1 |
| 20 dB | 100 | 26 | 6 | 3 | 1 | 3 |
| 25 dB | 100 | 52 | 6 | 6 | 4 | 4 |
| 30 dB | 100 | 44 | 10 | 4 | 2 | 4 |

Reading:
- **Noise only:** the minimum detectable fault falls steadily with SNR, from 32% loss at 10 dB to 2% at 30 dB; at 5 dB nothing is detectable.
- **Tolerance sets a floor that more SNR cannot lower:** about 16-18% loss at 0.25x and about 32% at 0.5x, flat from 20 to 30 dB, roughly 6-7 times the capacitor tolerance (two data points on a 2% grid, so only a rough guide). At the project's assumed bands (1x), nothing reaches 95% even at 30 dB.
- Once tolerance is present, raising SNR beyond about 15-20 dB buys almost nothing; reducing component uncertainty (measuring or calibrating parts) would help more.
- Adding Q (T2 OR Q) changes little: it lowers the 10 dB noise-only minimum from 32% to 28% loss and makes 40% loss detectable at 10 dB and 0.5x, with false alarms of at most 3.4%. All false alarm rates in the map stay below 5%.
- Resolution: 2% severity grid, 100 trials per severity, one seed. The (15 dB, no tolerance) cell reads 14% here versus 12.9% in the table above, most likely because of the coarser grid and the smaller sweep used to fit the PCA basis.

### Multi-point sensing (`multipoint_sensing.m`)

Everything above measures one signal, the voltage across the last capacitor. Here the same swept-PCA T2 detector is given more signals from the same simulated circuit. Each channel has its own noise at 15 dB and is scaled by its nominal healthy RMS before PCA. Protocol as in the baseline comparison (separate healthy training set, fresh test trials: 1000 healthy, 100 per severity). Two capacities are shown: 2 PCs (the project's choice) and 10 PCs.

Sensor sets: **S1** last capacitor V5 (what we had), **S2** V3 only (the faulted stage), **S3** V1+V3+V5, **S4** all five capacitor voltages, **S5** all five capacitor voltages + all five inductor currents.

Minimum detectable capacitance loss (>=95% detection, FAR <= 5%), 2 PCs / 10 PCs:

| Sensors | Noise only | Noise + tolerance |
|---|---|---|
| S1 V5 (baseline) | 11.5% / 14.2% | not reached |
| S2 V3 | 12.2% / 14.2% | not reached |
| S3 V1,V3,V5 | 8.8% / 10.8% | not reached |
| S4 all 5 caps | 6.8% / 8.8% | not reached |
| S5 caps + currents | 4.7% / 5.4% | **31.9% / 29.8%** |

Detection rate (%) under noise + tolerance, 2 PCs:

| Sensors | 70% cap (30% loss) | 60% cap (40% loss) | AUC at 90% cap |
|---|---|---|---|
| S1 | 21 | 54 | 0.491 |
| S2 | 48 | 84 | 0.609 |
| S3 | 38 | 73 | 0.573 |
| S4 | 36 | 72 | 0.561 |
| S5 | 92 | 100 | 0.650 |

Reading:
- **Noise only:** more sensors steadily lower the minimum detectable fault, from 11.5% to 4.7% loss, although each added channel brings its own noise.
- **With tolerance:** capacitor voltages alone do not help; S1-S4 still never reach 95%. Only S5, which adds the inductor currents, does, at about 30-32% loss. This is the first configuration in the project that clears the bar under the assumed tolerances.
- The improvement is for large faults only: at 90% capacitance S5 detects just 6% of faults and its AUC is 0.65, so a 10% fault is still not reliably separable from tolerance.
- Adding capacitor channels without currents (S2 vs S3 vs S4) slightly hurts at large faults. With 2 PCs the retained directions are dominated by tolerance spread, not the fault, which motivates nuisance-aware features.
- 10 PCs is worse than 2 PCs for the capacitor-only sets under tolerance, for the same reason.
- All false alarm rates are at most 1.8%. Resolution: one seed, 100 trials per severity.

## 6. Repository structure

```
.
├── build_rlc_ladder.m              % physics: state-space matrices for an N-stage ladder
├── get_params.m                    % single source of truth for all parameters
├── save_fig.m                      % exports figures as 300-dpi PNG
├── simulate_trials.m               % fresh noisy trials at one severity (optional tolerance)
├── main_simulation.m               % fault sweep + AWGN + PCA + eigenvalue migration
├── analyze_eigen_pca_link.m        % eigenvalue <-> PCA correlation + first threshold
├── epistemic_uncertainty_sweep.m   % tolerance Monte Carlo on top of AWGN
├── detection_metrics.m             % FAR, detection rate, ROC/AUC, minimum detectable fault
├── detection_metrics_q.m           % adds Q (SPE); reports T2, Q, T2 OR Q
├── baseline_comparison.m           % swept PCA vs healthy-only PCA vs raw/deviation energy thresholds
├── eigen_sensitivity.m             % analytical eigenvalue sensitivity vs simulated migration; tolerance spread
├── snr_tolerance_heatmap.m         % minimum detectable fault over an SNR x tolerance grid
├── multipoint_sensing.m            % detection with 1 to 10 measured signals (caps, inductor currents)
├── test_build_rlc_ladder.m         % model verification against closed form
├── run_all.m                       % runs everything in order
└── results/                        % created automatically (.mat, .csv, figures/)
```

## 7. Getting started

**Requirements:** MATLAB R2021a or newer, with **Control System Toolbox** (`ss`, `lsim`, `step`) and **Statistics and Machine Learning Toolbox** (`pca`, `chi2inv`, `perfcurve`, `prctile`).

```matlab
cd path/to/this/folder
run_all          % tests, then every script, in order (snr_tolerance_heatmap is the slow step)
```

Scripts read parameters from `get_params.m`; change values there, not inside scripts. Runs are deterministic (seed 42).

## 8. Limitations and scope

- **Simulation only.** The RLC ladder has not been built; planned hardware validation uses an RC ladder (no inductors).
- With the single output, detection under tolerance never reaches 95% within the tested range. Measuring all ten states (`multipoint_sensing.m`) reaches it only at about 30-32% loss; smaller faults remain undetectable under tolerance.
- The PCA basis is fit on swept data that includes faulty trials; metrics describe the method as defined. A healthy-only PCA basis was tested (`baseline_comparison.m`) and is blind with 2 components; a healthy-only basis with more components was not tried.
- A simple raw-energy threshold detects as well as the PCA T2 detector when only noise is present (see Section 5), so the PCA trajectory is justified by its physical interpretation, not by detection power.
- A single fault type (capacitance loss in one stage); the main results use a single measured output, and the multi-point study assumes up to ten independently noisy sensors.
- Noise level (15 dB) and tolerance bands are assumed values. The SNR x tolerance map sweeps both, but tolerance is modelled as independent uniform draws within each band; other distributions were not tried.
- The analytical eigenvalue check is done, but the Monte Carlo tolerance spread uses greedy nearest-neighbour mode matching, which is order-dependent for the -590 +/- 622i pair (Monte Carlo 40.0 rad/s vs linearised 26.2 rad/s). The linearised value is the reliable one there; those modes barely respond to the fault, so no conclusion depends on it.
- We do not address data heterogeneity, cascading-failure metrics, or the 2ⁿ state-explosion problem.

## 9. Roadmap

- [x] Reproducibility, saved results, model test
- [x] Detection metrics: FAR, detection rate, ROC/AUC, minimum detectable fault
- [x] Q statistic alongside T²: reported for T² alone, Q alone, and T² OR Q
- [x] Baseline comparison (plain PCA + T², raw-energy threshold)
- [x] Analytical eigenvalue sensitivity vs. simulated migration
- [x] SNR × tolerance heatmap
- [x] Multi-point sensing (several capacitor voltages, inductor currents)
- [ ] Self-baseline (per-unit reference) detection with sequential CUSUM/EWMA
- [ ] Kalman-filter joint state and capacitance estimation
- [ ] Nuisance-aware / supervised features
- [ ] RC-ladder hardware validation (Arduino)

