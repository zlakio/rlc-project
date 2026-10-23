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

**Detectors can lower the tolerance floor, but not remove it.** Three different upgrades (all ten signals, a supervised matched filter that knows the tolerance spread, an extended Kalman filter that estimates the faulted capacitance) all bring the smallest detectable loss under tolerance down from about 33% to about 22-23%, and none gets a 10% fault above about 41% detection with a fleet reference (Section 5). Comparing each unit with its own healthy record (self-baseline) removes the floor almost entirely in simulation: 12% loss on the last capacitor alone and 4% on ten signals, the same as with no tolerance, but it needs a per-unit reference and the fault shape. Under noise alone they all land within a percentage point or two of each other. Blindly removing the tolerance directions without using labels made things worse in a tried-and-dropped variant, because the faulted capacitor's own tolerance lies along the fault direction.

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

### Supervised matched filter (`nuisance_aware_features.m`)

`multipoint_sensing.m` showed a symptom: with 2 PCs the retained directions are dominated by tolerance spread. This script gives the detector the labels, on the same five sensor sets and protocol (one 60-component PCA basis fit on the swept data; separate healthy calibration set; fresh test trials). It compares:

- **D1** swept PCA T2, the project's detector (reference).
- **D2** Fisher matched filter (supervised): the fault signature is the slope of the data against capacity loss on the training sweep, whitened by the covariance of healthy trials (noise + tolerance); threshold from a separate healthy calibration set.

Minimum detectable capacitance loss (>=95% detection, FAR <= 5%), D1 / D2:

| Sensors | Noise only | Noise + tolerance |
|---|---|---|
| S1 V5 | 11.5% / 11.5% | n.r. / n.r. |
| S2 V3 | 12.9% / 12.2% | n.r. / n.r. |
| S3 V1,V3,V5 | 8.1% / 8.1% | n.r. / n.r. |
| S4 all caps | 6.8% / 6.1% | n.r. / n.r. |
| S5 caps + currents | 4.7% / 4.1% | 33.2% / **23.1%** |

Under noise + tolerance, detection rate (%) at 90% / 80% / 70% / 60% capacitance, and AUC at 90%:

| Sensors | D1 | D2 |
|---|---|---|
| S1 V5 | 7 / 10 / 27 / 48 (AUC 0.565) | 6 / 12 / 35 / 61 (0.690) |
| S2 V3 | 6 / 14 / 41 / 70 (0.545) | 10 / 24 / 56 / 86 (0.727) |
| S3 V1,V3,V5 | 3 / 10 / 37 / 63 (0.559) | 9 / 19 / 55 / 82 (0.731) |
| S4 all caps | 4 / 10 / 35 / 58 (0.571) | 8 / 14 / 45 / 82 (0.705) |
| S5 caps + currents | 11 / 37 / 90 / 99 (0.657) | 29 / 83 / 99 / 100 (0.863) |

Reading:
- **Noise only:** the two detectors are within about 1 point of each other; the detector is not the limit there.
- **With tolerance, the supervised filter helps:** AUC rises for every sensor set (for example 0.545 to 0.727 for V3), and on all ten signals the minimum detectable loss falls from 33% to 23% and detection of a 10% fault from 11% to 29%. Capacitor-only sets still do not reach 95%.
- **Labels matter, blind removal does not work.** An unsupervised variant that deleted the top 10 directions of the healthy covariance (the tolerance subspace) and then applied the 2-PC T2 was tried in an earlier version of this script and removed. It collapsed detection to roughly the false alarm rate on every sensor set (AUC about 0.5): the faulted capacitor's own tolerance lies along the fault direction, so projecting out the tolerance subspace removes the signal too. That result is described here but no longer reproduced by the script.
- **Interpretability is kept:** the cosine similarity between the regression fault signature and the model-based physical signature (derivative of the nominal response with respect to capacitance) is 0.94-0.96 for every sensor set, in both conditions. PC1 of the swept PCA also matches it (0.93-0.98 under tolerance), so the PCA direction is physically right even where its T2 distance fails to separate fault from tolerance. The matched filter itself is deliberately different from the signature under tolerance (cosine 0.13-0.38, against 0.91-0.94 with noise only), because it is shaped to avoid the tolerance directions: the supervised detector gains power but is less directly readable as "the mode moving".
- Caveats: the training sweep includes tolerance draws from the same bands as the test, so D2 knows the tolerance distribution and needs labelled faulty data that a real deployment would not have; `n.r.` = not reached; one seed, 100 trials per severity.

### Kalman-filter capacitance estimation (`kalman_capacitance.m`)

An extended Kalman filter tracks the ten circuit states plus one hidden state, theta = C(stage 3) / C(nominal), using the nominal model in which only theta may differ (the identity `A(theta) = A0 + (1/theta - 1)E` from `eigen_sensitivity.m`). Euler discretisation at 10 us; each measured channel is a scalar update with its own noise. Two versions: `qrel = 0` trusts the nominal model, `qrel = 1e-3` adds process noise of 0.1% of each state's RMS per step so the filter can absorb mismatch from unknown R, L, C. Decision rule: alarm if the final theta estimate falls below the 1st percentile of a separate healthy calibration set. 500 healthy and 100 trials per severity on a 2% grid.

Estimate quality at 90% true capacitance (mean estimate, standard deviation over trials, and the uncertainty the filter claims for itself), `qrel = 0`:

| Sensors | Noise only | Noise + tolerance |
|---|---|---|
| S1 V5 | 0.909, std 0.026, claims 0.025 | 0.887, std 0.138, claims 0.025 |
| S3 V1,V3,V5 | 0.908, std 0.018, claims 0.018 | 0.889, std 0.120, claims 0.018 |
| S5 caps + currents | 0.904, std 0.010, claims 0.009 | 0.891, std 0.086, claims 0.009 |

Minimum detectable capacitance loss, `qrel = 0` / `qrel = 1e-3`:

| Sensors | Noise only | Noise + tolerance |
|---|---|---|
| S1 V5 | 12% / 12% | n.r. / n.r. |
| S2 V3 | 12% / 14% | n.r. / n.r. |
| S3 V1,V3,V5 | 8% / 12% | n.r. / n.r. |
| S4 all caps | 6% / 12% | n.r. / 40% |
| S5 caps + currents | 4% / 6% | 34% / **22%** |

Under noise + tolerance, `qrel = 1e-3`: detection rate (%) at 90% / 80% / 70% / 60% capacitance, and AUC at 90%:

| Sensors | Detection | AUC |
|---|---|---|
| S1 V5 | 7 / 15 / 42 / 71 | 0.702 |
| S2 V3 | 9 / 32 / 69 / 92 | 0.727 |
| S3 V1,V3,V5 | 7 / 34 / 70 / 91 | 0.732 |
| S4 all caps | 9 / 38 / 75 / 96 | 0.746 |
| S5 caps + currents | 41 / 93 / 100 / 100 | **0.880** |

Reading:
- **Noise only: an efficient estimator.** `qrel = 0` is unbiased to within about 1 point and its standard deviation matches the uncertainty the filter claims (0.026 vs 0.025 on one signal, 0.010 vs 0.009 on ten). Capacitance is tracked to about 2.5% from the last node and about 1% from all ten signals. Detection is essentially the same as PCA's (minimum detectable loss 12% on one signal, 4% on ten), which suggests the PCA detector was already near the noise-limited optimum there; this is an inference from the estimator being efficient, not a proven bound.
- **Process noise costs accuracy when the model is exact:** with `qrel = 1e-3` and noise only, the estimate lags the truth (0.697 at a true 0.6 on one signal) and detection worsens.
- **With tolerance, the filter is overconfident:** the real scatter is 5 to 10 times the filter's own claim at `qrel = 0` (0.138 vs 0.025 on one signal, 0.086 vs 0.009 on ten) and still about 2 times at `qrel = 1e-3` on ten signals (0.056 vs 0.026). The filter's covariance describes noise only; epistemic uncertainty is invisible to it, which is exactly the problem named in the review paper.
- **Allowing for mismatch helps where it matters:** at `qrel = 1e-3` on all ten signals the minimum detectable loss is 22% (34% at `qrel = 0`), detection of a 10% fault is 41% and the AUC of 0.880 is the highest of any detector in this project. All-five-capacitor voltages reach 95% only at the 40% loss limit; the other capacitor-only sets do not.
- Caveats: `qrel = 1e-3` was picked after a small preliminary run on tolerance data, not cross-validated, and so encodes knowledge of the mismatch size; the filter assumes only the stage-3 capacitor differs from nominal; Euler discretisation; theta clipped to [0.3, 1.6]; one seed.

### Self-baseline detection (`selfbaseline_detection.m`)

Every detector above compares a measurement with a fleet reference, so a healthy unit with an unlucky capacitor looks faulty. Here each unit (one random draw of R, L, C within tolerance, fixed for its life) is compared with its own healthy response, recorded once when it was known to be good: its tolerance is in both signals and cancels. Each unit's reference is the average of m healthy recordings; the test is the average of K recordings at some capacity loss; the score is applied to the difference. Thresholds come from 1000 separate healthy calibration units (1% false alarm target); rates are measured on 500 fresh units per severity (2% grid, 0% to 40% loss).

- **SB-MF** self-baseline matched filter: score = <difference, g> / |g|, with g the model-based fault signature (derivative of the nominal response with respect to capacity loss). Needs the circuit equations, no fault examples.
- **SB-E** self-baseline difference energy: model-free.
- **Fleet-MF** control: the same matched filter against the nominal response (a fleet reference). Isolates what the self-baseline buys.
- Configurations: **A** m = 1, K = 1; **B** m = 10, K = 1; **C** m = 10, K = 5. Averaging K repeated measurements is what a CUSUM or EWMA would accumulate for a constant fault; a fault growing over time was not simulated.

Minimum detectable capacitance loss (>=95% detection, FAR <= 5%) under noise + tolerance, A / B / C:

| Sensors | SB-MF | SB-E (model-free) | Fleet-MF (control) |
|---|---|---|---|
| S1 V5 | 16% / 12% / **6%** | n.r. / 40% / 24% | n.r. / n.r. / n.r. |
| S2 V3 | 16% / 12% / 6% | n.r. / 38% / 26% | n.r. / n.r. / n.r. |
| S3 V1,V3,V5 | 10% / 8% / 6% | n.r. / 36% / 26% | n.r. / n.r. / n.r. |
| S4 all caps | 10% / 6% / 4% | n.r. / 34% / 24% | n.r. / n.r. / n.r. |
| S5 caps + currents | 6% / 4% / 4% | n.r. / 34% / 22% | 34% / 34% / 34% |

Under noise + tolerance, detection (%) of a 10% capacitance loss (AUC at 90% in brackets), A / B / C:

| Sensors | SB-MF | Fleet-MF (control) |
|---|---|---|
| S1 V5 | 66 / 91 / 100 (0.971 / 0.995 / 1.000) | 6 (0.684) |
| S3 V1,V3,V5 | 97 / 100 / 100 (0.997 / 1.000 / 1.000) | 7 (0.714) |
| S5 caps + currents | 100 / 100 / 100 (1.000) | 13 (0.784) |

Reading:
- **Tolerance is removed by calibration.** SB-MF under tolerance performs as it does with no tolerance at all (last capacitor, config A: 16% loss with tolerance, 16% without). The control, the same filter against a fleet reference, never reaches 95% on the last capacitor and needs 34% on all ten signals, so the gain comes from the per-unit reference and not from the filter.
- **A 10% fault becomes detectable** (66% to 100% on the last capacitor alone, depending on averaging), against 5-13% for the fleet reference and at most 41% for the best earlier detector.
- **Averaging helps as expected:** a single noisy reference costs about 4 points of loss against an averaged one (config A vs B, last capacitor: 16% vs 12%), and five repeated tests lower the floor to 6%.
- **The model-free version is much weaker.** SB-E needs averaged references and repeated tests: with config C about 22-26% loss under tolerance and 14-20% without it; with config B 34-40% under tolerance and 20-32% without. Most of the power comes from knowing the fault shape.
- Caveats: this is an idealised case. Every unit has a clean healthy reference; tolerance is fixed over the unit's life; only the stage-3 capacitor degrades (other parts ageing would look like a fault); the fault shape is the nominal one; the fault stage is known. Units that were already degraded when first measured cannot be handled this way. The first-order cancellation is not exact because the fault's effect depends slightly on the unit's other component values. One seed.

### What each fix buys under tolerance

| Detector | S1 V5 | S4 all caps | S5 caps + currents |
|---|---|---|---|
| Swept PCA T2 (project's, `multipoint_sensing.m`) | n.r. | n.r. | 31.9-33.2% loss |
| Fisher matched filter (supervised) | n.r. | n.r. | 23.1% |
| Kalman filter, `qrel = 0` | n.r. | n.r. | 34% |
| Kalman filter, `qrel = 1e-3` | n.r. | 40% | **22%** |
| Self-baseline matched filter, m = 10, K = 1 (own reference; different protocol, see above) | 12% | 6% | **4%** |

Detection of a 10% capacitance loss on the ten-signal set: 11% (PCA), 29% (Fisher), 41% (Kalman, `qrel = 1e-3`), 100% (self-baseline, which uses each unit's own healthy reference). Differences of a percentage point or two in these tables are within the 2% severity grid and trial-to-trial noise (the same PCA detector read 31.9% and 33.2% in two runs). The robust conclusions are: (1) every approach needs the inductor currents to reach 95% under tolerance; (2) using what is known about the tolerance (labels and whitening for the Fisher filter, a process-noise allowance for the Kalman filter) moves the floor from about 33% to about 22-23% loss; (3) a tried-and-dropped variant that removed the tolerance directions blindly made things worse; (4) with a fleet reference a 10% fault is still not reliably detectable under the assumed tolerances with any of these; (5) a per-unit healthy reference removes the tolerance problem almost entirely in simulation, but needs that reference and the fault shape.

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
├── selfbaseline_detection.m        % each unit against its own healthy reference (matched filter, energy, fleet control)
├── nuisance_aware_features.m       % supervised Fisher matched filter vs swept PCA, with physical-signature cosine check
├── kalman_capacitance.m            % extended Kalman filter estimating the faulted capacitance (slow, needs R2020b+)
├── test_build_rlc_ladder.m         % model verification against closed form
├── run_all.m                       % runs everything in order
└── results/                        % created automatically (.mat, .csv, figures/)
```

## 7. Getting started

**Requirements:** MATLAB R2021a or newer, with **Control System Toolbox** (`ss`, `lsim`, `step`) and **Statistics and Machine Learning Toolbox** (`pca`, `chi2inv`, `perfcurve`, `prctile`).

```matlab
cd path/to/this/folder
run_all          % tests, then every script, in order (snr_tolerance_heatmap and kalman_capacitance are the slow steps)
```

Scripts read parameters from `get_params.m`; change values there, not inside scripts. Runs are deterministic (seed 42).

## 8. Limitations and scope

- **Simulation only.** The RLC ladder has not been built; planned hardware validation uses an RC ladder (no inductors).
- With a fleet reference and the single output, detection under tolerance never reaches 95% within the tested range for any detector tried (PCA, supervised matched filter, Kalman filter). The self-baseline detector removes the floor but rests on a clean per-unit healthy reference, a fixed tolerance over the unit's life, a known fault shape and a known fault stage; ageing of other parts, or units that were already degraded at first measurement, would break it. A fault growing over time (true CUSUM/EWMA tracking) was not simulated. Measuring all ten states reaches it at about 33% loss with PCA and about 22-23% with the supervised filter or the Kalman filter; a 10% fault stays below about 41% detection.
- The PCA basis is fit on swept data that includes faulty trials; metrics describe the method as defined. A healthy-only PCA basis was tested (`baseline_comparison.m`) and is blind with 2 components; a healthy-only basis with more components was not tried.
- A simple raw-energy threshold detects as well as the PCA T2 detector when only noise is present (see Section 5), so the PCA trajectory is justified by its physical interpretation, not by detection power.
- The supervised matched filter is trained on sweeps that include tolerance draws from the same bands as the test, and the Kalman filter's process-noise setting was chosen from a small preliminary run on tolerance data; both therefore encode knowledge of the tolerance that a real deployment would have to obtain. The Kalman filter's own uncertainty ignores tolerance (it is overconfident by 2-10 times). A tried-and-dropped unsupervised removal of the tolerance subspace destroyed detection.
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
- [x] Self-baseline (per-unit reference) detection, repeat-averaged (true CUSUM/EWMA tracking of a drifting fault not simulated)
- [x] Kalman-filter joint state and capacitance estimation
- [x] Supervised features (Fisher matched filter)
- [ ] RC-ladder hardware validation (Arduino)

## 10. Relation to the literature

This project was built from the review by Wang and Geng (2025), "Modeling of complex electromechanical systems based on state-space, graph network and intelligent algorithms: a review" (*Journal of Reliability Science and Engineering* 1, 032002, doi:10.1088/3050-2454/ade71d). A review has no method to reproduce; it lists open problems in the reliability modelling of complex systems. This project works on some of them, on one fault mode in a small ladder, and leaves others untouched. One wording trap: the review's "state-space methods" mostly means discrete multi-state models (Markov, universal generating function); here "state-space" means the continuous-time electrical state equations of the circuit.

| Gap named in the review | Where it is addressed here | Status |
|---|---|---|
| Distinguishing aleatory from epistemic uncertainty (Sec. 2.2: "a persistent issue") | Noise vs component tolerance, separated everywhere; SNR x tolerance map; tolerance spread vs fault shift per mode (`eigen_sensitivity.m`) | Addressed for this circuit |
| Physical interpretability of data-driven methods (Sec. 2.2) | PC1 anchored to the dominant eigenvalue (r = 0.99); analytic eigenvalue sensitivity; cosine check of learned directions against the model-based fault signature | Addressed |
| Dimensionality reduction onto low-dimensional latent spaces that keep failure signatures (Sec. 4) | PCA latent signature trajectory | Addressed |
| Physics-data fusion and physics-informed digital twins with "correction and updating" (Sec. 2.2, 4, Fig. 5) | Extended Kalman filter that corrects a nominal model against the measurement and estimates the faulted capacitance with its own uncertainty (`kalman_capacitance.m`) | Partly: one parameter, nominal model |
| Field uncertainties such as sensor noise (Sec. 4) | Noise and tolerance both modelled; SNR swept 5-30 dB | Partly: assumed noise and tolerance models |
| Data heterogeneity across sensors (Sec. 2.2) | Multi-point sensing (`multipoint_sensing.m`): capacitor voltages and inductor currents | Partly: identical sampling, independent noise, up to ten channels |
| Multi-failure coupling, cascading and common-cause failure | Not addressed: one capacitor degrades at a time | Open |
| Scalability, state-space explosion (Sec. 2.2) | Not addressed: N = 5 stages; the continuous model grows linearly (2N states), unlike the 2^n of multi-state methods, but only N = 5 was studied | Open |
| Reliability metrics such as failure probability and mean time to failure | Not addressed: results are detection rates and minimum detectable fault | Open |
| Validation on physical or semi-physical platforms (Sec. 4) | Planned RC-ladder hardware check | Open |

The closest technical works cited in the review are Chen, Wang and Kang (2024), "Epistemic uncertainty propagation and reliability evaluation of feedback control system" (*IEEE Trans. Reliab.* 73, 521-535), and Xu, Bansal, Wang and Li (2025), "Physics-informed machine learning for system reliability analysis and design with partially observed information" (*Reliab. Eng. Syst. Saf.* 254, 110598). The Kalman and supervised-filter results (Section 5) follow the review's suggestion to combine a physics model with data-driven correction; the main finding, that component uncertainty rather than measurement noise sets the floor on detection, is a concrete instance of the epistemic-uncertainty problem it describes.

Reference: S. Wang and Y. Geng, "Modeling of complex electromechanical systems based on state-space, graph network and intelligent algorithms: a review," J. Reliab. Sci. Eng. 1, 032002 (2025).
