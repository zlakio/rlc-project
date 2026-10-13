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

**Bug found: eigenvalue ordering.** `eig()` does not guarantee a consistent order, and complex-conjugate pairs can swap. In the first run the "dominant" eigenvalue jumped from `-562 + 1188i` to `-560 - 1287i` (a jump of 2457 between consecutive steps), so the correlation of 0.43 was meaningless. Fix: track each mode by matching it to its nearest neighbour from the previous severity step. The largest jump dropped to 7.6 and the correlation became **r = 0.99**. Caveat: two smoothly trending quantities correlate easily, so this is supporting evidence; analytical sensitivity (below) is the stronger test.

**The first detector was optimistic.** The initial approach fit the healthy baseline on 20 healthy trials, scored those same trials, and flagged a severity when the *mean* distance crossed the threshold. That reported detection at 94.6% capacitance. Problems: the baseline was evaluated on its own training data, and a mean cannot give a false alarm rate, detection rate, or ROC. The metrics script (README item 2.1) fixes this: fit on separate healthy training trials, score fresh trials, and alarm per trial.

**Baseline size matters.** The first held-out evaluation fit the healthy baseline on only 200 trials and gave a false alarm rate of 2.9% and a minimum detectable fault of 89.2% capacitance. Refitting on 1000 trials dropped the false alarm rate to 0.3% and moved the minimum detectable fault to 87.1%. The earlier figure was flattered by a loose, noisy baseline. Detection numbers are only trustworthy once the baseline is large enough.

**Closed-form threshold.** For 2 PCs the Mahalanobis distance squared follows a chi-square distribution with 2 degrees of freedom, whose quantile is `-2·ln(1 − p)`. So the 99% threshold is `sqrt(-2·ln(0.01)) ≈ 3.03` with no toolbox needed.

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
├── test_build_rlc_ladder.m         % model verification against closed form
├── run_all.m                       % runs everything in order
└── results/                        % created automatically (.mat, .csv, figures/)
```

## 7. Getting started

**Requirements:** MATLAB R2021a or newer, with **Control System Toolbox** (`ss`, `lsim`, `step`) and **Statistics and Machine Learning Toolbox** (`pca`, `chi2inv`, `perfcurve`, `prctile`).

```matlab
cd path/to/this/folder
run_all          % tests, then every script, in order
```

Scripts read parameters from `get_params.m`; change values there, not inside scripts. Runs are deterministic (seed 42).

## 8. Limitations and scope

- **Simulation only.** The RLC ladder has not been built; planned hardware validation uses an RC ladder (no inductors).
- Detection under tolerance never reaches 95% within the tested range; more informative measurements (several nodes, not one output) have not been tried.
- The PCA basis is fit on swept data that includes faulty trials; metrics describe the method as defined. A healthy-only PCA baseline is planned.
- A single fault type (capacitance loss in one stage) and a single measured output.
- Noise level (15 dB) and tolerance bands are assumed values; sensitivity sweeps are planned.
- The eigenvalue correlation is supporting evidence; the analytical perturbation check is still to do.
- We do not address data heterogeneity, cascading-failure metrics, or the 2ⁿ state-explosion problem.

## 9. Roadmap

- [x] Reproducibility, saved results, model test
- [x] Detection metrics: FAR, detection rate, ROC/AUC, minimum detectable fault
- [x] Q statistic alongside T²: reported for T² alone, Q alone, and T² OR Q
- [ ] Baseline comparison (plain PCA + T², raw-energy threshold)
- [ ] Analytical eigenvalue sensitivity vs. simulated migration
- [ ] SNR × tolerance heatmap
- [ ] RC-ladder hardware validation (Arduino)
