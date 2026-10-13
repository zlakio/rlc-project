# Study Guide: Understanding the RLC Fault-Detection Project from the Basics

*This is for you, not for GitHub. It builds up from first principles to everything done so far, explains what each script does in the background, and ends with likely questions and honest answers.*

---

## 0. The one-minute story (say this first)

> "We simulate a ladder of resistors, inductors and capacitors and slowly damage one capacitor. We record the circuit's response to a step input, add realistic noise, and use PCA to compress each response into two numbers. We show that the way those numbers move matches how the circuit's natural frequencies, its eigenvalues, move. Then we use a statistical distance from the healthy cluster as a fault detector, and measure honestly how small a fault it can catch. With noise alone it catches about a 13% capacitance loss. When we also include manufacturing tolerance, small faults get hidden."

If you can say that and defend each sentence, you understand the project. The rest of this guide supports each sentence.

---

## 1. The circuit

- **Resistor (R)** resists current and burns energy as heat.
- **Inductor (L)** resists *changes* in current; it stores energy in a magnetic field.
- **Capacitor (C)** resists *changes* in voltage; it stores energy in an electric field.

Our circuit is a **ladder**: stage 1 is a series R and L, then a capacitor to ground; stage 2 repeats, and so on for N = 5 stages. A source with internal resistance `Rs` drives stage 1, and a load resistor `RL` sits at the end.

A **step input** means the source jumps from 0 V to 1 V at time zero and stays there. The circuit's output (voltage across the last capacitor) then rises and wobbles before settling. That rise-and-wobble shape is the circuit's **fingerprint**. Change a component and the fingerprint changes slightly. That is the whole idea: *a degrading component changes the fingerprint before the circuit actually fails.*

Default values: R = 10 Ω, L = 10 mH, C = 100 µF in every stage, RL = 50 Ω, Rs = 5 Ω.

---

## 2. State-space and eigenvalues (the physics core)

### 2.1 State-space

Instead of one big differential equation, we describe the circuit with a **state vector**: the minimum set of numbers that fully describes the circuit at an instant. For an inductor that is its current; for a capacitor, its voltage. With 5 stages we have 10 states:

```
x = [ iL1, vC1, iL2, vC2, ..., iL5, vC5 ]
```

The dynamics are written as:

```
dx/dt = A·x + B·u        (how the states change)
y     = Cm·x + D·u       (what we measure; here y = vC5 and D = 0)
```

`A` is a 10×10 matrix; each entry says "how strongly does state *j* push the rate of change of state *i*". Examples from `build_rlc_ladder.m`:

- `A(iL_k, iL_k) = -R_k/L_k` : an inductor's current decays through its own series resistor.
- `A(iL_k, vC_k) = -1/L_k` : the capacitor voltage opposes the inductor's current.
- `A(vC_k, iL_k) = +1/C_k` : current flowing into a capacitor raises its voltage.

Note `1/C_k` appears in `A`. **That is why changing a capacitor changes `A`, and therefore everything below.**

### 2.2 Eigenvalues = natural modes

The eigenvalues of `A` are numbers λ such that the circuit can "ring" as `e^(λt)`. They are usually complex pairs: `λ = σ ± jω`.

- **Imaginary part ω** is the oscillation frequency (rad/s). Divide by 2π for Hz.
- **Real part σ** is the decay rate. It must be negative for the circuit to be stable (left half of the complex plane). The time constant is `1/|σ|`.

Worked example from our run: the dominant eigenvalue was about `-508 − 1850j`. That is an oscillation at about 1850 rad/s, roughly 294 Hz, that dies out with a time constant of about 2 ms.

At 60% capacitance it became `-502 − 2114j`: the frequency rose to about 336 Hz while the decay barely changed. Intuition: a smaller capacitor makes the "spring" stiffer, so the circuit rings faster.

**Eigenvalue migration map** = the plot of where each eigenvalue sits in the complex plane as the fault grows. Eigenvalues drift, and none should cross into the right half (instability).

### 2.3 Verification

`test_build_rlc_ladder.m` checks the code against a pencil-and-paper result. For N = 1 with no load and no source resistance, the circuit is a plain series RLC with `λ = -R/(2L) ± sqrt((R/(2L))² − 1/(LC))`. The test compares this to `eig(A)` and also checks two other N = 1 cases and that the N = 5 system is stable. All four passed. This is your "is the model correct?" answer.

---

## 3. The fault and the data

- **Fault:** the capacitor in the middle stage (stage 3) is multiplied by a factor from 1.00 down to 0.60 in 60 steps. 1.0 is healthy; 0.6 is a 40% loss of capacitance (electrolytic capacitors really do lose capacitance as they age, which is why this is a realistic "pre-failure" mode).
- **Simulation:** for each severity, `lsim` computes the step response over 20 ms at 10 µs spacing, which is **2001 samples**.
- **Data matrix:** every noisy response becomes one row. 60 severities × 20 noisy repeats = **1200 rows × 2001 columns**. This `dataMatrix` is the input to PCA.

---

## 4. Uncertainty: aleatory vs. epistemic

| | Aleatory | Epistemic |
|---|---|---|
| Meaning | Randomness you cannot remove | Not knowing the true values |
| Example | Sensor noise | Resistors are labelled 10 Ω but are really 9.6 or 10.4 |
| Model | Add random noise to each response | Randomly perturb R, L, C within tolerance per trial |
| In reality | Measure the same thing twice, get different readings | Two "identical" circuits behave differently |

### Noise (AWGN)

AWGN = additive white Gaussian noise. We set it by **signal-to-noise ratio**:

```
noisePower = signalPower / 10^(SNR_dB / 10)
```

At 15 dB, `10^1.5 ≈ 31.6`, so the noise power is about 1/32 of the signal power. Noise is added independently to each of the 2001 samples of each trial.

### Tolerance

Each trial draws every component uniformly within its band: R ±5%, L ±10%, C ±10% (the fault is then applied on top). This is a **Monte Carlo** approach: run many random trials and look at the spread of outcomes.

### Random seed

`rng(42)` fixes the random number generator, so noise and tolerance draws are identical on every run. That is why your MATLAB Online and laptop results should match.

---

## 5. PCA from scratch

**Problem:** each response is a point in 2001-dimensional space. We cannot see that.

**PCA** finds the directions along which the 1200 points vary the most.

1. Subtract the mean response from every row ("centering").
2. Find orthogonal directions (principal components, PCs) ordered by how much variance each explains.
3. Project each row onto the first few PCs. These coordinates are the **scores**.

In MATLAB: `[coeff, score, ~, ~, explained] = pca(dataMatrix)`.

- `coeff` (2001 × k): the PC directions. Each is itself a "response shape".
- `score` (1200 × k): where each trial sits in PC space. Column 1 is PC1, column 2 is PC2.
- `explained`: percent of variance per PC.

**Why it works here:** the fault changes the response along a smooth, one-parameter path. The biggest source of variation across the dataset is this fault, so PC1 mostly captures "how degraded is the circuit". Plotting PC1 against PC2 coloured by severity gives the **latent signature trajectory**: a smooth curve from healthy to degraded, with noise scattering points around it.

---

## 6. "Eigenvalue-anchored": tying PCA to physics

PCA by itself is statistics. It would give a trajectory for any data. The "anchor" claim is that the PC1 movement is *explained by* the eigenvalue movement.

**How we test it** (`analyze_eigen_pca_link.m`):

1. Pick the eigenvalue that moves the most between healthy and most degraded.
2. Compute its drift `|λ(s) − λ(healthy)|` at each severity.
3. Compute the mean PC1 score at each severity.
4. Correlate the two. Result: **r = 0.99**.

**How to talk about it honestly:** two quantities that both change smoothly and in one direction will correlate strongly, almost regardless of mechanism. So r = 0.99 is supporting evidence. The stronger test, still to do, is README item 2.2: use perturbation theory to *predict* the eigenvalue movement analytically and overlay it on the simulated movement.

### The bug we found and fixed

First run: the correlation was only 0.43, and the "largest eigenvalue jump" between consecutive steps was 2457, enormous. The dominant eigenvalue appeared to go from `-562 + 1188j` to `-560 − 1287j`. The real part barely changed and the imaginary part flipped sign. Cause: `eig()` returns eigenvalues in no guaranteed order, and a complex-conjugate pair `a ± bj` can swap places between calls, so we were comparing two different modes.

Fix (in `main_simulation.m`): at each step, match every eigenvalue to its **nearest neighbour** from the previous step. Afterwards the largest jump was 7.6 and r became 0.99. This is a good story for Q&A: you found a silent bug by sanity-checking a number.

---

## 7. From plot to detector: Mahalanobis distance

### 7.1 Why not ordinary distance?

Euclidean distance treats all directions equally. But the healthy cluster is not round: it may be wider in one direction than another. **Mahalanobis distance** measures how many standard deviations away a point is, taking the cluster's shape into account:

```
d = sqrt( (x − μ) · Σ⁻¹ · (x − μ)ᵀ )
```

`μ` = mean of the healthy cluster, `Σ` = its covariance matrix. In code: `sqrt(d / Sigma * d')`.

This distance is also called the **T² statistic** (Hotelling's T², which is d²) in process monitoring.

### 7.2 The threshold

If the healthy scores are roughly Gaussian in 2 dimensions, then d² follows a chi-square distribution with 2 degrees of freedom. We choose 99% confidence: only 1% of healthy points should exceed it.

For 2 degrees of freedom the chi-square quantile has a simple closed form, `-2·ln(1 − p)`:

```
threshold = sqrt( -2 · ln(0.01) ) = sqrt(9.21) ≈ 3.03
```

That is why the code does not need `chi2inv`. A point with d > 3.03 raises an alarm.

---

## 8. Why the first detection result was optimistic, and how we measured properly

### The original approach (`analyze_eigen_pca_link.m`)

Fit the healthy baseline on the first 20 trials, compute the mean distance per severity, and report the first severity where the mean crosses 3.03. Result: **94.6% capacitance**.

Three problems:

1. **Training and testing on the same data.** The healthy trials define the baseline and are then scored against it, so they look unrealistically close to it.
2. **Averages hide individual trials.** If the *mean* distance crosses the line, many single trials still do not. You cannot compute a false alarm rate or a detection rate from means.
3. **PCA itself was fit on all trials including faulty ones**, so PC1 was handed the fault direction. (Still true; listed as a limitation.)

### The fix (`detection_metrics.m`)

- Fit the healthy baseline (`μ`, `Σ`) on separate healthy **training** trials.
- Generate **fresh** test trials the detector has never seen.
- Raise an alarm **per trial**.

### The metrics, in plain language

| Metric | Question it answers |
|---|---|
| **False alarm rate (FAR)** | Of healthy trials, what fraction wrongly raise an alarm? Ideal near 1% at 99% confidence. |
| **Detection rate** | At a given fault severity, what fraction of trials correctly raise an alarm? |
| **ROC curve** | Plot detection rate vs. false alarm rate as the threshold slides. Shows the trade-off. |
| **AUC** | Area under the ROC. 0.5 = coin flip, 1.0 = perfect separation. |
| **Minimum detectable fault** | The smallest fault with ≥95% detection (at ≤5% false alarms). |

---

## 9. The results and how to read them

Setup for these numbers: baseline fit on 1000 healthy training trials; tests on 1000 fresh healthy trials and 100 fresh trials per severity.

**Aleatory only (noise), detection rate in %:**

| Capacitance | T² | Q | T² OR Q |
|---|---|---|---|
| 98% | 1 | 2 | 3 |
| 95% | 12 | 9 | 21 |
| 90% | 75 | 10 | 78 |
| 85% | 100 | 28 | 100 |
| 80% | 100 | 39 | 100 |
| 70% | 100 | 74 | 100 |
| 60% | 100 | 96 | 100 |

- **Headline: minimum detectable fault = 87.1% capacitance (12.9% loss)** with T². "Minimum detectable" means the first severity from which detection stays at or above 95%.
- False alarm rate: T² 0.3%, Q 1.4%, T² OR Q 1.7%. All well under the 5% limit.
- T² is the strong statistic here. Q alone only reaches 95% at 40% loss.

**Aleatory + epistemic, detection rate in %:**

| Capacitance | T² | Q | T² OR Q |
|---|---|---|---|
| 98% | 1 | 2 | 3 |
| 95% | 1 | 3 | 3 |
| 90% | 2 | 3 | 5 |
| 85% | 10 | 9 | 15 |
| 80% | 8 | 10 | 16 |
| 70% | 21 | 22 | 29 |
| 60% | 52 | 42 | 61 |

- No method reaches 95% detection, even at 40% loss. Best case (T² OR Q) catches 61% of trials at 60% capacitance.
- False alarm rates: T² 1.4%, Q 1.2%, T² OR Q 2.3%. The detector is quiet on healthy circuits; it just misses most faults.
- **Q is better than T² under tolerance** by AUC at every severity (e.g. 0.608 vs 0.494 at 95% capacitance; 0.944 vs 0.917 at 60%). T² at 95% capacitance has AUC 0.494, which is chance level. Interpretation: tolerance moves healthy responses around in ways the first two PCs mostly absorb, while part of the fault signature sits outside them, which is what Q measures. This is a plausible explanation, not a proven one.

**Why tolerance hides faults (the key physical insight):** with +/-10% capacitor tolerance, a *healthy* circuit's capacitor can already be 10% off nominal, so a 10% capacitance fault looks like ordinary manufacturing variation. The detector cannot tell "this unit was built with a slightly small capacitor" from "this unit's capacitor has degraded". This is what epistemic uncertainty means in practice, and it is a result, not a bug.

**Lesson about baseline size (good for Q&A):** our first held-out evaluation used a 200-trial healthy baseline and gave false alarm rate 2.9% and minimum detectable fault 89.2%. With a 1000-trial baseline, false alarms fell to 0.3% and the minimum detectable fault became 87.1%. The small baseline had made the threshold effectively too loose, flattering detection. Quote the 1000-trial numbers.

**Reading the figure:** each point is 100 trials, so values wobble by roughly +/-5 percentage points. Look at the overall shape, not single points. Left panel: T² jumps from near 0 to 100% between about 95% and 85% capacitance, while Q climbs slowly. Right panel: all three curves stay low and rise only past about 80% capacitance.

---

## 10. T² vs. Q (README item 2.4)

PCA keeps only a few directions (we keep 2 out of up to 1200). Every response then splits into two parts:

- The part **inside** the kept PCs, monitored by **T²** (Mahalanobis).
- The part **outside** them (the residual), monitored by **Q**.

```
x_hat = reconstruction of x from the kept PCs
Q = || x - x_hat ||^2           (squared prediction error, SPE)
```

A fault that moves the response *along* the kept PCs is caught by T². A fault that pushes the response into a direction PCA discarded is caught only by Q. Standard practice is to use both and alarm if **either** fires: **T² OR Q**.

Thresholds in `detection_metrics_q.m`:

- T²: the chi-square threshold (3.03 for 2 PCs at 99%).
- Q: the empirical 99th percentile of Q on healthy training trials. Each statistic is divided by its own threshold, so "alarm" means value > 1, and "T² OR Q" means `max(T²/thrT, Q/thrQ) > 1`.

**What we found:** under noise alone T² does nearly all the work and Q adds little. Under tolerance Q is the more informative of the two by AUC, and combining them roughly doubles detection at mid severities, but neither gets near 95%. Note T² OR Q is not automatically better by AUC: in the noise-only case its AUC at 95% capacitance (0.722) is below T² alone (0.805), because merging in the noisier Q score adds some false alarms.

---

## 11. What every file does in the background

| File | Role | What it does step by step |
|---|---|---|
| `get_params.m` | Settings | Returns one struct `cfg` with N, R/L/C, fault range, SNR, tolerance, seed, thresholds, folders. Creates `results/`. |
| `build_rlc_ladder.m` | Physics | Fills the 10×10 `A`, input vector `B`, output row `Cm` from R, L, C, RL, Rs. |
| `test_build_rlc_ladder.m` | Verification | Compares `eig(A)` with closed-form roots for N = 1; checks N = 5 stability. |
| `main_simulation.m` | Experiment | Loop over 60 severities: build `A`, track eigenvalues (nearest-neighbour matching), simulate clean response, add noise 20 times, store rows. Then PCA, three figures, save `.mat`. |
| `analyze_eigen_pca_link.m` | Core claim | Loads saved results, finds the dominant eigenvalue, correlates its drift with mean PC1, builds the first Mahalanobis threshold. |
| `epistemic_uncertainty_sweep.m` | Robustness | Same sweep, but every trial draws random R, L, C first. Sorts eigenvalues so modes align across trials. Then PCA and detection. |
| `simulate_trials.m` | Helper | Generates fresh noisy trials at one severity, with or without tolerance. Used by the metrics scripts. |
| `detection_metrics.m` | Honest evaluation | Trains the baseline on separate healthy trials; scores 1000 fresh healthy + 100 per severity; reports FAR, detection rates, ROC/AUC, min detectable fault; both uncertainty cases. |
| `detection_metrics_q.m` | T², Q, OR | Same, with 1000 training trials, adding the Q statistic. |
| `save_fig.m` | Helper | Exports a figure as a 300-dpi PNG. |
| `run_all.m` | Runner | Calls everything in order. |

Order matters: `main_simulation` and `epistemic_uncertainty_sweep` must run before the metrics scripts, because those load their saved PCA basis from `results/`.

**Shapes to remember:** `dataMatrix` 1200×2001; `score` 1200×(many PCs); `eigTrack` 60×10; `coeff` 2001×(many PCs).

---

## 12. Likely questions and honest answers

**Why PCA instead of a neural network?** PCA is linear, fast, interpretable, and needs no fault labels. Its components can be tied back to the circuit's physics, which a neural network would not give us for free.

**What does "eigenvalue-anchored" actually mean?** The PCA trajectory is linked to the eigenvalue movement of the system matrix. Currently supported by a correlation of 0.99; the analytical sensitivity check is the stronger proof and is still to do.

**Isn't r = 0.99 trivial?** Partly, yes: both quantities trend monotonically. That is why we plan the perturbation-theory comparison. Say this before they do.

**Why is the 94.6% detection not your final number?** It used the mean distance and scored the training data. The honest figure is 87.1% capacitance (12.9% loss), from per-trial alarms on held-out trials with a 1000-trial baseline.

**What is the difference between aleatory and epistemic uncertainty?** Aleatory: irreducible randomness (noise). Epistemic: lack of knowledge (we don't know exact component values), which could in principle be reduced by measuring them.

**Why can't it detect small faults once tolerance is included?** Because tolerance variation of ±10% in C is as large as the fault. The healthy population and mildly faulty population overlap.

**Why 15 dB SNR and these tolerances?** They are assumed, typical values, not measured. A sensitivity sweep over SNR and tolerance (README 2.5) is planned.

**Why 99% confidence?** A design choice trading false alarms against sensitivity. The ROC curves show the whole trade-off.

**Is this validated on hardware?** No. It is simulation only. Hardware would use an RC ladder because we do not have inductors.

**Does this work for other faults or components?** Not tested. We only degrade one capacitor, in the middle stage.

**What is the biggest weakness?** PCA is fit on data that includes faulty trials, so the latent space is partly "told" about the fault. A healthy-only PCA baseline is planned for the comparison step.

**Does the result change between runs?** No. The random seed is fixed at 42.

---

## 13. Glossary

- **State-space model:** `dx/dt = Ax + Bu`, `y = Cx + Du`.
- **Eigenvalue:** a number describing one natural mode of the system.
- **Stable:** all eigenvalues have negative real part.
- **Step response:** output when the input jumps from 0 to a constant.
- **AWGN:** additive white Gaussian noise.
- **SNR:** signal power divided by noise power, usually in dB.
- **Monte Carlo:** repeat random trials to see the spread of outcomes.
- **PCA:** projection onto directions of maximum variance.
- **Score:** a trial's coordinates in PC space.
- **Mahalanobis distance / T²:** distance from the healthy cluster in units of its own spread.
- **Q / SPE:** squared size of what PCA discarded.
- **FAR:** false alarm rate.
- **ROC / AUC:** trade-off curve and area under it.
- **Aleatory / epistemic:** randomness vs. lack of knowledge.

---

## 14. Where you are and what is left

**Done:** reproducibility, saved results and figures, model verification, eigenvalue tracking fix, honest detection metrics (FAR, detection rate, ROC/AUC, minimum detectable fault), and the Q statistic reported for T², Q, and T² OR Q.

**To do next:** baseline comparison (README 2.6); analytical eigenvalue sensitivity (2.2); SNR × tolerance heatmap (2.5); then updating the GitHub README with final results.

**Optional, if time permits:** fault localisation (2.3), scalability with N (3.1), harder fault scenarios (3.2), hardware validation (Tier 4).
