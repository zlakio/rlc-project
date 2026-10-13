%% epistemic_uncertainty_sweep.m
%
% PURPOSE : Robustness under epistemic uncertainty (component tolerance) ON TOP
%           of aleatory noise (AWGN). Every trial draws fresh R, L, C values.
% INPUTS  : get_params.m, build_rlc_ladder.m
%           (optional overlay) results/rlc_analysis_extended.mat
% OUTPUTS : results/rlc_epistemic_results.mat
%           results/figures/fig_epi_pca_trajectory.png
%           results/figures/fig_epi_detection_comparison.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N = cfg.N; R_nom = cfg.R_nom; L_nom = cfg.L_nom; C_nom = cfg.C_nom;
RL = cfg.RL; Rs = cfg.Rs;
fault_stage = cfg.fault_stage; fault_range = cfg.fault_range;
snr_db = cfg.snr_db; trials_per_step = cfg.trials_per_step;
t = cfg.t; u = cfg.u;
R_tol = cfg.R_tol; L_tol = cfg.L_tol; C_tol = cfg.C_tol;

fprintf('Component tolerance bands: R +/-%.0f%%, L +/-%.0f%%, C +/-%.0f%%\n\n', ...
    R_tol*100, L_tol*100, C_tol*100);

numSteps  = length(fault_range);
totalRows = numSteps * trials_per_step;

eigTrackEpi   = zeros(totalRows, 2*N);
dataMatrix    = zeros(totalRows, length(t));
severityLabel = zeros(totalRows, 1);

row = 1;
tic;
for s = 1:numSteps
    for trIdx = 1:trials_per_step
        R_pert = R_nom .* (1 + R_tol * (2*rand(1,N) - 1));
        L_pert = L_nom .* (1 + L_tol * (2*rand(1,N) - 1));
        C_pert = C_nom .* (1 + C_tol * (2*rand(1,N) - 1));
        C_pert(fault_stage) = C_pert(fault_stage) * fault_range(s);

        [A, B, Cm, D] = build_rlc_ladder(R_pert, L_pert, C_pert, RL, Rs);
        sys = ss(A, B, Cm, D);

        % Sort eigenvalues (real part, then imag) so column k means the same
        % mode across randomly perturbed trials. Approximate: modes whose real
        % parts nearly coincide could still swap - treat per-mode stats as indicative.
        e = eig(A);
        [~, ix] = sortrows([real(e) imag(e)]);
        eigTrackEpi(row, :) = e(ix).';

        y_clean = lsim(sys, u, t);

        sigPower   = mean(y_clean.^2);
        noisePower = sigPower / (10^(snr_db/10));
        y_noisy    = y_clean + sqrt(noisePower) * randn(size(y_clean));

        dataMatrix(row, :) = y_noisy.';
        severityLabel(row) = fault_range(s);
        row = row + 1;
    end
    if mod(s,10) == 0
        fprintf('Completed severity step %d / %d (%.1f s elapsed)\n', s, numSteps, toc);
    end
end
fprintf('\nAll trials complete. Total simulation time: %.1f s\n\n', toc);

%% ---- Epistemic uncertainty band at the HEALTHY state ----
healthyRows = 1:trials_per_step;
[~, dominantIdx] = max(var(real(eigTrackEpi(healthyRows,:)), 0, 1));
healthyEigReal = real(eigTrackEpi(healthyRows, dominantIdx));
fprintf('Dominant eigenvalue (real part) spread from component tolerance ALONE (healthy):\n');
fprintf('  mean = %.3f, std = %.3f, range = [%.3f, %.3f]\n\n', ...
    mean(healthyEigReal), std(healthyEigReal), min(healthyEigReal), max(healthyEigReal));

%% ---- PCA on combined ensemble ----
[coeff, score, ~, ~, explained] = pca(dataMatrix);

f1 = figure;
scatter(score(:,1), score(:,2), 20, severityLabel, 'filled');
colormap(jet);
cb = colorbar;
cb.Label.String = 'Capacitance fraction remaining (1 = healthy)';
xlabel(sprintf('PC1 (%.1f%% variance)', explained(1)));
ylabel(sprintf('PC2 (%.1f%% variance)', explained(2)));
title('Latent Signature Trajectory - Aleatory + Epistemic Uncertainty');
grid on;
save_fig(f1, 'fig_epi_pca_trajectory', cfg);

%% ---- Detection threshold under combined uncertainty ----
baselineScores = score(healthyRows, 1:2);
mu    = mean(baselineScores, 1);
Sigma = cov(baselineScores);

mahalDist = zeros(totalRows, 1);
for i = 1:totalRows
    d = score(i, 1:2) - mu;
    mahalDist(i) = sqrt(d / Sigma * d');
end

confLevel = cfg.confLevel;
threshold = sqrt(-2 * log(1 - confLevel));   % = sqrt(chi2inv(confLevel, 2))

meanMahalEpi = zeros(numSteps, 1);
for s = 1:numSteps
    rows = (s-1)*trials_per_step + (1:trials_per_step);
    meanMahalEpi(s) = mean(mahalDist(rows));
end

detectedIdxEpi = find(meanMahalEpi > threshold, 1, 'first');
if isempty(detectedIdxEpi)
    fprintf('Under combined uncertainty, fault never crosses the %.0f%% threshold in range.\n', confLevel*100);
else
    fprintf('Under combined (aleatory + epistemic) uncertainty:\n');
    fprintf('  Fault first detected at capacitance = %.1f%% of nominal (step %d of %d)\n', ...
        fault_range(detectedIdxEpi)*100, detectedIdxEpi, numSteps);
end

%% ---- Comparison plot ----
f2 = figure; hold on;
plot(fault_range, meanMahalEpi, '-o', 'LineWidth', 1.5, 'DisplayName', 'Aleatory + Epistemic');
try
    prev = load(fullfile(cfg.resultsDir, 'rlc_analysis_extended.mat'), 'meanMahal', 'cfg');
    if isequal(prev.cfg.fault_range, fault_range)
        plot(fault_range, prev.meanMahal, '--s', 'LineWidth', 1.5, 'DisplayName', 'Aleatory only (baseline)');
    else
        warning('fault_range differs from the earlier run - skipping overlay. Re-run main_simulation + analyze_eigen_pca_link.');
    end
catch
    fprintf('(No prior aleatory-only results found - run main_simulation + analyze_eigen_pca_link first for the overlay.)\n');
end
yline(threshold, '--r', sprintf('%.0f%% detection threshold', confLevel*100), 'HandleVisibility', 'off');
set(gca, 'XDir', 'reverse');
xlabel('Capacitance fraction remaining (1 = healthy)');
ylabel('Mahalanobis distance (dimensionless)');
title('Detection Robustness: Aleatory-only vs. Aleatory+Epistemic');
legend('Location', 'best');
grid on;
hold off;
save_fig(f2, 'fig_epi_detection_comparison', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'rlc_epistemic_results.mat'), ...
    'dataMatrix', 'severityLabel', 'eigTrackEpi', 'fault_range', 't', ...
    'coeff', 'score', 'explained', 'mahalDist', 'meanMahalEpi', ...
    'threshold', 'detectedIdxEpi', 'R_tol', 'L_tol', 'C_tol', ...
    'healthyEigReal', 'dominantIdx', 'cfg');

fprintf('\nEpistemic analysis complete. Saved to %s\n', fullfile(cfg.resultsDir, 'rlc_epistemic_results.mat'));
