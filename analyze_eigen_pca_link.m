%% analyze_eigen_pca_link.m
%
% PURPOSE : (1) correlate dominant-eigenvalue drift with mean PC1 score,
%           (2) build a Mahalanobis detection threshold from healthy trials.
% INPUTS  : results/rlc_fault_simulation_results.mat (from main_simulation.m)
% OUTPUTS : results/rlc_analysis_extended.mat
%           results/figures/fig_link_eigen_vs_pc1.png
%           results/figures/fig_link_detection_curve.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');
load(fullfile(cfg.resultsDir, 'rlc_fault_simulation_results.mat'));

numSteps        = length(fault_range);
trials_per_step = size(dataMatrix, 1) / numSteps;
if mod(trials_per_step, 1) ~= 0
    error('dataMatrix row count is not evenly divisible by numSteps - results file does not match current parameters. Re-run main_simulation.m.');
end

%% ================= PART 1: Eigenvalue-PCA correlation =================
eigHealthy = eigTrack(1, :);

% Sanity check: consecutive steps should be close (eig() ordering assumption)
maxJump = max(max(abs(diff(eigTrack, 1, 1))));
fprintf('Largest eigenvalue jump between consecutive fault steps: %.4g\n', maxJump);
fprintf('(If unexpectedly large vs. eigTrack values, eig() may be reordering - flag before trusting Part 1.)\n\n');

[~, dominantIdx] = max(abs(eigTrack(numSteps, :) - eigHealthy));
eigDrift = abs(eigTrack(:, dominantIdx) - eigHealthy(dominantIdx));

fprintf('Dominant eigenvalue at healthy state:        %.4f %+.4fi\n', ...
    real(eigHealthy(dominantIdx)), imag(eigHealthy(dominantIdx)));
fprintf('Dominant eigenvalue at most degraded state:  %.4f %+.4fi\n\n', ...
    real(eigTrack(numSteps, dominantIdx)), imag(eigTrack(numSteps, dominantIdx)));

meanPC1 = zeros(numSteps, 1);
meanPC2 = zeros(numSteps, 1);
for s = 1:numSteps
    rows = (s-1)*trials_per_step + (1:trials_per_step);
    meanPC1(s) = mean(score(rows, 1));
    meanPC2(s) = mean(score(rows, 2));
end

Rmat = corrcoef(eigDrift, meanPC1);
r = Rmat(1, 2);
fprintf('Correlation between dominant-eigenvalue drift and mean PC1 score: r = %.4f\n', r);
if abs(r) < 0.5
    fprintf('WARNING: weak correlation. Try PC2, or a different eigenvalue.\n');
end

f1 = figure;
yyaxis left
plot(fault_range, eigDrift, '-o', 'LineWidth', 1.5);
ylabel('|\Delta\lambda_{dominant}| from healthy (rad/s)');
yyaxis right
plot(fault_range, meanPC1, '-s', 'LineWidth', 1.5);
ylabel('Mean PC1 score (a.u.)');
xlabel('Capacitance fraction remaining (1 = healthy)');
set(gca, 'XDir', 'reverse');
title(sprintf('Eigenvalue Drift vs. PCA Latent Trajectory  (r = %.3f)', r));
grid on;
legend('Eigenvalue drift', 'Mean PC1 score', 'Location', 'best');
save_fig(f1, 'fig_link_eigen_vs_pc1', cfg);

%% ================= PART 2: Detection threshold =================
healthyRows    = 1:trials_per_step;
baselineScores = score(healthyRows, 1:2);
mu    = mean(baselineScores, 1);
Sigma = cov(baselineScores);

mahalDist = zeros(size(score, 1), 1);
for i = 1:size(score, 1)
    d = score(i, 1:2) - mu;
    mahalDist(i) = sqrt(d / Sigma * d');
end

% chi2inv(p, 2) has the closed form -2*ln(1-p): no toolbox needed.
confLevel = cfg.confLevel;
threshold = sqrt(-2 * log(1 - confLevel));

meanMahal = zeros(numSteps, 1);
for s = 1:numSteps
    rows = (s-1)*trials_per_step + (1:trials_per_step);
    meanMahal(s) = mean(mahalDist(rows));
end

detectedIdx = find(meanMahal > threshold, 1, 'first');
if isempty(detectedIdx)
    fprintf('\nFault never crosses the %.0f%% threshold in the tested range.\n', confLevel*100);
else
    fprintf('\nFault first detected at capacitance = %.1f%% of nominal (step %d of %d)\n', ...
        fault_range(detectedIdx)*100, detectedIdx, numSteps);
    fprintf('Mahalanobis distance at detection = %.2f (threshold = %.2f)\n', ...
        meanMahal(detectedIdx), threshold);
end

f2 = figure;
plot(fault_range, meanMahal, '-o', 'LineWidth', 1.5, 'DisplayName', 'Mean Mahalanobis distance'); hold on;
yline(threshold, '--r', sprintf('%.0f%% detection threshold', confLevel*100));
if ~isempty(detectedIdx)
    plot(fault_range(detectedIdx), meanMahal(detectedIdx), 'ko', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'y', 'DisplayName', 'First detected');
end
set(gca, 'XDir', 'reverse');
xlabel('Capacitance fraction remaining (1 = healthy)');
ylabel('Mahalanobis distance (dimensionless)');
title('Fault Detection Curve');
legend('Location', 'best');
grid on;
save_fig(f2, 'fig_link_detection_curve', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'rlc_analysis_extended.mat'), ...
    'eigDrift', 'meanPC1', 'meanPC2', 'r', ...
    'mahalDist', 'meanMahal', 'threshold', 'detectedIdx', 'cfg');
fprintf('\nExtended analysis complete. Saved to %s\n', fullfile(cfg.resultsDir, 'rlc_analysis_extended.mat'));
