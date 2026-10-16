%% detection_metrics.m  (README 2.1)
%
% PURPOSE : Turn "the threshold is crossed" into numbers:
%             - false alarm rate (FAR) on FRESH healthy trials
%             - detection rate vs. fault severity (per-trial alarms)
%             - ROC curves + AUC at selected severities
%             - minimum detectable fault (>= 95% detection at <= 5% FAR)
%           Run for aleatory-only AND aleatory+epistemic uncertainty.
%
% DETECTOR: unchanged from the project so far - PCA basis (first 2 PCs) taken
%           from the saved sweep results, Mahalanobis distance, 99% chi-square
%           threshold. What is NEW: the healthy baseline (mean, covariance) is
%           fit on separate healthy TRAINING trials, and every number below is
%           measured on FRESH test trials the detector has never seen.
%
% KNOWN LIMITATION: the PCA basis itself was fit on the swept data (healthy +
%           faulty). Metrics therefore describe the method as currently defined.
%           A healthy-only PCA is tested later (README 2.6 baseline).
%
% INPUTS  : results/rlc_fault_simulation_results.mat
%           results/rlc_epistemic_results.mat
% OUTPUTS : results/detection_metrics.mat, results/detection_table.csv
%           results/figures/fig_metrics_detection_rate.png
%           results/figures/fig_metrics_roc.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

nTrain       = 1000;   % healthy trials used to fit the baseline (README: 200 was too few)
nTestHealthy = 1000;   % fresh healthy trials -> false alarm rate
nTestFault   = 100;    % fresh trials per severity -> detection rate
rocSev       = [0.95 0.90 0.80 0.60];   % severities shown as ROC curves
detTarget    = 0.95;   % required detection rate
farMax       = 0.05;   % allowed false alarm rate
tableSev     = [0.98 0.95 0.90 0.85 0.80 0.70 0.60];

fault_range = cfg.fault_range;
numSteps    = numel(fault_range);
threshold   = sqrt(-2 * log(1 - cfg.confLevel));   % = sqrt(chi2inv(conf, 2))

conds(1).name = 'Aleatory only';         conds(1).useTol = false; conds(1).file = 'rlc_fault_simulation_results.mat';
conds(2).name = 'Aleatory + Epistemic';  conds(2).useTol = true;  conds(2).file = 'rlc_epistemic_results.mat';

detRate = zeros(numSteps, 2);
FAR     = zeros(1, 2);
minIdx  = nan(1, 2);
aucAll  = zeros(numel(rocSev), 2);
rocFPR  = cell(numel(rocSev), 2);
rocTPR  = cell(numel(rocSev), 2);

for c = 1:2
    fprintf('\n===== %s =====\n', conds(c).name);
    S = load(fullfile(cfg.resultsDir, conds(c).file), 'dataMatrix', 'coeff');
    muData  = mean(S.dataMatrix, 1);
    W       = S.coeff(:, 1:2);
    project = @(Y) (Y - muData) * W;

    % ---- fit healthy baseline on TRAINING trials ----
    Ptrain = project(simulate_trials(cfg, 1.00, nTrain, conds(c).useTol));
    mu     = mean(Ptrain, 1);
    Sigma  = cov(Ptrain);
    mahal  = @(P) sqrt(sum(((P - mu) / Sigma) .* (P - mu), 2));

    % ---- false alarm rate on FRESH healthy trials ----
    dH     = mahal(project(simulate_trials(cfg, 1.00, nTestHealthy, conds(c).useTol)));
    FAR(c) = mean(dH > threshold);
    fprintf('False alarm rate (fresh healthy trials): %.1f%%  (nominal ~%.0f%%)\n', ...
        FAR(c)*100, (1-cfg.confLevel)*100);

    % ---- detection rate vs. severity ----
    dFault = cell(numSteps, 1);
    for s = 1:numSteps
        dFault{s}    = mahal(project(simulate_trials(cfg, fault_range(s), nTestFault, conds(c).useTol)));
        detRate(s,c) = mean(dFault{s} > threshold);
    end

    % ---- minimum detectable fault ----
    bad = find(detRate(:,c) < detTarget, 1, 'last');   % last severity that fails
    if FAR(c) > farMax
        fprintf('False alarm rate exceeds %.0f%% - minimum detectable fault not valid.\n', farMax*100);
    elseif isempty(bad)
        minIdx(c) = 1;
    elseif bad < numSteps
        minIdx(c) = bad + 1;
    end
    if ~isnan(minIdx(c))
        fprintf('Minimum detectable fault: capacitance = %.1f%% of nominal (%.1f%% loss)\n', ...
            fault_range(minIdx(c))*100, (1 - fault_range(minIdx(c)))*100);
    else
        fprintf('No severity in the tested range reaches %.0f%% detection.\n', detTarget*100);
    end

    % ---- ROC / AUC at selected severities ----
    for k = 1:numel(rocSev)
        [~, idx] = min(abs(fault_range - rocSev(k)));
        labels   = [zeros(nTestHealthy,1); ones(nTestFault,1)];
        scores   = [dH; dFault{idx}];
        [rocFPR{k,c}, rocTPR{k,c}, ~, aucAll(k,c)] = perfcurve(labels, scores, 1);
        fprintf('AUC at %.0f%% capacitance: %.3f\n', fault_range(idx)*100, aucAll(k,c));
    end
end

%% ---- Table: severity vs detection rate vs false alarm rate ----
sevCol = zeros(numel(tableSev), 1);
detA = zeros(numel(tableSev), 1);
detB = zeros(numel(tableSev), 1);
for k = 1:numel(tableSev)
    [~, idx] = min(abs(fault_range - tableSev(k)));
    sevCol(k) = fault_range(idx);
    detA(k) = detRate(idx, 1);
    detB(k) = detRate(idx, 2);
end
T = table(sevCol, detA*100, detB*100, ...
    'VariableNames', {'CapFraction', 'DetRate_AleatoryOnly_pct', 'DetRate_Combined_pct'});
fprintf('\nFalse alarm rate: %.1f%% (aleatory only), %.1f%% (combined)\n', FAR(1)*100, FAR(2)*100);
disp(T);
writetable(T, fullfile(cfg.resultsDir, 'detection_table.csv'));

%% ---- Figure: detection rate vs severity ----
f1 = figure; hold on;
plot(fault_range, detRate(:,1)*100, '-o', 'LineWidth', 1.5, 'DisplayName', 'Aleatory only');
plot(fault_range, detRate(:,2)*100, '-s', 'LineWidth', 1.5, 'DisplayName', 'Aleatory + Epistemic');
yline(detTarget*100, '--k', 'DisplayName', '95% target');
yline(FAR(1)*100, ':b', 'DisplayName', 'False alarm rate (aleatory)');
yline(FAR(2)*100, ':r', 'DisplayName', 'False alarm rate (combined)');
set(gca, 'XDir', 'reverse');
xlabel('Capacitance fraction remaining (1 = healthy)');
ylabel('Detection rate (%)');
title('Detection Rate vs. Fault Severity');
legend('Location', 'southeast'); grid on; ylim([0 100]);
save_fig(f1, 'fig_metrics_detection_rate', cfg);

%% ---- Figure: ROC curves ----
f2 = figure;
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
for c = 1:2
    nexttile; hold on;
    for k = 1:numel(rocSev)
        plot(rocFPR{k,c}, rocTPR{k,c}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%.0f%% cap. (AUC %.3f)', rocSev(k)*100, aucAll(k,c)));
    end
    plot([0 1], [0 1], '--k', 'DisplayName', 'Chance');
    xlabel('False positive rate'); ylabel('True positive rate');
    title(conds(c).name); legend('Location', 'southeast'); grid on; axis square;
end
title(tl, 'ROC Curves at Selected Fault Severities');
save_fig(f2, 'fig_metrics_roc', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'detection_metrics.mat'), ...
    'detRate', 'FAR', 'minIdx', 'aucAll', 'rocSev', 'fault_range', 'threshold', ...
    'nTrain', 'nTestHealthy', 'nTestFault', 'T', 'cfg');
fprintf('\nDetection metrics complete. Saved to %s\n', fullfile(cfg.resultsDir, 'detection_metrics.mat'));
