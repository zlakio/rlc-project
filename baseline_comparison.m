%% baseline_comparison.m  (README roadmap: baseline comparison)
%
% PURPOSE : Does the PCA detector earn its complexity? Compare it, on identical
%           data and with identical metrics, against simpler alternatives:
%
%     1  Swept PCA : T2         the project's current detector (PCA basis fit on
%                               the swept healthy + faulty data)
%     2  Healthy PCA : T2       "plain" PCA: basis fit on healthy trials ONLY,
%                               Mahalanobis distance in the first nPC scores
%     3  Healthy PCA : Q        squared residual outside those nPC components
%     4  Healthy PCA : T2 OR Q  either of the two alarms
%     5  Raw energy             mean(y.^2); alarm when it leaves the healthy band
%     6  Deviation energy       mean((y - healthy mean response).^2)
%
%   Methods 2-4 remove the known limitation of method 1 (its basis has already
%   seen faulty data). Methods 5-6 use no PCA at all: if they match PCA, then
%   PCA is not buying anything.
%
% PROTOCOL: identical to detection_metrics_q.m. Every detector is fit on healthy
%           TRAINING trials, thresholds for methods 3, 5, 6 come from a separate
%           healthy CALIBRATION set, and every reported number is measured on
%           FRESH trials. All scores are divided by their own threshold, so
%           "alarm" always means score > 1.
%
%   Thresholds  T2 (1, 2)  : sqrt(chi2inv(confLevel, nPC))
%               Q  (3)     : empirical confLevel percentile, calibration set
%               Raw energy : two-sided, |E - mean(E_train)|, empirical percentile
%               Dev. energy: empirical confLevel percentile, calibration set
%
% INPUTS  : results/rlc_fault_simulation_results.mat   (swept PCA basis, aleatory)
%           results/rlc_epistemic_results.mat          (swept PCA basis, combined)
%           simulate_trials.m, build_rlc_ladder.m, get_params.m, save_fig.m
% OUTPUTS : results/baseline_comparison.mat
%           results/baseline_summary.csv
%           results/baseline_table_aleatory.csv, results/baseline_table_combined.csv
%           results/figures/fig_baseline_detection_rate.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

nPC          = 2;      % retained principal components (same as the project)
nTrain       = 1000;   % healthy trials: fit PCA / baseline mean + covariance
nCal         = 1000;   % separate healthy trials: set empirical thresholds
nTestHealthy = 1000;   % fresh healthy trials -> false alarm rate
nTestFault   = 100;    % fresh trials per severity -> detection rate
rocSev       = [0.95 0.90 0.80 0.60];
detTarget    = 0.95;
farMax       = 0.05;
tableSev     = [0.98 0.95 0.90 0.85 0.80 0.70 0.60];

fault_range = cfg.fault_range;
numSteps    = numel(fault_range);
thrT        = sqrt(chi2inv(cfg.confLevel, nPC));
pctl        = 100 * cfg.confLevel;

methodNames = {'Swept PCA: T2', 'Healthy PCA: T2', 'Healthy PCA: Q', ...
               'Healthy PCA: T2 OR Q', 'Raw energy', 'Deviation energy'};
nM = numel(methodNames);

conds(1).name = 'Aleatory only';         conds(1).useTol = false; conds(1).file = 'rlc_fault_simulation_results.mat';
conds(2).name = 'Aleatory + Epistemic';  conds(2).useTol = true;  conds(2).file = 'rlc_epistemic_results.mat';

detRate = zeros(numSteps, nM, 2);
FAR     = zeros(nM, 2);
minIdx  = nan(nM, 2);
aucAll  = zeros(numel(rocSev), nM, 2);

for c = 1:2
    fprintf('\n===== %s =====\n', conds(c).name);

    % ---- (1) swept-data PCA basis: the project's existing detector ----
    S = load(fullfile(cfg.resultsDir, conds(c).file), 'dataMatrix', 'coeff');
    m.muS = mean(S.dataMatrix, 1);
    m.WS  = S.coeff(:, 1:nPC);
    clear S

    % ---- healthy training set: healthy-only PCA + every baseline ----
    Ytr    = simulate_trials(cfg, 1.00, nTrain, conds(c).useTol);
    m.muH  = mean(Ytr, 1);
    coeffH = pca(Ytr, 'NumComponents', nPC);
    m.WH   = coeffH(:, 1:nPC);

    P1 = (Ytr - m.muS) * m.WS;   m.mu1 = mean(P1, 1);   m.Sig1 = cov(P1);
    P2 = (Ytr - m.muH) * m.WH;   m.mu2 = mean(P2, 1);   m.Sig2 = cov(P2);
    m.muE  = mean(mean(Ytr.^2, 2));
    m.thr  = [thrT, thrT, 1, 1, 1];          % placeholders, filled from calibration

    % ---- thresholds from a separate healthy calibration set ----
    Ycal   = simulate_trials(cfg, 1.00, nCal, conds(c).useTol);
    rawCal = raw_stats(Ycal, m);
    m.thr(3) = prctile(rawCal(:, 3), pctl);
    m.thr(4) = prctile(rawCal(:, 4), pctl);
    m.thr(5) = prctile(rawCal(:, 5), pctl);
    clear Ytr Ycal rawCal P1 P2

    % ---- false alarm rate on FRESH healthy trials ----
    scH = score_matrix(simulate_trials(cfg, 1.00, nTestHealthy, conds(c).useTol), m);
    FAR(:, c) = mean(scH > 1, 1).';

    % ---- detection rate vs. severity ----
    scF = cell(numSteps, 1);
    for s = 1:numSteps
        scF{s} = score_matrix(simulate_trials(cfg, fault_range(s), nTestFault, conds(c).useTol), m);
        detRate(s, :, c) = mean(scF{s} > 1, 1);
    end

    % ---- minimum detectable fault + AUC ----
    for mm = 1:nM
        minIdx(mm, c) = min_detect(detRate(:, mm, c), FAR(mm, c), detTarget, farMax);
    end
    labels = [zeros(nTestHealthy, 1); ones(nTestFault, 1)];
    for k = 1:numel(rocSev)
        [~, idx] = min(abs(fault_range - rocSev(k)));
        for mm = 1:nM
            [~, ~, ~, aucAll(k, mm, c)] = perfcurve(labels, [scH(:, mm); scF{idx}(:, mm)], 1);
        end
    end

    % ---- console summary ----
    kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
    fprintf('%-24s %8s %24s %10s\n', 'Method', 'FAR %', 'Min detectable fault', 'AUC @90%');
    for mm = 1:nM
        if isnan(minIdx(mm, c))
            mdTxt = 'not reached';
        else
            mdTxt = sprintf('%.1f%% cap (%.1f%% loss)', fault_range(minIdx(mm, c))*100, ...
                            (1 - fault_range(minIdx(mm, c)))*100);
        end
        fprintf('%-24s %8.1f %24s %10.3f\n', methodNames{mm}, FAR(mm, c)*100, mdTxt, aucAll(kRef, mm, c));
    end
end

%% ---- Tables ----
rowsCond = {}; rowsMeth = {};
farCol = []; minCol = []; aucCol = [];
kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
for c = 1:2
    for mm = 1:nM
        rowsCond{end+1, 1} = conds(c).name;      %#ok<SAGROW>
        rowsMeth{end+1, 1} = methodNames{mm};    %#ok<SAGROW>
        farCol(end+1, 1)   = FAR(mm, c) * 100;   %#ok<SAGROW>
        if isnan(minIdx(mm, c)), minCol(end+1, 1) = NaN; else, minCol(end+1, 1) = fault_range(minIdx(mm, c)); end %#ok<SAGROW>
        aucCol(end+1, 1)   = aucAll(kRef, mm, c); %#ok<SAGROW>
    end
end
T_sum = table(rowsCond, rowsMeth, farCol, minCol, aucCol, ...
    'VariableNames', {'Condition', 'Method', 'FAR_pct', 'MinDetectable_CapFraction', 'AUC_at_90pct_cap'});
writetable(T_sum, fullfile(cfg.resultsDir, 'baseline_summary.csv'));

colNames = {'CapFraction', 'SweptPCA_T2', 'HealthyPCA_T2', 'HealthyPCA_Q', ...
            'HealthyPCA_T2orQ', 'RawEnergy', 'DeviationEnergy'};
csvNames = {'baseline_table_aleatory.csv', 'baseline_table_combined.csv'};
T_rate = cell(1, 2);
for c = 1:2
    M = zeros(numel(tableSev), 1 + nM);
    for k = 1:numel(tableSev)
        [~, idx] = min(abs(fault_range - tableSev(k)));
        M(k, 1)     = fault_range(idx);
        M(k, 2:end) = detRate(idx, :, c) * 100;
    end
    T_rate{c} = array2table(M, 'VariableNames', colNames);
    writetable(T_rate{c}, fullfile(cfg.resultsDir, csvNames{c}));
    fprintf('\nDetection rate (%%) by severity - %s:\n', conds(c).name);
    disp(T_rate{c});
end

%% ---- Figure ----
f1 = figure;
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
mk = {'-o', '-s', '-^', '-d', '--x', '--+'};
for c = 1:2
    nexttile; hold on;
    for mm = 1:nM
        plot(fault_range, detRate(:, mm, c)*100, mk{mm}, 'LineWidth', 1.2, ...
             'MarkerSize', 4, 'DisplayName', methodNames{mm});
    end
    yline(detTarget*100, '--', '95% target', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
    set(gca, 'XDir', 'reverse');
    xlabel('Capacitance fraction remaining (1 = healthy)');
    ylabel('Detection rate (%)');
    title(conds(c).name); legend('Location', 'southwest', 'FontSize', 7); grid on; ylim([0 100]);
end
title(tl, 'Baseline Comparison: PCA detectors vs. energy thresholds');
save_fig(f1, 'fig_baseline_detection_rate', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'baseline_comparison.mat'), ...
    'detRate', 'FAR', 'minIdx', 'aucAll', 'rocSev', 'fault_range', 'methodNames', ...
    'nPC', 'nTrain', 'nCal', 'nTestHealthy', 'nTestFault', 'T_sum', 'T_rate', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'baseline_comparison.mat'));

%% ================= local functions =================
function raw = raw_stats(Y, m)
% Unnormalised statistics, one row per trial. Columns:
%   1 T2 distance, swept PCA     2 T2 distance, healthy PCA
%   3 Q, healthy PCA             4 |raw energy - healthy mean energy|
%   5 deviation energy (distance to the healthy mean response)
    P1 = (Y - m.muS) * m.WS;
    d1 = sqrt(sum(((P1 - m.mu1) / m.Sig1) .* (P1 - m.mu1), 2));

    Xc = Y - m.muH;
    P2 = Xc * m.WH;
    d2 = sqrt(sum(((P2 - m.mu2) / m.Sig2) .* (P2 - m.mu2), 2));
    Q  = max(sum(Xc.^2, 2) - sum(P2.^2, 2), 0);

    dE = abs(mean(Y.^2, 2) - m.muE);
    D  = mean(Xc.^2, 2);
    raw = [d1, d2, Q, dE, D];
end

function sc = score_matrix(Y, m)
% Threshold-normalised scores (alarm when > 1), columns in methodNames order.
    s = raw_stats(Y, m) ./ m.thr;
    sc = [s(:, 1), s(:, 2), s(:, 3), max(s(:, 2), s(:, 3)), s(:, 4), s(:, 5)];
end

function idx = min_detect(rate, far, detTarget, farMax)
% First severity from which detection stays >= detTarget (and FAR is acceptable).
    idx = NaN;
    if far > farMax, return; end
    bad = find(rate < detTarget, 1, 'last');
    if isempty(bad)
        idx = 1;
    elseif bad < numel(rate)
        idx = bad + 1;
    end
end
