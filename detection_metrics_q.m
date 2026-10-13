%% detection_metrics_q.m  (README 2.4)
%
% PURPOSE : Add the Q (squared prediction error, SPE) statistic next to the
%           Mahalanobis / T^2 statistic and report detection for:
%             T^2 alone, Q alone, and T^2 OR Q.
%           Run for aleatory-only AND aleatory + epistemic uncertainty.
%
% WHY Q  : T^2 only sees movement INSIDE the retained principal components.
%          Q = ||x - xhat||^2 measures what falls OUTSIDE them (xhat = the
%          reconstruction from the retained PCs).
%
% THRESHOLDS
%   T^2 : sqrt(chi2inv(confLevel, nPC)), baseline fit on healthy training trials
%   Q   : empirical confLevel percentile of Q on the healthy training trials
%   Each statistic is divided by its own threshold, so "alarm" = value > 1.
%   T^2 OR Q uses max(T^2/thrT, Q/thrQ) > 1, which is exactly "either alarms".
%
% INPUTS  : results/rlc_fault_simulation_results.mat
%           results/rlc_epistemic_results.mat
%           simulate_trials.m, build_rlc_ladder.m, get_params.m, save_fig.m
% OUTPUTS : results/detection_metrics_q.mat, results/detection_table_q.csv
%           results/figures/fig_q_detection_rate.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

nPC          = 2;      % retained principal components (T^2 and Q both use these)
nTrain       = 1000;   % healthy trials used to fit baselines / thresholds
nTestHealthy = 1000;   % fresh healthy trials -> false alarm rate
nTestFault   = 100;    % fresh trials per severity -> detection rate
rocSev       = [0.95 0.90 0.80 0.60];
detTarget    = 0.95;
farMax       = 0.05;
tableSev     = [0.98 0.95 0.90 0.85 0.80 0.70 0.60];

fault_range = cfg.fault_range;
numSteps    = numel(fault_range);
thrT        = sqrt(chi2inv(cfg.confLevel, nPC));
methodNames = {'T2 only', 'Q only', 'T2 OR Q'};

conds(1).name = 'Aleatory only';         conds(1).useTol = false; conds(1).file = 'rlc_fault_simulation_results.mat';
conds(2).name = 'Aleatory + Epistemic';  conds(2).useTol = true;  conds(2).file = 'rlc_epistemic_results.mat';

detRate = zeros(numSteps, 3, 2);
FAR     = zeros(3, 2);
minIdx  = nan(3, 2);
aucAll  = zeros(numel(rocSev), 3, 2);

for c = 1:2
    fprintf('\n===== %s =====\n', conds(c).name);
    S = load(fullfile(cfg.resultsDir, conds(c).file), 'dataMatrix', 'coeff');
    muData = mean(S.dataMatrix, 1);
    W      = S.coeff(:, 1:nPC);

    % ---- fit baseline + thresholds on healthy TRAINING trials ----
    [Ptr, Qtr] = pc_stats(simulate_trials(cfg, 1.00, nTrain, conds(c).useTol), muData, W);
    mu    = mean(Ptr, 1);
    Sigma = cov(Ptr);
    thrQ  = prctile(Qtr, 100 * cfg.confLevel);

    % ---- false alarm rate on FRESH healthy trials ----
    [Ph, Qh] = pc_stats(simulate_trials(cfg, 1.00, nTestHealthy, conds(c).useTol), muData, W);
    scH = norm_scores(Ph, Qh, mu, Sigma, thrT, thrQ);
    FAR(:, c) = mean(scH > 1, 1).';
    fprintf('False alarm rate:  T2 %.1f%%   Q %.1f%%   T2 OR Q %.1f%%\n', FAR(:,c)*100);

    % ---- detection rate vs. severity ----
    scF = cell(numSteps, 1);
    for s = 1:numSteps
        [Pf, Qf] = pc_stats(simulate_trials(cfg, fault_range(s), nTestFault, conds(c).useTol), muData, W);
        scF{s} = norm_scores(Pf, Qf, mu, Sigma, thrT, thrQ);
        detRate(s, :, c) = mean(scF{s} > 1, 1);
    end

    % ---- minimum detectable fault + AUC ----
    for m = 1:3
        minIdx(m, c) = min_detect(detRate(:, m, c), FAR(m, c), detTarget, farMax);
        if isnan(minIdx(m, c))
            fprintf('Min detectable fault (%s): not reached within tested range\n', methodNames{m});
        else
            fprintf('Min detectable fault (%s): %.1f%% capacitance (%.1f%% loss)\n', methodNames{m}, ...
                fault_range(minIdx(m,c))*100, (1 - fault_range(minIdx(m,c)))*100);
        end
    end
    for k = 1:numel(rocSev)
        [~, idx] = min(abs(fault_range - rocSev(k)));
        labels = [zeros(nTestHealthy,1); ones(nTestFault,1)];
        for m = 1:3
            [~, ~, ~, aucAll(k, m, c)] = perfcurve(labels, [scH(:,m); scF{idx}(:,m)], 1);
        end
        fprintf('AUC at %.0f%% cap.:  T2 %.3f   Q %.3f   T2 OR Q %.3f\n', ...
            fault_range(idx)*100, aucAll(k,1,c), aucAll(k,2,c), aucAll(k,3,c));
    end
end

%% ---- Table: detection rate (%) by severity ----
M = zeros(numel(tableSev), 7);
for k = 1:numel(tableSev)
    [~, idx] = min(abs(fault_range - tableSev(k)));
    M(k, 1)   = fault_range(idx);
    M(k, 2:4) = detRate(idx, :, 1) * 100;
    M(k, 5:7) = detRate(idx, :, 2) * 100;
end
T = array2table(M, 'VariableNames', {'CapFraction', ...
    'Aleat_T2', 'Aleat_Q', 'Aleat_OR', 'Comb_T2', 'Comb_Q', 'Comb_OR'});
fprintf('\nDetection rate (%%) by severity:\n');
disp(T);
writetable(T, fullfile(cfg.resultsDir, 'detection_table_q.csv'));

%% ---- Figure ----
f1 = figure;
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
mk = {'-o', '-s', '-^'};
for c = 1:2
    nexttile; hold on;
    for m = 1:3
        plot(fault_range, detRate(:, m, c)*100, mk{m}, 'LineWidth', 1.3, 'DisplayName', methodNames{m});
    end
    yline(detTarget*100, '--', '95% target', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
    set(gca, 'XDir', 'reverse');
    xlabel('Capacitance fraction remaining (1 = healthy)');
    ylabel('Detection rate (%)');
    title(conds(c).name); legend('Location', 'best'); grid on; ylim([0 100]);
end
title(tl, 'Detection Rate vs. Severity: T^2, Q, and T^2 OR Q');
save_fig(f1, 'fig_q_detection_rate', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'detection_metrics_q.mat'), ...
    'detRate', 'FAR', 'minIdx', 'aucAll', 'rocSev', 'fault_range', ...
    'thrT', 'nPC', 'nTrain', 'nTestHealthy', 'nTestFault', 'T', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'detection_metrics_q.mat'));

%% ================= local functions =================
function [P, Q] = pc_stats(Y, muData, W)
% Scores in the retained PCs, and Q = squared residual outside them.
    Xc = Y - muData;
    P  = Xc * W;
    Q  = max(sum(Xc.^2, 2) - sum(P.^2, 2), 0);
end

function sc = norm_scores(P, Q, mu, Sigma, thrT, thrQ)
% Columns: T^2 distance / thrT, Q / thrQ, and max of the two. Alarm when > 1.
    d  = sqrt(sum(((P - mu) / Sigma) .* (P - mu), 2));
    sT = d / thrT;
    sQ = Q / thrQ;
    sc = [sT, sQ, max(sT, sQ)];
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