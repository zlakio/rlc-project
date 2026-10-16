%% snr_tolerance_heatmap.m  (README roadmap: SNR x tolerance heatmap)
%
% PURPOSE : Replace the two assumed operating points (15 dB, tolerance as given)
%           with a map: for every combination of measurement-noise level and
%           component-tolerance level, how small a fault can the detector find?
%
%   Axes   : SNR (dB)  x  tolerance multiplier
%            multiplier 0 = no epistemic uncertainty (noise only)
%            multiplier 1 = the tolerances in get_params.m (R 5%, L 10%, C 10%)
%            multiplier 2 = twice those bands, and so on.
%   Cells  : for each of T2 and T2 OR Q, with the same protocol as
%            detection_metrics_q.m (healthy-trained baseline, fresh test trials):
%              - minimum detectable fault (>= detTarget detection, FAR <= farMax)
%              - detection rate at refSev
%              - false alarm rate
%
% HOW EACH CELL IS BUILT (self-contained, so every cell is a full experiment):
%   1. simulate a sweep (nSweep trials at each severity in sevGrid) and fit the
%      PCA basis on it, as the project does (flip useHealthyPCA to use a basis
%      fit on healthy trials only instead);
%   2. fit the healthy baseline (mean, covariance, Q threshold) on nTrain trials;
%   3. measure FAR on nTestHealthy fresh trials and detection on nTestFault
%      fresh trials per severity.
%
% COST    : ~3,500 simulated trials per cell, 36 cells. Expect on the order of
%           10-20 minutes in MATLAB (cells with tolerance cost the most).
%           Lower nTestFault / nSweep / the grid sizes for a quick look.
% RESOLUTION: sevGrid steps are 2%, so "minimum detectable fault" is quantised
%           to 2% capacitance. Detection rates use nTestFault trials per point.
%
% INPUTS  : get_params.m, simulate_trials.m, build_rlc_ladder.m, save_fig.m
% OUTPUTS : results/snr_tolerance_heatmap.mat, results/snr_tolerance_heatmap.csv
%           results/figures/fig_heatmap_snr_tolerance.png

clear; clc; close all;
cfg0 = get_params();
rng(cfg0.seed, 'twister');

snrList  = [5 10 15 20 25 30];          % dB
tolScale = [0 0.25 0.5 1 1.5 2];        % multiples of the get_params tolerances
p.sevGrid       = 1 - (0:0.02:0.40);    % 1.00 ... 0.60
p.nPC           = 2;
p.nSweep        = 20;                   % sweep trials per severity (PCA basis)
p.nTrain        = 500;                  % healthy trials for the baseline
p.nTestHealthy  = 500;                  % fresh healthy trials -> FAR
p.nTestFault    = 100;                  % fresh trials per severity -> detection
p.detTarget     = 0.95;
p.farMax        = 0.05;
p.useHealthyPCA = false;                % false = project's swept basis
refSev          = 0.90;                 % severity for the detection-rate maps

methodNames = {'T2', 'T2 OR Q'};
nS = numel(snrList);  nT = numel(tolScale);

minLoss = nan(nS, nT, 2);      % minimum detectable capacitance LOSS (%), NaN = not reached
farMap  = zeros(nS, nT, 2);    % false alarm rate (%)
detRef  = zeros(nS, nT, 2);    % detection rate (%) at refSev

[~, iRef] = min(abs(p.sevGrid - refSev));
tAll = tic;
for i = 1:nS
    for j = 1:nT
        tCell = tic;
        cfg = cfg0;
        cfg.snr_db = snrList(i);
        cfg.R_tol  = cfg0.R_tol * tolScale(j);
        cfg.L_tol  = cfg0.L_tol * tolScale(j);
        cfg.C_tol  = cfg0.C_tol * tolScale(j);
        useTol = tolScale(j) > 0;

        res = run_cell(cfg, useTol, p, iRef);
        minLoss(i, j, :) = res.minLoss;
        farMap(i, j, :)  = res.far;
        detRef(i, j, :)  = res.detRef;

        fprintf('SNR %2d dB | tol x%-4g | FAR %4.1f%% / %4.1f%% | min detectable loss: T2 %s, T2 OR Q %s | %.0f s\n', ...
            snrList(i), tolScale(j), res.far(1), res.far(2), ...
            fmt_loss(res.minLoss(1)), fmt_loss(res.minLoss(2)), toc(tCell));
    end
end
fprintf('\nAll cells done in %.1f min.\n', toc(tAll)/60);

%% ---- Figure: 2 x 3 maps ----
f1 = figure('Position', [100 100 1250 700]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact');
nexttile; draw_map(minLoss(:,:,1), snrList, tolScale, 'Min detectable loss (%) - T^2', '%.0f');
nexttile; draw_map(minLoss(:,:,2), snrList, tolScale, 'Min detectable loss (%) - T^2 OR Q', '%.0f');
nexttile; draw_map(farMap(:,:,1),  snrList, tolScale, 'False alarm rate (%) - T^2', '%.1f');
nexttile; draw_map(detRef(:,:,1),  snrList, tolScale, sprintf('Detection rate (%%) at %.0f%% cap. - T^2', refSev*100), '%.0f');
nexttile; draw_map(detRef(:,:,2),  snrList, tolScale, sprintf('Detection rate (%%) at %.0f%% cap. - T^2 OR Q', refSev*100), '%.0f');
nexttile; draw_map(farMap(:,:,2),  snrList, tolScale, 'False alarm rate (%) - T^2 OR Q', '%.1f');
title(tl, 'Detectability vs. measurement noise and component tolerance  (n.r. = target not reached within 40% loss, or FAR above limit)');
save_fig(f1, 'fig_heatmap_snr_tolerance', cfg0);

%% ---- Long-form CSV ----
rows = nS * nT * 2;
snrCol = zeros(rows, 1); tolCol = zeros(rows, 1); methCol = cell(rows, 1);
farCol = zeros(rows, 1); minCol = zeros(rows, 1); refCol = zeros(rows, 1);
r = 0;
for i = 1:nS
    for j = 1:nT
        for m = 1:2
            r = r + 1;
            snrCol(r) = snrList(i);  tolCol(r) = tolScale(j);  methCol{r} = methodNames{m};
            farCol(r) = farMap(i, j, m);  minCol(r) = minLoss(i, j, m);  refCol(r) = detRef(i, j, m);
        end
    end
end
T = table(snrCol, tolCol, methCol, farCol, minCol, refCol, 'VariableNames', ...
    {'SNR_dB', 'TolMultiplier', 'Method', 'FAR_pct', 'MinDetectableLoss_pct', 'DetRateAtRefSev_pct'});
writetable(T, fullfile(cfg0.resultsDir, 'snr_tolerance_heatmap.csv'));

save(fullfile(cfg0.resultsDir, 'snr_tolerance_heatmap.mat'), ...
    'snrList', 'tolScale', 'minLoss', 'farMap', 'detRef', 'refSev', 'methodNames', 'p', 'T', 'cfg0');
fprintf('Saved to %s\n', fullfile(cfg0.resultsDir, 'snr_tolerance_heatmap.mat'));

%% ================= local functions =================
function res = run_cell(cfg, useTol, p, iRef)
% One full experiment at (cfg.snr_db, cfg tolerances). Returns, for [T2, T2 OR Q]:
%   res.minLoss (% capacitance loss, NaN = not reached), res.far (%), res.detRef (%).
    nSev = numel(p.sevGrid);
    thrT = sqrt(chi2inv(cfg.confLevel, p.nPC));

    Ytr = simulate_trials(cfg, 1.00, p.nTrain, useTol);
    if p.useHealthyPCA
        muData = mean(Ytr, 1);
        coeff  = pca(Ytr, 'NumComponents', p.nPC);
    else
        Ysw = zeros(nSev * p.nSweep, numel(cfg.t));
        for s = 1:nSev
            Ysw((s-1)*p.nSweep + (1:p.nSweep), :) = simulate_trials(cfg, p.sevGrid(s), p.nSweep, useTol);
        end
        muData = mean(Ysw, 1);
        coeff  = pca(Ysw, 'NumComponents', p.nPC);
    end
    W = coeff(:, 1:p.nPC);

    [Ptr, Qtr] = pc_stats(Ytr, muData, W);
    mu    = mean(Ptr, 1);
    Sigma = cov(Ptr);
    thrQ  = prctile(Qtr, 100 * cfg.confLevel);

    [Ph, Qh] = pc_stats(simulate_trials(cfg, 1.00, p.nTestHealthy, useTol), muData, W);
    scH = norm_scores(Ph, Qh, mu, Sigma, thrT, thrQ);        % columns: T2, T2 OR Q
    far = mean(scH > 1, 1);

    dr = zeros(nSev, 2);
    for s = 1:nSev
        [Pf, Qf] = pc_stats(simulate_trials(cfg, p.sevGrid(s), p.nTestFault, useTol), muData, W);
        dr(s, :) = mean(norm_scores(Pf, Qf, mu, Sigma, thrT, thrQ) > 1, 1);
    end

    res.far    = far * 100;
    res.detRef = dr(iRef, :) * 100;
    res.minLoss = nan(1, 2);
    for m = 1:2
        idx = min_detect(dr(:, m), far(m), p.detTarget, p.farMax);
        if ~isnan(idx), res.minLoss(m) = (1 - p.sevGrid(idx)) * 100; end
    end
end

function [P, Q] = pc_stats(Y, muData, W)
% Scores in the retained PCs, and Q = squared residual outside them.
    Xc = Y - muData;
    P  = Xc * W;
    Q  = max(sum(Xc.^2, 2) - sum(P.^2, 2), 0);
end

function sc = norm_scores(P, Q, mu, Sigma, thrT, thrQ)
% Columns: T2 distance / thrT, and max(T2/thrT, Q/thrQ) (= "T2 OR Q"). Alarm when > 1.
    d  = sqrt(sum(((P - mu) / Sigma) .* (P - mu), 2));
    sT = d / thrT;
    sQ = Q / thrQ;
    sc = [sT, max(sT, sQ)];
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

function s = fmt_loss(x)
    if isnan(x), s = 'n.r.'; else, s = sprintf('%.0f%%', x); end
end

function draw_map(Z, snrList, tolScale, ttl, fmt)
% Annotated heat map: rows = SNR (low at bottom), columns = tolerance multiplier.
    h = imagesc(Z);
    set(h, 'AlphaData', ~isnan(Z));
    set(gca, 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
    colormap(gca, parula);
    colorbar;
    set(gca, 'XTick', 1:numel(tolScale), 'XTickLabel', arrayfun(@(v) sprintf('%g', v), tolScale, 'UniformOutput', false));
    set(gca, 'YTick', 1:numel(snrList),  'YTickLabel', arrayfun(@(v) sprintf('%g', v), snrList,  'UniformOutput', false));
    xlabel('Tolerance multiplier (1 = R \pm5%, L \pm10%, C \pm10%)');
    ylabel('SNR (dB)');
    title(ttl);
    zmin = min(Z(:), [], 'omitnan');  zmax = max(Z(:), [], 'omitnan');
    if isempty(zmin) || isnan(zmin), zmin = 0; zmax = 1; end
    for r = 1:size(Z, 1)
        for c = 1:size(Z, 2)
            if isnan(Z(r, c))
                txt = 'n.r.';  col = 'k';
            else
                txt = sprintf(fmt, Z(r, c));
                frac = (Z(r, c) - zmin) / max(zmax - zmin, eps);
                if frac < 0.45, col = 'w'; else, col = 'k'; end
            end
            text(c, r, txt, 'HorizontalAlignment', 'center', 'Color', col, 'FontSize', 8);
        end
    end
end
