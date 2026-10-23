%% selfbaseline_detection.m  (idea: compare each unit with ITS OWN healthy record)
%
% PURPOSE : Every detector so far compares a measurement with a generic
%           reference built from many units (the fleet). Component tolerance then
%           looks like a fault, because a healthy unit with an unlucky capacitor
%           differs from the fleet reference by about as much as a 10% fault does.
%
%           If instead each unit is compared with its OWN healthy response,
%           recorded once when it was known to be good, its fixed tolerance is in
%           both signals and cancels. What is left is the fault plus noise.
%           This script tests how much of the tolerance problem that removes.
%
%   Protocol per unit (one "unit" = one random draw of R, L, C within tolerance):
%     1. Reference: record the unit's healthy step response m times, average.
%     2. Test: record the unit again at some capacity loss, K times, average.
%     3. Difference d = test - reference. Score it.
%
%   Detectors (all one-sided, thresholds from separate healthy calibration units,
%   target false alarm rate 1%; every rate measured on fresh units):
%     SB-MF   Self-baseline matched filter. Score = <d, g> / |g|, where g is the
%             model-based fault signature (derivative of the nominal response
%             with respect to capacity loss). Needs the circuit equations,
%             needs NO fault examples.
%     SB-E    Self-baseline difference energy. Score = mean(d^2). Model-free.
%     Fleet-MF  CONTROL. The same matched filter, but the reference is the
%             nominal healthy response (a fleet reference), not the unit's own.
%             Isolates exactly what the self-baseline buys.
%
%   Configurations (m reference repeats, K test repeats averaged):
%     A  m = 1,  K = 1   single shot, one noisy reference
%     B  m = 10, K = 1   averaged reference
%     C  m = 10, K = 5   averaged reference, five test measurements averaged
%   Averaging K repeated measurements is what a CUSUM/EWMA accumulates for a
%   constant fault. Tracking a fault that grows over time is not simulated here.
%
%   Sensor sets are the same five as multipoint_sensing.m; each channel has its
%   own noise at cfg.snr_db and is divided by its nominal healthy RMS.
%   Conditions: noise only (all units identical to nominal) and noise +
%   tolerance (each unit drawn from the tolerance bands).
%
% INPUTS  : build_rlc_ladder.m, get_params.m, save_fig.m
% OUTPUTS : results/selfbaseline_detection.mat
%           results/selfbaseline_summary.csv
%           results/figures/fig_selfbaseline_detection_rate.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N   = cfg.N;
nS  = 2*N;
nT  = numel(cfg.t);

cfgMK        = [1 1; 10 1; 10 5];       % [m K] for configurations A, B, C
cfgNames     = {'A: m=1, K=1', 'B: m=10, K=1', 'C: m=10, K=5'};
nCfg         = size(cfgMK, 1);
mMax         = max(cfgMK(:, 1));
KMax         = max(cfgMK(:, 2));
nCal         = 1000;                    % healthy calibration units (thresholds)
nUnits       = 500;                     % fresh test units
sevGrid      = 1 - (0:0.02:0.40);       % 21 severities, 1 = healthy
nSev         = numel(sevGrid);
detTarget    = 0.95;
farMax       = 0.05;
pctl         = 100 * cfg.confLevel;
physDelta    = 0.05;
tableSev     = [0.90 0.80 0.70 0.60];

sensors(1).name = 'V5 (last cap, baseline)';   sensors(1).idx = 2*N;
sensors(2).name = 'V3 (faulted stage)';        sensors(2).idx = 2*cfg.fault_stage;
sensors(3).name = 'V1,V3,V5';                  sensors(3).idx = 2*unique([1 cfg.fault_stage N]);
sensors(4).name = 'All 5 caps';                sensors(4).idx = 2*(1:N);
sensors(5).name = 'All caps + currents';       sensors(5).idx = 1:nS;
nSens = numel(sensors);

detNames = {'SB-MF (self-baseline matched filter)', 'SB-E (self-baseline energy)', ...
            'Fleet-MF (control: fleet reference)'};
nD = numel(detNames);

conds(1).name = 'Aleatory only';        conds(1).useTol = false;
conds(2).name = 'Aleatory + Epistemic'; conds(2).useTol = true;

% ---- nominal reference, fault signature and sensor-set masks ----
yNom   = sim_clean(cfg, cfg.R_nom, cfg.L_nom, cfg.C_nom);            % nT x nS
refRMS = sqrt(mean(yNom.^2, 1));
yNomN  = yNom ./ refRMS;
Cp = cfg.C_nom;  Cp(cfg.fault_stage) = Cp(cfg.fault_stage) * (1 - physDelta);
yPhys  = (sim_clean(cfg, cfg.R_nom, cfg.L_nom, Cp) - yNom) / physDelta ./ refRMS;   % per unit loss
Msk = zeros(nS, nSens);
Gm  = zeros(nT*nS, nSens);
for k = 1:nSens
    Msk(sensors(k).idx, k) = 1;
    g = yPhys .* Msk(:, k).';                      % zero outside the sensor set
    Gm(:, k) = g(:) / norm(g(:));
end
nChan = sum(Msk, 1).';                              % channels per sensor set

detRate = zeros(nSev, nSens, nCfg, nD, 2);
FAR     = zeros(nSens, nCfg, nD, 2);
minIdx  = nan(nSens, nCfg, nD, 2);
auc90   = zeros(nSens, nCfg, nD, 2);
[~, s90] = min(abs(sevGrid - 0.90));

for c = 1:2
    useTol = conds(c).useTol;
    fprintf('\n===== %s =====\n', conds(c).name);
    t0 = tic;

    % ---- cache of clean responses when all units are identical ----
    if ~useTol
        cacheClean = zeros(nT, nS, nSev);
        for s = 1:nSev
            cacheClean(:, :, s) = clean_unit(cfg, draw_unit(cfg, false), sevGrid(s));
        end
    end

    % ---- calibration: healthy units against their own references ----
    scoreCal = zeros(nCal, nSens, nCfg, nD);
    for u = 1:nCal
        un  = draw_unit(cfg, useTol);
        ycH = clean_unit(cfg, un, 1.0);
        scoreCal(u, :, :, :) = unit_scores(cfg, ycH, ycH, refRMS, yNomN, Gm, Msk, nChan, cfgMK, mMax, KMax);
    end
    thr = zeros(nSens, nCfg, nD);
    for k = 1:nSens
        for q = 1:nCfg
            for d = 1:nD
                thr(k, q, d) = prctile(scoreCal(:, k, q, d), pctl);
            end
        end
    end
    clear scoreCal

    % ---- fresh test units across all severities ----
    scoreTest = zeros(nUnits, nSev, nSens, nCfg, nD);
    for u = 1:nUnits
        un  = draw_unit(cfg, useTol);
        if useTol, ycH = clean_unit(cfg, un, 1.0); else, ycH = cacheClean(:, :, 1); end
        % one fixed set of reference repeats per unit, shared by all severities
        refRep = zeros(nT, nS, mMax);
        for r = 1:mMax, refRep(:, :, r) = noisy(ycH, cfg.snr_db, refRMS); end
        refM = zeros(nT, nS, nCfg);
        for q = 1:nCfg, refM(:, :, q) = mean(refRep(:, :, 1:cfgMK(q, 1)), 3); end
        for s = 1:nSev
            if useTol, yc = clean_unit(cfg, un, sevGrid(s)); else, yc = cacheClean(:, :, s); end
            tRep = zeros(nT, nS, KMax);
            for r = 1:KMax, tRep(:, :, r) = noisy(yc, cfg.snr_db, refRMS); end
            for q = 1:nCfg
                Ybar = mean(tRep(:, :, 1:cfgMK(q, 2)), 3);
                scoreTest(u, s, :, q, :) = score_meas(Ybar, refM(:, :, q), yNomN, Gm, Msk, nChan);
            end
        end
    end

    % ---- metrics ----
    for k = 1:nSens
        for q = 1:nCfg
            for d = 1:nD
                alarm = squeeze(scoreTest(:, :, k, q, d)) > thr(k, q, d);     % nUnits x nSev
                detRate(:, k, q, d, c) = mean(alarm, 1).';
                FAR(k, q, d, c) = detRate(1, k, q, d, c);
                minIdx(k, q, d, c) = min_detect(detRate(:, k, q, d, c), FAR(k, q, d, c), detTarget, farMax);
                labels = [zeros(nUnits, 1); ones(nUnits, 1)];
                [~, ~, ~, auc90(k, q, d, c)] = perfcurve(labels, ...
                    [scoreTest(:, 1, k, q, d); scoreTest(:, s90, k, q, d)], 1);
            end
        end
    end
    clear scoreTest
    fprintf('  (%.0f s)\n', toc(t0));

    % ---- console summary ----
    fprintf('%-38s %-26s %-14s %6s %24s %9s %9s\n', 'Detector', 'Sensors', 'Config', 'FAR %', ...
            'Min detectable fault', 'AUC@90%', 'Det@90%');
    for d = 1:nD
        for k = 1:nSens
            for q = 1:nCfg
                if isnan(minIdx(k, q, d, c))
                    mdTxt = 'not reached';
                else
                    mdTxt = sprintf('%.1f%% cap (%.1f%% loss)', sevGrid(minIdx(k, q, d, c))*100, ...
                                    (1 - sevGrid(minIdx(k, q, d, c)))*100);
                end
                fprintf('%-38s %-26s %-14s %6.1f %24s %9.3f %9.0f\n', detNames{d}, sensors(k).name, ...
                        cfgNames{q}, FAR(k, q, d, c)*100, mdTxt, auc90(k, q, d, c), ...
                        detRate(s90, k, q, d, c)*100);
            end
        end
    end
end

%% ---- Table ----
rowsCond = {}; rowsSens = {}; rowsCfg = {}; rowsDet = {};
farCol = []; minCol = []; aucCol = []; dCols = zeros(0, numel(tableSev));
for c = 1:2
    for d = 1:nD
        for k = 1:nSens
            for q = 1:nCfg
                rowsCond{end+1, 1} = conds(c).name;    %#ok<SAGROW>
                rowsSens{end+1, 1} = sensors(k).name;  %#ok<SAGROW>
                rowsCfg{end+1, 1}  = cfgNames{q};      %#ok<SAGROW>
                rowsDet{end+1, 1}  = detNames{d};      %#ok<SAGROW>
                farCol(end+1, 1)   = FAR(k, q, d, c)*100; %#ok<SAGROW>
                if isnan(minIdx(k, q, d, c)), minCol(end+1, 1) = NaN; else, minCol(end+1, 1) = sevGrid(minIdx(k, q, d, c)); end %#ok<SAGROW>
                aucCol(end+1, 1)   = auc90(k, q, d, c);   %#ok<SAGROW>
                row = zeros(1, numel(tableSev));
                for r = 1:numel(tableSev)
                    [~, ii] = min(abs(sevGrid - tableSev(r)));
                    row(r) = detRate(ii, k, q, d, c)*100;
                end
                dCols(end+1, :) = row; %#ok<SAGROW>
            end
        end
    end
end
T_sum = table(rowsCond, rowsSens, rowsCfg, rowsDet, farCol, minCol, aucCol, ...
    dCols(:,1), dCols(:,2), dCols(:,3), dCols(:,4), 'VariableNames', ...
    {'Condition', 'Sensors', 'Config', 'Detector', 'FAR_pct', 'MinDetectable_CapFraction', ...
     'AUC_at_90pct_cap', 'Det_at_90pct_cap', 'Det_at_80pct_cap', 'Det_at_70pct_cap', 'Det_at_60pct_cap'});
writetable(T_sum, fullfile(cfg.resultsDir, 'selfbaseline_summary.csv'));

%% ---- Figure: rows = condition, columns = configuration; S1 and S5 ----
f1 = figure; set(f1, 'Position', [100 100 1350 650]);
tl = tiledlayout(2, nCfg, 'TileSpacing', 'compact');
dcol = lines(nD);
showSens = [1 5];  lst = {'--', '-'};
for c = 1:2
    for q = 1:nCfg
        nexttile; hold on;
        for d = 1:nD
            for j = 1:numel(showSens)
                k = showSens(j);
                plot(sevGrid, detRate(:, k, q, d, c)*100, lst{j}, 'Color', dcol(d, :), ...
                     'LineWidth', 1.3, 'DisplayName', sprintf('%s | %s', detNames{d}, sensors(k).name));
            end
        end
        yline(detTarget*100, '--', '95% target', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
        set(gca, 'XDir', 'reverse'); grid on; ylim([0 100]);
        xlabel('Capacitance fraction remaining'); ylabel('Detection rate (%)');
        title(sprintf('%s | %s', conds(c).name, cfgNames{q}), 'FontSize', 9);
        if c == 1 && q == 1, legend('Location', 'southwest', 'FontSize', 6); end
    end
end
title(tl, 'Self-baseline detection (dashed: last capacitor only, solid: all ten signals)');
save_fig(f1, 'fig_selfbaseline_detection_rate', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'selfbaseline_detection.mat'), ...
    'detRate', 'FAR', 'minIdx', 'auc90', 'sevGrid', 'cfgMK', 'cfgNames', 'detNames', ...
    'sensors', 'nUnits', 'nCal', 'T_sum', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'selfbaseline_detection.mat'));

%% ================= local functions =================
function yc = sim_clean_states(cfg, R, L, C)
    [A, B, ~, ~] = build_rlc_ladder(R, L, C, cfg.RL, cfg.Rs);
    nS = size(A, 1);
    yc = lsim(ss(A, B, eye(nS), zeros(nS, 1)), cfg.u, cfg.t);
end

function yc = sim_clean(cfg, R, L, C)
    yc = sim_clean_states(cfg, R, L, C);
end

function un = draw_unit(cfg, useTol)
% One physical unit: component values fixed for the unit's whole life.
    N = cfg.N;
    if useTol
        un.R = cfg.R_nom .* (1 + cfg.R_tol * (2*rand(1,N) - 1));
        un.L = cfg.L_nom .* (1 + cfg.L_tol * (2*rand(1,N) - 1));
        un.C = cfg.C_nom .* (1 + cfg.C_tol * (2*rand(1,N) - 1));
    else
        un.R = cfg.R_nom;  un.L = cfg.L_nom;  un.C = cfg.C_nom;
    end
end

function yc = clean_unit(cfg, un, sev)
% Noise-free response of the unit with the faulted capacitor scaled by sev.
    C = un.C;  C(cfg.fault_stage) = C(cfg.fault_stage) * sev;
    yc = sim_clean_states(cfg, un.R, un.L, C);
end

function Y = noisy(yc, snr_db, refRMS)
% One noisy measurement of all channels, normalised by the nominal healthy RMS.
    noisePower = mean(yc.^2, 1) / (10^(snr_db/10));
    Y = (yc + sqrt(noisePower) .* randn(size(yc))) ./ refRMS;
end

function sc = score_meas(Ybar, ref, yNomN, Gm, Msk, nChan)
% Scores of one averaged test measurement against a unit reference.
% Returns nSens x 3: [self-baseline matched filter, self-baseline energy, fleet matched filter].
    nT = size(Ybar, 1);
    d  = Ybar - ref;
    sMF = (d(:).' * Gm).';
    sE  = (sum(d.^2, 1) * Msk).' ./ (nT * nChan);
    dF  = Ybar - yNomN;
    sF  = (dF(:).' * Gm).';
    sc  = [sMF, sE, sF];
end

function out = unit_scores(cfg, ycRef, ycTest, refRMS, yNomN, Gm, Msk, nChan, cfgMK, mMax, KMax)
% Healthy calibration unit: scores of healthy test measurements against the unit's own
% reference, for every configuration. out: 1 x nSens x nCfg x 3.
    [nT, nS] = size(ycRef);
    refRep = zeros(nT, nS, mMax);   tRep = zeros(nT, nS, KMax);
    for r = 1:mMax, refRep(:, :, r) = noisy(ycRef,  cfg.snr_db, refRMS); end
    for r = 1:KMax, tRep(:, :, r)   = noisy(ycTest, cfg.snr_db, refRMS); end
    nCfg = size(cfgMK, 1);  nSens = size(Gm, 2);
    out = zeros(1, nSens, nCfg, 3);
    for q = 1:nCfg
        ref  = mean(refRep(:, :, 1:cfgMK(q, 1)), 3);
        Ybar = mean(tRep(:, :, 1:cfgMK(q, 2)), 3);
        out(1, :, q, :) = reshape(score_meas(Ybar, ref, yNomN, Gm, Msk, nChan), 1, nSens, 1, 3);
    end
end

function idx = min_detect(rate, far, detTarget, farMax)
    idx = NaN;
    if far > farMax, return; end
    bad = find(rate < detTarget, 1, 'last');
    if isempty(bad)
        idx = 1;
    elseif bad < numel(rate)
        idx = bad + 1;
    end
end
