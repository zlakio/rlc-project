%% nuisance_aware_features.m  (next-step fix 4: stop wasting detection power on tolerance)
%
% PURPOSE : With component tolerance, the project's swept-PCA T2 detector keeps
%           the two directions of largest VARIANCE, and under tolerance those
%           mostly describe part-to-part spread (the "nuisance"), not the fault.
%           multipoint_sensing.m showed the symptom: adding capacitor channels
%           did not help, and 10 PCs did worse than 2.
%
%           This script asks: if the detector is told about the nuisance, how
%           much of the lost performance comes back, on the SAME sensor sets as
%           multipoint_sensing.m?
%
%   Detectors (both share one 60-component PCA basis fit on the swept data):
%     D1  Swept PCA T2            the project's detector (first 2 PCs). Reference.
%     D2  Fisher matched filter   SUPERVISED. Fault signature g = slope of the
%                                 data against capacity loss (regression on the
%                                 training sweep); weight w = Sw^-1 g, where Sw
%                                 is the covariance of HEALTHY trials, i.e. of
%                                 noise + tolerance. Directions where healthy
%                                 trials already vary a lot are down-weighted.
%                                 Threshold: 99th percentile of a separate
%                                 healthy calibration set.
%
%   Note: an unsupervised variant that deleted the top 10 tolerance directions
%   (no labels) was tried in an earlier version of this script and removed. It
%   collapsed detection to the false alarm rate, because the faulted capacitor's
%   own tolerance lies along the fault direction and the projection removes the
%   signal with it. Knowing the fault direction (labels) is what makes D2 work.
%
%   Interpretability check (is the learned direction still physical?):
%     physics signature = d(response)/d(capacity loss) of the nominal circuit,
%     from the model, no data. Reported as cosine similarity with
%       - the regression signature g (what D2 learned)            cosG
%       - the filter itself, w mapped back to the signal domain   cosW
%       - PC1 of the swept PCA (what D1 uses; sign-agnostic)      cosPC1
%     1 = same shape as the physical fault signature, 0 = unrelated.
%
% PROTOCOL: as in multipoint_sensing.m. Healthy training set for score
%           covariances, separate healthy calibration set (D2 threshold), every
%           reported rate on FRESH trials (1000 healthy, 100 per severity).
%
% INPUTS  : build_rlc_ladder.m, get_params.m, save_fig.m
% OUTPUTS : results/nuisance_aware_features.mat
%           results/nuisance_summary.csv
%           results/nuisance_physics_cosine.csv
%           results/figures/fig_nuisance_detection_rate.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N   = cfg.N;
nS  = 2*N;
nT  = numel(cfg.t);

mComp        = 60;      % PCA components kept for D2 (score space)
ridgeFrac    = 1e-3;    % ridge added to healthy covariance, as a fraction of its mean eigenvalue
nPC          = 2;       % T2 capacity of D1 (same as the project)
nTrainHealthy = 1000;
nCal         = 1000;
nTestHealthy = 1000;
nTestFault   = 100;
rocSev       = [0.95 0.90 0.80 0.60];
detTarget    = 0.95;
farMax       = 0.05;
tableSev     = [0.95 0.90 0.80 0.70 0.60];
physDelta    = 0.05;    % capacity loss used for the model-based fault signature

fault_range = cfg.fault_range;
numSteps    = numel(fault_range);
thrT        = sqrt(chi2inv(cfg.confLevel, nPC));
pctl        = 100 * cfg.confLevel;

sensors(1).name = 'V5 (last cap, baseline)';   sensors(1).idx = 2*N;
sensors(2).name = 'V3 (faulted stage)';        sensors(2).idx = 2*cfg.fault_stage;
sensors(3).name = 'V1,V3,V5';                  sensors(3).idx = 2*unique([1 cfg.fault_stage N]);
sensors(4).name = 'All 5 caps';                sensors(4).idx = 2*(1:N);
sensors(5).name = 'All caps + currents';       sensors(5).idx = 1:nS;
nSens = numel(sensors);

detNames = {'D1 Swept PCA T2', 'D2 Fisher matched filter'};
nD = numel(detNames);

conds(1).name = 'Aleatory only';        conds(1).useTol = false;
conds(2).name = 'Aleatory + Epistemic'; conds(2).useTol = true;

% ---- nominal reference RMS and the model-based fault signature ----
yNom   = sim_clean(cfg, cfg.R_nom, cfg.L_nom, cfg.C_nom);          % nT x nS
refRMS = sqrt(mean(yNom.^2, 1));
Cp = cfg.C_nom;  Cp(cfg.fault_stage) = Cp(cfg.fault_stage) * (1 - physDelta);
yPhys  = (sim_clean(cfg, cfg.R_nom, cfg.L_nom, Cp) - yNom) / physDelta ./ refRMS;   % nT x nS

detRate = zeros(numSteps, nSens, nD, 2);
FAR     = zeros(nSens, nD, 2);
minIdx  = nan(nSens, nD, 2);
aucAll  = zeros(numel(rocSev), nSens, nD, 2);
cosG    = zeros(nSens, 2);  cosW = zeros(nSens, 2);  cosPC1 = zeros(nSens, 2);

for c = 1:2
    useTol = conds(c).useTol;
    fprintf('\n===== %s =====\n', conds(c).name);

    tps    = cfg.trials_per_step;
    nSw    = numSteps*tps;
    Zsweep = zeros(nSw, nS, nT, 'single');
    delta  = zeros(nSw, 1);
    for s = 1:numSteps
        Zsweep((s-1)*tps+1 : s*tps, :, :) = gen_states(cfg, fault_range(s), tps, useTol, refRMS);
        delta((s-1)*tps+1 : s*tps) = 1 - fault_range(s);
    end
    Zh = gen_states(cfg, 1.00, nTrainHealthy, useTol, refRMS);

    % ---- fit the shared basis and the two detectors, per sensor set ----
    models = cell(nSens, 1);
    for k = 1:nSens
        idx = sensors(k).idx;
        B   = fit_basis(Zsweep, idx, mComp);
        Ps  = project(Zsweep, idx, B);
        Ph  = project(Zh,     idx, B);
        mdl.B = B;
        % D1: first nPC PCs, healthy score statistics
        mdl.mu1  = mean(Ph(:, 1:nPC), 1);  mdl.Sig1 = cov(Ph(:, 1:nPC));
        % D2: regression signature in score space, whitened by healthy covariance
        dc  = delta - mean(delta);
        gS  = ((Ps - mean(Ps, 1)).' * dc) / sum(dc.^2);              % mComp x 1
        Sw  = cov(Ph);
        Sw  = Sw + ridgeFrac * mean(diag(Sw)) * eye(mComp);
        mdl.w   = Sw \ gS;
        mdl.muH = mean(Ph, 1);
        mdl.thr2 = 1;                                                % filled below
        models{k} = mdl;

        % interpretability: cosine with the model-based physical signature
        phys = vec_stack(yPhys(:, idx).');                           % channels x nT -> vector
        gFeat = zeros(numel(idx), nT);
        wFeat = zeros(numel(idx), nT);
        for q = 1:numel(idx)
            Xq = double(reshape(Zsweep(:, idx(q), :), nSw, nT));
            gFeat(q, :) = (dc.' * (Xq - mean(Xq, 1))) / sum(dc.^2);
            wFeat(q, :) = (B.W{q} * mdl.w).';
        end
        pc1 = zeros(numel(idx), nT);
        for q = 1:numel(idx), pc1(q, :) = B.W{q}(:, 1).'; end
        cosG(k, c)   = cosine(gFeat(:), phys);
        cosW(k, c)   = cosine(wFeat(:), phys);
        cosPC1(k, c) = abs(cosine(pc1(:), phys));
        clear Xq gFeat wFeat pc1
    end
    clear Zsweep Zh

    % ---- thresholds for D2 from a separate healthy calibration set ----
    Zc = gen_states(cfg, 1.00, nCal, useTol, refRMS);
    for k = 1:nSens
        Pc = project(Zc, sensors(k).idx, models{k}.B);
        sc2 = (Pc - models{k}.muH) * models{k}.w;
        models{k}.thr2 = prctile(sc2, pctl);
    end
    clear Zc

    % ---- false alarm rate on FRESH healthy trials ----
    Zt  = gen_states(cfg, 1.00, nTestHealthy, useTol, refRMS);
    scH = cell(nSens, 1);
    for k = 1:nSens
        scH{k} = score_all(project(Zt, sensors(k).idx, models{k}.B), models{k}, thrT);
        FAR(k, :, c) = mean(scH{k} > 1, 1);
    end
    clear Zt

    % ---- detection rate vs. severity ----
    scF = cell(numel(rocSev), 1);
    rocIdx = zeros(1, numel(rocSev));
    for r = 1:numel(rocSev), [~, rocIdx(r)] = min(abs(fault_range - rocSev(r))); end
    for s = 1:numSteps
        Zf = gen_states(cfg, fault_range(s), nTestFault, useTol, refRMS);
        sc = cell(nSens, 1);
        for k = 1:nSens
            sc{k} = score_all(project(Zf, sensors(k).idx, models{k}.B), models{k}, thrT);
            detRate(s, k, :, c) = reshape(mean(sc{k} > 1, 1), 1, 1, nD);
        end
        kr = find(rocIdx == s);
        if ~isempty(kr), scF{kr(1)} = sc; end
        clear Zf
    end

    % ---- minimum detectable fault + AUC ----
    labels = [zeros(nTestHealthy, 1); ones(nTestFault, 1)];
    for k = 1:nSens
        for d = 1:nD
            minIdx(k, d, c) = min_detect(detRate(:, k, d, c), FAR(k, d, c), detTarget, farMax);
            for r = 1:numel(rocSev)
                [~, ~, ~, aucAll(r, k, d, c)] = perfcurve(labels, [scH{k}(:, d); scF{r}{k}(:, d)], 1);
            end
        end
    end

    % ---- console summary ----
    kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
    [~, sRef] = min(abs(fault_range - 0.90));
    fprintf('%-24s %-26s %6s %24s %9s %9s\n', 'Detector', 'Sensors', 'FAR %', ...
            'Min detectable fault', 'AUC@90%', 'Det@90%');
    for d = 1:nD
        for k = 1:nSens
            if isnan(minIdx(k, d, c))
                mdTxt = 'not reached';
            else
                mdTxt = sprintf('%.1f%% cap (%.1f%% loss)', fault_range(minIdx(k, d, c))*100, ...
                                (1 - fault_range(minIdx(k, d, c)))*100);
            end
            fprintf('%-24s %-26s %6.1f %24s %9.3f %9.0f\n', detNames{d}, sensors(k).name, ...
                    FAR(k, d, c)*100, mdTxt, aucAll(kRef, k, d, c), detRate(sRef, k, d, c)*100);
        end
    end
    fprintf('\nCosine with the physical fault signature (1 = same shape):\n');
    fprintf('%-26s %8s %8s %8s\n', 'Sensors', 'cosG', 'cosW', 'cosPC1');
    for k = 1:nSens
        fprintf('%-26s %8.3f %8.3f %8.3f\n', sensors(k).name, cosG(k, c), cosW(k, c), cosPC1(k, c));
    end
end

%% ---- Tables ----
rowsCond = {}; rowsSens = {}; rowsDet = {};
farCol = []; minCol = []; aucCol = []; dCols = zeros(0, numel(tableSev));
kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
for c = 1:2
    for d = 1:nD
        for k = 1:nSens
            rowsCond{end+1, 1} = conds(c).name;    %#ok<SAGROW>
            rowsSens{end+1, 1} = sensors(k).name;  %#ok<SAGROW>
            rowsDet{end+1, 1}  = detNames{d};      %#ok<SAGROW>
            farCol(end+1, 1)   = FAR(k, d, c)*100; %#ok<SAGROW>
            if isnan(minIdx(k, d, c)), minCol(end+1, 1) = NaN; else, minCol(end+1, 1) = fault_range(minIdx(k, d, c)); end %#ok<SAGROW>
            aucCol(end+1, 1)   = aucAll(kRef, k, d, c); %#ok<SAGROW>
            row = zeros(1, numel(tableSev));
            for r = 1:numel(tableSev)
                [~, ii] = min(abs(fault_range - tableSev(r)));
                row(r) = detRate(ii, k, d, c)*100;
            end
            dCols(end+1, :) = row; %#ok<SAGROW>
        end
    end
end
T_sum = table(rowsCond, rowsSens, rowsDet, farCol, minCol, aucCol, ...
    dCols(:,1), dCols(:,2), dCols(:,3), dCols(:,4), dCols(:,5), 'VariableNames', ...
    {'Condition', 'Sensors', 'Detector', 'FAR_pct', 'MinDetectable_CapFraction', 'AUC_at_90pct_cap', ...
     'Det_at_95pct_cap', 'Det_at_90pct_cap', 'Det_at_80pct_cap', 'Det_at_70pct_cap', 'Det_at_60pct_cap'});
writetable(T_sum, fullfile(cfg.resultsDir, 'nuisance_summary.csv'));

rowsCond = {}; rowsSens = {}; cg = []; cw = []; cp = [];
for c = 1:2
    for k = 1:nSens
        rowsCond{end+1, 1} = conds(c).name;   %#ok<SAGROW>
        rowsSens{end+1, 1} = sensors(k).name; %#ok<SAGROW>
        cg(end+1, 1) = cosG(k, c); cw(end+1, 1) = cosW(k, c); cp(end+1, 1) = cosPC1(k, c); %#ok<SAGROW>
    end
end
T_cos = table(rowsCond, rowsSens, cg, cw, cp, 'VariableNames', ...
    {'Condition', 'Sensors', 'cosG_regression_signature', 'cosW_matched_filter', 'cosPC1'});
writetable(T_cos, fullfile(cfg.resultsDir, 'nuisance_physics_cosine.csv'));

%% ---- Figure: rows = condition, columns = detector ----
f1 = figure; set(f1, 'Position', [100 100 950 650]);
tl = tiledlayout(2, nD, 'TileSpacing', 'compact');
mk = {'-o', '-s', '-^', '-d', '-p'};
for c = 1:2
    for d = 1:nD
        nexttile; hold on;
        for k = 1:nSens
            plot(fault_range, detRate(:, k, d, c)*100, mk{k}, 'LineWidth', 1.1, ...
                 'MarkerSize', 3.5, 'DisplayName', sensors(k).name);
        end
        yline(detTarget*100, '--', '95% target', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
        set(gca, 'XDir', 'reverse'); grid on; ylim([0 100]);
        xlabel('Capacitance fraction remaining'); ylabel('Detection rate (%)');
        title(sprintf('%s | %s', conds(c).name, detNames{d}), 'FontSize', 9);
        if c == 1 && d == 1, legend('Location', 'southwest', 'FontSize', 7); end
    end
end
title(tl, 'Nuisance-aware features: detection rate on the multi-point sensor sets');
save_fig(f1, 'fig_nuisance_detection_rate', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'nuisance_aware_features.mat'), ...
    'detRate', 'FAR', 'minIdx', 'aucAll', 'rocSev', 'fault_range', 'detNames', 'sensors', ...
    'cosG', 'cosW', 'cosPC1', 'mComp', 'nPC', 'T_sum', 'T_cos', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'nuisance_aware_features.mat'));

%% ================= local functions =================
function yc = sim_clean(cfg, R, L, C)
    [A, B, ~, ~] = build_rlc_ladder(R, L, C, cfg.RL, cfg.Rs);
    nS = size(A, 1);
    yc = lsim(ss(A, B, eye(nS), zeros(nS, 1)), cfg.u, cfg.t);
end

function Z = gen_states(cfg, sev, nTrials, useTol, refRMS)
% Fresh noisy responses of all states: nTrials x nS x nT (single), each channel
% with its own noise at cfg.snr_db, divided by its nominal healthy RMS.
    N = cfg.N; nS = 2*N; nT = numel(cfg.t);
    Z = zeros(nTrials, nS, nT, 'single');
    if ~useTol
        C = cfg.C_nom;  C(cfg.fault_stage) = C(cfg.fault_stage) * sev;
        yc = sim_clean(cfg, cfg.R_nom, cfg.L_nom, C);
    end
    for i = 1:nTrials
        if useTol
            R = cfg.R_nom .* (1 + cfg.R_tol * (2*rand(1,N) - 1));
            L = cfg.L_nom .* (1 + cfg.L_tol * (2*rand(1,N) - 1));
            C = cfg.C_nom .* (1 + cfg.C_tol * (2*rand(1,N) - 1));
            C(cfg.fault_stage) = C(cfg.fault_stage) * sev;
            yc = sim_clean(cfg, R, L, C);
        end
        noisePower = mean(yc.^2, 1) / (10^(cfg.snr_db/10));
        yn = (yc + sqrt(noisePower) .* randn(size(yc))) ./ refRMS;
        Z(i, :, :) = single(yn.');
    end
end

function B = fit_basis(Z, idx, m)
% PCA basis from the swept data, computed one channel at a time (low memory):
% mean per channel, Gram matrix accumulated over channels, then loadings per channel.
    n = size(Z, 1);  nT = size(Z, 3);  nCh = numel(idx);
    B.mu = zeros(nCh, nT);
    G = zeros(n);
    for q = 1:nCh
        X = double(reshape(Z(:, idx(q), :), n, nT));
        B.mu(q, :) = mean(X, 1);
        X = X - B.mu(q, :);
        G = G + X * X.';
    end
    [U, D] = eig((G + G.') / 2);
    [d, ord] = sort(diag(D), 'descend');
    U = U(:, ord(1:m));  d = d(1:m);
    B.W = cell(nCh, 1);
    for q = 1:nCh
        X = double(reshape(Z(:, idx(q), :), n, nT)) - B.mu(q, :);
        B.W{q} = (X.' * U) ./ sqrt(d).';                  % nT x m
    end
end

function P = project(Z, idx, B)
% Scores of new trials on the basis: nTrials x m.
    n = size(Z, 1);  nT = size(Z, 3);
    P = zeros(n, size(B.W{1}, 2));
    for q = 1:numel(idx)
        X = double(reshape(Z(:, idx(q), :), n, nT)) - B.mu(q, :);
        P = P + X * B.W{q};
    end
end

function sc = score_all(P, m, thrT)
% Threshold-normalised scores for D1 and D2 (alarm when > 1).
    nPC = numel(m.mu1);
    d1 = P(:, 1:nPC) - m.mu1;
    s1 = sqrt(sum((d1 / m.Sig1) .* d1, 2)) / thrT;
    s2 = ((P - m.muH) * m.w) / m.thr2;
    sc = [s1, s2];
end

function v = vec_stack(M)
% Stack a (channels x time) matrix into one column, channel by channel.
    v = M(:);
end

function r = cosine(a, b)
    r = (a(:).' * b(:)) / (norm(a) * norm(b));
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
