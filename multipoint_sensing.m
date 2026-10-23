%% multipoint_sensing.m  (next-step fix 2: measure at more than one place)
%
% PURPOSE : Everything so far looks at ONE signal, the voltage across the last
%           capacitor. The SNR x tolerance map showed that, with one signal,
%           component tolerance puts a floor under the smallest detectable
%           fault, because a 10% capacitance loss moves the dominant mode by
%           about as much as normal part-to-part spread does.
%
%           Different sensors see the fault and the tolerance in different
%           proportions. This script asks: if the detector is given MORE
%           signals from the same circuit, does the floor come down?
%
%   Sensor sets (nested, all from the SAME simulated circuit and noise draw):
%     1  V5            last capacitor voltage        <- the project so far
%     2  V3            capacitor voltage at the FAULTED stage only
%     3  V1,V3,V5      three capacitor voltages
%     4  All caps      all five capacitor voltages
%     5  All states    all five capacitor voltages + all five inductor currents
%
%   Every sensor has its own noise at the project's SNR (15 dB relative to that
%   channel's own signal power), so more sensors = more independent noise
%   sources AND more information. Each channel is divided by its nominal healthy
%   RMS so voltage and current channels are on the same footing before PCA.
%
%   Detector (same idea as the project's swept-PCA T2):
%     PCA basis fit on the swept healthy+faulty training data (60 severities x
%     cfg.trials_per_step trials), Mahalanobis distance in the first nPC scores
%     using the healthy training trials, threshold sqrt(chi2inv(confLevel, nPC)).
%   Two capacities are reported: nPC = 2 (the project's choice) and nPC = 10,
%   because with several channels two components may simply be too few.
%
% PROTOCOL: identical to baseline_comparison.m. Separate healthy training set,
%           every reported rate measured on FRESH trials.
%
% INPUTS  : build_rlc_ladder.m, get_params.m, save_fig.m
% OUTPUTS : results/multipoint_sensing.mat
%           results/multipoint_summary.csv
%           results/multipoint_table_aleatory.csv, results/multipoint_table_combined.csv
%           results/figures/fig_multipoint_detection_rate.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N   = cfg.N;
nS  = 2*N;                        % states: [iL1; vC1; ...; iLN; vCN]
nT  = numel(cfg.t);

nPCList      = [2 10];            % PCA capacity levels compared
nTrainHealthy = 1000;             % healthy trials for mean / covariance of scores
nTestHealthy = 1000;              % fresh healthy trials -> false alarm rate
nTestFault   = 100;               % fresh trials per severity -> detection rate
rocSev       = [0.95 0.90 0.80 0.60];
detTarget    = 0.95;
farMax       = 0.05;
tableSev     = [0.98 0.95 0.90 0.85 0.80 0.70 0.60];

fault_range = cfg.fault_range;
numSteps    = numel(fault_range);

sensors(1).name = 'V5 (last cap, baseline)';   sensors(1).idx = 2*N;
sensors(2).name = 'V3 (faulted stage)';        sensors(2).idx = 2*cfg.fault_stage;
sensors(3).name = 'V1,V3,V5';                  sensors(3).idx = 2*unique([1 cfg.fault_stage N]);
sensors(4).name = 'All 5 caps';                sensors(4).idx = 2*(1:N);
sensors(5).name = 'All caps + currents';       sensors(5).idx = 1:nS;
nSens = numel(sensors);
nPCn  = numel(nPCList);

conds(1).name = 'Aleatory only';        conds(1).useTol = false;
conds(2).name = 'Aleatory + Epistemic'; conds(2).useTol = true;

% ---- per-channel reference RMS: the nominal healthy waveform ----
yNom   = sim_clean(cfg, cfg.R_nom, cfg.L_nom, cfg.C_nom);      % nT x nS
refRMS = sqrt(mean(yNom.^2, 1));                               % 1 x nS

detRate = zeros(numSteps, nSens, nPCn, 2);
FAR     = zeros(nSens, nPCn, 2);
minIdx  = nan(nSens, nPCn, 2);
aucAll  = zeros(numel(rocSev), nSens, nPCn, 2);
dispPC  = zeros(nSens, 2, 2);          % explained-variance share of PC1, PC2 (information only)

for c = 1:2
    useTol = conds(c).useTol;
    fprintf('\n===== %s =====\n', conds(c).name);

    % ---- training sweep (healthy + faulty) and healthy training set ----
    tps    = cfg.trials_per_step;
    Zsweep = zeros(numSteps*tps, nS, nT, 'single');
    for s = 1:numSteps
        Zsweep((s-1)*tps+1 : s*tps, :, :) = gen_states(cfg, fault_range(s), tps, useTol, refRMS);
    end
    Zh = gen_states(cfg, 1.00, nTrainHealthy, useTol, refRMS);

    % ---- fit one model per sensor set ----
    models = cell(nSens, 1);
    for k = 1:nSens
        Xs = feat(Zsweep, sensors(k).idx);
        Xh = feat(Zh,     sensors(k).idx);
        models{k} = fit_model(Xs, Xh, max(nPCList), nPCList, cfg.confLevel);
        dispPC(k, :, c) = models{k}.varShare(1:2).';
        clear Xs Xh
    end
    clear Zsweep Zh

    % ---- false alarm rate on FRESH healthy trials ----
    Zt = gen_states(cfg, 1.00, nTestHealthy, useTol, refRMS);
    scH = cell(nSens, 1);
    for k = 1:nSens
        scH{k} = score_model(feat(Zt, sensors(k).idx), models{k});
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
            sc{k} = score_model(feat(Zf, sensors(k).idx), models{k});
            detRate(s, k, :, c) = reshape(mean(sc{k} > 1, 1), 1, 1, nPCn);
        end
        kr = find(rocIdx == s);
        if ~isempty(kr), scF{kr(1)} = sc; end
        clear Zf
    end

    % ---- minimum detectable fault + AUC ----
    labels = [zeros(nTestHealthy, 1); ones(nTestFault, 1)];
    for k = 1:nSens
        for p = 1:nPCn
            minIdx(k, p, c) = min_detect(detRate(:, k, p, c), FAR(k, p, c), detTarget, farMax);
            for r = 1:numel(rocSev)
                [~, ~, ~, aucAll(r, k, p, c)] = perfcurve(labels, ...
                    [scH{k}(:, p); scF{r}{k}(:, p)], 1);
            end
        end
    end

    % ---- console summary ----
    kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
    [~, sRef] = min(abs(fault_range - 0.90));
    fprintf('%-26s %4s %7s %26s %9s %12s\n', 'Sensors', 'nPC', 'FAR %', ...
            'Min detectable fault', 'AUC@90%', 'Det@90% (%)');
    for k = 1:nSens
        for p = 1:nPCn
            if isnan(minIdx(k, p, c))
                mdTxt = 'not reached';
            else
                mdTxt = sprintf('%.1f%% cap (%.1f%% loss)', fault_range(minIdx(k, p, c))*100, ...
                                (1 - fault_range(minIdx(k, p, c)))*100);
            end
            fprintf('%-26s %4d %7.1f %26s %9.3f %12.0f\n', sensors(k).name, nPCList(p), ...
                    FAR(k, p, c)*100, mdTxt, aucAll(kRef, k, p, c), detRate(sRef, k, p, c)*100);
        end
    end
end

%% ---- Tables ----
rowsCond = {}; rowsSens = {}; colNPC = []; farCol = []; minCol = []; aucCol = []; d90Col = [];
kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
[~, sRef] = min(abs(fault_range - 0.90));
for c = 1:2
    for k = 1:nSens
        for p = 1:nPCn
            rowsCond{end+1, 1} = conds(c).name;     %#ok<SAGROW>
            rowsSens{end+1, 1} = sensors(k).name;   %#ok<SAGROW>
            colNPC(end+1, 1)   = nPCList(p);        %#ok<SAGROW>
            farCol(end+1, 1)   = FAR(k, p, c)*100;  %#ok<SAGROW>
            if isnan(minIdx(k, p, c)), minCol(end+1, 1) = NaN; else, minCol(end+1, 1) = fault_range(minIdx(k, p, c)); end %#ok<SAGROW>
            aucCol(end+1, 1)   = aucAll(kRef, k, p, c);        %#ok<SAGROW>
            d90Col(end+1, 1)   = detRate(sRef, k, p, c)*100;   %#ok<SAGROW>
        end
    end
end
T_sum = table(rowsCond, rowsSens, colNPC, farCol, minCol, aucCol, d90Col, 'VariableNames', ...
    {'Condition', 'Sensors', 'nPC', 'FAR_pct', 'MinDetectable_CapFraction', 'AUC_at_90pct_cap', 'DetRateAt90pct_cap'});
writetable(T_sum, fullfile(cfg.resultsDir, 'multipoint_summary.csv'));

csvNames = {'multipoint_table_aleatory.csv', 'multipoint_table_combined.csv'};
T_rate = cell(1, 2);
for c = 1:2
    varNames = {'CapFraction'};
    M = zeros(numel(tableSev), 1 + nSens*nPCn);
    for k = 1:nSens
        for p = 1:nPCn
            varNames{end+1} = sprintf('S%d_PC%d', k, nPCList(p)); %#ok<SAGROW>
        end
    end
    for r = 1:numel(tableSev)
        [~, idx] = min(abs(fault_range - tableSev(r)));
        M(r, 1) = fault_range(idx);
        col = 2;
        for k = 1:nSens
            for p = 1:nPCn
                M(r, col) = detRate(idx, k, p, c)*100; col = col + 1;
            end
        end
    end
    T_rate{c} = array2table(M, 'VariableNames', varNames);
    writetable(T_rate{c}, fullfile(cfg.resultsDir, csvNames{c}));
    fprintf('\nDetection rate (%%) by severity - %s  (S1..S5 = sensor sets 1..5):\n', conds(c).name);
    disp(T_rate{c});
end
fprintf('\nSensor sets:\n');
for k = 1:nSens, fprintf('  S%d = %s\n', k, sensors(k).name); end

%% ---- Figure: rows = PCA capacity, columns = uncertainty condition ----
f1 = figure; set(f1, 'Position', [100 100 1100 720]);
tl = tiledlayout(nPCn, 2, 'TileSpacing', 'compact');
mk = {'-o', '-s', '-^', '-d', '-p'};
for p = 1:nPCn
    for c = 1:2
        nexttile; hold on;
        for k = 1:nSens
            plot(fault_range, detRate(:, k, p, c)*100, mk{k}, 'LineWidth', 1.2, ...
                 'MarkerSize', 4, 'DisplayName', sensors(k).name);
        end
        yline(detTarget*100, '--', '95% target', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
        set(gca, 'XDir', 'reverse');
        xlabel('Capacitance fraction remaining (1 = healthy)');
        ylabel('Detection rate (%)');
        title(sprintf('%s, nPC = %d', conds(c).name, nPCList(p)));
        grid on; ylim([0 100]);
        if p == 1 && c == 1, legend('Location', 'southwest', 'FontSize', 7); end
    end
end
title(tl, 'Multi-point sensing: detection rate vs. number of measured signals');
save_fig(f1, 'fig_multipoint_detection_rate', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'multipoint_sensing.mat'), ...
    'detRate', 'FAR', 'minIdx', 'aucAll', 'rocSev', 'fault_range', 'nPCList', ...
    'sensors', 'dispPC', 'nTrainHealthy', 'nTestHealthy', 'nTestFault', ...
    'T_sum', 'T_rate', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'multipoint_sensing.mat'));

%% ================= local functions =================
function yc = sim_clean(cfg, R, L, C)
% Noise-free unit-step response of ALL states: nT x nS matrix.
    [A, B, ~, ~] = build_rlc_ladder(R, L, C, cfg.RL, cfg.Rs);
    nS = size(A, 1);
    yc = lsim(ss(A, B, eye(nS), zeros(nS, 1)), cfg.u, cfg.t);
end

function Z = gen_states(cfg, sev, nTrials, useTol, refRMS)
% Fresh noisy responses of all states at one severity: nTrials x nS x nT (single).
% Each channel gets its own noise at cfg.snr_db relative to that channel's own
% signal power (as simulate_trials does), then is divided by its nominal RMS.
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
        noisePower = mean(yc.^2, 1) / (10^(cfg.snr_db/10));      % 1 x nS
        yn = (yc + sqrt(noisePower) .* randn(size(yc))) ./ refRMS;
        Z(i, :, :) = single(yn.');
    end
end

function X = feat(Z, idx)
% Stack the chosen channels of every trial into one feature row: trials x (nT*numel(idx)).
    X = reshape(permute(double(Z(:, idx, :)), [1 3 2]), size(Z, 1), []);
end

function m = fit_model(Xs, Xh, kMax, nPCList, confLevel)
% PCA basis from the swept data (Gram-matrix trick: trials << features), then
% healthy-score mean / covariance for each capacity level in nPCList.
    m.mu = mean(Xs, 1);
    Xc   = Xs - m.mu;
    G    = Xc * Xc.';
    [U, D] = eig((G + G.') / 2);
    [d, ord] = sort(diag(D), 'descend');
    m.varShare = d(1:min(kMax, numel(d))) / sum(max(d, 0));
    U = U(:, ord(1:kMax));  d = d(1:kMax);
    m.W = (Xc.' * U) ./ sqrt(d).';                  % features x kMax, orthonormal columns
    Ph = (Xh - m.mu) * m.W;
    m.nPCList = nPCList;
    m.muH  = cell(1, numel(nPCList));
    m.SigH = cell(1, numel(nPCList));
    m.thr  = zeros(1, numel(nPCList));
    for p = 1:numel(nPCList)
        n = nPCList(p);
        m.muH{p}  = mean(Ph(:, 1:n), 1);
        m.SigH{p} = cov(Ph(:, 1:n));
        m.thr(p)  = sqrt(chi2inv(confLevel, n));
    end
end

function sc = score_model(X, m)
% Threshold-normalised T2 distance per capacity level (alarm when > 1).
    P  = (X - m.mu) * m.W;
    sc = zeros(size(X, 1), numel(m.nPCList));
    for p = 1:numel(m.nPCList)
        n  = m.nPCList(p);
        dP = P(:, 1:n) - m.muH{p};
        sc(:, p) = sqrt(sum((dP / m.SigH{p}) .* dP, 2)) / m.thr(p);
    end
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
