%% kalman_capacitance.m  (next-step fix 3: estimate the faulty capacitance directly)
%
% PURPOSE : Instead of asking "does this response look unusual?" (PCA), estimate
%           the quantity that actually changes: the capacitance of the stage that
%           degrades. An extended Kalman filter (EKF) tracks the 10 circuit
%           states AND one extra hidden state, theta = C_stage3 / C_nominal.
%           This is the "physics-informed digital twin, corrected against the
%           measurement" scheme of the review paper (Fig. 5), and it returns an
%           estimate WITH its own uncertainty.
%
%   Model used by the filter (nominal circuit; only theta is allowed to differ):
%       x_{k+1} = x_k + dt * ( A(theta) x_k + B u ),   A(theta) = A0 + (1/theta - 1) E,
%   where E has the single non-zero row of the faulted capacitor (the same
%   identity checked in eigen_sensitivity.m). Each measured channel is processed
%   as a scalar update (channels have independent noise).
%
%   Two versions of the filter, to separate noise from tolerance:
%       qrel = 0      trusts the nominal model completely (no model-mismatch allowance)
%       qrel = 1e-3   adds process noise of 0.1% of each state's RMS per step, so
%                     the filter can absorb mismatch from unknown R, L, C
%                     (epistemic uncertainty) at the price of a slower, noisier estimate
%
%   Decision rule: after the 20 ms window, alarm if the final estimate
%   theta_hat < tau, where tau is the 1st percentile of theta_hat over a separate
%   healthy calibration set (so the nominal false alarm rate is 1%).
%
%   Sensor sets are the same five as multipoint_sensing.m. Channels are measured
%   in physical units with noise at cfg.snr_db relative to each channel's own power.
%
% PROTOCOL: separate healthy calibration set, every reported rate on FRESH trials
%           (500 healthy, 100 per severity on a 2% grid from 0% to 40% loss).
%           Needs pagemtimes (MATLAB R2020b or newer). Runtime: roughly 10-20 min.
%
% INPUTS  : build_rlc_ladder.m, get_params.m, save_fig.m
% OUTPUTS : results/kalman_capacitance.mat
%           results/kalman_summary.csv, results/kalman_estimates.csv
%           results/figures/fig_kalman_estimates.png
%           results/figures/fig_kalman_detection_rate.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N   = cfg.N;
nS  = 2*N;
nT  = numel(cfg.t);
dt  = cfg.t(2) - cfg.t(1);
k3  = 2*cfg.fault_stage;          % state index of the faulted capacitor voltage

qrelList     = [0 1e-3];
qTheta       = 0;                 % random-walk noise on theta per step (0: theta is constant)
p0Theta      = 0.1^2;             % initial variance of theta (nominal 1)
nCal         = 500;
nTestHealthy = 500;
nTestFault   = 100;
sevGrid      = 1 - (0:0.02:0.40); % 21 severities, 1 = healthy
rocSev       = [0.94 0.90 0.80 0.60];
tableSev     = [1.00 0.90 0.80 0.70 0.60];
detTarget    = 0.95;
farMax       = 0.05;
pctlLow      = 100 * (1 - cfg.confLevel);

sensors(1).name = 'V5 (last cap, baseline)';   sensors(1).idx = 2*N;
sensors(2).name = 'V3 (faulted stage)';        sensors(2).idx = 2*cfg.fault_stage;
sensors(3).name = 'V1,V3,V5';                  sensors(3).idx = 2*unique([1 cfg.fault_stage N]);
sensors(4).name = 'All 5 caps';                sensors(4).idx = 2*(1:N);
sensors(5).name = 'All caps + currents';       sensors(5).idx = 1:nS;
nSens = numel(sensors);
nQ    = numel(qrelList);
nSev  = numel(sevGrid);

conds(1).name = 'Aleatory only';        conds(1).useTol = false;
conds(2).name = 'Aleatory + Epistemic'; conds(2).useTol = true;

% ---- nominal model and per-channel scales ----
[A0, B0, ~, ~] = build_rlc_ladder(cfg.R_nom, cfg.L_nom, cfg.C_nom, cfg.RL, cfg.Rs);
yNom = sim_clean(cfg, cfg.R_nom, cfg.L_nom, cfg.C_nom);          % nT x nS
refRMS = sqrt(mean(yNom.^2, 1));                                 % 1 x nS
M.dt = dt;  M.A0 = A0;  M.B0 = B0;  M.k3 = k3;
M.xscale = refRMS;                                               % typical state size
M.Rv = refRMS.^2 / 10^(cfg.snr_db/10);                           % nominal noise variance per channel
M.qTheta = qTheta;  M.p0Theta = p0Theta;

detRate = zeros(nSev, nSens, nQ, 2);
FAR     = zeros(nSens, nQ, 2);
minIdx  = nan(nSens, nQ, 2);
aucAll  = zeros(numel(rocSev), nSens, nQ, 2);
estMean = zeros(nSev, nSens, nQ, 2);
estStd  = zeros(nSev, nSens, nQ, 2);
estSig  = zeros(nSev, nSens, nQ, 2);     % the filter's own claimed 1-sigma

for c = 1:2
    useTol = conds(c).useTol;
    fprintf('\n===== %s =====\n', conds(c).name);

    % ---- data: calibration, fresh healthy, fresh faulty (all severities in one batch) ----
    Ycal = gen_phys(cfg, 1.00, nCal, useTol, refRMS);
    Yh   = gen_phys(cfg, 1.00, nTestHealthy, useTol, refRMS);
    Yf   = zeros(nSev*nTestFault, nS, nT, 'single');
    for s = 1:nSev
        Yf((s-1)*nTestFault+1 : s*nTestFault, :, :) = gen_phys(cfg, sevGrid(s), nTestFault, useTol, refRMS);
    end

    for k = 1:nSens
        for q = 1:nQ
            t0 = tic;
            [thCal, ~]  = run_ekf(Ycal, sensors(k).idx, qrelList(q), M);
            [thH,   ~]  = run_ekf(Yh,   sensors(k).idx, qrelList(q), M);
            [thF, sgF]  = run_ekf(Yf,   sensors(k).idx, qrelList(q), M);

            tau = prctile(thCal, pctlLow);                 % alarm when theta_hat < tau
            scal = @(th) (1 - th) / (1 - tau);             % >1  <=>  theta_hat < tau
            scH  = scal(thH);
            FAR(k, q, c) = mean(scH > 1);
            for s = 1:nSev
                rows = (s-1)*nTestFault+1 : s*nTestFault;
                detRate(s, k, q, c) = mean(scal(thF(rows)) > 1);
                estMean(s, k, q, c) = mean(thF(rows));
                estStd(s, k, q, c)  = std(thF(rows));
                estSig(s, k, q, c)  = mean(sgF(rows));
            end
            minIdx(k, q, c) = min_detect(detRate(:, k, q, c), FAR(k, q, c), detTarget, farMax);
            labels = [zeros(nTestHealthy, 1); ones(nTestFault, 1)];
            for r = 1:numel(rocSev)
                [~, si] = min(abs(sevGrid - rocSev(r)));
                rows = (si-1)*nTestFault+1 : si*nTestFault;
                [~, ~, ~, aucAll(r, k, q, c)] = perfcurve(labels, [scH; scal(thF(rows))], 1);
            end
            fprintf('  %-26s qrel=%-6g done in %.0f s\n', sensors(k).name, qrelList(q), toc(t0));
        end
    end
    clear Ycal Yh Yf

    % ---- console summary ----
    kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
    [~, sRef] = min(abs(sevGrid - 0.90));
    fprintf('\n%-26s %7s %6s %24s %9s %9s %12s\n', 'Sensors', 'qrel', 'FAR %', ...
            'Min detectable fault', 'AUC@90%', 'Det@90%', 'std(theta)@90');
    for k = 1:nSens
        for q = 1:nQ
            if isnan(minIdx(k, q, c))
                mdTxt = 'not reached';
            else
                mdTxt = sprintf('%.1f%% cap (%.1f%% loss)', sevGrid(minIdx(k, q, c))*100, ...
                                (1 - sevGrid(minIdx(k, q, c)))*100);
            end
            fprintf('%-26s %7g %6.1f %24s %9.3f %9.0f %12.3f\n', sensors(k).name, qrelList(q), ...
                    FAR(k, q, c)*100, mdTxt, aucAll(kRef, k, q, c), detRate(sRef, k, q, c)*100, ...
                    estStd(sRef, k, q, c));
        end
    end
end

%% ---- Tables ----
rowsCond = {}; rowsSens = {}; qCol = []; farCol = []; minCol = []; aucCol = []; dCols = zeros(0, numel(tableSev));
kRef = find(abs(rocSev - 0.90) < 1e-9, 1);
for c = 1:2
    for k = 1:nSens
        for q = 1:nQ
            rowsCond{end+1, 1} = conds(c).name;     %#ok<SAGROW>
            rowsSens{end+1, 1} = sensors(k).name;   %#ok<SAGROW>
            qCol(end+1, 1)   = qrelList(q);         %#ok<SAGROW>
            farCol(end+1, 1) = FAR(k, q, c)*100;    %#ok<SAGROW>
            if isnan(minIdx(k, q, c)), minCol(end+1, 1) = NaN; else, minCol(end+1, 1) = sevGrid(minIdx(k, q, c)); end %#ok<SAGROW>
            aucCol(end+1, 1) = aucAll(kRef, k, q, c); %#ok<SAGROW>
            row = zeros(1, numel(tableSev));
            for r = 1:numel(tableSev)
                [~, ii] = min(abs(sevGrid - tableSev(r)));
                row(r) = detRate(ii, k, q, c)*100;
            end
            dCols(end+1, :) = row; %#ok<SAGROW>
        end
    end
end
T_sum = table(rowsCond, rowsSens, qCol, farCol, minCol, aucCol, ...
    dCols(:,1), dCols(:,2), dCols(:,3), dCols(:,4), dCols(:,5), 'VariableNames', ...
    {'Condition', 'Sensors', 'qrel', 'FAR_pct', 'MinDetectable_CapFraction', 'AUC_at_90pct_cap', ...
     'Det_at_100pct_cap', 'Det_at_90pct_cap', 'Det_at_80pct_cap', 'Det_at_70pct_cap', 'Det_at_60pct_cap'});
writetable(T_sum, fullfile(cfg.resultsDir, 'kalman_summary.csv'));

rowsCond = {}; rowsSens = {}; qCol = []; tr = []; em = []; es = []; eg = [];
for c = 1:2
    for k = 1:nSens
        for q = 1:nQ
            for r = 1:numel(tableSev)
                [~, ii] = min(abs(sevGrid - tableSev(r)));
                rowsCond{end+1, 1} = conds(c).name;    %#ok<SAGROW>
                rowsSens{end+1, 1} = sensors(k).name;  %#ok<SAGROW>
                qCol(end+1, 1) = qrelList(q);          %#ok<SAGROW>
                tr(end+1, 1) = sevGrid(ii);            %#ok<SAGROW>
                em(end+1, 1) = estMean(ii, k, q, c);   %#ok<SAGROW>
                es(end+1, 1) = estStd(ii, k, q, c);    %#ok<SAGROW>
                eg(end+1, 1) = estSig(ii, k, q, c);    %#ok<SAGROW>
            end
        end
    end
end
T_est = table(rowsCond, rowsSens, qCol, tr, em, es, eg, 'VariableNames', ...
    {'Condition', 'Sensors', 'qrel', 'TrueCapFraction', 'MeanEstimate', 'StdEstimate', 'FilterClaimedSigma'});
writetable(T_est, fullfile(cfg.resultsDir, 'kalman_estimates.csv'));

%% ---- Figure 1: estimate vs truth (mean +/- 1 std over trials) ----
showSens = [1 3 5];
f1 = figure; set(f1, 'Position', [100 100 1250 650]);
tl = tiledlayout(2, nQ, 'TileSpacing', 'compact');
col = lines(nSens);
for c = 1:2
    for q = 1:nQ
        nexttile; hold on;
        plot([0.6 1], [0.6 1], 'k--', 'HandleVisibility', 'off');
        for k = showSens
            errorbar(sevGrid, estMean(:, k, q, c), estStd(:, k, q, c), '-o', 'Color', col(k, :), ...
                     'MarkerSize', 3, 'CapSize', 3, 'DisplayName', sensors(k).name);
        end
        set(gca, 'XDir', 'reverse'); grid on;
        xlabel('True capacitance fraction (stage 3)'); ylabel('Estimated fraction (mean \pm 1 std)');
        title(sprintf('%s, qrel = %g', conds(c).name, qrelList(q)), 'FontSize', 9);
        if c == 1 && q == 1, legend('Location', 'northeast', 'FontSize', 7); end
    end
end
title(tl, 'Kalman filter estimate of the faulted capacitance');
save_fig(f1, 'fig_kalman_estimates', cfg);

%% ---- Figure 2: detection rate ----
f2 = figure; set(f2, 'Position', [100 100 1250 650]);
tl = tiledlayout(2, nQ, 'TileSpacing', 'compact');
mk = {'-o', '-s', '-^', '-d', '-p'};
for c = 1:2
    for q = 1:nQ
        nexttile; hold on;
        for k = 1:nSens
            plot(sevGrid, detRate(:, k, q, c)*100, mk{k}, 'LineWidth', 1.1, 'MarkerSize', 3.5, ...
                 'DisplayName', sensors(k).name);
        end
        yline(detTarget*100, '--', '95% target', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
        set(gca, 'XDir', 'reverse'); grid on; ylim([0 100]);
        xlabel('Capacitance fraction remaining'); ylabel('Detection rate (%)');
        title(sprintf('%s, qrel = %g', conds(c).name, qrelList(q)), 'FontSize', 9);
        if c == 1 && q == 1, legend('Location', 'southwest', 'FontSize', 7); end
    end
end
title(tl, 'Kalman-filter detector: detection rate on the multi-point sensor sets');
save_fig(f2, 'fig_kalman_detection_rate', cfg);

%% ---- Save ----
save(fullfile(cfg.resultsDir, 'kalman_capacitance.mat'), ...
    'detRate', 'FAR', 'minIdx', 'aucAll', 'estMean', 'estStd', 'estSig', 'sevGrid', 'rocSev', ...
    'qrelList', 'sensors', 'T_sum', 'T_est', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'kalman_capacitance.mat'));

%% ================= local functions =================
function yc = sim_clean(cfg, R, L, C)
    [A, B, ~, ~] = build_rlc_ladder(R, L, C, cfg.RL, cfg.Rs);
    nS = size(A, 1);
    yc = lsim(ss(A, B, eye(nS), zeros(nS, 1)), cfg.u, cfg.t);
end

function Y = gen_phys(cfg, sev, nTrials, useTol, refRMS)
% Fresh noisy responses of all states in PHYSICAL units: nTrials x nS x nT (single).
% Each channel has its own noise at cfg.snr_db relative to its own signal power.
    N = cfg.N; nS = 2*N; nT = numel(cfg.t);
    Y = zeros(nTrials, nS, nT, 'single');
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
        Y(i, :, :) = single((yc + sqrt(noisePower) .* randn(size(yc))).');
    end
end

function [th, sig] = run_ekf(Y, idx, qrel, M)
% Extended Kalman filter for the 10 states plus theta = C_stage3/C_nominal,
% run on all trials of Y in parallel (page-wise 11x11 algebra).
%   th  : final theta estimate per trial (nTrials x 1)
%   sig : the filter's own final 1-sigma for theta (nTrials x 1)
    n  = size(Y, 1);  nT = size(Y, 3);
    dt = M.dt;  A0 = M.A0;  B0 = M.B0;  k3 = M.k3;
    a3 = A0(k3, :);                                   % 1 x 10, capacitor row
    X  = zeros(11, n);  X(11, :) = 1;                 % states at rest, theta = 1
    P  = zeros(11, 11, n);
    for i = 1:10, P(i, i, :) = 1e-12 * M.xscale(i)^2; end
    P(11, 11, :) = M.p0Theta;
    Q  = zeros(11);  Q(1:10, 1:10) = diag((qrel * M.xscale).^2);  Q(11, 11) = M.qTheta^2;
    F0 = eye(11);    F0(1:10, 1:10) = F0(1:10, 1:10) + dt * A0;
    Ek = zeros(11);  Ek(k3, 1:10) = dt * a3;
    for k = 2:nT
        thk = X(11, :);
        x   = X(1:10, :);
        r   = a3 * x;                                 % 1 x n
        Ax  = A0 * x;
        Ax(k3, :) = Ax(k3, :) + (1 ./ thk - 1) .* r;
        F   = repmat(F0, 1, 1, n) + reshape(1 ./ thk - 1, 1, 1, n) .* Ek;
        F(k3, 11, :) = reshape(-dt ./ thk.^2 .* r, 1, 1, n);
        X(1:10, :) = x + dt * (Ax + B0);              % unit step input
        P = pagemtimes(pagemtimes(F, P), 'none', F, 'transpose') + Q;
        for ch = idx
            yk = double(reshape(Y(:, ch, k), 1, n));
            Pc = reshape(P(:, ch, :), 11, n);         % column ch of P (= row, by symmetry)
            s  = Pc(ch, :) + M.Rv(ch);
            K  = Pc ./ s;
            X  = X + K .* (yk - X(ch, :));
            P  = P - reshape(K, 11, 1, n) .* reshape(Pc, 1, 11, n);
        end
        P = (P + permute(P, [2 1 3])) / 2;
        X(11, :) = min(max(X(11, :), 0.3), 1.6);
    end
    th  = X(11, :).';
    sig = sqrt(reshape(P(11, 11, :), n, 1));
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
