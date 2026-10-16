%% eigen_sensitivity.m  (README roadmap: analytical eigenvalue sensitivity)
%
% PURPOSE : Replace "two trends correlate" with a derivation. For a simple
%           eigenvalue lambda_i of A with right/left eigenvectors v_i, w_i:
%
%               d(lambda_i)/dp = (w_i' * (dA/dp) * v_i) / (w_i' * v_i)
%
%           Capacitor C_k enters A ONLY through row 2k (the dvC_k/dt equation),
%           and scales that whole row by 1/C_k. So for C_k -> s*C_k,
%
%               A(s) = A0 + (1/s - 1) * E_k ,     dA/ds = -E_k / s^2 ,
%
%           where E_k is A0 with every row except 2k zeroed. The script checks
%           this identity numerically, then:
%
%     Part 1  Predicts the whole eigenvalue migration two ways and compares both
%             with the simulated migration from main_simulation.m
%               first-order : lambda(s) = lambda(1) + (s - 1) * dlambda/ds|s=1
%               integrated  : dlambda/ds re-evaluated along the path (Heun/RK2)
%     Part 2  Correlates the ANALYTIC drift of the dominant mode with mean PC1.
%     Part 3  Which stage's capacitor does the dominant mode feel most?
%     Part 4  Same sensitivities for every R, L, C give the eigenvalue spread
%             caused by manufacturing tolerance. Compared with the fault shift,
%             this explains analytically why detection collapses under
%             epistemic uncertainty. A Monte Carlo run checks the linearisation.
%
% INPUTS  : results/rlc_fault_simulation_results.mat   (eigTrack, fault_range)
%           results/rlc_analysis_extended.mat          (meanPC1)
% OUTPUTS : results/eigen_sensitivity.mat, results/eigen_sensitivity_modes.csv
%           results/figures/fig_sens_migration.png
%           results/figures/fig_sens_stage.png
%           results/figures/fig_sens_tolerance.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N = cfg.N;  k = cfg.fault_stage;
refSev = 0.90;     % severity used for the fault-vs-tolerance comparison
nMC    = 2000;     % Monte Carlo draws to check the tolerance linearisation

Sim = load(fullfile(cfg.resultsDir, 'rlc_fault_simulation_results.mat'), 'eigTrack', 'fault_range');
Ext = load(fullfile(cfg.resultsDir, 'rlc_analysis_extended.mat'),        'meanPC1');
sev      = Sim.fault_range(:);       % numSteps x 1, 1.00 -> 0.60
eigSim   = Sim.eigTrack;             % numSteps x 2N, tracked mode by mode
numSteps = numel(sev);
n        = 2 * N;

[A0, ~, ~, ~] = build_rlc_ladder(cfg.R_nom, cfg.L_nom, cfg.C_nom, cfg.RL, cfg.Rs);
Ek = row_pattern(A0, 2*k);

%% ---- Check the structural identity A(s) = A0 + (1/s - 1) E_k ----
sTest = 0.60;
Cts = cfg.C_nom;  Cts(k) = Cts(k) * sTest;
[At, ~, ~, ~] = build_rlc_ladder(cfg.R_nom, cfg.L_nom, Cts, cfg.RL, cfg.Rs);
idErr = norm(At - (A0 + (1/sTest - 1) * Ek), 'fro') / norm(At, 'fro');
fprintf('Identity check  A(0.6) vs A0 + (1/s - 1)E_k : relative error %.2e\n', idErr);
if idErr > 1e-12
    error('Structural identity failed - build_rlc_ladder may have changed.');
end

%% ---- Nominal eigenstructure, ordered like the simulated tracks ----
[V0, D0] = eig(A0);
lam0 = diag(D0);
W0   = inv(V0);                      % rows = left eigenvectors, W0 * V0 = I

gaps = abs(lam0 - lam0.');  gaps(1:n+1:end) = Inf;
minGap = min(gaps(:));
fprintf('Smallest gap between nominal eigenvalues: %.3g rad/s (modes must be simple)\n', minGap);
if minGap < 1e-6 * max(abs(lam0))
    warning('Repeated eigenvalues: first-order sensitivity formula is not valid.');
end

perm = zeros(1, n);  free = true(n, 1);        % column j of eigSim <-> mode perm(j)
for j = 1:n
    d = abs(lam0 - eigSim(1, j));  d(~free) = Inf;
    [dmin, i] = min(d);
    perm(j) = i;  free(i) = false;
    if dmin > 1e-6 * max(abs(lam0))
        warning('Saved eigTrack(1,%d) does not match eig(A0); re-run main_simulation.m.', j);
    end
end

[~, dom] = max(abs(eigSim(end, :) - eigSim(1, :)));   % same rule as analyze_eigen_pca_link.m
fprintf('Dominant mode (largest drift): column %d, lambda(healthy) = %.2f %+.2fi\n\n', ...
    dom, real(eigSim(1, dom)), imag(eigSim(1, dom)));

%% ================= PART 1: predicted vs simulated migration =================
% ---- first-order, linearised once at the healthy point ----
dl1 = zeros(n, 1);
for i = 1:n
    dl1(i) = W0(i, :) * (-Ek) * V0(:, i);            % dlambda_i/ds at s = 1
end
lam1 = zeros(numSteps, n);
for j = 1:n
    lam1(:, j) = lam0(perm(j)) + (sev - 1) * dl1(perm(j));
end

% ---- integrated: re-linearise along the path (Heun / RK2) ----
lamI = zeros(numSteps, n);
for j = 1:n
    lamI(1, j) = lam0(perm(j));
    for s = 1:numSteps-1
        h    = sev(s+1) - sev(s);                    % negative: capacitance falls
        g1   = sens_at(A0, Ek, sev(s),   lamI(s, j));
        pred = lamI(s, j) + h * g1;
        g2   = sens_at(A0, Ek, sev(s+1), pred);
        lamI(s+1, j) = lamI(s, j) + 0.5 * h * (g1 + g2);
    end
end

err1 = abs(lam1 - eigSim);
errI = abs(lamI - eigSim);

fprintf('Prediction error vs simulated migration (rad/s)\n');
fprintf('%-8s %14s %16s %16s %18s %18s\n', 'Cap.', '|dlam| (dom.)', 'err 1st (dom.)', ...
    'err integ (dom.)', 'max err 1st (all)', 'max err integ (all)');
for sv = [0.95 0.90 0.80 0.60]
    [~, ix] = min(abs(sev - sv));
    fprintf('%-8.2f %14.3f %16.3f %16.3g %18.3f %18.3g\n', sev(ix), ...
        abs(eigSim(ix, dom) - eigSim(1, dom)), err1(ix, dom), errI(ix, dom), ...
        max(err1(ix, :)), max(errI(ix, :)));
end

f1 = figure;
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
plot(real(eigSim(:, dom)), imag(eigSim(:, dom)), 'ko', 'DisplayName', 'Simulated (eig)');
plot(real(lam1(:, dom)),   imag(lam1(:, dom)),   '-',  'LineWidth', 1.5, 'DisplayName', 'First-order');
plot(real(lamI(:, dom)),   imag(lamI(:, dom)),   '--', 'LineWidth', 1.5, 'DisplayName', 'Integrated');
xlabel('Real(\lambda) (rad/s)'); ylabel('Imag(\lambda) (rad/s)');
title('Dominant mode: analytic vs simulated'); legend('Location', 'best'); grid on;
nexttile;
semilogy(sev, max(err1(:, dom), 1e-12), '-o', 'LineWidth', 1.3, 'DisplayName', 'First-order'); hold on;
semilogy(sev, max(errI(:, dom), 1e-12), '-s', 'LineWidth', 1.3, 'DisplayName', 'Integrated');
set(gca, 'XDir', 'reverse');
xlabel('Capacitance fraction remaining (1 = healthy)'); ylabel('|prediction error| (rad/s)');
title('Dominant mode: prediction error'); legend('Location', 'best'); grid on;
title(tl, 'Analytical eigenvalue sensitivity vs. simulated migration');
save_fig(f1, 'fig_sens_migration', cfg);

%% ================= PART 2: analytic drift vs PCA latent trajectory =================
driftSim = abs(eigSim(:, dom) - eigSim(1, dom));
drift1   = abs(lam1(:, dom)   - lam1(1, dom));
driftI   = abs(lamI(:, dom)   - lamI(1, dom));
pc1      = Ext.meanPC1(:);
rSim = corr_xy(driftSim, pc1);
r1   = corr_xy(drift1,   pc1);
rI   = corr_xy(driftI,   pc1);
fprintf('\nCorrelation of dominant-mode drift with mean PC1:\n');
fprintf('  simulated eigenvalues : r = %.4f\n', rSim);
fprintf('  analytic, integrated  : r = %.4f\n', rI);
fprintf('  analytic, first-order : r = %.4f  (linear in s, so this mainly measures linearity of PC1 in s)\n', r1);

%% ================= PART 3: which stage does the dominant mode feel? =================
idom  = perm(dom);
sensC = zeros(N, 1);                  % dlambda_dom / d ln(C_kk)   (complex)
for kk = 1:N
    sensC(kk) = W0(idom, :) * (-row_pattern(A0, 2*kk)) * V0(:, idom);
end
per1pct = 0.01 * abs(sensC);          % |shift| in rad/s for a 1% capacitance change
fprintf('\nDominant-mode shift per 1%% capacitance change, by stage:\n');
for kk = 1:N
    mark = ''; if kk == k, mark = '   <- fault stage used in the simulations'; end
    fprintf('  stage %d : %.3f rad/s%s\n', kk, per1pct(kk), mark);
end

f2 = figure;
b = bar(per1pct); b.FaceColor = 'flat';
b.CData(:) = repmat([0.6 0.6 0.6], N, 1);  b.CData(k, :) = [0.85 0.33 0.10];
xlabel('Stage whose capacitor changes'); ylabel('|\Delta\lambda_{dominant}| per 1% change (rad/s)');
title('Which capacitor does the dominant mode feel? (orange = simulated fault)'); grid on;
save_fig(f2, 'fig_sens_stage', cfg);

%% ================= PART 4: tolerance spread vs fault shift =================
% Parameters in the same order as simulate_trials: R_1..R_N, L_1..L_N, C_1..C_N
pNom   = [cfg.R_nom, cfg.L_nom, cfg.C_nom];
tolRel = [cfg.R_tol * ones(1, N), cfg.L_tol * ones(1, N), cfg.C_tol * ones(1, N)];
nP     = 3 * N;
hRel   = 1e-6;
dLam   = zeros(n, nP);                % dlambda_i / d ln(p_j)
for ip = 1:nP
    pp = pNom;  pm = pNom;
    pp(ip) = pNom(ip) * (1 + hRel);
    pm(ip) = pNom(ip) * (1 - hRel);
    [Ap, ~, ~, ~] = build_rlc_ladder(pp(1:N), pp(N+1:2*N), pp(2*N+1:3*N), cfg.RL, cfg.Rs);
    [Am, ~, ~, ~] = build_rlc_ladder(pm(1:N), pm(N+1:2*N), pm(2*N+1:3*N), cfg.RL, cfg.Rs);
    dA = (Ap - Am) / (2 * hRel);
    for i = 1:n
        dLam(i, ip) = W0(i, :) * dA * V0(:, i);
    end
end

sig   = tolRel / sqrt(3);             % std of a uniform +/- tol relative error
varRe = (real(dLam).^2) * (sig.^2).';
varIm = (imag(dLam).^2) * (sig.^2).';
sigLin = sqrt(varRe + varIm);         % linearised tolerance spread of each mode (rad/s)
idxC   = 2*N + k;                     % column of the faulted capacitor itself
shareC = ((real(dLam(:, idxC))*sig(idxC)).^2 + (imag(dLam(:, idxC))*sig(idxC)).^2) ./ sigLin.^2;

% ---- Monte Carlo check of the linearisation ----
reMC = zeros(nMC, n);  imMC = zeros(nMC, n);
for q = 1:nMC
    R = cfg.R_nom .* (1 + cfg.R_tol * (2*rand(1, N) - 1));
    L = cfg.L_nom .* (1 + cfg.L_tol * (2*rand(1, N) - 1));
    C = cfg.C_nom .* (1 + cfg.C_tol * (2*rand(1, N) - 1));
    [Aq, ~, ~, ~] = build_rlc_ladder(R, L, C, cfg.RL, cfg.Rs);
    eq = eig(Aq);  freeq = true(n, 1);
    for i = 1:n                                    % match to nominal modes
        d = abs(eq - lam0(i));  d(~freeq) = Inf;
        [~, jq] = min(d);  freeq(jq) = false;
        reMC(q, i) = real(eq(jq));  imMC(q, i) = imag(eq(jq));
    end
end
sigMC = sqrt(var(reMC, 0, 1) + var(imMC, 0, 1)).';
fprintf('\nMonte Carlo vs linearised tolerance spread (largest relative gap: %.1f%%)\n', ...
    100 * max(abs(sigMC - sigLin) ./ sigMC));
gapRatio = minGap / (3 * max(sigMC));
if gapRatio < 1
    fprintf('NOTE: 3-sigma spread exceeds the smallest nominal gap - nearest-neighbour matching may swap modes.\n');
end

% ---- fault shift at refSev vs tolerance spread, mode by mode (sim ordering) ----
[~, ir] = min(abs(sev - refSev));
shiftRef = abs(lamI(ir, :) - lamI(1, :)).';                % n x 1, sim column order
sigLinS  = sigLin(perm);  sigMCS = sigMC(perm);  shareCS = shareC(perm);
ratioMC  = shiftRef ./ sigMCS;

fprintf('\nFault shift at %.0f%% capacitance vs. tolerance spread, per mode (rad/s)\n', sev(ir)*100);
fprintf('%-5s %20s %12s %12s %12s %10s %14s\n', 'Mode', 'lambda (healthy)', 'fault shift', ...
    'tol sigma lin', 'tol sigma MC', 'shift/sig', 'tol share C_k');
for j = 1:n
    fprintf('%-5d %20s %12.3f %12.3f %12.3f %10.2f %13.0f%%\n', j, ...
        sprintf('%.0f%+.0fi', real(eigSim(1, j)), imag(eigSim(1, j))), ...
        shiftRef(j), sigLinS(j), sigMCS(j), ratioMC(j), 100 * shareCS(j));
end
fprintf('(tol share C_k = fraction of that mode''s tolerance variance caused by the faulted capacitor''s own tolerance;\n');
fprintf(' that part moves the mode along the SAME direction as the fault and cannot be separated from it.)\n');

f3 = figure;
bar([shiftRef, sigMCS]);
set(gca, 'YScale', 'log');
xlabel('Mode (column of eigTrack)'); ylabel('Eigenvalue shift (rad/s)');
legend({sprintf('Fault shift at %.0f%% capacitance', sev(ir)*100), 'Tolerance spread (1 sigma, Monte Carlo)'}, ...
       'Location', 'best');
title('Fault signal vs. manufacturing-tolerance spread, per eigenvalue mode'); grid on;
save_fig(f3, 'fig_sens_tolerance', cfg);

%% ---- Save ----
modeIdx = (1:n).';
T_modes = table(modeIdx, real(eigSim(1, :)).', imag(eigSim(1, :)).', shiftRef, sigLinS, sigMCS, ratioMC, shareCS, ...
    'VariableNames', {'Mode', 'Re_lambda', 'Im_lambda', 'FaultShift_refSev', 'TolSigma_linear', ...
                      'TolSigma_MC', 'Shift_over_Sigma', 'TolShare_FaultedC'});
writetable(T_modes, fullfile(cfg.resultsDir, 'eigen_sensitivity_modes.csv'));
save(fullfile(cfg.resultsDir, 'eigen_sensitivity.mat'), ...
    'sev', 'eigSim', 'lam1', 'lamI', 'err1', 'errI', 'dom', 'rSim', 'rI', 'r1', ...
    'sensC', 'per1pct', 'dLam', 'sigLin', 'sigMC', 'shareC', 'shiftRef', 'refSev', 'T_modes', 'cfg');
fprintf('\nDone. Saved to %s\n', fullfile(cfg.resultsDir, 'eigen_sensitivity.mat'));

%% ================= local functions =================
function E = row_pattern(A, r)
% Copy of A with every row except r set to zero.
    E = zeros(size(A));
    E(r, :) = A(r, :);
end

function dl = sens_at(A0, Ek, s, lamGuess)
% dlambda/ds of the eigenvalue nearest lamGuess, for A(s) = A0 + (1/s - 1) Ek.
    As = A0 + (1/s - 1) * Ek;
    [V, D] = eig(As);
    lam = diag(D);
    [~, i] = min(abs(lam - lamGuess));
    W = inv(V);
    dl = -(W(i, :) * Ek * V(:, i)) / s^2;
end

function r = corr_xy(x, y)
    R = corrcoef(x(:), y(:));
    r = R(1, 2);
end
