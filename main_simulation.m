%% main_simulation.m
% Open-loop RLC N-stage ladder: fault sweep + AWGN + PCA + eigenvalue migration
%
% PURPOSE : Sweep one capacitor from healthy -> degraded, simulate noisy step
%           responses, track eigenvalues, run PCA on the noisy outputs.
% INPUTS  : get_params.m (all parameters), build_rlc_ladder.m
% OUTPUTS : results/rlc_fault_simulation_results.mat
%           results/figures/fig_main_healthy_step.png
%           results/figures/fig_main_pca_trajectory.png
%           results/figures/fig_main_eigen_migration.png

clear; clc; close all;
cfg = get_params();
rng(cfg.seed, 'twister');

N = cfg.N; R_nom = cfg.R_nom; L_nom = cfg.L_nom; C_nom = cfg.C_nom;
RL = cfg.RL; Rs = cfg.Rs;
fault_stage = cfg.fault_stage; fault_range = cfg.fault_range;
snr_db = cfg.snr_db; trials_per_step = cfg.trials_per_step;
t = cfg.t; u = cfg.u;

%% ---- Sanity check: healthy system ----
[A0, B0, C0, D0] = build_rlc_ladder(R_nom, L_nom, C_nom, RL, Rs);
sys0 = ss(A0, B0, C0, D0);

fprintf('Healthy system eigenvalues:\n');
disp(eig(A0));
if any(real(eig(A0)) >= 0)
    warning('Nominal system is NOT stable - check component values before proceeding.');
else
    fprintf('Nominal system confirmed stable (all eigenvalues in left half-plane).\n');
end

f1 = figure;
step(sys0);
title('Healthy system step response (sanity check)');
grid on;
save_fig(f1, 'fig_main_healthy_step', cfg);

%% ---- Storage ----
numSteps      = length(fault_range);
eigTrack      = zeros(numSteps, 2*N);
dataMatrix    = zeros(numSteps * trials_per_step, length(t));
severityLabel = zeros(numSteps * trials_per_step, 1);

row = 1;
for s = 1:numSteps
    C_this = C_nom;
    C_this(fault_stage) = C_nom(fault_stage) * fault_range(s);

    [A, B, Cm, D] = build_rlc_ladder(R_nom, L_nom, C_this, RL, Rs);
    sys = ss(A, B, Cm, D);

    e = eig(A);
    if s == 1
        [~, ix] = sortrows([real(e) imag(e)]);   % fixed starting order
        e = e(ix);
    else
        prev = eigTrack(s-1, :).';
        matched = zeros(size(e));
        free = true(size(e));
        for k = 1:numel(prev)                    % match each mode to its nearest
            d = abs(e - prev(k));                % neighbour from the previous step
            d(~free) = Inf;
            [~, j] = min(d);
            matched(k) = e(j);
            free(j) = false;
        end
        e = matched;
    end
    eigTrack(s, :) = e.';
    y_clean = lsim(sys, u, t);

    sigPower   = mean(y_clean.^2);
    noisePower = sigPower / (10^(snr_db/10));

    for tr = 1:trials_per_step
        noise   = sqrt(noisePower) * randn(size(y_clean));
        dataMatrix(row, :) = (y_clean + noise).';
        severityLabel(row) = fault_range(s);
        row = row + 1;
    end
end

%% ---- PCA ----
[coeff, score, ~, ~, explained] = pca(dataMatrix);

f2 = figure;
scatter(score(:,1), score(:,2), 20, severityLabel, 'filled');
colormap(jet);
cb = colorbar;
cb.Label.String = 'Capacitance fraction remaining (1 = healthy)';
xlabel(sprintf('PC1 (%.1f%% variance explained)', explained(1)));
ylabel(sprintf('PC2 (%.1f%% variance explained)', explained(2)));
title('Latent Signature Trajectory');
grid on;
save_fig(f2, 'fig_main_pca_trajectory', cfg);

%% ---- Eigenvalue Migration Map ----
f3 = figure;
hold on;
for k = 1:2*N
    scatter(real(eigTrack(:,k)), imag(eigTrack(:,k)), 15, fault_range, 'filled');
end
xline(0, '--k', 'Stability boundary');
colormap(jet);
cb2 = colorbar;
cb2.Label.String = 'Capacitance fraction remaining (1 = healthy)';
xlabel('Real(\lambda) (rad/s)');
ylabel('Imag(\lambda) (rad/s)');
title('Eigenvalue Migration Map');
grid on; axis equal;
hold off;
save_fig(f3, 'fig_main_eigen_migration', cfg);

%% ---- Save (parameters saved alongside results) ----
save(fullfile(cfg.resultsDir, 'rlc_fault_simulation_results.mat'), ...
     'dataMatrix', 'severityLabel', 'eigTrack', 'fault_range', 't', ...
     'coeff', 'score', 'explained', 'R_nom', 'L_nom', 'C_nom', 'RL', 'Rs', ...
     'fault_stage', 'cfg');

fprintf('\nDone. Results saved to %s\n', fullfile(cfg.resultsDir, 'rlc_fault_simulation_results.mat'));
