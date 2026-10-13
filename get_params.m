function cfg = get_params()
% GET_PARAMS  Single source of truth for every shared parameter (README 1.1).
%
%   Every script calls:   cfg = get_params();
%   Change a value HERE, never inside an individual script.
%
%   OUTPUT
%     cfg : struct with circuit values, fault sweep, noise, tolerance,
%           Monte Carlo count, RNG seed, and output folders.

    % ---- reproducibility ----
    cfg.seed = 42;

    % ---- circuit ----
    cfg.N      = 5;                          % ladder stages
    cfg.R_nom  = 10    * ones(1, cfg.N);     % Ohms
    cfg.L_nom  = 10e-3 * ones(1, cfg.N);     % Henries
    cfg.C_nom  = 100e-6* ones(1, cfg.N);     % Farads
    cfg.RL     = 50;                         % load resistance, Ohms
    cfg.Rs     = 5;                          % source resistance, Ohms

    % ---- fault sweep ----
    cfg.fault_stage = ceil(cfg.N/2);               % which capacitor degrades
    cfg.fault_range = linspace(1.00, 0.60, 60);    % 1.00 = healthy, 0.60 = 40% loss

    % ---- aleatory uncertainty (measurement noise) ----
    cfg.snr_db          = 15;    % dB
    cfg.trials_per_step = 20;    % Monte Carlo repeats per severity level

    % ---- epistemic uncertainty (manufacturing tolerance) ----
    cfg.R_tol = 0.05;            % +/-5%
    cfg.L_tol = 0.10;            % +/-10%
    cfg.C_tol = 0.10;            % +/-10%

    % ---- time axis / input ----
    cfg.t = 0:1e-5:0.02;         % 20 ms window, 10 us resolution
    cfg.u = ones(size(cfg.t));   % unit step

    % ---- detection ----
    cfg.confLevel = 0.99;

    % ---- output folders (created if missing) ----
    cfg.resultsDir = 'results';
    cfg.figDir     = fullfile('results', 'figures');
    if ~exist(cfg.resultsDir, 'dir'), mkdir(cfg.resultsDir); end
    if ~exist(cfg.figDir,     'dir'), mkdir(cfg.figDir);     end
end
