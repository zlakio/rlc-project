function Y = simulate_trials(cfg, sev, nTrials, useTol)
% SIMULATE_TRIALS  Fresh noisy step responses at one fault severity.
%
%   Y = simulate_trials(cfg, sev, nTrials, useTol)
%
%   INPUTS
%     cfg     : struct from get_params()
%     sev     : capacitance fraction remaining at the fault stage (1 = healthy)
%     nTrials : number of noisy trials to generate
%     useTol  : false = aleatory noise only; true = also random R, L, C within
%               tolerance for every trial (epistemic uncertainty)
%   OUTPUT
%     Y       : nTrials x numel(cfg.t) matrix, one trial per row
    N = cfg.N; t = cfg.t; u = cfg.u;
    Y = zeros(nTrials, numel(t));
    if ~useTol
        C = cfg.C_nom;  C(cfg.fault_stage) = C(cfg.fault_stage) * sev;
        [A, B, Cm, D] = build_rlc_ladder(cfg.R_nom, cfg.L_nom, C, cfg.RL, cfg.Rs);
        yc = lsim(ss(A, B, Cm, D), u, t);
    end
    for i = 1:nTrials
        if useTol
            R = cfg.R_nom .* (1 + cfg.R_tol * (2*rand(1,N) - 1));
            L = cfg.L_nom .* (1 + cfg.L_tol * (2*rand(1,N) - 1));
            C = cfg.C_nom .* (1 + cfg.C_tol * (2*rand(1,N) - 1));
            C(cfg.fault_stage) = C(cfg.fault_stage) * sev;
            [A, B, Cm, D] = build_rlc_ladder(R, L, C, cfg.RL, cfg.Rs);
            yc = lsim(ss(A, B, Cm, D), u, t);
        end
        noisePower = mean(yc.^2) / (10^(cfg.snr_db/10));
        Y(i,:) = (yc + sqrt(noisePower) * randn(size(yc))).';
    end
end
