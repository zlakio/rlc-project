%% test_build_rlc_ladder.m  (README 1.3)
%
% Verifies build_rlc_ladder.m against closed-form eigenvalues for N = 1.
%
%   Test A: pure series RLC (Rs = 0, RL = Inf)
%       lambda = -R/(2L) +/- sqrt((R/(2L))^2 - 1/(LC))
%   Test B: same, but with source resistance (R -> R + Rs)
%   Test C: N = 1 with finite load RL. States [i; v]:
%       s^2 + (R/L + 1/(RL*C)) s + (1/(LC))*(1 + R/RL) = 0
%   Test D: N = 5 nominal system is stable.
%
% PASS criterion: max relative eigenvalue error < 1e-10.

clear; clc;
tol = 1e-10;
R = 10; L = 10e-3; C = 100e-6;

% ---- Test A ----
[A,~,~,~] = build_rlc_ladder(R, L, C, Inf, 0);
lam_ref = roots([1, R/L, 1/(L*C)]);
errA = relerr(eig(A), lam_ref);

% ---- Test B ----
Rs = 5;
[A,~,~,~] = build_rlc_ladder(R, L, C, Inf, Rs);
lam_ref = roots([1, (R+Rs)/L, 1/(L*C)]);
errB = relerr(eig(A), lam_ref);

% ---- Test C ----
RL = 50;
[A,~,~,~] = build_rlc_ladder(R, L, C, RL, 0);
lam_ref = roots([1, R/L + 1/(RL*C), (1/(L*C))*(1 + R/RL)]);
errC = relerr(eig(A), lam_ref);

% ---- Test D ----
cfg = get_params();
[A5,~,~,~] = build_rlc_ladder(cfg.R_nom, cfg.L_nom, cfg.C_nom, cfg.RL, cfg.Rs);
stableD = all(real(eig(A5)) < 0);

fprintf('Test A (series RLC)          : max rel. error = %.2e  %s\n', errA, passfail(errA < tol));
fprintf('Test B (with Rs)             : max rel. error = %.2e  %s\n', errB, passfail(errB < tol));
fprintf('Test C (finite RL)           : max rel. error = %.2e  %s\n', errC, passfail(errC < tol));
fprintf('Test D (N=5 nominal stable)  :                         %s\n', passfail(stableD));

if errA < tol && errB < tol && errC < tol && stableD
    fprintf('\nALL TESTS PASSED - build_rlc_ladder.m verified.\n');
else
    error('One or more tests FAILED - do not trust simulation results yet.');
end

function e = relerr(a, b)
    a = sortrows([real(a) imag(a)]);
    b = sortrows([real(b) imag(b)]);
    e = max(abs((a(:,1)+1i*a(:,2)) - (b(:,1)+1i*b(:,2))) ./ abs(b(:,1)+1i*b(:,2)));
end
function s = passfail(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
