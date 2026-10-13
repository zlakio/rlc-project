function [A, B, Cm, D] = build_rlc_ladder(R, L, C, RL, Rs)
% BUILD_RLC_LADDER  Generate state-space matrices for an N-stage
% series-RL / shunt-C ladder network, driven by a Thevenin source
% with internal resistance Rs, terminated by a load resistor RL.
%
%   Topology (N stages):
%
%   u(t) --Rs-- R1,L1 --+-- R2,L2 --+-- ... --+-- RN,LN --+
%                        |          |         |           |
%                       C1(gnd)    C2(gnd)   ...        CN(gnd) -- RL(gnd)
%
%   State vector: x = [iL1; vC1; iL2; vC2; ...; iLN; vCN]   (size 2N)
%   Input:        u  = source voltage
%   Output:       y  = vCN (voltage across the load)
%
%   INPUTS
%     R  : 1xN vector, series resistance of each stage (Ohms)
%     L  : 1xN vector, series inductance of each stage (Henries)
%     C  : 1xN vector, shunt capacitance of each stage (Farads)
%     RL : scalar, load resistance at the final node (Ohms)
%     Rs : scalar, source (internal) resistance (Ohms)
%
%   OUTPUTS
%     A, B, Cm, D : state-space matrices (D is 0, Cm avoids shadowing
%                   the built-in C)

    N = length(R);
    if length(L) ~= N || length(C) ~= N
        error('R, L, and C must all be vectors of the same length N.');
    end

    n = 2*N;
    A  = zeros(n, n);
    B  = zeros(n, 1);
    Cm = zeros(1, n);
    D  = 0;

    for k = 1:N
        iL_idx = 2*k - 1;
        vC_idx = 2*k;

        % ---- inductor current equation: diLk/dt ----
        if k == 1
            % first stage sees source resistance in series
            A(iL_idx, iL_idx) = -(R(1) + Rs) / L(1);
            B(iL_idx)         = 1 / L(1);
        else
            vCprev_idx = 2*(k-1);
            A(iL_idx, iL_idx)     = -R(k) / L(k);
            A(iL_idx, vCprev_idx) = 1 / L(k);
        end
        A(iL_idx, vC_idx) = -1 / L(k);

        % ---- capacitor voltage equation: dvCk/dt ----
        if k < N
            iLnext_idx = 2*(k+1) - 1;
            A(vC_idx, iL_idx)     = 1 / C(k);
            A(vC_idx, iLnext_idx) = -1 / C(k);
        else
            % last stage also discharges through the load resistor
            A(vC_idx, iL_idx) = 1 / C(k);
            A(vC_idx, vC_idx) = -1 / (RL * C(k));
        end
    end

    Cm(2*N) = 1;  % output = voltage across the final capacitor / load

end