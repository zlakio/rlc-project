%% run_all.m - runs the full pipeline in README order
addpath(pwd);
test_build_rlc_ladder       % 0. verify the model first
main_simulation             % 1. fault sweep, PCA, eigenvalue migration
analyze_eigen_pca_link      % 2. correlation + detection threshold
epistemic_uncertainty_sweep % 3. Monte Carlo tolerance study
