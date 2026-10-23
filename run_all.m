%% run_all.m - runs the full pipeline in README order
addpath(pwd);
test_build_rlc_ladder       % 0. verify the model first
main_simulation             % 1. fault sweep, PCA, eigenvalue migration
analyze_eigen_pca_link      % 2. correlation + first detection threshold
epistemic_uncertainty_sweep % 3. Monte Carlo tolerance study
detection_metrics           % 4. FAR, detection rate, ROC/AUC, minimum detectable fault
detection_metrics_q         % 5. adds the Q statistic (T2, Q, T2 OR Q)
baseline_comparison         % 6. PCA vs. healthy-only PCA vs. energy thresholds
eigen_sensitivity           % 7. analytical eigenvalue sensitivity vs. simulation
multipoint_sensing          % 8. detection with 1-10 measured signals
nuisance_aware_features     % 9. supervised matched filter vs swept PCA
kalman_capacitance          % 10. EKF capacitance estimator (slow: roughly 10-20 min)
selfbaseline_detection      % 11. each unit against its own healthy reference
snr_tolerance_heatmap       % 12. SNR x tolerance map (slow: roughly 10-20 min)
