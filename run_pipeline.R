# run_pipeline.R
# Master execution script to run all analysis modules and generate output figures

cat(">> Running clinical analytics pipeline...\n\n")

# Source all modular scripts
source("R/00_utils_and_simulation.R")
source("R/01_icu_early_warning_model.R")
source("R/02_wearable_biosignal_analytics.R")
source("R/03_health_economics_outcomes_research.R")

# 1. Generate & Save Inpatient Telemetry Cohort
cat(">> Step 1: Synthesizing Multi-Parameter Telemetry Cohort (N = 2,500)...\n")
cohort_data <- generate_clinical_cohort(n_patients = 2500, seed = 2026)
if (!dir.exists("data")) dir.create("data")
write.csv(cohort_data, "data/clinical_telemetry_cohort.csv", row.names = FALSE)
cat("   Saved: data/clinical_telemetry_cohort.csv\n\n")

# 2. Run Module 1: ICU Early Warning & Decision Curve Analysis
cat(">> Step 2: Executing Module 1 (ICU Early Warning & Clinical Utility)...\n")
mod1_res <- run_early_warning_pipeline(cohort_data, output_dir = "figures")
cat("\n")

# 3. Generate Raw Telemetry & Run Module 2: Wearable DSP & HRV
cat(">> Step 3: Synthesizing High-Frequency PPG Telemetry Waveform...\n")
raw_stream <- generate_raw_ppg_stream(duration_sec = 60, sampling_rate = 100, base_hr = 76, noise_level = 0.08)
write.csv(raw_stream, "data/wearable_raw_ppg_stream.csv", row.names = FALSE)
cat("   Saved: data/wearable_raw_ppg_stream.csv\n")
cat("   Executing Module 2 (Digital Signal Processing & HRV Analytics)...\n")
mod2_res <- run_biosignal_pipeline(raw_stream, output_dir = "figures")
cat("\n")

# 4. Run Module 3: Health Economics & Outcomes Research (HEOR)
cat(">> Step 4: Executing Module 3 (Markov Cost-Effectiveness & PSA)...\n")
mod3_res <- run_psa_simulation(n_sim = 1000, seed = 2026, output_dir = "figures")
cat("\n")

cat("=======================================================================\n")
cat("  PIPELINE EXECUTION COMPLETED SUCCESSFULLY!\n")
cat("  Generated Figures in 'figures/':\n")
cat("   - 01_roc_curve.png\n")
cat("   - 01_calibration_curve.png\n")
cat("   - 01_decision_curve_analysis.png\n")
cat("   - 01_survival_kaplan_meier.png\n")
cat("   - 02_ppg_filtering_and_peak_detection.png\n")
cat("   - 02_hrv_power_spectral_density.png\n")
cat("   - 03_cost_effectiveness_plane.png\n")
cat("   - 03_ceac_acceptability_curve.png\n")
cat("=======================================================================\n")
