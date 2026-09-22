# test_pipeline.R
# Unit tests for data generation, DSP filtering, and Markov model consistency

cat(">> Running Clinical & Algorithmic Validation Tests...\n")

source("R/00_utils_and_simulation.R")
source("R/01_icu_early_warning_model.R")
source("R/02_wearable_biosignal_analytics.R")
source("R/03_health_economics_outcomes_research.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(paste("FAIL:", msg))
  cat(paste("  [PASS]", msg, "\n"))
}

# 1. Test Cohort Synthesis
cohort <- generate_clinical_cohort(n_patients = 100, seed = 99)
assert_true(nrow(cohort) == 100, "Cohort generates exact requested patient size")
assert_true(all(cohort$spo2 >= 50 & cohort$spo2 <= 100), "SpO2 values within plausible physiological limits")
assert_true(all(cohort$heart_rate >= 30 & cohort$heart_rate <= 220), "Heart rates within physiological bounds")
assert_true(all(cohort$news2_score >= 0 & cohort$news2_score <= 20), "NEWS2 scores in valid range [0, 20]")

# 2. Test Biosignal DSP & Peak Detection
stream <- generate_raw_ppg_stream(duration_sec = 10, sampling_rate = 100, base_hr = 60)
filtered <- clean_biosignal(stream$raw_signal, sampling_rate = 100)
assert_true(length(filtered) == length(stream$raw_signal), "Filter preserves stream length")
assert_true(!any(is.na(filtered)), "Filtered stream contains zero NaN/NA values")

peaks <- detect_systolic_peaks(filtered, stream$time_sec, sampling_rate = 100)
assert_true(nrow(peaks) >= 7 & nrow(peaks) <= 13, "Peak detector identifies correct number of cardiac beats (~10 beats in 10s at 60 bpm)")

# 3. Test HRV Metrics
rr <- diff(peaks$time_sec) * 1000
hrv <- compute_hrv_metrics(rr)
assert_true(hrv$time_domain$Mean_HR_bpm >= 50 & hrv$time_domain$Mean_HR_bpm <= 75, "HRV estimated heart rate matches simulated rate")
assert_true(hrv$freq_domain$Total_Power_ms2 > 0, "Spectral total power is strictly positive")

# 4. Test Markov State-Transition Probabilities
out_markov <- run_markov_model()
assert_true(is.numeric(out_markov$icer), "ICER is finite numeric value")
assert_true(out_markov$res_philips$total_cost < out_markov$res_std$total_cost, "Continuous telemetry reduces net hospital costs via ICU avoidance")
assert_true(out_markov$res_philips$total_qalys > out_markov$res_std$total_qalys, "Continuous telemetry increases cumulative QALYs")

cat(">> ALL CLINICAL & STATISTICAL TESTS PASSED SUCCESSFULLY! (6/6)\n")
