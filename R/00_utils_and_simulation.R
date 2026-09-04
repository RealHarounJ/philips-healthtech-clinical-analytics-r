# 00_utils_and_simulation.R
# Patient cohort simulation and wearable telemetry stream generation
# Based on MIMIC-IV vital sign distributions and NHS NEWS2 standards

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

#' Generate Synthetic Longitudinal ICU / Ward Clinical Cohort
#'
#' Simulates a realistic inpatient cohort mimicking MIMIC-IV / Philips IntelliVue telemetry,
#' incorporating baseline demographics, multi-parameter vital signs, laboratory markers,
#' and clinical deterioration endpoints (acute septic shock / unexpected ICU escalation).
#'
#' @param n_patients Integer, number of unique patients (default: 2500)
#' @param seed Integer, pseudo-random generator seed for strict reproducibility
#' @return A tibble with clinical variables and deterioration status
generate_clinical_cohort <- function(n_patients = 2500, seed = 42) {
  set.seed(seed)
  
  # Patient demographics
  patient_id <- sprintf("PT-%05d", 1:n_patients)
  age <- round(pmax(18, pmin(95, rnorm(n_patients, mean = 66, sd = 14))))
  sex <- sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.54, 0.46))
  charlson_index <- pmax(0, rpois(n_patients, lambda = 2.4))
  
  # Baseline latent risk index (underlying physiological frailty)
  latent_frailty <- 0.03 * (age - 60) + 0.35 * charlson_index + rnorm(n_patients, 0, 0.8)
  
  # Vital Signs (Continuous telemetry spot values)
  # Heart Rate (bpm): Elevated with frailty/stress
  heart_rate <- round(pmax(40, pmin(180, rnorm(n_patients, mean = 78 + 4 * latent_frailty, sd = 16))))
  
  # Mean Arterial Pressure (mmHg): Decreases in shock
  map <- round(pmax(45, pmin(140, rnorm(n_patients, mean = 85 - 5 * latent_frailty, sd = 14))))
  
  # Respiration Rate (breaths/min): Tachypnea is prime indicator of deterioration
  resp_rate <- round(pmax(8, pmin(45, rnorm(n_patients, mean = 17 + 2.5 * latent_frailty, sd = 5))))
  
  # Oxygen Saturation SpO2 (%): Decreases in respiratory compromise
  spo2 <- round(pmax(75, pmin(100, 98 - rbeta(n_patients, shape1 = 1.2, shape2 = 12) * 20 - pmax(0, latent_frailty * 1.5))))
  
  # Body Temperature (°C): Hypothermia (<36) or fever (>38) indicate systemic inflammatory response
  temperature <- round(rnorm(n_patients, mean = 37.1 + 0.2 * latent_frailty, sd = 0.8), 1)
  
  # Glasgow Coma Scale (GCS): 3 (deep coma) to 15 (fully conscious)
  gcs_prob_drop <- plogis(latent_frailty - 1.5)
  gcs <- ifelse(runif(n_patients) < gcs_prob_drop, 
                sample(9:14, n_patients, replace = TRUE), 15)
  
  # Supplemental Oxygen requirement (Binary)
  on_o2 <- ifelse(spo2 < 94 | (runif(n_patients) < plogis(latent_frailty - 1)), 1, 0)
  
  # Laboratory biomarkers (ICU admission panel)
  serum_lactate <- round(pmax(0.5, rgamma(n_patients, shape = 2.5, rate = 1.4) + pmax(0, latent_frailty * 0.8)), 2)
  creatinine <- round(pmax(0.4, rnorm(n_patients, mean = 1.0 + 0.15 * charlson_index, sd = 0.4)), 2)
  wbc_count <- round(pmax(2.0, rnorm(n_patients, mean = 8.5 + 1.2 * latent_frailty, sd = 3.5)), 1)
  
  # Calculate National Early Warning Score 2 (NEWS2) - Clinical Standard
  compute_news2 <- function(hr, bp_map, rr, sat, temp, gcs_val, o2_supp) {
    score <- 0
    # Respiration Rate
    score <- score + ifelse(rr <= 8, 3, ifelse(rr <= 11, 1, ifelse(rr <= 20, 0, ifelse(rr <= 24, 2, 3))))
    # SpO2 (Scale 1)
    score <- score + ifelse(sat <= 91, 3, ifelse(sat <= 93, 2, ifelse(sat <= 95, 1, 0)))
    # Supplemental O2
    score <- score + ifelse(o2_supp == 1, 2, 0)
    # Systolic / MAP proxy
    score <- score + ifelse(bp_map <= 65, 3, ifelse(bp_map <= 75, 2, ifelse(bp_map <= 85, 1, 0)))
    # Heart rate
    score <- score + ifelse(hr <= 40, 3, ifelse(hr <= 50, 1, ifelse(hr <= 90, 0, ifelse(hr <= 110, 1, ifelse(hr <= 130, 2, 3)))))
    # Consciousness (GCS < 15)
    score <- score + ifelse(gcs_val < 15, 3, 0)
    # Temperature
    score <- score + ifelse(temp <= 35.0, 3, ifelse(temp <= 36.0, 1, ifelse(temp <= 38.0, 0, ifelse(temp <= 39.0, 1, 2))))
    return(score)
  }
  
  news2_score <- mapply(compute_news2, heart_rate, map, resp_rate, spo2, temperature, gcs, on_o2)
  
  # True outcome: Clinical Deterioration within 24 hours (Adverse event: ICU transfer / Cardiac arrest)
  # Simulated via logistic sigmoid with realistic clinical coefficients
  log_odds <- -4.0 + 
    0.35 * news2_score + 
    0.65 * (serum_lactate - 2.0) + 
    0.02 * (age - 65) + 
    0.25 * charlson_index +
    0.03 * (heart_rate - 80) - 
    0.04 * (map - 80)
  
  prob_deterioration <- plogis(log_odds)
  deterioration_event <- rbinom(n_patients, size = 1, prob = prob_deterioration)
  
  # Time-to-event for Survival Analysis (hours to event or censoring up to 72 hours)
  time_to_event <- ifelse(
    deterioration_event == 1,
    pmax(1, round(rexp(n_patients, rate = 0.05 + 0.08 * prob_deterioration))),
    72
  )
  time_to_event <- pmin(72, time_to_event)
  
  tibble(
    patient_id = patient_id,
    age = age,
    sex = sex,
    charlson_index = charlson_index,
    heart_rate = heart_rate,
    map = map,
    resp_rate = resp_rate,
    spo2 = spo2,
    temperature = temperature,
    gcs = gcs,
    on_o2 = on_o2,
    serum_lactate = serum_lactate,
    creatinine = creatinine,
    wbc_count = wbc_count,
    news2_score = news2_score,
    prob_true = round(prob_deterioration, 4),
    deterioration_event = deterioration_event,
    time_to_event_hrs = time_to_event
  )
}


#' Generate Synthetic Raw High-Frequency PPG / ECG Telemetry Biosignals
#'
#' Generates raw continuous photoplethysmogram (PPG) waveform at 100 Hz sampling frequency,
#' containing realistic baseline respiratory modulation, dicrotic notches, high-frequency
#' sensor noise, and movement artifacts typical of wearable hospital monitoring devices.
#'
#' @param duration_sec Numeric, length of telemetry window in seconds (default: 60)
#' @param sampling_rate Integer, sampling rate in Hz (default: 100)
#' @param base_hr Numeric, underlying resting heart rate in bpm (default: 75)
#' @param noise_level Numeric, Gaussian sensor noise variance (default: 0.08)
#' @param seed Integer, pseudo-random generator seed
#' @return A tibble with time (sec) and raw amplitude (mV)
generate_raw_ppg_stream <- function(duration_sec = 60, sampling_rate = 100, 
                                    base_hr = 75, noise_level = 0.08, seed = 42) {
  set.seed(seed)
  n_samples <- duration_sec * sampling_rate
  t <- seq(0, duration_sec - 1/sampling_rate, by = 1/sampling_rate)
  
  # Fundamental cardiac pulse frequency
  f_cardiac <- base_hr / 60 # Hz
  
  # Heart Rate Variability (HRV) modulation: Low Frequency (0.1 Hz) + High Frequency (0.25 Hz)
  hrv_modulation <- 0.04 * sin(2 * pi * 0.25 * t) + 0.03 * sin(2 * pi * 0.10 * t)
  instantaneous_phase <- 2 * pi * cumsum((f_cardiac + hrv_modulation) / sampling_rate)
  
  # PPG Pulse morphology (Systolic peak + secondary Diastolic wave/notch)
  systolic_component <- (sin(instantaneous_phase) + 1) / 2
  diastolic_component <- 0.35 * (sin(instantaneous_phase * 2 - 0.8) + 1) / 2
  clean_ppg <- systolic_component^3 + diastolic_component^4
  
  # Baseline wander (Respiration induced thoracic impedance variation, ~0.25 Hz)
  respiratory_baseline <- 0.25 * sin(2 * pi * 0.22 * t + 0.5)
  
  # High frequency electronic noise
  electronic_noise <- rnorm(n_samples, mean = 0, sd = noise_level)
  
  # Random transient motion artifact (e.g., patient repositioning at t=35s)
  motion_artifact <- exp(-((t - 35)^2) / (2 * 1.2^2)) * 0.7 * sin(2 * pi * 3.5 * t)
  
  # Composite raw biosignal
  raw_ppg <- clean_ppg + respiratory_baseline + electronic_noise + motion_artifact
  
  tibble(
    time_sec = t,
    clean_signal = clean_ppg,
    raw_signal = raw_ppg
  )
}
