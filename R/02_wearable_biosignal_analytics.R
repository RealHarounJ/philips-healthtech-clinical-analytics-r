# 02_wearable_biosignal_analytics.R
# Wearable PPG biosignal processing, baseline detrending, and HRV spectral analysis

suppressPackageStartupMessages({
  library(tidyverse)
})

#' Apply Moving Average Filter (Pure R)
moving_average <- function(x, n = 5) {
  stats::filter(x, rep(1 / n, n), sides = 2) %>% as.numeric()
}

#' Digital Bandpass & Baseline Wander Removal Filter
#'
#' Removes respiratory baseline wander (0.1-0.3 Hz) and high-frequency noise
#' using a zero-phase cascade of high-pass baseline subtraction and low-pass smoothing.
clean_biosignal <- function(raw_signal, sampling_rate = 100) {
  # Baseline drift estimation (moving median / broad moving average of 1.2 seconds)
  w_baseline <- round(sampling_rate * 1.2)
  baseline <- moving_average(raw_signal, n = w_baseline)
  # Handle NA boundaries
  baseline[is.na(baseline)] <- mean(raw_signal, na.rm = TRUE)
  
  # Detrended signal
  detrended <- raw_signal - baseline
  
  # High-frequency smoothing (moving average of 50 ms = 5 samples at 100 Hz)
  w_smooth <- round(sampling_rate * 0.05)
  smoothed <- moving_average(detrended, n = w_smooth)
  smoothed[is.na(smoothed)] <- detrended[is.na(smoothed)]
  
  return(smoothed)
}

#' Adaptive Systolic Peak Detection (R-Peak / Systolic Pulse)
detect_systolic_peaks <- function(signal, time_vec, sampling_rate = 100, min_distance_sec = 0.45) {
  min_samples <- round(min_distance_sec * sampling_rate)
  n <- length(signal)
  peaks_idx <- integer(0)
  
  # Adaptive threshold: upper 60th percentile
  thresh <- stats::quantile(signal, probs = 0.65)
  
  i <- 2
  while (i < (n - 1)) {
    if (signal[i] > thresh && signal[i] > signal[i - 1] && signal[i] >= signal[i + 1]) {
      peaks_idx <- c(peaks_idx, i)
      i <- i + min_samples # refractory period
    } else {
      i <- i + 1
    }
  }
  
  tibble(
    peak_idx = peaks_idx,
    time_sec = time_vec[peaks_idx],
    amplitude = signal[peaks_idx]
  )
}

#' Compute Comprehensive Time-Domain and Frequency-Domain HRV Metrics
compute_hrv_metrics <- function(rr_intervals_ms) {
  # Time-Domain Metrics
  mean_rr <- mean(rr_intervals_ms)
  mean_hr <- 60000 / mean_rr
  sdnn    <- sd(rr_intervals_ms)
  
  diff_rr <- diff(rr_intervals_ms)
  rmssd   <- sqrt(mean(diff_rr^2))
  pnn50   <- (sum(abs(diff_rr) > 50) / length(diff_rr)) * 100
  
  # Frequency-Domain Spectral Analysis (FFT on interpolated RR series)
  # Resample RR time series to 4 Hz grid
  t_rr <- cumsum(rr_intervals_ms) / 1000 # seconds
  fs_resample <- 4 # Hz
  t_grid <- seq(min(t_rr), max(t_rr), by = 1 / fs_resample)
  rr_interp <- stats::approx(t_rr, rr_intervals_ms, xout = t_grid, rule = 2)$y
  
  # Remove mean & apply Hanning window
  rr_zero_mean <- rr_interp - mean(rr_interp)
  N <- length(rr_zero_mean)
  hanning <- 0.5 * (1 - cos(2 * pi * (0:(N - 1)) / (N - 1)))
  windowed_signal <- rr_zero_mean * hanning
  
  # Fast Fourier Transform
  fft_vals <- stats::fft(windowed_signal)
  psd <- (Mod(fft_vals)^2) / (N * fs_resample)
  freqs <- (0:(N - 1)) * (fs_resample / N)
  
  # Take positive spectrum
  half_n <- floor(N / 2)
  freqs <- freqs[1:half_n]
  psd   <- psd[1:half_n] * 2
  
  # Band powers (Trapezoidal integration)
  integrate_band <- function(f, p, f_min, f_max) {
    idx <- which(f >= f_min & f <= f_max)
    if (length(idx) < 2) return(0)
    sum(diff(f[idx]) * (p[idx[-length(idx)]] + p[idx[-1]]) / 2)
  }
  
  vlf_power <- integrate_band(freqs, psd, 0.0033, 0.04)
  lf_power  <- integrate_band(freqs, psd, 0.04, 0.15)
  hf_power  <- integrate_band(freqs, psd, 0.15, 0.40)
  total_power <- vlf_power + lf_power + hf_power
  lf_hf_ratio <- ifelse(hf_power > 0, lf_power / hf_power, NA)
  
  list(
    time_domain = tibble(
      Mean_RR_ms = round(mean_rr, 1),
      Mean_HR_bpm = round(mean_hr, 1),
      SDNN_ms = round(sdnn, 2),
      RMSSD_ms = round(rmssd, 2),
      pNN50_pct = round(pnn50, 2)
    ),
    freq_domain = tibble(
      VLF_ms2 = round(vlf_power, 1),
      LF_ms2  = round(lf_power, 1),
      HF_ms2  = round(hf_power, 1),
      Total_Power_ms2 = round(total_power, 1),
      LF_HF_Ratio = round(lf_hf_ratio, 2)
    ),
    spectrum = tibble(frequency_hz = freqs, psd = psd)
  )
}

#' Run Wearable Biosignal Processing & Telemetry Pipeline
run_biosignal_pipeline <- function(raw_stream, output_dir = "figures") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("--- [Module 2] Processing High-Frequency Wearable PPG Telemetry ---\n")
  
  # 1. Filter raw PPG waveform
  filtered <- clean_biosignal(raw_stream$raw_signal, sampling_rate = 100)
  raw_stream$filtered_signal <- filtered
  
  # 2. Systolic Peak Detection
  peaks <- detect_systolic_peaks(filtered, raw_stream$time_sec, sampling_rate = 100)
  
  # 3. Inter-Beat Intervals (RR / Peak-to-Peak in milliseconds)
  rr_intervals <- diff(peaks$time_sec) * 1000
  
  # 4. HRV Analytics
  hrv_results <- compute_hrv_metrics(rr_intervals)
  
  cat("  Wearable HRV Metrics Computed:\n")
  print(hrv_results$time_domain)
  cat("  Spectral Autonomic Balance:\n")
  print(hrv_results$freq_domain)
  
  # 5. Figure 2A: Filtering and Peak Detection Waveform (15-second window for clarity)
  zoom_stream <- raw_stream %>% filter(time_sec >= 10 & time_sec <= 25)
  zoom_peaks  <- peaks %>% filter(time_sec >= 10 & time_sec <= 25)
  
  p_waveform <- ggplot() +
    geom_line(data = zoom_stream, aes(x = time_sec, y = raw_signal, color = "Raw Telemetry (Noise + Drift)"), alpha = 0.5) +
    geom_line(data = zoom_stream, aes(x = time_sec, y = filtered_signal, color = "Zero-Phase Filtered PPG"), linewidth = 1) +
    geom_point(data = zoom_peaks, aes(x = time_sec, y = amplitude, shape = "Detected Systolic Peak"), color = "red", size = 3) +
    scale_color_manual(values = c("Raw Telemetry (Noise + Drift)" = "gray40", "Zero-Phase Filtered PPG" = "navyblue")) +
    scale_shape_manual(values = c("Detected Systolic Peak" = 17)) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Real-Time Telemetry Digital Signal Processing (DSP)",
      subtitle = "Baseline wander removal, high-frequency filtering and adaptive peak detection",
      x = "Telemetry Time (Seconds)",
      y = "PPG Amplitude (Arbitrary Units / Volts)",
      color = "Signal Stream",
      shape = "Feature"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "02_ppg_filtering_and_peak_detection.png"), 
         plot = p_waveform, width = 10, height = 5.5, dpi = 300)
  
  # 6. Figure 2B: HRV Power Spectral Density (FFT)
  spectrum_df <- hrv_results$spectrum %>% filter(frequency_hz <= 0.45)
  
  p_spectral <- ggplot(spectrum_df, aes(x = frequency_hz, y = psd)) +
    geom_area(data = filter(spectrum_df, frequency_hz >= 0.04 & frequency_hz <= 0.15),
              aes(fill = "LF Band (0.04 - 0.15 Hz): Sympathetic & Baroreflex"), alpha = 0.4) +
    geom_area(data = filter(spectrum_df, frequency_hz >= 0.15 & frequency_hz <= 0.40),
              aes(fill = "HF Band (0.15 - 0.40 Hz): Vagal / Parasympathetic"), alpha = 0.4) +
    geom_line(color = "black", linewidth = 0.9) +
    scale_fill_manual(values = c("orange", "skyblue")) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Autonomic Spectral Analysis (Heart Rate Variability)",
      subtitle = sprintf("LF/HF Ratio: %.2f | Vagal Tone (HF): %.1f ms² | Sympathetic Tone (LF): %.1f ms²",
                         hrv_results$freq_domain$LF_HF_Ratio,
                         hrv_results$freq_domain$HF_ms2,
                         hrv_results$freq_domain$LF_ms2),
      x = "Frequency (Hz)",
      y = "Power Spectral Density (ms² / Hz)",
      fill = "Autonomic Frequency Band"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "02_hrv_power_spectral_density.png"), 
         plot = p_spectral, width = 9, height = 5.5, dpi = 300)
  
  cat("  [Module 2] Saved Figures: 02_ppg_filtering_and_peak_detection.png, 02_hrv_power_spectral_density.png\n")
  
  return(list(
    peaks = peaks,
    rr_intervals = rr_intervals,
    hrv_metrics = hrv_results
  ))
}
