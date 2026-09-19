# 03_health_economics_outcomes_research.R
# Markov cohort state-transition model & Monte Carlo Probabilistic Sensitivity Analysis (PSA)

suppressPackageStartupMessages({
  library(tidyverse)
})

#' Run Deterministic Markov Cohort Model for Inpatient Monitoring Strategies
#'
#' Evaluates Standard Intermittent Care (every 4-6h) vs Philips Connected Care continuous telemetry
#' across a 30-day time horizon for a cohort of 1,000 acute ward inpatients.
run_markov_model <- function(
  n_cohort = 1000,
  horizon_days = 30,
  cost_ward = 650,          # Cost per day in General Ward (€)
  cost_icu = 2800,          # Cost per day in ICU (€)
  cost_telemetry_std = 15,  # Consumables / spot-check monitor amortized (€/day)
  cost_telemetry_philips = 85, # Philips Continuous Telemetry & analytics license (€/day)
  utility_ward = 0.85,      # Daily QALY utility weight
  utility_icu = 0.40,       # ICU QALY utility weight
  utility_post = 0.80,      # Discharged post-recovery utility
  early_rescue_rr = 0.62    # Relative risk of severe ICU deterioration with Philips telemetry (38% reduction)
) {
  # 4 States: 1=Ward, 2=Deteriorating, 3=ICU, 4=Recovered/Discharged
  # State transitions per day:
  # Standard Care Transition Matrix
  P_std <- matrix(c(
    0.78, 0.12, 0.02, 0.08,  # From Ward
    0.05, 0.45, 0.40, 0.10,  # From Deteriorating (unrecognized -> 40% escalate to ICU)
    0.00, 0.00, 0.82, 0.18,  # From ICU (mean stay ~ 5.5 days)
    0.00, 0.00, 0.00, 1.00   # Absorbing / Discharged
  ), nrow = 4, byrow = TRUE)
  
  # Philips Connected Care Transition Matrix (Early detection rescues deterioration early)
  P_philips <- matrix(c(
    0.80, 0.08, 0.01, 0.11,  # From Ward
    0.35, 0.35, 0.40 * early_rescue_rr, 0.15, # Early rescue pushes 35% back to stable ward!
    0.00, 0.00, 0.75, 0.25,  # Faster ICU discharge due to earlier intervention (mean stay 4 days)
    0.00, 0.00, 0.00, 1.00   # Absorbing
  ), nrow = 4, byrow = TRUE)
  
  # Normalize row probabilities
  P_philips <- P_philips / rowSums(P_philips)
  
  # Initial distribution: 100% in General Ward
  init_state <- c(1, 0, 0, 0)
  
  # Simulation function
  simulate_cohort <- function(P, c_monitor, early_stay_red = 1.0) {
    state_trace <- matrix(0, nrow = horizon_days, ncol = 4)
    curr <- init_state
    
    total_cost <- 0
    total_qalys <- 0
    icu_days <- 0
    
    for (t in 1:horizon_days) {
      curr <- curr %*% P
      state_trace[t, ] <- curr
      
      # Daily patient counts
      n_ward <- curr[1] * n_cohort
      n_det  <- curr[2] * n_cohort
      n_icu  <- curr[3] * n_cohort
      n_disc <- curr[4] * n_cohort
      
      icu_days <- icu_days + n_icu
      
      # Costs (€)
      cost_t <- (n_ward + n_det) * (cost_ward + c_monitor) + 
                n_icu * (cost_icu + c_monitor)
      total_cost <- total_cost + cost_t
      
      # QALYs (daily QALY weight / 365)
      qaly_t <- ((n_ward + n_det) * utility_ward + 
                 n_icu * utility_icu + 
                 n_disc * utility_post) / 365
      total_qalys <- total_qalys + qaly_t
    }
    
    list(
      trace = state_trace,
      total_cost = total_cost,
      cost_per_patient = total_cost / n_cohort,
      total_qalys = total_qalys,
      qalys_per_patient = total_qalys / n_cohort,
      total_icu_days = icu_days,
      icu_days_per_patient = icu_days / n_cohort
    )
  }
  
  res_std <- simulate_cohort(P_std, cost_telemetry_std)
  res_philips <- simulate_cohort(P_philips, cost_telemetry_philips)
  
  # Incremental Cost Effectiveness Ratio (ICER)
  delta_cost <- res_philips$cost_per_patient - res_std$cost_per_patient
  delta_qaly <- res_philips$qalys_per_patient - res_std$qalys_per_patient
  icer <- delta_cost / delta_qaly
  
  list(
    res_std = res_std,
    res_philips = res_philips,
    delta_cost = delta_cost,
    delta_qaly = delta_qaly,
    icer = icer
  )
}


#' Run Probabilistic Sensitivity Analysis (PSA) Monte Carlo Simulation
run_psa_simulation <- function(n_sim = 1000, seed = 42, output_dir = "figures") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("--- [Module 3] Running Health Economics & Outcomes Research (HEOR) ---\n")
  cat(sprintf("  Running %d-iteration Monte Carlo Probabilistic Sensitivity Analysis...\n", n_sim))
  set.seed(seed)
  
  # Sample distributions
  rr_samples <- rbeta(n_sim, shape1 = 15, shape2 = 10) * 0.4 + 0.45 # RR around 0.60
  cost_icu_samples <- rgamma(n_sim, shape = 28, scale = 100) # mean ~ 2800
  cost_ward_samples <- rgamma(n_sim, shape = 26, scale = 25) # mean ~ 650
  philips_lic_samples <- rnorm(n_sim, mean = 85, sd = 10)
  
  psa_results <- vector("list", n_sim)
  for (i in 1:n_sim) {
    out <- run_markov_model(
      cost_icu = cost_icu_samples[i],
      cost_ward = cost_ward_samples[i],
      cost_telemetry_philips = philips_lic_samples[i],
      early_rescue_rr = rr_samples[i]
    )
    psa_results[[i]] <- tibble(
      sim = i,
      delta_cost = out$delta_cost,
      delta_qaly = out$delta_qaly,
      icer = out$icer,
      icu_days_saved = out$res_std$icu_days_per_patient - out$res_philips$icu_days_per_patient
    )
  }
  
  psa_df <- bind_rows(psa_results)
  
  mean_dcost <- mean(psa_df$delta_cost)
  mean_dqaly <- mean(psa_df$delta_qaly)
  mean_icu_saved <- mean(psa_df$icu_days_saved)
  prob_cost_saving <- mean(psa_df$delta_cost < 0) * 100
  
  cat(sprintf("  Mean Incremental Cost: €%.1f per patient\n", mean_dcost))
  cat(sprintf("  Mean QALYs Gained: +%.4f per patient\n", mean_dqaly))
  cat(sprintf("  Mean ICU Days Averted: %.2f days/patient\n", mean_icu_saved))
  cat(sprintf("  Probability of Cost-Saving (Dominant): %.1f%%\n", prob_cost_saving))
  
  # Figure 3A: Cost-Effectiveness Plane
  wtp_threshold <- 30000 # €30,000 / QALY (NICE / EU standard)
  
  p_cep <- ggplot(psa_df, aes(x = delta_qaly, y = delta_cost)) +
    geom_point(alpha = 0.35, color = "navyblue", size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_abline(slope = wtp_threshold, intercept = 0, color = "forestgreen", linewidth = 1.1, linetype = "dotdash") +
    geom_point(data = tibble(x = mean_dqaly, y = mean_dcost), aes(x = x, y = y), color = "red", size = 4, shape = 18) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Cost-Effectiveness Plane (1,000 Monte Carlo PSA Iterations)",
      subtitle = sprintf("Red diamond = Base case | Green line = €30,000/QALY Willingness-to-Pay | %.1f%% Cost-Saving", prob_cost_saving),
      x = "Incremental Effectiveness (Δ QALYs Gained)",
      y = "Incremental Cost per Inpatient (Δ Cost €)"
    ) +
    annotate("text", x = max(psa_df$delta_qaly) * 0.8, y = min(psa_df$delta_cost) * 0.9, 
             label = "DOMINANT QUADRANT\n(Higher QALYs, Lower Costs)", color = "forestgreen", fontface = "bold") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "03_cost_effectiveness_plane.png"), 
         plot = p_cep, width = 9, height = 6, dpi = 300)
  
  # Figure 3B: Cost-Effectiveness Acceptability Curve (CEAC)
  wtp_range <- seq(0, 60000, by = 2000)
  ceac_df <- tibble(
    wtp = wtp_range,
    prob_ce = sapply(wtp_range, function(lambda) {
      mean((lambda * psa_df$delta_qaly - psa_df$delta_cost) > 0)
    })
  )
  
  p_ceac <- ggplot(ceac_df, aes(x = wtp, y = prob_ce)) +
    geom_line(color = "navyblue", linewidth = 1.3) +
    geom_vline(xintercept = 30000, linetype = "dashed", color = "forestgreen") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_x_continuous(labels = scales::dollar_format(prefix = "€", big.mark = ",")) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Cost-Effectiveness Acceptability Curve (CEAC)",
      subtitle = "Probability that Philips Connected Care is Cost-Effective as a function of Payer WTP Threshold",
      x = "Willingness-to-Pay Threshold (€ / QALY)",
      y = "Probability of Cost-Effectiveness"
    ) +
    annotate("text", x = 32000, y = 0.5, label = "Reference Standard\n€30.000 / QALY", hjust = 0, color = "forestgreen") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "03_ceac_acceptability_curve.png"), 
         plot = p_ceac, width = 9, height = 5.5, dpi = 300)
  
  cat("  [Module 3] Saved Figures: 03_cost_effectiveness_plane.png, 03_ceac_acceptability_curve.png\n")
  
  return(list(
    psa_df = psa_df,
    summary = tibble(
      Mean_Delta_Cost_EUR = round(mean_dcost, 2),
      Mean_Delta_QALY = round(mean_dqaly, 4),
      Mean_ICU_Days_Saved = round(mean_icu_saved, 2),
      Pct_Cost_Saving_Dominant = round(prob_cost_saving, 1)
    )
  ))
}
