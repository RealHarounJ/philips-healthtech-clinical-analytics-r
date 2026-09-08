# 01_icu_early_warning_model.R
# Predictive deterioration modeling, ROC validation and Decision Curve Analysis (DCA)
# Note: Net benefit computation follows Vickers & Elkin (Med Decis Making, 2006)

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

#' Compute Area Under ROC Curve (Trapezoidal Rule)
calc_roc_auc <- function(actual, predicted) {
  ord <- order(predicted, decreasing = TRUE)
  actual <- actual[ord]
  tp <- cumsum(actual)
  fp <- cumsum(1 - actual)
  tpr <- tp / sum(actual)
  fpr <- fp / sum(1 - actual)
  
  # Trapezoidal integration
  auc <- sum((fpr[-1] - fpr[-length(fpr)]) * (tpr[-1] + tpr[-length(tpr)]) / 2)
  return(list(auc = auc, tpr = tpr, fpr = fpr))
}

#' Run Comprehensive ICU Early Warning Modeling Pipeline
#'
#' Fits baseline clinical NEWS2 model, multivariable predictive telemetry model,
#' Kaplan-Meier survival curves, and evaluates Decision Curve Analysis (Net Benefit).
#'
#' @param cohort_data Tibble from generate_clinical_cohort()
#' @param output_dir Character, directory to save generated plots
#' @return A list of model objects and performance metrics
run_early_warning_pipeline <- function(cohort_data, output_dir = "figures") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("--- [Module 1] Training Predictive Clinical Early Warning Models ---\n")
  
  # 1. Train/Test Split (70/30 stratified)
  set.seed(123)
  n <- nrow(cohort_data)
  train_idx <- sample(1:n, size = floor(0.70 * n))
  train_data <- cohort_data[train_idx, ]
  test_data  <- cohort_data[-train_idx, ]
  
  # 2. Benchmark Model: Standard NEWS2 Logistic Regression
  fit_news2 <- glm(deterioration_event ~ news2_score, 
                   data = train_data, family = binomial)
  
  # 3. Advanced Connected Care Model: Multivariable Physiological Telemetry
  fit_telemetry <- glm(
    deterioration_event ~ news2_score + heart_rate + map + resp_rate + 
      spo2 + serum_lactate + creatinine + age + charlson_index,
    data = train_data, family = binomial
  )
  
  # Predict test probabilities
  test_data$pred_news2 <- predict(fit_news2, newdata = test_data, type = "response")
  test_data$pred_telemetry <- predict(fit_telemetry, newdata = test_data, type = "response")
  
  # 4. Performance Metrics: ROC-AUC and Brier Score
  roc_news2 <- calc_roc_auc(test_data$deterioration_event, test_data$pred_news2)
  roc_telemetry <- calc_roc_auc(test_data$deterioration_event, test_data$pred_telemetry)
  
  brier_news2 <- mean((test_data$pred_news2 - test_data$deterioration_event)^2)
  brier_telemetry <- mean((test_data$pred_telemetry - test_data$deterioration_event)^2)
  
  cat(sprintf("  Standard NEWS2 Benchmark   -> AUC: %.3f | Brier: %.4f\n", 
              roc_news2$auc, brier_news2))
  cat(sprintf("  Connected Care Full Model  -> AUC: %.3f | Brier: %.4f\n", 
              roc_telemetry$auc, brier_telemetry))
  
  # 5. Figure 1: Comparative ROC Curve & Calibration Plot
  # ROC Plot Data
  roc_df <- bind_rows(
    tibble(fpr = roc_news2$fpr, tpr = roc_news2$tpr, Model = sprintf("NEWS2 Benchmark (AUC = %.2f)", roc_news2$auc)),
    tibble(fpr = roc_telemetry$fpr, tpr = roc_telemetry$tpr, Model = sprintf("Connected Care Model (AUC = %.2f)", roc_telemetry$auc))
  )
  
  p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, color = Model)) +
    geom_line(linewidth = 1.1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("firebrick", "navyblue")) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Receiver Operating Characteristic (ROC) Validation",
      subtitle = "Philips Connected Care Telemetry vs Standard Hospital NEWS2",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Clinical Strategy"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  # Calibration Plot (Deciles)
  test_data <- test_data %>%
    mutate(decile = ntile(pred_telemetry, 10))
  
  cal_df <- test_data %>%
    group_by(decile) %>%
    summarise(
      mean_pred = mean(pred_telemetry),
      obs_rate = mean(deterioration_event),
      n = n(),
      .groups = "drop"
    )
  
  p_cal <- ggplot(cal_df, aes(x = mean_pred, y = obs_rate)) +
    geom_point(size = 3, color = "navyblue") +
    geom_line(color = "navyblue", linetype = "solid") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "forestgreen") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Clinical Calibration Plot (Deciles of Risk)",
      subtitle = "Observed Deterioration Frequency vs Predicted Telemetry Risk",
      x = "Mean Predicted Probability",
      y = "Observed Deterioration Proportion"
    ) +
    theme(plot.title = element_text(face = "bold"))
  
  # Save ROC and Calibration plots
  ggsave(file.path(output_dir, "01_roc_curve.png"), plot = p_roc, width = 7, height = 5.5, dpi = 300)
  ggsave(file.path(output_dir, "01_calibration_curve.png"), plot = p_cal, width = 7, height = 5.5, dpi = 300)
  
  # 6. Decision Curve Analysis (DCA) - Addressing Alarm Fatigue
  # Net Benefit = (TP / N) - (FP / N) * (p_t / (1 - p_t))
  thresholds <- seq(0.02, 0.60, by = 0.02)
  N_test <- nrow(test_data)
  prevalence <- mean(test_data$deterioration_event)
  
  dca_list <- lapply(thresholds, function(pt) {
    w <- pt / (1 - pt)
    
    # Model: Telemetry
    tp_tel <- sum(test_data$pred_telemetry >= pt & test_data$deterioration_event == 1)
    fp_tel <- sum(test_data$pred_telemetry >= pt & test_data$deterioration_event == 0)
    nb_tel <- (tp_tel / N_test) - (fp_tel / N_test) * w
    
    # Model: NEWS2
    tp_n2 <- sum(test_data$pred_news2 >= pt & test_data$deterioration_event == 1)
    fp_n2 <- sum(test_data$pred_news2 >= pt & test_data$deterioration_event == 0)
    nb_n2 <- (tp_n2 / N_test) - (fp_n2 / N_test) * w
    
    # Treat All Strategy
    nb_all <- prevalence - (1 - prevalence) * w
    
    # Treat None Strategy
    nb_none <- 0
    
    tibble(
      threshold = pt,
      nb_telemetry = nb_tel,
      nb_news2 = nb_n2,
      nb_all = nb_all,
      nb_none = nb_none
    )
  })
  
  dca_df <- bind_rows(dca_list) %>%
    pivot_longer(cols = starts_with("nb_"), names_to = "strategy", values_to = "net_benefit") %>%
    mutate(
      Strategy = recode(strategy,
        "nb_telemetry" = "Philips Smart Telemetry",
        "nb_news2" = "Standard NEWS2 Benchmark",
        "nb_all" = "Escalate All Patients",
        "nb_none" = "Escalate No Patients"
      )
    )
  
  p_dca <- ggplot(dca_df, aes(x = threshold, y = net_benefit, color = Strategy, linetype = Strategy)) +
    geom_line(linewidth = 1.1) +
    coord_cartesian(ylim = c(-0.05, max(dca_df$net_benefit) * 1.1)) +
    scale_color_manual(values = c("forestgreen", "black", "firebrick", "navyblue")) +
    scale_linetype_manual(values = c("dotdash", "dotted", "dashed", "solid")) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Clinical Decision Curve Analysis (DCA): Mitigating Alarm Fatigue",
      subtitle = "Standardized Net Benefit across clinical intervention thresholds (Vickers & Elkin)",
      x = "Threshold Probability for Clinical Escalation (pt)",
      y = "Standardized Net Benefit"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "01_decision_curve_analysis.png"), 
         plot = p_dca, width = 9, height = 6, dpi = 300)
  
  # 7. Kaplan-Meier Survival Analysis by Risk Tier
  test_data <- test_data %>%
    mutate(
      risk_tier = case_when(
        pred_telemetry < 0.15 ~ "Low Risk (<15%)",
        pred_telemetry < 0.40 ~ "Moderate Risk (15-40%)",
        TRUE ~ "High Risk (>40%)"
      ),
      risk_tier = factor(risk_tier, levels = c("Low Risk (<15%)", "Moderate Risk (15-40%)", "High Risk (>40%)"))
    )
  
  surv_obj <- Surv(time = test_data$time_to_event_hrs, event = test_data$deterioration_event)
  km_fit <- survfit(surv_obj ~ risk_tier, data = test_data)
  
  # Base survival plot with ggplot
  km_summary <- summary(km_fit)
  km_df <- tibble(
    time = km_summary$time,
    surv = km_summary$surv,
    strata = km_summary$strata
  )
  
  p_km <- ggplot(km_df, aes(x = time, y = surv, color = strata)) +
    geom_step(linewidth = 1.1) +
    scale_color_manual(values = c("darkgreen", "orange3", "firebrick"),
                       labels = c("Low Risk", "Moderate Risk", "High Risk")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Time-to-Deterioration Kaplan-Meier Survival Analysis",
      subtitle = "Stratified by Philips Connected Care Predicted Telemetry Risk Tiers",
      x = "Inpatient Observation Time (Hours)",
      y = "Deterioration-Free Survival Probability",
      color = "Risk Stratification"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "01_survival_kaplan_meier.png"), 
         plot = p_km, width = 8, height = 5.5, dpi = 300)
  
  cat("  [Module 1] Saved Figures: 01_roc_curve.png, 01_calibration_curve.png, 01_decision_curve_analysis.png, 01_survival_kaplan_meier.png\n")
  
  return(list(
    fit_news2 = fit_news2,
    fit_telemetry = fit_telemetry,
    auc_news2 = roc_news2$auc,
    auc_telemetry = roc_telemetry$auc,
    brier_news2 = brier_news2,
    brier_telemetry = brier_telemetry
  ))
}
