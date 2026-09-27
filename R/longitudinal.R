## =============================================================================
## Longitudinal biomarker analysis (from R_Code_Logintudinal_data_set_Final_code.R)
## Thresholds come from the Control group: control mean ± z * residual SD of
## lm(Biomarker ~ days + ClinicalNumber) fitted on Control samples only.
## Every measurement is then flagged as Significantly High / Low / Normal.
## =============================================================================

UVEAL_REQUIRED <- c("SampleID", "Class", "ClinicalNumber", "days", "Biomarker")

compute_thresholds <- function(df, control = "Control", z = 1.96) {
  ctrl <- df[df$Class == control, ]
  if (nrow(ctrl) < 4) stop("Need at least 4 Control measurements to estimate thresholds.")
  fit <- lm(Biomarker ~ days + ClinicalNumber, data = ctrl)
  sd_res <- sigma(fit)
  avg <- mean(ctrl$Biomarker)
  list(lower = avg - z * sd_res, upper = avg + z * sd_res, mean = avg,
       sd = sd_res, z = z, model = fit)
}

flag_points <- function(df, th) {
  df$Status <- ifelse(df$Biomarker > th$upper, "Significantly High",
                      ifelse(df$Biomarker < th$lower, "Significantly Low", "Normal"))
  df
}

patient_summary <- function(df) {
  out <- do.call(rbind, lapply(split(df, df$SampleID), function(d) {
    d <- d[order(d$days), ]
    slope <- if (nrow(d) >= 2) coef(lm(Biomarker ~ days, data = d))[["days"]] else NA
    hi <- d$days[d$Status == "Significantly High"]
    data.frame(
      Patient = d$SampleID[1], Class = d$Class[1], Visits = nrow(d),
      `Follow-up (days)` = max(d$days) - min(d$days),
      `First value` = round(d$Biomarker[1], 2), `Last value` = round(tail(d$Biomarker, 1), 2),
      `Slope per 100 days` = round(100 * slope, 3),
      `High points` = sum(d$Status == "Significantly High"),
      `Low points` = sum(d$Status == "Significantly Low"),
      `First high (day)` = if (length(hi)) min(hi) else NA,
      check.names = FALSE)
  }))
  rownames(out) <- NULL
  out[order(out$Class, -out$`High points`), ]
}

plot_longitudinal <- function(df, th, show_labels = TRUE, highlight = NULL) {
  df$ClinicalNumber_Fact <- factor(df$ClinicalNumber)
  df$Status <- factor(df$Status, levels = c("Significantly High", "Significantly Low", "Normal"))
  if (!"Label" %in% names(df)) df$Label <- df$SampleID
  last <- do.call(rbind, lapply(split(df, df$SampleID), function(d) d[which.max(d$days), ]))
  n_vis <- nlevels(df$ClinicalNumber_Fact)
  pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
  pal <- setNames(rep(pal, length.out = n_vis), levels(df$ClinicalNumber_Fact))

  p <- ggplot2::ggplot(df, ggplot2::aes(days, Biomarker)) +
    ggplot2::facet_wrap(~Class) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = th$lower, ymax = th$upper,
                      fill = "grey70", alpha = 0.12) +
    ggplot2::geom_line(ggplot2::aes(group = SampleID), color = "grey55", alpha = 0.45) +
    ggplot2::geom_hline(yintercept = c(th$lower, th$upper), linetype = "dashed",
                        color = "black", alpha = 0.35)
  if (length(highlight)) {
    p <- p + ggplot2::geom_line(data = df[df$SampleID %in% highlight, ],
                                ggplot2::aes(group = SampleID), color = "black", linewidth = 1.1)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = ClinicalNumber_Fact, fill = ClinicalNumber_Fact,
                                     shape = Status, size = Status), stroke = 1) +
    ggplot2::scale_shape_manual(values = c("Significantly High" = 17, "Significantly Low" = 25,
                                           "Normal" = 16), drop = FALSE) +
    ggplot2::scale_size_manual(values = c("Significantly High" = 4, "Significantly Low" = 4,
                                          "Normal" = 2), drop = FALSE) +
    ggplot2::scale_color_manual(values = pal) + ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(color = "black"),
                   strip.text = ggplot2::element_text(size = 13, face = "bold"),
                   legend.position = "right") +
    ggplot2::labs(title = "Longitudinal Biomarker Analysis",
                  subtitle = sprintf("Control-derived thresholds (mean +/- %.2f x residual SD): [%.2f, %.2f]",
                                     th$z, th$lower, th$upper),
                  x = "Days", y = "Biomarker value", color = "Clinical visit",
                  fill = "Clinical visit", shape = "Significance") +
    ggplot2::guides(size = "none")
  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(data = last, ggplot2::aes(label = Label), size = 3,
                                      fontface = "bold", box.padding = 0.5, seed = 1)
  }
  p
}
