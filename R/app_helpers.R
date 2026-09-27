## =============================================================================
## Shared UI / plotting helpers
## =============================================================================

APP_PRIMARY <- "#1f6f78"

CLASSIFIER_COLORS <- c(
  "SVM (linear, class-weighted)" = "#1f77b4", "SVM (radial)" = "#17becf",
  "Random forest" = "#2ca02c", "MLP (nnet)" = "#9467bd", "KNN" = "#ff7f0e",
  "Naive Bayes" = "#8c564b", "GBM" = "#d62728")

theme_app <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

## Data-source picker: bundled dummy data or a CSV upload
data_source_ui <- function(id, dummy_label) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::radioButtons(ns("source"), "Data", inline = FALSE,
                        choices = c(setNames("dummy", dummy_label), "Upload my own CSV" = "upload")),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'upload'", ns("source")),
      shiny::fileInput(ns("file"), NULL, accept = ".csv", placeholder = "Choose a .csv file")),
    shiny::downloadLink(ns("dl_dummy"), shiny::tagList(shiny::icon("download"),
                        "Download the dummy CSV (shows the format)"), class = "small")
  )
}

## Call inside a module server whose UI used data_source_ui(ns("data"), ...)
data_source_server <- function(input, output, dummy_path) {
  output[["data-dl_dummy"]] <- shiny::downloadHandler(
    basename(dummy_path), function(f) file.copy(dummy_path, f))
  shiny::reactive({
    if (identical(input[["data-source"]], "upload")) {
      shiny::req(input[["data-file"]])
      d <- tryCatch(read.csv(input[["data-file"]]$datapath, check.names = FALSE,
                             stringsAsFactors = FALSE), error = function(e) NULL)
      shiny::validate(shiny::need(!is.null(d) && nrow(d) > 0, "Could not read the uploaded CSV."))
      d
    } else {
      read.csv(dummy_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  })
}

dt_table <- function(df, page = 10, digits = 3) {
  num <- names(df)[vapply(df, is.double, logical(1))]
  t <- DT::datatable(df, rownames = FALSE, filter = "none",
                     options = list(pageLength = page, scrollX = TRUE, dom = "ftip"))
  if (length(num)) t <- DT::formatRound(t, num, digits)
  t
}

SET_LEVELS <- c("Training (cross-validation)", "Test (held-out)")

## ROC curves coloured by classifier, one panel per set, AUCs printed in each panel
plot_roc_curves <- function(roc_df, metrics, title = NULL) {
  roc_df$set <- factor(roc_df$set, levels = SET_LEVELS)
  roc_df <- roc_df[order(roc_df$classifier, roc_df$set, roc_df$fpr, roc_df$tpr), ]
  key <- rbind(
    data.frame(classifier = metrics$classifier, set = SET_LEVELS[1], auc = metrics$cv_pooled_auc),
    data.frame(classifier = metrics$classifier, set = SET_LEVELS[2], auc = metrics$test_auc))
  key <- key[key$set %in% as.character(unique(roc_df$set)), ]
  key$set <- factor(key$set, levels = SET_LEVELS)
  key <- key[order(key$set, -key$auc), ]
  key$y <- ave(seq_len(nrow(key)), key$set, FUN = function(i) 0.04 + 0.075 * (length(i) - seq_along(i)))
  key$label <- sprintf("%s  %.2f", key$classifier, key$auc)
  ggplot2::ggplot(roc_df, ggplot2::aes(fpr, tpr, color = classifier)) +
    ggplot2::geom_abline(linetype = "dashed", color = "grey60") +
    ggplot2::geom_path(linewidth = 0.9) +
    ggplot2::geom_text(data = key, ggplot2::aes(x = 1, y = y, label = label), hjust = 1,
                       size = 3.4, fontface = "bold", show.legend = FALSE) +
    ggplot2::annotate("text", x = 1, y = 0.04 + 0.075 * length(unique(key$classifier)),
                      label = "AUC", hjust = 1, size = 3.4, color = "grey40") +
    ggplot2::facet_wrap(~set, drop = TRUE) +
    ggplot2::scale_color_manual(values = CLASSIFIER_COLORS) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "False positive rate (1 - specificity)",
                  y = "True positive rate (sensitivity)", color = NULL, title = title) +
    theme_app() + ggplot2::theme(legend.position = "none")
}

plot_auc_box <- function(metrics) {
  long <- rbind(
    data.frame(classifier = metrics$classifier, set = "Training (CV)", auc = metrics$cv_auc),
    data.frame(classifier = metrics$classifier, set = "Test (held-out)", auc = metrics$test_auc))
  # reversed so that, after coord_flip, Training sits above Test for each classifier
  long$set <- factor(long$set, levels = c("Test (held-out)", "Training (CV)"))
  ord <- names(sort(tapply(metrics$test_auc, metrics$classifier, mean)))
  long$classifier <- factor(long$classifier, levels = ord)
  ggplot2::ggplot(long, ggplot2::aes(classifier, auc, fill = set)) +
    ggplot2::geom_boxplot(alpha = 0.55, outlier.shape = NA, position = ggplot2::position_dodge(0.8)) +
    ggplot2::geom_point(ggplot2::aes(color = set), position = ggplot2::position_jitterdodge(0.1, dodge.width = 0.8),
                        size = 1.8, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
    ggplot2::scale_fill_manual(values = c("Training (CV)" = "#9ecae1", "Test (held-out)" = APP_PRIMARY)) +
    ggplot2::scale_color_manual(values = c("Training (CV)" = "#3182bd", "Test (held-out)" = "#0b3d42")) +
    ggplot2::coord_flip(ylim = c(min(0.4, min(long$auc, na.rm = TRUE)), 1)) +
    ggplot2::labs(x = NULL, y = "ROC AUC", fill = NULL,
                  title = "AUC across random train/test splits") + theme_app() +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
}

method_note <- function(...) {
  shiny::div(class = "method-note", shiny::HTML(paste0(...)))
}
