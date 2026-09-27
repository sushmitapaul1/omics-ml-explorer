## =============================================================================
## Module 1 - Prostate cancer: protease panel classifiers (BPH vs PCa)
## Workflow of R_code_Splitting_Data.R + RCode_CaretPackage_Main_Final.R:
##   repeated random 80/20 splits -> caret classifiers tuned by repeated CV
##   (metric = ROC) -> training CV AUC and held-out test AUC per split.
## =============================================================================

prostate_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      data_source_ui(ns("data"), "Dummy protease data (BPH vs PCa, n = 120)"),
      shiny::hr(),
      shiny::selectInput(ns("class_col"), "Class column", choices = NULL),
      shiny::selectInput(ns("positive"), "Positive class", choices = NULL),
      shiny::radioButtons(ns("preset"), "Features", choices = c(
        "All proteases + PSA" = "all", "Proteases only" = "no_psa",
        "PSA only" = "psa", "Custom selection" = "custom")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'custom'", ns("preset")),
        shiny::selectizeInput(ns("features"), NULL, choices = NULL, multiple = TRUE)),
      shiny::checkboxGroupInput(ns("methods"), "Classifiers",
        choices = names(CLASSIFIERS),
        selected = c("SVM (linear, class-weighted)", "Random forest", "MLP (nnet)",
                     "KNN", "Naive Bayes", "GBM")),
      bslib::accordion(open = FALSE, bslib::accordion_panel(
        "Validation settings", icon = shiny::icon("sliders"),
        shiny::sliderInput(ns("n_splits"), "Random train/test splits", 1, 10, 3, step = 1),
        shiny::sliderInput(ns("p_train"), "Training fraction", 0.6, 0.9, 0.8, step = 0.05),
        shiny::checkboxInput(ns("stratified"), "Stratify splits by class", TRUE),
        shiny::sliderInput(ns("folds"), "CV folds", 3, 10, 5, step = 1),
        shiny::sliderInput(ns("repeats"), "CV repeats", 1, 10, 2, step = 1),
        shiny::sliderInput(ns("tune"), "Tuning values per parameter", 1, 5, 3, step = 1),
        shiny::numericInput(ns("seed"), "Random seed", 123)
      )),
      shiny::actionButton(ns("run"), "Train classifiers", icon = shiny::icon("play"),
                          class = "btn-primary w-100"),
      shiny::helpText("The original analysis used 10 splits and 10 CV repeats. ",
                      "Fewer are used here by default so a run takes seconds.")
    ),
    bslib::navset_card_tab(
      id = ns("tabs"),
      bslib::nav_panel("Data", icon = shiny::icon("table"),
        bslib::layout_columns(col_widths = c(5, 7),
          shiny::plotOutput(ns("class_plot"), height = 260, fill = FALSE),
          shiny::plotOutput(ns("feature_plot"), height = 260, fill = FALSE)),
        DT::DTOutput(ns("data_table"), fill = FALSE)),
      bslib::nav_panel("Performance", icon = shiny::icon("chart-simple"), value = "perf",
        shiny::uiOutput(ns("perf_hint")),
        DT::DTOutput(ns("summary_table"), fill = FALSE),
        shiny::plotOutput(ns("auc_plot"), height = 380, fill = FALSE)),
      bslib::nav_panel("ROC curves", icon = shiny::icon("chart-line"),
        shiny::uiOutput(ns("split_picker")),
        shiny::plotOutput(ns("roc_plot"), height = 480, fill = FALSE)),
      bslib::nav_panel("All results", icon = shiny::icon("list"),
        shiny::downloadButton(ns("dl_metrics"), "Download per-split metrics (CSV)", class = "btn-sm mb-2 align-self-start"),
        DT::DTOutput(ns("metrics_table"), fill = FALSE)),
      bslib::nav_panel("Method", icon = shiny::icon("book"), method_note(
        "<h5>What this module does</h5>",
        "<ol><li>Splits the samples into training and test sets at random (default 80/20). ",
        "This is repeated several times, as in <code>R_code_Splitting_Data.R</code>.</li>",
        "<li>On each training set, every selected classifier is tuned with repeated k-fold ",
        "cross-validation using <code>caret::train</code> with <code>metric = \"ROC\"</code>.</li>",
        "<li>The tuned model predicts the untouched test set. Test AUC, accuracy, sensitivity ",
        "and specificity are reported, together with the cross-validated training AUC.</li></ol>",
        "<h5>Differences from the original script</h5><ul>",
        "<li>ROC curves and AUC use <code>pROC</code> instead of <code>MLeval</code>.</li>",
        "<li>The MLP uses <code>nnet</code> (one hidden layer) instead of <code>RSNNS::mlp</code>.</li>",
        "<li>Splits are stratified by class by default; untick the option to use plain random ",
        "splits as in the original.</li></ul>",
        "<p><b>Input format:</b> one row per sample, numeric feature columns, and a two-level ",
        "class column (for example <code>BPH</code> / <code>PCa</code>). An ID column is ignored ",
        "if it is not numeric.</p>"))
    )
  )
}

prostate_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dummy <- "data/prostate_protease_dummy.csv"
    dat <- data_source_server(input, output, dummy)

    num_cols <- shiny::reactive(names(dat())[vapply(dat(), is.numeric, logical(1))])

    shiny::observe({
      d <- dat()
      cand <- names(d)[vapply(d, function(v) length(unique(v)) == 2, logical(1))]
      sel <- if ("Class" %in% cand) "Class" else cand[1]
      shiny::updateSelectInput(session, "class_col", choices = cand, selected = sel)
      shiny::updateSelectizeInput(session, "features", choices = num_cols(), selected = num_cols())
    })
    shiny::observe({
      shiny::req(input$class_col %in% names(dat()))
      lv <- sort(unique(as.character(dat()[[input$class_col]])))
      sel <- if ("PCa" %in% lv) "PCa" else lv[length(lv)]
      shiny::updateSelectInput(session, "positive", choices = lv, selected = sel)
    })

    features <- shiny::reactive({
      nc <- setdiff(num_cols(), input$class_col)
      psa <- grep("^PSA$", nc, ignore.case = TRUE, value = TRUE)
      switch(input$preset,
             all = nc,
             no_psa = setdiff(nc, psa),
             psa = if (length(psa)) psa else character(0),
             custom = intersect(input$features, nc))
    })

    output$class_plot <- shiny::renderPlot({
      d <- dat(); shiny::req(input$class_col %in% names(d))
      ggplot2::ggplot(d, ggplot2::aes(x = .data[[input$class_col]], fill = .data[[input$class_col]])) +
        ggplot2::geom_bar(width = 0.6, show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = c("#9ecae1", APP_PRIMARY, "#fdae6b", "#bcbddc")) +
        ggplot2::labs(x = NULL, y = "Samples", title = "Class balance") + theme_app()
    }, res = 96)
    output$feature_plot <- shiny::renderPlot({
      d <- dat(); shiny::req(input$class_col %in% names(d), length(features()) > 0)
      f <- head(features(), 12)
      long <- do.call(rbind, lapply(f, function(v) data.frame(
        feature = v, value = as.numeric(scale(d[[v]])), class = d[[input$class_col]])))
      ggplot2::ggplot(long, ggplot2::aes(feature, value, fill = class)) +
        ggplot2::geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
        ggplot2::scale_fill_manual(values = c("#9ecae1", APP_PRIMARY)) +
        ggplot2::labs(x = NULL, y = "z-score", fill = NULL, title = "Selected features by class") +
        theme_app(11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    }, res = 96)
    output$data_table <- DT::renderDT(dt_table(dat(), page = 8))

    res <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run, {
      if (length(input$methods) == 0) {
        shiny::showNotification("Select at least one classifier.", type = "warning"); return()
      }
      if (!isTRUE(input$class_col %in% names(dat()))) return()
      if (length(features()) == 0) {
        shiny::showNotification("No features selected (is there a PSA column?).", type = "warning"); return()
      }
      d <- dat()
      out <- shiny::withProgress(message = "Training classifiers", value = 0, {
        tryCatch(run_classifiers(
          d, features(), input$class_col, input$positive, input$methods,
          n_splits = input$n_splits, p_train = input$p_train, stratified = input$stratified,
          cv_folds = input$folds, cv_repeats = input$repeats, tune_length = input$tune,
          seed = input$seed, progress = function(v, detail) shiny::setProgress(v, detail = detail)),
          error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); NULL })
      })
      res(out)
      bslib::nav_select("tabs", "perf")
    })

    output$perf_hint <- shiny::renderUI({
      if (is.null(res())) return(shiny::div(class = "text-muted p-3",
        shiny::icon("circle-info"), " Choose features and classifiers, then press ",
        shiny::strong("Train classifiers"), "."))
      r <- res()
      shiny::p(class = "small text-muted", sprintf(
        "%d feature(s): %s. Mean ± SD over %d split(s); positive class = %s.",
        length(r$features), paste(head(r$features, 12), collapse = ", "),
        max(r$metrics$split), r$positive))
    })
    output$summary_table <- DT::renderDT({
      shiny::req(res()); DT::datatable(summarise_metrics(res()$metrics), rownames = FALSE,
                                       options = list(dom = "t", pageLength = 20))
    })
    output$auc_plot <- shiny::renderPlot({ shiny::req(res()); plot_auc_box(res()$metrics) }, res = 96)
    output$split_picker <- shiny::renderUI({
      shiny::req(res())
      shiny::selectInput(ns("roc_split"), "Show split", choices = unique(res()$metrics$split), width = 160)
    })
    output$roc_plot <- shiny::renderPlot({
      shiny::req(res(), input$roc_split)
      r <- res(); s <- as.integer(input$roc_split)
      plot_roc_curves(r$roc[r$roc$split == s, ], r$metrics[r$metrics$split == s, ],
                      sprintf("ROC curves, split %d", s))
    }, res = 96)
    output$metrics_table <- DT::renderDT({ shiny::req(res()); dt_table(res()$metrics, 15) })
    output$dl_metrics <- shiny::downloadHandler("prostate_classifier_metrics.csv",
      function(f) write.csv(res()$metrics, f, row.names = FALSE))
  })
}
