## =============================================================================
## Module 2 - PDAC: ML feature selection + classification on the consensus panel
## Workflow of R_Code_ML_FeatureSelection.R + R_Code_For_Classifier.R
## =============================================================================

pdac_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      data_source_ui(ns("data"), "Dummy Olink NPX data (PDAC, 25 vs 25, 80 proteins)"),
      shiny::hr(),
      shiny::selectInput(ns("group_col"), "Group column", choices = NULL),
      shiny::selectInput(ns("positive"), "Positive class", choices = NULL),
      shiny::selectizeInput(ns("exclude"), "Exclude columns (IDs, grouping markers)",
                            choices = NULL, multiple = TRUE),
      shiny::radioButtons(ns("fs_scope"), "Run feature selection on", choices = c(
        "Training split only (recommended)" = "train",
        "All samples (as in the original script)" = "all")),
      bslib::accordion(open = FALSE,
        bslib::accordion_panel("Step 1 - Feature selection", icon = shiny::icon("filter"),
          shiny::checkboxGroupInput(ns("fs_methods"), NULL, choices = FS_METHODS, selected = FS_METHODS),
          shiny::sliderInput(ns("n_boot"), "Stability: resamples", 20, 300, 100, step = 10),
          shiny::sliderInput(ns("alpha"), "Stability: elastic-net alpha", 0.1, 1, 0.7, step = 0.1),
          shiny::sliderInput(ns("stab_thr"), "Stability: selection frequency cut-off", 0.3, 0.9, 0.6, step = 0.05),
          shiny::sliderInput(ns("boruta_runs"), "Boruta: max runs", 30, 300, 100, step = 10),
          shiny::sliderInput(ns("rfe_repeats"), "RFE: CV repeats (5-fold)", 1, 5, 2, step = 1),
          shiny::sliderInput(ns("k_mrmr"), "mRMR: top k", 5, 50, 20, step = 1)),
        bslib::accordion_panel("Step 2 - Classification", icon = shiny::icon("sliders"),
          shiny::checkboxGroupInput(ns("methods"), "Classifiers", choices = names(CLASSIFIERS),
            selected = c("SVM (radial)", "Random forest", "MLP (nnet)", "KNN", "Naive Bayes", "GBM")),
          shiny::sliderInput(ns("folds"), "CV folds", 3, 10, 5, step = 1),
          shiny::sliderInput(ns("repeats"), "CV repeats", 1, 10, 3, step = 1),
          shiny::numericInput(ns("seed"), "Random seed", 42))
      ),
      shiny::actionButton(ns("run_fs"), "1. Run feature selection", icon = shiny::icon("filter"),
                          class = "btn-primary w-100"),
      shiny::uiOutput(ns("panel_ui")),
      shiny::actionButton(ns("run_clf"), "2. Classify with this panel", icon = shiny::icon("play"),
                          class = "btn-outline-primary w-100")
    ),
    bslib::navset_card_tab(
      id = ns("tabs"),
      bslib::nav_panel("Data", icon = shiny::icon("table"),
        shiny::uiOutput(ns("data_summary")), DT::DTOutput(ns("data_table"), fill = FALSE)),
      bslib::nav_panel("Feature selection", icon = shiny::icon("filter"), value = "fs",
        shiny::uiOutput(ns("fs_hint")),
        bslib::layout_columns(col_widths = c(6, 6),
          shiny::plotOutput(ns("votes_plot"), height = 560, fill = FALSE),
          shiny::div(shiny::plotOutput(ns("stab_plot"), height = 270, fill = FALSE),
                     shiny::plotOutput(ns("rfe_plot"), height = 270, fill = FALSE)))),
      bslib::nav_panel("Consensus table", icon = shiny::icon("list-check"),
        shiny::downloadButton(ns("dl_cons"), "Download consensus table (CSV)", class = "btn-sm mb-2 align-self-start"),
        DT::DTOutput(ns("cons_table"), fill = FALSE)),
      bslib::nav_panel("Classification", icon = shiny::icon("chart-line"), value = "clf",
        shiny::uiOutput(ns("clf_hint")),
        DT::DTOutput(ns("clf_table"), fill = FALSE),
        bslib::layout_columns(col_widths = c(5, 7),
          shiny::plotOutput(ns("cv_plot"), height = 420, fill = FALSE),
          shiny::plotOutput(ns("roc_plot"), height = 420, fill = FALSE)),
        shiny::downloadButton(ns("dl_clf"), "Download test-set metrics (CSV)", class = "btn-sm w-auto align-self-start")),
      bslib::nav_panel("Method", icon = shiny::icon("book"), method_note(
        "<h5>Step 1: feature selection (four complementary methods)</h5><ol>",
        "<li><b>Stability selection</b>: an elastic-net logistic regression (<code>glmnet</code>, ",
        "lambda.1se) is refit on many 80% subsamples. A protein is kept if it is selected in at ",
        "least the chosen fraction of resamples.</li>",
        "<li><b>Boruta</b>: an all-relevant random-forest wrapper that compares each protein ",
        "with shuffled \"shadow\" copies.</li>",
        "<li><b>RFE-CV</b>: recursive feature elimination with random forests and repeated ",
        "5-fold CV (<code>caret::rfe</code>).</li>",
        "<li><b>mRMR</b>: minimum redundancy, maximum relevance (<code>praznik</code>), top k.</li></ol>",
        "<p>The <b>consensus</b> counts how many methods chose each protein and ranks proteins ",
        "with Robust Rank Aggregation (<code>RobustRankAggreg</code>). The panel is every protein ",
        "chosen by at least the number of methods you set. You can also edit it by hand.</p>",
        "<h5>Step 2: classification</h5><p>Six classifiers are tuned with repeated k-fold CV on ",
        "the training split (metric = ROC). They are then evaluated once on the held-out 20%, ",
        "as in <code>R_Code_For_Classifier.R</code>.</p>",
        "<h5>A note on selection bias</h5><p>If proteins are selected using <i>all</i> samples, ",
        "the test samples have already influenced which proteins were chosen, so the test AUC ",
        "comes out too optimistic. The default here runs feature selection on the training ",
        "split only. The held-out samples are then genuinely unseen. Choose \"All samples\" ",
        "to reproduce the original workflow.</p>",
        "<p><b>Input format:</b> one row per patient, a group column with two levels, and one ",
        "numeric column per protein (for example Olink NPX values). Exclude ID columns and any ",
        "marker that defines the groups (here <code>CTSG</code>), because including it would be circular.</p>"))
    )
  )
}

pdac_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dummy <- "data/pdac_olink_dummy.csv"
    dat <- data_source_server(input, output, dummy)

    shiny::observe({
      d <- dat()
      cand <- names(d)[vapply(d, function(v) length(unique(v)) == 2, logical(1))]
      shiny::updateSelectInput(session, "group_col", choices = cand,
                               selected = if ("Group" %in% cand) "Group" else cand[1])
      default_ex <- intersect(c("Patient ID", "SampleID", "ID", "CTSG"), names(d))
      shiny::updateSelectizeInput(session, "exclude", choices = names(d), selected = default_ex)
    })
    shiny::observe({
      shiny::req(input$group_col %in% names(dat()))
      lv <- sort(unique(as.character(dat()[[input$group_col]])))
      hi <- grep("high|short|pos|case", lv, ignore.case = TRUE, value = TRUE)
      shiny::updateSelectInput(session, "positive", choices = lv,
                               selected = if (length(hi)) hi[1] else lv[1])
    })

    proteins <- shiny::reactive({
      d <- dat(); shiny::req(input$group_col %in% names(d))
      num <- names(d)[vapply(d, is.numeric, logical(1))]
      setdiff(num, c(input$group_col, input$exclude))
    })

    output$data_summary <- shiny::renderUI({
      d <- dat(); shiny::req(input$group_col %in% names(d))
      tb <- table(d[[input$group_col]])
      bslib::layout_columns(col_widths = c(4, 4, 4),
        bslib::value_box("Samples", nrow(d), showcase = shiny::icon("users"), theme = "primary"),
        bslib::value_box("Protein features", length(proteins()), showcase = shiny::icon("dna"),
                         theme = "secondary"),
        bslib::value_box("Groups", shiny::HTML(paste(names(tb), tb, sep = ": ", collapse = "<br>")),
                         showcase = shiny::icon("scale-balanced"), theme = "light"))
    })
    output$data_table <- DT::renderDT(dt_table(dat()[, 1:min(20, ncol(dat()))], page = 8))

    ## ---- Step 1: feature selection -------------------------------------------
    fs <- shiny::reactiveVal(NULL)
    split_idx <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run_fs, {
      if (length(input$fs_methods) < 1) {
        shiny::showNotification("Select at least one feature-selection method.", type = "warning"); return()
      }
      d <- dat(); y <- prep_class(d[[input$group_col]])
      set.seed(input$seed)
      idx <- split_indices(y, 0.8, stratified = TRUE)
      split_idx(idx)
      rows <- if (input$fs_scope == "train") idx else seq_len(nrow(d))
      out <- shiny::withProgress(message = "Feature selection", value = 0, {
        tryCatch(run_feature_selection(
          d[rows, proteins(), drop = FALSE], d[[input$group_col]][rows],
          methods = input$fs_methods, n_boot = input$n_boot, alpha = input$alpha,
          stab_threshold = input$stab_thr, boruta_max_runs = input$boruta_runs,
          rfe_folds = 5, rfe_repeats = input$rfe_repeats, k_mrmr = input$k_mrmr,
          seed = input$seed, progress = function(v, detail) shiny::setProgress(v, detail = detail)),
          error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); NULL })
      })
      fs(out); clf(NULL)
      bslib::nav_select("tabs", "fs")
    })

    output$panel_ui <- shiny::renderUI({
      f <- fs()
      if (is.null(f)) return(shiny::helpText("Run step 1 to build the consensus panel."))
      n_m <- length(f$methods)
      shiny::tagList(
        shiny::sliderInput(ns("min_votes"), "Panel: proteins chosen by at least",
                           1, n_m, min(3, n_m), step = 1, post = " method(s)"),
        shiny::selectizeInput(ns("panel"), "Panel (editable)", choices = f$consensus$protein,
                              selected = NULL, multiple = TRUE))
    })
    shiny::observeEvent(list(input$min_votes, fs()), {
      f <- fs(); shiny::req(f, input$min_votes)
      sel <- f$consensus$protein[f$consensus$n_methods_selected >= input$min_votes]
      shiny::updateSelectizeInput(session, "panel", choices = f$consensus$protein, selected = sel)
    })

    output$fs_hint <- shiny::renderUI({
      f <- fs()
      if (is.null(f)) return(shiny::div(class = "text-muted p-3", shiny::icon("circle-info"),
        " Press ", shiny::strong("1. Run feature selection"), ". With the default settings this takes under a minute."))
      n_sel <- vapply(f$selected, length, integer(1))
      shiny::p(class = "small text-muted",
        paste0("Selected per method: ", paste(names(n_sel), n_sel, sep = " = ", collapse = "; "),
               ". Scope: ", if (input$fs_scope == "train") "training split only" else "all samples", "."))
    })
    output$votes_plot <- shiny::renderPlot({
      f <- fs(); shiny::req(f)
      top <- head(f$consensus[f$consensus$n_methods_selected > 0, ], 25)
      top$protein <- factor(top$protein, levels = rev(top$protein))
      ggplot2::ggplot(top, ggplot2::aes(protein, n_methods_selected,
                                        fill = n_methods_selected >= (input$min_votes %||% 3))) +
        ggplot2::geom_col(width = 0.75) + ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(`TRUE` = APP_PRIMARY, `FALSE` = "#b8c7c9"),
                                   labels = c(`TRUE` = "In panel", `FALSE` = "Not in panel"), name = NULL) +
        ggplot2::scale_y_continuous(breaks = 0:4) +
        ggplot2::labs(x = NULL, y = "Number of methods selecting the protein",
                      title = "Consensus votes (top 25, ordered by votes then RRA rank)") + theme_app(12)
    }, res = 96)
    output$stab_plot <- shiny::renderPlot({
      f <- fs(); shiny::req(f, f$stability)
      s <- head(f$stability, 20); s$protein <- factor(s$protein, levels = rev(s$protein))
      ggplot2::ggplot(s, ggplot2::aes(protein, stability_freq)) +
        ggplot2::geom_segment(ggplot2::aes(xend = protein, y = 0, yend = stability_freq), color = "grey70") +
        ggplot2::geom_point(color = APP_PRIMARY, size = 2.5) +
        ggplot2::geom_hline(yintercept = input$stab_thr, linetype = "dashed", color = "#d73027") +
        ggplot2::coord_flip() + ggplot2::ylim(0, 1) +
        ggplot2::labs(x = NULL, y = "Selection frequency", title = "Stability selection (top 20)") + theme_app(11)
    }, res = 96)
    output$rfe_plot <- shiny::renderPlot({
      f <- fs(); shiny::req(f, f$rfe_profile)
      ggplot2::ggplot(f$rfe_profile, ggplot2::aes(Variables, Accuracy)) +
        ggplot2::geom_line(color = "#4575b4") + ggplot2::geom_point(color = "#4575b4") +
        ggplot2::geom_vline(xintercept = f$rfe_optsize, linetype = "dashed", color = "#d73027") +
        ggplot2::labs(title = "RFE-CV: accuracy vs number of proteins",
                      subtitle = paste("Optimal size =", f$rfe_optsize),
                      x = "Number of proteins", y = "CV accuracy") + theme_app(11)
    }, res = 96)
    output$cons_table <- DT::renderDT({
      f <- fs(); shiny::req(f)
      t <- dt_table(f$consensus, 15, digits = 4)
      DT::formatSignif(t, "rra_score", 3)
    })
    output$dl_cons <- shiny::downloadHandler("pdac_feature_selection_consensus.csv",
      function(file) write.csv(fs()$consensus, file, row.names = FALSE))

    ## ---- Step 2: classification on the panel ------------------------------------
    clf <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$run_clf, {
      if (is.null(fs())) { shiny::showNotification("Run step 1 first.", type = "warning"); return() }
      if (length(input$panel) < 1) { shiny::showNotification("The panel is empty: lower the vote threshold or add proteins.", type = "warning"); return() }
      if (length(input$methods) < 1) { shiny::showNotification("Select at least one classifier.", type = "warning"); return() }
      d <- dat()
      out <- shiny::withProgress(message = "Training classifiers", value = 0, {
        tryCatch(run_classifiers(
          d, input$panel, input$group_col, input$positive, input$methods,
          cv_folds = input$folds, cv_repeats = input$repeats, seed = input$seed,
          train_idx = split_idx(),
          progress = function(v, detail) shiny::setProgress(v, detail = detail)),
          error = function(e) { shiny::showNotification(conditionMessage(e), type = "error"); NULL })
      })
      clf(out)
      bslib::nav_select("tabs", "clf")
    })
    output$clf_hint <- shiny::renderUI({
      r <- clf()
      if (is.null(r)) return(shiny::div(class = "text-muted p-3", shiny::icon("circle-info"),
        " After step 1, press ", shiny::strong("2. Classify with this panel"), "."))
      shiny::p(class = "small text-muted", sprintf(
        "Panel (%d proteins): %s. Train n = %d, held-out test n = %d. Positive class: %s.",
        length(r$features), paste(r$features, collapse = ", "), r$metrics$n_train[1],
        r$metrics$n_test[1], r$positive))
    })
    output$clf_table <- DT::renderDT({
      r <- clf(); shiny::req(r)
      m <- r$metrics[order(-r$metrics$test_auc),
                     c("classifier", "cv_auc", "test_auc", "accuracy", "sensitivity", "specificity", "f1", "kappa", "best_tune")]
      names(m) <- c("Classifier", "CV AUC (train)", "Test AUC", "Accuracy", "Sensitivity",
                    "Specificity", "F1", "Kappa", "Best tuning")
      DT::formatRound(DT::datatable(m, rownames = FALSE, options = list(dom = "t", scrollX = TRUE)),
                      2:8, 3)
    })
    output$cv_plot <- shiny::renderPlot({
      r <- clf(); shiny::req(r)
      rs <- r$resamples
      rs$classifier <- stats::reorder(rs$classifier, rs$ROC, FUN = median)
      ggplot2::ggplot(rs, ggplot2::aes(classifier, ROC, fill = classifier)) +
        ggplot2::geom_boxplot(alpha = 0.7, show.legend = FALSE, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.12, size = 1, alpha = 0.5, show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = CLASSIFIER_COLORS) + ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "ROC AUC per CV fold",
                      title = "CV AUC per fold (training)") + theme_app(12)
    }, res = 96)
    output$roc_plot <- shiny::renderPlot({
      r <- clf(); shiny::req(r)
      plot_roc_curves(r$roc[r$roc$set == "Test (held-out)", ], r$metrics,
                      "Held-out test set ROC curves")
    }, res = 96)
    output$dl_clf <- shiny::downloadHandler("pdac_classifier_test_performance.csv",
      function(file) write.csv(clf()$metrics, file, row.names = FALSE))
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
