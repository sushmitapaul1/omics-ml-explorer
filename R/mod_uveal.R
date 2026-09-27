## =============================================================================
## Module 3 - Uveal melanoma: longitudinal biomarker monitoring
## Based on R_Code_Logintudinal_data_set_Final_code.R
## =============================================================================

uveal_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      data_source_ui(ns("data"), "Dummy data (10 Control + 10 Metastasized, 4 visits)"),
      shiny::hr(),
      shiny::selectInput(ns("control"), "Reference (control) group", choices = NULL),
      shiny::selectizeInput(ns("groups"), "Groups to show", choices = NULL, multiple = TRUE),
      shiny::sliderInput(ns("z"), "Threshold width (z x residual SD)", 1, 3, 1.96, step = 0.01),
      shiny::helpText("1.96 is roughly a 95% reference interval, 2.58 roughly 99%."),
      shiny::selectizeInput(ns("highlight"), "Highlight patients", choices = NULL, multiple = TRUE),
      shiny::checkboxInput(ns("labels"), "Label each patient's last visit", TRUE),
      shiny::downloadButton(ns("dl_plot"), "Download plot (PNG)", class = "btn-sm w-100 mb-1"),
      shiny::downloadButton(ns("dl_flags"), "Download flagged data (CSV)", class = "btn-sm w-100")
    ),
    bslib::navset_card_tab(
      bslib::nav_panel("Trajectories", icon = shiny::icon("chart-line"),
        shiny::uiOutput(ns("boxes")),
        shiny::plotOutput(ns("plot"), height = 540, fill = FALSE)),
      bslib::nav_panel("Patient summary", icon = shiny::icon("user"),
        shiny::p(class = "small text-muted",
          "One row per patient: number of measurements outside the control-derived thresholds, ",
          "the day of the first high value, and the patient's own linear trend."),
        DT::DTOutput(ns("summary"), fill = FALSE)),
      bslib::nav_panel("All measurements", icon = shiny::icon("table"), DT::DTOutput(ns("flags"), fill = FALSE)),
      bslib::nav_panel("Method", icon = shiny::icon("book"), method_note(
        "<h5>How the thresholds are set</h5>",
        "<p>Only the control group is used. A linear model <code>Biomarker ~ days + ClinicalNumber</code> ",
        "is fitted to the control measurements, and its residual standard deviation is taken ",
        "as the normal visit-to-visit variation. The reference band is:</p>",
        "<p style='text-align:center'><code>control mean ± z × residual SD</code></p>",
        "<p>Every measurement from every group is then labelled <b>Significantly High</b> ",
        "(above the band, shown as an upward triangle), <b>Significantly Low</b> (below it, ",
        "a downward triangle) or <b>Normal</b> (a dot). Point colour shows the clinical visit number.</p>",
        "<p><b>Input format (long format, one row per measurement):</b> <code>SampleID</code> ",
        "(patient), <code>Class</code> (for example Control / Metastasized), <code>ClinicalNumber</code> ",
        "(visit 1, 2, 3, ...), <code>days</code> (time since the first visit), <code>Biomarker</code> ",
        "(value), and an optional <code>Label</code> for the plot.</p>"))
    )
  )
}

uveal_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    dummy <- "data/uveal_melanoma_longitudinal_dummy.csv"
    raw <- data_source_server(input, output, dummy)
    dat <- shiny::reactive({
      d <- raw()
      miss <- setdiff(UVEAL_REQUIRED, names(d))
      shiny::validate(shiny::need(!length(miss), paste("Missing column(s):", paste(miss, collapse = ", "))))
      d$Class <- as.character(d$Class); d$SampleID <- as.character(d$SampleID)
      if (!"Label" %in% names(d)) d$Label <- d$SampleID
      d
    })
    shiny::observe({
      cl <- unique(dat()$Class)
      shiny::updateSelectInput(session, "control", choices = cl,
                               selected = if ("Control" %in% cl) "Control" else cl[1])
      shiny::updateSelectizeInput(session, "groups", choices = cl, selected = cl)
      shiny::updateSelectizeInput(session, "highlight", choices = sort(unique(dat()$SampleID)))
    })

    th <- shiny::reactive({
      shiny::req(input$control %in% dat()$Class)
      tryCatch(compute_thresholds(dat(), control = input$control, z = input$z),
               error = function(e) shiny::validate(shiny::need(FALSE, conditionMessage(e))))
    })
    flagged <- shiny::reactive({
      shiny::req(input$groups)
      d <- dat(); d <- d[d$Class %in% input$groups, ]
      shiny::req(nrow(d) > 0)
      flag_points(d, th())
    })
    the_plot <- shiny::reactive(plot_longitudinal(flagged(), th(), input$labels, input$highlight))

    output$boxes <- shiny::renderUI({
      f <- flagged(); t <- th()
      pts <- tapply(f$Status == "Significantly High", f$Class, sum)
      pat <- tapply(f$Status == "Significantly High", f$SampleID, any)
      pat_class <- tapply(f$Class, f$SampleID, `[`, 1)
      hi_pat <- tapply(pat, pat_class, sum)
      n_pat <- table(pat_class)
      bslib::layout_columns(col_widths = c(4, 4, 4), fill = FALSE,
        bslib::value_box("Reference band", sprintf("%.2f to %.2f", t$lower, t$upper),
                         sprintf("mean %.2f, residual SD %.2f", t$mean, t$sd),
                         showcase = shiny::icon("ruler-vertical"), theme = "primary"),
        bslib::value_box("Patients with a high value",
                         shiny::HTML(paste0(names(hi_pat), ": ", hi_pat, " / ", n_pat[names(hi_pat)], collapse = "<br>")),
                         showcase = shiny::icon("triangle-exclamation"), theme = "secondary"),
        bslib::value_box("High measurements",
                         shiny::HTML(paste0(names(pts), ": ", pts, collapse = "<br>")),
                         showcase = shiny::icon("arrow-up"), theme = "light"))
    })
    output$plot <- shiny::renderPlot(the_plot(), res = 96)
    output$summary <- DT::renderDT(dt_table(patient_summary(flagged()), 20))
    output$flags <- DT::renderDT({
      f <- flagged()[, c("SampleID", "Class", "ClinicalNumber", "days", "Biomarker", "Status")]
      DT::formatStyle(dt_table(f, 20), "Status", fontWeight = DT::styleEqual(
        c("Significantly High", "Significantly Low"), c("bold", "bold")),
        color = DT::styleEqual(c("Significantly High", "Significantly Low"), c("#b2182b", "#2166ac")))
    })
    output$dl_plot <- shiny::downloadHandler("uveal_longitudinal_plot.png", function(file)
      ggplot2::ggsave(file, the_plot(), width = 12, height = 6, dpi = 300, bg = "white"))
    output$dl_flags <- shiny::downloadHandler("uveal_longitudinal_flagged.csv", function(file)
      write.csv(flagged(), file, row.names = FALSE))
  })
}
