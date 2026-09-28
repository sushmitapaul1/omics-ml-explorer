## =============================================================================
## Omics ML Explorer
## Shiny interface for three analyses from the Laboratory of Systems Tumor
## Immunology (Dept. of Dermatology, University Hospital Erlangen):
##   1. Prostate cancer  - classifiers on a protease panel (BPH vs PCa)
##   2. PDAC             - ML feature selection + classification (Olink NPX)
##   3. Uveal melanoma   - longitudinal biomarker monitoring
## Run locally:  shiny::runApp()   (files in R/ are loaded automatically)
## =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
})

theme <- bs_theme(
  version = 5, primary = APP_PRIMARY, secondary = "#3d5a80",
  base_font = font_collection("Inter", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif"),
  heading_font = font_collection("Inter", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif")
)

module_card <- function(title, icon_name, what, data, target) {
  card(
    card_header(class = "d-flex align-items-center gap-2", icon(icon_name), strong(title)),
    p(what), p(class = "small text-muted mb-2", icon("database"), " ", data),
    actionButton(paste0("go_", target), "Open", class = "btn-outline-primary btn-sm mt-auto")
  )
}

# bslib >= 0.9 groups navbar colours in navbar_options(); older versions take bg/inverse directly
navbar_args <- if ("navbar_options" %in% getNamespaceExports("bslib")) {
  list(navbar_options = bslib::navbar_options(bg = "#123f44", theme = "dark"))
} else {
  list(bg = "#123f44", inverse = TRUE)
}

ui <- do.call(page_navbar, c(list(
  id = "main",
  title = span(icon("microscope"), " Omics ML Explorer"),
  theme = theme,
  fillable = FALSE,
  header = tags$head(tags$style(HTML("
    .method-note { max-width: 820px; line-height: 1.6; }
    .method-note code { color: #0b5560; }
    .hero { padding: 2rem 0 1rem; }
    .hero h1 { font-weight: 700; }
    .bslib-value-box .value-box-value { font-size: 1.15rem; }
    .demo-banner { background:#fff6e0; border-radius:8px; padding:.55rem .9rem; font-size:.9rem; }
  "))),
  nav_panel("Overview", icon = icon("house"), value = "home",
    div(class = "container-lg",
      div(class = "hero",
        h1("Omics ML Explorer"),
        p(class = "lead", "Interactive versions of three cancer-biomarker analyses: machine-learning ",
          "classifiers, robust feature selection, and longitudinal monitoring.")),
      div(class = "demo-banner mb-4", icon("flask"), " ",
          strong("Demo mode:"), " every module starts with simulated (dummy) data that has the same ",
          "structure as the real studies. No patient data is included. You can upload your own CSV in each module."),
      layout_columns(col_widths = c(4, 4, 4),
        module_card("Prostate cancer: protease classifiers", "venus-mars",
          what = paste0("Train and compare SVM, random forest, MLP, KNN, naive Bayes and GBM on a ",
                        "protease panel with and without PSA, over repeated random train/test splits."),
          data = "Dummy: 120 samples (60 BPH / 60 PCa), 10 proteases + PSA", target = "prostate"),
        module_card("PDAC: feature selection + classification", "filter",
          what = paste0("Select proteins with stability selection, Boruta, RFE-CV and mRMR, build a ",
                        "consensus panel, then test six classifiers on held-out patients."),
          data = "Dummy: 50 patients (25 vs 25), 80 Olink-style proteins", target = "pdac"),
        module_card("Uveal melanoma: longitudinal analysis", "chart-line",
          what = paste0("Follow a blood biomarker over repeated visits and flag values outside a ",
                        "reference band derived from control patients."),
          data = "Dummy: 10 Control + 10 Metastasized patients, 4 visits each", target = "uveal")
      ),
      p(class = "small text-muted mt-4",
        "Built with R Shiny, caret, pROC, glmnet, Boruta, praznik and RobustRankAggreg. ",
        "Source code: ", a("github.com/sushmitapaul1/omics-ml-explorer",
                            href = "https://github.com/sushmitapaul1/omics-ml-explorer", target = "_blank"))
    )),
  nav_panel("Prostate cancer", icon = icon("venus-mars"), value = "prostate", prostate_ui("prostate")),
  nav_panel("PDAC", icon = icon("filter"), value = "pdac", pdac_ui("pdac")),
  nav_panel("Uveal melanoma", icon = icon("chart-line"), value = "uveal", uveal_ui("uveal")),
  nav_spacer(),
  nav_item(a(icon("github"), href = "https://github.com/sushmitapaul1/omics-ml-explorer", target = "_blank"))
), navbar_args))

server <- function(input, output, session) {
  for (t in c("prostate", "pdac", "uveal")) local({
    target <- t
    observeEvent(input[[paste0("go_", target)]], nav_select("main", target))
  })
  prostate_server("prostate")
  pdac_server("pdac")
  uveal_server("uveal")
}

shinyApp(ui, server)
