# Omics ML Explorer

An R Shiny interface for three cancer-biomarker analyses:

| Module | Question | Methods |
|---|---|---|
| **Prostate cancer** | Can a blood protease panel separate prostate cancer (PCa) from benign prostatic hyperplasia (BPH), with or without PSA? | SVM, random forest, MLP, KNN, naive Bayes and GBM (`caret`), tuned by repeated k-fold CV over repeated random 80/20 splits |
| **PDAC** | Which proteins separate *high CTSG, short survival* from *low CTSG, long survival* patients, and how well do they classify? | Stability selection (elastic net), Boruta, RFE-CV, mRMR, consensus with Robust Rank Aggregation, then six classifiers on held-out patients |
| **Uveal melanoma** | Does a biomarker leave the normal range during follow-up in metastasized patients? | Control-derived reference band (mean ± z × residual SD of `lm(Biomarker ~ days + ClinicalNumber)`), per-visit flags |

**Project page:** https://sushmitapaul1.github.io/omics-ml-explorer/

> **Demo data only.** Every file in `data/` is simulated. It has the same structure as the real
> studies but contains no patient information. You can upload your own CSV in each module.

## Run the app

You need R version 4.1 or later.

```r
# Option 1: run directly from GitHub
install.packages("shiny")
source("https://raw.githubusercontent.com/sushmitapaul1/omics-ml-explorer/main/install.R")
shiny::runGitHub("omics-ml-explorer", "sushmitapaul1")

# Option 2: from a local copy
# git clone https://github.com/sushmitapaul1/omics-ml-explorer.git
# then, inside the folder:
source("install.R")
shiny::runApp()
```

## Input formats

**Prostate** (wide format): one row per sample, numeric feature columns (for example `Cath.G`,
`ADAM10`, …, `PSA`), and a two-level class column (`Class`: `BPH` / `PCa`).

**PDAC** (wide format): one row per patient, a two-level group column (`Group`), an ID column, and
one numeric column per protein (for example Olink NPX). Exclude the ID column and any marker that
defines the groups (`CTSG`) in the sidebar.

**Uveal melanoma** (long format): one row per measurement, with the columns `SampleID`, `Class`
(`Control` / `Metastasized`), `ClinicalNumber` (visit number), `days`, `Biomarker`, and an
optional `Label`.

Download the matching dummy CSV from each module to see the exact layout.

## Repository layout

```
app.R                 Shiny UI + server (modules are loaded from R/ automatically)
R/classifiers.R       shared caret training / evaluation engine
R/feature_selection.R stability selection, Boruta, RFE-CV, mRMR, consensus
R/longitudinal.R      thresholds, flags, patient summary, trajectory plot
R/mod_*.R             one Shiny module per analysis
R/dummy_data.R        generators for the simulated datasets (Rscript R/dummy_data.R)
data/                 simulated CSVs
docs/                 GitHub Pages project website
install.R             installs all required packages
```

## Notes on the implementation

- ROC curves and AUCs are computed with `pROC`. The original scripts used `MLeval`.
- The MLP is `caret` method `nnet` (one hidden layer). The original prostate script used `mlp` from `RSNNS`.
- By default, PDAC feature selection uses the **training split only**, so the held-out test AUC is
  not inflated by selection bias. Choose *All samples* in the sidebar to reproduce the original
  workflow.
- Default run settings are lighter than in the original analyses (for example 3 splits × 2 CV
  repeats instead of 10 × 10) so each run takes seconds. You can raise them in the *Validation
  settings* panel.

## Hosting the live app

GitHub Pages only serves static files, so it hosts the project website (`docs/`), not the R
server. To give the app a public link, deploy this folder to [shinyapps.io](https://www.shinyapps.io)
(free tier) or a lab Shiny Server:

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
rsconnect::deployApp(appName = "omics-ml-explorer")
```

## Author

Dr. Sushmita Paul, Laboratory of Systems Tumor Immunology, Department of Dermatology,
University Hospital Erlangen (FAU Erlangen-Nürnberg).
