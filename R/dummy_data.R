## =============================================================================
## Dummy (simulated) datasets used by the app.
## All values are synthetic. They mimic the structure of the real datasets
## (column names, group sizes, visit layout) but contain no patient data.
## Re-create the CSVs in data/ with:  Rscript R/dummy_data.R
## =============================================================================

## ---- 1. Prostate cancer: protease panel + PSA, BPH vs PCa ---------------------
make_prostate_dummy <- function(n_per_group = 60, seed = 123) {
  set.seed(seed)
  proteases <- c("Cath.G", "ADAM10", "ADAM17", "Aggrecanase.1", "MMP8",
                 "MMP10", "Neprilysin.1", "Neprilysin.2", "Renin", "Hepsin")
  # Group shift (in SD units) for PCa relative to BPH; most proteases carry a
  # modest signal so that multi-marker models can beat PSA alone.
  shift <- c(Cath.G = 0.5, ADAM10 = 0.3, ADAM17 = 0.6, Aggrecanase.1 = 0.1,
             MMP8 = 0.4, MMP10 = 0.8, Neprilysin.1 = -0.4, Neprilysin.2 = -0.2,
             Renin = 0.0, Hepsin = 0.9)
  n <- 2 * n_per_group
  cls <- rep(c("BPH", "PCa"), each = n_per_group)
  # Shared latent factor gives realistic correlation between proteases
  latent <- rnorm(n)
  X <- sapply(proteases, function(p) {
    base <- 0.5 * latent + rnorm(n, sd = 0.85)
    round(5 + base + ifelse(cls == "PCa", shift[[p]], 0), 3)
  })
  # PSA: log-normal, overlapping groups (the "grey zone" problem)
  psa <- round(exp(rnorm(n, mean = ifelse(cls == "PCa", 2.0, 1.6), sd = 0.45)), 2)
  df <- data.frame(SampleID = sprintf("PR%03d", seq_len(n)), X, PSA = psa,
                   Class = cls, check.names = FALSE)
  df[sample(n), ]
}

## ---- 2. PDAC: Olink-style NPX proteomics, high vs low CTSG ---------------------
make_pdac_dummy <- function(n_per_group = 25, n_proteins = 80, seed = 42) {
  set.seed(seed)
  informative <- c("ADA", "VSNL1", "MME", "CNP", "CALML4", "IL6", "CXCL13", "LAMP3")
  effect <- c(ADA = 1.3, VSNL1 = 1.1, MME = -1.2, CNP = 1.0, CALML4 = -1.0,
              IL6 = 0.6, CXCL13 = 0.5, LAMP3 = 0.4)
  filler_pool <- c("CD40", "CCL19", "CXCL9", "CXCL10", "GZMA", "GZMB", "PDCD1",
    "CD274", "TNFRSF9", "LAG3", "IL10", "IL18", "IFNG", "TNF", "CCL2", "CCL3",
    "CCL4", "CSF1", "VEGFA", "ANGPT2", "MMP7", "MMP12", "TGFB1", "HGF", "EGF",
    "FGF2", "PGF", "KDR", "TEK", "CD70", "CD27", "ICOSLG", "NCR1", "KLRD1",
    "CD8A", "CD4", "CD5", "CD83", "ARG1", "NOS3", "HMOX1", "GAL1", "LGALS9",
    "MUC16", "CEACAM5", "MSLN", "EPCAM", "KRT19", "SPP1", "TFF1", "TFF3",
    "REG4", "OLFM4", "AGR2", "CLDN18", "S100A4", "S100P", "LCN2", "PIGR",
    "CA9", "PLAU", "PLAUR", "SERPINE1", "THBS2", "COL1A1", "POSTN", "FAP",
    "ACTA2", "TIMP1", "TIMP3", "ICAM1", "VCAM1", "SELE", "ADM", "GDF15")
  fillers <- head(filler_pool, n_proteins - length(informative))
  proteins <- c(informative, fillers)
  n <- 2 * n_per_group
  grp <- rep(c("high CTSG, short survival", "Low CTSG, long survival"), each = n_per_group)
  is_high <- grp == "high CTSG, short survival"
  latent <- matrix(rnorm(n * 3), n, 3)             # 3 correlated protein modules
  X <- sapply(seq_along(proteins), function(j) {
    p <- proteins[j]
    base <- 0.6 * latent[, (j %% 3) + 1] + rnorm(n, sd = 0.8)
    eff <- if (p %in% names(effect)) effect[[p]] else 0
    round(4 + base + ifelse(is_high, eff / 2, -eff / 2), 3)
  })
  colnames(X) <- proteins
  ctsg <- round(4 + ifelse(is_high, 1.2, -1.2) + rnorm(n, sd = 0.5), 3)
  df <- data.frame(`Patient ID` = sprintf("PDAC-%02d", seq_len(n)), Group = grp,
                   CTSG = ctsg, X, check.names = FALSE)
  df[sample(n), ]
}

## ---- 3. Uveal melanoma: longitudinal biomarker, 10 Control + 10 Metastasized ---
make_uveal_dummy <- function(n_per_group = 10, n_visits = 4, seed = 7) {
  set.seed(seed)
  rows <- list()
  for (cl in c("Control", "Metastasized")) {
    for (i in seq_len(n_per_group)) {
      id <- sprintf("%s%02d", ifelse(cl == "Control", "C", "M"), i)
      days <- cumsum(c(0, round(runif(n_visits - 1, 60, 120))))
      patient_level <- rnorm(1, 0, 0.35)
      if (cl == "Control") {
        biom <- 3 + patient_level + rnorm(n_visits, 0, 0.3)
      } else {
        # Metastasized patients: rising trajectory in most, a few stay flat/low
        slope <- sample(c(rnorm(1, 0.006, 0.002), rnorm(1, -0.003, 0.001), 0),
                        1, prob = c(0.7, 0.15, 0.15))
        biom <- 3.2 + patient_level + slope * days + rnorm(n_visits, 0, 0.3)
      }
      rows[[length(rows) + 1]] <- data.frame(
        SampleID = id, Label = id, Class = cl, ClinicalNumber = seq_len(n_visits),
        days = days, Biomarker = round(biom, 3))
    }
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0) {
  dir.create("data", showWarnings = FALSE)
  write.csv(make_prostate_dummy(), "data/prostate_protease_dummy.csv", row.names = FALSE)
  write.csv(make_pdac_dummy(), "data/pdac_olink_dummy.csv", row.names = FALSE)
  write.csv(make_uveal_dummy(), "data/uveal_melanoma_longitudinal_dummy.csv", row.names = FALSE)
  cat("Dummy datasets written to data/\n")
}
