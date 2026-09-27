## =============================================================================
## Robust ML feature selection (from R_Code_ML_FeatureSelection.R)
##   1. Stability selection (elastic net, bootstrap/subsample resampling)
##   2. Boruta (all-relevant random-forest wrapper)
##   3. RFE with repeated CV (random-forest based, caret)
##   4. mRMR (praznik)
##   5. Consensus: votes across methods + Robust Rank Aggregation
## =============================================================================

FS_METHODS <- c("Stability selection", "Boruta", "RFE-CV", "mRMR")

run_feature_selection <- function(X, y, methods = FS_METHODS,
                                  n_boot = 100, alpha = 0.7, stab_threshold = 0.6,
                                  boruta_max_runs = 100,
                                  rfe_folds = 5, rfe_repeats = 2,
                                  k_mrmr = 20, seed = 42, progress = NULL) {
  X <- as.matrix(X); y <- factor(y)
  feats <- colnames(X)
  res <- list(selected = list(), ranked = list())
  tick <- function(v, d) if (!is.null(progress)) progress(v, d)
  n_m <- length(methods); k <- 0

  if ("Stability selection" %in% methods) {
    k <- k + 1; set.seed(seed)
    freq <- setNames(numeric(ncol(X)), feats)
    for (b in seq_len(n_boot)) {
      if (b %% 10 == 1) tick((k - 1 + b / n_boot) / n_m,
                             sprintf("Stability selection: resample %d/%d", b, n_boot))
      idx <- sample(nrow(X), floor(0.8 * nrow(X)))
      cv <- tryCatch(glmnet::cv.glmnet(X[idx, ], y[idx], family = "binomial",
                                       alpha = alpha, nfolds = 5),
                     error = function(e) NULL)
      if (is.null(cv)) next
      cf <- as.matrix(coef(cv, s = "lambda.1se"))
      sel <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
      freq[sel] <- freq[sel] + 1
    }
    freq <- freq / n_boot
    res$stability <- data.frame(protein = feats, stability_freq = freq)
    res$stability <- res$stability[order(-res$stability$stability_freq), ]
    res$selected[["Stability selection"]] <- res$stability$protein[res$stability$stability_freq >= stab_threshold]
    res$ranked[["Stability selection"]] <- res$stability$protein
  }

  if ("Boruta" %in% methods) {
    k <- k + 1; tick((k - 0.5) / n_m, "Boruta (random-forest shadow features)")
    set.seed(seed)
    bf <- Boruta::Boruta(x = as.data.frame(X), y = y, doTrace = 0, maxRuns = boruta_max_runs)
    bf <- Boruta::TentativeRoughFix(bf)
    imp <- bf$ImpHistory[, names(bf$finalDecision), drop = FALSE]
    imp[!is.finite(imp)] <- NA
    res$boruta <- data.frame(protein = names(bf$finalDecision),
                             boruta_decision = as.character(bf$finalDecision),
                             boruta_meanImp = colMeans(imp, na.rm = TRUE))
    res$boruta <- res$boruta[order(-res$boruta$boruta_meanImp), ]
    res$selected[["Boruta"]] <- res$boruta$protein[res$boruta$boruta_decision == "Confirmed"]
    res$ranked[["Boruta"]] <- res$boruta$protein
  }

  if ("RFE-CV" %in% methods) {
    k <- k + 1; tick((k - 0.5) / n_m, "Recursive feature elimination with repeated CV")
    set.seed(seed)
    sizes <- unique(c(2, 5, 10, 15, 20, 30, 50))
    sizes <- sizes[sizes < ncol(X)]
    ctrl <- caret::rfeControl(functions = caret::rfFuncs, method = "repeatedcv",
                              number = rfe_folds, repeats = rfe_repeats, verbose = FALSE)
    rf <- caret::rfe(x = as.data.frame(X), y = y, sizes = sizes, rfeControl = ctrl)
    sel <- caret::predictors(rf)
    res$rfe_profile <- rf$results
    res$rfe_optsize <- rf$optsize
    res$selected[["RFE-CV"]] <- sel
    res$ranked[["RFE-CV"]] <- c(sel, setdiff(feats, sel))
  }

  if ("mRMR" %in% methods) {
    k <- k + 1; tick((k - 0.5) / n_m, "mRMR (minimum redundancy, maximum relevance)")
    mr <- praznik::MRMR(as.data.frame(X), y, k = min(k_mrmr, ncol(X)))
    res$mrmr <- data.frame(protein = names(mr$score), mrmr_score = unname(mr$score))
    res$selected[["mRMR"]] <- res$mrmr$protein
    res$ranked[["mRMR"]] <- c(res$mrmr$protein, setdiff(feats, res$mrmr$protein))
  }

  ## ---- Consensus -------------------------------------------------------------
  tick(1, "Building consensus panel")
  rra <- RobustRankAggreg::aggregateRanks(glist = res$ranked)
  colnames(rra) <- c("protein", "rra_score")
  all_sel <- unique(unlist(res$selected))
  votes <- vapply(feats, function(p) sum(vapply(res$selected, function(s) p %in% s, logical(1))), numeric(1))
  cons <- data.frame(protein = feats, n_methods_selected = votes)
  for (m in names(res$selected)) cons[[paste0("in_", make.names(m))]] <- cons$protein %in% res$selected[[m]]
  cons <- merge(cons, rra, by = "protein", all.x = TRUE)
  if (!is.null(res$stability)) cons <- merge(cons, res$stability, by = "protein", all.x = TRUE)
  if (!is.null(res$boruta)) cons <- merge(cons, res$boruta[, c("protein", "boruta_meanImp")], by = "protein", all.x = TRUE)
  cons <- cons[order(-cons$n_methods_selected, cons$rra_score), ]
  rownames(cons) <- NULL
  res$consensus <- cons
  res$methods <- names(res$selected)
  res
}
