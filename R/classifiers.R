## =============================================================================
## Shared classifier engine (used by the Prostate and PDAC modules)
## Based on RCode_CaretPackage_Main_Final.R and R_Code_For_Classifier.R:
##   caret::train with repeated k-fold CV, metric = ROC, evaluation on a
##   held-out test set. MLeval is replaced by pROC for ROC curves / AUC.
## =============================================================================

# Label -> caret method + extra arguments.
# Note: the original prostate script used caret "mlp" (RSNNS package); here the
# MLP is caret "nnet" (single hidden layer), as in the PDAC classifier script.
CLASSIFIERS <- list(
  "SVM (linear, class-weighted)" = list(method = "svmLinearWeights", scale = TRUE),
  "SVM (radial)"                 = list(method = "svmRadial",        scale = TRUE),
  "Random forest"                = list(method = "rf",               scale = FALSE),
  "MLP (nnet)"                   = list(method = "nnet",             scale = TRUE,
                                        extra = list(trace = FALSE, maxit = 300)),
  "KNN"                          = list(method = "knn",              scale = TRUE),
  "Naive Bayes"                  = list(method = "nb",               scale = FALSE),
  "GBM"                          = list(method = "gbm",              scale = FALSE,
                                        extra = list(verbose = FALSE),
                                        # small-n friendly grid (from R_Code_For_Classifier.R)
                                        grid = expand.grid(n.trees = c(50, 100, 150),
                                                           interaction.depth = 1:3,
                                                           shrinkage = 0.1,
                                                           n.minobsinnode = 5))
)

# Make class labels safe for caret (classProbs needs syntactic names)
prep_class <- function(y) {
  y <- factor(y)
  lv <- levels(y)
  levels(y) <- make.names(lv)
  attr(y, "original_levels") <- setNames(lv, levels(y))
  y
}

split_indices <- function(y, p = 0.8, stratified = TRUE) {
  if (stratified) as.integer(caret::createDataPartition(y, p = p, list = FALSE))
  else sample(seq_along(y), size = floor(p * length(y)))
}

fit_classifier <- function(label, x, y, ctrl, tune_length = 3) {
  spec <- CLASSIFIERS[[label]]
  args <- c(list(x = x, y = y, method = spec$method, trControl = ctrl, metric = "ROC"),
            if (is.null(spec$grid)) list(tuneLength = tune_length) else list(tuneGrid = spec$grid),
            if (spec$scale) list(preProcess = c("center", "scale")),
            spec$extra)
  suppressWarnings(do.call(caret::train, args))
}

roc_obj <- function(truth, prob_pos, positive) {
  negative <- setdiff(levels(truth), positive)
  tryCatch(pROC::roc(response = truth, predictor = prob_pos,
                     levels = c(negative, positive), direction = "<", quiet = TRUE),
           error = function(e) NULL)
}

roc_to_df <- function(r, ...) {
  if (is.null(r)) return(NULL)
  data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities, ...)
}

## Train several classifiers over one or more random train/test splits.
##   data      : data.frame with feature columns + class column
##   features  : character vector of feature column names
##   positive  : the (original) name of the positive class
##   progress  : optional function(value, detail) for Shiny progress bars
run_classifiers <- function(data, features, class_col, positive, methods,
                            n_splits = 3, p_train = 0.8, stratified = TRUE,
                            cv_folds = 5, cv_repeats = 2, tune_length = 3,
                            seed = 123, progress = NULL, train_idx = NULL) {
  # train_idx: optional fixed training rows (used by the PDAC module so that
  # feature selection and classification share the same train/test split)
  if (!is.null(train_idx)) n_splits <- 1
  y_all <- prep_class(data[[class_col]])
  pos <- make.names(positive)
  stopifnot(pos %in% levels(y_all), length(levels(y_all)) == 2)
  X_all <- data[, features, drop = FALSE]

  ctrl <- caret::trainControl(method = "repeatedcv", number = cv_folds,
                              repeats = cv_repeats, classProbs = TRUE,
                              summaryFunction = caret::twoClassSummary,
                              savePredictions = "final")
  metrics <- list(); rocs <- list(); resamp <- list(); models <- list()
  total <- n_splits * length(methods); step <- 0
  set.seed(seed)
  for (s in seq_len(n_splits)) {
    idx <- if (!is.null(train_idx)) train_idx else split_indices(y_all, p_train, stratified)
    x_tr <- X_all[idx, , drop = FALSE]; y_tr <- y_all[idx]
    x_te <- X_all[-idx, , drop = FALSE]; y_te <- y_all[-idx]
    for (m in methods) {
      step <- step + 1
      if (!is.null(progress)) progress(step / total, sprintf("Split %d/%d: %s", s, n_splits, m))
      set.seed(seed + s)
      fit <- tryCatch(fit_classifier(m, x_tr, y_tr, ctrl, tune_length),
                      error = function(e) { message(m, " failed: ", conditionMessage(e)); NULL })
      if (is.null(fit)) next

      # Cross-validated AUC of the best tuning (training set)
      cv_auc <- max(fit$results$ROC, na.rm = TRUE)
      cv_pred <- fit$pred
      cv_roc <- roc_obj(cv_pred$obs, cv_pred[[pos]], pos)
      # Held-out test set
      prob <- suppressWarnings(predict(fit, newdata = x_te, type = "prob"))[[pos]]
      cls  <- suppressWarnings(predict(fit, newdata = x_te))
      te_roc <- roc_obj(y_te, prob, pos)
      cm <- caret::confusionMatrix(cls, y_te, positive = pos)

      metrics[[length(metrics) + 1]] <- data.frame(
        split = s, classifier = m,
        cv_auc = cv_auc,
        cv_pooled_auc = if (is.null(cv_roc)) NA else as.numeric(pROC::auc(cv_roc)),
        test_auc = if (is.null(te_roc)) NA else as.numeric(pROC::auc(te_roc)),
        accuracy = unname(cm$overall["Accuracy"]),
        kappa = unname(cm$overall["Kappa"]),
        sensitivity = unname(cm$byClass["Sensitivity"]),
        specificity = unname(cm$byClass["Specificity"]),
        f1 = unname(cm$byClass["F1"]),
        n_train = length(idx), n_test = length(y_te),
        best_tune = paste(names(fit$bestTune), vapply(fit$bestTune, function(v) if (is.numeric(v)) format(signif(v, 3)) else as.character(v), character(1)), sep = "=", collapse = "; "),
        stringsAsFactors = FALSE)
      rocs[[length(rocs) + 1]] <- rbind(
        roc_to_df(cv_roc, split = s, classifier = m, set = "Training (cross-validation)"),
        roc_to_df(te_roc, split = s, classifier = m, set = "Test (held-out)"))
      rs <- fit$resample
      resamp[[length(resamp) + 1]] <- data.frame(split = s, classifier = m, ROC = rs$ROC,
                                                 Sens = rs$Sens, Spec = rs$Spec)
      if (s == n_splits) models[[m]] <- fit
    }
  }
  list(metrics = do.call(rbind, metrics), roc = do.call(rbind, rocs),
       resamples = do.call(rbind, resamp), models = models, positive = positive,
       features = features)
}

summarise_metrics <- function(metrics) {
  if (is.null(metrics) || !nrow(metrics)) return(NULL)
  agg <- function(v) sprintf("%.3f ± %.3f",mean(v, na.rm = TRUE),
                             ifelse(sum(!is.na(v)) > 1, sd(v, na.rm = TRUE), 0))
  out <- do.call(rbind, lapply(split(metrics, metrics$classifier), function(d) data.frame(
    Classifier = d$classifier[1], Splits = nrow(d),
    `CV AUC (train)` = agg(d$cv_auc), `Test AUC` = agg(d$test_auc),
    Accuracy = agg(d$accuracy), Sensitivity = agg(d$sensitivity),
    Specificity = agg(d$specificity), mean_test = mean(d$test_auc, na.rm = TRUE),
    check.names = FALSE)))
  out <- out[order(-out$mean_test), ]
  out$mean_test <- NULL
  rownames(out) <- NULL
  out
}
