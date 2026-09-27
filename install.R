## Install every package the app needs (run once):  Rscript install.R
pkgs <- c("shiny", "bslib", "DT", "ggplot2", "ggrepel", "dplyr",
          "caret", "pROC", "randomForest", "gbm", "kernlab", "nnet", "klaR", "e1071",
          "glmnet", "Boruta", "praznik", "RobustRankAggreg")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("All packages installed. Start the app with:  shiny::runApp()\n")
