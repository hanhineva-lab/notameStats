# ----------- Random Forest ----------------

#' Fit Random Forest
#'
#' Fits a random forest, where given response column in pheno data is predicted 
#' using the features. Can be used both for classification and regression. For 
#' more information, see the documentation of 
#' \code{\link[randomForest]{randomForest}}.
#' After fitting the random forest, use \code{\link{importance_rf}} as a 
#' shortcut for getting the feature importance in random forest prediction.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param y character, column name of pheno data giving the dependent variable 
#' of the model
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates character, column names of pheno data
#' to use as covariates in the model, in addition to molecular features
#' @param importance Should importance of features be assessed?
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters passed to 
#' \code{\link[randomForest]{randomForest}}
#'
#' @return An object of class randomForest.
#'
#' @seealso \code{\link[randomForest]{randomForest}}, 
#' \code{\link{importance_rf}}
#'
#' @examples
#' data(toy_notame_set, package = "notame")
#' rf <- fit_rf(toy_notame_set, y = "Group")
#' rf
#' importance_rf(rf)
#'
#' @export
fit_rf <- function(object, y, all_features = FALSE, 
                   covariates = NULL, importance = TRUE, 
                   assay.type = NULL, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package \"randomForest\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("randomForest package was used to fit random forest models:",
                citation("randomForest"))
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates),
                         assay.type = from)
  object <- drop_flagged(object, all_features = all_features)

  x <- combined_data(object, assay.type = from)
  x <- x[, c(rownames(object), covariates)]
  rf <- randomForest::randomForest(x = x, y = colData(object)[, y], 
                                   importance = importance, ...)

  rf
}


#' Feature importance in random forest
#'
#' Extracts feature importance in random forest in a nice format.
#'
#' @param rf An object of class randomForest
#'
#' @return A data frame of feature importance.
#'
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{fit_rf}}
#'
#' @examples
#' data(toy_notame_set, package = "notame")
#' rf <- fit_rf(toy_notame_set, y = "Group")
#' rf
#' importance_rf(rf)
#'
#' @export
importance_rf <- function(rf) {
  # Extract metrics and feature ID
  df <- data.frame(Feature_ID = rownames(rf$importance),
                   as.data.frame(rf$importance),
                   stringsAsFactors = FALSE, check.names = FALSE)
  df
}


# ------------------ mixOmics PLS ---------------------------

#' A helper function for extracting predictor matrix with covariates
#'
#' @param object a SummarizedExperiment object
#' @param covariates character, column names of pheno data to use as covariates 
#' in the model, in addition to molecular features
#' @return A data frame with predictors, including covariates.
#' @noRd
.get_x <- function(object, covariates, assay.type = NULL) {
  # Convert covariates to numeric
  if (any(!vapply(colData(object)[, covariates], .looks_numeric, logical(1)))) {
    stop("All covariates should be convertable to numeric.")
  }
  colData(object)[covariates] <- lapply(colData(object)[covariates], as.numeric)

  # Extract X
  x <- combined_data(object, assay.type)[, c(rownames(object), covariates)]
  x
}

#' PLS
#'
#' Simple wrappers for fitting a PLS model using mixOmics package. The result 
#' can then be passed to many of the mixOmics functions for prediction, 
#' performance evaluation etc. 
#' 
#' \itemize{
#' \item{\code{mixomics_pls} A simple PLS model with set number of components 
#' and all features}
#' \item{\code{mixomics_pls_optimize} Test different numbers of components}
#' \item{\code{mixomics_spls_optimize} sPLS model: Test different numbers of 
#' components and features}
#' }
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param y character vector, column names of the grouping variable to predict
#' @param ncomp number of X components
#' @param plot_perf plot performance of models in cross-validation
#' @param folds the number of folds to use in k-fold cross validation
#' @param nrepeat the number of times to repeat the cross validation. Lower 
#' this for faster testing.
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates character, column names of pheno datato use as covariates 
#' in the model, in addition to molecular features
#' @param n_features the number of features to try for each component
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... any parameters passed to 
#' \code{\link[mixOmics]{pls}} or \code{\link[mixOmics]{spls}}
#'
#' @return An object of class "mixo_pls" or "mixo_spls". For the optimized and 
#' sparse models, a list with object of class "mixo_plsda" and a performance 
#' plot.
#'
#' @examples
#' data(toy_notame_set, package = "notame")
#' pls_res <- mixomics_pls(toy_notame_set, y = "Injection_order", ncomp = 3)
#' # Cross-validation repeated only 5 times for quick run time 
#' pls_opt <- mixomics_pls_optimize(toy_notame_set, 
#'   y = "Injection_order", ncomp = 3, nrepeat = 5)
#' spls_opt <- mixomics_spls_optimize(toy_notame_set,
#'   y = "Injection_order", ncomp = 3,
#'   n_features = c(1:10, 12, 15, 20), nrepeat = 5
#' )
#' # Plot score plot of any final model
#' mixOmics::plotIndiv(pls_res,  
#'   comp = seq_len(2), group = toy_notame_set$Group, 
#'   ind.names = FALSE, title = "PLS scores plot", legend = TRUE)
#' 
#' # Proportion of variance explained
#' pls_res$prop_expl_var$X[seq_len(2)] |> round(digits = 3) * 100
#' 
#' @name pls
#' @seealso \code{\link[mixOmics]{pls}}, \code{\link[mixOmics]{perf}},
#' \code{\link[mixOmics]{spls}}, \code{\link[mixOmics]{tune.spls}}
NULL

#' @rdname pls
#' @export
mixomics_pls <- function(object, y, ncomp,
                         all_features = FALSE, covariates = NULL, 
                         assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))

  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[y]

  log_text("Fitting PLS")
  pls_model <- mixOmics::pls(predictors, outcome, ncomp = ncomp, ...)

  pls_model
}

#' @rdname pls
#'
#' @export
mixomics_pls_optimize <- function(object, y, ncomp, plot_perf = FALSE, 
                                  folds = 5, nrepeat = 50,
                                  all_features = FALSE, covariates = NULL, 
                                  assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)
  pls_res <- mixomics_pls(object = object, y = y, ncomp = ncomp, 
                          all_features = all_features,
                          covariates = covariates, assay.type = from, ...)
    
  log_text("Evaluating PLS performance")
  perf_pls <- mixOmics::perf(pls_res, validation = "Mfold", 
                             folds = folds, nrepeat = nrepeat)
  
  # Find the optimal number of components
  ncomp_opt <- which(perf_pls$measure$MSEP$summary$mean ==
    min(perf_pls$measure$MSEP$summary$mean))[1]
  log_text(paste0("Choosing a PLS model with ", ncomp_opt, 
                  " component(s) based on the minimal MSE\n",
                  "Take a look at the performance plot and make sure this", 
                  " is the correct number of components"
  ))

  pls_final <- mixomics_pls(object = object, y = y, ncomp = ncomp_opt,
                            all_features = all_features,
                            covariates = covariates, assay.type = from, ...)
  if (plot_perf) {
    p <- plot(perf_pls, measure = "MSEP")
    return(list(model = pls_final, plot_perf = p))
  }
  pls_final
}

#' @rdname pls
#'
#' @export
mixomics_spls_optimize <- function(object, y, ncomp, plot_perf = FALSE,
                                   n_features =
                                   c(seq_len(10), seq(20, 300, 10)), 
                                   folds = 5, nrepeat = 50,
                                   all_features = FALSE, covariates = NULL,
                                   assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[y]

  # Test different number of components and features with cross validation
  log_text("Tuning sPLS")
  tuned_spls <- mixOmics::tune.spls(predictors, outcome, ncomp = ncomp,
                                    test.keepX = n_features,
                                    validation = "Mfold", folds = folds,
                                    nrepeat = nrepeat, measure = "MAE")

  # Choose optimal numbers of components and features
  ncomp_opt <- tuned_spls$choice.ncomp$ncomp
  keep_x <- tuned_spls$choice.keepX[seq_len(ncomp_opt)]
  log_text(paste("Final model has", ncomp_opt, 
                 "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))
  # Fit the final model
  spls_final <- mixOmics::spls(predictors, outcome, ncomp = ncomp_opt, 
                               keepX = keep_x, ...)

  if (plot_perf) {
    # Plot error for each component with different number of features
    p < plot(tuned_spls, measure = "MSEP")
    return(list(model = spls_final, plot_perf = p))
  }
  spls_final
}

#' PLS-DA
#'
#' A simple wrapper for fitting a PLS-DA model using mixOmics package. 
#' The object can then be passed to many of the mixOmics 
#' functions for prediction, performance evaluation etc.
#' \itemize{
#' \item{\code{mixomics_plsda} A simple PLS-DA model with set number of 
#' components and all features}
#' \item{\code{mixomics_plsda_optimize} Test different numbers of components,
#' choose the one with minimal balanced error rate}
#' \item{\code{mixomics_splsda_optimize} Test different numbers 
#' of components and features, choose the one with minimal balanced error rate}
#' }
#'
#' @param object a SummarizedExperiment object
#' @param y character, column name of the grouping variable to predict
#' @param ncomp the number of X components
#' @param plot_perf plot performance of models in cross-validation
#' @param folds the number of folds to use in k-fold cross validation
#' @param nrepeat the number of times to repeat the cross validation. Lower 
#' this for faster testing.
#' @param n_features the number of features to try for each component
#' @param dist the distance metric to use, one of "max.dist", 
#' "mahalanobis.dist", "centroids.dist".
#' use \code{\link{mixomics_plsda_optimize}} to find the best distance metric
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates character, column names of pheno data to use as covariates 
#' in the model, in addition to molecular features
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... any parameters passed to 
#' \code{\link[mixOmics]{plsda}}
#'
#' @return An object of class "mixo_plsda" or for the optimized and sparse 
#' models, a list with object of class "mixo_plsda" and a performance plot.
#'
#' @examples
#' data(toy_notame_set, package = "notame")
#' noqc <- notame::drop_qcs(toy_notame_set)
#' plsda_res <- mixomics_plsda(noqc, y = "Group", ncomp = 2)
#' # Cross-validation repeated only 5 times for quick run time 
#' set.seed(38)
#' plsda_opt <- mixomics_plsda_optimize(noqc, 
#'   y = "Group", ncomp = 3, nrepeat = 5
#' )
#' set.seed(38)
#' splsda_opt <- mixomics_splsda_optimize(noqc,
#'   y = "Group", dist = "max.dist", ncomp = 2,
#'   n_features = c(1:10, 12, 15, 20), nrepeat = 5
#' )
#' # Plot PLS-DA scores
#' mixOmics::plotIndiv(plsda_res, 
#'   comp = seq_len(2), group = notame::drop_qcs(toy_notame_set)$Group, 
#'   ind.names = FALSE, title = "PLS-DA scores plot", legend = TRUE, 
#'   ellipse = TRUE)
#'
#' # Plot prediction areas
#' background <- mixOmics::background.predict(plsda_res, 
#'   comp.predicted = 2, dist = "max.dist")
#' mixOmics::plotIndiv(plsda_res,
#'   comp = seq_len(2), group = notame::drop_qcs(toy_notame_set)$Group, 
#'   ind.names = FALSE, 
#'   title = "prediction areas", legend = TRUE, background = background)
#' 
#' @name pls_da
#' @seealso \code{\link[mixOmics]{plsda}}, \code{\link[mixOmics]{perf}},
#' \code{\link[mixOmics]{splsda}}, \code{\link[mixOmics]{tune.splsda}}
NULL

#' @rdname pls_da
#' @export
mixomics_plsda <- function(object, y, ncomp, 
                           all_features = FALSE, covariates = NULL,
                           assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         "Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)  
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                          assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[, y]
  # outcome needs to be a factor, this ensures the levels are right
  if (!is(outcome, "factor")) {
    outcome <- as.factor(outcome)
    warning(y, " is not encoded as a factor, converted to factor with levels: ",
            paste(levels(outcome), collapse = ", "))
  }
  log_text("Fitting PLS-DA")
  plsda_model <- mixOmics::plsda(predictors, outcome, ncomp = ncomp, ...)

  plsda_model
}


#' @rdname pls_da
#' @export
mixomics_plsda_optimize <- function(object, y, ncomp, plot_perf = FALSE,
                                    folds = 5, nrepeat = 50,
                                    all_features = FALSE, covariates = NULL, 
                                    assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)  
  object <- .check_object(object, pheno_cols = c(y, covariates), 
                          assay.type = from)

  plsda_res <- mixomics_plsda(object = object, y = y, ncomp = ncomp,
                              all_features = all_features,
                              covariates = covariates, 
                              assay.type = from, ...)

  log_text("Evaluating PLS-DA performance")
  perf_plsda <- mixOmics::perf(plsda_res, validation = "Mfold", folds = folds,
                               auc = TRUE, nrepeat = nrepeat)

  # Find the distance metric with minimum BER
  ber <- perf_plsda$error.rate$BER
  inds <- which(ber == min(ber), arr.ind = TRUE)[1, ]
  dist_met <- colnames(ber)[inds[2]]
  # Find the optimal number of components
  ncomp_opt <- perf_plsda$choice.ncomp["BER", dist_met]
  log_text(paste("Choosing a PLS-DA model with", ncomp_opt, 
                 "components using", dist_met))

  plsda_final <- mixomics_plsda(object = object, y = y, ncomp = ncomp_opt, 
                                all_features = all_features,
                                covariates = covariates, assay.type = from, ...)

  if (plot_perf) {
    plot(perf_plsda, col = mixOmics::color.mixo(seq_len(3)), 
         sd = TRUE, legend.position = "horizontal")
    graphics::title("Performance of PLS-DA models")
    p <- recordPlot()
    return(list(model = plsda_final, plot_perf = p))
  }

  plsda_final
}


#' @rdname pls_da
#' @export
mixomics_splsda_optimize <- function(object, y, ncomp, dist, plot_perf = FALSE,
                                     n_features = c(seq_len(10), 
                                                    seq(20, 300, 10)),
                                     folds = 5, nrepeat = 50,
                                     all_features = FALSE, covariates = NULL, assay.type = NULL, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package \"mixOmics\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation("mixOmics package was used to fit PLS models:",
                citation("mixOmics"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object,  pheno_cols = c(y, covariates),
                         assay.type = from)

  predictors <- .get_x(object, covariates, from)
  outcome <- colData(object)[, y]
  # outcome needs to be a factor, this ensures the levels are right
  if (!is(outcome, "factor")) {
    outcome <- as.factor(outcome)
    warning(y, " is not encoded as a factor, converted to factor with levels: ",
            paste(levels(outcome), collapse = ", "))
  }
  # Test different components and numbers of features with cross validation
  log_text("Tuning sPLS-DA")
  tuned_splsda <- mixOmics::tune.splsda(predictors, outcome, ncomp = ncomp,
                                        validation = "Mfold", folds = folds,
                                        dist = dist, measure = "BER", 
                                        nrepeat = nrepeat,
                                        test.keepX = n_features)
  
  # Choose optimal numbers of components and features
  ncomp_opt <- tuned_splsda$choice.ncomp$ncomp
  keep_x <- tuned_splsda$choice.keepX[seq_len(ncomp_opt)]
  log_text(paste("Final model has", ncomp_opt, 
                 "components with the numbers of features:",
                 paste(keep_x, collapse = ", ")))
  # Fit the final model
  splsda_final <- mixOmics::splsda(predictors, outcome, ncomp = ncomp_opt,
                                   keepX = keep_x)
  if (plot_perf) {
    # Plot error rate of different components as a function of number of features
    p <- plot(tuned_splsda)
    return(list(model = splsda_final, plot_perf = p))
  }
  
  splsda_final
}

# ------------------- MUVR --------------------------------

#' Multivariate modelling with minimally biased variable selection (MUVR)
#'
#' A wrapper around \code{\link[MUVR2]{MUVR2}} (random forest,
#' PLS(-DA)) and \code{\link[MUVR2]{MUVR2_EN}} (elastic net)
#' functions from the MUVR2 package. 
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param y character, column name in pheno data of the target variable
#' @param id character, column name in pheno data of the subject ID variable in 
#' case of repeated measurements
#' @param multi_level logical, whether multi-level modeling should be applied, 
#' see Details
#' @param multi_level_var character, column name in pheno data of the variable 
#' for splitting the data in multi-level modeling
#' @param all_features logical, should all features be included in the model? 
#' if FALSE, flagged features are left out
#' @param covariates,static_covariates character, column names of 
#' pheno data to use as covariates in the model, in addition to molecular 
#' features. 
#' \code{static_covariates} are ignored for non-multi-level models.
#' For multi-level models, the change in \code{covariates} is computed, while
#' \code{static_covariates} are taken from the first time point. 
#' @param nRep Number of repetitions of double CV, parameter of MUVR
#' @param nOuter Number of outer CV loop segments, parameter of MUVR
#' @param nInner Number of inner CV loop segments, parameter of MUVR
#' @param varRatio Ratio of variables to include in subsequent inner loop 
#' iteration, parameter of MUVR
#' @param method Multivariate method. Supports 'PLS', 'RF' and 'EN'
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters to \code{\link[MUVR2]{MUVR2}} or 
#' \code{\link[MUVR2]{MUVR2_EN}} and 
#' \code{\link[MUVR2]{getVar}} (when method == "EN")
#'
#' @return A MUVR object.
#'
#' @details This function is now using the MUVR2 package, characterized as an 
#' upgrade extending the original MUVR package by the inclusion of elastic net 
#' regression (EN) and some functionality not covered by this wrapper. Elastic 
#' net regression supports covariate adjustment by suppressing regularization 
#' of specified features from the regularization procedure. Note that this is 
#' different from simply including covariates such as sex. EN also differs from 
#' PLS and RF in that no recursive variable elimination is performed, so an 
#' additional scheme is used to obtain the 'min', 'mid' and 'max' models using
#' \code{\link[MUVR2]{getVar}}. 
#' 
#' Sex would be entered as a static covariate, since the change in sex is zero 
#' for all individuals, so computing the change and using that as a covariate 
#' does not make sense.
#'
#' Note that there are several more plots available in MUVR2 for inspecting the 
#' results, notably \code{\link[MUVR2]{plotMV}}, 
#' \code{\link[MUVR2]{plotStability}} and \code{\link[MUVR2]{plotVIRank}}
#' Many of these return different plots depending on the model specification.
#'
#' @examples
#' data(toy_notame_set, package = "notame")
#' ex_set <- notame::drop_qcs(toy_notame_set)[1:10, ]
#' ex_set$Injection_order <- as.numeric(ex_set$Injection_order)
#' # Simple PLS regression model
#' pls_model <- muvr_analysis(ex_set, 
#'   y = "Injection_order", nRep = 2, method = "PLS")
#'
#' # RF classification with covariate and repeated measures (not longitudinal)
#' rf_model <- muvr_analysis(ex_set, y = "Group", id = "Subject_ID", 
#'   nRep = 2, method = "RF", covariates = "Injection_order")
#'
#' # RF classification on multilevel variable comparing levels of y
#' rf_model_ <- muvr_analysis(ex_set, 
#'   y = "Group", multi_level = TRUE, id = "Subject_ID", 
#'   multi_level_var = "Time", method = "RF", nRep = 2)
#'
#' # EN regression on multilevel variable with covariate and static covariate
#' ex_set$Group <- as.numeric(ex_set$Group)
#' en_model <- muvr_analysis(ex_set, id = "Subject_ID", 
#'  multi_level = TRUE, multi_level_var = "Time", 
#'  covariates = "Injection_order", static_covariates = "Group", 
#'  method = "EN", nRep = 2)
#'
#' @seealso \code{\link[MUVR2]{MUVR2}} \code{\link[MUVR2]{MUVR2_EN}} 
#' \code{\link[MUVR2]{getVar}} \code{\link[MUVR2]{plotMV}} 
#' \code{\link[MUVR2]{plotStability}} \code{\link[MUVR2]{plotVIRank}}
#' \code{\link[MUVR2]{plotVAL}}
#' @export
muvr_analysis <- function(object, y = NULL, id = NULL, multi_level = FALSE,
                          multi_level_var = NULL, covariates = NULL,
                          static_covariates = NULL, all_features = FALSE,
                          nRep = 50, nOuter = 6, nInner = nOuter - 1,
                          varRatio = 0.75, method = c("PLS", "RF"),
                          assay.type = NULL, ...) {
  if (!requireNamespace("MUVR2", quietly = TRUE)) {
    stop("Package \"MUVR2\" needed for this function to work.",
         " Please install it.", call. = FALSE)
  }
  .add_citation(paste("MUVR2 package was used to fit multivariate models",
                      "with variable selection:"),
                citation("MUVR2"))
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, 
                         pheno_cols = c(y, id, covariates, static_covariates,
                                        multi_level_var, covariates), 
                         assay.type = from)

  # MUVR2 can only use numeric predictors
  classes <- vapply(colData(object)[, c(covariates, static_covariates)], 
                    class, character(1))
  if (length(classes) && any(classes != "numeric")) {
    stop("MUVR2 can only deal with numeric inputs,", 
         " please transform all covariates to numeric", call. = FALSE)
  }

  if (any(!vapply(colData(object)[, covariates], .looks_numeric, logical(1)))) {
    stop("All covariates should be convertable to numeric.")
  }
  colData(object)[covariates] <- lapply(colData(object)[covariates], as.numeric)
  
  # Do do.call with MUVR2::MUVR2_EN if method == "EN", to avoid nesting
  if (method == "EN") {
    func <- MUVR2::MUVR2_EN
  } else {
    func <- MUVR2::MUVR2
  }
  
  # Substitute additional arguments to formal arguments of MUVR2 or MUVR2_EN
  add_args <- substitute(...())
  add_args <- add_args[names(add_args) %in% formalArgs(func)]
  
  # Classic MUVR
  if (!multi_level) {
    if (is.null(y)) {
      stop("y variable needs to be defined unless doing multi-level modeling.")
    }
    predictors <- combined_data(object, from)[, c(rownames(object), covariates)]
    outcome <- colData(object)[, y]

    # Independent samples
    if (is.null(id)) {
      # Make list of arguments for do.call, match additional arguments to 
      # formal arguments of MUVR2 or MUVR2_EN
      args <- list(X = predictors, Y = outcome, nRep = nRep, nOuter = nOuter, 
                   nInner = nInner, varRatio = varRatio, method = method)
      muvr_model <- do.call(func, c(args, add_args))
    } else {
      # Multiple measurements
      ID <- as.numeric(colData(object)[, id])
      args <- list(X = predictors, Y = outcome, ID = ID, nRep = nRep, 
                   nOuter = nOuter, nInner = nInner, varRatio = varRatio, 
                   method = method)

      muvr_model <- do.call(func, c(args, add_args))
    }
  } else { # Multi-level analysis
    if (is.null(id) || is.null(multi_level_var)) {
      stop("id and multi_level_var needed for multi-level modeling.")
    }
    # Check that multi_level_var has only 2 unique values
    ml_var <- colData(object)[, multi_level_var] <-
      as.factor(colData(object)[, multi_level_var])
    if (length(levels(ml_var)) != 2) {
      stop("The multilevel variable should have exactly 2 unique values.")
    } else {
      message("Computing effect matrix according to ", multi_level_var, " : ",
              levels(ml_var)[2], " - ", levels(ml_var)[1])
    }

    # Compute effect matrix with covariates
    cd <- combined_data(object, from)
    cd <- cd[order(cd[, id]), ]
    x1 <- cd[cd[, multi_level_var] == levels(ml_var)[1],
             c(rownames(object), covariates)]
    x2 <- cd[cd[, multi_level_var] == levels(ml_var)[2],
             c(rownames(object), covariates)]
    predictors <- x2 - x1
    # Add static covariates, where we don't want to compute change, such as sex
    predictors[, static_covariates] <- 
      cd[cd[, multi_level_var] == levels(ml_var)[1], static_covariates]
    rownames(predictors) <- unique(cd[, id])

    # Modeling
    if (!is.null(y)) { # Compare change of multi_level_var between levels of y
      outcome <- cd[cd[, multi_level_var] == levels(ml_var)[1], y]
      args <- list(X = predictors, Y = outcome, nRep = nRep, nOuter = nOuter,
                   nInner = nInner, varRatio = varRatio, method = method)
      muvr_model <- do.call(func, c(args, add_args))
    } else { # Compare levels of multi_level_var
      args <- list(X = predictors, ML = TRUE, nRep = nRep, nOuter = nOuter,
                   nInner = nInner, varRatio = varRatio, method = method)
      muvr_model <- do.call(func, c(args, add_args))
    }
  }
  
  if (method == "EN") {
    add_args <- substitute(...())
    add_args <- add_args[names(add_args) %in% formalArgs(MUVR2::getVar)]    
    muvr_model <- do.call(MUVR2::getVar, 
                          c(list(rdCVnetObject = muvr_model), add_args))
  }
  
  muvr_model
}

# PERMANOVA ----

#' PERMANOVA
#'
#' Performs permutational multivariate analysis of variance. Uses package 
#' called PERMANOVA.
#'
#' @param object a \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object
#' @param group character, name of the column to compare
#' @param all_features should all features be included?
#' @param transform Transformation to use in 
#' \code{\link[PERMANOVA]{IniTransform}}. By default uses "Standardize columns".
#' @param coef Coefficient to calculate continuous distances in 
#' \code{\link[PERMANOVA]{IniTransform}}.
#' By default uses Pythagorean distances.
#' @param assay.type character, assay to be used in case of multiple assays
#' @param ... other parameters to \code{\link[PERMANOVA]{PERMANOVA}}
#'
#' @return A PERMANOVA object.
#'
#' @examples
#' data(toy_notame_set, package = "notame")
#' permanova_res <- perform_permanova(
#'   notame::drop_qcs(toy_notame_set), 
#'   group = "Group")
#'
#' @export
perform_permanova <- function(object, group, all_features = FALSE,
                              transform = "Standardize columns",
                              coef = "Pythagorean", assay.type = NULL, ...) {
  if (!requireNamespace("PERMANOVA", quietly = TRUE)) {
    stop("Package \"PERMANOVA\" needed for this function to work.", 
         " Please install it.", call. = FALSE)
  }
  
  .add_citation(paste("PERMANOVA was used for permutational multivariate", 
                      "analysis of variance:"),
                citation("PERMANOVA"))
                
  object <- drop_flagged(object, all_features = all_features)
  from <- .get_from_name(object, assay.type)
  object <- .check_object(object, pheno_factors = group,
                         assay.type = from)

  log_text("Starting PERMANOVA tests")
  data <- t(assay(object, from))
  data <- PERMANOVA::IniTransform(data, transform = transform)
  initialized <- PERMANOVA::DistContinuous(data, coef = coef)
  res <- PERMANOVA::PERMANOVA(initialized, colData(object)[, group], ...)
  log_text("PERMANOVA performed")

  res
}