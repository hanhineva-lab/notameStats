#' @import BiocGenerics
#' @importFrom utils citation
#' @importFrom dplyr .data
#' @import methods
#' @importFrom notame drop_flagged drop_qcs combined_data flag merge_notame_sets
#' mark_nas "flag<-" flag_quality log_text init_log finish_log citations
#' @import SummarizedExperiment
NULL

utils::globalVariables(c('i', '.'))

# Get internal notame functions
.add_citation <- notame:::.add_citation
.check_object <- notame:::.check_object
.get_from_name <- notame:::.get_from_name
.looks_numeric <- notame:::.looks_numeric
finite_mean <- notame:::finite_mean
finite_sd <- notame:::finite_sd
finite_median <- notame:::finite_median
finite_mad <- notame:::finite_mad
finite_quantile <- notame:::finite_quantile
finite_max <- notame:::finite_max
finite_min <- notame:::finite_min