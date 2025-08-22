#' @import BiocGenerics
#' @importFrom utils citation
#' @import methods
#' @importFrom notame drop_qcs combined_data flag merge_objects mark_nas
#' "flag<-" flag_quality log_text init_log finish_log citations
#' @import SummarizedExperiment
NULL

utils::globalVariables(c('i', '.'))

# Get internal notame functions
.add_citation <- notame:::.add_citation
.get_from_name <- notame:::.get_from_name
.check_object <- notame:::.check_object
finite_mean <- notame:::finite_mean
finite_sd <- notame:::finite_sd
finite_median <- notame:::finite_median
finite_mad <- notame:::finite_mad
finite_quantile <- notame:::finite_quantile
finite_max <- notame:::finite_max
finite_min <- notame:::finite_min