#' @title Mixed inheritance model simulations of disease within family studies
#'
#' @description A package for simulating data from family studies of disease. A mixed inheritance model (MIM) is used to model disease. This is a liability threshold approach where the risk of disease is jointly influenced by a polygenic component (\eqn{A}), a major disease locus (\eqn{G}) and an environmental component (\eqn{E}). Simulated pedigrees are bi-linear and span 3 or 4 generations. The ascertainment of families can be influenced by the family history of disease, allowing users of \code{MIMsim} to explore the effect of ascertainment bias within family studies of disease risk and penetrance.
#'
#' @details \code{\link{mimsim}} is the simulation function. Users specify inputs for the disease model, family structure and ascertainment procedure. 
#' A range of summary functions exist within \code{MIMsim} to aid exploration of the outputted families. These are: \itemize{
#' \item \code{\link{polygenic_A_summary}}; summarises the distribution of the polygenic component \eqn{A} for relatives of the proband, 
#' \item \code{\link{major_locus_G_summary}}; summarises the distribution of the major disease locus \eqn{G} for a user-specified relative of the proband, stratified by the major locus genotype of the proband and, if applicable, sibship size,
#' \item \code{\link{environment_E_summary}}; summarises the distribution of the environmental component \eqn{E} for relatives of the proband, and,
#' \item \code{\link{simulation_summary}}; summarises the disease risk for relatives of the proband, stratified by major locus genotype.
#' }
#' A vignette with sample code is provided for calculating empirical 95\% confidence intervals for summary statistics outputted by the summary functions of \code{MIMsim}.
#'
#' @seealso The mixed inheritance model is a liability threshold approach where disease can be jointly influenced by a combination of: 1. a polygenic component, 2. a major disease locus and 3. an environmental component. For a brief description of the mixed inheritance model see \code{\link{mimsim}}. For a more details see 'References', in particular; Morton and MacLean (1974) for the mixed inheritance model, and Falconer (1965) for the liability threshold model.
#'
#' Output from summary functions within \code{MIMsim} come in the form of a \code{data.table}; please see \code{\link[data.table]{data.table}} for details.
#' 
#' @references 
#' Morton, N.E. and MacLean C.J. (1974) Analysis of family resemblance. 3. Complex segregation of quantitative traits. \emph{American Journal of Human Genetics (26)}, 489 - 503.
#'
#' Falconer, D.S. (1960) \emph{Introduction to Quantitative Genetics.} New York, Ronald Press Co.
#'  
#' Falconer, D.S. (1965) The inheritance of liability to certain diseases, estimated from the incidence among relatives. \emph{Annuals of Human Genetics (29)}, 51 - 76.
#'  
#' @docType package
#' @name MIMsim-package
#' 
NULL