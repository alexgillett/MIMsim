#' An S4 class representing output from environment_E_summary
#'
#' This is an internal S4 class within \code{\link{environment_E_summary}}.
#'
#' @slot E data.table. Summary measures for the distribution of environmental component \eqn{E} by relative type.
#' @slot E_G0 data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 0}, by relative type. That is, the distribution of \eqn{E | G = 0}, by relative type.
#' @slot E_G1 data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 1}, by relative type. That is, the distribution of \eqn{E | G = 1}, by relative type.
#' @slot E_G2 data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 2}, by relative type. That is, the distribution of \eqn{E | G = 2}, by relative type.
#' @slot E_Carrier data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals carrying a risk genotype at the major locus, by relative type. That is, the distribution of \eqn{E | G \geq 1}{E | G \ge 1}, by relative type.
#'
#' @import data.table
#' @name s4c_environment_summary
#' @rdname s4c_environment_summary
#' @aliases s4c_environment_summary-class
#' @exportClass s4c_environment_summary

setClass("s4c_environment_summary", slots=c(E = "data.table", E_G0 = "data.table", E_G1 = "data.table", E_G2 = "data.table", E_Carrier = "data.table"))