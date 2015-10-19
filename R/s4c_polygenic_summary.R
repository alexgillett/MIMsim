#' An S4 class for \code{polygenic_A_summary} output
#'
#' This is an internal S4 class within \code{\link{polygenic_A_summary}}.
#'
#' @slot A data.table. Summary measures for the distribution of polygenic component \eqn{A} by relative type. 
#' @slot A_G0 data.table. Summary measures for the distribution of polygenic component \eqn{A} for individuals with \eqn{G = 0}, by relative type. That is, the distribution of \eqn{A | G = 0}, by relative type.
#' @slot A_G1 data.table. Summary measures for the distribution of polygenic component \eqn{A} for individuals with \eqn{G = 1}, by relative type. That is, the distribution of \eqn{A | G = 1}, by relative type.
#' @slot A_G2 data.table. Summary measures for the distribution of polygenic component \eqn{A} for individuals with \eqn{G = 2}, by relative type. That is, the distribution of \eqn{A | G = 2}, by relative type. 
#' @slot A_Carrier data.table. Summary measures for the distribution of polygenic component \eqn{A} for individuals carrying a risk genotype at the major locus, by relative type. That is, the distribution of \eqn{A | G \geq 1}{A | G \ge 1}, by relative type.
#'
#' @import data.table
#' @name s4c_polygenic_summary
#' @rdname s4c_polygenic_summary
#' @aliases s4c_polygenic_summary-class
#' @exportClass s4c_polygenic_summary

setClass("s4c_polygenic_summary", slots=c(A = "data.table", A_G0 = "data.table", A_G1 = "data.table", A_G2 = "data.table", A_Carrier = "data.table"))