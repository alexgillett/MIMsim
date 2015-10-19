#' An S4 class to represent disease risk and count summaries within ascertained families
#'
#' This is an internal S4 class used within \code{\link{simulation_summary}}
#'
#' @slot RRR data.table. For all family members, estimated probability of disease (column \code{K}) and recurrence risk ratio (column \code{lambda}) by relative type.
#' @slot RRR_G0 data.table. For family members with major locus genotype \eqn{G = 0}, estimated probability of disease (column \code{K_G0}) and recurrence risk ratio (column \code{lambda_G0}) by relative type.
#' @slot RRR_G1 data.table. For family members with major locus genotype \eqn{G = 1}, estimated probability of disease (column \code{K_G1}) and recurrence risk ratio (column \code{lambda_G1}) by relative type.
#' @slot RRR_G2 data.table. For family members with major locus genotype \eqn{G = 2}, estimated probability of disease (column \code{K_G2}) and recurrence risk ratio (column \code{lambda_G2}) by relative type.
#' @slot RRR_carrier data.table. For family members with major locus genotype \eqn{G \geq 1}{G \ge 1}, estimated probability of disease (column \code{K_carrier}) and recurrence risk ratio (column \code{lambda_carrier}) by relative type.
#' @slot Count_summary data.table. For each relative type, \code{Count_summary} provides a count of the number of families containing \eqn{\geq 1}{\ge 1} relative of type \eqn{R} (column \code{N_peds}), and the number of families containing \eqn{\geq 1}{\ge 1} affected relative of type \eqn{R} (column \code{N_peds_ge1_D}).
#'
#' @import data.table
#' @name s4c_simulation_summary
#' @rdname s4c_simulation_summary
#' @aliases s4c_simulation_summary-class
#' @exportClass s4c_simulation_summary

setClass("s4c_simulation_summary",  slots=c(RRR = "data.table", RRR_G0 = "data.table", RRR_G1 = "data.table", RRR_G2 = "data.table", RRR_carrier = "data.table", Count_summary = "data.table"))