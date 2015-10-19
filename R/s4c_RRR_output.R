#' An S4 class to represent disease risk within ascertained families
#'
#' This is an internal S4 class used within \code{\link{simulation_summary}}
#'
#' @slot RRR data.table. For all family members, estimated probability of disease (column \code{K}) and recurrence risk ratio (column \code{lambda}) by relative type.
#' @slot RRR_G0 data.table. For family members with major locus genotype \eqn{G = 0}, estimated probability of disease (column \code{K_G0}) and recurrence risk ratio (column \code{lambda_G0}) by relative type.
#' @slot RRR_G1 data.table. For family members with major locus genotype \eqn{G = 1}, estimated probability of disease (column \code{K_G1}) and recurrence risk ratio (column \code{lambda_G1}) by relative type.
#' @slot RRR_G2 data.table. For family members with major locus genotype \eqn{G = 2}, estimated probability of disease (column \code{K_G2}) and recurrence risk ratio (column \code{lambda_G2}) by relative type.
#' @slot RRR_carrier data.table. For family members with major locus genotype \eqn{G \geq 1}{G \ge 1}, estimated probability of disease (column \code{K_carrier}) and recurrence risk ratio (column \code{lambda_carrier}) by relative type.
#'
#' @import data.table
#' @name s4c_RRR_output
#' @rdname s4c_RRR_output
#' @aliases s4c_RRR_output-class
#' @exportClass s4c_RRR_output
setClass(Class="s4c_RRR_output", slots = c(RRR = "data.table", RRR_G0 = "data.table", RRR_G1 = "data.table", RRR_G2 = "data.table", RRR_carrier = "data.table"))