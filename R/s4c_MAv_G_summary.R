#' An S4 class for \code{major_locus_G_summary} output: aunts and uncles (by marriage)
#'
#' This is an internal S4 class within \code{\link{major_locus_G_summary}}.
#'
#' @slot Description character. Brief overview of output.
#' @slot G_summary1 data.table. Counts and proportions of observed major locus genotypes for aunts and uncles (by marriage).
#' @slot G_summary2 data.table. Counts and proportions of observed major locus genotypes for aunts and uncles (by marriage), stratified by the major locus genotype of the proband.
#'
#' @import data.table
#' @name s4c_MAv_G_summary
#' @rdname s4c_MAv_G_summary
#' @aliases s4c_MAv_G_summary-class
#' @exportClass s4c_MAv_G_summary

setClass("s4c_MAv_G_summary", slots = c(Description="character", G_summary1="data.table", G_summary2="data.table"))