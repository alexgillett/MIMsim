#' An S4 class for \code{major_locus_G_summary} output: cousins
#'
#' This is an internal S4 class within \code{\link{major_locus_G_summary}}.
#'
#' @slot Description character. Brief overview of output.
#' @slot G_summary1 data.table. Counts and proportions of observed major locus genotypes for cousins.
#' @slot G_summary2 data.table. Counts and proportions of observed major locus genotypes for cousins, stratified by the major locus genotype of the proband.
#' @slot G_summary3 data.table. Counts and proportions of observed major locus genotypes for cousins, stratified by sibship size.
#' @slot G_summary4 data.table. Counts and proportions of observed major locus genotypes for cousins, stratified by the major locus genotype of the proband and sibship size. 
#'
#' @import data.table
#' @name s4c_Cousin_G_summary
#' @rdname s4c_Cousin_G_summary
#' @aliases s4c_Cousin_G_summary-class
#' @exportClass s4c_Cousin_G_summary
setClass("s4c_Cousin_G_summary", slots = c(Description="character", G_summary1="data.table", G_summary2="data.table", G_summary3="data.table", G_summary4="data.table"))