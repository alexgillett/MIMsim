#' An S4 class representing output from function \code{mimsim}
#'
#' This is an internal S4 class within \code{\link{mimsim}}, used to store \code{mimsim} output. Inherits slots from S4 class \code{input_normal}.
#'
#' @slot pedigree data.frame. Contains ascertained families.
#' @slot proband_info data.frame. Contains additional information for the proband.
#' @slot pedigree_unascertained data.frame. Contains eligible but unascertained probands.
#' @slot proband_info_unascertained data.frame. Contains additional information for the unascertained but eligible probands.
#'
#' @include s4c_mimsim_input.R
#' @name s4c_mimsim_output
#' @rdname s4c_mimsim_output
#' @aliases s4c_mimsim_output-class
#' @exportClass s4c_mimsim_output
setClass("s4c_mimsim_output", slots = c(pedigree="data.frame", proband_info = "data.frame", pedigree_unascertained = "data.frame", proband_info_unascertained = "data.frame"), contains="s4c_mimsim_input")