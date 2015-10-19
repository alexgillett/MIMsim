#' An S4 class representing disease model parameters and the distribution of liability to disease
#'
#' This is an internal S4 class within \code{\link{mimsim}}.
#'
#' @slot model_parameters matrix. Contains the following disease model parameters: \code{D_threshold} (\eqn{t}), \code{beta_0} (\eqn{beta_0}{beta(0)}), \code{beta_1} (\eqn{beta_1}{beta(1)}) and \code{beta_2} (\eqn{beta_2}{beta(2)}).
#' @slot marginal_distn matrix. Contains the \code{Mean} and \code{Variance} of liability to disease \eqn{L}.
#'
#' @name s4c_marginal_output
#' @rdname s4c_marginal_output
#' @aliases s4c_marginal_output-class
#' @exportClass s4c_marginal_output
setClass("s4c_marginal_output", slots = c(model_parameters = "matrix", marginal_distn = "matrix"))