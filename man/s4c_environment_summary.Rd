% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/s4c_environment_summary.R
\docType{class}
\name{s4c_environment_summary}
\alias{s4c_environment_summary}
\alias{s4c_environment_summary-class}
\title{An S4 class for \code{environment_E_summary} output}
\description{
This is an internal S4 class within \code{\link{environment_E_summary}}.
}
\section{Slots}{

\describe{
\item{\code{E}}{data.table. Summary measures for the distribution of environmental component \eqn{E} by relative type.}

\item{\code{E_G0}}{data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 0}, by relative type. That is, the distribution of \eqn{E | G = 0}, by relative type.}

\item{\code{E_G1}}{data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 1}, by relative type. That is, the distribution of \eqn{E | G = 1}, by relative type.}

\item{\code{E_G2}}{data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 2}, by relative type. That is, the distribution of \eqn{E | G = 2}, by relative type.}

\item{\code{E_Carrier}}{data.table. Summary measures for the distribution of environmental component \eqn{E} for individuals carrying a risk genotype at the major locus, by relative type. That is, the distribution of \eqn{E | G \geq 1}{E | G \ge 1}, by relative type.}
}}

