#' An S4 class representing inputs used within function mimsim
#'
#' This is an internal S4 class within \code{\link{mimsim}}.
#'
#' @slot n numeric. The number of pedigrees to be generated. User specified input.
#' @slot penetrance matrix.Lifetime penetrance of a bi-allelic major locus \eqn{G}; [\eqn{K_{0}, K_{1}, K_{2}}{K(0), K(1), K(2)}] \eqn{=} [\eqn{p(D = 1 | G = 0), p(D = 1 | G = 1), p(D = 1 | G = 2)}]. Order dependent input. User specified input.
#' @slot maf numeric. Risk allele frequency at the major locus, \eqn{G \sim Binom(2, maf)}{G ~ Binom(2, maf)}. User specified input.
#' @slot h2 numeric. Proportion of variability in liability to disease, \eqn{L}, attributable to a polygenic component, \eqn{A \sim N(0, h2)}{A ~ N(0, h2)}. User specified input.
#' @slot r2 numeric. Proportion of variability in liability to disease, \eqn{L}, attributable to an environmental component, \eqn{E \sim N(0, r2)}{E ~ N(0, r2)}. User specified input.
#' @slot eligible character. Eligibility criteria for ascertainment of a family. 2 options: 1. \code{"Disease"} = family is eligible if there exists \eqn{\geq 1}{\ge 1} individual in Generation III who has the disease, 2. \code{"Genotype"} = family is eligible if there exists \eqn{\geq 1}{\ge 1} individual in Generation III who has the disease \emph{and} carries a risk genotype at the major locus. User specified input.
#' @slot FH character. Family history of disease definition. 2 options: 1. \code{"binary"} = family history variable is binary and indicates if any relative has the disease, 2. \code{"count"} = family history variable counts the number of affected relatives. Relatives included in the definition of family history are specified using input variable \code{relatives}. User specified input.
#' @slot relatives vector. Character vector. Relative types used in to calculate family history of disease (see \code{FH}). Can be any combination of \code{"P"} (parents), \code{"S"} (siblings) and \code{"GP"} (grandparents), or \code{NULL}. Not order dependent. User specified input.
#' @slot pConsent vector. Defines the conditional probability that an eligible individual will consent to the study, given family history of disease. 2 options depending on \code{FH} input. 1. If \code{FH = "binary"} (default), \code{pConsent} is a vector of length \eqn{2} and contains the probability for consent given eligible and \eqn{FH = i}, (\eqn{i = 0, 1}). 2. If \code{FH = "count"}, \code{pConsent} is a vector of length \eqn{4} and contains model parameters for a logistic model of the probability of consent given eligible and the number of affected relatives. Order dependent input. User specified input.
#' @slot sibship character. Pedigree structure input. 2 options: 1. \code{"Fixed"} = sibship size is fixed within each generation, 2. \code{"Variable"} = sibship size is variable within and between generations. User specified input.
#' @slot mean_sibship vector. Average sibship size in Generation II, III and IV. Default is \code{c(2.25, 2.25, 2.25)}. Order dependent input. User specified input.
#' @slot offspring character. \code{Yes} means that offspring of the proband will be generated. \code{No} means that offspring of the proband will not be generated. User specified input.
#' @slot p numeric.
#' @slot store character. \code{Yes} means that unascertained but eligible pedigrees will be stored in addition to ascertained pedigrees. \code{No} means otherwise. User specified input.
#' @slot maxit numeric.
#' @slot D_threshold numeric. The disease threshold: \eqn{D = 1} if \eqn{L > T} and \eqn{D = 0} otherwise, where \eqn{D} is disease status. Calculated within \code{mimsim}.
#' @slot beta matrix. The mean shifts in liability to disease \eqn{L} for having major locus genotype \eqn{= \{0, 1, 2\}} respectively. See \code{\link{mimsim}} for details. Calculated within \code{mimsim}.
#' @slot G_MOI character. The mode of inheritance of the major locus \eqn{G}. Calculated within \code{mimsim}.
#' @slot L_given_G_0 matrix. Contains the distribution parameters for the conditional distribution of liability to disease given major locus genotype \eqn{G = 0}; \eqn{L | G = 0 \sim N(\beta_{0}, h2 + r2)}{L | G = 0 ~ N(\beta(0), h2 + r2)}. Calculated within \code{mimsim}.
#' @slot L_given_G_1 matrix. Contains the distribution parameters for the conditional distribution of liability to disease given major locus genotype \eqn{G = 1}; \eqn{L | G = 1 \sim N(\beta_{1}, h2 + r2)}{L | G = 1 ~ N(\beta(1), h2 + r2)}. Calculated within \code{mimsim}.
#' @slot L_given_G_2 matrix. Contains the distribution parameters for the conditional distribution of liability to disease given major locus genotype \eqn{G = 2}; \eqn{L | G = 2 \sim N(\beta_{2}, h2 + r2)}{L | G = 2 ~ N(\beta(2), h2 + r2)}. Calculated within \code{mimsim}.
#' @slot L matrix. Contains the distribution parameters for the marginal distribution of liability to disease \eqn{L}, where \eqn{L} is a mixture of the normal distributions defined in \code{L_given_G_0}, \code{L_given_G_1} and \code{L_given_G_2}. Calculated within \code{mimsim}.
#' @slot prevalence numeric. The lifetime risk of disease in the population: \eqn{K = p(D = 1)}. Calculated within \code{mimsim}.
#'
#' @name s4c_mimsim_input
#' @rdname s4c_mimsim_input
#' @aliases s4c_mimsim_input-class
#' @exportClass s4c_mimsim_input
setClass("s4c_mimsim_input", slots = c(n="numeric", penetrance="matrix", maf="numeric", h2="numeric", r2="numeric", eligible="character", FH="character", relatives="vector", pConsent="vector", sibship="character", mean_sibship="vector", offspring="character", p="numeric", store="character", maxit="numeric", D_threshold="numeric", beta="matrix", G_MOI="character", L_given_G_0="matrix", L_given_G_1="matrix", L_given_G_2="matrix", L="matrix", prevalence="numeric"))
