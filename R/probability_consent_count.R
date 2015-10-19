#' Calculates the probability of consent when the family history of disease is a count variable
#'
#' \code{probability_consent_count} calculates the probability of consent for an eligible individual with a user-specified number of affected parents, siblings and grandparents. This is useful to understand the probability of consent used by \code{mimsim} when \code{FH = "count"}.
#'
#'  The output from this function is the probability of consent for an eligible individual with \code{n_D_P} affected parents, \code{n_D_S} affected siblings and \code{n_D_GP} grandparents. It calculates the probability of consent when: 1. the probability of consent is a function of the family history of disease, and, 2. Family history of disease is defined using the number of affected parents, siblings and grandparents. 
#'
#' The probability of consent is calculated using the inverse logit of a linear combination of \code{n_D_P}, \code{n_D_S} and \code{n_D_GP}. That is \eqn{pC = p(Consent | Eligible, n\_D\_P, n\_D\_S, n\_D\_GP)}{pC = p(Consent | Eligible, n_D_P, n_D_S, n_D_GP)} is the inverse of:
#' \deqn{log(\frac{pC}{1 - pC}) = a + b*n\_D\_P + c*n\_D\_S + d*n\_D\_GP}{log(pC/(1 - pC)) = a + b*n_D_P + c*n_D_S + d*n_D_GP} 
#' \code{pConsent = c(a, b, c, d)}.
#'
#' Note, the probability of consent for an ineligible individual is always \eqn{0}.
#'
#' See \code{help(mimsim)} for details on what makes an individual eligible. 
#'
#' @param n_D_P numeric. Number of affected parents. Should be \eqn{\leq 2}{\le 2}.
#' @param n_D_S numeric. Number of affected siblings.
#' @param n_D_GP numeric. Number of affected grandparents. Should be \eqn{\leq 4}{\le 4}.
#' @param pConsent vector. Specifies the coefficients used to calculate the probability of consent. It contains the intercept, and the coefficients for \code{n_D_P}, \code{n_D_S} and \code{n_D_GP}. Order dependent. Must be of length 4. See 'Details' for how consent probabilities are calculated.
#'
#' @return \code{probability_consent_count} returns an interger. This is the probability of consent for an eligible individual with \code{n_D_P} affected parents, \code{n_D_S} affected siblings and \code{n_D_GP} grandparents. See 'Details' for how this is calculated.
#'
#' @examples
#' # Example 1. 
#' # Eligible individual has 1 affected parent, 0 affected siblings and 1 affected grandparent.
#' probability_consent_count(1, 0, 1, pConsent = c(-0.5, 2, 2, 1))
#' # Eligible individual has no affected parents, siblings or grandparents.
#' probability_consent_count(0, 0, 0, pConsent = c(-0.5, 2, 2, 1))
#' # Eligible individual has 1 affected parent and no affected siblings and grandparents.
#' probability_consent_count(1, 0, 0, pConsent = c(-0.5, 2, 2, 1))
#'
#'  \dontrun{probability_consent_count(3, 0, 1, pConsent=c(-0.5, 2, 2, 1))}
#'
#' # Example 2. 
#' # Eligible individual has 1 affected parent, 1 affected siblings and 1 affected grandparent.
#' probability_consent_count(1, 1, 1, pConsent = c(0.6, 1, 2, 0.5))
#' # Eligible individual has no affected parents, siblings or grandparents.
#' probability_consent_count(0, 0, 0, pConsent = c(0.6, 1, 2, 0.5))
#' # Eligible individual has 1 affected parent and no affected siblings and grandparents.
#' probability_consent_count(1, 0, 0,  pConsent = c(0.6, 1, 2, 0.5))
#'
#'  \dontrun{probability_consent_count(0, 0, 0, pConsent=c(1, 2, 0.5))}
#'
#' @export
probability_consent_count <- function(n_D_P, n_D_S, n_D_GP, pConsent){
	if(n_D_P > 2){stop("Number of affected parents (n_D_P) should be less that or equal to 2", call. = FALSE)}
	if(n_D_GP > 4){stop("Number of affected grandparents (n_D_GP) should be less that or equal to 4", call. = FALSE)}
	if(length(pConsent) != 4){stop("pConsent must be of length 4. See help(probability_consent_count) for details", call.=FALSE)}
	sum_linear_combo <- sum(c(1, n_D_P, n_D_S, n_D_GP)*pConsent)
	exp(sum_linear_combo)/(1+exp(sum_linear_combo))
}
