% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/environment_E_summary.R
\name{environment_E_summary}
\alias{environment_E_summary}
\title{Summary measures for environmental component, \eqn{E}, by relative type}
\usage{
environment_E_summary(mimsim.object, include = "Ascertained")
}
\arguments{
\item{mimsim.object}{S4 class object resulting from a call to \code{\link{mimsim}}.}

\item{include}{character. 2 options: 1. \code{"Ascertained"} = include families containing an eligible and consenting individual (ascertained families). These families are stored in \code{mimsim.object@pedigree}. 2. \code{"Complete"} = include families containing an eligible individual. This combines families stored in \code{mimsim.object@pedigree} (consenting, and so ascertained, families) and \code{mimsim.object@pedigree_unascertained} (non-consenting, and so unascertained, families). Default is \code{include = "Ascertained"}.}
}
\value{
\code{environment_E_summary} returns an S4 class object with 5 slots. Each slot contains a data table providing summary measures for the distribution of \eqn{E}, by relative type, for the generated families in \code{mimsim.object}.

Slots are:
\item{E}{data table. Summary measures for the distribution of environmental component \eqn{E} by relative type.}
\item{E_G0}{data table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 0}, by relative type. That is, the distribution of \eqn{E | G = 0}, by relative type.}
\item{E_G1}{data table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 1}, by relative type. That is, the distribution of \eqn{E | G = 1}, by relative type.}
\item{E_G2}{data table. Summary measures for the distribution of environmental component \eqn{E} for individuals with \eqn{G = 2}, by relative type. That is, the distribution of \eqn{E | G = 2}, by relative type.}
\item{E_Carrier}{data table. Summary measures for the distribution of environmental component \eqn{E} for individuals carrying a risk genotype at the major locus, by relative type. That is, the distribution of \eqn{E | G \geq 1}{E | G \ge 1}, by relative type.}

Summary measures are calculated using 1: all relatives of type \eqn{R} (available in slot \code{E}), and 2. relatives of type \eqn{R} stratified by major locus genotype (available in slots \code{E_G0}, \code{E_G1}, \code{E_G2} and \code{E_Carrier}).

Each data table has 7 columns. These are: 1. \code{Relationship}; relative type with respect to the proband, 2. \code{N}; the total number of relatives of type \eqn{R}, 3. \code{Mean}; mean of \eqn{E}, 4. \code{Variance}; observed variance of \eqn{E}, 5. \code{LowerIQR}; the lower interquartile range (\eqn{25^{th}}{25-th} percentile) of \eqn{E}, 6. \code{Median}; the median of \eqn{E}, and, 7. \code{UpperIQR}; the upper interquartile range (\eqn{75^{th}}{75-th} percentile) of \eqn{E}.

Relative types used are: \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood),\code{"C"} (cousins), \code{"MAv"} (aunts and uncles by marriage) and \code{"M"} (mates).
}
\description{
\code{environment_E_summary} produces summary measures for the distribution of the environmental component, \eqn{E}, by relative type. Relative specific summary measures stratified by major locus genotype, \eqn{G}, are also outputted.
}
\details{
Individuals outputted by \code{\link{mimsim}} belong to families where \eqn{\geq 1}{\ge 1} individual, the proband, is known to be affected with the disease of interest. That is, the outputted individuals have a family history of disease. In the population \eqn{E \sim N(0, r2)}{E ~ N(0, r2)}, where \eqn{r2} is the proportion of variability in liability to disease attributable to an environmental component \eqn{E}. \code{environment_E_summary} summarises the impact of a positive family history of disease on the distribution of \eqn{E} under the user-specified ascertainment scheme.

Summary measures of \eqn{E} are stratified by relative type because, for example, the effect on the distribution of \eqn{E} for offspring of the proband is expected to be different to the effect for cousins of the proband. Relative types used are: \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood),\code{"C"} (cousins), \code{"MAv"} (aunts and uncles by marriage) and \code{"M"} (mates).

Summary measures provided do not take into account correlation between individuals from the same family.

Note: If in \code{mimsim.object} non-consenting but eligible families are not stored and/ or complete ascertainment has been specified then \code{include = "Ascertained"} is equivalent to \code{include = "Complete"}.
}
\examples{
# 1. Using Example 1 from \\code{\\link{mimsim}} help page.

ped1 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4)

# Because ascertainment is complete include= "Ascertained" (default so not typed).

Esummary1 <- environment_E_summary(ped1)

# Get slot names
slotNames(Esummary1)
# Note: to access a slot type, for example:
# Esummary1@E

# 2. Using Example 4 from mimsim help page.

ped4 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4,
pConsent=c(0,1), store=TRUE)

# Summary of environmental component, just using ascertained families:

Esummary2a <- environment_E_summary(ped4)

# Summary of environmental component using ascertained and 'eligible but non-consenting' families:

Esummary2b <- environment_E_summary(ped4, include = "Complete")
}
\seealso{
\code{\link[data.table]{data.table}} for information on using data tables.

\code{\link{polygenic_A_summary}} for summary measures of polygenic component \eqn{A} and \code{\link{major_locus_G_summary}} for summary measures of major locus \eqn{G}.

\code{\link{simulation_summary}} for estimates of disease risk by relative type.
}

