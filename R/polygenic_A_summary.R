#' Summary measures for polygenic component, \eqn{A}, by relative type
#'
#' \code{polygenic_A_summary} produces summary measures for the distribution of the polygenic component, \eqn{A}, by relative type. Relative specific summary measures stratified by major locus genotype, \eqn{G}, are also outputted.
#'
#' Individuals outputted by \code{\link{mimsim}} belong to families where \eqn{\geq 1}{\ge 1} individual, the proband, is known to be affected with the disease of interest. That is, the outputted individuals have a family history of disease. In the population \eqn{A \sim N(0, h2)}{A ~ N(0, h2)}, where \eqn{h2} is the proportion of variability in liability to disease attributable to polygenic component \eqn{A}. \code{polygenic_A_summary} summarises the impact of a positive family history of disease on the distribution of \eqn{A} under the user-specified ascertainment scheme. 
#'
#' Summary measures of \eqn{A} are stratified by relative type because, for example, the effect on the distribution of \eqn{A} for offspring of the proband is expected to be different to the effect for cousins of the proband. Relative types used are: \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood),\code{"C"} (cousins), \code{"MAv"} (aunts and uncles by marriage) and \code{"M"} (mates).
#'
#' Summary measures provided do not take into account correlation between individuals from the same family.
#'
#' Note: If in \code{mimsim.object} non-consenting but eligible families are not stored and/ or complete ascertainment has been specified then \code{include = "Ascertained"} is equivalent to \code{include = "Complete"}.
#'
#' @param mimsim.object S4 class object resulting from a call to \code{mimsim}.
#' @param include character. 2 options: 1. \code{"Ascertained"} = include families containing an eligible and consenting individual (ascertained families). These families are stored in \code{mimsim.object@@pedigree}. 2. \code{"Complete"} = include families containing an eligible individual. This combines families stored in \code{mimsim.object@@pedigree} (consenting, and so ascertained, families) and \code{mimsim.object@@pedigree_unascertained} (non-consenting, and so unascertained, families). Default is \code{include = "Ascertained"}.
#'
#' @return \code{polygenic_A_summary} returns an S4 class object with 5 slots. Each slot contains a data table providing summary measures for the distribution of \eqn{A}, by relative type, for the generated families in \code{mimsim.object}. 
#'
#' Slots are:
#' \item{A}{data table. Summary measures for the distribution of polygenic component \eqn{A} by relative type.}
#' \item{A_G0}{data table. Summary measures for the distribution of polygenic component \eqn{A} for individuals with \eqn{G = 0}, by relative type. That is, the distribution of \eqn{A | G = 0}, by relative type.}
#' \item{A_G1}{data table. Summary measures for the distribution of polygenic component \eqn{A} for individuals with \eqn{G = 1}, by relative type. That is, the distribution of \eqn{A | G = 1}, by relative type.}
#' \item{A_G2}{data table. Summary measures for the distribution of polygenic component \eqn{A} for individuals with \eqn{G = 2}, by relative type. That is, the distribution of \eqn{A | G = 2}, by relative type.}
#' \item{A_Carrier}{data table. Summary measures for the distribution of polygenic component \eqn{A} for individuals carrying a risk genotype at the major locus, by relative type. That is, the distribution of \eqn{A | G \geq 1}{A | G \ge 1}, by relative type.}
#'
#' Summary measures are calculated using 1: all relatives of type \eqn{R} (available in slot \code{A}), and 2. relatives of type \eqn{R} stratified by major locus genotype (available in slots \code{A_G0}, \code{A_G1}, \code{A_G2} and \code{A_Carrier}). 
#'
#' Each data table has 7 columns. These are: 1. \code{Relationship}; relative type with respect to the proband, 2. \code{N}; the total number of relatives of type \eqn{R}, 3. \code{Mean}; mean of \eqn{A}, 4. \code{Variance}; observed variance of \eqn{A}, 5. \code{LowerIQR}; the lower interquartile range (\eqn{25^{th}}{25-th} percentile) of \eqn{A}, 6. \code{Median}; the median of \eqn{A}, and, 7. \code{UpperIQR}; the upper interquartile range (\eqn{75^{th}}{75-th} percentile) of \eqn{A}.
#'
#' Relative types used are: \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood),\code{"C"} (cousins), \code{"MAv"} (aunts and uncles by marriage) and \code{"M"} (mates).
#'
#' @examples
#' # 1. Using Example 1 from mimsim help page.
#'
#' ped1 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4)
#'
#' # Because ascertainment is complete include= "Ascertained" (default so not typed).
#'
#' Asummary1 <- polygenic_A_summary(ped1)
#'
#' # Get slot names
#' slotNames(Asummary1)
#' # Note: to access a slot type, for example:
#' # Asummary1@@A
#'
#' # 2. Using Example 4 from mimsim help page.
#' 
#' ped4 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4, 
#' pConsent=c(0,1), store=TRUE)
#'
#' # Summary of polygenic component, just using ascertained families:
#'
#' Asummary2a <- polygenic_A_summary(ped4)
#' 
#' # Summary of polygenic component using ascertained and 'eligible but non-consenting' families:
#'
#' Asummary2b <- polygenic_A_summary(ped4, include = "Complete")
#'
#' @seealso \code{\link[data.table]{data.table}} for information on using data tables.
#'	
#' \code{\link{environment_E_summary}} for summary measures of environmental component \eqn{E} and \code{\link{major_locus_G_summary}} for summary measures of major locus \eqn{G}.
#'
#' \code{\link{simulation_summary}} for estimates of disease risk by relative type.
#'
#' @import data.table
#' @import methods
#'
#' @export
polygenic_A_summary <-
function(mimsim.object, include="Ascertained"){	
	Relationship <- G <- A <- Mean <- Carrier <- NULL
	if(include == "Ascertained"){ped <- as.data.table(mimsim.object@pedigree)}else{if(dim(mimsim.object@pedigree_unascertained)[1]>1){ped <- as.data.table(rbind(mimsim.object@pedigree, mimsim.object@pedigree_unascertained))}else{ped <- as.data.table(mimsim.object@pedigree)}}
	setkey(ped, Relationship, G)		
	summA <- ped[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M")), .(N=.N, Mean=if(.N > 1){mean(A, na.rm=T)}else{-999}, Variance=var(A, na.rm=T), LowerIQR = quantile(A, na.rm=T)[2], Median = quantile(A, na.rm=T)[3], UpperIQR = quantile(A, na.rm=T)[4]), by=.EACHI][,"Relationship" := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	summA[Mean==-999, "Mean":=NA]
	## By genotype?
	G0summA <- ped[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 0), .(N=.N, Mean=if(.N > 0){mean(A, na.rm=T)}else{-999}, Variance=var(A, na.rm=T), LowerIQR = quantile(A, na.rm=T)[2], Median = quantile(A, na.rm=T)[3], UpperIQR = quantile(A, na.rm=T)[4]), by=.EACHI][,"Relationship" := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	G0summA[Mean==-999, "Mean":=NA]
	G0summA[,"G":=NULL]
	G1summA <- ped[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 1), .(N=.N, Mean=if(.N > 0){mean(A, na.rm=T)}else{-999}, Variance=var(A, na.rm=T), LowerIQR = quantile(A, na.rm=T)[2], Median = quantile(A, na.rm=T)[3], UpperIQR = quantile(A, na.rm=T)[4]), by=.EACHI][,"Relationship" := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	G1summA[Mean==-999, "Mean":=NA]
	G1summA[,"G":=NULL]
	G2summA <- ped[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 2), .(N=.N, Mean=if(.N > 0){mean(A, na.rm=T)}else{-999}, Variance=var(A, na.rm=T), LowerIQR = quantile(A, na.rm=T)[2], Median = quantile(A, na.rm=T)[3], UpperIQR = quantile(A, na.rm=T)[4]), by=.EACHI][,"Relationship" := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	G2summA[Mean==-999, "Mean":=NA]	
	G2summA[,"G":=NULL]
	ped[, "Carrier" := as.numeric(G > 0)]
	setkey(ped, Relationship, Carrier)
	CarriersummA <- ped[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 1), .(N= .N, Mean=if(.N > 0){mean(A, na.rm=T)}else{-999}, Variance=var(A, na.rm=T), LowerIQR = quantile(A, na.rm=T)[2], Median = quantile(A, na.rm=T)[3], UpperIQR = quantile(A, na.rm=T)[4]), by=.EACHI][,"Relationship" := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	CarriersummA[Mean==-999, "Mean":=NA]
	CarriersummA[,"Carrier":=NULL]

	new(Class = "s4c_polygenic_summary", A=summA, A_G0=G0summA, A_G1=G1summA, A_G2=G2summA, A_Carrier=CarriersummA)
}
