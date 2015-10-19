#' Summary of mimsim simulation output
#'
#' \code{simulation_summary} produces estimates of the probability of disease and recurrence risk ratio for relatives of a proband using output from \code{mimsim}. Estimates of disease risk and recurrence risk ratios stratified by major locus genotype, \eqn{G}, are also given.
#' 
#' Note: For all summaries, relative types are \eqn{R} = \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood), \code{"C"} (cousins), \code{"MAv"} (aunts and unclues by marriage) and \code{"M"} (mates).
#' 
#' \code{simulation_summary} produces simple estimates of the probability of disease for relatives of type \eqn{R}, \eqn{K_{R}}{K(R)}, for families generated using \code{\link{mimsim}}. Estimates of \eqn{K_{R}}{K(R)} are calculated using:
#' \deqn{E[K_{R}] = E[p(D_{R} = 1 | ascertained)] = \frac{\sum_{i = 1}^{N_{R}}d_{R, i}}{N_{R}} ~ ,}{E[K(R)] = E[p(D(R) = 1 | ascertained)] = \sum d(R, i)/N(R) ,}
#' where: 
#' \itemize{
#' \item \eqn{D_{R}}{D(R)} is the disease status variable for a relative of type \eqn{R}, 
#' \item \eqn{d_{R, i}}{d(R, i)} is the observed disease status for the \eqn{i^{th}}{i-th} relative of type \eqn{R} in \code{mimsim.object}, and,
#' \item \eqn{N_{R}}{N(R)} is the total number of relatives of type \eqn{R} in \code{mimsim.object}.
#' }
#' Estimates of the recurrence risk ratio for relatives of type \eqn{R}, \eqn{\lambda_{R}}{\lambda(R)} are calculated using:
#' \deqn{E[\lambda_{R}] = \frac{E[K_{R}]}{K} ~ ,}{E[\lambda(R)] = E[K(R)]/K ,}
#' where \eqn{K = p(D = 1)} is the population prevalence of disease.
#' 
#' Estimates of disease risk to relatives of a proband, stratified by major locus genotype (\eqn{G}) are calculated using:
#' \deqn{E[K_{R | G_{R} = g}] = E[p(D_{R} = 1 | G_{R} = g, ascertained)] = \frac{\sum_{i = 1}^{N_{R | G_{R} = g}} d_{R | G_{R} = g, i}}{N_{R | G_{R} = g}} ~ ,}{E[K(R | G(R) = g)] = E[p(D(R) = 1 | G(R) = g, ascertained)] = \sum d(R | G(R) = g, i)/N(R | G(R) = g) ,}
#' for \eqn{g = 0, 1, 2}, and where:
#' \itemize{
#' \item \eqn{D_{R}}{D(R)} is the disease status variable for a relative of type \eqn{R},
#' \item \eqn{G_{R}}{G(R)} is the major locus genotype variable for a relative of type \eqn{R},
#' \item \eqn{d_{R | G_{R} = g, i}}{d(R | G(R) = g, i)} is the observed disease status for the \eqn{i^{th}}{i-th} relative of type \eqn{R} who has a major locus genotype equal to \eqn{g} in \code{mimsim.object}, and,
#' \item \eqn{N_{R | G_{R} = g}}{N(R | G(R) = g)} is the total number fo relatives of type \eqn{R} who have a major locus genotype equal to \eqn{g} in \code{mimsim.object}.
#' }
#' Estimates of recurrence risk ratio for relatives of a proband, stratified by major locus genotype (\eqn{G}) are calculated using:
#' \deqn{E[\lambda_{R | G_{R} = g}] = \frac{E[K_{R | G_{R} = g}]}{K_{g}} ~ ,}{E[\lambda(R | G(R) = g)] = E[K(R | G(R) = g)]/K(g) ,}
#' for \eqn{g = 0, 1, 2}, and where \eqn{K_{g} = p(D = 1 | G = g)}{K(g) = p(D = 1 | G = g)}.
#'
#' Summary measures provided do not take into account correlation between individuals from the same family.
#'
#' Note: If in \code{mimsim.object} non-consenting but eligible families are not stored and/ or complete ascertainment has been specified then \code{include = "Ascertained"} is equivalent to \code{include = "Complete"}.
#'
#' @param mimsim.object S4 class object resulting from a call to \code{mimsim}.
#' @param include character. 2 options: 1. \code{"Ascertained"} = include families containing an eligible and consenting individual (ascertained families). These families are stored in \code{mimsim.object@@pedigree}. 2. \code{"Complete"} = include families containing an eligible individual. This combines families stored in \code{mimsim.object@@pedigree} (consenting, and so ascertained, families) and \code{mimsim.object@@pedigree_unascertained} (non-consenting, and so unascertained, families). Default is \code{include = "Ascertained"}.
#' 
#' @return \code{simulation_summary} returns an S4 class object. To view a slot type: \code{mimsim.summary.object@@slot_name}. 
#' 
#' Slots are:
#' \item{RRR}{data table. For all family members, estimated probability of disease (column \code{K}) and recurrence risk ratio (column \code{lambda}) by relative type.}
#' \item{RRR_G0}{data table. For family members with major locus genotype \eqn{G = 0}, estimated probability of disease (column \code{K_G0}) and recurrence risk ratio (column \code{lambda_G0}) by relative type.} 
#' \item{RRR_G1}{data table. For family members with major locus genotype \eqn{G = 1}, estimated probability of disease (column \code{K_G1}) and recurrence risk ratio (column \code{lambda_G1}) by relative type.}
#' \item{RRR_G2}{data table. For family members with major locus genotype \eqn{G = 2}, estimated probability of disease (column \code{K_G2}) and recurrence risk ratio (column \code{lambda_G2}) by relative type.}
#' \item{RRR_carrier}{data table. For family members with major locus genotype \eqn{G \geq 1}{G \ge 1}, estimated probability of disease (column \code{K_carrier}) and recurrence risk ratio (column \code{lambda_carrier}) by relative type.}
#' \item{Count_summary}{data table. For each relative type, \code{Count_summary} provides a count of the number of families containing \eqn{\geq 1}{\ge 1} relative of type \eqn{R} (column \code{N_peds}), and the number of families containing \eqn{\geq 1}{\ge 1} affected relative of type \eqn{R} (column \code{N_peds_ge1_D}).}
#'
#' See 'Details' section for relative types included.
#'
#' @examples
#' # 1. Using Example 1 from mimsim help page.
#'
#' ped1 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4)
#'
#' summary1 <- simulation_summary(ped1)
#'
#' # Get slot names
#' slotNames(summary1)
#' # Note: to access a slot type, for example:
#' # summary1@@RRR
#'
#' # 2. Using Example 4 from mimsim help page.
#' 
#' ped4 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4, 
#' pConsent=c(0,1), store=TRUE)
#'
#' # Estimates of disease risk and recurrence risk ratio for ascertained families
#'
#' summary2a <- simulation_summary(ped4)
#'
#' # Estimates of disease risk and recurrence risk ratio for ascertained families and
#' # eligible but non-consenting families
#'
#' summary2b <- simulation_summary(ped4, include = "Complete")
#'
#' @import data.table
#' @import methods
#'
#' @seealso \code{\link[data.table]{data.table}} for information on using data tables.
#'	
#' \code{\link{polygenic_A_summary}}, \code{\link{major_locus_G_summary}} and \code{\link{environment_E_summary}} produce summary measures for the polygenic component \eqn{A}, the major locus \eqn{G} and the environmental component \eqn{E} respectively.
#'  
#' @export 
simulation_summary <-
function(mimsim.object, include = "Ascertained"){
	Relationship <- G <- NULL
## Other include option == "Complete"
if(include == "Ascertained"){ped <- as.data.table(mimsim.object@pedigree)}else{if(dim(mimsim.object@pedigree_unascertained)[1]>1){ped <- as.data.table(rbind(mimsim.object@pedigree, mimsim.object@pedigree_unascertained))}else{ped <- as.data.table(mimsim.object@pedigree)}}

setkey(ped, Relationship, G)

RRR_function <- function(pedigree, prevalence, penetrance, maf){
	Relationship <- N_G1 <- N_G2 <- K_G1 <- K_G2 <- K_carrier <- lambda_carrier <- NULL
	RRR_out <- pedigree[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M")), .(N = .N, K=sum(D)/.N, lambda=(sum(D)/.N)/prevalence), by=.EACHI][,Relationship := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	
	RRR_0 <- pedigree[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 0), .(N_G0 = .N, K_G0=sum(D)/.N, lambda_G0 = (sum(D)/.N)/penetrance[1,1]), by=.EACHI][,Relationship := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	RRR_1 <- pedigree[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 1), .(N_G1 = .N, K_G1=sum(D)/.N, lambda_G1 = (sum(D)/.N)/penetrance[1,2]), by=.EACHI][,Relationship := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	RRR_2 <- pedigree[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), 2), .(N_G2 = .N, K_G2=sum(D)/.N, lambda_G2 = (sum(D)/.N)/penetrance[1,3]), by=.EACHI][,Relationship := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
		
	RRR_carrier <- data.table(Relationship = c("S", "P", "O", "GP", "Av", "C", "MAv", "M"), N_carrier = RRR_1[,N_G1] + RRR_2[,N_G2], K_carrier = rowSums(cbind(RRR_1[,N_G1*K_G1], RRR_2[,N_G2*K_G2]), na.rm=T)/(RRR_1[,N_G1] + RRR_2[,N_G2]))[,K_carrier:=(K_carrier <= 1)*K_carrier][,lambda_carrier := K_carrier/((prevalence - (penetrance[1,1]*((1-maf)^2)))/(1 - ((1-maf)^2)))]

	output <- new(Class = "s4c_RRR_output", RRR = RRR_out, RRR_G0 = RRR_0, RRR_G1 = RRR_1, RRR_G2 = RRR_2, RRR_carrier = RRR_carrier)
	output
}

count_summary_function <- function(pedigree){
		Pedigree <- Relationship <- N_peds <- Proportion <- N_peds_ge1_D <- NULL
		count_out <- pedigree[.(c("S", "P", "O", "GP", "Av", "C", "MAv", "M")), .(N_peds = length(unique(Pedigree)), N_peds_ge1_D = length(unique(Pedigree[D==1]))), by=.EACHI][,N_peds := (1-is.na(Relationship))*N_peds][, N_peds_ge1_D := (1-is.na(Relationship))*N_peds_ge1_D][, Proportion := N_peds_ge1_D/N_peds][, Proportion := (Proportion <= 1)*Proportion][,Relationship := c("S", "P", "O", "GP", "Av", "C", "MAv", "M")]
	count_out
}

RRR_summary <- RRR_function(pedigree=ped, prevalence=mimsim.object@prevalence, penetrance=mimsim.object@penetrance, maf=mimsim.object@maf)
Count_summary <- count_summary_function(pedigree=ped)
	
	output <- new(Class = "s4c_simulation_summary", RRR=RRR_summary@RRR, RRR_G0=RRR_summary@RRR_G0, RRR_G1=RRR_summary@RRR_G1, RRR_G2=RRR_summary@RRR_G2, RRR_carrier=RRR_summary@RRR_carrier, Count_summary=Count_summary)
	output
}
