#' Mixed inheritance model simulation of disease within families
#'
#' \code{mimsim} uses a mixed inheritance model to simulate disease within 3 or 4 generation families. Families are generated forwards-in-time. Outputted families contain a \emph{proband}; an eligible and consenting individual in Generation III. The probability that an eligible individual consents can be a function of the family history of disease, allowing the user to generate family data when there is ascertainment bias. 
#' 
#' \code{mimsim} inputs can broadly be categorised into 3 classes:
#' 
#' \bold{1.} Disease model inputs: \code{penetrance}, \code{maf}, \code{h2} and \code{r2}. 
#' The disease model is a mixed inheritance model (MIM); a liability threshold approach where covariates are linked to disease status, \eqn{D}, through a continuous variable called liability to disease, \eqn{L}. If \eqn{L > T} then \eqn{D = 1}, else if \eqn{L \leq T}{L \le T} then \eqn{D = 0}, where \eqn{T} is the disease threshold. Within \code{mimsim} \eqn{L} is a mixture of the following normal distributions: \deqn{L | \{G = g\} = \beta_{g} + A + E ~ \sim N(0, h2 + r2)~,}{L | (G = g) = \beta(g) + A + E  ~ N(0, h2 + r2) ,} for \eqn{g = 0, 1, 2}. Where polygenic component \eqn{A \sim N(0, h2)}{A ~ N(0, h2)}, environmental component \eqn{E \sim N(0, r2)}{E ~ N(0, r2)} and \eqn{\beta_{0} = 0}{\beta(0) = 0}.
#' Disease model inputs are used to: \itemize{
#' \item calculate \eqn{T}, \eqn{\beta_{1}}{\beta(1)} and \eqn{\beta_{2}}{\beta(2)}, and,
#' \item simulate \eqn{A}, \eqn{E} and \eqn{G} for each family member, and thus calculate \eqn{L} and \eqn{D}.
#' }
#'
#' \bold{2.} Family structure inputs: \code{sibship}, \code{mean_sibship} and \code{offspring}. 
#' Together \code{sibship} and \code{mean_sibship} specify sibship sizes across a family. \code{offspring} deterimines whether a family will span \eqn{3} or \eqn{4} generations. 
#' If \code{sibship = "Fixed"} then sibship size, \eqn{S}, is fixed \emph{within} each generation. \eqn{S} may still vary \emph{between} generations, for example, if \code{mean_sibship = c(3, 2, 1)} then \eqn{S_{II} = 3}{S(II) = 3} for all sibships in Generation II, \eqn{S_{III} = 2}{S(III) = 2} for all sibships in Generation III, and \eqn{S_{IV} = 1}{S(IV) = 1} for all sibships in Generation IV. When \code{sibship = "Fixed"} and a non-integer value is specified within \code{mean_sibship} then the ceiling function is used, i.e. \code{mean_sibship = c(2.25, 2.25, 2.25)} will become \code{mean_sibship = c(3, 3, 3)}. 
#' When \code{sibship = "Variable"}, the size of each sibship within a family is a random draw from a Poisson distribution, with a common mean within a generation specified by \code{mean_sibship}. For example, if \code{sibship = "Variable"} and \code{mean_sibship = c(3, 2, 1)} then: \eqn{S_{II} \sim Pois(3)}{S(II) ~ Pois(3)}, \eqn{S_{III} \sim Pois(2)}{S(III) ~ Pois(2)} and \eqn{S_{IV} \sim Pois(1)}{S(IV) ~ Pois(1)}.
#'
#' \bold{3.} Ascertainment inputs: \code{eligible}, \code{FH}, \code{relatives} and \code{pConsent}. 
#' These inputs set criteria for which families should be outputted. \code{mimsim} mimics some aspects of ascertainment for family studies investigating lifetime risk of disease. Users set eligibility criteria, applied in Generation III, as either: 1. \code{eligible = "Disease"}; an eligible individual has \eqn{D = 1}, or 2. \code{eligible = "Genotype"}; an eligible individual has \eqn{\{D = 1, G = risk genotype\}}{{D = 1, G = risk genotype}} (a penetrance study). Not all eligible individuals will \emph{consent} to the study, and the decision to consent may be influenced by family history of disease. That is, there may be ascertainment bias. This can be specified using \code{FH}, \code{relatives} and \code{pConsent}. \code{FH} specifies whether family history is defined as a binary variable (\code{FH = "binary"}) or a count variable (\code{FH = "count"}). \code{relatives} specifies which relatives are used in the definition of family history. This can be: 1. any combination of parents (\code{"P"}), siblings (\code{"S"}) and/ or grandparents (\code{"GP"}), or, 2. \code{NULL}. 
#' \code{pConsent} is a vector used to define the probability of consent. It's length and meaning differ depending on \code{FH}. When:
#'
#' \emph{(i)}  \code{FH = "binary"} then \code{pConsent = c(a, b)}; a vector of length 2 where \code{a} is the probability of consent given eligible and no family history of disease (\eqn{a = p(Consent | Eligible, FH = 0)}) and \code{b} is the probability of consent given eligible and there is a family history of disease (\eqn{a = p(Consent | Eligible, FH = 1)}).
#'
#' \emph{(ii)} \code{FH = "count"} then \code{pConsent = c(a, b, c, d)}; a vector of length 4 used to define the probability of consent of an individual given: 1. they are eligible and 2. the disease status of \code{relatives}. The probability of consent given eligible and the disease status of \code{relative}, is calculated using the following logistic model: 
#' \deqn{logit\big(p(Consent | Eligible, D(relatives))\big) = a + b\sum_{i=1}^{2}D_{"P", i} + c\sum_{j=1}^{S_{e} - 1}D_{"S", j} + d\sum_{k=1}^{4}D_{"GP", k} ,}{logit(p(Consent | Eligible, D(relatives))) = a + b\sumD("P", i) + c\sumD("S", j), d\sumD("GP", k) ,}
#'  where: 
#' \itemize{
#' \item \eqn{D_{"P", i}}{D("P", i)} is a binary variable, which equals \eqn{1} if the \eqn{i^{th}}{i-th} parent of the proband is affected with the disease, and \eqn{0} otherwise; for \eqn{i = 1, 2},
#' \item \eqn{D_{"S", j}}{D("S", j)} is a binary variable, which equals \eqn{1} if the \eqn{j^{th}}{j-th} sibling of the proband is affected with the disease, and \eqn{0} otherwise; for \eqn{j = 1, 2, ..., S_{e} - 1}{j = 1, 2, ..., S(e) - 1}, and where \eqn{S_{e}}{S(e)} is sibship size of the sibship containing the eligible individual under investigation,
#' \item \eqn{D_{"GP", k}}{D("GP", k)} is a binary variable, which equals \eqn{1} if the \eqn{k^{th}}{k-th} grandparent of the proband is affected with the disease, and \eqn{0} otherwise; for \eqn{k = 1, 2, 3, 4},
#' \item \code{a} defines the probability of consent for eligible individuals with no family history of disease,
#' \item \code{b} defines the increase in \eqn{logit} consent probability for each additional parent (\code{"P"}) affected with the disease,
#' \item \code{c} defines the increase in \eqn{logit} consent probability for each additional sibling (\code{"S"}) affected with the disease, and,
#' \item \code{d} defines the increase in \eqn{logit} consent probability for each additional grandparent (\code{"GP"}) affected with the disease.
#'}
#' If user is not confident with using \code{FH = "count"} due to uncertainty over what consent probabilities are being used, please see \code{\link{probability_consent_count}}. This function allows the user to explore consent probabilities for user-specified \code{pConsent} and \code{n_D_P} \eqn{= \sum_{i=1}^{2}D_{"P", i}}{= \sumD("P", i)}, \code{n_D_S} \eqn{= \sum_{j=1}^{S_{e} - 1}D_{"S", j}}{= \sumD("S", j)} and \code{n_D_GP} \eqn{= \sum_{k=1}^{4}D_{"GP", k}}{= \sumD("GP", k)}.
#'
#' Note: If \code{relatives = NULL} then everyone is assumed to have no family history of disease. In this instance to save computational time it is recommended to use complete ascertainment. The simplest way to do this is to use the default settings: \code{FH = "binary"} and \code{pConsent = c(1, 1)}.
#'
#' @param n numeric. Number of families to be generated.
#' @param penetrance vector. Lifetime penetrance of a bi-allelic major locus \eqn{G}; [\eqn{K_{0}, K_{1}, K_{2}}{K(0), K(1), K(2)}] \eqn{=} [\eqn{p(D = 1 | G = 0), p(D = 1 | G = 1), p(D = 1 | G = 2)}]. Order dependent input.
#' @param maf numeric. Risk allele frequency at the major locus, \eqn{G \sim Binom(2, maf)}{G ~ Binom(2, maf)}.
#' @param h2 numeric. Proportion of variability in liability to disease, \eqn{L}, attributable to a polygenic component, \eqn{A \sim N(0, h2)}{A ~ N(0, h2)}.
#' @param r2 numeric. Proportion of variability in liability to disease, \eqn{L}, attributable to an environmental component, \eqn{E \sim N(0, r2)}{E ~ N(0, r2)}.
#' @param eligible character. Eligibility criteria for ascertainment of a family. 2 options: 1. \code{"Disease"} = family is eligible if there exists \eqn{\geq 1}{\ge 1} individual in Generation III who has the disease, 2. \code{"Genotype"} = family is eligible if there exists \eqn{\geq 1}{\ge 1} individual in Generation III who has the disease \emph{and} carries a risk genotype at the major locus. The default is \code{"Genotype"}.
#' @param FH character. Family history of disease definition. 2 options: 1. \code{"binary"} = family history variable is binary and indicates if any relative has the disease, 2. \code{"count"} = family history variable counts the number of affected relatives. Relatives included in the definition of family history are specified using input variable \code{relatives}.
#' @param relatives character vector. Relative types used in to calculate family history of disease (see \code{FH}). Can be any combination of \code{"P"} (parents), \code{"S"} (siblings) and \code{"GP"} (grandparents), or \code{NULL}. Default is \code{c("P", "S", "GP")}. Not order dependent.
#' @param pConsent vector. Defines the conditional probability that an eligible individual will consent to the study, given family history of disease. 2 options depending on \code{FH} input. 1. If \code{FH = "binary"} (default), \code{pConsent} is a vector of length \eqn{2} and contains the probability for consent given eligible and \eqn{FH = i}, (\eqn{i = 0, 1}). 2. If \code{FH = "count"}, \code{pConsent} is a vector of length \eqn{4} and contains model parameters for a logistic model of the probability of consent given eligible and the number of affected relatives. A detailed description is given under 'Details'. Order dependent input.
#' @param sibship character. Pedigree structure input. 2 options: 1. \code{"Fixed"} = sibship size is fixed within each generation, 2. \code{"Variable"} = sibship size is variable within and between generations. Default is \code{"Variable"}.
#' @param mean_sibship vector. Average sibship size in Generation II, III and IV. Default is \code{c(2.25, 2.25, 2.25)}. Order dependent input.
#' @param offspring logical. Pedigree structure input. If \code{TRUE} then the offspring of the proband are generated and outputted families will span \eqn{4} generations. Default is \code{TRUE}.
#' @param p numeric. Probability that gene dropping is used. Default is 1; gene dropping is always used.
#' @param store logical. If \code{TRUE} eligible but non-consenting families will be stored and outputted. Default is \code{FALSE}.
#' @param maxit numeric. Simulation stops when: a. \code{n} families have been generated or b. iteration number\eqn{=} \code{maxit}. Default is \code{1e+06}. Substantially reducing this number may lead to \eqn{<} \code{n} families in \code{mimsim} output.
#'
#' @return \code{mimsim} returns an S4 class object. To view a slot type: \code{mimsim.object@@slot_name}.
#' 
#' Slots are:
#' \item{pedigree}{data frame. Pedigree data. Has \eqn{14} columns which includes variables such as individual ID (\code{ID}), ID of parents (\code{Father} and \code{Mother}), liability to disease (\code{L}) and disease status (\code{D}).}
#' \item{proband_info}{data frame. Information about probands. Has \eqn{21} columns which includes variables such as a binary variable indicating whether any 1st degree relatives (defined as parents and siblings) are affected (\code{FH1}), total number of relatives of type \code{R} (\code{n_R}) and total number affected relatives of type \code{R} (\code{n_D_R}), for \code{R} = \code{P} (parents), \code{S} (siblings), \code{O} (offspring), \code{GP} (grandparents), \code{Av} (aunts and uncles by blood), \code{MAv} (aunts and uncles by marriage), \code{C} (cousins) and \code{M} (mate).} 
#' \item{pedigree_unascertained}{data frame. Pedigree data for eligible but non-consenting families. If \code{store = FALSE} (default) or \code{store = TRUE} but no such families exist then \code{NA} is returned. Else, a data frame of same structure as \code{pedigree} is returned.}
#' \item{proband_info_unascertained}{data frame. Proband information for eligible but non-consenting families. If \code{store = FALSE} (default) or \code{store = TRUE} but no such families exist then \code{NA} is returned. Else, a data frame of same structure as \code{proband_info} is returned.}
#' \item{D_threshold}{numeric. Disease threshold, \eqn{T}; \eqn{T = \sqrt{h2 + r2}\Phi^{-1}\big(1 - K_{0}\big)}{T = \sqrt(h2 + r2)\Phi^(-1)[1 - K(0)]}, where \eqn{K_{0} = p(D = 1 | G = 0)}{K(0) = p(D = 1 | G = 0)} and is the first element in vector \code{penetrance}.}
#' \item{beta}{matrix. Mean shift in the distribution of liability to disease, \eqn{L}, when [\eqn{G = 0, G = 1, G = 2}]. \eqn{G=0} is the reference category; \eqn{\beta_{0} = 0}{\beta(0) = 0}. \eqn{\beta_{i} = T - \sqrt{h2 + r2}\Phi^{-1}\big(1 - K_{i}\big)}{\beta(i) = T - sqrt(h2 + r2)\Phi^(-1)[1 - K(i)]}, where \eqn{K_{i} = p(D = 1 | G = i)}{K(i) = p(D = 1 | G = i)} for \eqn{i = 1, 2} are specified within \code{penetrance}.}
#' \item{G_MOI}{character. Mode of inheritance for the major locus \eqn{G}. E.g. if \code{penetrance = c(0.1, 0.6, 0.6)}, then \code{G_MOI = "Dominant"}.}
#' \item{L_given_G_0}{matrix. Distribution parameters for conditional liability to disease, \eqn{L |G = 0}; \eqn{L | \{G = 0\} \sim N(0, h2 + r2)}{L | (G = 0) ~ N(0, h2 + r2)}.}
#' \item{L_given_G_1}{matrix. Distribution parameters for conditional liability to disease, \eqn{L | G = 1}; \eqn{L | \{G = 1\} \sim N(\beta_{1}, h2 + r2)}{L | (G = 1) ~ N(\beta(1), h2 + r2)}.}
#' \item{L_given_G_2}{matrix. Distribution parameters for conditional liability to disease, \eqn{L | G = 2}; \eqn{L | \{G = 2\} \sim N(\beta_{2}, h2 + r2)}{L | (G = 2) ~ N(\beta(2), h2 + r2)}.}
#' \item{L}{matrix. Distribution parameters for liability to disease, \eqn{L}. \code{L} is a mixture of the conditional distributions specified in \code{L_given_G_0}, \code{L_given_G_1} and \code{L_given_G_2}.}
#' \item{prevalence}{numeric. Lifetime population prevalence of disease.}
#' \item{...}{Input arguments. E.g. To view input variable \code{eligible} type: \code{mimsim_object_name@@eligible}.
#' Note: logical variables are changed such that \code{TRUE} = \code{Yes} and \code{FALSE} = \code{No}.}
#'
#'@section Additional information:
#' 1. Error given when inputs \code{penetrance}, \code{maf}, \code{h2} and \code{r2} result in \eqn{var(L) > 1}. 
#'
#' 2. When \code{maf = 0} regardless of specification in inputs, \code{eligible = "Disease"} and \code{p = 0}.
#'
#'
#' @examples 
#' # 1. Penetrance study: eligible = "Genotype" (default so not typed).
#' #     Complete ascertainment: no ascertainment bias.
#' #     100 families. 
#' #     Dominant major locus. prevalence = 0.01688.
#' #     Variable sibship size (sibship = "Variable") with average sibship size = 2.25 
#' #     pedigree-wide.
#' 
#' ped1 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4)
#'
#'
#' # 2. Family study of lifetime disease risk in families. eligible = "Disease".
#' #     No major locus carriers: maf = 0. 
#' #     Note 1: when maf = 0, regardless of p, no gene dropping is used.
#' #     mean_sibship=c(2.2, 2, 1): average sibship size differ across generations.
#' #     Note 2: ped2a has penetrance = c(0.01, 0.2, 0.2) and ped2b has 
#' #     penetrance = c(0.01, 0.01, 0.01). This does not matter as long as 1st element
#' #     in penetrance = p(D = 1 | G = 0) is the same.
#'
#' ped2a <- mimsim(n=100, penetrance=c(0.01, 0.2, 0.2), maf=0, h2=0.7, r2=0.2 , 
#' eligible="Disease", mean_sibship=c(2.2,2,1))
#'
#' ped2b <- mimsim(n=100, penetrance=c(0.01, 0.01, 0.01), maf=0, h2=0.7, r2=0.2 , 
#' eligible="Disease", mean_sibship=c(2.2,2,1))
#'
#'
#' # 3. Penetrance study: eligible = "Genotype".
#' #     Sibship = "Fixed", mean_sibship = c(2, 2, 2).
#' #     Recessive major locus (see penetrance).
#' ped3 <- mimsim(n = 50, penetrance = c(0.005, 0.005, 0.8), maf=0.001, h2=0.5, r2=0.1, 
#' sibship="Fixed", mean_sibship=c(2,2,2))
#'
#'
#' # 4. Example 1 with ascertainment bias: 
#' #    FH = "binary" (default so not typed),
#' #    relatives = c("P", "S", "GP") (default so not typed).
#' #    pConsent = c(0.8, 1) : eligible individuals with no family history of 
#' #    disease consent 80% of the time, eligible individuals with a family
#' #    history of disease always consent.
#' #    store = TRUE (stored eligible but non-consenting families).
#'
#' ped4 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4, 
#' pConsent=c(0.8,1), store=TRUE)
#'
#'
#' # 5. Penetrance study: eligible = "Genotype".
#' #    with ascertainment bias: 
#' #    FH = "count",
#' #    pConsent = c(0.5, 1, 1, 0.5), this means, for example:
#' # - p(Consent | Eligible, no FH) = 0.62,
#' # - p(Consent | Eligible, {1 affected parent OR 1 affected sibling}) = 0.82,
#' # - p(Consent | Eligible, {2 affected parent OR 1 affected parent, 1 affected sibling}) = 0.92,
#' # - p(Consent | Eligible, 1 affected grandparent) = 0.73,
#' # - p(Consent | Eligible, 1 affected parents, 1 affected sibling, 1 affected grandparent) = 0.95.
#' #     store = TRUE (stored eligible but non-consenting families).
#' #     Dominant major locus.
#'
#' ped5 <- mimsim(n = 100, penetrance = c(0.001, 0.55, 0.55), maf=0.005, h2 = 0.25, r2 = 0.2, 
#' FH="count", pConsent=c(0.5, 1, 1, 0.5), store=TRUE)
#'
#'
#' # 6. Example of combination of polygenic component, environmental component 
#' #     and major locus giving var(liability) > 1
#' \dontrun{pederror <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.05, h2=0.3, r2=0.4)}
#'
#' @seealso \code{\link{polygenic_A_summary}}, \code{\link{major_locus_G_summary}} and \code{\link{environment_E_summary}} produce summary measures for the polygenic component \eqn{A}, the major locus \eqn{G} and the environmental component \eqn{E} respectively.
#'
#' \code{simulation_summary} for estimates of disease risk by relative type.
#'
#' \code{\link{probability_consent_count}} allows users to explore the probability of consent used by \code{mimsim} when \code{FH = "count"}.
#'
#' @references 
#' Morton, N.E. and MacLean C.J. (1974) Analysis of family resemblance. 3. Complex segregation of quantitative traits. \emph{American Journal of Human Genetics (26)}, 489 - 503.
#' 
#' Falconer, D.S. (1960) \emph{Introduction to Quantitative Genetics.} New York, Ronald Press Co.
#'  
#' Falconer, D.S. (1965) The inheritance of liability to certain diseases, estimated from the incidence among relatives. \emph{Annuals of Human Genetics (29)}, 51 - 76.
#'
#' @import methods
#'
#' @export
mimsim <-
function(n, penetrance, maf, h2, r2, eligible = "Genotype", FH = "binary", relatives = c("P","S","GP"), pConsent=c(1,1), sibship = "Variable", mean_sibship=c(2.25, 2.25, 2.25), offspring = TRUE, p = 1, store=FALSE, maxit = 1E+06){
# Error messages to correct for 1. too few/many probabilities specified within penetrance and 2. mismatch between maf and penetrance dimensions.
	if(!is.matrix(penetrance)){penetrance <- as.matrix(t(penetrance))}
	if(dim(penetrance)[2] !=3){
		stop("penetrance requires 3 probabilities for biallelic major locus", call. = FALSE)
	}	
	if(length(maf) != dim(penetrance)[1]){
		stop("unequal number of loci specified by maf and penetrance", call. = FALSE)
	}
	
## Check that variance of conditional liability (given major locus genotype) != 0:
if((r2+h2) == 0){
	stop("h2 + r2 = 0; please specify non-zero h2 and/ or r2", call. = FALSE)
}

### FH, pConsent and relatives:
if(length(relatives) > 3){
	stop("relatives maximum length = 3; relatives is some combination of P (parents), S (siblings) and GP (grandparents) or is NULL")
}
if(FH == "binary"){
	if(length(pConsent) != 2){
		stop("when FH = binary, length(pConsent) = 2; pConsent = [p(Consent | Eligible, FH = 0), p(Consent | Eligible, FH = 1)]^T")
	}
}else{
	if(length(pConsent) != 4){
		stop("when FH = count, length(pConsent) = 4; pConsent = [beta_0, beta_P, beta_S, beta_GP]^T, betas used in a logit framework to calculate p(Consent | Eligible, number of affected relatives)")
	}
	if(!any(relatives == "GP")){pConsent[4] <- 0}
	if(!any(relatives == "P")){pConsent[2] <- 0}
	if(!any(relatives == "S")){pConsent[3] <- 0}
}
if(maf == 0){p <- 0}
### Internal functions ###
## Re-number pedigrees so that ascertained pedigrees are a continuous sequence of natural numbers.
	re_number <- function(i, pedigree_number){
		pedigree_num_i <- unique(pedigree_number)[i]
		pedigree_num_new <- rep(i, sum(pedigree_number == pedigree_num_i))
		pedigree_num_new
	}
# Function to update id names to include pedigree number info
	id_names<- function(j, ped_num, id_col){
		if(is.numeric(id_col)){
			char_ped_j <- as.character(ped_num[j])
			char_id_j <- as.character(id_col[j])
			if(as.numeric(char_id_j) <10){
				output <- paste(char_ped_j, "", sep="_")
				output <- paste(output, char_id_j, sep="")
			}
			if(as.numeric(char_id_j) == 0){
				output <- ""
			}
			if(as.numeric(char_id_j)>=10){
				output <- paste(char_ped_j, char_id_j, sep="_")
			}
		}else{
			char_id_j <- as.character(id_col[j])
			char_ped_j <- as.character(ped_num[j])
			if(id_col[j] == ""){
				output <- ""
			}else{
				output <- paste(char_ped_j, char_id_j, sep="_")
			}
		}
		output			
	}
	
## Checking the marginal distribution of liability | parameter inputs... should not be > 1...
marginal_distn_NORMAL_f <- function(h2, r2, penetrance, maf){
	D_threshold <- qnorm(1-penetrance[1,1], sd=sqrt(h2 +r2))
	beta <- c(0, (D_threshold - qnorm(1-penetrance[1,2],  sd=sqrt(h2 +r2))), (D_threshold - qnorm(1-penetrance[1,3],  sd=sqrt(h2 + r2))))
	mu <- (beta[2]*2*maf*(1-maf)) + (beta[3]*(maf^2))
	var <- h2 + r2 + ((beta[2]^2)*2*maf*(1-maf)*(1-(2*maf*(1-maf)))) + ((beta[3]^2)*(maf^2)*(1-(maf^2))) - (2*beta[2]*beta[3]*2*maf*(1-maf)*(maf^2))
	model_parameters <- matrix(c(D_threshold, beta))
	colnames(model_parameters) <- "Parameters"
	rownames(model_parameters) <- c("D_threshold", "beta0", "beta1", "beta2")
	
	marginal_distn <- matrix(c(mu, var))
	rownames(marginal_distn) <- c("Mean", "Variance")
	colnames(marginal_distn) <- c("Parameters")
	output <- new(Class = "s4c_marginal_output", model_parameters=model_parameters, marginal_distn=marginal_distn)
	output
}	

## Function to define the mode of inheritance for the major locus variable...
major_locus_heritability_f <- function(penetrance){
	if(penetrance[1,1] == penetrance[1,2]){
		if(penetrance[1,1] == penetrance[1,3]){
			output <- "NO_EFFECT"
		}else{
			output <- "RECESSIVE"
		}
	}else{
		if(penetrance[1,1] == penetrance[1,3]){
			output <- "CENTRE"
		}else{
			if(penetrance[1,2] == penetrance[1,3]){
				output <- "DOMINANT"
			}else{
			output <- "ALL_VARY"
			}
		}
	}
	output
}

## Single locus prevalence...
prevalence_1locus_f <- function(maf, penetrance){
	prevalence <- penetrance[1,1]*((1-maf)^2) + penetrance[1,2]*(2*maf*(1-maf)) + penetrance[1,3]*(maf^2)
	prevalence
}

## Parameters of liability:
params_marginal <- marginal_distn_NORMAL_f(h2, r2, penetrance, maf)
if(params_marginal@marginal_distn[2,1] > 1){
		stop("This combination of maf, penetrance, r2 and h2 is not possible within this model structure. Please change one or more of these parameters", call.=FALSE)
}

##Calculate missing input params for S4 class input object
	D_threshold <- params_marginal@model_parameters[1,1]
	beta <- params_marginal@model_parameters[2:4,1]
	beta <- t(matrix(beta))
	colnames(beta) <- c("beta_0", "beta_1", "beta_2")
	L_given_G_0 <- matrix(c(beta[1,1], (h2 + r2)))
	L_given_G_1 <- matrix(c(beta[1,2], (h2 + r2)))
	L_given_G_2 <- matrix(c(beta[1,3], (h2 + r2)))
	prevalence <- prevalence_1locus_f(maf, penetrance)
	row.names(L_given_G_0) <- c("Mean", "Variance")
	row.names(L_given_G_1) <- c("Mean", "Variance")
	row.names(L_given_G_2) <- c("Mean", "Variance")
	colnames(L_given_G_0) <- c("Distn parameters")
	colnames(L_given_G_1) <- c("Distn parameters")
	colnames(L_given_G_2) <- c("Distn parameters")
	row.names(penetrance) <- c("p(D given G = j)")
	colnames(penetrance) <- c("j=0", "j=1", "j=2")
	marginal_distn <- params_marginal@marginal_distn
	MOI_G <- major_locus_heritability_f(penetrance)
	if(store){store_in <- "Yes"}else{store_in <- "No"}
	if(offspring){off_in <- "Yes"}else{off_in <- "No"}
	
	input <- new(Class = "s4c_mimsim_input", n=n, penetrance=penetrance, maf=maf, h2=h2, r2=r2, eligible=eligible, FH=FH, relatives=relatives, pConsent=pConsent, sibship=sibship, mean_sibship=mean_sibship, offspring=off_in, p=p, store=store_in, maxit=maxit, D_threshold=D_threshold, beta=beta, G_MOI=MOI_G, L_given_G_0=L_given_G_0, L_given_G_1=L_given_G_1, L_given_G_2=L_given_G_2, L=marginal_distn, prevalence=prevalence)
		
	pedigree <- NULL
	Proband_info <- NULL	
	pb_info_unasc <- NULL
	ped_unasc <- NULL

	if(FH=="binary"){
		for(i in 1:maxit){
			gene_drop <- as.logical(rbinom(1,1,p))
			for(j in 1:maxit){
				ped_i <- mim_single_ped_binaryFH(i, penetrance=penetrance, maf=maf, h2=h2, r2=r2, eligible=eligible, sibship=sibship, mean_sibship=mean_sibship, gene_drop=gene_drop, offspring=offspring, relatives=relatives, pConsent=pConsent, store=store, maxit=maxit)
				if(!is.null(ped_i)){break}
			}	

			if(!is.null(ped_i)){
			## Pedigree 'ascertained to study'?
			## Pedigree number...
				if(ped_i[1,15] == 1){
					if(is.null(pedigree)){
						ped_num <- 1
					}else{
						ped_num <- max(pedigree[,1]) + 1
					}
				}else{
					if(is.null(ped_unasc)){
						ped_num <- n + 1
					}else{
						ped_num <- max(ped_unasc[,1]) + 1
					}
				}
			
				ped_i[,1] <- ped_num		
				pb_id <- ped_i[1, 2]
						
				## Proband_info work
				#summ_ped <- ddply(ped_i, .(Relationship), summarize, N=length(Affected), N_D=sum(Affected))
				## Siblings
				if(any(ped_i[,14] == "S")){
					sibs <- ped_i[ped_i[,14] == "S",]
					n_S <- dim(sibs)[1]
					n_S_D <- sum(sibs[,12])			
				}else{
					n_S <- 0
					n_S_D <- NA
				}
				## Parents
				parents <- ped_i[ped_i[,14] == "P",]
				n_P <- dim(parents)[1]
				n_P_D <- sum(parents[,12])
				## Aunts/ Uncles:
				if(any(ped_i[,14] == "Av")){
					Av <- ped_i[ped_i[,14] == "Av",]
					n_Av <- dim(Av)[1]
					n_Av_D <- sum(Av[,12])
					MAv <- ped_i[ped_i[,14] == "MAv",]
					n_MAv <- dim(MAv)[1]
					n_MAv_D <- sum(MAv[,12])
				}else{
					n_Av <- 0
					n_Av_D <- NA
					n_MAv <- 0
					n_MAv_D <- NA
				}
				## Grandparents:			
				gp <- ped_i[ped_i[,14] == "GP",]
				n_GP <- dim(gp)[1]
				n_GP_D <- sum(gp[,12])
				## Cousins:
				if(any(ped_i[,14] == "C")){
					cousins <- ped_i[ped_i[,14] == "C",]
					n_C <- dim(cousins)[1]
					n_C_D <- sum(cousins[,12])
				}else{
					n_C <- 0
					n_C_D <- NA
				}
				## Offspring:
				if(any(ped_i[,14] == "O")){
					off <- ped_i[ped_i[,14] == "O",]
					n_O <- dim(off)[1]
					n_O_D <- sum(off[,12])
				}else{
					n_O <- 0
					n_O_D <- NA
				}
				## Mates:
				if(any(ped_i[,14] == "M")){
					mates <- ped_i[ped_i[,14] == "M", ]
					n_M <- 1
					n_M_D <- mates[1,12]
				}else{
					n_M <- 0
					n_M_D <- NA
				}
			
				FH1 <- any(c(n_S_D,n_P_D) > 0, na.rm=T)
				FH2 <- any(c(n_GP_D, n_Av_D) > 0, na.rm=T)
				FH3 <- any(n_C_D > 0, na.rm=T)
			
				pb_info_i <- data.frame(ped_num, pb_id, FH1, FH2, FH3, n_S, n_S_D, n_P, n_P_D, n_O, n_O_D, n_GP, n_GP_D, n_Av, n_Av_D, n_MAv, n_MAv_D, n_C, n_C_D, n_M, n_M_D)

				if(ped_i[1,15] == 1){
					pedigree <- rbind(pedigree, ped_i)
					Proband_info <- rbind(Proband_info, pb_info_i)
				}else{
					ped_unasc <- rbind(ped_unasc, ped_i)
					pb_info_unasc <- rbind(pb_info_unasc, pb_info_i)
				}	
			}
			if(length(unique(pedigree[,1])) == n){break}
		}
		pedigree <- pedigree[,1:14]
		Proband_info[,3:5] <- apply(Proband_info[,3:5], 2, as.numeric)
	
		colnames(pedigree) <- c("Pedigree", "ID", "Generation", "Father", "Mother", "Mate", "Sex", "A", "G", "E", "L", "D", "child_num", "Relationship")
		pedigree[,2] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,2]))
			pedigree[,4] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,4]))
		pedigree[,5] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,5]))
		pedigree[,6] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,6]))

		Proband_info[,2] <- pedigree[pedigree[,14] == "PB", 2]
		colnames(Proband_info) <- c("Pedigree_ID", "Proband_ID", "FH1", "FH2", "FH3", "n_S", "n_D_S", "n_P", "n_D_P", "n_O", "n_D_O", "n_GP", "n_D_GP", "n_Av", "n_D_Av", "n_MAv", "n_D_MAv", "n_C", "n_D_C", "n_M", "n_D_M")
	
		if(!is.null(ped_unasc)){
			ped_unasc <- ped_unasc[,1:14]
			pb_info_unasc[,3:5] <- apply(pb_info_unasc[,3:5], 2, as.numeric)
			colnames(ped_unasc) <- colnames(pedigree)
			ped_unasc[,2] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,2]))
			ped_unasc[,4] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,4]))
			ped_unasc[,5] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,5]))
			ped_unasc[,6] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,6]))
			pb_info_unasc[,2] <- ped_unasc[ped_unasc[,14] == "PB", 2]
			colnames(pb_info_unasc) <- colnames(Proband_info)
		}else{
			ped_unasc <- data.frame(NA)
			pb_info_unasc <- data.frame(NA)
		}
		output <- new(Class = "s4c_mimsim_output", n=n, penetrance=penetrance, maf=maf, h2=h2, r2=r2, eligible=eligible, FH=FH, relatives=relatives, pConsent=pConsent, sibship=sibship, mean_sibship=mean_sibship, offspring=off_in, p=p, store=store_in, maxit=maxit, D_threshold=D_threshold, beta=beta, G_MOI=MOI_G, L_given_G_0=L_given_G_0, L_given_G_1=L_given_G_1, L_given_G_2=L_given_G_2, L=marginal_distn, prevalence=prevalence, pedigree=pedigree, proband_info=Proband_info, pedigree_unascertained=ped_unasc, proband_info_unascertained=pb_info_unasc)	
	}else{
		for(i in 1:maxit){
			gene_drop <- as.logical(rbinom(1,1,p))
			for(j in 1:maxit){
				ped_i <- mim_single_ped_countFH(i, penetrance=penetrance, maf=maf, h2=h2, r2=r2, eligible=eligible, sibship=sibship, mean_sibship=mean_sibship, gene_drop=gene_drop, offspring=offspring, pConsent=pConsent, store=store, maxit=maxit)
				if(!is.null(ped_i)){break}
			}	

			if(!is.null(ped_i)){
			## Pedigree 'ascertained to study'?
			## Pedigree number...
				if(ped_i[1,15] == 1){
					if(is.null(pedigree)){
						ped_num <- 1
					}else{
						ped_num <- max(pedigree[,1]) + 1
					}
				}else{
					if(is.null(ped_unasc)){
						ped_num <- n + 1
					}else{
						ped_num <- max(ped_unasc[,1]) + 1
					}
				}
			
				ped_i[,1] <- ped_num		
				pb_id <- ped_i[1, 2]
						
				## Proband_info work
				#summ_ped <- ddply(ped_i, .(Relationship), summarize, N=length(Affected), N_D=sum(Affected))
				## Siblings
				if(any(ped_i[,14] == "S")){
					sibs <- ped_i[ped_i[,14] == "S",]
					n_S <- dim(sibs)[1]
					n_S_D <- sum(sibs[,12])			
				}else{
					n_S <- 0
					n_S_D <- NA
				}
				## Parents
				parents <- ped_i[ped_i[,14] == "P",]
				n_P <- dim(parents)[1]
				n_P_D <- sum(parents[,12])
				## Aunts/ Uncles:
				if(any(ped_i[,14] == "Av")){
					Av <- ped_i[ped_i[,14] == "Av",]
					n_Av <- dim(Av)[1]
					n_Av_D <- sum(Av[,12])
					MAv <- ped_i[ped_i[,14] == "MAv",]
					n_MAv <- dim(MAv)[1]
					n_MAv_D <- sum(MAv[,12])
				}else{
					n_Av <- 0
					n_Av_D <- NA
					n_MAv <- 0
					n_MAv_D <- NA
				}
				## Grandparents:			
				gp <- ped_i[ped_i[,14] == "GP",]
				n_GP <- dim(gp)[1]
				n_GP_D <- sum(gp[,12])
				## Cousins:
				if(any(ped_i[,14] == "C")){
					cousins <- ped_i[ped_i[,14] == "C",]
					n_C <- dim(cousins)[1]
					n_C_D <- sum(cousins[,12])
				}else{
					n_C <- 0
					n_C_D <- NA
				}
				## Offspring:
				if(any(ped_i[,14] == "O")){
					off <- ped_i[ped_i[,14] == "O",]
					n_O <- dim(off)[1]
					n_O_D <- sum(off[,12])
				}else{
					n_O <- 0
					n_O_D <- NA
				}
				## Mates:
				if(any(ped_i[,14] == "M")){
					mates <- ped_i[ped_i[,14] == "M", ]
					n_M <- 1
					n_M_D <- mates[1,12]
				}else{
					n_M <- 0
					n_M_D <- NA
				}
			
				FH1 <- any(c(n_S_D,n_P_D) > 0, na.rm=T)
				FH2 <- any(c(n_GP_D, n_Av_D) > 0, na.rm=T)
				FH3 <- any(n_C_D > 0, na.rm=T)
			
				pb_info_i <- data.frame(ped_num, pb_id, FH1, FH2, FH3, n_S, n_S_D, n_P, n_P_D, n_O, n_O_D, n_GP, n_GP_D, n_Av, n_Av_D, n_MAv, n_MAv_D, n_C, n_C_D, n_M, n_M_D)

				if(ped_i[1,15] == 1){
					pedigree <- rbind(pedigree, ped_i)
					Proband_info <- rbind(Proband_info, pb_info_i)
				}else{
					ped_unasc <- rbind(ped_unasc, ped_i)
					pb_info_unasc <- rbind(pb_info_unasc, pb_info_i)
				}	
			}
			if(length(unique(pedigree[,1])) == n){break}
		}
		pedigree <- pedigree[,1:14]
		Proband_info[,3:5] <- apply(Proband_info[,3:5], 2, as.numeric)

		colnames(pedigree) <- c("Pedigree", "ID", "Generation", "Father", "Mother", "Mate", "Sex", "A", "G", "E", "L", "D", "child_num", "Relationship")
		pedigree[,2] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,2]))
			pedigree[,4] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,4]))
		pedigree[,5] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,5]))
		pedigree[,6] <- mapply(id_names, j=1:dim(pedigree)[1], ped_num=list(pedigree[,1]), id_col=list(pedigree[,6]))
		#Proband_info[,1] <- 1:n
		Proband_info[,2] <- pedigree[pedigree[,14] == "PB", 2]
		colnames(Proband_info) <- c("Pedigree_ID", "Proband_ID", "FH1", "FH2", "FH3", "n_S", "n_D_S", "n_P", "n_D_P", "n_O", "n_D_O", "n_GP", "n_D_GP", "n_Av", "n_D_Av", "n_MAv", "n_D_MAv", "n_C", "n_D_C", "n_M", "n_D_M")
	
		if(!is.null(ped_unasc)){
			ped_unasc <- ped_unasc[,1:14]
			pb_info_unasc[,3:5] <- apply(pb_info_unasc[,3:5], 2, as.numeric)
			colnames(ped_unasc) <- colnames(pedigree)
			ped_unasc[,2] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,2]))
			ped_unasc[,4] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,4]))
			ped_unasc[,5] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,5]))
			ped_unasc[,6] <- mapply(id_names, j=1:dim(ped_unasc)[1], ped_num=list(ped_unasc[,1]), id_col=list(ped_unasc[,6]))
			pb_info_unasc[,2] <- ped_unasc[ped_unasc[,14] == "PB", 2]
			colnames(pb_info_unasc) <- colnames(Proband_info)
		}else{
			ped_unasc <- data.frame(NA)
			pb_info_unasc <- data.frame(NA)
		}
		output <- new(Class = "s4c_mimsim_output", n=n, penetrance=penetrance, maf=maf, h2=h2, r2=r2, eligible=eligible, FH=FH, relatives=relatives, pConsent=pConsent, sibship=sibship, mean_sibship=mean_sibship, offspring=off_in, p=p, store=store_in, maxit=maxit, D_threshold=D_threshold, beta=beta, G_MOI=MOI_G, L_given_G_0=L_given_G_0, L_given_G_1=L_given_G_1, L_given_G_2=L_given_G_2, L=marginal_distn, prevalence=prevalence, pedigree=pedigree, proband_info=Proband_info, pedigree_unascertained=ped_unasc, proband_info_unascertained=pb_info_unasc)	
	}
	output
}
