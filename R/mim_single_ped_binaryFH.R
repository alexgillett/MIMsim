#' Single pedigree mixed inheritance model simulations of disease (1)
#'
#' \code{mim_single_ped_binaryFH} is a workhorse function for \code{\link{mimsim}}. Users will not typically call \code{mim_single_ped_binaryFH} directly. It can generate a single, 3 or 4 generation family (pedigree) where \eqn{\geq 1}{\ge 1} individual in Generation III has the disease of interest. Ascertainment of a family can be influenced by disease family history, where family history is a \emph{binary} variable. That is, family history will equal 1 if there exists \eqn{\geq1}{\ge 1} affected relative, and 0 otherwise.
#'
#' When \code{mimsim} input \code{FH} is set to \code{"binary"}, \code{mim_single_ped_binaryFH} is repeatedly called until the number of families in output equals the required number \code{n}.
#'
#' Since \code{mim_single_ped_binaryFH} can produce a \code{NULL} output, it is recommended to use \code{\link{mimsim}} rather than call \code{mim_single_ped_binaryFH} directly.
#'
#' @param i numeric. Family ID.
#' @param penetrance vector. Lifetime penetrance of a bi-allelic major locus \eqn{G}; [\eqn{K_{0}, K_{1}, K_{2}}{K(0), K(1), K(2)}] \eqn{=} [\eqn{p(D = 1 | G = 0), p(D = 1 | G = 1), p(D = 1 | G = 2)}]. Order dependent input.
#' @param maf numeric. Risk allele frequency at the major locus, \eqn{G \sim Binom(2, maf)}{G ~ Binom(2, maf)}.
#' @param h2 numeric. Proportion of variability in liability to disease, \eqn{L}, attributable to a polygenic component, \eqn{A \sim N(0, h2)}{A ~ N(0, h2)}.
#' @param r2 numeric. Proportion of variability in liability to disease, \eqn{L}, attributable to an environmental component, \eqn{E \sim N(0, r2)}{E ~ N(0, r2)}.
#' @param eligible character. Eligibility criteria for ascertainment of a family. 2 options: 1. \code{"Disease"} = family is eligible if there exists \eqn{\geq 1}{\ge 1} individual in Generation III who has the disease, 2. \code{"Genotype"} = family is eligible if there exists \eqn{\geq 1}{\ge 1} individual in Generation III who has the disease \emph{and} carries a risk genotype at the major locus. The default is \code{"Genotype"}.
#' @param relatives character vector. Relative types used to calculate family history of disease (see \code{FH}). Can be any combination of \code{"P"} (parents), \code{"S"} (siblings) and \code{"GP"} (grandparents), or \code{NULL}. Default is \code{c("P", "S", "GP")}. Not order dependent.
#' @param pConsent vector. The conditional probability that an eligible individual will consent to the study, given family history of disease. \code{pConsent} is a vector of length \eqn{2} and contains the probability for consent given eligible and \eqn{FH = i}, (\eqn{i = 0, 1}). Order dependent input. Default is \code{c(1, 1)}; complete ascertainment.
#' @param sibship character. Pedigree structure input. 2 options: 1. \code{"Fixed"} = sibship size is fixed within each generation, 2. \code{"Variable"} = sibship size is variable within and between generations. Default is \code{"Variable"}.
#' @param mean_sibship vector. Average sibship size in Generation II, III and IV. Default is \code{c(2.25, 2.25, 2.25)}. Order dependent input.
#' @param offspring logical. Pedigree structure input. If \code{TRUE} then the offspring of the proband are generated and outputted families will span \eqn{4} generations. Default is \code{TRUE}.
#' @param gene_drop logical. If \code{TRUE} a risk genotype at the major locus is always assigned to 1 individual in Generation I. All other individuals in Generation I are assigned a genotype using the population distribution \eqn{G \sim Bi(2, maf)}{G ~ Bi(2, maf)}. If \code{FALSE} all individuals in Generation I are assigned a major locus genotype using the population distribution. The default is \code{TRUE}.
#' @param store logical. If \code{TRUE} eligible but non-consenting families will be stored and outputted. Default is \code{FALSE}.
#' @param maxit numeric. Maximum number of iterations. Used in internal functions. Default is \code{1e+06}.
#'
#' @return \code{mim_single_ped_binaryFH} returns either 1. a data frame with 15 columns containing information about the simulated family, or 2. \code{NULL}, meaning there was no individual in Generation III who matched the eligibility criteria specified in \code{eligible}.
#'
#' \item{Pedigree}{data frame. If there exists \eqn{\geq 1}{\ge 1} individual in Generation III matching the eligibility criteria set in \code{eligible} then a data frame containing pedigree information is outputted. Otherwise output is \code{NULL}.}
#'
#' If output is a data frame, the information it contains is:
#' \itemize{
#' \item Pedigree; the pedigree ID (\code{i}),
#' \item ID; individual ID,
#' \item Generation; the generation an individual belongs to (1 - 4);
#' \item Father; ID of father,
#' \item Mother; ID of mother,
#' \item Mate; ID of mate,
#' \item Sex; gender (1 = Male, 2 = Female),
#' \item A; polygenic component,
#' \item G; major locus genotype,
#' \item E; environmental component,
#' \item L; liability to disease,
#' \item D; disease status,
#' \item child_num; number of offspring,
#' \item Relationship; relationship to the proband,
#' \item Ascertained; given there is an eligible individual in the family, if they consent then 1 is assigned (ascertained), else 0 is assigned (unascertained).
#'}
#'
#' Because \code{mim_single_ped_binaryFH} can have \code{NULL} output, users are encouraged to use \code{\link{mimsim}} rather than call \code{mim_single_ped_binaryFH} directly. If required, a single family can be generated within \code{mimsim} by setting input \code{n = 1}. See \code{help(mimsim)} for details.
#'
#' @examples
#' 
#' single_ped <- mim_single_ped_binaryFH(i=1, penetrance = c(0.001, 0.55, 0.55), maf=0.005, 
#' h2 = 0.25, r2 = 0.2)
#'
#' # Alternatively use mimsim to guarantee a family is outputted.
#'
#' single_ped <- mimsim(n=1, penetrance = c(0.001, 0.55, 0.55), maf=0.005, h2 = 0.25, r2 = 0.2)
#'
#' # Pedigree data accessed using:
#' single_ped@@pedigree
#'
#' @import plyr
#'
#' @seealso \code{\link{mimsim}}, which is the recommended function to generate a single pedigree.
#' \code{\link{mim_single_ped_countFH}} for the generation of a single pedigree using an alternative definition of disease family history.
#'
#' @export
mim_single_ped_binaryFH <-
function(i, penetrance, maf, h2, r2, eligible = "Genotype", relatives = c("P", "S", "GP"), pConsent=c(1,1), sibship = "Variable", mean_sibship=c(2.25, 2.25, 2.25), offspring = TRUE, gene_drop = TRUE, store = FALSE, maxit = 1E+06){

###########################################	
###### SET UP #######
###########################################
## NOTE: Error messages are in wrapper function 'mim_sim'
## Internal corrections.
# If there is no major locus (maf == 0), then ascertain method must be disease
if(maf == 0){eligible = "Disease"}

if(!is.matrix(penetrance)){penetrance <- t(matrix(penetrance))}
##dim(penetrance)[2] : need this to be 3...
if(dim(penetrance)[1] == 3){penetrance <- t(penetrance)}

###########################################
##### INTERNAL FUNCTIONS#####
############################################

# Major locus type. Identifed through the penetrance function/ matrix.
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

### Assigning genotypes to an individual based on parental genotypes (Major Locus G related functions)
genotype_assignment_single <- function(x){
	output <- 0
	parent1_genotype <- x[1]
	parent2_genotype <- x[2]
	
	if(parent1_genotype == 0){
		if(parent2_genotype == 1){
			output <- rbinom(1,1,0.5)
		}
		if(parent2_genotype == 2){
			output <- 1
		}
	}
	if(parent1_genotype == 2){
		if(parent2_genotype == 0){
			output <- 1
		}
		if(parent2_genotype == 2){
			output <- 2
		}
		if(parent2_genotype == 1){
			output <- 1 + rbinom(1,1,0.5)
		}
	}
	if(parent1_genotype == 1){
		if(parent2_genotype == 0){
			output <- rbinom(1,1,0.5)
		}
		if(parent2_genotype == 2){
			output <- 1 + rbinom(1,1,0.5)
		}
		if(parent2_genotype == 1){
			output <- rbinom(1,1,0.5) + rbinom(1,1,0.5)
		}
	}
	output
}

## In Gen. 2. we want to get to the next generation but want to allow zero offspring -> sum(offspring) > 0
child_num_variable_f2 <- function(S, mean_sibship, generation, maxit){
	for(s in 1:maxit){
		out <- rpois(S, mean_sibship[generation])
		if(any(out > 0)){break}
	}
	out
}

########################################################################################
############ INITIAL VARIABLE SET UP #####################################################
########################################################################################
# Variance of the conditional liability distn, when only normally distributed environmental component is included....
varL_G <- r2 + h2
# What is the disease threshold? For discrete env may change...
D_threshold <- qnorm(1-penetrance[1,1], mean = 0, sd=sqrt(varL_G))
#Mean shift/ Effect of being G = {0,1,2}:
beta <- c(0, D_threshold - qnorm(1-penetrance[1,2], mean = 0, sd=sqrt(varL_G)), D_threshold - qnorm(1-penetrance[1,3], mean = 0, sd=sqrt(varL_G)))
## What is the effect of the major locus
MOI_G <- major_locus_heritability_f(penetrance)

########################################################################################
############ PEDIGREE SET UP ###########################################################
########################################################################################
# set up pedigree output
pedigree_i <- NULL
# 1st mating pair in Gen I:
parent_1 <- data.frame(matrix(0, nrow=1, ncol=13))
parent_2 <- data.frame(matrix(0, nrow=1, ncol=13))
names_ped <- c("Pedigree", "ID", "Generation", "Father", "Mother", "Mate", "Sex", "A", "G", "E","L","D", "child_num", "Relationship")
# Assign A (polygenic component) for founding mating pair.
parent_1[1,8] <- rnorm(1,mean=0, sd=sqrt(h2))
parent_2[1,8] <- rnorm(1,mean=0, sd=sqrt(h2))
# Assign G for founding mating pair
parent_1[1,9] <- rbinom(1,2,maf)
if(gene_drop){
	if(maf == 0){parent_2[1,9] <- rbinom(1,2,maf)}else{parent_2[1,9] <- rbinom(1,1, ((maf^2)/((maf^2) + (2*maf*(1-maf))))) + 1}
}else{
	parent_2[1,9] <- rbinom(1,2,maf)
}
# Assign unique environmental variable
parent_1[1,10] <- rnorm(1,mean=0, sd=sqrt(r2))
parent_2[1,10] <- rnorm(1,mean=0, sd=sqrt(r2))
# Assign founding parent IDs
parent_1[1,2] <- "1_01"
parent_2[1,2] <- "1_02"
# Assign founding parent generation number - 1
parent_1[1,3] <- 1		
parent_2[1,3] <- 1
# Define founding parents as each others mate
parent_1[1,6] <- parent_2[1,2]
parent_2[1,6] <- parent_1[1,2]
# Assign sex to patent 1 and assign opposite sex to parent 2- so they can be mates in the reproduction sense.
parent_1[1,7] <- rbinom(1,1,0.5)
parent_2[1,7] <- abs(1-parent_1[1,7])
# how many children do a coupling have. As founding pair need > 0 to continue:
if(sibship == "Fixed"){
	parent_1[1,13] <- ceiling(mean_sibship[1])
}else{
	for(s in 1:maxit){
		parent_1[1,13] <- rpois(1, mean_sibship[1])
		if(parent_1[1,13] > 0){break}
	}
}
parent_2[1,13] <- parent_1[1,13]
# Liability + disease status:
parent_1[,11] <- parent_1[,8] + parent_1[,10] + beta[parent_1[,9] + 1] 
parent_2[,11] <- parent_2[,8] + parent_2[,10] + beta[parent_2[,9] + 1] 
parent_1[,12] <- as.numeric(parent_1[,11] > D_threshold)
parent_2[,12] <- as.numeric(parent_2[,11] > D_threshold)

###################################################################################
# Set up of parents of mates of above gen 1 children.
###################################################################################
## Need to generate mating pairs == {parents of mates of parent1, parent2's children}
gen1_mates <- data.frame(matrix(0, nrow=(parent_1[1,13]*2), ncol = 13))
## They belong to generation 1:
gen1_mates[,3] <- 1
# IDs of parents of mates of generation 2: col 2.
if(2*parent_1[,13]+2 < 10){gen1_mates[,2] <- paste(rep(1, 2*parent_1[,13]), paste(rep("0", 2*parent_1[,13]), 2 + 1:(2*parent_1[,13]), sep=""), sep="_")}else{gen1_mates[,2] <- c(paste(rep(1, sum((2 + 1:(2*parent_1[,13])) <10)), paste(rep("0", sum((2 + 1:(2*parent_1[,13])) <10)), 2 + 1:sum((2 + 1:(2*parent_1[,13])) <10), sep=""), sep="_"), paste(rep(1, (2*parent_1[,13] - sum((2 + 1:(2*parent_1[,13])) <10))), 2 + (sum((2 + 1:(2*parent_1[,13])) <10)+1):(2*parent_1[,13]), sep="_"))}
#Assign G genotypes
gen1_mates[,9] <- rbinom((parent_1[1,13]*2),2,maf)
# Assign A
gen1_mates[,8] <- rnorm((parent_1[1,13]*2), 0, sd=sqrt(h2))
#Assign E
gen1_mates[,10] <- rnorm((parent_1[1,13]*2), mean=0, sd=sqrt(r2))
# Assign Sex. Col 7.
Sex_work <- rbinom(parent_1[, 13], 1, 0.5)
gen1_mates[,7] <- as.vector(rbind(Sex_work, abs(1 - Sex_work)))
## Mate ID. Col 6.
gen1_mates[seq(2, 2*parent_1[,13], 2), 6] <- gen1_mates[seq(1, 2*parent_1[,13], 2), 2]
gen1_mates[seq(1, 2*parent_1[,13], 2), 6] <- gen1_mates[seq(2, 2*parent_1[,13], 2), 2]
# child_num
if(sibship == "Fixed"){
	Child_work  <-  rep(ceiling(mean_sibship[1]), parent_1[, 13])
}else{
	for(s in 1: maxit){
		Child_work  <-  rpois(parent_1[, 13], mean_sibship[1])
		if(all(Child_work) > 0){break}
	}
}
gen1_mates[, 13] <- as.vector(rbind(Child_work, Child_work))
# Liab + D:
gen1_mates[,11] <- gen1_mates[,8] + gen1_mates[,10] + beta[gen1_mates[,9] + 1]
gen1_mates[,12] <- as.numeric(gen1_mates[,11] > D_threshold)
#########################################################################################
# For each subsequent generation k = 2, 3
#########################################################################################
for(k in 2:3){
	n_off <- sum(parent_1[, 13])
	## Combine parents. Create Father and Mother datasets. Make sure mates have same row number.
	parent_work <- rbind(parent_1, parent_2)
	Father_d <- parent_work[parent_work[,7] == 1 & parent_work[,13] > 0, ]
	Mother_d <- parent_work[parent_work[,7] == 0 & parent_work[,13] > 0, ]
	Father_d <- Father_d[order(Father_d[,2]), ]
	Mother_d <- Mother_d[order(Mother_d[,6]), ]

	# set up object to input children and their respective mate info into.
	children_k <- data.frame(matrix(0, nrow=n_off, ncol=13))
	mates_k <- children_k
	# ID of children. Col 2:
	if(n_off < 10){children_k[,2] <- paste(rep(k, n_off), paste(rep("0", n_off), 1:n_off, sep=""), sep="_")}else{children_k[,2] <- c(paste(rep(k, 9), paste(rep("0", 9), 1:9, sep=""), sep="_"), paste(rep(k, n_off-9), 10:n_off, sep="_"))}
	## Generation. Col 3:
	children_k[,3] <- rep(k, n_off)
	## Father ID. Col 4. And Mother ID. Col 5.
	children_k[,4] <- rep(Father_d[,2], Father_d[,13])
	children_k[,5] <- rep(Mother_d[,2], Mother_d[,13])
	## Sex. Col 7.
	children_k[,7] <- rbinom(n_off, 1, 0.5)
	## A. Col 8:
	children_k[,8] <- rep(0.5*(Father_d[,8] + Mother_d[,8]), Father_d[,13]) + rnorm(n_off, sd=sqrt(h2/2))
	## G. Col 9:
	children_k[,9] <- apply(cbind(rep(Father_d[,9], Father_d[,13]), rep(Mother_d[,9], Mother_d[,13])), 1, FUN=genotype_assignment_single)
	## E. Col 10:
	children_k[,10] <- rnorm(n_off, sd=sqrt(r2))
	## Liability + D. Col 11 + Col 12:
	children_k[,11] <- children_k[,8] + children_k[,10] + beta[children_k[,9] + 1]
	children_k[,12] <- as.numeric(children_k[,11] > D_threshold)
		## When Generation == 2:
		if(k==2){
			## Mate ID. Col 6.
			if(2*n_off < 10){children_k[,6] <- paste(rep(k, n_off), paste(rep("0", n_off), (1+n_off):(n_off*2), sep=""), sep="_")}else{if(all((n_off + 1):(2*n_off)>9)){children_k[,6] <- paste(rep(k, n_off), (n_off+1):(2*n_off),sep="_")}else{children_k[,6] <- c(paste(rep(k, sum((1+n_off):(2*n_off) <10)), paste(rep("0", sum((1+n_off):(2*n_off) <10)), (1+n_off):9, sep=""), sep="_"), paste(rep(k, n_off - sum((1+n_off):(2*n_off) <10)), 10:(2*n_off), sep="_"))}}
		## number of offspring col 13:
		children_k[,13] <- if(sibship == "Fixed"){rep(ceiling(mean_sibship[2]), sum(parent_1[, 13]))}else{child_num_variable_f2(S=sum(parent_1[, 13]), mean_sibship=mean_sibship, generation=k, maxit=maxit)}
		
		## Mate. 1. ID (col2), 2. Generation (col3), 3. Mate (col6), 4. Sex (col7), 5. child_num (col13)
		mates_k[,2] <- children_k[,6]
		mates_k[,3]<- children_k[,3]
		mates_k[,6] <- children_k[,2]
		mates_k[,7] <- abs(1-children_k[,7])
		mates_k[,13]<- children_k[,13]
		## Set up Father and Mother datasets, ordered such that Mating pairs have same row #.
		F_gen1_mates <- gen1_mates[gen1_mates[,7] == 1,]
		M_gen1_mates <- gen1_mates[gen1_mates[,7] == 0,]
		F_gen1_mates <- F_gen1_mates[order(F_gen1_mates[,2]),]
		M_gen1_mates <- M_gen1_mates[order(M_gen1_mates[,6]),]
		## 1. Father ID (col4), 2. Mother ID (col5), 3. A (col8), 4. G (col9), 5. E (col10), 6. L (col11), 7. D (col12):
		mates_k[,4] <- F_gen1_mates[,2]
		mates_k[,5] <- M_gen1_mates[,2]
		mates_k[,8] <- 0.5*(F_gen1_mates[,8] + M_gen1_mates[,8]) + rnorm(sum(parent_1[,13]), sd=sqrt(h2/2))
		mates_k[,9] <- apply(cbind(F_gen1_mates[,9], M_gen1_mates[,9]), 1, FUN=genotype_assignment_single)
		mates_k[,10] <- rnorm(sum(parent_1[, 13]), sd=sqrt(r2))
		mates_k[,11] <- mates_k[,8] + mates_k[,10] + beta[mates_k[,9] + 1]
		mates_k[,12] <- as.numeric(mates_k[,11] > D_threshold)
	}else{
		mates_k <- NULL
		# Change Gen 3 offspring to 0
		children_k[, 13] <- 0
	}

	# move current parents into the pedigree matrix as they are not needed for the next generation
	pedigree_i <- rbind(pedigree_i, parent_1, parent_2)
	# current children become the next generations parents
	parent_1 <- children_k
	parent_2 <- mates_k
	# And repeat...
}
parent_1[,13] <- NA
### function for ordering gen3 randomly later on:
sample_function <- function(x){
	sample(1:x)
}

#rm(parent_2, parent_work, F_gen1_mates, M_gen1_mates, Father_d, Mother_d)

#############################################################################################
## Create an eligibility binary variable, under ascertainment criteria (using MOI of G where applicable):
if(eligible == "Disease"){
	## When ascertained by disease Eligibility == Affected status
	Eligible <- parent_1[, 12]
}else{
	if(MOI_G == "RECESSIVE"){
		Eligible <- parent_1[, 12]*as.numeric(parent_1[,9]==2)
	}else{
		Eligible <- parent_1[, 12]*as.numeric(parent_1[,9]>0)
	}	
}

## If no one is eligible -> pedigree_out <- NULL
if(all(Eligible == 0)){
	pedigree_out <- NULL
}else{	
	Gen3_specific <- cbind(as.numeric(rep(Father_d[,12], Father_d[,13]) + rep(Mother_d[,12], Mother_d[,13]) > 0), rep(as.vector(by(parent_1[,12], parent_1[,4], sum)), Father_d[,13]), rep(mapply(function(tmp_MF, tmp_GP, z){tmp_GP[tmp_GP[,2] == tmp_MF[z,1],12] + tmp_GP[tmp_GP[,2] == tmp_MF[z,2],12] + tmp_GP[tmp_GP[,2] == tmp_MF[z,3],12] + tmp_GP[tmp_GP[,2] == tmp_MF[z,4],12]}, z=1:dim(Father_d)[1], tmp_MF=list(cbind(Father_d[,4:5], Mother_d[,4:5])), tmp_GP=list(rbind(pedigree_i[1:2, ], gen1_mates))), Father_d[,13]))
	colnames(Gen3_specific) <- c("FH_P", "FH_S", "FH_GP")
#	FH_P <- as.numeric(rep(Father_d[,12], Father_d[,13]) + rep(Mother_d[,12], Mother_d[,13]) > 0)
#	FH_S <- as.numeric(rep(as.vector(by(parent_1[,12], parent_1[,4], sum)), Father_d[,13]) > 0)
#	FH_GP <- as.numeric(rep(mapply(function(tmp_MF, tmp_GP, z){tmp_GP[tmp_GP[,2] == tmp_MF[z,1],12] + tmp_GP[tmp_GP[,2] == tmp_MF[z,2],12] + tmp_GP[tmp_GP[,2] == tmp_MF[z,3],12] + tmp_GP[tmp_GP[,2] == tmp_MF[z,4],12]}, z=1:dim(Father_d)[1], tmp_MF=list(cbind(Father_d[,4:5], Mother_d[,4:5])), tmp_GP=list(rbind(pedigree_i[1:2, ], gen1_mates))), Father_d[,13]) > 0)

	eligible_tmp <- cbind(parent_1, Gen3_specific)[Eligible == 1,]
	## Orders table by 1. Sibship ID, and 2. Individual ID.
	FH_consent <- rep(0, dim(eligible_tmp)[1])
	if(any(relatives == "P")){FH_consent <- FH_consent + eligible_tmp$FH_P}
	if(any(relatives == "S")){FH_consent <- FH_consent + as.numeric({eligible_tmp$FH_S - 1} > 0)}
	if(any(relatives == "GP")){FH_consent <- FH_consent + as.numeric(eligible_tmp$FH_GP > 0)}
	FH_consent <- as.numeric(FH_consent > 0) + 1
	Consent <- apply(matrix(FH_consent), 1, function(x) rbinom(1,1,pConsent[x]))
	eligible_tmp <- cbind(eligible_tmp, Consent)
	if(all(Consent == 0)){
		pb_c <- NULL
		pb_nc <- eligible_tmp[sample(dim(eligible_tmp)[1], 1), 2]
	}else{
		pb_c <- eligible_tmp[sample(dim(eligible_tmp)[1], 1), 2]
		pb_nc <- NULL
	}	
	#	rm(eligible_tmp)
	### Now have pb_nc and pb_c.
	### If pb_c <- NULL will only proceed if store == T
	pb_ID <- pb_c
	if(is.null(pb_c)){
		if(store){
			pb_ID <- pb_nc
		}else{
			pb_ID <- NULL
		}
	}
	if(is.null(pb_ID)){
		pedigree_out <- NULL
	}else{
		## Create Relationship Variable. Start by saying everyone == "C" (for cousin):
		Relationship <- rep("C", dim(parent_1)[1])
		## Then update S = sibs
		Relationship[parent_1[,4] == parent_1[parent_1[,2] == pb_ID,4]] <- "S"
		## Then update PB = proband
		Relationship[parent_1[,2] == pb_ID] <- "PB"
		## Already generated Aunts and Uncles (in pedigree_i); need to split by MAv and Av.
		## Already generated Av + MAv:
		ped_rel <- rep("", dim(pedigree_i)[1])
		ped_rel[as.numeric(pedigree_i[,4] == c("1_01")) + as.numeric(pedigree_i[,4] == c("1_02")) == 1] <- "Av"
		ped_rel[as.numeric(pedigree_i[,4] == c("1_01")) + as.numeric(pedigree_i[,4] == c("1_02")) == 0] <- "MAv"
		## Identify the parents + GP (they are in pedigree_i)
		ped_rel[pedigree_i[,2] == parent_1[parent_1[,2] == pb_ID, 4]] <- "P"
		ped_rel[pedigree_i[,2] == parent_1[parent_1[,2] == pb_ID, 5]] <- "P"
		ped_rel[1:2] <- "GP"
#		setkey(gen1_mates, ID)
		gen1_rel <- rep(NA, dim(gen1_mates)[1])	
		gen1_rel[gen1_mates[,2] == pedigree_i[ped_rel == "P", ][2,4]] <- "GP"
		gen1_rel[gen1_mates[,2] == pedigree_i[ped_rel == "P", ][2,5]] <- "GP"

		## 'Missing' Aunts/ Uncles/ Cousin
		if(all(gen1_mates[!is.na(gen1_rel), ][1,13] == 1)){
			## There are no missing Av/MAv/C:
			av_missing <- NULL
			c_missing <- NULL
		}else{
			n_av <- (gen1_mates[!is.na(gen1_rel), ][1,13] - 1)*2
			n2 <- dim(pedigree_i)[1] - 2

			sex_work <- rbinom(n_av/2, 1, 0.5)
			
			av_missing <- data.frame(matrix(0, nrow=n_av, ncol=14))
			## av ID, col 2:
			if(n_av+n2 < 10){av_missing[,2] <- paste(rep(2, n_av), paste(rep("0", n_av), (1+n2):(n_av+n2), sep=""), sep="_")}else{if(any((1+n2):(n2+n_av) < 10)){av_missing[,2] <- c(paste(rep(2, sum((n2+1):(n2+n_av) <10)), paste(rep("0", sum((1+n2):(n2+n_av) <10)), (1+n2):(n2+sum(((1+n2):(n2+n_av)) <10)), sep=""), sep="_"), paste(rep(2, n_av - sum(((1+n2):(n_av+n2)) <10)), (sum((1+n2):(n2+n_av) <10)+n2+1):(n2+n_av), sep="_"))}else{av_missing[,2] <- paste(rep(2, n_av), (n2+1):(n2+n_av), sep="_")}}
			## Generation, col 3:
			av_missing[,3] <- rep(2, n_av)
			## Father ID + Mother ID + mate ID. Col 4 + Col 5 + col 6:
			av_missing[,4] <- c(rep(gen1_mates[!is.na(gen1_rel), ][order(gen1_mates[!is.na(gen1_rel), 7], decreasing=T),][1,2], n_av/2), rep("",n_av/2))
			av_missing[,5] <- c(rep(gen1_mates[!is.na(gen1_rel), ][order(gen1_mates[!is.na(gen1_rel), 7], decreasing=F),][1,2], n_av/2), rep("",n_av/2))
			av_missing[,6] <- c(av_missing[(1 + n_av/2):(n_av), 2], av_missing[1:(n_av/2), 2])
			#Sex, col 7:
			av_missing[,7] <- c(sex_work, abs(1-sex_work))
			# A, col 8:
			av_missing[,8] <- c(0.5*sum(gen1_mates[!is.na(gen1_rel), 8]) + rnorm(n_av/2,mean=0,sd=sqrt(0.5*h2)), rnorm(n_av/2,mean=0,sd=sqrt(h2))) 
			# G, col 9:
			av_missing[,9] <- c(apply(t(matrix(rep(gen1_mates[!is.na(gen1_rel), 9],n_av/2), nrow=2)),1,genotype_assignment_single),rbinom(n_av/2, 2, maf))
			# E, col 10:
			av_missing[,10] <- rnorm(n_av,mean=0,sd=sqrt(r2))
			# L + d, cols 11 + 12:
			av_missing[,11] <- av_missing[,8] + av_missing[,10] + beta[av_missing[,9] + 1]
			av_missing[,12] <- as.numeric(av_missing[,11] > D_threshold)
			# Child_num:
			if(sibship=="Fixed"){av_missing[,13] <- rep(ceiling(mean_sibship[2]), n_av)}else{av_missing[,13] <- rep(rpois(n_av/2, mean_sibship[2]),2)}
			# Relationship col 14:
			av_missing[,14] <- c(rep("Av", n_av/2),rep("MAv",n_av/2))
			
			if(all(av_missing[,13]==0)){
				c_missing <- NULL
			}else{
				n_c <- sum(av_missing[,13])/2
				n3 <- dim(parent_1)[1]
				c_missing <- data.frame(matrix(0, nrow=n_c, ncol=14))
				# ID + Generation, cols 2 + 3
				if(n_c+n3 < 10){c_missing[,2] <- paste(rep(3, n_c), paste(rep("0", n_c), (1+n3):(n_c+n3), sep=""), sep="_")}else{if(any((1+n3):(n3+n_c) < 10)){c_missing[,2] <- c(paste(rep(3, sum((n3+1):(n3+n_c) <10)), paste(rep("0", sum((1+n3):(n3+n_c) <10)), (1+n3):(n3+sum(((1+n3):(n3+n_c)) <10)), sep=""), sep="_"), paste(rep(3, n_c - sum(((1+n3):(n_c+n3)) <10)), (sum((1+n3):(n3+n_c) <10)+n3+1):(n3+n_c), sep="_"))}else{c_missing[,2] <- paste(rep(3, n_c), (n3+1):(n3+n_c), sep="_")}}
				c_missing[,3] <- rep(3,n_c)
				## Father + Mother ID. Cols 4 + 5. 
				F_av <- av_missing[av_missing[,7] == 1,][order(av_missing[av_missing[,7] == 1,2]),]
				M_av <- av_missing[av_missing[,7] == 0,][order(av_missing[av_missing[,7] == 0,6]),]
				c_missing[,4] <- rep(F_av[,2],F_av[,13]) 
				c_missing[,5] <- rep(M_av[,2],M_av[,13]) 
				## 1. Sex, col 7, 2. A col 8, 3. G col 9, 4. E col 10, 5. L col 11, 6. D col 12, 7. Relationship, col 14:
				c_missing[,7] <- rbinom(n_c,1,0.5)
				c_missing[,8] <- 0.5*(rep(F_av[,8],F_av[,13]) + rep(M_av[,8],M_av[,13])) + rnorm(n_c,mean=0,sd=sqrt(h2/2))
				c_missing[,9] <- apply(cbind(rep(F_av[,9], F_av[,13]), rep(M_av[,9], M_av[,13])),1,genotype_assignment_single)
				c_missing[,10] <- rnorm(n_c,mean=0,sd=sqrt(r2))
				c_missing[,11] <- c_missing[,8] + c_missing[,10] + beta[c_missing[,9] + 1]
				c_missing[,12] <- as.numeric(c_missing[,11] > D_threshold)							
				c_missing[,14] <- rep("C",n_c)
				c_missing[,13] <- rep(NA, n_c)
			}
		}
		if(offspring){
			if(sibship == "Variable"){parent_1[parent_1[,2] ==pb_ID, 13] <- rpois(1, mean_sibship[3])}else{parent_1[parent_1[,2] ==pb_ID, 13] <- ceiling(mean_sibship[3])}
			if(is.null(c_missing)){n3 <- dim(parent_1)[1]}else{n3 <- dim(parent_1)[1] + dim(c_missing)[1]}			
			Mate1 <- data.frame(matrix(0, nrow=1, ncol=14))
			## 1. ID col 2, 2. Generation col 3, 3. Mate ID col 6, 4. Sex col 7, 5. A col 8, 6. G col 9, 7. E col 10, 8. L col 11, 9. D col 12, 10. Child num col 13, 11. Relatiopship col 14
			if(n3+1 <10){Mate1[,2] <- paste("3_0", n3+1,sep="")}else{Mate1[,2] <- paste("3", n3+1, sep="_")}
			Mate1[,3] <- 3
			Mate1[,6] <- pb_ID
			Mate1[,7] <- abs(1 - parent_1[parent_1[,2] == pb_ID, 7])
			Mate1[,8] <- rnorm(1,mean=0,sd=sqrt(h2))
			Mate1[,9] <- rbinom(1,2,maf)
			Mate1[,10] <- rnorm(1,mean=0,sd=sqrt(r2))
			Mate1[,11] <- Mate1[,8] + Mate1[,10] + beta[Mate1[,9] + 1]
			Mate1[,12] <- as.numeric(Mate1[,11] > D_threshold)
			Mate1[,13] <- parent_1[parent_1[,2] ==pb_ID, 13]
			Mate1[,14] <- "M"
			
			parent_1[parent_1[,2] == pb_ID, 6] <- Mate1[,2]			
			n4 <- Mate1[,13]
			if(all(n4==0)){
				off_1 <- NULL
			}else{
				parent_work <- rbind(Mate1[,1:13], parent_1[parent_1[,2] == pb_ID,])
				off_1 <- data.frame(matrix(0, nrow=n4, ncol=14))
				# 1. ID col 2, 2. Generation col 3, 3. Father ID col 4, 4. Mother ID col 5, 5. Sex col 7, 6. A col 8, 7. G col 9, 8. E col 10, 9. L col 11, 10. D col 12, 11. child num col 13, 12. Relationship col 14:
				if(n4 < 10){off_1[,2] <- paste(rep(4, n4), paste(rep("0", n4), 1:n4, sep=""), sep="_")}else{off_1[,2] <- c(paste(rep(4, 9), paste(rep("0", 9), 1:9, sep=""), sep="_"), paste(rep(4, n4 - 9), 10:n4, sep="_"))}
				off_1[,3] <- rep(4,n4)
				off_1[,4] <- parent_work[parent_work[,7] == 1, 2]
				off_1[,5] <- parent_work[parent_work[,7] == 0, 2]
				off_1[,7] <- rbinom(n4, 1, 0.5)
				off_1[,8] <- 0.5*sum(parent_work[,8]) + rnorm(n4, mean=0, sd=sqrt(h2/2))
				off_1[,9] <- apply(t(matrix(rep(parent_work[,9], n4), nrow=2)), 1, genotype_assignment_single)
				off_1[,10] <- rnorm(n4, mean=0, sd=sqrt(r2))
				off_1[,11] <- off_1[,8] + off_1[,10] + beta[off_1[,9] + 1]
				off_1[,12] <- as.numeric(off_1[,11] > D_threshold)
				off_1[,13] <- rep(NA, n4)
				off_1[,14] <- rep("O", n4)
			}				
		}else{
			Mate1 <- NULL
			off_1 <- NULL
		}
		parent_1 <- cbind(parent_1, Relationship)
		colnames(parent_1) <- names_ped
		pedigree_i <- cbind(pedigree_i, ped_rel)
		colnames(pedigree_i) <- names_ped
		gen1_mates <- cbind(gen1_mates[!is.na(gen1_rel), ], gen1_rel[!is.na(gen1_rel)])
		colnames(gen1_mates) <- names_ped
		if(!is.null(av_missing)){colnames(av_missing) <- names_ped}
		if(!is.null(c_missing)){colnames(c_missing) <- names_ped}
		if(!is.null(Mate1)){colnames(Mate1) <- names_ped}
		if(!is.null(off_1)){colnames(off_1) <- names_ped}
		
		pedigree_out <- rbind(parent_1, pedigree_i, gen1_mates, av_missing, c_missing, Mate1, off_1)
		order_ped <- rep(0, dim(pedigree_out)[1])
		order_ped[pedigree_out[,14] == "PB"] <- 1
		order_ped[pedigree_out[,14] == "S"] <- 2
		order_ped[pedigree_out[,14] == "P"] <- 3
		order_ped[pedigree_out[,14] == "GP"] <- 4
		order_ped[pedigree_out[,14] == "Av"] <- 5
		order_ped[pedigree_out[,14] == "MAv"] <- 6
		order_ped[pedigree_out[,14] == "C"] <- 7
		order_ped[pedigree_out[,14] == "M"] <- 8
		order_ped[pedigree_out[,14] == "O"] <- 9	
		pedigree_out <- pedigree_out[order(order_ped), ]	
		Ascertained <- rep(as.numeric(is.null(pb_nc)), dim(pedigree_out)[1])
		pedigree_out <- cbind(pedigree_out, Ascertained)
		pedigree_out[pedigree_out[,4] == 0,4:5] <- ""
		pedigree_out[pedigree_out[,6] == 0, 6] <- ""
		pedigree_out[,1] <- i
		row.names(pedigree_out) <- NULL
		pedigree_out$Sex[pedigree_out$Sex == 0] <- 2
	}
}
pedigree_out
}
