#' Reformat a pedigree for plotting
#'
#' Takes a single pedigree from \code{mimsim} output and creates a pedigree structure format for plotting with \code{\link[kinship2]{plot.pedigree}}, from the \code{kinship2} package. 
#'
#' Individual IDs within \code{mimsim.object} have the format "\code{\{ped.id\}_}\code{\{Generation\}_}\code{\{Number within a generation\}}". For plotting, we remove \code{ped.id} so that inidividuals have ID of the format "\code{\{Generation\}_}\code{\{Number within a generation\}}". This is to save space as some pedigrees may be very large. 
#'
#' If IDs overlap when plotting the complete pedigree, we recommend using \code{nuclear = TRUE}. This plots the proband and their 1-st degree relatives only.
#'
#' If \code{G = TRUE} then a 2-nd binary outcome variable is created which equals \eqn{1} if an individual carries a risk variant at the major locus \eqn{G}, and \eqn{0} otherwise. What values of \eqn{G} are considered risk variants depends on the mode of inheritance of \eqn{G} with respect to the disease of interest. This is found in \code{mimsim} output \code{mimsim.object@@G_MOI}. If \code{mimsim.object@@G_MOI = "RECESSIVE"}, the 2-nd outcome will be \eqn{1} if \eqn{G = 2} and \eqn{0} otherwise. If \code{mimsim.object@@G_MOI = "DOMINANT"} or \code{mimsim.object@@G_MOI = "ALL_VARY"}, the 2-nd outcome will be \eqn{1} if \eqn{G \geq 1}{G \ge 1} and \eqn{0} otherwise. 
#'
#' Note. \code{mimsim.object@@G_MOI = "ALL_VARY"} means \code{penetrance = c(a, b, c)}, where \code{a}, \code{b} and \code{c} all differ, but \code{a} \eqn{<} \code{b}, and \code{a} \eqn{<} \code{c}.
#' 
#' If \code{G = TRUE}, but the major locus genotype has no effect, that is; \code{penetrance = c(K, K, K)}, then only disease status is considered as an outcome. This is because there are no risk variants at the major locus.
#'
#' @param  mimsim.object S4 class object resulting from a call to \code{\link{mimsim}}.
#' @param ped.id numeric. ID of the family, contained within \code{mimsim.object}, to be plotted.
#' @param nuclear logical. If \code{TRUE} then only the proband and their 1-st degree relatives (parents, siblings, offspring) will be reformatted for plotting. If \code{FALSE} then all of the family is reformatted for plotting. Default is \code{FALSE}.
#' @param G logical. If \code{TRUE} carrying a risk variant at the major locus genotype \eqn{G} will be an outcome to be plotted. This is in addition to disease status. If \code{FALSE} only disease status will be be an outcome to be plotted. Default is \code{FALSE}.
#'
#' @return \code{reformat4plotting} returns an object of class \code{pedigree}. Please see \code{\link[kinship2]{pedigree}}, in the \code{kinship2} package, for details.
#' @import kinship2
#' @examples 
#' # Using Example 1 from \code{\link{mimsim}} help page.
#'
#' ped1 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4)
#'
#' # Get a data set of class \code{pedigree} and plot it for all individuals in Pedigree 15:
#' # 1. Plot all individuals within family. 
#' # Disease status is the only outcome.
#' ped1.15 <- reformat4plotting(ped1, 15)
#' kinship2::plot.pedigree(ped1.15)
#' 
#' # 2. Plot 'nuclear' family members, with respect to the proband. 
#' # Disease status is the only outcome.
#' ped1.15 <- reformat4plotting(ped1, 15, nuclear = TRUE)
#' kinship2::plot.pedigree(ped1.15)
#' 
#' # 3. Plot 'nuclear' family members, with respect to the proband. 
#' # 2 outcomes: disease status and risk variant carrier at the major locus.
#' ped1.15 <- reformat4plotting(ped1, 15, nuclear = TRUE, G = TRUE)
#' kinship2::plot.pedigree(ped1.15)
#'
#' @export
reformat4plotting <- function(mimsim.object, ped.id, nuclear = FALSE, G = FALSE){
	if(is.na(mimsim.object@pedigree_unascertained)){
		ped.all <- mimsim.object@pedigree
	}else{
		ped.all <- rbind(mimsim.object@pedigree, mimsim.object@pedigree_unascertained)
	}
	ped.s <- ped.all[ped.all$Pedigree == ped.id, ]
	if(dim(ped.s)[1] == 0){
		stop("ped.id does not match any Pedigree IDs in mimsim.object", call.=FALSE)
	} 
	if(nuclear){
		keep.bin <- as.numeric(ped.s$Relationship == "PB")
		keep.bin[ped.s$Relationship == "S"] <- 1
		keep.bin[ped.s$Relationship == "P"] <- 1
		keep.bin[ped.s$Relationship == "O"] <- 1
		keep.bin[ped.s$Relationship == "M"] <- 1
		ped.s <- ped.s[keep.bin == 1, ]
		ped.s$Father[ped.s$Relationship == "P"] <- ""
		ped.s$Mother[ped.s$Relationship == "P"] <- ""		
	}
	id <- ped.s$ID
	sex <- ped.s$Sex
	dadid <- ped.s$Father
	dadid[dadid == ""] <- NA
	momid <- ped.s$Mother
	momid[momid == ""] <- NA
	relation <- cbind(ped.s$ID, ped.s$Mate, rep(4, dim(ped.s)[1]))
	relation <- relation[relation[,2] != "", ]
	colnames(relation) <- c("id1", "id2", "code")
	rownames(relation) <- NULL
	relation <- data.frame(relation)
	relation$code <- as.numeric(as.character(relation$code))
	dadid.unique <- unique(dadid)
	momid.unique <- unique(momid)
	bin_dad <- apply(matrix(dadid.unique), 1, function(x, id_reference) 1 - as.numeric(any(x == id_reference)), id_reference=id)
	bin_mom <- apply(matrix(momid.unique), 1, function(x, id_reference) 1 - as.numeric(any(x == id_reference)), id_reference=id)
	bin_dad[is.na(bin_dad)] <- 0
	bin_mom[is.na(bin_mom)] <- 0
	n_dad <- sum(bin_dad)
	n_mom <- sum(bin_mom)
	id <- c(id, dadid.unique[bin_dad == 1], momid.unique[bin_mom == 1])
	sex <- c(sex, rep(1, n_dad), rep(2, n_mom))
	dadid <- c(dadid, rep(NA, n_dad+n_mom))
	momid <- c(momid, rep(NA, n_dad+n_mom))
	affected<- c(ped.s$D, rep(NA, n_dad+n_mom))	
	G.risk <- NULL
	if(G){
		if(mimsim.object@G_MOI == "RECESSIVE"){G.risk <- as.numeric(ped.s$G == 2)}
		if(mimsim.object@G_MOI == "DOMINANT"){G.risk <- as.numeric(ped.s$G >= 1)}
		if(mimsim.object@G_MOI == "ALL_VARY"){G.risk <- as.numeric(ped.s$G >= 1)}	
	}
	if(!is.null(G.risk)){G.risk <- c(G.risk, rep(NA, n_dad+n_mom))}
	affected <- cbind(affected, G.risk)
	dadid <- apply(matrix(dadid), 1, function(x) if(is.na(x)){x}else{paste(strsplit(x[1], split = "_")[[1]][2], strsplit(x[1], split = "_")[[1]][3], sep="_")})
	momid <- apply(matrix(momid), 1, function(x) if(is.na(x)){x}else{paste(strsplit(x[1], split = "_")[[1]][2], strsplit(x[1], split = "_")[[1]][3], sep="_")})
	id <- apply(matrix(id), 1, function(x) if(is.na(x)){x}else{paste(strsplit(x[1], split = "_")[[1]][2], strsplit(x[1], split = "_")[[1]][3], sep="_")})
	relation[,1] <- apply(matrix(relation[,1]), 1, function(x) if(is.na(x)){x}else{paste(strsplit(x[1], split = "_")[[1]][2], strsplit(x[1], split = "_")[[1]][3], sep="_")})
	relation[,2] <- apply(matrix(relation[,2]), 1, function(x) if(is.na(x)){x}else{paste(strsplit(x[1], split = "_")[[1]][2], strsplit(x[1], split = "_")[[1]][3], sep="_")})
	dat.ped.plot <- kinship2::pedigree(id=id, dadid=dadid, momid=momid, sex=sex, affected<- affected, missid=rep("", length(id)), relation=relation)
	dat.ped.plot
}