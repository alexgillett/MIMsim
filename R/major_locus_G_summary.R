#' Summary measures for the major locus, \eqn{G}, by relative type
#' 
#' \code{major_locus_G_summary} produces summary tables of counts and proportions for the major locus genotypes of a user-specified relative type. Relative type is defined with respect to the proband, e.g. a sibling means 'a sibling of a proband'.
#'
#' For a user-specified relative type, summary measures produced are counts and proportions for the major locus genotypes. Options for the relative type input (\code{relative_type}) are: \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood),\code{"C"} (cousins), \code{"MAv"} (aunts and uncles by marriage) or \code{"M"} (mates). Only one relative type is specified.
#' 
#' The number of tables produced depends upon the relative type specifed. For all types the following 2 tables are outputted:
#' \itemize{
#' \item \code{G_summary1}; table containing counts and proportions of observed major locus genotypes for all relatives of type \code{relative_type},
#' \item \code{G_summary2}; table containing counts and proportions of observed major locus genotypes for all relatives of type \code{relative_type}, stratified by the major locus genotype of the proband.
#' } 
#' For relative types \code{"S"}, \code{"O"} and \code{"C"}, 2 additional tables are generated:
#' \itemize{
#' \item \code{G_summary3}; table containing counts and proportions of observed major locus genotypes for all relatives of type \code{relative_type}, stratified by sibship size.
#' \item \code{G_summary4}; table containing counts and proportions of observed major locus genotypes for all relatives of type \code{relative_type}, stratified by the major locus genotype of the proband and sibship size.
#' }
#' 
#' Summary measures provided do not take into account correlation between individuals from the same family.
#'
#' Note: If in \code{mimsim.object} non-consenting but eligible families are not stored and/ or complete ascertainment has been specified then \code{include = "Ascertained"} is equivalent to \code{include = "Complete"}.
#'
#' @param mimsim.object S4 class object resulting from a call to \code{\link{mimsim}}.
#' @param relative_type Relative type. One of: \code{"P"} (parents), \code{"S"} (siblings), \code{"O"} (offspring), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood),\code{"C"} (cousins), \code{"MAv"} (aunts and uncles by marriage) or \code{"M"} (mates).
#' @param include character. 2 options: 1. \code{"Ascertained"} = include families containing an eligible and consenting individual (ascertained families). These families are stored in \code{mimsim.object@@pedigree}. 2. \code{"Complete"} = include families containing an eligible individual. This combines families stored in \code{mimsim.object@@pedigree} (consenting, and so ascertained, families) and \code{mimsim.object@@pedigree_unascertained} (non-consenting, and so unascertained, families). Default is \code{include = "Ascertained"}.
#'
#' @return \code{major_locus_G_summary} returns an S4 class object. The number of slots depends on input variable \code{relative_type}. For \code{relative_type =} \code{"P"} (parents), \code{"GP"} (grandparents), \code{"Av"} (aunts and uncles by blood), \code{"MAv"} (aunts and uncles by marriage) or \code{"M"} (mates) there are 3 slots outputted. For \code{relative_type =} \code{"S"} (siblings), \code{"O"} (offspring) or \code{"C"} (cousins) there are 5 slots outputted.
#'
#' Slots are:
#' \item{Description}{character. Brief overview of output.}
#' \item{G_summary1}{data table. Counts and proportions of observed major locus genotypes for \code{relative_type}.}
#' \item{G_summary2}{data table. Counts and proportions of observed major locus genotypes for \code{relative_type}, stratified by the major locus genotype of the proband.}
#' \item{G_summary3}{data table. Counts and proportions of observed major locus genotypes for \code{relative_type}, stratified by sibship size. Only available for \code{relative_type =} \code{"S"} (siblings), \code{"O"} (offspring) or \code{"C"} (cousins).}
#' \item{G_summary4}{data table. Counts and proportions of observed major locus genotypes for \code{relative_type}, stratified by the major locus genotype of the proband and sibship size. Only available for \code{relative_type =} \code{"S"} (siblings), \code{"O"} (offspring) or \code{"C"} (cousins).}
#' 
#' Column headings: 
#' \itemize{
#' \item G; observed genotype of the relatives of type \code{relative_type},
#' \item Count; number of relatives of type \code{relative_type} with observed G,
#' \item Proportion; proportion of relatives of type \code{relative_type} with observed G,
#' \item Genotype_PB; the observed major locus genotypes of Probands (PB), and,
#' \item S; observed sibship size (when \code{relative_type} is \code{"S"}, \code{"O"} or \code{"C"}).
#'}
#'
#' @examples
#' # 1. Using Example 1 from mimsim help page.
#'
#' ped1 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4)
#'
#' # Summary of G for siblings of probands in ped1
#' Gsummary1S <- major_locus_G_summary(ped1, relative_type = "S")
#'
#' # Get slot names
#' slotNames(Gsummary1S)
#' # Note: to access a slot type, for example:
#' # Gsummary1S@@G_summary1
#' 
#' # Summary of G for parents of probands in ped1
#' Gsummary1P <- major_locus_G_summary(ped1, relative_type = "P")
#'
#' # 2. Using Example 4 from mimsim help page.
#' 
#' ped4 <- mimsim(n=100, penetrance=c(0.01, 0.7, 0.7), maf=0.005, h2=0.3, r2=0.4, 
#' pConsent=c(0,1), store=TRUE)
#'
#' # Summary of G for offspring of probands in ped4
#' 
#' Gsummary2aO <- major_locus_G_summary(ped1, relative_type = "O")
#' 
#' # Summary of G for offspring of probands in: 1. pedigree slot and 2. pedigree_unascertained slot
#'
#' Gsummary2bO <- major_locus_G_summary(ped1, relative_type = "O", include = "Complete")
#'
#' @import data.table
#' @import methods
#'
#' @seealso \code{\link[data.table]{data.table}} for information on using data tables.
#'	
#' \code{\link{polygenic_A_summary}} for summary measures of polygenic component \eqn{A} and \code{\link{environment_E_summary}} for summary measures of environmental component \eqn{E}.
#'
#' \code{\link{simulation_summary}} for estimates of disease risk by relative type.
#'
#' @export
major_locus_G_summary <-
function(mimsim.object, relative_type, include="Ascertained"){
	Genotype_PB <- G <- Pedigree <- G_PB <- N <- Relationship <- NULL
# relative.type options are "S", "P", "GP", "Av", "MAv", "C".
	if(!is.character(relative_type)){stop("relative_type must be a character", call. = FALSE)}
	if(include == "Ascertained"){ped <- as.data.table(mimsim.object@pedigree)}else{if(dim(mimsim.object@pedigree_unascertained)[1]>1){ped <- as.data.table(rbind(mimsim.object@pedigree, mimsim.object@pedigree_unascertained))}else{ped <- as.data.table(mimsim.object@pedigree)}}
	ped[, Genotype_PB := ped[, .(G_PB = .SD[1, G], N = .N), by=Pedigree][,rep(G_PB, N)]]
	setkey(ped, Relationship, Genotype_PB, G)
			
	## INTERNAL FUNCTIONS ##
	## Sibling genotype:
sibling_genotype_summary <- function(ped){	
	G <- Proportion <- Count <- Genotype_PB <- S <- Pedigree <- NULL
	## Overall summary
	Gsumm <- ped["S", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]

	## By PB genotype
	PB_Gsumm <- ped["S", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]
	
	## Count sibship size:
	sibs <- ped["S"]
	sibs[,S:= (.N + 1), by=Pedigree]
	setkey(sibs, Genotype_PB, S, G)	
	## By sibship size:	
	SibshipGsumm <- sibs[, .(Count = .N), by=.(S,G)][, .(G, Count, Proportion=Count/sum(Count)), by=S]
	
	## By G_PB and Sibship size?
	PBSibshipGsumm <- sibs[, .(Count = .N), by=.(Genotype_PB,S,G)][, .(G, Count, Proportion=Count/sum(Count)), by=.(Genotype_PB, S)]

  new(Class="s4c_Sibling_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype, Summary3 = Stratified by sibship size, Summary4= Stratified by PB genotype & sibship size", G_summary1 = Gsumm, G_summary2=PB_Gsumm, G_summary3=SibshipGsumm, G_summary4=PBSibshipGsumm)
}

	## Offspring genotype:
offspring_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- S <- Pedigree <- NULL
	## Overall summary
	Gsumm <- ped["O", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)]

	## By PB genotype
	PB_Gsumm <- ped["O", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]
	
	## Count sibship size:
	offs <- ped["O"]
	offs[,S:= (.N), by=Pedigree]
	setkey(offs, Genotype_PB, S, G)	
	## By sibship size:	
	SibshipGsumm <- offs[, .(Count = .N), by=.(S,G)][, .(G, Count, Proportion=Count/sum(Count)), by=S]
	
	## By G_PB and Sibship size?
	PBSibshipGsumm <- offs[, .(Count = .N), by=.(Genotype_PB,S,G)][, .(G, Count, Proportion=Count/sum(Count)), by=.(Genotype_PB, S)]

  new(Class="s4c_Offspring_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype, Summary3 = Stratified by sibship size, Summary4= Stratified by PB genotype & sibship size", G_summary1 = Gsumm, G_summary2=PB_Gsumm, G_summary3=SibshipGsumm, G_summary4=PBSibshipGsumm)
}

## parents genotype:
parent_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- NULL
	## Overall summary
	Gsumm <- ped["P", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]
	## By PB genotype
	PB_Gsumm <- ped["P", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]
	
  new(Class = "s4c_Parent_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype", G_summary1 = Gsumm, G_summary2=PB_Gsumm)
}

## aunts and uncles genotype:
aunts_uncles_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- NULL
	## Overall summary
	Gsumm <- ped["Av", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]
	## By PB genotype
	PB_Gsumm <- ped["Av", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]

  new(Class="s4c_Av_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype. NB. Av = Aunts and/or uncles by blood", G_summary1 = Gsumm, G_summary2=PB_Gsumm)
}

## married-in aunts and uncles genotype:
aunts_uncles_by_marriage_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- NULL
	## Overall summary
	Gsumm <- ped["MAv", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]
	## By PB genotype
	PB_Gsumm <- ped["MAv", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]

  new(Class = "s4c_MAv_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype. NB. MAv = Aunts and/or uncles by marriage", G_summary1 = Gsumm, G_summary2=PB_Gsumm)	
}

## grandparents genotype:
grandparent_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- NULL
	## Overall summary
	Gsumm <- ped["GP", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]
	## By PB genotype
	PB_Gsumm <- ped["GP", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]

  new(Class = "s4c_Grandparent_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype", G_summary1 = Gsumm, G_summary2=PB_Gsumm)	
}

## mates genotype:
mates_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- NULL
	## Overall summary
	Gsumm <- ped["M", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]
	## By PB genotype
	PB_Gsumm <- ped["M", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]

  new(Class = "s4c_Mate_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype", G_summary1 = Gsumm, G_summary2=PB_Gsumm)
}

## cousin genotype:
cousin_genotype_summary <- function(ped){
	G <- Proportion <- Count <- Genotype_PB <- S <- Pedigree <- Father <- NULL
	## Overall summary
	Gsumm <- ped["C", .(Count=.N), by=G, nomatch = 0][,Proportion:=Count/sum(Count)][order(G)]

	## By PB genotype
	PB_Gsumm <- ped["C", .(Count = .N), by=.(Genotype_PB, G)][order(Genotype_PB), .(G, Count, Proportion=Count/sum(Count)), by=Genotype_PB]
	
	## Count sibship size:
	cousins <- ped["C"]
	cousins[,S:=.N, by=.(Pedigree, Father)]
	setkey(cousins, Genotype_PB, S, G)	
	## By sibship size:	
	SibshipGsumm <- cousins[, .(Count = .N), by=.(S,G)][, .(G, Count, Proportion=Count/sum(Count)), by=S]
	
	## By G_PB and Sibship size?
	PBSibshipGsumm <- cousins[, .(Count = .N), by=.(Genotype_PB,S,G)][, .(G, Count, Proportion=Count/sum(Count)), by=.(Genotype_PB, S)]
		
  new(Class = "s4c_Cousin_G_summary", Description="Summary1 = Overall, Summary2 = Stratified by PB genotype, Summary3 = Stratified by sibship size, Summary4= Stratified by PB genotype & sibship size", G_summary1 = Gsumm, G_summary2=PB_Gsumm, G_summary3=SibshipGsumm, G_summary4=PBSibshipGsumm)
}
	
	if(relative_type == "S"){
		out <- sibling_genotype_summary(ped)
	}
	if(relative_type == "P"){
		out <- parent_genotype_summary(ped)
	}
	if(relative_type == "Av"){
		out <- aunts_uncles_genotype_summary(ped)
	}	
	if(relative_type == "MAv"){
		out <- aunts_uncles_by_marriage_genotype_summary(ped)
	}	
	if(relative_type == "GP"){
		out <- grandparent_genotype_summary(ped)
	}	
	if(relative_type == "C"){
		out <- cousin_genotype_summary(ped)
	}	
	if(relative_type == "O"){
		out <- offspring_genotype_summary(ped)
	}	
	if(relative_type == "M"){
		out <- mates_genotype_summary(ped)
	}	
	out
}
