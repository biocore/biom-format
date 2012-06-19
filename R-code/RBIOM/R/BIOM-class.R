################################################################################
#' Build and return an instance of the BIOM-class.
#'
#' \code{BIOM()} is a constructor method, This is the main method
#' suggested for constructing an experiment-level (\code{\link{BIOM-class}})
#' object from its component data 
#' (component data classes: \code{\link{otuTable-class}}, \code{\link{sampleData-class}}, 
#'  \code{\link{taxonomyTable-class}}, \code{\link{phylo-class}}).
#'
#' @usage BIOM(abundance, header=list(), taxonomy=NULL, sampleData=NULL, tree=NULL)
#'
#' @param abundance (REQUIRED). A \code{\link{Matrix}}-class object of abundance values.
#'  By convention, rows should represent taxa and columns different samples.
#'
#' @param header (OPTIONAL). A list or \code{NULL}. The header information describing
#'  this BIOM file, it's format, and other meta-data.
#'  Basically all the elements above \code{"rows"} in the BIOM format description.
#'  \url{http://biom-format.org/documentation/format_versions/biom-1.0.html}.
#'  Default is \code{NULL}.
#'
#' @param taxonomy (OPTIONAL). A (usually character-based) \code{link{matrix}} of
#'  taxonomic ranks as columns and taxa as rows, or \code{NULL}.
#'  Default is \code{NULL}.
#'
#' @param sampleData (OPTIONAL). A \code{link{data.frame}} representing Sample Data, or \code{NULL}.
#'  Default is \code{NULL}
#'
#' @param tree (OPTIONAL). A \code{\link[ape]{phylo}}-class phylogenetic tree, or \code{NULL}.
#'  Default is \code{NULL}
#'
#' @return An instance of the \code{\link{BIOM-class}}. 
#'
#' @export
#' @examples #
#' # # Examples to go here.
BIOM <- function(abundance, header=list(), taxonomy=NULL, sampleData=NULL, tree=NULL){
	
	# # Some instantiation checks should go here...
	# # Some can also be wrapped into validity methods.
	biom <- new("BIOM",
		header     = header,
		abundance  = abundance,
		taxonomy   = taxonomy,
		sampleData = sampleData,
		tree       = tree
	)
	
	return(biom)
}
################################################################################
################################################################################
#' Universal slot accessor for BIOM-class.
#'
#' This function is used by slot-specific accessors. 
#'
#' @usage access(physeq, slot, errorIfNULL=FALSE)
#'
#' @param biom (Required). \code{\link{BIOM-class}}.
#'
#' @param slot (Required). A character string indicating the slot (not data class)
#'  of the component data type that is desired.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{FALSE}. 
#'
#' @return Returns the component object specified by the argument \code{slot}. 
#'  Returns NULL if slot does not exist. Returns \code{biom} as-is 
#'  if it is a component class that already matches the slot name.
#'
#' @export
#' @examples #
#' # # Examples to go here.
access <- function(biom, slot, errorIfNULL=FALSE){
	
	slots <- getSlots("BIOM")
	
	# Check that slot is actually a slot in this class.
	 if(!slot %in% slotNames(biom) ){
		out <- NULL
	} else {
		out <- eval(parse(text=paste("biom@", slot, sep=""))) 
	}
	
	# Test if you should error upon the emptiness of the slot being accessed
	if( errorIfNULL & is.null(out) ){
		stop(slot, " slot is empty.")
	}
	return(out)
}
################################################################################
header <- function(biom, errorIfNULL=TRUE){
	access(biom, "header", errorIfNULL)
}
################################################################################
abundance <- function(biom, errorIfNULL=TRUE){
	access(biom, "abundance", errorIfNULL)
}
################################################################################
taxonomy <- function(biom, errorIfNULL=TRUE){
	access(biom, "taxonomy", errorIfNULL)
}
################################################################################
sampleData <- function(biom, errorIfNULL=TRUE){
	access(biom, "sampleData", errorIfNULL)
}
################################################################################
tree <- function(biom, errorIfNULL=TRUE){
	access(biom, "tree", errorIfNULL)
}
################################################################################