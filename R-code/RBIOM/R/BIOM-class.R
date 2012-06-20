################################################################################
#' Build and return an instance of the BIOM-class.
#'
#' This is for instantiating a BIOM object within R (\code{\link{BIOM-class}}),
#' and assumes relevant data is already available in R. 
#' This is different than reading a BIOM file into R.
#' If you are instead interested in importing a BIOM file into R,
#' you should use the \code{\link{read_biom}} function. 
#' It is probably worth noting that \code{\link{read_biom}} uses this function
#' to construct its \code{\link{BIOM-class}} instance after parsing the raw
#' data from the BIOM file into R. 
#' This function is made available (exported) so that other users/developers
#' can easily represent analogous data in this structure if needed.
#'
#' \code{BIOM()} is a constructor method. This is the main method
#' suggested for constructing an experiment-level (\code{\link{BIOM-class}})
#' object from its component data.
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
#' Universal slot accessor for the \code{\link{BIOM-class}}.
#'
#' This function is used by slot-specific accessors. 
#'
#' @usage access(biom, slot, errorIfNULL=FALSE) 
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
#' @keywords internal
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
#' Accessors for the \code{\link{BIOM-class}}. 
#' 
#' Convenience functions for accessing particular components of
#' the \code{\link{BIOM-class}}. 
#'
#' Some features of the BIOM-format can be essentially empty, and are 
#' represented as \code{NULL} in R once imported. 
#' By default, these accessors will throw an error if
#' a \code{NULL} is encountered, under the assumption that you
#' are using this accessor because you expect the component to be non-empty.
#' Simply setting the second argument to \code{FALSE} relieves this and allows
#' \code{NULL} to be returned silently.
#'
#' Internally these acessors wrap the not-exported \code{\link{access}} function.
#' This helps ensure consistent behavior from accessors.
#'
#' @usage header(biom, errorIfNULL=TRUE)
#' @usage abundance(biom, errorIfNULL=TRUE)
#' @usage taxonomy(biom, errorIfNULL=TRUE)
#' @usage sampleData(biom, errorIfNULL=TRUE)
#' @usage tree(biom, errorIfNULL=TRUE)
#'
#' @param biom (Required). \code{\link{BIOM-class}}.
#'
#' @param errorIfNULL (Optional). Logical. Should the accessor stop with 
#'  an error if the slot is empty (\code{NULL})? Default \code{FALSE}. 
#'
#' @aliases accessors
#' @aliases header
#' @aliases taxonomy
#' @aliases sampleData
#' @aliases tree
#' @rdname accessor-functions
#' @export
#' @examples
#' ## Examples here.
header <- function(biom, errorIfNULL=TRUE){
	access(biom, "header", errorIfNULL)
}
################################################################################
#' @export
#' @aliases header
#' @aliases abundance
#' @rdname accessor-functions
abundance <- function(biom, errorIfNULL=TRUE){
	access(biom, "abundance", errorIfNULL)
}
################################################################################
#' @export
#' @aliases header
#' @aliases taxonomy
#' @rdname accessor-functions
taxonomy <- function(biom, errorIfNULL=TRUE){
	access(biom, "taxonomy", errorIfNULL)
}
################################################################################
#' @export
#' @aliases header
#' @aliases sampleData
#' @rdname accessor-functions
sampleData <- function(biom, errorIfNULL=TRUE){
	access(biom, "sampleData", errorIfNULL)
}
################################################################################
#' @export
#' @aliases header
#' @aliases tree
#' @rdname accessor-functions
tree <- function(biom, errorIfNULL=TRUE){
	access(biom, "tree", errorIfNULL)
}
################################################################################