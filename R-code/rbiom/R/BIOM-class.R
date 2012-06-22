################################################################################
#' Build and return an instance of the biom-class.
#'
#' This is for instantiating a biom object within R (\code{\link{biom-class}}),
#' and assumes relevant data is already available in R. 
#' This is different than reading a biom file into R.
#' If you are instead interested in importing a biom file into R,
#' you should use the \code{\link{read_biom}} function. 
#' It is probably worth noting that \code{\link{read_biom}} uses this function
#' to construct its \code{\link{biom-class}} instance after parsing the raw
#' data from the biom file into R. 
#' This function is made available (exported) so that other users/developers
#' can easily represent analogous data in this structure if needed.
#'
#' \code{biom()} is a constructor method. This is the main method
#' suggested for constructing an experiment-level (\code{\link{biom-class}})
#' object from its component data.
#'
#' @usage biom(abundance, header=list(), taxonomy=NULL, sampleData=NULL, tree=NULL)
#'
#' @param abundance (REQUIRED). A \code{\link{Matrix}}-class object of abundance values.
#'  By convention, rows should represent taxa and columns different samples.
#'
#' @param header (OPTIONAL). A list or \code{NULL}. The header information describing
#'  this biom file, it's format, and other meta-data.
#'  Basically all the elements above \code{"rows"} in the biom format description.
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
#' @return An instance of the \code{\link{biom-class}}. 
#'
#' @seealso 
#' The \code{\link{read_biom}} import function.
#'
#' Accessor functions like \code{\link{header}}.
#'
#' @export
#' @examples #
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' x <- read_biom(biom_file)
#' show(x)
#' print(x)
#' header(x)
#' abundance(x)
#' taxonomy(x)
#' sampleData(x)
#' tree(x, FALSE)
biom <- function(abundance, header=list(), taxonomy=NULL, sampleData=NULL, tree=NULL){
	
	# # Some instantiation checks should go here...
	# # Some can also be wrapped into validity methods.
	biom <- new("biom",
		header     = header,
		abundance  = abundance,
		taxonomy   = taxonomy,
		sampleData = sampleData,
		tree       = tree
	)
	
	return(biom)
}
################################################################################
#' Method extensions to show for biom objects.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior. 
#'
#' @seealso \code{\link[methods]{show}}
#' 
#' @export
#' @aliases show,biom-method
#' @docType methods
#' @rdname show-methods
#' @examples
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' (x <- read_biom(biom_file) )
#' show(x)
setMethod("show", "biom", function(object){
	# cat("biom-class experiment-level object", fill=TRUE)
	# print header
	for( i in names(header(object)) ){
		# Only show non-null elements
		if( !is.null(header(object)[[i]]) ){
			cat(paste(i, header(object)[[i]][[1]], sep=": "), fill=TRUE)
		}
	}

	# print otuTable (always there).
	cat(paste("OTU Table:          [", nrow(abundance(object)), " taxa by ", 
        ncol(abundance(object)), " samples]", sep = ""), fill = TRUE
	)	

	# print Sample Data if there
	if(!is.null(sampleData(object, FALSE))){
        cat(paste("Sample Data:         [", dim(sampleData(object))[1], " samples by ", 
	        dim(sampleData(object))[2], 
            " sample variables]:", sep = ""), fill = TRUE)
	}

	# print taxonomy if present
	if(!is.null(taxonomy(object, FALSE))){
        cat(paste("Taxa Data:     [", dim(taxonomy(object))[1], " taxa by ", 
	        dim(taxonomy(object))[2], 
            " taxa variables]:", sep = ""), fill = TRUE)
	}
	
	# print tree if there
	if(!is.null(tree(object, FALSE))){
        cat(paste("Phylogenetic Tree:  [", length(tree(object)$tip.label), " tips and ", 
	        tree(object)$Nnode,
            " internal nodes]", sep = ""),
        	fill = TRUE
        )
		if( is.rooted(tree(object)) ){
			cat("                     rooted", fill=TRUE)
		} else {
			cat("                     unrooted", fill=TRUE)
		}        
	}
})
################################################################################
################################################################################
#' Universal slot accessor for the \code{\link{biom-class}}.
#'
#' This function is used by slot-specific accessors. 
#'
#' @usage access(biom, slot, errorIfNULL=FALSE) 
#'
#' @param biom (Required). \code{\link{biom-class}}.
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
	
	slots <- getSlots("biom")
	
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
#' Accessors for the \code{\link{biom-class}}. 
#' 
#' Convenience functions for accessing particular components of
#' the \code{\link{biom-class}}. 
#'
#' Some features of the biom-format can be essentially empty, and are 
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
#' @param biom (Required). \code{\link{biom-class}}.
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