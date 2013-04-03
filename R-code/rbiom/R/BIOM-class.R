################################################################################
#' Build and return an instance of the biom-class.
#'
#' This is for instantiating a biom object within R (\code{\link{biom-class}}),
#' and assumes relevant data is already available in R. 
#' This is different than reading a biom file into R.
#' If you are instead interested in importing a biom file into R,
#' you should use the \code{\link{read_biom}} function. 
#' This function is made available (exported) so that 
#' advanced-users/developers
#' can easily represent analogous data in this structure if needed.
#' However, most users are expected to instead rely on the
#' \code{\link{read_biom}} function for data import, followed by 
#' \code{\link{accessors}} (accessor functions) that extract R-friendly
#'  subsets of the data stored in the biom-format derived list.
#'
#' \code{biom()} is a constructor method. This is the main method
#' suggested for constructing an experiment-level (\code{\link{biom-class}})
#' object from its component data.
#'
#' @usage biom(abundance, header=list(), taxonomy=NULL, sampleData=NULL, tree=NULL)
#'
#' @param x (REQUIRED). A named list conforming to conventions arising from 
#'  the \code{\link{fromJSON}} function reading a biom-format file with 
#'  default settings. See \code{\link{read_biom}} for more details, and 
#'  \code{\link{accessors}} (accessor functions) that extract R-friendly
#'  subsets of the data stored in the biom-format derived list.
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
#' # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' x <- read_biom(biom_file)
#' show(x)
#' print(x)
#' header(x)
#' abundance(x)
#' taxonomy(x)
#' sampleData(x)
biom <- function(x){
	
	# Some instantiation checks chould go here, or wrap them in validity methods.
	# Depends on how strict the check should be.
	biom <- new("biom", x)
	
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
	# print header
	for( i in names(header(object)) ){
		# Only show non-null elements
		if( !is.null(header(object)[[i]]) ){
			cat(paste(i, header(object)[[i]][[1]], sep=": "), fill=TRUE)
		}
	}

})
################################################################################
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
#' @usage header(x)
#' @usage abundance(x, parallel=FALSE)
#' @usage taxonomy(x, parallel=FALSE)
#' @usage sampleData(x, parallel=FALSE)
#'
#' @param x (Required). An instance of the \code{\link{biom-class}}.
#' @param parallel (Optional). Logical. Whether to perform the accession parsing
#'  using a parallel-computing backend supported by the \code{\link{plyr-package}}
#'  via the \code{\link{foreach-package}}. Note: At the moment, the header
#'  accessor does not need nor does it support parallel-computed parsing.
#'
#' @aliases accessors
#' @aliases header
#' @aliases taxonomy
#' @aliases sampleData
#' @rdname accessor-functions
#' @export
#' @examples
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' x <- read_biom(biom_file)
#' header(x)
#' abundance(x)
#' taxonomy(x)
#' sampleData(x)
header <- function(x){
	# Store header as everything up to the "rows" key.
	# Protect against missing header keys causing an error
	header <- list() # Initialize header
	try(header <- x[1:(charmatch("rows", names(x)) - 1)], TRUE)
	if( identical(header, list()) ){
		warning("Header information missing from this file.")
	}
}
################################################################################
#' @export
#' @aliases header
#' @aliases abundance
#' @rdname accessor-functions
#' @importFrom plyr laply
abundance <- function(x, parallel=FALSE){
	# Check if sparse. Must parse differently than dense
	if( x$matrix_type == "sparse" ){
		otumat <- Matrix(0, nrow=x$shape[1], ncol=x$shape[2])
		# Loop through each sparse line and assign to relevant position in otumat.
		for( i in x$data ){
			otumat[(i[1]+1), (i[2]+1)] <- i[3]
		}
	} else if( x$matrix_type == "dense" ){ 
		# parse the dense matrix instead using plyr's laply
		otumat <- laply(x$data, function(i) i, .parallel=parallel)
	} 
	
	# Get row (OTU) and col (sample) names
	rownames(otumat) <- sapply(x$rows, function(i) i$id )
	colnames(otumat) <- sapply(x$columns, function(i) i$id )
	
	# Instantiates a "Matrix" daughter class, usually sparse,
	# but precise daughter class is chosen dynamically based
	# on properties of the data (if actually dense, it stays dense).
	return(Matrix(otumat))
}
################################################################################
#' @export
#' @aliases header
#' @aliases taxonomy
#' @rdname accessor-functions
#' @importFrom plyr ldply
taxonomy <- function(x, parallel=FALSE){
	# Need to check if taxonomy information is empty (minimal biom file)
	if(  all( sapply(sapply(x$rows, function(i) i$metadata), is.null) )  ){
		taxdf <- NULL
	} else {
		taxdf <- ldply(x$rows, function(i) i$metadata$taxonomy, .parallel=parallel)
		# Add rownames to taxonomy data.frame.
		rownames(taxdf) <- sapply(x$rows, function(i) i$id )
	}
	return(taxdf)
}
################################################################################
#' @export
#' @aliases header
#' @aliases sampleData
#' @rdname accessor-functions
#' @importFrom plyr ldply
sampleData <- function(x, parallel=FALSE){
	# Sample Data ("columns" in biom)
	if(  all( sapply(sapply(x$columns, function(i) i$metadata), is.null) )  ){
		# If there is no metadata (all NULL),
		# then set samdata to NULL, representing empty.
		samdata <- NULL
	} else {
		# Otherwise, parse it.
		samdata <- ldply(x$columns, function(i){
			if( class(i$metadata) == "list"){
				return(i$metadata[[1]])
			} else {
				return(i$metadata)				
			}
		}, .parallel=parallel)
		rownames(samdata) <- sapply(x$columns, function(i) i$id)
	}
}
################################################################################
################################################################################