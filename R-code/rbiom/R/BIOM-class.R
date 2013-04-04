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
#' @usage biom(x)
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
#' biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' x = read_biom(biom_file)
#' show(x)
#' print(x)
#' header(x)
#' biom_table(x)
#' observ_meta(x)
#' sample_meta(x)
biom = function(x){
	# Some instantiation checks chould go here, or wrap them in validity methods.
	# Depends on how strict the check should be.
	biom = new("biom", x)
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
#' biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' (x = read_biom(biom_file) )
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
#' Some features of the biom-format can be essentially empty,
#' represented by the string \code{"null"} in the file.
#' These fields are returned as \code{NULL} when accessed 
#' by an accessor function.
#'
#' @usage header(x)
#' @usage biom_table(x, parallel=FALSE)
#' @usage observ_meta(x, parallel=FALSE)
#' @usage sample_meta(x, parallel=FALSE)
#' @usage biomnull()
#'
#' @param x (Required). An instance of the \code{\link{biom-class}}.
#' @param key (Optional). Character string. The key for the metadata type
#'  that you are attempting to access. Default value depends upon whether
#'  observation or sample metadata. This argument only applies to 
#'  metadat accessors.
#' @param parallel (Optional). Logical. Whether to perform the accession parsing
#'  using a parallel-computing backend supported by the \code{\link{plyr-package}}
#'  via the \code{\link[foreach]{foreach-package}}. Note: At the moment, the header
#'  accessor does not need nor does it support parallel-computed parsing.
#'
#' @aliases accessors
#' @aliases header
#' @aliases biomnull
#' @aliases biomclass
#' @aliases biomshape
#' @aliases biom_table
#' @aliases observ_meta
#' @aliases sample_meta
#' @rdname accessor-functions
#' @export
#' @examples
#' biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' x = read_biom(biom_file)
#' header(x)
#' biomshape(x)
#' biomclass(x)
#' biom_table(x)
#' observ_meta(x)
#' sample_meta(x)
#' biomnull()
header = function(x){
	# Store header as everything up to the "rows" key.
	# Protect against missing header keys causing an error
	biomheader = list() # Initialize biomheader
	try(biomheader <- x[1:(charmatch("rows", names(x)) - 1)], TRUE)
	if( identical(biomheader, list()) ){
		stop("biom header information missing")
	}
  return(header)
}
################################################################################
#' @export
#' @aliases header
#' @aliases biomnull
#' @rdname accessor-functions
biomnull = function(){
  # No-argument function that defines the NULL/missing value string in biom-format
  return("null")
}
################################################################################
#' @export
#' @aliases header
#' @aliases biomshape
#' @rdname accessor-functions
biomshape = function(x){
  bsv = as.integer(c(nrow=x$shape[1], ncol=x$shape[2]))
  if(!inherits(bsv, "integer")){stop("problem with biom shape value type")}
  if(!identical(length(bsv), 2L)){stop("problem with biom shape value length")}
  if(any(bsv < 0)){stop("problem with biom shape value: negative value")}
  return(bsv)
}
################################################################################
#' @export
#' @aliases header
#' @aliases biomclass
#' @rdname accessor-functions
biomclass = function(x){
  bcl = x$matrix_element_type
  if(!identical(length(bcl), 1L)){
    stop("biom matrix_element_type field should have only 1 element")
  }
  if( !bcl %in% c("int", "float", "unicode") ){
    stop("biom matrix_element_type value is unsupported.")
  }
  return(bcl)
}
################################################################################
#' @export
#' @aliases header
#' @aliases biom_table
#' @rdname accessor-functions
#' @importFrom plyr laply
#' @import Matrix
biom_table = function(x, parallel=FALSE){
	if( x$matrix_type == "sparse" ){
	  # If sparse, must parse accordingly. Dense below.
    biomshape = biomshape(x)
    if(biomclass(x) %in% c("int", "float")){
      # If data is numeric, initialize with Matrix (numeric data only)
      otumat = Matrix(0, nrow=biomshape["nrow"], ncol=biomshape["ncol"]) 
    } else {
      # Else, biomclass must be "unicode" for a unicode string.
      # Use a standard R character matrix
      otumat = matrix(NA_character_, nrow=biomshape["nrow"], ncol=biomshape["ncol"]) 
    }
		# Loop through each sparse line and assign to relevant position in otumat.
		for( i in x$data ){
			otumat[(i[1]+1), (i[2]+1)] <- i[3]
		}
	} else if( x$matrix_type == "dense" ){ 
		# parse the dense matrix instead using plyr's laply
		otumat = laply(x$data, function(i) i, .parallel=parallel)
		if(biomclass(x) %in% c("int", "float")){
		  # If data is numeric, attempt to coerce to Matrix-inherited class
      # Mainly because it might still be sparse and this is a good way
      # to handle it in R. 
		  otumat = Matrix(otumat)
		} else if(!inherits(otumat, "matrix")){
      # Else, check that results of laply above are in fact a matrix
      stop("biom: problem parsing dense character matrix")
		}
	} 
	
	# Define row and col indices
	rownames(otumat) <- sapply(x$rows, function(i) i$id )
	colnames(otumat) <- sapply(x$columns, function(i) i$id )
	
	return(otumat)
}
################################################################################
#' @export
#' @aliases header
#' @aliases observ_meta
#' @rdname accessor-functions
#' @importFrom plyr ldply
observ_meta = function(x, key="taxonomy", parallel=FALSE){
	# Need to check if observ_meta data is empty (minimal biom file)
	if(  all( sapply(sapply(x$rows, function(i) i$metadata), is.null) )  ){
		obsmetadf = NULL
	} else {
		obsmetadf = ldply(x$rows, function(i, key) i$metadata[[key]],
                       key=key, .parallel=parallel)
		# Add rownames to observ_meta data.frame.
		rownames(obsmetadf) <- sapply(x$rows, function(i) i$id )
	}
	return(obsmetadf)
}
################################################################################
#' @export
#' @aliases header
#' @aliases sample_meta
#' @rdname accessor-functions
#' @importFrom plyr ldply
sample_meta = function(x, parallel=FALSE){
	# Sample Data ("columns" in biom)
	if(  all( sapply(sapply(x$columns, function(i) i$metadata), is.null) )  ){
		# If there is no metadata (all NULL),
		# then set samdata to NULL, representing empty.
		samdata = NULL
	} else {
		# Otherwise, parse it.
		samdata = ldply(x$columns, function(i){
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