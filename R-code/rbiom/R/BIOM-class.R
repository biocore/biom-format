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
	cat("biom object. \n")
	cat("type:", object$type, "\n")
	cat("matrix_type:", object$matrix_type, "\n")
	cat(biomshape(object)[1], "rows and", biomshape(object)[2], "columns \n")
})
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
header = function(x){
	biomheadkeys = c("id", "format", "format_url", "type", "generated_by", "date",
									 "matrix_type", "matrix_element_type", "shape")
  return(x[biomheadkeys])
}
################################################################################
#' @export
#' @aliases header
#' @aliases biomshape
#' @rdname accessor-functions
biomshape = function(x){
  bsv = as(c(nrow=x$shape[1], ncol=x$shape[2]), "integer")
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
# NEED SPECIAL DOCUMENTATION FOR THIS ONE
# AND NEED TO ADD ROXYGEN HEADERS TO THE METHODS BELOW
# ALSO NEED TO ADD TESTS FOR PROPER BEHAVIOR OF THIS NEW S4 Generic accessor, biom_table
#' @export
#' @aliases header
#' @aliases biom_table
#' @rdname accessor-functions
#' @importFrom plyr laply
#' @import Matrix
setGeneric("biom_table", function(x, rows, columns, parallel=FALSE){
  standardGeneric("biom_table")
})
# All methods funnel toward signature biom,numeric,numeric
setMethod("biom_table", c("biom", "missing", "missing"), function(x, rows, columns, parallel){
  # Dispatch with full rows and cols
  biom_table(x, 1:biomshape(x)["nrow"], 1:biomshape(x)["ncol"], parallel)
})
setMethod("biom_table", c("biom", "character", "ANY"), function(x, rows, columns, parallel){
  rows = which(sapply(x$rows, function(s) s$id) %in% rows)
  # Dispatch with specified numeric rows and pass cols
  biom_table(x, rows, columns)
})
setMethod("biom_table", c("biom", "ANY", "character"), function(x, rows, columns, parallel){
  columns = which(sapply(x$columns, function(s) s$id) %in% columns)
  # Dispatch with specified numeric columns and pass rows
  biom_table(x, rows, columns)
})
setMethod("biom_table", c("biom", "numeric", "missing"), function(x, rows, columns, parallel){
  # Dispatch with specified rows and full cols
  biom_table(x, rows, 1:biomshape(x)["ncol"], parallel)
})
setMethod("biom_table", c("biom", "missing", "numeric"), function(x, rows, columns, parallel){
  # Dispatch with full rows and specified cols
  biom_table(x, 1:biomshape(x)["nrow"], columns, parallel)
})
setMethod("biom_table", c("biom", "numeric", "numeric"), function(x, rows, columns, parallel){
  if( identical(x$matrix_type, "dense") ){
    # Begin dense section
    # If matrix is stored as dense, create "vanilla" R matrix, m
    m = laply(x$data[rows], function(i) i[columns], .parallel=parallel)
    if(biomclass(x) %in% c("int", "float")){
      # If data is numeric, attempt to coerce to Matrix-inherited class
      # Mainly because it might still be sparse and this is a good way
      # to handle it in R. 
      m = Matrix(m)
    } else if(!inherits(m, "matrix")){
      # Else, check that results of laply above are in fact a matrix
      stop("biom: problem parsing dense character matrix")
    }    
  } else {
    # Begin sparse section
    ## Initialize sparse matrix as either Matrix or matrix, depending on data class
    biomshape = biomshape(x)
    if(biomclass(x) %in% c("int", "float")){
      # If data is numeric, initialize with Matrix (numeric data only)
      m = Matrix(0, nrow=biomshape["nrow"], ncol=biomshape["ncol"]) 
    } else {
      # Else, biomclass must be "unicode" for a unicode string.
      # Use a standard R character matrix
      m = matrix(NA_character_, nrow=biomshape["nrow"], ncol=biomshape["ncol"]) 
    }
    # Create an assignment data.frame
    adf = ldply(x$data, function(x){
      data.frame(r=x[[1]], c=x[[2]], data=x[[3]], stringsAsFactors=FALSE)
    })
    # indices start at 0 in biom sparse format
    adf$r <- adf$r + 1L
    adf$c <- adf$c + 1L
    # Subset to just indices that are in both arguments `rows` and `columns`
    adf = subset(adf, r %in% rows & c %in% columns)
    # Add dummy index for plyr .variable and
    # iterate and assign to matrix with double-arrow, m[r, c] <<-
    adf$iter <- 1:nrow(adf)
    d_ply(adf, "iter", function(x) m[x$r, x$c] <<- x$data)
    # Subset this biggest-size m to just `rows` and `columns`
    m = m[rows, columns]
  # End sparse section
  }   
  # Add row and column names
  rownames(m) <- sapply(x$rows[rows], function(i) i$id )
  colnames(m) <- sapply(x$columns[columns], function(i) i$id )
  return(m)
})
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
	if( all(sapply(sapply(x$columns, function(i) i$metadata), is.null)) ){
		# If there is no metadata (all NULL),
		# then set samdata to NULL, representing empty.
		samdata = NULL
	} else {
		# Otherwise, parse it.
		samdata = ldply(x$columns, function(i) i$metadata, .parallel=parallel)
		rownames(samdata) <- sapply(x$columns, function(i) i$id)
	}
	return(samdata)
}
################################################################################
################################################################################