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
#' @usage biomclass(x)
#' @usage biomshape(x)
#'
#' @param x (Required). An instance of the \code{\link{biom-class}}.
#'
#' @return A list or vector appropriate to the type of data being accessed.
#' \code{header(x)} - returns a list with all the required elements that are
#' not the main data or index metadata;
#' \code{biomclass(x)} - returns a character string indicating the class of the
#' data stored in the main observation matrix, with expected values
#' \code{"int"}, \code{"float"}, \code{"unicode"};
#' \code{biomshape(x)} - a two-element numeric vector indicating the respective
#' row and column dimensions of the observation matrix.
#'
#' @aliases accessors
#' @aliases header
#' @aliases biomclass
#' @aliases biomshape
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
#' Access main data observation matrix data from \code{\link{biom-class}}. 
#' 
#' Retrieve and organize main data from \code{\link{biom-class}},
#' represented as a matrix with index names.
#'
#' @usage biom_table(x, rows, columns, parallel=FALSE)
#'
#' @param x (Required). An instance of the \code{\link{biom-class}}.
#' @param rows (Optional). The subset of row indices described in the
#'  returned object. For large datasets, specifying the row subset here,
#'  rather than after creating the whole matrix first,
#'  can improve speed/efficiency.
#'  Can be vector of index numbers (\code{\link{numeric-class}}) or 
#'  index names (\code{\link{character-class}}).
#' @param columns (Optional). The subset of column indices described in the
#'  returned object. For large datasets, specifying the column subset here,
#'  rather than after creating the whole matrix first,
#'  can improve speed/efficiency.
#'  Can be vector of index numbers (\code{\link{numeric-class}}) or 
#'  index names (\code{\link{character-class}}).
#' @param parallel (Optional). Logical. Whether to perform the accession parsing
#'  using a parallel-computing backend supported by the \code{\link{plyr-package}}
#'  via the \code{\link[foreach]{foreach-package}}. Note: At the moment, the header
#'  accessor does not need nor does it support parallel-computed parsing.
#'  
#'  @return A matrix containing the main observation data, with index names.
#'   The type of data (numeric or character) 
#'   will depend on the results of \code{\link{biomclass}(x)}.
#'   The class of the matrix returned will depend on the sparsity of the data,
#'   and whether it has numeric or character data.
#'   For now, only numeric data can be stored in a \code{\link{Matrix-class}},
#'   which will be stored sparsely, if possible.
#'   Character data will be returned as a vanilla \code{\link{matrix-class}}.
#'  
#' @aliases biom_table
#' @rdname biom_table-methods
#' @export
#' @examples 
#' min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
#' min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
#' rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
#' rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
#' rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "rbiom")
#' rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "rbiom")
#' # Read the biom-format files
#' x1 = read_biom(min_dense_file)
#' x2 = read_biom(min_sparse_file)
#' x3 = read_biom(rich_dense_file)
#' x4 = read_biom(rich_sparse_file)
#' x5 = read_biom(rich_dense_char)
#' x6 = read_biom(rich_sparse_char)
#' # Extract the data matrices
#' biom_table(x1)
#' biom_table(x2)
#' biom_table(x3)
#' biom_table(x4)
#' biom_table(x5)
#' biom_table(x6)
setGeneric("biom_table", function(x, rows, columns, parallel=FALSE){
  standardGeneric("biom_table")
})
# All methods funnel toward signature biom,numeric,numeric
#' @aliases biom_table,biom,missing,missing-method
#' @rdname biom_table-methods
setMethod("biom_table", c("biom", "missing", "missing"), function(x, rows, columns, parallel){
  # Dispatch with full rows and cols
  biom_table(x, 1:biomshape(x)["nrow"], 1:biomshape(x)["ncol"], parallel)
})
#' @aliases biom_table,biom,character,ANY-method
#' @rdname biom_table-methods
setMethod("biom_table", c("biom", "character", "ANY"), function(x, rows, columns, parallel){
  rows = which(sapply(x$rows, function(s) s$id) %in% rows)
  # Dispatch with specified numeric rows and pass cols
  biom_table(x, rows, columns)
})
#' @aliases biom_table,biom,ANY,character-method
#' @rdname biom_table-methods
setMethod("biom_table", c("biom", "ANY", "character"), function(x, rows, columns, parallel){
  columns = which(sapply(x$columns, function(s) s$id) %in% columns)
  # Dispatch with specified numeric columns and pass rows
  biom_table(x, rows, columns)
})
#' @aliases biom_table,biom,numeric,missing-method
#' @rdname biom_table-methods
setMethod("biom_table", c("biom", "numeric", "missing"), function(x, rows, columns, parallel){
  # Dispatch with specified rows and full cols
  biom_table(x, rows, 1:biomshape(x)["ncol"], parallel)
})
#' @aliases biom_table,biom,missing,numeric-method
#' @rdname biom_table-methods
setMethod("biom_table", c("biom", "missing", "numeric"), function(x, rows, columns, parallel){
  # Dispatch with full rows and specified cols
  biom_table(x, 1:biomshape(x)["nrow"], columns, parallel)
})
#' @aliases biom_table,biom,numeric,numeric-method
#' @rdname biom_table-methods
#' @import Matrix
#' @importFrom plyr d_ply
#' @importFrom plyr ldply
#' @importFrom plyr laply
setMethod("biom_table", c("biom", "numeric", "numeric"), function(x, rows, columns, parallel){
  if( identical(length(rows), 0) ){
    stop("argument `rows` must have non-zero length.")
  }
  if( identical(length(columns), 0) ){
    stop("argument `columns` must have non-zero length.")
  }    
  # Begin matrix section
  if( identical(x$matrix_type, "dense") ){
    # Begin dense section
    # If matrix is stored as dense, create "vanilla" R matrix, m
    m = laply(x$data[rows], function(i) i[columns], .parallel=parallel) 
    if( length(rows) > 1L & length(columns) > 1L & biomclass(x) %in% c("int", "float")){
      # If either dimension is length-one, don't call coerce to "Matrix"
      # Note that laply() does still work in this case.
      # If both dimension lengths > 1 & data is numeric, 
      # attempt to coerce to Matrix-inherited class,
      # Mainly because it might still be sparse and this is a good way
      # to handle it in R. 
      m = Matrix(m)
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
  if( identical(length(rows), 1L) | identical(length(columns), 1L) ){
    # If either dimension is length-one
    # Try naming by colnames first, then rownames
    if( identical(length(rows), 1L) ){
      names(m) <- sapply(x$columns[columns], function(i) i$id )
    } else {
      names(m) <- sapply(x$rows[rows], function(i) i$id )
    }
  } else {
    # Else, both dimensions are longer than 1,
    # can assume is a matrix and assign names to both dimensions
    rownames(m) <- sapply(x$rows[rows], function(i) i$id )
    colnames(m) <- sapply(x$columns[columns], function(i) i$id )
  }
  return(m)
})
################################################################################
#' Access meta data from \code{\link{biom-class}}. 
#' 
#' Retrieve and organize meta data from \code{\link{biom-class}},
#' represented as a \code{\link{data.frame}} (if possible, or a list) 
#' with proper index names.
#'
#' @usage sample_meta(x, columns, parallel=FALSE)
#'
#' @param x (Required). An instance of the \code{\link{biom-class}}.
#' @param columns (Optional). The subset of column indices described in the
#'  returned object. For large datasets, specifying the column subset here,
#'  rather than after creating the whole matrix first,
#'  can improve speed/efficiency.
#'  Can be vector of index numbers (\code{\link{numeric-class}}) or 
#'  index names (\code{\link{character-class}}).
#' @param parallel (Optional). Logical. Whether to perform the accession parsing
#'  using a parallel-computing backend supported by the \code{\link{plyr-package}}
#'  via the \code{\link[foreach]{foreach-package}}. 
#'  
#' @return A \code{\link{data.frame}} or \code{\link{list}} containing 
#'  the meta data, with index names. The precise form of the object returned
#'  depends on the metadata stored in \code{x}. A \code{data.frame} is
#'  created if possible.
#'  
#' @aliases sample_meta
#' @rdname sample_meta-methods
#' @export
#' @examples 
#' min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
#' min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
#' rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
#' rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
#' rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "rbiom")
#' rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "rbiom")
#' # Read the biom-format files
#' x1 = read_biom(min_dense_file)
#' x2 = read_biom(min_sparse_file)
#' x3 = read_biom(rich_dense_file)
#' x4 = read_biom(rich_sparse_file)
#' x5 = read_biom(rich_dense_char)
#' x6 = read_biom(rich_sparse_char)
#' # Extract metadata
#' sample_meta(x1)
#' sample_meta(x2)
#' sample_meta(x3)
#' sample_meta(x3, 1:4)
#' sample_meta(x4)
#' sample_meta(x5)
#' sample_meta(x6)
setGeneric("sample_meta", function(x, columns, parallel=FALSE){
	standardGeneric("sample_meta")
})
# All methods funnel toward signature biom,numeric
#' @aliases sample_meta,biom,missing-method
#' @rdname sample_meta-methods
setMethod("sample_meta", c("biom", "missing"), function(x, columns, parallel=FALSE){
	# Dispatch with full rows and cols
	sample_meta(x, 1:biomshape(x)["ncol"], parallel)
})
#' @aliases sample_meta,biom,character-method
#' @rdname sample_meta-methods
setMethod("sample_meta", c("biom", "character"), function(x, columns, parallel=FALSE){
	columns = which(sapply(x$columns, function(j) j$id) %in% columns)
	if( length(columns)==0 ){
		stop("The column ID names you provided do not match the column IDs in x")
	}
	# Dispatch with specified numeric columns
	sample_meta(x, columns, parallel)
})
#' @rdname sample_meta-methods
#' @aliases sample_meta,biom,numeric-method
setMethod("sample_meta", c("biom", "numeric"), function(x, columns, parallel=FALSE){
	if( any(columns > biomshape(x)["ncol"]) ){
		warning(paste0("column indices ",
									 paste0(columns[columns > biomshape(x)["ncol"]], collapse=" "),
									 " are greater than available columns in data. They were ignored."))
		columns = columns[columns <= biomshape(x)["ncol"]]
	}
	return(extract_metadata(x, "columns", columns, parallel))
})
################################################################################
#' Access observation (row) meta data from \code{\link{biom-class}}. 
#' 
#' Retrieve and organize meta data from \code{\link{biom-class}},
#' represented as a \code{\link{data.frame}} (if possible)
#' or a list, with proper index names.
#'
#' @usage observ_meta(x, rows, parallel=FALSE)
#'
#' @param x (Required). An instance of the \code{\link{biom-class}}.
#' @param rows (Optional). The subset of row indices described in the
#'  returned object. For large datasets, specifying the row subset here,
#'  -- rather than first creating the complete data object --
#'  can improve speed/efficiency.
#'  This parameter can be vector of index numbers (\code{\link{numeric-class}}) or 
#'  index names (\code{\link{character-class}}).
#' @param parallel (Optional). Logical. Whether to perform the accession parsing
#'  using a parallel-computing backend supported by the \code{\link{plyr-package}}
#'  via the \code{\link[foreach]{foreach-package}}. 
#'  
#' @return A \code{\link{data.frame}} or \code{\link{list}} containing 
#'  the meta data, with index names. The precise form of the object returned
#'  depends on the metadata stored in \code{x}. A \code{data.frame} is
#'  created if possible.
#'  
#' @aliases observ_meta
#' @rdname observ_meta-methods
#' @export
#' @examples 
#' min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
#' min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
#' rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
#' rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
#' rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "rbiom")
#' rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "rbiom")
#' # Read the biom-format files
#' x1 = read_biom(min_dense_file)
#' x2 = read_biom(min_sparse_file)
#' x3 = read_biom(rich_dense_file)
#' x4 = read_biom(rich_sparse_file)
#' x5 = read_biom(rich_dense_char)
#' x6 = read_biom(rich_sparse_char)
#' # Extract metadata
#' observ_meta(x1)
#' observ_meta(x2)
#' observ_meta(x3)
#' observ_meta(x3, 2:4)
#' observ_meta(x3, 2)
#' observ_meta(x3, c("GG_OTU_3", "GG_OTU_4", "whoops"))
#' observ_meta(x4)
#' observ_meta(x5)
#' observ_meta(x6)
setGeneric("observ_meta", function(x, rows, parallel=FALSE){
	standardGeneric("observ_meta")
})
# All methods funnel toward signature biom,numeric
#' @aliases observ_meta,biom,missing-method
#' @rdname observ_meta-methods
setMethod("observ_meta", c("biom", "missing"), function(x, rows, parallel=FALSE){
	# Dispatch with full rows and cols
	observ_meta(x, 1:biomshape(x)["nrow"], parallel)
})
#' @aliases observ_meta,biom,character-method
#' @rdname observ_meta-methods
setMethod("observ_meta", c("biom", "character"), function(x, rows, parallel=FALSE){
	rows = which(sapply(x$rows, function(j) j$id) %in% rows)
	if( length(rows)==0 ){
		stop("The row ID names you provided do not match the row IDs in x")
	}
	# Dispatch with specified numeric rows
	observ_meta(x, rows, parallel)
})
#' @rdname observ_meta-methods
#' @aliases observ_meta,biom,numeric-method
setMethod("observ_meta", c("biom", "numeric"), function(x, rows, parallel=FALSE){
	if( any(rows > biomshape(x)["nrow"]) ){
		warning(paste0("Row indices ",
									 paste0(rows[rows > biomshape(x)["nrow"]], collapse=" "),
									 " are greater than available rows in data. They were ignored."))
		rows = rows[rows <= biomshape(x)["nrow"]]
	}
	return(extract_metadata(x, "rows", rows, parallel))
})
################################################################################
# Generic internal function for extracting metadata from either rows or columns
#' @importFrom plyr ldply
#' @importFrom plyr llply
#' @keywords internal
extract_metadata = function(x, indextype, indices, parallel=FALSE){
	# Immediately extract just those index indicated by `index` argument
	metalist = x[[indextype]][indices]
	# Extract metadata elements as a list, for checking dimensions, NULL, etc.
	rx = llply(metalist, function(i) unlist(i$metadata), .parallel=parallel)
	if( all(sapply(rx, is.null)) ){
		# If there is no metadata (all NULL),
		# then set metadata to NULL, representing empty.
		metadata = NULL
	} else {
		# Else, extract names and metadata (both required)
		# Extract names
		metaids = sapply(metalist, function(i) i$id)
		# Test if length of metadata entries is same for all indices.
		rxlengths = sapply(rx, length)
		if( all( rxlengths == rxlengths[1]) ){
			# If so, can parse it as data.frame with ldply
			# return a data.frame with colnames
			metadata = ldply(rx, .parallel=parallel)
			rownames(metadata) <- metaids
		} else {
			# Else, should keep it as a list. But should name the entries
			metadata = rx
			names(metadata) <- metaids
		}
	}
	return(metadata)
}
################################################################################