################################################################################
# Eventually, class definitions should go here. Most expected is a class
# definition for the "biom" class that represents generally the biom format
# structure within R.
################################################################################
#' A native R-representation of the biom-format.
#'
#' This class inherits from the \code{\link{list-class}},
#' with validity checks specific to the definition to the biom-format.
#' For documentation of the text file format see \url{http://biom-format.org/}. 
#' 
#' @seealso
#' The constructor, \code{\link{biom}}
#' 
#' Accessor functions: \code{\link{accessors}}.
#' 
#' @name biom-class
#' @rdname biom-class
#' @exportClass biom
#'
#' @examples #
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
#' x <- read_biom(biom_file)
#' show(x)
#' print(x)
#' header(x)
#' biom_table(x)
#' observ_meta(x)
#' sample_meta(x)
setClass("biom", contains="list")
################################################################################