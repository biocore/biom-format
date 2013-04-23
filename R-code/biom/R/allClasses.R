################################################################################
#' The biom format data class.
#'
#' This class inherits from the \code{\link{list-class}},
#' with validity checks specific to the definition to the biom-format.
#' Effectively this means the list must have certain index names,
#' some elements of which must have a specific structure or value.
#' For further details see
#' \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}.
#' Importantly, this means other special properties of lists,
#' like operations with \code{$} and single- or double-square-braces
#' are also supported; as-is the \code{apply}-family function
#' that can operate on lists. 
#' Note that some features of the biom-format can be essentially empty,
#' represented by the string \code{"null"} in the file.
#' These fields are returned as \code{\link{NULL}} when accessed 
#' by an accessor function.
#' 
#' @seealso
#' The constructor, \code{\link{biom}}
#' 
#' Accessor functions:
#' 
#' \code{\link{header}}, 
#' \code{\link{biom_shape}},
#' \code{\link{nrow}},
#' \code{\link{ncol}},
#' \code{\link{matrix_element_type}},
#' \code{\link{biom_data}},
#' \code{\link{observation_metadata}},
#' \code{\link{sample_metadata}}
#' 
#' @name biom-class
#' @rdname biom-class
#' @exportClass biom
#'
#' @examples 
#' biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
#' x = read_biom(biom_file)
#' header(x)
#' biom_shape(x)
#' nrow(x)
#' ncol(x)
#' rownames(x)
#' colnames(x)
#' matrix_element_type(x)
#' biom_data(x)
#' observation_metadata(x)
#' sample_metadata(x)
setClass("biom", contains="list")
################################################################################
