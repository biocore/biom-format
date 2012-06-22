################################################################################
#' An S4 copy of the main phylogenetic tree class from the ape package.
#'
#' See the \code{\link[ape]{ape}} package for details about this type of
#' representation of a phylogenetic tree. It is used throught ape.
#'
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#'
#' @import ape
#' @name phylo-class
#' @rdname phylo-class
#' @exportClass phylo
setOldClass("phylo")
################################################################################
# Use setClassUnion to define the unholy NULL-data union as a virtual class.
# This is a way of dealing with the expected scenarios in which one or more of
# the component data classes is not available, in which case NULL will be used
# instead.
################################################################################
#' @keywords internal
setClassUnion("dataFrameOrNull", c("data.frame", "NULL"))
#' @keywords internal
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
################################################################################
# Eventually, class definitions should go here. Most expected is a class
# definition for the "BIOM" class that represents generally the BIOM format
# structure within R.
################################################################################
#' A native R-representation of the BIOM format.
#'
#' For documentation of the text file format see \url{http://biom-format.org/}. 
#' 
#' slots:
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{samData}{ a single object of class sampleData.}
#'    \item{taxTab}{ a single object of class taxonomyTable.}
#'    \item{tre}{ a single object of class phylo, from the package ape}
#' }
#'
#' @seealso
#' The constructor, \code{\link{BIOM}}
#' 
#' The accessors, \code{\link{header}}, \code{\link{abundance}},
#'  \code{\link{taxonomy}}, \code{\link{sampleData}}, \code{\link{tree}}.
#' 
#' @import ape
#' @import Matrix
#' @name BIOM-class
#' @rdname BIOM-class
#' @exportClass BIOM
#'
#' @examples #
#' # # # import with default parameters, specify a file
#' biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "RBIOM")
#' x <- read_biom(biom_file)
#' show(x)
#' print(x)
#' header(x)
#' abundance(x)
#' taxonomy(x)
#' sampleData(x)
#' tree(x, FALSE)
setClass("BIOM",
	representation(
	header     = "list",
	abundance  = "Matrix",
	taxonomy   = "dataFrameOrNull",
	sampleData = "dataFrameOrNull",
	tree       = "phyloOrNULL")
)
################################################################################