################################################################################
#' An S4 copy of the main phylogenetic tree class from the ape package.
#'
#' See the \code{\link[ape]{ape}} package for details about this type of
#' representation of a phylogenetic tree. It is used throught ape.
#'
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#'
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
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
#' @keywords internal
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
################################################################################
# Eventually, class definitions should go here. Most expected is a class
# definition for the "BIOM" class that represents generally the BIOM format
# structure within R.
################################################################################
#' A native R-representation of the BIOM format.
#'
#' For documentation of the format see \url{http://biom-format.org/}. 
#' 
#' slots:
#' \describe{
#'    \item{otuTable}{a single object of class otuTable.}
#'    \item{samData}{ a single object of class sampleData.}
#'    \item{taxTab}{ a single object of class taxonomyTable.}
#'    \item{tre}{ a single object of class phylo, from the package ape}
#' }
#'
#' @seealso The constructor, \code{\link{phyloseq}}, 
#'  the merger \code{\link{merge_phyloseq}}, and also the component 
#'  constructor/accessors \code{\link{otuTable}}, \code{\link{sampleData}},
#'  \code{\link{taxTab}}, and \code{\link{tre}}.
#' 
#' @name BIOM-class
#' @rdname BIOM-class
#' @exportClass BIOM
setClass("BIOM",
	representation(
	header     = "list",
	abundance  = "Matrix",
	taxonomy   = "matrixOrNULL",
	sampleData = "dataFrameOrNull",
	tree       = "phyloOrNULL")
)
################################################################################