\name{BIOM-class}
\alias{BIOM-class}
\title{A native R-representation of the BIOM format.}
\description{
  For documentation of the text file format see
  \url{http://biom-format.org/}.
}
\details{
  slots: \describe{ \item{otuTable}{a single object of
  class otuTable.} \item{samData}{ a single object of class
  sampleData.} \item{taxTab}{ a single object of class
  taxonomyTable.} \item{tre}{ a single object of class
  phylo, from the package ape} }
}
\examples{
#
# # # import with default parameters, specify a file
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "RBIOM")
x <- read_biom(biom_file)
show(x)
print(x)
header(x)
abundance(x)
taxonomy(x)
sampleData(x)
tree(x, FALSE)
}
\seealso{
  The constructor, \code{\link{BIOM}}

  The accessors, \code{\link{header}},
  \code{\link{abundance}}, \code{\link{taxonomy}},
  \code{\link{sampleData}}, \code{\link{tree}}.
}

