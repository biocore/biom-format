\docType{methods}
\name{show}
\alias{show}
\alias{show,biom-method}
\title{Method extensions to show for biom objects.}
\description{
  See the general documentation of
  \code{\link[methods]{show}} method for expected behavior.
}
\examples{
# # # import with default parameters, specify a file
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
(x = read_biom(biom_file) )
show(x)
}
\seealso{
  \code{\link[methods]{show}}
}

