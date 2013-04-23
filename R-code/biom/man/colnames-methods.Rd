\docType{methods}
\name{colnames}
\alias{colnames}
\alias{colnames,biom-method}
\title{Method extensions to \code{\link[base]{colnames}}
for \code{\link{biom-class}} objects.}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}
}
\value{
  The number of columns in \code{x}. A length 1
  \code{\link{integer-class}}.
}
\description{
  See the general documentation of
  \code{\link[base]{colnames}} method for expected
  behavior.
}
\examples{
# # # import with default parameters, specify a file
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
(x = read_biom(biom_file) )
colnames(x)
}
\seealso{
  \code{\link{nrow}}

  \code{\link[base]{colnames}}

  \code{\link{biom_shape}}
}

