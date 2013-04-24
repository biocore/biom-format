\docType{methods}
\name{nrow}
\alias{nrow}
\alias{nrow,biom-method}
\title{Method extensions to \code{\link[base]{nrow}}
for \code{\link{biom-class}} objects.}
\usage{
  nrow(x)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}
}
\value{
  The number of rows in \code{x}. A length 1
  \code{\link{integer-class}}.
}
\description{
  See the general documentation of \code{\link[base]{nrow}}
  method for expected behavior.
}
\examples{
# # # import with default parameters, specify a file
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
(x = read_biom(biom_file) )
nrow(x)
}
\seealso{
  \code{\link{ncol}}

  \code{\link[base]{nrow}}

  \code{\link{biom_shape}}
}

