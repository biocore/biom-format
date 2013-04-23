\docType{methods}
\name{biom_shape}
\alias{biom_shape}
\alias{biom_shape,biom-method}
\title{The matrix dimensions
of a \code{\link{biom-class}} object.}
\usage{
  biom_shape(x)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}
}
\value{
  A length two \code{\link{integer-class}} vector
  indicating the \code{\link{nrow}} and \code{\link{ncol}}
  of the main data matrix stored in \code{x}.
}
\description{
  The matrix dimensions of a \code{\link{biom-class}}
  object.
}
\examples{
# # # import with default parameters, specify a file
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
(x = read_biom(biom_file) )
biom_shape(x)
}
\seealso{
  \code{\link{biom-class}}
}

