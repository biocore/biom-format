\docType{methods}
\name{matrix_element_type}
\alias{matrix_element_type}
\alias{matrix_element_type,biom-method}
\title{Access class of data in the matrix elements
of a \code{\link{biom-class}} object}
\usage{
  matrix_element_type(x)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}
}
\value{
  A \code{\link{character-class}} string indicating the
  class of the data stored in the main observation matrix
  of \code{x}, with expected values \code{"int"},
  \code{"float"}, \code{"unicode"}.
}
\description{
  Access class of data in the matrix elements of a
  \code{\link{biom-class}} object
}
\examples{
# # # import with default parameters, specify a file
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
(x = read_biom(biom_file) )
matrix_element_type(x)
}
\seealso{
  \code{\link{biom-class}}
}

