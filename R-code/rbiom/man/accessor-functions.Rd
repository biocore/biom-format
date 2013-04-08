\name{header}
\alias{accessors}
\alias{biomclass}
\alias{biomshape}
\alias{header}
\title{Accessors for the \code{\link{biom-class}}.}
\usage{
  header(x)

  biomshape(x)

  biomclass(x)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}
}
\value{
  A list or vector appropriate to the type of data being
  accessed. \code{header(x)} - returns a list with all the
  required elements that are not the main data or index
  metadata; \code{biomclass(x)} - returns a character
  string indicating the class of the data stored in the
  main observation matrix, with expected values
  \code{"int"}, \code{"float"}, \code{"unicode"};
  \code{biomshape(x)} - a two-element numeric vector
  indicating the respective row and column dimensions of
  the observation matrix.
}
\description{
  Convenience functions for accessing particular components
  of the \code{\link{biom-class}}.
}
\details{
  Some features of the biom-format can be essentially
  empty, represented by the string \code{"null"} in the
  file. These fields are returned as \code{NULL} when
  accessed by an accessor function.
}
\examples{
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
x = read_biom(biom_file)
header(x)
biomshape(x)
biomclass(x)
biom_table(x)
observ_meta(x)
sample_meta(x)
}

