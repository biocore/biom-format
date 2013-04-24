\docType{methods}
\name{header}
\alias{header}
\alias{header,biom-method}
\title{Extract the header from a \code{\link{biom-class}} object as a list.}
\usage{
  header(x)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}
}
\value{
  A list containing the header data. That is, all the
  required elements that are not the main data or index
  metadata.
}
\description{
  Extract the header from a \code{\link{biom-class}} object
  as a list.
}
\examples{
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
x = read_biom(biom_file)
header(x)
}

