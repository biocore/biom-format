\name{header}
\alias{accessors}
\alias{biom_table}
\alias{biomclass}
\alias{biomnull}
\alias{biomshape}
\alias{header}
\alias{observ_meta}
\alias{sample_meta}
\title{Accessors for the \code{\link{biom-class}}.}
\usage{
  header(x)

  biomshape(x)

  biomclass(x)

  biom_table(x, parallel = FALSE)

  observ_meta(x, key = "taxonomy", parallel = FALSE)

  sample_meta(x, parallel = FALSE)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}

  \item{key}{(Optional). Character string. The key for the
  metadata type that you are attempting to access. Default
  value depends upon whether observation or sample
  metadata. This argument only applies to metadat
  accessors.}

  \item{parallel}{(Optional). Logical. Whether to perform
  the accession parsing using a parallel-computing backend
  supported by the \code{\link{plyr-package}} via the
  \code{\link[foreach]{foreach-package}}. Note: At the
  moment, the header accessor does not need nor does it
  support parallel-computed parsing.}
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
biomnull()
}

