\name{observ_meta}
\alias{observ_meta}
\alias{sample_meta}
\title{Access meta data from \code{\link{biom-class}}.}
\usage{
  observ_meta(x, rows, key="taxonomy", parallel=FALSE)

  sample_meta(x, columns, parallel = FALSE)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}

  \item{rows}{(Optional). The subset of row indices
  described in the returned object. For large datasets,
  specifying the row subset here, rather than after
  creating the whole matrix first, can improve
  speed/efficiency. Can be vector of index numbers
  (\code{\link{numeric-class}}) or index names
  (\code{\link{character-class}}).}

  \item{columns}{(Optional). The subset of column indices
  described in the returned object. For large datasets,
  specifying the column subset here, rather than after
  creating the whole matrix first, can improve
  speed/efficiency. Can be vector of index numbers
  (\code{\link{numeric-class}}) or index names
  (\code{\link{character-class}}).}

  \item{key}{(Optional). Character string. The key for the
  metadata type that you are attempting to access. Default
  value depends upon whether observation or sample
  metadata. This argument only applies to sample metadata.}

  \item{parallel}{(Optional). Logical. Whether to perform
  the accession parsing using a parallel-computing backend
  supported by the \code{\link{plyr-package}} via the
  \code{\link[foreach]{foreach-package}}. Note: At the
  moment, the header accessor does not need nor does it
  support parallel-computed parsing.}
}
\value{
  A \code{data.frame} containing the meta data, with index
  names.
}
\description{
  Retrieve and organize meta data from
  \code{\link{biom-class}}, represented as a
  \code{\link{data.frame}} with index names.
}
\examples{
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "rbiom")
rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "rbiom")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "rbiom")
rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "rbiom")
rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "rbiom")
# Read the biom-format files
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
x5 = read_biom(rich_dense_char)
x6 = read_biom(rich_sparse_char)
# Extract metadata
observ_meta(x1)
observ_meta(x2)
observ_meta(x3)
observ_meta(x4)
observ_meta(x5)
observ_meta(x6)
sample_meta(x1)
sample_meta(x2)
sample_meta(x3)
sample_meta(x4)
sample_meta(x5)
sample_meta(x6)
}

