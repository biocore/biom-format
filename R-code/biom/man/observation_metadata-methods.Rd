\name{observation_metadata}
\alias{observation_metadata}
\alias{observation_metadata,biom,character-method}
\alias{observation_metadata,biom,missing-method}
\alias{observation_metadata,biom,numeric-method}
\title{Access observation (row) meta data from \code{\link{biom-class}}.}
\usage{
  observation_metadata(x, rows, parallel=FALSE)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}

  \item{rows}{(Optional). The subset of row indices
  described in the returned object. For large datasets,
  specifying the row subset here, -- rather than first
  creating the complete data object -- can improve
  speed/efficiency. This parameter can be vector of index
  numbers (\code{\link{numeric-class}}) or index names
  (\code{\link{character-class}}).}

  \item{parallel}{(Optional). Logical. Whether to perform
  the accession parsing using a parallel-computing backend
  supported by the \code{\link{plyr-package}} via the
  \code{\link[foreach]{foreach-package}}.}
}
\value{
  A \code{\link{data.frame}} or \code{\link{list}}
  containing the meta data, with index names. The precise
  form of the object returned depends on the metadata
  stored in \code{x}. A \code{data.frame} is created if
  possible.
}
\description{
  Retrieve and organize meta data from
  \code{\link{biom-class}}, represented as a
  \code{\link{data.frame}} (if possible) or a list, with
  proper index names.
}
\examples{
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "biom")
min_sparse_file  = system.file("extdata", "min_sparse_otu_table.biom", package = "biom")
rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom", package = "biom")
rich_sparse_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
min_dense_file   = system.file("extdata", "min_dense_otu_table.biom", package = "biom")
rich_dense_char  = system.file("extdata", "rich_dense_char.biom", package = "biom")
rich_sparse_char  = system.file("extdata", "rich_sparse_char.biom", package = "biom")
# Read the biom-format files
x1 = read_biom(min_dense_file)
x2 = read_biom(min_sparse_file)
x3 = read_biom(rich_dense_file)
x4 = read_biom(rich_sparse_file)
x5 = read_biom(rich_dense_char)
x6 = read_biom(rich_sparse_char)
# Extract metadata
observation_metadata(x1)
observation_metadata(x2)
observation_metadata(x3)
observation_metadata(x3, 2:4)
observation_metadata(x3, 2)
observation_metadata(x3, c("GG_OTU_3", "GG_OTU_4", "whoops"))
observation_metadata(x4)
observation_metadata(x5)
observation_metadata(x6)
}

