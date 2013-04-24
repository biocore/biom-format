\name{sample_metadata}
\alias{sample_metadata}
\alias{sample_metadata,biom,character-method}
\alias{sample_metadata,biom,missing-method}
\alias{sample_metadata,biom,numeric-method}
\title{Access meta data from \code{\link{biom-class}}.}
\usage{
  sample_metadata(x, columns, parallel=FALSE)
}
\arguments{
  \item{x}{(Required). An instance of the
  \code{\link{biom-class}}.}

  \item{columns}{(Optional). The subset of column indices
  described in the returned object. For large datasets,
  specifying the column subset here, rather than after
  creating the whole matrix first, can improve
  speed/efficiency. Can be vector of index numbers
  (\code{\link{numeric-class}}) or index names
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
  \code{\link{data.frame}} (if possible, or a list) with
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
sample_metadata(x1)
sample_metadata(x2)
sample_metadata(x3)
sample_metadata(x3, 1:4)
sample_metadata(x4)
sample_metadata(x5)
sample_metadata(x6)
}

