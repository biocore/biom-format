\name{biom_data}
\alias{biom_data}
\alias{biom_data,biom,ANY,character-method}
\alias{biom_data,biom,character,ANY-method}
\alias{biom_data,biom,missing,missing-method}
\alias{biom_data,biom,missing,numeric-method}
\alias{biom_data,biom,numeric,missing-method}
\alias{biom_data,biom,numeric,numeric-method}
\title{Access main data observation matrix data from \code{\link{biom-class}}.}
\usage{
  biom_data(x, rows, columns, parallel=FALSE)
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

  \item{parallel}{(Optional). Logical. Whether to perform
  the accession parsing using a parallel-computing backend
  supported by the \code{\link{plyr-package}} via the
  \code{\link[foreach]{foreach-package}}. Note: At the
  moment, the header accessor does not need nor does it
  support parallel-computed parsing.}
}
\value{
  A matrix containing the main observation data, with index
  names. The type of data (numeric or character) will
  depend on the results of
  \code{\link{matrix_element_type}(x)}. The class of the
  matrix returned will depend on the sparsity of the data,
  and whether it has numeric or character data. For now,
  only numeric data can be stored in a
  \code{\link{Matrix-class}}, which will be stored
  sparsely, if possible. Character data will be returned as
  a vanilla \code{\link{matrix-class}}.
}
\description{
  Retrieve and organize main data from
  \code{\link{biom-class}}, represented as a matrix with
  index names.
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
# Extract the data matrices
biom_data(x1)
biom_data(x2)
biom_data(x3)
biom_data(x4)
biom_data(x5)
biom_data(x6)
}

