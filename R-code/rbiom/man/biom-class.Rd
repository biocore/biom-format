\name{biom-class}
\alias{biom-class}
\title{A native R-representation of the biom-format.}
\description{
  This class inherits from the \code{\link{list-class}},
  with validity checks specific to the definition to the
  biom-format. For documentation of the text file format
  see \url{http://biom-format.org/}.
}
\examples{
#
# # # import with default parameters, specify a file
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
x <- read_biom(biom_file)
show(x)
print(x)
header(x)
biom_table(x)
observ_meta(x)
sample_meta(x)
}
\seealso{
  The constructor, \code{\link{biom}}

  Accessor functions: \code{\link{accessors}}.
}

