\name{biom-class}
\alias{biom-class}
\title{The biom format data class.}
\description{
  This class inherits from the \code{\link{list-class}},
  with validity checks specific to the definition to the
  biom-format. Effectively this means the list must have
  certain index names, some elements of which must have a
  specific structure or value. For further details see
  \href{http://biom-format.org/documentation/biom_format.html}{the
  biom-format definition}.
}
\examples{
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

