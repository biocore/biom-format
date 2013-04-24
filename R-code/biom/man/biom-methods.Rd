\docType{methods}
\name{biom}
\alias{biom}
\alias{biom,list-method}
\title{Build and return an instance of the biom-class.}
\usage{
  biom(x)
}
\arguments{
  \item{x}{(REQUIRED). A named list conforming to
  conventions arising from the \code{\link{fromJSON}}
  function reading a biom-format file with default
  settings. See \code{\link{read_biom}} for more details
  about data import and \code{\link{biom-class}} for more
  details about accessor functions that extract R-friendly
  subsets of the data and metadata stored in \code{x}.}
}
\value{
  An instance of the \code{\link{biom-class}}.
}
\description{
  This is for instantiating a biom object within R
  (\code{\link{biom-class}}), and assumes relevant data is
  already available in R. This is different than reading a
  biom file into R. If you are instead interested in
  importing a biom file into R, you should use the
  \code{\link{read_biom}} function. This function is made
  available (exported) so that advanced-users/developers
  can easily represent analogous data in this structure if
  needed. However, most users are expected to instead rely
  on the \code{\link{read_biom}} function for data import,
  followed by accessor functions that extract R-friendly
  subsets of the data stored in the biom-format derived
  list.
}
\details{
  \code{biom()} is a constructor method. This is the main
  method suggested for constructing an experiment-level
  (\code{\link{biom-class}}) object from its component
  data.
}
\examples{
#
# import with default parameters, specify a file
biom_file = system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
x = read_biom(biom_file)
show(x)
print(x)
header(x)
biom_data(x)
biom_shape(x)
nrow(x)
ncol(x)
observation_metadata(x)
sample_metadata(x)
}
\seealso{
  The \code{\link{read_biom}} import function.

  Accessor functions like \code{\link{header}}.
}

