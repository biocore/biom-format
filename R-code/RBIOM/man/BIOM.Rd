\name{BIOM}
\alias{BIOM}
\title{Build and return an instance of the BIOM-class.}
\usage{
  BIOM(abundance, header=list(), taxonomy=NULL,
  sampleData=NULL, tree=NULL)
}
\arguments{
  \item{abundance}{(REQUIRED). A \code{\link{Matrix}}-class
  object of abundance values. By convention, rows should
  represent taxa and columns different samples.}

  \item{header}{(OPTIONAL). A list or \code{NULL}. The
  header information describing this BIOM file, it's
  format, and other meta-data. Basically all the elements
  above \code{"rows"} in the BIOM format description.
  \url{http://biom-format.org/documentation/format_versions/biom-1.0.html}.
  Default is \code{NULL}.}

  \item{taxonomy}{(OPTIONAL). A (usually character-based)
  \code{link{matrix}} of taxonomic ranks as columns and
  taxa as rows, or \code{NULL}. Default is \code{NULL}.}

  \item{sampleData}{(OPTIONAL). A \code{link{data.frame}}
  representing Sample Data, or \code{NULL}. Default is
  \code{NULL}}

  \item{tree}{(OPTIONAL). A \code{\link[ape]{phylo}}-class
  phylogenetic tree, or \code{NULL}. Default is
  \code{NULL}}
}
\value{
  An instance of the \code{\link{BIOM-class}}.
}
\description{
  This is for instantiating a BIOM object within R
  (\code{\link{BIOM-class}}), and assumes relevant data is
  already available in R. This is different than reading a
  BIOM file into R. If you are instead interested in
  importing a BIOM file into R, you should use the
  \code{\link{read_biom}} function. It is probably worth
  noting that \code{\link{read_biom}} uses this function to
  construct its \code{\link{BIOM-class}} instance after
  parsing the raw data from the BIOM file into R. This
  function is made available (exported) so that other
  users/developers can easily represent analogous data in
  this structure if needed.
}
\details{
  \code{BIOM()} is a constructor method. This is the main
  method suggested for constructing an experiment-level
  (\code{\link{BIOM-class}}) object from its component
  data.
}
\examples{
#
# # # import with default parameters, specify a file
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "RBIOM")
x <- read_biom(biom_file)
show(x)
print(x)
header(x)
abundance(x)
taxonomy(x)
sampleData(x)
tree(x, FALSE)
}
\seealso{
  The \code{\link{read_biom}} import function.

  Accessor functions like \code{\link{header}}.
}

