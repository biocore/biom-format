\name{read_biom}
\alias{read_biom}
\title{Read a BIOM-format file, returning a \code{biom-class}.}
\usage{
  read_biom(biom_file, taxaPrefix=NULL, parallel=FALSE,
  version=1.0)
}
\arguments{
  \item{biom_file}{(Required). A character string
  indicating the file location of the BIOM formatted file.
  This is a JSON formatted file, specific to biological
  datasets, as described in

  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}}

  \item{taxaPrefix}{(Optional). Character string. What
  category of prefix precedes the taxonomic label at each
  taxonomic rank. Currently only ``greengenes'' is a
  supported option, and implies that the first letter
  indicates the taxonomic rank, followed by two underscores
  and then the actual taxonomic assignment at that rank.
  The default value is \code{NULL}, meaning that no prefix
  or rank identifier will be interpreted.}

  \item{parallel}{(Optional). Logical. Wrapper option for
  \code{.parallel} parameter in \code{plyr-package}
  functions. If \code{TRUE}, apply parsing functions in
  parallel using parallel backend declared with (probably)
  \code{foreach}-package and its supporting backend
  packages. This is implemented via plyr, so backend
  options could change depending on how plyr-package
  evolves. A second caveat, plyr-parallelization currently
  works most-cleanly with \code{multicore}-like backends
  (Mac OS X, Unix?), and may throw warnings for SNOW-like
  backends. See the example below for code invoking
  multicore-style backend within the \code{doParallel}
  package.}

  \item{version}{(Optional). Numeric. The expected version
  number of the file. As the BIOM format evolves,
  version-specific importers will be available by adjusting
  the version value. Default is \code{1.0}. Not
  implemented. Has no effect (yet).}
}
\value{
  An instance of the \code{biom-class}.
}
\description{
  New versions of QIIME produce a more-comprehensive and
  formally-defined JSON file format. From the QIIME
  website:
}
\details{
  ``The biom file format (canonically pronounced `biome')
  is designed to be a general-use format for representing
  counts of observations in one or more biological
  samples.''

  \url{http://biom-format.org/}
}
\examples{
# # # import with default parameters, specify a file
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "RBIOM")
read_biom(biom_file)
x <- read_biom(biom_file)
show(x)
print(x)
header(x)
abundance(x)
taxonomy(x)
sampleData(x)
tree(x, FALSE)
## The previous example uses system.file() because of constraints in specifying a fixed
##   path within a reproducible example in a package.
## In practice, however, you should simply provide "hard-link"
## character string path to your file:
# myBIOMfile <- "path/to/my/biomfile.biom"
# read_biom(myBIOMfile)
}
\references{
  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
}
\seealso{
  The \code{\link{BIOM}} constructor function.

  Accessor functions like \code{\link{header}}.
}

