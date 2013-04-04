\name{read_biom}
\alias{read_biom}
\title{Read a biom-format file, returning a \code{biom-class}.}
\usage{
  read_biom(biom_file, version=1.0)
}
\arguments{
  \item{biom_file}{(Required). A character string
  indicating the file location of the biom formatted file.
  This is a JSON formatted file specific to biological
  datasets. The format is formally defined at
  \href{http://biom-format.org/documentation/biom_format.html}{the
  biom-format definition}

  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}}

  \item{version}{(Optional). Numeric. The expected version
  number of the file-format. As the biom format evolves,
  version-specific importers will be available by adjusting
  the version value. Default is \code{1.0}. Not
  implemented. Has no effect (yet).}
}
\value{
  An instance of the \code{biom-class}.
}
\description{
  Import the data from a biom-format file into R,
  represented as an instance of the
  \code{\link{biom-class}}; essentially a
  \code{\link{list}} with special constraints that map to
  \href{http://biom-format.org/documentation/biom_format.html}{the
  biom-format definition}.
}
\details{
  The BIOM file format (canonically pronounced biome) is
  designed to be a general-use format for representing
  biological sample by observation contingency tables. BIOM
  is a recognized standard for the
  \href{http://www.earthmicrobiome.org/}{Earth Microbiome
  Project} and is a \href{http://gensc.org/}{Genomics
  Standards Consortium} candidate project. Please see
  \href{http://biom-format.org/}{the biom-format home page}
  for more details.
}
\examples{
# # # import with default parameters, specify a file
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "rbiom")
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
# mybiomfile <- "path/to/my/biomfile.biom"
# read_biom(mybiomfile)
}
\references{
  \url{http://www.qiime.org/svn_documentation/documentation/biom_format.html}
}
\seealso{
  The \code{\link{biom}} constructor function.

  Accessor functions like \code{\link{header}}.
}

