\name{read_biom}
\alias{read_biom}
\title{Read a biom-format file, returning a \code{biom-class}.}
\usage{
  read_biom(biom_file)
}
\arguments{
  \item{biom_file}{(Required). A character string
  indicating the file location of the biom formatted file.
  This is a JSON formatted file specific to biological
  datasets. The format is formally defined at
  \href{http://biom-format.org/documentation/biom_format.html}{the
  biom-format definition}}
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

  It is tempting to include an argument identifying the
  biom-format version number of the data file being
  imported. However, the biom-format version number is a
  required field in the biom-format definition. Rather than
  duplicate this formal specification and allow the
  possibility of a conflict, the version number of the biom
  format will be referred to only by the "format" field in
  the biom formatted data, or its representation in R.
}
\examples{
# # # import with default parameters, specify a file
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
biom_file
read_biom(biom_file)
biom_file <- system.file("extdata", "min_sparse_otu_table.biom", package = "biom")
biom_file
read_biom(biom_file)
## The previous examples use system.file() because of constraints in specifying a fixed
##   path within a reproducible example in a package.
## In practice, however, you can simply provide "hard-link"
## character string path to your file:
# mybiomfile <- "path/to/my/biomfile.biom"
# read_biom(mybiomfile)
}
\references{
  \url{http://biom-format.org/}
}
\seealso{
  The \code{\link{biom}} constructor function.

  Accessor functions like \code{\link{header}}.
}

