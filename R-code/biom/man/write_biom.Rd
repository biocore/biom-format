\name{write_biom}
\alias{write_biom}
\title{Write a biom-format file, returning a \code{biom-class}.}
\usage{
  write_biom(x, biom_file)
}
\arguments{
  \item{x}{(Required). A biom object that is going to be
  written to file as a proper biom formatted file, adhering
  to
  \href{http://biom-format.org/documentation/biom_format.html}{the
  biom-format definition}.}

  \item{biom_file}{(Required). A character string
  indicating the file location of the biom formatted file.
  This is a JSON formatted file specific to biological
  datasets. The format is formally defined at
  \href{http://biom-format.org/documentation/biom_format.html}{the
  biom-format definition}}
}
\value{
  Nothing. The first argument, \code{x}, is written to a
  file.
}
\description{
  Write a biom-format file, returning a \code{biom-class}.
}
\examples{
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biom")
x = read_biom(biom_file)
outfile = tempfile()
write_biom(x, outfile)
y = read_biom(outfile)
identical(x, y)
}
\references{
  \url{http://biom-format.org/}
}
\seealso{
  Accessor functions like \code{\link{header}}.
  \href{http://biom-format.org/}{the biom-format home page}
}

