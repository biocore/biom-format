\name{BIOM-class}
\alias{BIOM-class}
\title{A native R-representation of the BIOM format.}
\description{
  For documentation of the text file format see
  \url{http://biom-format.org/}.
}
\details{
  slots: \describe{ \item{otuTable}{a single object of
  class otuTable.} \item{samData}{ a single object of class
  sampleData.} \item{taxTab}{ a single object of class
  taxonomyTable.} \item{tre}{ a single object of class
  phylo, from the package ape} }
}
\seealso{
  The constructor, \code{\link{BIOM}}

  The accessors, \code{\link{header}},
  \code{\link{abundance}}, \code{\link{taxonomy}},
  \code{\link{sampleData}}, \code{\link{tree}}.
}

