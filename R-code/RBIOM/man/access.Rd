\name{access}
\alias{access}
\title{Universal slot accessor for the \code{\link{BIOM-class}}.}
\usage{
  access(biom, slot, errorIfNULL=FALSE)
}
\arguments{
  \item{biom}{(Required). \code{\link{BIOM-class}}.}

  \item{slot}{(Required). A character string indicating the
  slot (not data class) of the component data type that is
  desired.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{FALSE}.}
}
\value{
  Returns the component object specified by the argument
  \code{slot}. Returns NULL if slot does not exist. Returns
  \code{biom} as-is if it is a component class that already
  matches the slot name.
}
\description{
  This function is used by slot-specific accessors.
}
\keyword{internal}

