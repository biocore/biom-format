\name{header}
\alias{abundance}
\alias{accessors}
\alias{header}
\alias{sampleData}
\alias{taxonomy}
\alias{tree}
\title{Accessors for the \code{\link{biom-class}}.}
\usage{
  header(biom, errorIfNULL=TRUE)

  abundance(biom, errorIfNULL = TRUE)

  taxonomy(biom, errorIfNULL = TRUE)

  sampleData(biom, errorIfNULL = TRUE)

  tree(biom, errorIfNULL = TRUE)
}
\arguments{
  \item{biom}{(Required). \code{\link{biom-class}}.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{FALSE}.}
}
\description{
  Convenience functions for accessing particular components
  of the \code{\link{biom-class}}.
}
\details{
  Some features of the biom-format can be essentially
  empty, and are represented as \code{NULL} in R once
  imported. By default, these accessors will throw an error
  if a \code{NULL} is encountered, under the assumption
  that you are using this accessor because you expect the
  component to be non-empty. Simply setting the second
  argument to \code{FALSE} relieves this and allows
  \code{NULL} to be returned silently.

  Internally these acessors wrap the not-exported
  \code{\link{access}} function. This helps ensure
  consistent behavior from accessors.
}
\examples{
## Examples here.
}

