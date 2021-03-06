% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-MLGWAS.R
\docType{methods}
\name{trainClassifiers}
\alias{trainClassifiers}
\alias{trainClassifiers,MDT,data.frame-method}
\title{Calls trainClassifier for different combination of parameters.}
\usage{
trainClassifiers(mdt, params, rbind = TRUE, verbose = TRUE, seed = NULL,
  .export = NULL)

\S4method{trainClassifiers}{MDT,data.frame}(mdt, params, rbind = TRUE,
  verbose = TRUE, seed = NULL, .export = NULL)
}
\arguments{
\item{mdt}{\code{\link{MDT}} object.}

\item{params}{\code{data.frame} of different combinations of parameters.
See \code{\link[base]{expand.grid}} to see how to create such a data.frame.
To specify a parameter of a length bigger than one,
and for which combinations of values are not relevant
(such as \code{phen_vars}), specify it as a list.
For example, \code{expand.grid(p = 0.5, phen_vars = list('Sex', 'Age'))}.
Arguments with value NA will not be passed on to \code{\link{trainClassifier}}.}

\item{rbind}{\code{logical} specifying if results should be bound at the end with
\code{rbind} and return just one \code{\link{MLGWAS}} object, or, if results
should be returned as a \code{list} of \code{\link{MLGWAS}} objects}

\item{verbose}{\code{logical} specifying if function should be run in
verbose mode.}

\item{seed}{\code{integer} specifying random number generation seed.}

\item{.export}{character vector of variables to export.
    This can be useful when accessing a variable that isn't defined in the
    current environment.
    The default value in \code{NULL}.}
}
\value{
An \code{\link{MLGWAS}} object or a \code{list} of \code{\link{MLGWAS}}
objects according to \code{rbind} argument.
}
\description{
Calls trainClassifier for different combination of parameters as specified
by a \code{data.frame}.
}
\seealso{
\code{\link{trainClassifier}}
}

