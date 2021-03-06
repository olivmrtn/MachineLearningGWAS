% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Documentation.R, R/methods-MLGWAS.R
\docType{methods}
\name{plotPerf}
\alias{compareFeatimp}
\alias{plotErrors}
\alias{plotErrors,MLGWAS-method}
\alias{plotFeatimp}
\alias{plotFeatimp,MLGWAS-method}
\alias{plotMLGWAS}
\alias{plotPerf}
\alias{plotPerf,MLGWAS-method}
\title{MLGWAS Plots}
\usage{
plotPerf(x, measure, group)

plotFeatimp(x, sort_by_classifiers, features, max_features = 20)

plotErrors(x)

\S4method{plotPerf}{MLGWAS}(x, measure, group)

\S4method{plotErrors}{MLGWAS}(x)

\S4method{plotFeatimp}{MLGWAS}(x, sort_by_classifiers, features,
  max_features = 20)
}
\arguments{
\item{x}{\code{\link{MLGWAS}} object.}

\item{measure}{\code{character} specifying performance measure to be plotted.
Can be any column name in \code{performance(x)} that is not an identifier.}

\item{group}{Facultative \code{list} containing \code{ClassifierID}s.
\code{ClassifierID}s in the same list entry will be plotted side by side.
By default, all \code{ClassifierID} are plotted separately.}

\item{sort_by_classifiers}{Classifier identifiers
\code{classifierNames(x)} to be used to sort feature importance.}

\item{features}{Subset of \code{FeatureID}s.}

\item{max_features}{\code{integer} specifying maximum number of features in plot.}
}
\description{
Plots performance, feature importance and stability, or errors of
\code{\link{MLGWAS}} objects.
}
\seealso{
\code{\link[ggplot2]{ggplot}}
}

