% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator.R
\name{cMLEimat}
\alias{cMLEimat}
\title{Internal function: cMLEimat}
\usage{
cMLEimat(Fk, data, s, wSave = FALSE, S = NULL, onlylogLike = !wSave)
}
\arguments{
\item{Fk}{An \emph{n} by \emph{K} matrix of basis function values with
each column being a basis function taken values at \code{loc}.}

\item{data}{An \emph{n} by \emph{T} data matrix (NA allowed) with
\eqn{z[t]} as the \emph{t}-th column.}

\item{s}{A positive numeric.}

\item{wSave}{A logic.}

\item{S}{An \emph{n} by \emph{n} matrix}

\item{onlylogLike}{A logic.}
}
\value{
A list.
}
\description{
Internal function: cMLEimat
}
\keyword{internal.}
