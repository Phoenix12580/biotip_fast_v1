% Generated manually for legacy-compatible package cleanup.
\name{getTopMCI}
\alias{getTopMCI}
\title{Get the top MCI score per state}
\usage{
getTopMCI(modulesL, MCI1, membersL, min, n = 1)
}
\arguments{
\item{modulesL}{A list of module membership vectors, typically the \code{members} element returned by \code{\link{getMCI}}.}

\item{MCI1}{A list of module-level MCI score vectors, typically the \code{MCI} element returned by \code{\link{getMCI}}.}

\item{membersL}{A list of numeric vectors from which the selected module statistics are extracted.}

\item{min}{Minimum module size passed to \code{\link{getMaxMCImember}}.}

\item{n}{Number of top states to return.}
}
\value{
A named numeric vector of top MCI scores.
}
\description{
Returns the largest MCI scores after filtering modules by minimum size.
}
