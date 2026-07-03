% Generated manually for legacy-compatible package cleanup.
\name{getNextMaxStats}
\alias{getNextMaxStats}
\title{Get the next-ranked MCI statistic}
\usage{
getNextMaxStats(membersL, idL = NULL, whoisnext, which.next = 2)
}
\arguments{
\item{membersL}{A list of per-state numeric vectors, typically the \code{MCI} element returned by \code{\link{getMCI}}.}

\item{idL}{A list of ranked module indices, usually the \code{idx} element returned by \code{\link{getMaxMCImember}}.}

\item{whoisnext}{Character vector naming the states to evaluate.}

\item{which.next}{Rank to extract from \code{idL}; default is the second-ranked module.}
}
\value{
A named numeric vector of selected module scores.
}
\description{
Retrieves lower-ranked module scores for states with multiple candidate modules.
}
