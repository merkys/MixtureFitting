\name{rsimplex_start}
\alias{rsimplex_start}
\docType{methods}
\title{Nelder-Mead's Simplex Method for Function Minimization.}
\description{
    Generate initial simplices for simplex().
}
\usage{
    rsimplex_start( seed, n, lower, upper )
}
\arguments{
    \item{seed}{seed for random number generator}
    \item{n}{number of simplices}
    \item{lower}{vector with lower bounds of each dimension}
    \item{upper}{vector with upper bounds of each dimension}
}
\value{A list with \emph{n} simplices.}
\author{Andrius Merkys}
