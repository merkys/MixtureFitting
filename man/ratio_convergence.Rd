\name{ratio_convergence}
\alias{ratio_convergence}
\docType{methods}
\title{Ratio Convergence Check.}
\description{
    Compare two values to tell whether an optimization process has converged.
    The absolute difference between values of two iteractions is divided by the value of previous iteration and compared to the epsilon value.
}
\usage{
    ratio_convergence( p_now, p_prev, epsilon = 1e-6 )
}
\arguments{
    \item{p_now}{function value of \emph{i}-th iteration.}
    \item{p_prev}{function value of \emph{i-1}-th iteration.}
    \item{epsilon}{convergence criterion}
}
\value{TRUE if deemed to have converged, FALSE otherwise}
\author{Andrius Merkys}
