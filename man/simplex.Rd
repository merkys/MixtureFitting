\name{simplex}
\alias{simplex}
\docType{methods}
\title{Nelder-Mead's Simplex Method for Function Minimization.}
\description{
    Nelder-Mead's Simplex Method for Function Minimization.
}
\usage{
    simplex( fn, start, ..., epsilon = 0.000001, alpha = 1,
             gamma = 2, rho = 0.5, delta = 0.5, trace = FALSE )
}
\arguments{
    \item{fn}{minimized function, has to accept the argmin vector as first parameter}
    \item{start}{start vector}
    \item{...}{other parameters passed to the minimized function}
    \item{epsilon}{convergence criterion}
    \item{alpha}{reflection coefficient}
    \item{gamma}{expansion coefficient}
    \item{rho}{contraction coefficient}
    \item{delta}{shrink coefficient}
    \item{trace}{should debug trace be printed?}
}
\value{Vector yielding the minimum value of the minimized function}
\references{
    Nelder, J. A. & Mead, R.
    A Simplex Method For Function Minimization.
    The Computer Journal, 1965, 308-313.

    Users of Wikipedia.
    Nelder-Mead method.
    \url{https://en.wikipedia.org/w/index.php?title=Nelder\%E2\%80\%93Mead_method&oldid=1287347131}
}
\author{Andrius Merkys}
