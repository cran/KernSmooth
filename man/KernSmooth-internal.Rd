\name{KernSmooth-internal}
\alias{blkest}
\alias{cpblock}
\alias{linbin}
\alias{linbin2D}
\alias{rlbin}
\alias{sdiag}
\alias{sstdiag}
\title{Internal KernSmooth functions}
\description{
 Internal KernSmooth functions
}
\usage{
blkest(x, y, Nval, q)
cpblock(X, Y, Nmax, q)
linbin(X, gpoints, truncate=TRUE)
linbin2D(X, gpoints1, gpoints2)
rlbin(X, Y, gpoints, truncate=TRUE)
sdiag(x, drv=0, degree=1, kernel="normal",
      bandwidth, gridsize=401, bwdisc=25, range.x,
      binned=FALSE, truncate=TRUE)
sstdiag(x, drv=0, degree=1, kernel="normal",
        bandwidth, gridsize=401, bwdisc=25, range.x,
        binned=FALSE, truncate=TRUE)
}
\details{
  These are not to be called by the user.

  \code{blkest} is used to obtain preliminary estimates of
  quantities required for the ``direct plug-in''
  regression bandwidth selector based on
  blocked qth degree polynomial fits.

  \code{cpblock} chooses the number of blocks for the preliminary
  step of a plug-in rule using Mallows' \eqn{C_p}.

  \code{linbin} applies linear binning to a univariate data set.
  
  \code{linbin2D} applies linear binning to a bivariate data set.

  \code{rlbin} applies linear binning to a regression data set.

  \code{sdiag} computes the binned diagonal entries of a smoother
  matrix for local polynomial kernel regression.

  \code{sstdiag} computes the binned diagonal entries of \eqn{SS^T}
  where \eqn{S} is a smoother matrix for local polynomial 
  kernel regression.
}
\keyword{internal}
