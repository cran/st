\name{shrinkt.stat}
\alias{shrinkt.stat}
\alias{shrinkt.fun}
\title{The Shrinkage t Statistic}

\usage{
shrinkt.stat(X, L, var.equal=TRUE, verbose=TRUE)
shrinkt.fun(L, var.equal=TRUE, verbose=TRUE)
}

\description{
  \code{shrinkt.stat} and \code{shrinkt.fun} compute the ``shrinkage t'' statistic
   of Opgen-Rhein and Strimmer (2006).
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{group indicator vector.  Samples belonging to the first group  are
           assigned a `1', and those belonging to the second group a `2'. }
  \item{var.equal}{assume equal (default) or unequal variances in each group.}
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
  The ``shrinkage t'' statistic is similar to the usual t statistic, with the
  replacement of the sample variances by corresponding shrinkage estimates.
  These are derived in a distribution-free fashion and with little a priori
  assumptions.   Despite its simplicity,
  it can be shown that ``shrinkage t'' procduces highly accurate rankings -
  see Opgen-Rhein and Strimmer (2006).
}
\value{
  \code{shrinkt.stat} returns a vector containing the ``shrinkage t'' 
  statistic for each variable/gene.

  The corresponding \code{shrinkt.fun} functions return a function that
  produces the ``shrinkage t'' statistics when applied to a data matrix
  (this is very useful for simulations).
}

\references{
 Opgen-Rhein, R., and K. Strimmer. 2006. Accurate ranking of 
 differentially expressed genes by a distribution-free shrinkage 
 approach.  
 
 A preprint is available at 
 \url{http://strimmerlab.org/publications/shrinkt2006.pdf}.
}
\author{
  Rainer Opgen-Rhein (\url{http://opgen-rhein.de})
  and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{studentt.stat}},
         \code{\link{diffmean.stat}},
	 \code{\link{efront.stat}},
	 \code{\link{sam.stat}},
         \code{\link{samL1.stat}},
	 \code{\link{modt.stat}}.}

\examples{
# load st library 
library("st")

# load Choe et al. (2005) data
data(choedata)
X <- choe2.mat
dim(X) # 6 11475  
L <- choe2.L
L

# shrinkage t statistic (equal variances)
score = shrinkt.stat(X, L)
order(abs(score), decreasing=TRUE)[1:10]

# [1] 10979 11068    50  1022   724  5762    43  4790 10936  9939
#  lambda.var (variance vector):  0.3882


# shrinkage t statistic (unequal variances)
score = shrinkt.stat(X, L, var.equal=FALSE)
order(abs(score), decreasing=TRUE)[1:10]

# [1] 11068    50 10979   724    43  1022  5762 10936  9939  9769
#  lambda.var (variance vector):  0.3673   0.3362

}
\keyword{univar}