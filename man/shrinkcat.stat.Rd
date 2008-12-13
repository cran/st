\name{shrinkcat.stat}
\alias{shrinkcat.stat}
\alias{shrinkcat.fun}
\title{Correlation-Adjusted t Score}

\usage{
shrinkcat.stat(X, L, verbose=TRUE)
shrinkcat.fun(L, verbose=TRUE)
}

\description{
  \code{shrinkcat.stat} and \code{shrinkcat.fun} computes a shrinkage
   estimate of the ``correlation-adjusted t score''
   of Zuber and Strimmer (2008).
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{vector with class labels for the two groups. }
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
  The cat (``correlation-adjusted t'') score is the product of the square root of the
  inverse correlation matrix with a vector of t scores.  In Zuber and Strimmer (2008)
  it is shown that the cat score is
  a natural criterion to rank genes according to their ability to seperate two classes
  in the presence of correlation among genes.
  If there is no correlation, the cat score reduces to the usual t score
  (hence in this case the estimate from \code{shrinkcat.stat} equals that from \code{\link{shrinkt.stat}}).
}
\value{
  \code{shrinkcat.stat} returns a vector containing a shrinkage estimate of the
 ``cat score''  for each variable/gene.

  The corresponding \code{shrinkcat.fun} functions return a function that
  computes the cat score when applied to a data matrix
  (this is very useful for simulations).
}

\references{
Zuber, V., and K. Strimmer. 2008.  Gene ranking and differential expression
under correlation. In preparation.
}
\author{
  Verena Zuber 
  and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{shrinkt.stat}}.}

\examples{
# load st library 
library("st")


# load full Khan et al (2001) data set
data(khan2001)

# create data set containing only the RMS and EWS samples
idx = which( khan2001$y == "RMS" | khan2001$y == "EWS")
X = khan2001$x[idx,]
L = factor(khan2001$y[idx])
dim(X)
L


# shrinkage cat statistic
score = shrinkcat.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# [1] 1389 1955  509 1003  246  187 2050 2046  545 1954

# compute q-values and local false discovery rates
library("fdrtool")
fdr.out = fdrtool(as.vector(score))
sum(fdr.out$qval < 0.05)
sum(fdr.out$lfdr < 0.2)


# compared with:

# shrinkage t statistic 
score = shrinkt.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# [1] 1389 1955  187  246 1003 2046 2050  509  545 1799

# student t statistic
score = studentt.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# 1] 1389 1955 1003  187  246 2050 2046 1799  509  545

# difference of means ("Fold Change")
score = diffmean.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# [1]  509  187 1372 1955  246 1954  430 1645  545  129
}
\keyword{univar}
