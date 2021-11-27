\name{studentt.stat}
\alias{studentt.stat}
\alias{studentt.fun}


\title{(Paired) Student t Statistic}

\usage{
studentt.stat(X, L, var.equal=TRUE, paired=FALSE)
studentt.fun(L, var.equal=TRUE)
}

\description{
  These functions provide a simple interface to compute (paired) Student t statistics
  in the analysis of high-dimensional case-control studies.
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
\item{L}{factor with class labels for the two groups. If only a single label is given then a one-sample t score against 0 is computed.}
  \item{var.equal}{assume equal (default) or unequal variances in each group.}
\item{paired}{compute paired t-score (default is to use unpaired t-score).}
d}

\value{
  The studentt.stat function returns a vector containing the t-statistic for each variable.  It can be specified whether the variances in the two groups are equal. A paired t-score can also be computed.
  
  The studentt.fun function returns a function that computes the t-score
  statistics when applied to a data matrix (useful for simulations).
}

\author{
  Verena Zuber and
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{shrinkt.stat}}, \code{\link{shrinkcat.stat}}.}

\examples{
# load st library 
library("st")

# load Choe et al. (2005) data
data(choedata)
X <- choe2.mat
dim(X) # 6 11475  
L <- choe2.L
L

# L may also contain some real labels
L = c("group 1", "group 1", "group 1", "group 2", "group 2", "group 2")

# student t statistic (equal variances)
score = studentt.stat(X, L)
order(score^2, decreasing=TRUE)[1:10]
# [1] 11068   724  9990 11387 11310  9985  9996 11046    43    50

# the same computed with standard R methods (slower!)
#score2 = apply(X, 2, function(x) t.test(x ~ L, var.equal=TRUE)$statistic)
#sum((score-score2)^2)

# student t statistic (unequal variances)
score = studentt.stat(X, L, var.equal=FALSE)
order(score^2, decreasing=TRUE)[1:10]
# [1] 11068   724  9990 11387 11310  9985  9996 11046    43    50

# the same computed with standard R methods (slower!)
#score2 = apply(X, 2, function(x) t.test(x ~ L, var.equal=FALSE)$statistic)
#sum((score-score2)^2)

# paired student t statistic
score = studentt.stat(X, L, paired=TRUE)
order(score^2, decreasing=TRUE)[1:10]
# [1] 9985  7239  5393 11387 11310  9942 10238  9996 11015 11276

# the same computed with standard R methods (slower!)
#score2 = apply(X, 2, function(x) t.test(x ~ L, paired=TRUE)$statistic)
#sum((score-score2)^2)
}
\keyword{univar}
