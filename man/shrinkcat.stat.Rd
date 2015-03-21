\name{shrinkcat.stat}
\alias{shrinkcat.stat}
\alias{shrinkcat.fun}
\title{Correlation-Adjusted t Score (CAT score)}

\usage{
shrinkcat.stat(X, L, lambda, lambda.var, lambda.freqs, var.equal=TRUE, 
   paired=FALSE, verbose=TRUE)
shrinkcat.fun(L, lambda, lambda.var, lambda.freqs, var.equal=TRUE, 
   verbose=TRUE)
}

\description{
  \code{shrinkcat.stat} and \code{shrinkcat.fun} compute a shrinkage
   estimate of the ``correlation-adjusted t score''
   of Zuber and Strimmer (2009).
}

\arguments{
  \item{X}{data matrix. Note that the \emph{columns} correspond to variables (``genes'')
           and the \emph{rows} to samples.}
  \item{L}{factor with class labels for the two groups. If only a single label is given then a one-sample CAT score against 0 is computed. }
  \item{lambda}{Shrinkage intensity for the correlation matrix. If not specified it is 
    estimated from the data. \code{lambda=0} implies no shrinkage
    and \code{lambda=1} complete shrinkage. }
  \item{lambda.var}{Shrinkage intensity for the variances. If not specified it is 
    estimated from the data. \code{lambda.var=0} implies no shrinkage
    and \code{lambda.var=1} complete shrinkage. }
  \item{lambda.freqs}{Shrinkage intensity for the frequencies. If not specified it is 
    estimated from the data. \code{lambda.freqs=0} implies no shrinkage (i.e. empirical frequencies).}

  \item{var.equal}{assume equal (default) or unequal variances in each group.}
  \item{paired}{compute paired CAT score (default is to use unpaired CAT score).}
  \item{verbose}{print out some (more or less useful) information during computation.}
}
\details{
  The CAT (``correlation-adjusted t'') score is the product of the square root of the
  inverse correlation matrix with a vector of t scores. The CAT score thus describes the 
  contribution of each individual feature in separating the two groups, 
  after removing the effect of all other features.

   In Zuber and Strimmer (2009)
  it is shown that the CAT score is
  a natural criterion to rank features in the presence of correlation.
  If there is no correlation, the CAT score reduces to the usual t score
  (hence in this case the estimate from \code{shrinkcat.stat} equals that from \code{\link{shrinkt.stat}}).

  The function \code{\link[sda]{catscore}} implements multi-class CAT scores.
}
\value{
  \code{shrinkcat.stat} returns a vector containing a shrinkage estimate of the
 ``CAT score''  for each variable/gene.

  The corresponding \code{shrinkcat.fun} functions return a function that
  computes the cat score when applied to a data matrix
  (this is very useful for simulations).

 The scale factor in the ''shrinkage CAT'' statistic is computed from the estimated frequencies
 (to use the standard empirical scale factor set  \code{lambda.freqs=0}).
}

\references{
Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
Bioinformatics 25: 2700-2707.
\url{http://dx.doi.org/10.1093/bioinformatics/btp460}.
}
\author{
  Verena Zuber 
  and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link[sda]{catscore}}, \code{\link{shrinkt.stat}}, \code{\link{cst.stat}}, \code{\link{lait.stat}}.}

\examples{
# load st library 
library("st")

# prostate data set
data(singh2002)
X = singh2002$x
L = singh2002$y

dim(X)      # 102 6033 
length(L)   # 102


# shrinkage cat statistic
score = shrinkcat.stat(X, L)
idx = order(score^2, decreasing=TRUE)
idx[1:10]
# 610  364 1720 3647 3375  332 3282 3991 1557  914

# compute q-values and local false discovery rates
library("fdrtool")
fdr.out = fdrtool(as.vector(score))
sum(fdr.out$qval < 0.05)
sum(fdr.out$lfdr < 0.2)


# compared with:

# shrinkage t statistic 
score = shrinkt.stat(X, L)
idx = order(score^2, decreasing=TRUE)
idx[1:10]
# 610 1720 3940  914  364  332 3647 4331  579 1068

# shrinkage CAT score with zero correlation among predictors
# is the same as shrinkage t
score2 = shrinkcat.stat(X, L, lambda=1)
sum((score2-score)^2)


# Student t statistic
score = studentt.stat(X, L)
idx = order(score^2, decreasing=TRUE)
idx[1:10]
# 610 1720  364  332  914 3940 4546 1068  579 4331

# shrinkage CAT score with zero correlation and no shrinkage 
# is the same as student t
score2 = shrinkcat.stat(X, L, lambda=1, lambda.var=0, lambda.freqs=0,
  verbose=FALSE)
sum((score2-score)^2)


# difference of means ("Fold Change")
score = diffmean.stat(X, L)
idx = order(abs(score), decreasing=TRUE)
idx[1:10]
# 735  610  694  298  698  292  739 3940  702  721


## paired CAT score

# we drop two cancer cases to make samples size equal in
# the two groups to allow to compute paired statistic
X = X[1:100,]
L = L[1:100]
sum(L=="cancer") # 50
sum(L=="healthy") # 50

# paired shrinkage CAT score
scat.paired = shrinkcat.stat(X, L, paired=TRUE)

# for zero correlation the paired shrinkage CAT score
# reduces to the paired shrinkage t score
score = shrinkt.stat(X, L, paired=TRUE, verbose=FALSE)
score2 = shrinkcat.stat(X, L, lambda=1, paired=TRUE, verbose=FALSE)
sum((score-score2)^2)

}
\keyword{multivar}
