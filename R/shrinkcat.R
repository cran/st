### shrinkcat.R  (2012-12-02)
###
###    Shrinkage Estimation of Correlation-Adjusted t Statistic
###
### Copyright 2008-2012 Verena Zuber and Korbinian Strimmer
###
###
### This file is part of the `st' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


shrinkcat.stat = function (X, L, lambda, lambda.var, var.equal=TRUE, paired=FALSE, verbose=TRUE)
{
  if (paired)
  {
    X = pvt.pairX(X, L)
    L = rep("paired", dim(X)[1])
  }

  FUN = shrinkcat.fun(L=L, lambda=lambda, lambda.var=lambda.var, 
           var.equal=var.equal, verbose=verbose)
  score = FUN(X)
  
  return( score )
}


shrinkcat.fun = function (L, lambda, lambda.var, var.equal=TRUE, verbose=TRUE)
{
    if (missing(L)) stop("Class labels are missing!")
    if (missing(lambda)) auto.lambda=TRUE
    else auto.lambda=FALSE    
    if (missing(lambda.var)) auto.lambda.var=TRUE
    else auto.lambda.var=FALSE    

    function(X)
    {
      if (auto.lambda==FALSE)
      {
        if(lambda == 1) compute.cor=FALSE
        else compute.cor=TRUE
      }
      else
      {
        compute.cor=TRUE
      }     

      if (auto.lambda.var)
        tmp = centroids(X, L, var.groups=(!var.equal), 
               centered.data=compute.cor, verbose=verbose)
      else
        tmp = centroids(X, L, lambda.var=lambda.var, var.groups=(!var.equal),
               centered.data=compute.cor, verbose=verbose)       

      numClass = length(tmp$samples)
      if(numClass == 1) # one-sample t-score
      {
         diff = tmp$means[,1]
         v1 = tmp$variances[,"(pooled)"]/tmp$samples[1] # variance/n
         v2 = 0
      }
      else if(numClass == 2) # two-sample tscore
      {
        diff = tmp$means[,1]-tmp$means[,2]
        if(var.equal)
        {
          v1 = tmp$variances[,"(pooled)"]/tmp$samples[1] # pooled variance/n1
          v2 = tmp$variances[,"(pooled)"]/tmp$samples[2] # pooled variance/n2
        }
        else
        {
          v1 = tmp$variances[,1]/tmp$samples[1] # # group 1 variance/n1
          v2 = tmp$variances[,2]/tmp$samples[2] # # group 2 variance/n2
        }
      }
      else
      {
         stop("Incorrect class labels!")
      }
          
      # t statistic
      cat = diff/sqrt( v1 + v2 )

      if (compute.cor) # turn t score into cat score
      {
        if(verbose) cat("Computing the square root of the inverse pooled correlation matrix\n")     
        if (auto.lambda)
          cat = crossprod.powcor.shrink(tmp$centered.data, cat, alpha=-1/2, verbose=FALSE)
        else
          cat = crossprod.powcor.shrink(tmp$centered.data, cat, alpha=-1/2, lambda=lambda, verbose=FALSE)
        if(verbose)
        {
          if(attr(cat, "lambda.estimated") )
            cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
            round(attr(cat, "lambda"), 4), "\n")
          else
            cat("Specified shrinkage intensity lambda (correlation matrix):", 
            round(attr(cat, "lambda"), 4), "\n")
        }

      }

      return( as.double(cat) )
    }
}


# private function

# create data matrix for paired (ca)t score
pvt.pairX = function(X, L)
{
    lev = levels(factor(L)) 
    if(length(lev) != 2) stop("You need to specify labels for two groups!")
    g1 = (L == lev[1] )
    g2 = (L == lev[2] )
    if(sum(g1) != sum(g2)) stop("Sample sizes in the two groups must be equal!")

    X.new = X[g1,, drop=FALSE]-X[g2,, drop=FALSE]

    return(X.new)
}


