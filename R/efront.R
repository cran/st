### efront.R  (2008-10-27)
###
###    Efron t Statistic (2001)
###
### Copyright 2006-2008 Rainer Opgen-Rhein and Korbinian Strimmer
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


efront.stat = function (X, L, verbose=TRUE)
{
  FUN = efront.fun(L=L, verbose=verbose)
  score = FUN(X)
  
  return( score )
}

efront.fun <- function (L, verbose=TRUE)
{
    if (missing(L)) stop("class labels are missing!")
    L = factor(L)
    cl = levels(L)
    if (length(cl) != 2) stop("class labels must be specified for two groups, not more or less!")
    idx1 = (L == cl[1])
    idx2 = (L == cl[2])

    function(X)
    {
      tmp = pvt.group.moments(X, idx1, idx2, variances=TRUE)
      
      # differences between the two groups
      diff = tmp$mu1-tmp$mu2
      
      # standard error of diff
      n1 = sum(idx1); n2 = sum(idx2)      
      sd = sqrt( (1/n1 + 1/n2)*tmp$v.pooled )
      
      # tuning parameter
      a0 <- quantile(sd, probs=c(0.9))
      
      if (verbose) cat("Fudge factor a0 =", a0, "\n")
      
      # t statistic
      t = diff/(sd+a0)
                 
      return(t)
    }
}

