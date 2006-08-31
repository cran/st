### shrinkt.R  (2006-08-30)
###
###    Shrinkage t Statistic
###
### Copyright 2006 Rainer Opgen-Rhein and Korbinian Strimmer
###
###
### This file is part of the `st' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
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


shrinkt.stat = function (X, L, var.equal=TRUE, verbose=TRUE)
{
  FUN = shrinkt.fun(L=L, var.equal=var.equal, verbose=verbose)
  score = FUN(X)
  
  return( score )
}


shrinkt.fun <- function (L, var.equal=TRUE, verbose=TRUE)
{
    if (missing(L))
      stop("Please specify to which group (1 or 2) each sample is assigned!")
    
    function(X)
    {
      p <- ncol(X)
      n <- nrow(X)   
      n1 = sum(L==1)
      n2 = sum(L==2)  
      
      tmp = pvt.group.moments(X, L, variances=FALSE)
      
      # differences between the two groups
      diff = tmp$mu1-tmp$mu2
      
      #adiff = abs(diff)
      #cutoff = quantile(adiff, probs=c(0.5))
      #diff[ (adiff < cutoff) ] = 0   # hard thresholding
      
      if (var.equal) # compute pooled variance
      {	
	# center data
        xc1 = sweep(X[L==1,], 2, tmp$mu1)
        xc2 = sweep(X[L==2,], 2, tmp$mu2)
	
        v <- var.shrink(rbind(xc1, xc2), verbose=verbose)*(n-1)/(n-2)
      
        # standard error of diff    
        sd = sqrt( (1/n1 + 1/n2)*v )
      }
      else # allow different variances in each class
      {
        X1 = X[L==1,]
        X2 = X[L==2,]
        v1 = var.shrink(X1, verbose=verbose)
        v2 = var.shrink(X2, verbose=verbose)
   
        # standard error of diff 
        sd = sqrt( v1/n1 + v2/n2 )
      }
          
      # t statistic
      t <- diff/sd

      return(t)
    }
}
