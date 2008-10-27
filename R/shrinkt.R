### shrinkt.R  (2008-10-27)
###
###    Shrinkage t Statistic
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


shrinkt.stat = function (X, L, var.equal=TRUE, verbose=TRUE)
{
  FUN = shrinkt.fun(L=L, var.equal=var.equal, verbose=verbose)
  score = FUN(X)
  
  return( score )
}


shrinkt.fun = function (L, var.equal=TRUE, verbose=TRUE)
{
    if (missing(L)) stop("class labels are missing!")
    L = factor(L)
    cl = levels(L)
    if (length(cl) != 2) stop("class labels must be specified for two groups, not more or less!")
    idx1 = (L == cl[1])
    idx2 = (L == cl[2])
 
    function(X)
    {
      p = ncol(X)
      n = nrow(X)   
      n1 = sum(idx1)
      n2 = sum(idx2)  
      
      tmp = pvt.group.moments(X, idx1, idx2, variances=FALSE)
      
      # differences between the two groups
      diff = tmp$mu1-tmp$mu2
      
      #adiff = abs(diff)
      #cutoff = quantile(adiff, probs=c(0.5))
      #diff[ (adiff < cutoff) ] = 0   # hard thresholding
      
      if (var.equal) # compute pooled variance
      {	
	# center data
        xc1 = sweep(X[idx1,], 2, tmp$mu1)
        xc2 = sweep(X[idx2,], 2, tmp$mu2)
	
        v = as.vector( var.shrink(rbind(xc1, xc2), verbose=verbose)*(n-1)/(n-2) )
      
        # standard error of diff    
        sd = sqrt( (1/n1 + 1/n2)*v )
      }
      else # allow different variances in each class
      {
        X1 = X[idx1,]
        X2 = X[idx2,]
        v1 = as.vector(var.shrink(X1, verbose=verbose))
        v2 = as.vector(var.shrink(X2, verbose=verbose))
   
        # standard error of diff 
        sd = sqrt( v1/n1 + v2/n2 )
      }
          
      # t statistic
      t = diff/sd

      return(t)
    }
}
