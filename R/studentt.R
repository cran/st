### studentt.R  (2008-10-27)
###
###    Student t statistic and related stuff
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


# difference of means  ("fold change")

diffmean.stat = function (X, L)
{
  FUN = diffmean.fun(L=L)
  score = FUN(X)
  
  return( score )
}

diffmean.fun = function (L)
{
    if (missing(L)) stop("class labels are missing!")
    L = factor(L)
    cl = levels(L)
    if (length(cl) != 2) stop("class labels must be specified for two groups, not more or less!")
    idx1 = (L == cl[1])
    idx2 = (L == cl[2])
    
    function(X)
    { 
      tmp = pvt.group.moments(X, idx1, idx2, variances=FALSE)
      
      # differences between the two groups
      diff = tmp$mu1-tmp$mu2
      
      return(diff)
    }
}



# student t statistic

studentt.stat = function (X, L)
{
  FUN = studentt.fun(L=L)
  score = FUN(X)
  
  return( score )
}

studentt.fun = function (L)
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
      
      # t statistic
      t = diff/sd
      
      return(t)
    }
}


### private utility function

# compute group means and pooled variances
pvt.group.moments = function(X, idx1, idx2, variances=TRUE)
{
     # two groups
     X1 = X[idx1,]; n1 = sum(idx1)
     X2 = X[idx2,]; n2 = sum(idx2)

     if (n1 == 0 || n2 == 0 || n1+n2 < 3)
     {
       stop("increase sample size!")
     }
      
     # means of each group
     mu1 = colMeans(X1)
     mu2 = colMeans(X2)
     
     if (variances)
     {
       # pooled variance
       r1 = colSums(X1^2)-n1*mu1^2
       r2 = colSums(X2^2)-n2*mu2^2
       # r1 = (n1-1)*apply(X1, 2, var)   # equivalent but much slower code
       # r2 = (n2-1)*apply(X2, 2, var)
       
       v.pooled  = (r1+r2)/(n1+n2-2) 
       v1 = r1/(n1-1)
       v2 = r2/(n2-1)
     }
     else
     {
       v.pooled = NULL
       v1 = NULL
       v2 = NULL
     }
          
     return( list(mu1=mu1, mu2=mu2, v.pooled=v.pooled, v1=v1, v2=v2) )
}

