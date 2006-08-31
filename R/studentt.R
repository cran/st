### studentt.R  (2006-08-30)
###
###    Student t statistic and related stuff
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


# difference of means  ("fold change")

diffmean.stat = function (X, L)
{
  FUN = diffmean.fun(L=L)
  score = FUN(X)
  
  return( score )
}

diffmean.fun <- function (L)
{
    if (missing(L))
      stop("Please specify to which group (1 or 2) each sample is assigned!")
    
    function(X)
    { 
      tmp = pvt.group.moments(X, L, variances=FALSE)
      
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

studentt.fun <- function (L)
{
    if (missing(L))
      stop("Please specify to which group (1 or 2) each sample is assigned!")

    function(X)
    { 
      tmp = pvt.group.moments(X, L, variances=TRUE)
      
      # differences between the two groups
      diff = tmp$mu1-tmp$mu2
      
      # standard error of diff
      n1 = sum(L==1); n2 = sum(L==2)      
      sd = sqrt( (1/n1 + 1/n2)*tmp$v.pooled )
      
      # t statistic
      t = diff/sd
      
      return(t)
    }
}


### private utility function

# compute group means and pooled variances
pvt.group.moments = function(X, L, variances=TRUE)
{
     # two groups
     X1 = X[L==1,]; n1 = sum(L==1)
     X2 = X[L==2,]; n2 = sum(L==2)

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

