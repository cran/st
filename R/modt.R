### modt.R  (2006-08-30)
###
###    Moderated t Statistic
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


# Note: these function2 require the "limma" library


modt.stat = function (X, L)
{
  FUN = modt.fun(L=L)
  score = FUN(X)
  
  return( score )
}

modt.fun <- function (L)
{
    require("limma")
    
    if (missing(L))
      stop("Please specify to which group (1 or 2) each sample is assigned!")

    function(X)
    {
      d <- cbind(rep(1, length(L)), L)
      fit <- lmFit(t(X), design=d)
      eb.out <- ebayes(fit)
      modt <- -eb.out$t[,2]
  
      return(modt) 
    } 
}
