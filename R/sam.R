### sam.R  (2006-08-30)
###
###    SAM t Statistic
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


# Note: these function2 require the "samr" library

sam.stat = function (X, L)
{
  FUN = sam.fun(L=L)
  score = FUN(X)
  
  return( score )
}

sam.fun <- function(L)
{
   require("samr")

   if (missing(L))
      stop("Please specify to which group (1 or 2) each sample is assigned!")
  
   function(X)
    {
      dd = list(x=t(X),y=L, logged2=TRUE)
      out = samr(dd, resp.type="Two class unpaired", nperms=1)
  
      return(out$tt) # SAM test statistic
    } 
}
