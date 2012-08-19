### shrinkcat.R  (2012-08-19)
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


shrinkcat.stat = function (X, L, verbose=TRUE)
{
  FUN = shrinkcat.fun(L=L, verbose=verbose)
  score = FUN(X)
  
  return( score )
}


shrinkcat.fun = function (L, verbose=TRUE)
{
    if (missing(L)) stop("Class labels are missing!")

    function(X)
    {
      cat =  catscore(X, L, diagonal=FALSE, verbose=verbose)[,1]
      return(cat)
    }
}

