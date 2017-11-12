[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN version](http://www.r-pkg.org/badges/version/ragt2ridges)](http://cran.r-project.org/package=ragt2ridges)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/ragt2ridges)](http://cran.r-project.org/package=ragt2ridges/index.html)
[![Total CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/ragt2ridges)](http://www.r-pkg.org/pkg/ragt2ridges)


**ragt2ridges**
---------------

The R-package ```ragt2ridges``` performs ridge maximum likelihood estimation of vector auto-regressive processes: the VAR(1) model (more to be added). Prior knowledge may be incorporated in the estimation through *a)* specification of the edges believed to be absent in the time series chain graph, and *b)* a shrinkage target towards which the parameter estimate is shrunken for large penalty parameter values. Estimation functionality is accompanied by methodology for penalty parameter selection.

In addition, the package offers supporting functionality for the exploitation of estimated models. Among others, *i)* a procedure to infer the support of the non-sparse ridge estimate (and thereby of the time series chain graph) is implemented, *ii)* a table of node-wise network summary statistics, *iii)* mutual information analysis, and *iv)* impulse response analysis.  

### Installation

The `` ragt2ridges``  package is available via
[CRAN](http://cran.r-project.org/package=ragt2ridges) (Comprehensive R Archive Network) and can be installed from within R through:

```R
install.packages("ragt2ridges")
```

After installation run `news(package="ragt2ridges")` for the latest changes to the ```ragt2ridges``  package.

Previous versions of **ragt2ridges** are available via the [CRAN archive.](http://cran.r-project.org/src/contrib/Archive/ragt2ridges/)


### References

Publications related to **ragt2ridges** include:

 1. Miok, V., Wilting, S.M., & van Wieringen, W.N. (2017)
    *"Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data"*.
     _Biometrical Journal_, 59(1): 172-191
    ([doi:10.1002/bimj.201500269](http://onlinelibrary.wiley.com/doi/10.1002/bimj.201500269/abstract)). 
2. Miok, V., Wilting, S.M., & van Wieringen, W.N. (2017)
    *"Ridge estimation of network models from time-course omics data"*.
     _submitted_.    
 2. van Wieringen, W.N. (2017). 
    *"ragt2ridges: Ridge Estimation of Vector Auto-Regressive (VAR) Processes"*. 
    _R package_, version 0.2.4
 
Please cite the relevant publications if you use **ragt2ridges**.
