\name{NEWS}
\title{News for Package \pkg{ragt2ridges}}
\newcommand{\cpkg}{\href{http://cran.r-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.3.0}{
 \itemize{
    \item More extensive description added to the description field of the description file. 
    \item Correct package name inserted in registerDynamicSymbol.c. 
    \item R code aligned for nicer viewing on Github.	
    \item The \code{hpvP53mir} data object is replaced by the correct version.
    \item Several typo's in manuals corrected.
    \item More consistent wording across manuals introduced.
    \item Citation updated. 
    \item Functions \code{momentS}, \code{ridgePchordal}, \code{optPchordal}, and \code{support4ridgep} are now imported from the \code{rags2ridges}-package.
    \item Three new modules are added: fused ridge estimation of multiple VAR(1) models, and ridge estimation of the VAR(2) and VARX(1) models.		
  }
}

\section{Changes in version 0.2.4}{
 \itemize{
    \item Bug repaired in the \code{graphVAR1}-function. Input combination \code{type="globalPC"} and \code{prune=TRUE} produced a warning. 
    \item Some code from Rcpp-functions replaced by slightly more efficient code.     
    \item Meaningless code from Rcpp-functions removed. 
    \item Citation updated. 
    \item More data from the same experiment have been added. Documentation modified accordingly. Of note: the original data have been renamed to \code{hpvP53rna}.
    \item A misnomer in the output and documentation of the \code{nodeStatsVAR1}-function fixed.
    \item The formula of the VAR(1) model is added to documentation of the \code{ridgeVAR1}-function.  
  }
}

\section{Changes in version 0.2.3}{
 \itemize{
    \item Some re-organization of the src-directory. 
    \item Rcpp-functions made exportable to other packages.     
    \item Input checks of the \code{ridgeVAR1}-function contained inconsistencies: repaired.
    \item Function \code{momentS} added: evaluates the moments of the sample (inverse) covariance matrix.    
  }
}

\section{Changes in version 0.2.2}{
 \itemize{
    \item Change to the description and namespace-file. Importing now from methods, for compatibility with Bioconductor 3.4.
    \item A bug repaired in the \code{ridgeVAR1}-function. The option-combination \code{fitA="ss"} and a non-empty matrix being passed to the argument \code{zerosA} gave an error (an argument was not supplied to an internal C++ function).
    \item A bug repaired in the \code{ridgeVAR1}-function. When specifying the option-combination \code{fitA="ss"} and non-empty matrices being passed to the arguments \code{zerosA} and \code{zerosP} the wrong order of arguments was supplied to an internal C++ function.
  }
}

\section{Changes in version 0.2.1}{
  \itemize{
    \item Extra argument (\code{TburnIn}) for \code{dataVAR1}-function. It specifies number of burn-in time points before data is sampled from the process.
    \item The \code{impulseResponseVAR1}-function and its arguments changed slightly. Previously, the function aggregated all impulse responses between a variate at the current time point and those at the \code{T} future time points. The current implemention evaluates the impulse response in each variate at the current time point \code{t} on those at time point \code{t+T}. Apologies for inconvenience caused. 
    \item The \code{mutualInfoVAR1}-function and its arguments changed slightly. Previously, the function aggregated all mutual informations between a variate at the current time point and those at the \code{T} future time points. The current implemention evaluates the mutual information between each variate at the current time point \code{t} and those at time point \code{t+T}. Apologies for inconvenience caused. 
    \item The \code{nodeStatsVar1}-function changed internally to accomodate the changes in the \code{impulseResponseVAR1} and \code{mutualInfoVAR1}-functions.
    \item Extra argument (\code{zeros}) for the \code{sparsifyVAR1}-function. It specifies which elements in the auto-regression coefficient matrix are known to zero and should not be taken along in the sparsification.  
  }
}

\section{Changes in version 0.2.0}{
  \itemize{
    \item Added a \code{NEWS} file.
    \item Modified the example code of the \code{evaluateVAR1fit}. In this version, the \code{evaluateVAR1fit} function is no longer executed in the example code. This avoids saving .eps files during the package checking. 
  }
}

