\name{NEWS}
\title{News for Package \pkg{ragt2ridges}}
\newcommand{\cpkg}{\href{http://cran.r-project.org/package=#1}{\pkg{#1}}}

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
