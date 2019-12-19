\name{mdDelay}
\alias{mdDelay}

\title{Find optimal delay from a multi-dimensional dataset.}

\description{Estimates time delay for embedding of a 
multi-dimensional dataset.}

\usage{ mdDelay(data, nbins, maxlag, criterion, threshold) }

\arguments{
  \item{data}{The matrix containing all variables}
  \item{nbins}{The number of bins considered to estimate mutual information}
  \item{maxlag}{Number of lags considered}
  \item{criterion}{A string to indicate what delay optizes
  mutual information: 'firstBelow' uses the lowest delay 
  at which the AMI function drops below the value set 
  by the threshold parameter. 'localMin' uses the position 
  of the first local minimum of the AMI function.
  The categorical state on which phi is calculated}
   \item{threshold}{Value to select the delay when AMI drops below it.}
}

\value{

  It returns the recurrence phi-coefficient profile for
  state k for all delays considered

}

\references{Wallot, S., and Moenster, D. (2018). Calculation of 
average mutual information (AMI) and false-nearest neighbors (FNN) 
for the estimation of embedding parameters of multidimensional time
-series in Matlab. Front. Psychol. - Quantitative Psychology and
Measurement}

\author{Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics
     Dan Moenster, Aarhus University, 
     Moreno I. Coco, University of East London}

\seealso{\code{\link{mdFnn}}, \code{\link{optimizeParam}}} 

\examples{

## TODO


}

\keyword{ts}
