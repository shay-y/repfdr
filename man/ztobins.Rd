\name{ztobins}
\alias{ztobins}
\title{
Binning of z-scores and estimation of the probabilities in each bin for the null and non-null states.
}
\description{
For each study, the function discretizes the z-scores into bins and estimates the probabilities in each bin for the null and non-null states using the function \code{locfdr}.
}
\usage{
ztobins(zmat, n.association.status = 3, n.bins = 120,nulltype = 0, ...)
}
\arguments{
  \item{zmat}{
Matrix of z-scores of the features (in rows) in each study (columns).
}
  \item{n.association.status}{
either 2 for no-association\\association or 3 for no-associtation\\negative-association\\positive-association.
}
  \item{n.bins}{
Number of bins in the discretization of the z-score axis.
}
  \item{nulltype}{
Type of null hypothesis assumed in estimating \code{f0(z)}. \code{0} (default) is the theoretical null N(0,1). This is an argument to pass to \code{locfdr}. For other values, see documentation of \code{locfdr}.
}
  \item{\dots}{
Arguments to pass to \code{locfdr}. See \code{locfdr} package for details.
}
}
\details{
This utility function outputs the first two arguments to be input in the main function \code{\link{repfdr}}.
}
\value{
A list with:
  \item{pdf.binned.z}{
A 3-dimensional array which contains for each study (first dimension), the probabilities of a z-score to fall in the bin (second dimension), under each hypothesis status (third dimension). The third dimension can be of size 2 or 3, depending on the number of association states: if the association can be either null or only in one direction, the dimension is 2; if the association can be either null, or positive, or negative, the dimension is 3.  
}
  \item{binned.z.mat}{
A matrix of the bin numbers for each the z-scores (rows) in each study (columns).
}
}
\seealso{
\code{\link{repfdr}}
}
\examples{
\dontrun{
# download.file(url = "http://cran.r-project.org/src/contrib/Archive/locfdr/locfdr_1.1-7.tar.gz",
#           destfile = "locfdr_1.1-7.tar.gz")
# install.packages("locfdr_1.1-7.tar.gz", type="source", repos=NULL)
# unlink("locfdr_1.1-7.tar.gz")
library(locfdr)

# three association states case (H in {-1,0,1}):
download.file('http://www.math.tau.ac.il/~ruheller/repfdr_RData/zmat.RData',destfile = "zmat.RData")
load(file = "zmat.RData")

input.to.repfdr3 <- ztobins(zmat, 3, plot = TRUE, df = 15)
pbz    <- input.to.repfdr3$pdf.binned.z
bz     <- input.to.repfdr3$binned.z.mat

# two association states case (H in {0,1}):
data(zmat_sim)

input.to.repfdr <- ztobins(zmat_sim, 2, plot = TRUE, df = 15)
pbz_sim    <- input.to.repfdr$pdf.binned.z
bz_sim     <- input.to.repfdr$binned.z.mat
}
}