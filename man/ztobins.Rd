\name{ztobins}
\alias{ztobins}
\title{
Binning of z-scores and estimation of the probabilities in each bin for the null and non-null states.
}
\description{
For each study, the function discretizes the z-scores into bins and estimates the probabilities in each bin for the null and non-null states.
}
\usage{
ztobins(zmat, n.association.status = 3, n.bins = 120, type = 0,
        df = 7, central.prop = 0.5, pi0=NULL, plot.diagnostics = F,
        nr.plots.in.window=1, trim.z=F,trim.z.upper = 8,trim.z.lower = -8,
        force.bin.number = F)
}
\arguments{
  \item{zmat}{
Matrix of z-scores of the features (in rows) in each study (columns).
}
  \item{n.association.status}{
either 2 for no-association\\association or 3 for no-associtation\\negative-association\\positive-association.
}
  \item{n.bins}{
Number of bins in the discretization of the z-score axis (the number of bins is \code{n.bins - 1}). If the number of z-scores per study is small, we set \code{n.bins} to a number lower than the default of 120 (about equals to the square root of the number of z-scores). To override the bin number cap (and create a discretization of the data that is sparse), use the \code{force.bin.number = TRUE} argument.
}
  \item{type}{
Type of fitting used for f; 0 is a natural spline, 1 is a polynomial, in either case with degrees of freedom \code{df} (so total degrees of freedom including the intercept is \code{df+1}).
}
  \item{df}{
Degrees of freedom for fitting the estimated density f(z).
}
  \item{central.prop}{
Central proportion of the z-scores used like the area of zero-assumption to estimate pi0.
}

\item{pi0}{
Sets argument for estimation of proportion of null hypotheses. Default value is NULL (automatic estimation of pi0) for every study. Second option is to supply vector of values between 0 and 1 (with length of the number of studies/ columns of \code{zmat}. These values will be used for pi0.
}

\item{plot.diagnostics}{
  If set to \code{TRUE}, will show disgnostics plot for density estimation for each study. Plot is histogram of counts for each bin (Displayed as white bars), along with fitted density in green. Pink bars represent the observed number of counts in each bins, minus the expected number of null hypotheses by the model (truncated at zero). Red and Orange dashed lines represent the estimated densities for non null distributions fitted by the spline. A blue dashed line represents the density component of Z scores for null SNPS, N(0,1). Default value is \code{False}
}

\item{nr.plots.in.window}{Set the number of plots to be displayed in a single window, if plot.diagnostics is set to \code{TRUE} (usefull for visualizing many studies at once).
}

\item{trim.z}{ If set to \code{TRUE}, Z scores above \code{trim.z.upper} or below \code{trim.z.lower} will be trimmed at their respective limits. Default value if \code{FALSE}
}

\item{trim.z.upper}{
Upper bound for trimming Z scores. Default value is 8
}

\item{trim.z.lower}{
Lower bound for trimming Z scores. Default value is -8
}

\item{force.bin.number}{Set to \code{T} to be able to create a discretization with \code{n.bins>sqrt(nrow(zmat))}.
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

\item{breaks.matrix}{
A matrix with \code{n.bins + 1} rows and {ncol(zmat)} columns, representing for each study the discretization chosed. Values are the between bin breaks. First and last values are the edges of the outmost bins.
}

\item{df}{
Number of degrees of freedom, used for spline fitting of density.
}

\item{proportions}{
Matrix with \code{n.association.status} rows, and \code{ncol(zmat)} columns, giving the estimated proportion of each component, for each study.
}

\item{PlotWarnings}{
Vector of size \code{ncol{zmat}}, keeping the warnings given for each study (available here, in the plots for each study and printed to console). With no warnings given for study, value is \code{NA}
}

}
\seealso{
\code{\link{repfdr}}
}
\examples{
\dontrun{

# three association states case (H in {-1,0,1}):
download.file('http://www.math.tau.ac.il/~ruheller/repfdr_RData/zmat.RData',destfile = "zmat.RData")
load(file = "zmat.RData")

input.to.repfdr3 <- ztobins(zmat, 3, df = 15)
pbz    <- input.to.repfdr3$pdf.binned.z
bz     <- input.to.repfdr3$binned.z.mat

# two association states case (H in {0,1}):
data(zmat_sim)

input.to.repfdr <- ztobins(zmat_sim, 2, n.bins = 100 ,plot.diagnostics = T)
pbz_sim    <- input.to.repfdr$pdf.binned.z
bz_sim     <- input.to.repfdr$binned.z.mat
}
}