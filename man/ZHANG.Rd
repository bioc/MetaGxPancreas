\name{ZHANG}
\alias{ZHANG}
\docType{data}
\title{
 ZHANG Study Pancreatic Cancer Dataset
}
\description{
An expresion set object containing pancreatic cancer data
}
\format{
    SummarizedExperiment object
}
\details{
    more details can be found in the experimentData section of the object
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735
}
\references{
Zhang et al, PLoS One 2012
}
\examples{
ZHANGEset = loadPancreasEsets()$esets$ZHANG
experimentData(ZHANGEset)

}
\keyword{datasets}
