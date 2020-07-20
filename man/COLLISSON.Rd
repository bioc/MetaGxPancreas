\name{COLLISSON}
\alias{COLLISSON}
\docType{data}
\title{
 COLLISON Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17891
}
\references{
Collisson et al.,Nat Med 2011
}
\examples{
COLLISSONEset = loadPancreasEsets()$esets$COLLISSON
experimentData(COLLISSONEset)

}
\keyword{datasets}
