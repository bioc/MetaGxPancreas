\name{OUH}
\alias{OUH}
\docType{data}
\title{
 OUH Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60980
}
\references{
Sandhu et al, Mol Onc, 2015
}
\examples{
OUHEset = loadPancreasEsets()$esets$OUH
experimentData(OUHEset)

}
\keyword{datasets}
