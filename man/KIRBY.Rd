\name{KIRBY}
\alias{KIRBY}
\docType{data}
\title{
 KIRBY Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79670
}
\references{
Kirby et al., Mol Oncol 2016
}
\examples{
KIRBYEset = loadPancreasEsets()$esets$KIRBY
experimentData(KIRBYEset)

}
\keyword{datasets}
