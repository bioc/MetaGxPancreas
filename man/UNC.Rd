\name{UNC}
\alias{UNC}
\docType{data}
\title{
 UNC Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729
}
\references{
Moffitt et al, Nat Genet 2015
}
\examples{
UNCEset = loadPancreasEsets()$esets$UNC
 experimentData(UNCEset)

}
\keyword{datasets}
