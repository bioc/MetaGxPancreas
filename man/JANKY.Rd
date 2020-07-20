\name{JANKY}
\alias{JANKY}
\docType{data}
\title{
 JANKY Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62165
}
\references{
Janky et al, BMC Cancer 2016
}
\examples{
JANKYEset = loadPancreasEsets()$esets$JANKY
experimentData(JANKYEset) 
}
\keyword{datasets}
