\name{HAIDER}
\alias{HAIDER}
\docType{data}
\title{
 HAIDER Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56560
}
\references{
haider et al, Genome medicine, 2014
}
\examples{
HAIDEREset = loadPancreasEsets()$esets$HAIDER
experimentData(HAIDEREset) 
}
\keyword{datasets}
