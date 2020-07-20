\name{YANG}
\alias{YANG}
\docType{data}
\title{
 YANG Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452
}
\references{
Yang et al, 2016, Cancer Research
}
\examples{
YANGEset = loadPancreasEsets()$esets$YANG
experimentData(YANGEset) 
}
\keyword{datasets}
