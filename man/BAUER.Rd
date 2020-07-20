\name{BAUER}
\alias{BAUER}
\docType{data}
\title{
 BAUER Study Pancreatic Cancer Dataset
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
https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1791/
}
\references{
Bauer et al, 2016, Gastroenterology
}
\examples{
BAUEREset = loadPancreasEsets()$esets$BAUER
experimentData(BAUEREset) 
}
\keyword{datasets}
