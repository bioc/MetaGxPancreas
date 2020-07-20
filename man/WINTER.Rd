\name{WINTER}
\alias{WINTER}
\docType{data}
\title{
 WINTER Study Pancreatic Cancer Dataset
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
http://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-2780/
}
\references{
Winter et al, PLoS Comput Biol, 2012
}
\examples{
WINTEREset = loadPancreasEsets()$esets$WINTER
experimentData(WINTEREset)

}
\keyword{datasets}
