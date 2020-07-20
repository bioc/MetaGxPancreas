\name{CHEN}
\alias{CHEN}
\docType{data}
\title{
 CHEN Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57495
}
\references{
Chen et al, PLoS One 2015
}
\examples{
CHENEset = loadPancreasEsets()$esets$CHEN
  experimentData(CHENEset)

}
\keyword{datasets}
