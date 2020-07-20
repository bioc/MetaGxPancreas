\name{BADEA}
\alias{BADEA}
\docType{data}
\title{
 BADEA Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse15471
}
\references{
Badea et al, Hepatogastroenterology 2008
}
\examples{
BADEAEset = loadPancreasEsets()$esets$BADEA
experimentData(BADEAEset) 
}
\keyword{datasets}
