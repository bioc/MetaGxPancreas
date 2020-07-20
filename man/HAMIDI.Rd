\name{HAMIDI}
\alias{HAMIDI}
\docType{data}
\title{
 HAMIDI Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77858
}
\references{
NA
}
\examples{
HAMIDIEset = loadPancreasEsets()$esets$HAMIDI
experimentData(HAMIDIEset) 
}
\keyword{datasets}
