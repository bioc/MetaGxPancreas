\name{LUNARDI}
\alias{LUNARDI}
\docType{data}
\title{
 LUNARDI Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55643
}
\references{
Lunardi S et al, 2014, Oncotarget
}
\examples{
LUNARDIEset = loadPancreasEsets()$esets$LUNARDI
experimentData(LUNARDIEset) 
}
\keyword{datasets}
