\name{PEI}
\alias{PEI}
\docType{data}
\title{
 PEI Study Pancreatic Cancer Dataset
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
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse16515
}
\references{
Pei et al, Cancer Cell 2009
}
\examples{
PEIEset = loadPancreasEsets()$esets$PEI
  experimentData(PEIEset)

}
\keyword{datasets}
