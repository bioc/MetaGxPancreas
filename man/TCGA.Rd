\name{TCGA}
\alias{TCGA}
\docType{data}
\title{
 TCGA Study Pancreatic Cancer Dataset
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
https://portal.gdc.cancer.gov/projects/TCGA-PAAD
}
\references{
TCGA Research Network, Cancer Cell 2017
}
\examples{
TCGAEset = loadPancreasEsets()$esets$TCGA
experimentData(TCGAEset)

}
\keyword{datasets}
