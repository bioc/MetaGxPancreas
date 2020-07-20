\name{ICGCSEQ}
\alias{ICGCSEQ}
\docType{data}
\title{
 ICGCSEQ Study Pancreatic Cancer Dataset
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
http://icgc.org/icgc/cgp/68/304/798
}
\references{
Bailey et al, Nature, 2016
}
\examples{
ICGCSEQEset = loadPancreasEsets(removeSeqSubset = FALSE)$esets$ICGCSEQ
experimentData(ICGCSEQEset)

}
\keyword{datasets}
