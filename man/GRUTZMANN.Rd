\name{GRUTZMANN}
\alias{GRUTZMANN}
\docType{data}
\title{
 GRUTZMANN Study Pancreatic Cancer Dataset
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
https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-950/?query=pilarsky&s_page=1&s_pagesize=50
}
\references{
Grutzmann et al, Neoplasia, 2004
}
\examples{
GRUTZMANNEset = loadPancreasEsets()$esets$GRUTZMANN
 experimentData(GRUTZMANNEset)

}
\keyword{datasets}
