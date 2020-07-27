#' Function to load pancreas cancer expression profiles from the Experiment Hub
#'
#' This function returns pancreas cancer patient cohorts in SummarizedExperiment object
#' from the hub and a vector of patients from the datasets that are duplicates
#'
#' @param removeDuplicates remove patients with a Spearman correlation greater
#'        than or equal to 0.98 with other patient expression profiles (default TRUE)
#' @param quantileCutoff A numeric between 0 and 1 specifying to remove genes with
#'        standard deviation below the required quantile (default 0)
#' @param rescale apply centering and scaling to the expression sets (default FALSE)
#' @param minNumberGenes an integer specifying to remove expression sets with less genes
#'        than this number (default 0)
#' @param minSampleSize an integer specifying the minimum number of patients required in an SE (default 0)
#' @param minNumberEvents an integer specifying how man survival events must be in the dataset to keep the dataset
#'   (default 0)
#' @param removeSeqSubset currently only removes the ICGSSEQ dataset as it contains the same patients as
#'        the ICGS microarray dataset (default TRUE, currently just ICGSSEQ)
#' @param keepCommonOnly remove probes not common to all datasets (default FALSE)
#' @param imputeMissing impute missing expression value via knn
#'
#' @return a list with two elements. The First element named SummarizedExperiments
#'         contains the datasets as Bioconductor SummarizedExperiment objects.
#'         The second element named duplicates contains a vector with patient IDs for the
#'         duplicate patients (those with Spearman correlation greater than or equal to
#'         0.98 with other patient expression profiles).
#'
#' @examples
#' sumExptsAndDuplicates <- loadPancreasDatasets()
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays colData rowData assay
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom S4Vectors metadata
#' @importFrom AnnotationHub query
#' @importFrom stats quantile
#' @importFrom impute impute.knn
#' @importFrom stats sd complete.cases
#'
#' @export
loadPancreasDatasets <- function(removeDuplicates = TRUE, quantileCutoff = 0, rescale = FALSE,
                             minNumberGenes = NA, minSampleSize = NA, minNumberEvents = NA,
                             removeSeqSubset = FALSE, keepCommonOnly = FALSE,
                             imputeMissing = FALSE)
{
    ## -----------------------------------------------------------------------------
    ## Helper function
    ## -----------------------------------------------------------------------------
    .filterQuantile <- function(SE, quantile){
        message(metadata(SE)$name)
        exprs <- assay(SE, 'exprs')

        geneSd <- apply(exprs, 1, sd, na.rm=TRUE)
        geneQu <- stats::quantile(geneSd, probs = quantile, na.rm=TRUE)
        cutoff <- sum(geneSd < geneQu, na.rm=TRUE) / length(geneSd)

        if ( abs(quantile - cutoff) <= 0.01 ){
        exprs <- exprs[intersect(geneSd > geneQu, !is.na(geneSd)), ]

            return(SummarizedExperiment(metadata=list(name=metadata(SE)$name),
                    assays=list(exprs=exprs),
                    rowData=rowData(SE)[which(rowData(SE)$gene %in% rownames(exprs)),],
                    colData=colData(SE)))

            message()
        } else {
            return(SE)
        }
    }

    `%notin%` <- Negate(`%in%`)

    ## -----------------------------------------------------------------------------
    ## Load the SEs
    ## -----------------------------------------------------------------------------

    hub = ExperimentHub::ExperimentHub()
    pancreasData = query(hub, c("MetaGxPancreas", "SummarizedExperiment"))
    suppressMessages({ SEs = lapply(pancreasData, function(x) x[[names(x)]]) })
    names(SEs) = pancreasData$title

    duplicates = c("ICGC_0400", "ICGC_0402", "GSM388116", "GSM388118",
        "GSM388120", "GSM388145", "GSM299238", "GSM299239", "GSM299240")

    if (removeDuplicates) {
        SEs = lapply(SEs, function(SE) {
            keepSamples <- rownames(colData(SE)) %notin% duplicates
            SE[, keepSamples]
        })
        message( sprintf("Filtered out duplicated samples: %s",
            paste(duplicates, collapse=", ")) )
    }

    if(quantileCutoff > 0 && quantileCutoff < 1){
        SEs <- lapply(SEs, function(SE) .filterQuantile(SE,
                                                        quantile=quantileCutoff)
                      )
    }

    if (rescale) {
        SEs <- lapply(SEs, function(SE) {
            return(SummarizedExperiment(metadata=list(name=metadata(SE)$name),
                assays=list(exprs=t(scale(t(assays(SE)$exprs)))),
                rowData=rowData(SE),
                colData=colData(SE)))
        })
        message(paste0("Done: Expression profiles rescaled (rescale = ", rescale, ")"))
    }

    if (!is.na(minNumberGenes)) {
        SEs <- lapply(SEs, function(SE) {
            if (length(rownames(SE)) > minNumberGenes) {
                return(SE)
            } else {
                message(paste0("Filtered out: ", metadata(SE)$name, " (minNumberGenes < ",
                               minNumberGenes, ")"))
            }
        })
        SEs <- SEs[!sapply(SEs, is.null)]
    }

    if (!is.na(minSampleSize)) {
        SEs <- lapply(SEs, function(SE) {
            if (length(colnames(SE)) > minSampleSize) {
               return(SE)
            } else {
                message(paste0("Filtered out: ", metadata(SE)$name,
                " (minSampleSize < ", minSampleSize, ")"))
            }
        })
        SEs <- SEs[!sapply(SEs, is.null)]
    }

    if (!is.na(minNumberEvents)) {
        SEs <- lapply(SE, function(SE) {
            if (length(which(!is.na(colData(SE)$vital_status))) > minNumberEvents) {
                return(SE)
            } else {
                message(paste0("Filtered out: ", metadata(SE)$name,
                    " (minNumberEvents < ", minNumberEvents, ")"))
            }
        })
        SEs <- SEs[!sapply(SEs, is.null)]
    }

    if (removeSeqSubset) {
        SEs$ICGCSEQ = NULL
        SEs <- SEs[!sapply(SEs, is.null)]
    }

    if (keepCommonOnly) {
        featuresL <- lapply(SEs, rownames)
        commonFeatures <- Reduce(intersect, featuresL)
        SEs <- lapply(SEs, function(SE){
            return(SummarizedExperiment(metadata=list(name=metadata(SE)$name),
                assays=list(exprs=assays(SE)$exprs[which(rownames(assays(SE)$exprs)
                    %in% commonFeatures), , drop=FALSE]),
                rowData=rowData(SE)[which(rownames(rowData(SE)) %in% commonFeatures),],
                colData=colData(SE)))
        })
    }

    if (imputeMissing) {
        missingData <- which(vapply(SEs, function(SE)
            sum(!complete.cases(assays(SE)$exprs)) > 0, numeric(1)) == 1)
        message( sprintf("Missing data in : %s",
            paste(names(missingData), collapse=", ")) )
        if (length(missingData) > 0) {
            for (idx in missingData) {
              SE <- SEs[[idx]]
              SEs[[idx]] <- SummarizedExperiment(metadata=list(name=metadata(SE)$name),
                                                 assays=list(exprs=impute::impute.knn(assays(SE)$exprs)$data),
                                                 rowData=rowData(SE), colData=colData(SE))
            }
        }
    }

    metaGxPancreas = list(SEs, duplicates)
    names(metaGxPancreas) = c("SummarizedExperiments", "duplicates")

    return(metaGxPancreas)
}

