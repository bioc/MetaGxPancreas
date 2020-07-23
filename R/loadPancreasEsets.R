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
#' @param minSampleSize an integer specifying the minimum number of patients required in an eset (default 0)
#' @param minNumberEvents an integer specifying how man survival events must be in the dataset to keep the             dataset (default 0)
#' @param removeSeqSubset currently only removes the ICGSSEQ dataset as it contains the same patients as
#'        the ICGS microarray dataset (defeault TRUE, currently just ICGSSEQ)
#' @param keepCommonOnly remove probes not common to all datasets (default FALSE)
#' @param imputeMissing impute missing expression value via knn
#'
#' @return a list with two elements. The First element named esets contains the datasets.
#'         The second element named duplicates contains a vector with patient IDs for the
#'         duplicate patients (those with Spearman correlation greater than or equal to
#'         0.98 with other patient expression profiles).
#'
#' @examples
#' esetsAndDups = loadPancreasEsets()
#'
#' @export
#' @import SummarizedExperiment
#' @import ExperimentHub
#' @importFrom stats quantile
#' @importFrom impute impute.knn
#'
loadPancreasEsets = function(removeDuplicates = TRUE, quantileCutoff = 0, rescale = FALSE,
                             minNumberGenes = NA, minSampleSize = NA, minNumberEvents = NA,
                             removeSeqSubset = FALSE, keepCommonOnly = FALSE,
                             imputeMissing = FALSE)
{
    ## -----------------------------------------------------------------------------
    ## Helper function
    ## -----------------------------------------------------------------------------
    .filterQuantile <- function(eset, quantile){
        message(metadata(eset)$name)
        exprs <- assays(eset)$exprs

        geneSd <- apply(exprs, 1, sd, na.rm=TRUE)
        geneQu <- stats::quantile(geneSd, probs = quantile, na.rm=TRUE)
        cutoff <- sum(geneSd < geneQu, na.rm=TRUE) / length(geneSd)

        if ( abs(quantile - cutoff) <= 0.01 ){
        exprs <- exprs[intersect(geneSd > geneQu, !is.na(geneSd)), ]

            return(SummarizedExperiment(metadata=list(name=metadata(eset)$name),
                    assays=list(exprs=exprs),
                    rowData=rowData(eset)[which(rowData(eset)$gene %in% rownames(exprs)),],
                    colData=colData(eset)))

            message()
        } else {
            return(eset)
        }
    }

    `%notin%` <- Negate(`%in%`)
    
    ## -----------------------------------------------------------------------------
    ## Load the esets
    ## -----------------------------------------------------------------------------

    hub = ExperimentHub::ExperimentHub()
    pancreasData = query(hub, c("MetaGxPancreas", "SummarizedExperiment"))
    esets = lapply(pancreasData, function(x) x[[names(x)]])
    names(esets) = pancreasData$title

    duplicates = c("ICGC_0400", "ICGC_0402", "GSM388116", "GSM388118",
        "GSM388120", "GSM388145", "GSM299238", "GSM299239", "GSM299240")
        
    if (removeDuplicates) {
        esets = lapply(esets, function(eset) {
            return(SummarizedExperiment(metadata=list(name=metadata(eset)$name),
                assays=list(exprs=assays(eset)$exprs[,which(colnames(assays(eset)$exprs)
                    %notin% duplicates),drop=FALSE]),
                rowData=rowData(eset),
                colData=colData(eset)[which(rownames(colData(eset)) %notin% duplicates),]))
        })
        message( sprintf("Filtered out duplicated samples: %s",
            paste(duplicates, collapse=", ")) )
    }
    
    if(quantileCutoff > 0 && quantileCutoff < 1){
        esets <- lapply(esets, function(eset) .filterQuantile(eset, q=quantileCutoff))
    }
    
    if (rescale) {
        esets <- lapply(esets, function(eset) {
            return(SummarizedExperiment(metadata=list(name=metadata(eset)$name),
                assays=list(exprs=t(scale(t(assays(eset)$exprs)))),
                rowData=rowData(eset),
                colData=colData(eset)))
        })
        message(paste0("Done: Expression profiles rescaled (rescale = ", rescale, ")"))
    }

    if (!is.na(minNumberGenes)) {
        esets <- lapply(esets, function(eset) {
            if (length(rownames(eset)) > minNumberGenes) {
                return(eset)
            } else {
                message(paste0("Filtered out: ", metadata(eset)$name, " (minNumberGenes < ",
                    minNumberGenes, ")"))
            }
        })
        esets <- esets[!sapply(esets, is.null)]
    }
    
    if (!is.na(minSampleSize)) {
        esets <- lapply(esets, function(eset) {
            if (length(colnames(eset)) > minSampleSize) {
               return(eset)
            } else {
                message(paste0("Filtered out: ", metadata(eset)$name,
                " (minSampleSize < ", minSampleSize, ")"))
            }
        })
        esets <- esets[!sapply(esets, is.null)]
    }

    if (!is.na(minNumberEvents)) {
        eset <- lapply(eset, function(eset) {
            if (length(which(!is.na(colData(eset)$vital_status))) > minNumberEvents) {
                return(eset)
            } else {
                message(paste0("Filtered out: ", metadata(eset)$name,
                    " (minNumberEvents < ", minNumberEvents, ")"))
            }
        })
        esets <- esets[!sapply(esets, is.null)]
    }
                   
    if (removeSeqSubset) {
        esets$ICGCSEQ = NULL
        esets <- esets[!sapply(esets, is.null)]
    }
        
    if (keepCommonOnly) {
        featuresL <- lapply(esets, rownames)
        commonFeatures <- Reduce(intersect, featuresL)
        esets <- lapply(esets, function(eset){
            return(SummarizedExperiment(metadata=list(name=metadata(eset)$name),
                assays=list(exprs=assays(eset)$exprs[which(rownames(assays(eset)$exprs)
                    %in% commonFeatures), , drop=FALSE]),
                rowData=rowData(eset)[which(rownames(rowData(eset)) %in% commonFeatures),],
                colData=colData(eset)))
        })
    }

    if (imputeMissing) {
        missingData <- which(vapply(esets, function(eset)
            sum(!complete.cases(assays(eset)$exprs)) > 0, numeric(1)) == 1)
        message( sprintf("Missing data in : %s",
            paste(names(missingData), collapse=", ")) )
        if (length(missingData) > 0) {
            for (idx in missingData) {
              eset <- esets[[idx]]
              esets[[idx]] <- SummarizedExperiment(metadata=list(name=metadata(eset)$name),
                                  assays=list(exprs=impute::impute.knn(assays(eset)$exprs)$data),
                                  rowData=rowData(eset), colData=colData(eset))
            }
        }
    }
    
    metaGxPancreas = list(esets, duplicates)
    names(metaGxPancreas) = c("esets", "duplicates")
  
    return(metaGxPancreas)
}

