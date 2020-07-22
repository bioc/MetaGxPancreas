#' Function to load pancreas cancer expression profiles from the Experiment Hub
#'
#' This function returns pancreas cancer patient cohorts in SummarizedExperiment object
#' from the hub and a vector of patients from the datasets that are duplicates
#'
#' @param removeDuplicates remove patients with a Spearman correlation greater
#'        than or equal to 0.98 with other patient expression profiles (default TRUE)
#' @param rescale apply centering and scaling to the expression sets (default FALSE)
#' @param minNumberGenes an integer specifying to remove expression sets with less genes
#'        than this number (default 0)
#' @param minSampleSize an integer specifying the minimum number of patients required in an eset (default 0)
#' @param minNumberEvents an integer specifying how man survival events must be in the dataset to keep the             dataset (default 0)
#' @param removeSeqSubset currently only removes the ICGSSEQ dataset as it contains the same patients as
#'        the ICGS microarray dataset (defeault TRUE, currently just ICGSSEQ)
#' @param keepCommonOnly remove probes not common to all datasets (default FALSE)
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
loadPancreasEsets = function(removeDuplicates = TRUE, rescale = FALSE,
                             minNumberGenes = NA, minSampleSize = NA, minNumberEvents = NA,
                             removeSeqSubset = FALSE, keepCommonOnly = FALSE)
{
    ## -----------------------------------------------------------------------------
    ## Load the esets
    ## -----------------------------------------------------------------------------

    hub = ExperimentHub::ExperimentHub()
    pancreasData = query(hub, c("MetaGxPancreas", "SummarizedExperiment"))
    esets = lapply(pancreasData, function(x) x[[names(x)]])
    names(esets) = pancreasData$title
  
    `%notin%` <- Negate(`%in%`)

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
                if (!is.na(minNumberEvents)) {
                    if (length(which(!is.na(colData(eset)$vital_status))) > minNumberEvents) {
                        return(eset)
                    } else {
                        message(paste0("Filtered out: ", metadata(eset)$name,
                            " (minNumberEvents < ", minNumberEvents, ")"))
                    }
                } else {
                    return(eset)
                }
            } else {
                message(paste0("Filtered out: ", metadata(eset)$name,
                " (minSampleSize < ", minSampleSize, ")"))
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

    metaGxPancreas = list(esets, duplicates)
    names(metaGxPancreas) = c("esets", "duplicates")
  
    return(metaGxPancreas)
}

