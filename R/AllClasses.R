################################################################################
## All Classes
################################################################################

#################
## MDT ----
#################

checkValidityMDT <- function(object) {

  errors <- base::character()

  ## Retrieve data
  mtable <- MachineLearningGWAS::mtable(object)
  annotations <- MachineLearningGWAS::annotations(object)
  phenotype <- MachineLearningGWAS::phenotype(object)
  info <- MachineLearningGWAS::info(object)

  ## Empty object
  if (base::nrow(mtable) == 0) {
    return("Can not create MDT object with empty mtable.")

    ## Check mtable
  } else if (any(! c("FeatureID", "SampleID", "VALUE") %in% base::colnames(mtable))) {
    errors <- c(errors, "mtable(x) must have columns: 'FeatureID', 'SampleID', 'VALUE'.")
  }

  ## Check annotations
  if (base::nrow(annotations) != 0) {
    if (any(! c("FeatureID") %in% base::colnames(annotations))) {
      errors <- c(errors, "annotations(x) must have columns: 'FeatureID'.")
    }
  }

  ## Check info
  if (base::nrow(info) != 0) {
    if (any(! c("FeatureID", "SampleID") %in% base::colnames(info))) {
      errors <- c(errors, "info(x) must have columns: 'FeatureID', 'SampleID'.")
    }
  }

  ## Check phenotype
  if (base::nrow(phenotype) != 0) {
    if (any(! c("SampleID", "RESPONSE", "SEX", "AGE") %in% base::colnames(phenotype))) {
      errors <- c(errors, "phenotype(x) must have columns: 'SampleID', 'RESPONSE', 'SEX', 'AGE'.")
    }
  }

  if (length(errors) == 0) { TRUE
  } else { errors }
}

#' Mutation Data Table
#'
#' The MDT class can contain information about sequence variation, annotations
#' and phenotype for a collection of samples. It is meant to simplify
#' machine learning in genome-wide association studies.
#'
#' @name MDT-class
#' @rdname MDT-class
#' @exportClass MDT
#' @aliases MDT
#'
#' dim,MDT-method
#' length,MDT-method
#' nrow,MDT-method
#' ncol,MDT-method
#' show,MDT-method
#'
#' merge,MDT-method
#'
#' [,MDT-method
#' [,MDT,character,character,ANY-method
#' [,MDT,character,missing,ANY-method
#' [,MDT,missing,character,ANY-method
#'
#' samples
#' samples,MDT-method
#' samples<-
#' samples<-,MDT-method
#' features
#' features,MDT-method
#' features<-
#' features<-,MDT-method
#' mtable
#' mtable,MDT-method
#' mtable<-
#' mtable<-,MDT-method
#' annotations
#' annotations,MDT-method
#' annotations<-
#' annotations<-,MDT-method
#' phenotype
#' phenotype,MDT-method
#' phenotype<-
#' phenotype<-,MDT-method
#' info
#' info,MDT-method
#' info<-
#' info<-,MDT-method
#' response
#' response,MDT-method
#' response<-
#' response<-,MDT-method
#'
#' @param x \code{MDT} object.
#' @param value \code{character} vector.
#'
#' @slot mtable Variation matrix.
#' Must contain the following columns:
#' \code{FeatureID},
#' \code{SampleID},
#' \code{VALUE}.
#'
#' @slot annotations Feature annotations.
#' Must contain the following columns:
#' \code{FeatureID}.
#'
#' @slot phenotype Sample annotations.
#' Must contain the following columns:
#' \code{SampleID},
#' \code{RESPONSE} (variable to be predicted),
#' \code{SEX},
#' \code{AGE}.
#'
#' @slot info Feature and Sample annotations. Annotations that are dependant on
#' both sample and feature such as quality measures or observed alleles.
#' Must contain the following columns:
#' \code{FeatureID},
#' \code{SampleID}.
#'
#' @slot features \code{character} \code{vector} of unique \code{FeatureID}s.
#' Each combination of chromosome and position makes one unique \code{FeatureID}.
#'
#' @slot samples \code{character} \code{vector} of unique \code{SampleID}s.
#'
#' @section Constructors:
#'
#' \code{MDT} contructs a MDT object. Only the mtable argument is required.
#'
#' \code{\link{vcfsToMDT}} converts a list of VCF objects to MDT without
#' any phenotypic data, that is, with an empty \code{phenotype} slot.
#'
#' \code{\link{importPhenotype}} adds phenotypic data to MDT object.
#'
#' \code{\link{aggregateMDT}} aggregate FeatureIDs Into Higher Level IDs.
#'
#' @section Accessors:
#'
#' \code{mtable(x)}, \code{mtable(x) <- value} gets or sets \code{mtable}.
#'
#' \code{annotations(x)}, \code{annotations(x) <- value} gets or sets
#' \code{annotations}.
#'
#' \code{phenotype(x)}, \code{phenotype(x) <- value} gets or sets
#' \code{phenotype}.
#'
#' \code{info(x)}, \code{info(x) <- value} gets or sets \code{info}.
#'
#' \code{response(x)}, \code{response(x) <- value} gets or sets \code{Response}
#' column in \code{phenotype} slot.
#'
#' \code{features(x)}, \code{features(x) <- value} gets or sets
#' \code{FeatureID}s.
#'
#' \code{samples(x)}, \code{samples(x) <- value} gets or sets \code{SampleID}s.
#'
#' \code{asMatrixMDT(x, miss)} converts \code{mtable} tp a
#' \code{\link[base]{matrix}}
#' with \code{SampleID} as rows and
#' \code{FeatureID}s as columns.
#'
#' @section Subset:
#'
#' \code{x[SampleIDs, FeatureIDs]} subsets \code{MDT} objects using \code{SampleID}s
#' and \code{FeatureID}s. Must be \code{character}.
#'
#' \code{\link{filterMDT}} filter samples and features in \code{MDT}
#' objects using annotations in \code{annotations}, \code{phenotype}
#' and \code{info}.
#'
#' @section Plots:
#'
#' \code{\link{plotMDT}} allows to generate various plots such as PCA, heatmaps,
#' histograms.
#'
#' \code{\link{manhattanPlot}} creates a Manhattan plot.
#'
#' @section Add Annotations:
#'
#' See \code{\link{annotateMDT}} to add feature annotations to annotations(x).
#'
#' See \code{\link{importPhenotype}} to add sample annotations to phenotype(x).
#'
#' Use \code{\link{selectMDT}} to retrieve annotations from annotations(x),
#' phenotyope(x) and info(x).
#'
#' @section Train Classifier:
#'
#' \code{\link{trainClassifier}} trains a classifier, using an
#' \code{MDT} object, that that predicts phenotypic \code{response(x)}
#' using \code{mtable(x)} matrix. Creates a \code{\link{MLGWAS}} object.
#'
#' @seealso
#' \code{\link{MLGWAS}}
#' \code{\link[data.table]{data.table}}
#'
setClass(Class = "MDT",
         representation = representation(mtable = "data.table",
                                         annotations = "data.table",
                                         phenotype = "data.table",
                                         info = "data.table",
                                         features = "character",
                                         samples = "character"),
         validity = checkValidityMDT
)

#################
## MLGWAS ----
#################

checkValidityMLGWAS <- function(object) {

  errors <- character()

  ## Retrieve data
  perf <- MachineLearningGWAS::performances(object)
  featimp <- MachineLearningGWAS::featimp(object)
  err <- MachineLearningGWAS::errors(object)
  clas <- MachineLearningGWAS::classifiers(object)
  pred <- MachineLearningGWAS::predictions(object)
  part <- MachineLearningGWAS::partitions(object)
  call <- object@call

  ## Check if IDs are present
  if (any(! c("ClassifierID", "PartitionID") %in% base::colnames(perf))) {
    errors <- c(errors, "Missing ClassifierID or PartitionID column in performances slot.")
  }

  if (any(! c("ClassifierID", "PartitionID", "FeatureID") %in% base::colnames(featimp))) {
    errors <- c(errors, "Missing ClassifierID, PartitionID or FeatureID column in featimp slot")
  }

  if (any(! c("SampleID", "PartitionID") %in% base::colnames(err))) {
    errors <- c(errors, "Missing SampleID or PartitionID column in errors slot.")
  }

  ## Check if classifiers is a list of list of caret train objects
  if (any(base::sapply(clas, function(x) {
    any(base::sapply(x, function(y) class(y) != "train"))
  }))) {
    errors <- c(errors, "Classifiers must be caret train objects.")
  }

  ## Check if ClassifierIDs match
  perf_id <- base::sort(base::unique(perf$ClassifierID))
  featimp_id <- base::sort(base::unique(featimp$ClassifierID))
  err_id <- base::sort(base::unique(err$ClassifierID))
  clas_id <- base::sort(base::names(clas))
  pred_id <- base::sort(base::names(pred))
  part_id <- base::sort(base::names(part))
  call_id <- base::sort(base::names(call))

  if (any(! base::sapply(list(clas_id, featimp_id, err_id,
                              pred_id, part_id, call_id), identical, perf_id))) {
    errors <- c(errors, "ClassifierIDs do not match.")
  }

  ## Check if PartitionID match
  perf_id <- base::sort(base::unique(perf$PartitionID))
  featimp_id <- base::sort(base::unique(featimp$PartitionID))
  err_id <- base::sort(base::unique(err$PartitionID))
  clas_id <- base::sort(base::unique(base::unlist(base::lapply(clas, names))))
  pred_id <- base::sort(base::unique(base::unlist(base::lapply(pred, names))))
  part_id <- base::sort(base::unique(base::unlist(base::lapply(part, names))))

  if (any(! base::sapply(list(clas_id, featimp_id, pred_id, part_id, err_id), identical, perf_id))) {
    errors <- c(errors, "PartitionIDs do not match.")
  }

  if (length(errors) == 0) TRUE
  else errors
}

#' Machine Learning for Genome Wide Association Studies (MLGWAS)
#'
#' The MLGWAS contains one or more classifiers with information about their
#' global performance, local errors and feature importance.
#' It is meant to simplify machine learning in genome-wide association studies.
#'
#' @name MLGWAS-class
#' @rdname MLGWAS-class
#' @exportClass MLGWAS
#' @aliases MLGWAS
#'
#' dim,MLGWAS-method
#' length,MLGWAS-method
#' nrow,MLGWAS-method
#' ncol,MLGWAS-method
#' show,MLGWAS-method
#'
#' bindML
#' bindML,MLGWAS-method
#' [,MLGWAS-method
#' [,MLGWAS,character,ANY,ANY-method
#' [,MLGWAS,logical,ANY,ANY-method
#' [,MLGWAS,numeric,ANY,ANY-method
#'
#' classifiers
#' classifiers,MLGWAS-method
#' classifiers<-
#' classifiers<-,MLGWAS-method
#' featimp
#' featimp,MLGWAS-method
#' featimp<-
#' featimp<-,MLGWAS-method
#' performances
#' performances,MLGWAS-method
#' performances<-
#' performances<-,MLGWAS-method
#' errors
#' errors,MLGWAS-method
#' errors<-
#' errors<-,MLGWAS-method
#' partitions
#' partitions,MLGWAS-method
#' partitions<-
#' partitions<-,MLGWAS-method
#' predictions
#' predictions,MLGWAS-method
#' predictions<-
#' predictions<-,MLGWAS-method
#'
#' classifierNames
#' classifierNames,MLGWAS-method
#' classifierNames<-
#' classifierNames<-,MLGWAS-method
#' partitionNames
#' partitionNames,MLGWAS-method
#' partitionNames<-
#' partitionNames<-,MLGWAS-method
#' features,MLGWAS-method
#' features<-,MLGWAS-method
#' samples,MLGWAS-method
#' samples<-,MLGWAS-method
#'
#' @param x \code{MLGWAS} object.
#'
#' @slot performances Global performance measures of classifiers.
#' First two columns are \code{ClassifierID}s and \code{PartitionID}s.
#' Other columns are performance measures.
#' For categorical responses, these are:
#' sensitivity, specificity,
#' positive predictive value (PPV), negative predictive value (NPV),
#' ROC plot AUC and F1 score.
#' For continuous respones, these are:
#' root-mean-square deviation (RMSD) and coefficient of determination (r2).
#'
#' @slot featimp Variable importance of \code{FeatureID}
#' as estimated by \code{\link[caret]{varImp}}.
#' Columns are \code{ClassifierID}, \code{PartitionID},
#' \code{FeatureID} and \code{Importance}.
#'
#' @slot errors Errors per partition and per sample.
#' Columns are \code{ClassifierID}, \code{PartitionID},
#' \code{SampleID} and \code{Error}.
#' For categorical variables, values are integers: 1 correspond to a mistake,
#' 0 to a success.
#' For continuous, values the square of the difference between the prediction
#' and the actual response.
#' In both cases, NA means that this sample was not part of the test set.
#'
#' @slot classifiers Classifiers.
#' Coded as a \code{list} of \code{list}s of \code{\link[caret]{train}}
#' objects.
#' First level corresponds to \code{ClassifierID}
#' and second level to \code{PartitionID}.
#'
#' @slot partitions Train and test set partitions
#' Coded as a \code{list} of \code{list}s.
#' First level corresponds to \code{ClassifierID}s and second level to
#' \code{PartitionID}s.
#'
#' @slot predictions Predictions of \code{response(mdt)} for test set.
#' Coded as a \code{list} of \code{list}s.
#' First level corresponds to \code{ClassifierID}s and second level to
#' \code{PartitionID}s.
#'
#' @slot call List of calls. Names of list are \code{ClassifierIDs}
#'
#' @section Constructors:
#'
#' \code{\link{trainClassifier}} trains a classifier,
#' using an \code{\link{MDT}} object,
#' that that predicts \code{response(mdt)} affection using \code{mtable(mdt)} or
#' \code{phenotype(mdt)} or both.
#'
#' \code{\link{trainClassifiers}} allows to call \code{\link{trainClassifier}}
#' using different combination of parameters specified as a \code{data.frame}.
#'
#' @section Identifiers:
#'
#' Four identifiers are present in this class:
#' \itemize{
#' \item \code{ClassifierID} as specified by the \code{id} argument
#' in \code{\link{trainClassifier}}.
#' \item \code{PartitionID} as specified by the \code{names} of the
#' output of the \code{partition} argument in \code{\link{trainClassifier}}.
#' \item \code{FeatureID} as specified by the \code{FeatureID}s
#' of the \code{MDT} object.
#' \item \code{SampleID} as specified by the \code{Sample}s
#' of the \code{MDT} object.
#' }
#'
#' @section Subset and Combine:
#'
#' \code{x[ClassifierIDs]} subsets \code{MLGWAS} objects. Must be a \code{character}.
#'
#' \code{bindML(...)} combines \code{MLGWAS} objects.
#'
#' @section Summaries:
#'
#' The following functions summarize results throughout partitions:
#'
#' \code{\link{summaryPerf}} summarizes performance per classifier.
#'
#' \code{\link{summaryFeatimp}} summarizes feature importance per classifier.
#'
#' \code{\link{summaryErrors}} summarizes errors per classifier and per sample.
#'
#' @section Plots:
#'
#' \code{\link{plotPerf}} plots performance as boxplots.
#'
#' \code{\link{plotFeatimp}} plots feature importance as boxplots.
#'
#' \code{\link{plotErrors}} plots sample error as boxplots.
#'
#' \code{\link{manhattanPlot}} creates a Manhattan plot using the variable
#' importance of features as estimated by \code{\link[caret]{varImp}} instead
#' of -log10 pvalues.
#'
#' @section Accessors:
#'
#' \code{performances(x)} and \code{performances(x) <- value} gets or sets
#' \code{performances}.
#'
#' \code{featimp(x)} and \code{featimp(x) <- value} gets or sets
#' \code{featimp}.
#'
#' \code{errors(x)} and \code{errors(x) <- value} gets or sets
#' \code{errors}.
#'
#' \code{classifiers(x)} and \code{classifiers(x) <- value} gets or sets
#' \code{classifiers}.
#'
#' \code{partitions(x)} and \code{partitions(x) <- value} gets or sets
#' \code{partitions}.
#'
#' \code{predictions(x)} and \code{predictions(x) <- value} gets or sets
#' \code{predictions}.
#'
#' \code{classifierNames(x)} and \code{classifierNames(x) <- value} gets or sets
#' \code{ClassifiersID}s.
#'
#' \code{partitionNames(x)} and \code{partitionNames(x) <- value}
#' get or sets \code{PartitionID}s.
#'
#' \code{features(x)} and \code{features(x) <- value} get or sets
#' \code{FeatureID}.
#'
#' \code{samples(x)} and \code{samples(x) <- value} get or sets
#' \code{SampleID}s
#'
#' @seealso \code{\link{trainClassifier}}
#' \code{\link{trainClassifiers}}
#' \code{\link{MDT}}
#' \code{\link[caret]{train}}
#' \code{\link[data.table]{data.table}}
#'
setClass(Class = "MLGWAS",
         representation = representation(performances = "data.table",
                                         featimp = "data.table",
                                         errors = "data.table",
                                         classifiers = "list",
                                         partitions = "list",
                                         predictions = "list",
                                         call = "list"),
         validity = checkValidityMLGWAS
)
