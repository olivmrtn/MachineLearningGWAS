# annotateMDT ----

#' @title Annotate MDT Objects
#'
#' @description Appends annotations for each
#' \code{FeatureID} to \code{annotations(mdt)}.
#'
#' @name annotateMDT
#' @rdname annotateMDT
#' @aliases
#' annotateStat
#' annotateLocation
#' annotateBioconductor
#' annotateGOSim
#' annotateRSID
#' annotateCoding
#' annotateMAF
#' merge-MDT-data.frame-method
#'
#' @inheritParams AnnotationDbi::AnnotationDb
#' @inheritParams GOSemSim::godata
#' @inheritParams GOSemSim::mgoSim
#' @inheritParams VariantAnnotation::locateVariants
#' @inheritParams VariantAnnotation::predictCoding
#' @inheritParams data.table::merge
#'
#' @param x \code{\link{MDT}} object.
#'
#' @param y \code{data.frame} with a column named \code{FeatureID}.
#'
#' @param annotations_db \code{\link[AnnotationDbi]{AnnotationDb}} object.
#' Can also be \code{\link[VariantAnnotation]{PolyPhenDb}} or
#' \code{\link[VariantAnnotation]{SIFTDb}}.
#'
#' @param columns Column names to be appended to \code{annotations(mdt)}.
#' For \code{annotateLocation}, can me any column returned by
#' \code{\link[VariantAnnotation]{locateVariants}}.
#' For \code{annotateCoding}, can be any column returned by
#' \code{\link[VariantAnnotation]{predictCoding}}.
#'
#' @param function_list \code{list} of \code{function}s that take two arguments,
#' a dropped column of \code{mtable(mdt)} and \code{response(mdt)}, and returns
#' one \code{numeric} value.
#' \code{list} \code{names} are used as \code{colnames} in
#' \code{annotations(mdt)}.
#'
#' @param keys_colname Name of column to be used as \code{key} argument in
#' \code{\link[AnnotationDbi]{select}}.
#'
#' @param name Name of appended column to \code{annotations(mdt)}.
#'
#' @param seqSource \code{\link[BSgenome]{BSgenome-class}} instance or
#' a FaFile to be used for sequence extraction.
#'
#' @param snp_loc \code{\link[BSgenome]{SNPlocs}} object that contains
#' known SNP locations for a given organism.
#'
#' @param target_go Target GO list.
#'
#' @param txdb \code{\link[GenomicFeatures]{TxDb}} or
#' \code{\link[GenomicRanges]{GRangesList}} object that
#' serves as the annotation.
#' GFF files can be converted to \code{\link[GenomicFeatures]{TxDb}} objects
#' with \code{makeTxDbFromGF} in the \code{GenomicFeatures} package.
#'
#' @param verbose \code{logical} specifying if function should be run in verbose mode.
#'
#' @param fill Value to replace missing values in \code{mtable}.
#'
#' @section merge:
#' Merges a \code{data.frame} to annotations using \code{FeatureID}s.
#'
#' @section annotateStat:
#' Estimates association between \code{mtable(mdt)} feature
#' and \code{response(mdt)}. By default, appends Fisher's exact test p-value and
#' Hardy-Weinberg equilibrium exact test p-values.
#' Results can be visualized using \code{\link{manhattanPlot}}.
#'
#' @section annotateMAF:
#' Estimates minor allele frequency.
#'
#' @section annotateLocation:
#' Annotates genomic positions by locating variants
#' with respect to gene function using
#' \code{\link[VariantAnnotation]{locateVariants}}.
#' Please note that \code{GENEID} is renamed as \code{ENTREZID}.
#'
#' @section annotateBioconductor:
#' Annotates variants with Bioconductor
#' \code{\link[AnnotationDbi]{AnnotationDb}} objects using
#' \code{\link[AnnotationDbi]{select}}.
#'
#' @section annotateGOSim:
#'
#' Estimates GO semantic similarity scores between a gene vector
#' and a target GO terms vector.
#'
#' Uses \code{\link[GOSemSim]{mgoSim}} to estimate GO scores.
#' \code{annotation(mdt)} must contain \code{GO} and \code{ONTOLOGY} columns
#' before running function. Please call \code{annotateBioconductor} to do so.
#'
#' @section annotateRSID:
#' Annotates dbSNP RS identification numbers
#' using \code{\link[GenomicRanges]{findOverlaps}}.
#'
#' @section annotateCoding:
#' Predict coding changes (e.g. synonymous, frameshift) for variants by using
#' \code{\link[VariantAnnotation]{predictCoding}}.
#'
#' @seealso
#' \code{\link[AnnotationDbi]{AnnotationDb}}
#' \code{\link[VariantAnnotation]{locateVariants}}
#' \code{\link[GenomicFeatures]{TxDb}}
#' \code{\link[VariantAnnotation]{predictCoding}}
#' \code{\link[BSgenome]{BSgenome}}
#' \code{\link[VariantAnnotation]{PolyPhenDb}}
#' \code{\link[VariantAnnotation]{SIFTDb}}
#' \code{\link[BSgenome]{SNPlocs}}
#' \code{\link[GOSemSim]{mgoSim}}
#'
NULL

# plotMLGWAS ----

#' @title MLGWAS Plots
#'
#' @description Plots performance, feature importance and stability, or errors of
#' \code{\link{MLGWAS}} objects.
#'
#' @name plotMLGWAS
#' @rdname plotMLGWAS
#' @aliases
#' plotPerf
#' plotFeatimp
#' plotErrors
#' compareFeatimp
#'
#' @inheritParams ggplot2::geom_histogram
#'
#' @param x \code{\link{MLGWAS}} object.
#'
#' @param measure \code{character} specifying performance measure to be plotted.
#' Can be any column name in \code{performance(x)} that is not an identifier.
#'
#' @param group Facultative \code{list} containing \code{ClassifierID}s.
#' \code{ClassifierID}s in the same list entry will be plotted side by side.
#' By default, all \code{ClassifierID} are plotted separately.
#'
#' @param features Subset of \code{FeatureID}s.
#'
#' @param sort_by_classifiers Classifier identifiers
#' \code{classifierNames(x)} to be used to sort feature importance.
#'
#' @param max_features \code{integer} specifying maximum number of features in plot.
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#'
NULL

# summaryMLGWAS ----
#' @title Summary of MLGWAS Performance, Feature Importance and Errors
#'
#' @description Summarizes performance, feature importance or errors of
#' \code{\link{MLGWAS}} objects throughout partitions.
#'
#' @name summaryMLGWAS
#' @rdname summaryMLGWAS
#' @aliases
#' summaryPerf
#' summaryFeatimp
#' summaryErrors
#'
#' @param x \code{\link{MLGWAS}} object.
#'
#' @param summarize \code{list} of \code{function}s to summarize performance.
#' Each function must take a vector of \code{numeric}s
#' and return one \code{numeric} value.
#'
#' @return Summarized \code{\link{data.table}}
#'
NULL
