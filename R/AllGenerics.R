################################################################################
## All Generics
################################################################################

##
## 1. MDT Main Methods ---------------------------------------------------------
##

## vcfsToMDT ----
#' Convert a list of VCF objects into an MDT object
#'
#' Convert a \code{list} of
#' \code{\link[VariantAnnotation]{VCF}} objects, each representing an individual
#' patient, into a \code{\link{MDT}} object.
#'
#' Please note that chromosome X, Y, XY and M/MT are replaced by integers
#' 23, 24, 25 and 26. All other chromosome names are replaced by 99.
#' Please use \code{\link{importPhenotype}} to fill in \code{phenotype} slot.
#'
#' @param vcfs List of \code{\link[VariantAnnotation]{VCF}} objects.
#'
#' @param verbose \code{logical} specifying if function should be run in
#'  verbose mode.
#'
#' @param fill Value for entries absent in VCF file. By default, sets to "0/0"
#' corresponding to homozygote for REF allele.
#'
#' @return \code{\link{MDT}} object with \code{phenotype} slot.
#'
#' @seealso \code{\link[VariantAnnotation]{VCF}}
#' \code{\link[VariantAnnotation]{readVcf}}
#' \code{\link{readVcfDir}}
#' \code{\link{importPhenotype}}
#'
#' @import data.table
#'
#' @exportMethod vcfsToMDT
#'
setGeneric("vcfsToMDT", function(vcfs,
                                 fill = 0L,
                                 verbose = TRUE) {
  standardGeneric("vcfsToMDT")
})

## importPhenotype ----
#' Adds phenotype data to MDT object
#'
#' Adds phenotype data to \code{\link{MDT}} object.
#'
#' Uses \code{\link[data.table]{fread}} to read file.
#'
#' @param x \code{\link{MDT}} object.
#'
#' @param file Filename containing phenotypic information.
#' Must be tab separated, and contain a column named SampleID corresponding
#' to Sample names, one column named RESPONSE corresponding to the
#' response variable, and SEX and AGE columns.
#' Warning: the function is CAPS sensitive.
#'
#' @param response_type Class of RESPONSE. Choose between: \code{'numeric'},
#' \code{'integer'}, \code{'character'} or \code{'factor'}.
#'
#' @param ... Further arguments to pass on to \code{\link[data.table]{fread}}.
#'
#' @return Updates \code{\link{MDT}} object with phenotype data.
#'
#' @seealso \code{\link[data.table]{fread}}
#'
#' @exportMethod importPhenotype
#'
setGeneric("importPhenotype", function(x,
                                       file,
                                       response_type = "factor", ...) {
  standardGeneric("importPhenotype")
})

## aggregateMDT ----
#' Aggregate FeatureIDs Into Higher Level IDs
#'
#' Aggregate \code{FeatureID}s into Higher Level IDs.
#'
#' Melts \code{mtable} using \code{\link[data.table]{melt.data.table}} and
#' then recasts new \code{mtable} using \code{\link[data.table]{dcast}}.
#'
#' @inheritParams data.table::dcast
#'
#' @param x \code{\link{MDT}} object.
#'
#' @param group \code{data.frame} specifying how to aggregate \code{FeatureID}s.
#' The first column is the old \code{FeatureID},
#' the second is the new \code{FeatureID}.
#'
#' @param ... Further arguments passed on to \code{\link[data.table]{dcast}}
#'
#' @return \code{\link{MDT}} object with new \code{FeatureID}s.
#' All annotations are lost.
#'
#' @exportMethod aggregateMDT
#'
setGeneric("aggregateMDT", function(x,
                                    group,
                                    fun.aggregate = function(value) mean(value, na.rm = TRUE), ...) {
  standardGeneric("aggregateMDT")
})

# mdtToGRanges ----
#' Convert MDT to Genomic Ranges
#'
#' Converts \code{\link{MDT}} to \code{\link[GenomicRanges]{GRanges}}
#'
#' @param x \code{\link{MDT}} object.
#'
#' @return \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @exportMethod mdtToGRanges
#'
setGeneric("mdtToGRanges", function(x) {
  standardGeneric("mdtToGRanges")
})

# filterMDT ----

#' Filter MDT Objects
#'
#' Filter out samples or features on the basis of the entries in annotations,
#' phenotype or info.
#'
#' @param x \code{\link{MDT}} object.
#' @param condition \code{logical} expression indicating elements.
#' Missing values are taken as \code{FALSE}.
#' @param with \code{character(1)} specifying which annotation
#' \code{\link[data.table]{data.table}} should be used for filtering.
#'
#' @return Filtered \code{\link{MDT}} object.
#'
#' @exportMethod filterMDT
#'
setGeneric("filterMDT", function(x,
                                 condition,
                                 with = c("annotations", "phenotype", "info")) {
  standardGeneric("filterMDT")
})

#' Retrieve annotations from MDT objects.
#'
#' Retrieve annotations from annotations(x), phenotyope(x) and info(x) in MDT objects.
#'
#' @param x \code{\link{MDT}} object.
#' @param keys \code{character} or \code{\link{data.table}} with SampleID
#' and/or FeatureID key values to select records.
#' Can only be a \code{character} when \code{keytype} is not \code{info}.
#' @param columns The columns to be retrieved from the annotation data.table.
#' @param keytype \code{character} specifying which annotation data.table
#' should be used. Choose from \code{"annotations"}, \code{"phenotype"} or
#' \code{"info"}.
#' @param na.rm Should NA values be removed from results.
#'
#' @return \code{\link[data.table]{data.table}} with relevant key and annotations.
#'
#' @exportMethod selectMDT
setGeneric("selectMDT", function(x, keys, columns, keytype = "annotations", na.rm = TRUE) {
  standardGeneric("selectMDT")
})

## fillMDT
#' Fill missing values for MDT objects.
#'
#' Fill in missing values in \code{mtable} slot for \code{\link{MDT}} objects.
#' Identifies missing combinations of \code{SampleID} and \code{FeatureID} and
#' gives it value specified in \code{fill} argument. Maybe necessary after
#' applying \code{\link{filterMDT}}.
#'
#' @param x \code{\link{MDT}} object.
#' @param fill Value to replace missing values.
#'
#' @return Filled \code{\link{MDT}} object.
#'
#' @exportMethod fillMDT
#'
setGeneric("fillMDT", function(x, fill) {
  standardGeneric("fillMDT")
})

## asMatrixMDT
#' Convert mtable data.table to matrix
#'
#' Convert mtable \code{\link[data.table]{data.table}} to \code{matrix}
#'
#' @inheritParams data.table::dcast.data.table
#' @param x \code{\link{MDT}} object.
#'
#' @return Mutation \code{matrix}
#'
#' @exportMethod asMatrixMDT
#'
setGeneric("asMatrixMDT", function(x, fill = NA, fun.aggregate = NULL, ...) {
  standardGeneric("asMatrixMDT")
})

##
## 2. MLGWAS Main Methods ------------------------------------------------------
##

## bindML ----
#' @exportMethod bindML
setGeneric("bindML", function(...) {
  standardGeneric("bindML")
})

## trainClassifier ----
#' Train Classifier for MDT
#'
#' Trains a classifier that predicts phenotypic \code{response} for
#' \code{\link{MDT}} objects.
#'
#' @inheritParams caret::train
#' @inheritParams foreach::foreach
#'
#' @param mdt \code{\link{MDT}} object.
#' \code{FeatureID}s and \code{SampleID}s will correspond between the
#' \code{\link{MDT}} input and the \code{\link{MLGWAS}} output object.
#'
#' @param id Identification name of classifier: \code{ClassifierID}. By default,
#' uses the name of the method appended with '.permuted' if applicable.
#'
#' @param partitions Some \code{function} that takes a \code{\link{MDT}}
#' object as input and returns a \code{list} of \code{integer}s specifying
#' the training set indexes.
#' See \code{\link[caret]{createDataPartition}} for such functions.
#' Names of \code{list} will define the \code{PartitionID}.
#'
#' @param preprocess Some \code{function} that takes the training matrix,
#' with features as columns and samples as rows,
#' and the response (response(mdt)) as input and returns a
#' \code{function} that modifies a
#' a matrix with the same number of columns (features).
#' Can be used for feature selection or preprocessing.
#' The function will be created internally only using the training partition,
#' thus allowing to apply the feature selection to the test set without using
#' any of its information.
#'
#' @param permute \code{logical} specifying if \code{\link{MDT}} should be
#' permuted. Distribution of performance scores will correspond to a null
#' distribution and can be used to test the significance of the a classifier.
#'
#' @param validation \code{logical} specifying if test set should be data
#' samples not in \code{partition}.
#'
#' @param parallel \code{logical} specifying if parallel interface should be used.
#' Depends on \code{\link[foreach]{foreach}}.
#' See \url{https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf}
#'
#' @param verbose \code{logical} specifying if function should be run in verbose mode.
#'
#' @param mtable_vars Variables in \code{mtable} to be included in matrix.
#' Choose between \code{'none'}, \code{'all'}
#' or a subset of values in \code{features(x)}.
#'
#' @param phen_vars Variables in \code{phenotype} to be included in matrix.
#' Choose between \code{'none'}, \code{'all'} or a subset of values in
#' \code{colnames(phenotype(x))}. \code{SampleID} and \code{Response}
#'  are excluded.
#'
#' @param fill How should missing values in \code{\link{MDT}} be handled.
#' By default, does nothing: leaves then as \code{NA}.
#'
#' @param ... Further arguments to pass on to \code{\link[caret]{train}}.
#'
#' @return \code{\link{MLGWAS}} object.
#'
#' @seealso
#' \code{\link{trainClassifiers}}
#' \code{\link[caret]{train}}
#' \code{\link[caret]{createDataPartition}}
#'
#' @import caret
#'
#' @exportMethod trainClassifier
#'
setGeneric("trainClassifier", function(mdt,
                                       id,
                                       mtable_vars = "all",
                                       phen_vars = "none",
                                       fill = NA,
                                       partitions = caret::createMultiFolds(y = response(mdt), k = 10L, times = 1L),
                                       preprocess = function(x, y) function(x) x,
                                       method = "rf",
                                       permute = FALSE,
                                       validation = TRUE,
                                       parallel = FALSE,
                                       verbose = TRUE,
                                       .export = NULL,
                                       ...) {
  standardGeneric("trainClassifier")
})

## trainClassifiers ----
#' Calls trainClassifier for different combination of parameters.
#'
#' Calls trainClassifier for different combination of parameters as specified
#' by a \code{data.frame}.
#'
#' @inheritParams foreach::foreach
#'
#' @param mdt \code{\link{MDT}} object.
#'
#' @param params \code{data.frame} of different combinations of parameters.

#' See \code{\link[base]{expand.grid}} to see how to create such a data.frame.
#' To specify a parameter of a length bigger than one,
#' and for which combinations of values are not relevant
#' (such as \code{phen_vars}), specify it as a list.
#' For example, \code{expand.grid(p = 0.5, phen_vars = list('Sex', 'Age'))}.
#' Arguments with value NA will not be passed on to \code{\link{trainClassifier}}.
#'
#' @param rbind \code{logical} specifying if results should be bound at the end with
#' \code{rbind} and return just one \code{\link{MLGWAS}} object, or, if results
#' should be returned as a \code{list} of \code{\link{MLGWAS}} objects
#'
#' @param seed \code{integer} specifying random number generation seed.
#'
#' @param verbose \code{logical} specifying if function should be run in
#'  verbose mode.
#'
#' @return An \code{\link{MLGWAS}} object or a \code{list} of \code{\link{MLGWAS}}
#' objects according to \code{rbind} argument.
#'
#' @seealso \code{\link{trainClassifier}}
#'
#' @exportMethod trainClassifiers
#'
setGeneric("trainClassifiers", function(mdt,
                                        params,
                                        rbind = TRUE,
                                        verbose = TRUE,
                                        seed = NULL,
                                        .export = NULL) {
  standardGeneric("trainClassifiers")
})

##
## 3. Annotations --------------------------------------------------------------
##

## annotateStat ----
#' @rdname annotateMDT
#' @exportMethod annotateStat
setGeneric("annotateStat", function(x,
                                    function_list = list(
                                      FISHER = function(x, y) {
                                        stats::fisher.test(table(factor(x, c(0L, 1L, 2L)),
                                                                 factor(y, c(0L, 1L))))$p.value
                                      },
                                      HW = function(x, y) {
                                        HardyWeinberg::HWExact(table(factor(x, c(0L, 1L, 2L))),
                                                               verbose = FALSE)$pval
                                      }),
                                    fill = NA,
                                    verbose = TRUE) {
  standardGeneric("annotateStat")
})


## annotateLocation ----
#' @rdname annotateMDT
#' @exportMethod annotateLocation
setGeneric("annotateLocation", function(x,
                                        txdb,
                                        region = VariantAnnotation::AllVariants(),
                                        columns = c("LOCATION", "ENTREZID"),
                                        verbose = TRUE) {
  standardGeneric("annotateLocation")
})

## annotateBioconductor ----
#' @rdname annotateMDT
#' @exportMethod annotateBioconductor
setGeneric("annotateBioconductor", function(x,
                                            annotations_db,
                                            keys_colname,
                                            columns,
                                            keytype = keys_colname,
                                            verbose = TRUE) {
  standardGeneric("annotateBioconductor")
})

## annotateCoding ----
#' @rdname annotateMDT
#' @exportMethod annotateCoding
setGeneric("annotateCoding", function(x,
                                      seqSource,
                                      txdb,
                                      columns = c("CONSEQUENCE"),
                                      verbose = TRUE) {
  standardGeneric("annotateCoding")
})

#' @rdname annotateMDT
#' @exportMethod annotateMAF
setGeneric("annotateMAF", function(x) {
  standardGeneric("annotateMAF")
})

## annotateRSID ----
#' @rdname annotateMDT
#' @exportMethod annotateRSID
setGeneric("annotateRSID", function(x,
                                    snp_loc,
                                    verbose = TRUE) {
  standardGeneric("annotateRSID")
})

##
## 4. Plots --------------------------------------------------------------------
##

#' Plot MDT objects
#'
#' Various plots for \code{\link{MDT}} objects.
#'
#' @param x \code{\link{MDT}} object.
#'
#' @param type Type of plot
#' \itemize{
#' \item{scatter: }{Scatter plot using \code{\link[ggplot2]{geom_point}}}
#' \item{heatmap: }{Heatmap using \code{\link[gplots]{heatmap.2}}}
#' }
#'
#' @param preprocess Preprocessing transformation. Choose from 'scale',
#' 'center' and 'pca'.
#'
#' @param space For \code{y='scatter'}, specifies
#' if should plot feature space or sample space.
#'
#' @param dims For \code{y='scatter'}, specifies which features or samples
#' should be plotted. Can be integers specifying column index or
#' characters.
#'
#' @param fill Value to replace missing values in \code{mtable}.
#'
#' @param ... Further arguments to be passed on to
#' \code{\link[gplots]{heatmap.2}} or \code{\link[ggplot2]{geom_point}}
#'
#' @exportMethod plotMDT
#'
setGeneric("plotMDT", function(x,
                               type = c("scatter", "heatmap"),
                               preprocess = c("scale", "center", "pca"),
                               space = c("samples", "features"),
                               dims = c(1L, 2L),
                               fill = NA,
                               ...) {
  standardGeneric("plotMDT")
})

# manhattanPlot ----
#' Manhattan Plot for MDT and MLGWAS objects.
#'
#' Plot a Manhattan Plots for \code{\link{MDT}} and \code{\link{MLGWAS}}
#' objects. Can use any \code{numeric}
#' value in \code{annotations(x)} as specified by \code{value} argument.
#' Can also use feature importance as found in \code{featimp(mlgwas)} if
#' \code{value} is an \code{\link{MLGWAS}} object.
#'
#' @param x \code{\link{MDT}} object.
#'
#' @param group \code{character} specifying which column of \code{annotations(x)}
#' specifying grouping of 'x' axis such as chromosome number or
#' \code{\link{MLGWAS}} object.
#'
#' @param pos \code{character} specifying which column of \code{annotations(x)}
#' specifying position on 'x' axis.
#'
#' @param value \code{character} specifying which column of \code{annotations(x)}
#' should be used as values. For example, call \code{annotateStat} and then
#' use Fisher's exact test p-values. Can also be \code{\link{MLGWAS}}; plots
#' feature importance.
#'
#' @param f Some function that modifies a vector of integers.
#' For example, \code{-log10} for p-values.
#'
#' @param col Plot colors.
#'
#' @return A \code{\link{ggplot}} Manhattan plot.
#'
#' @exportMethod manhattanPlot
#'
setGeneric("manhattanPlot", function(x,
                                     value,
                                     group = "CHROM",
                                     pos = "POS",
                                     f = function(value) -log10(value),
                                     col = c("blue", "orange")) {
  standardGeneric("manhattanPlot")
})

# plotPerf ----
#' @rdname plotMLGWAS
#' @exportMethod plotPerf
setGeneric("plotPerf", function(x, measure, group) {
  standardGeneric("plotPerf")
})

# plotFeatimp ----
#' @rdname plotMLGWAS
#' @exportMethod plotFeatimp
setGeneric("plotFeatimp", function(x, sort_by_classifiers, features, max_features = 20) {
  standardGeneric("plotFeatimp")
})


# plotErrors ----
#' @rdname plotMLGWAS
#' @exportMethod plotErrors
setGeneric("plotErrors", function(x) {
  standardGeneric("plotErrors")
})

##
## 5. Summary ------------------------------------------------------------------
##

## summaryPerf ----
#' @rdname summaryMLGWAS
#' @exportMethod summaryPerf
setGeneric("summaryPerf", function(x,
                                   summarize = function(x) c(mean = mean(x, na.rm = TRUE),
                                                             sd = stats::sd(x, na.rm = TRUE),
                                                             len = length(stats::na.omit(x)))) {
  standardGeneric("summaryPerf")
})

## summaryFeatimp ----
#' @rdname  summaryMLGWAS
#' @exportMethod summaryFeatimp
setGeneric("summaryFeatimp", function(x,
                                      summarize = function(x) c(mean = mean(x, na.rm = TRUE),
                                                                sd = stats::sd(x, na.rm = TRUE),
                                                                len = length(stats::na.omit(x)))) {
  standardGeneric("summaryFeatimp")
})

## summaryErrors ----
#' @rdname summaryMLGWAS
#' @exportMethod summaryErrors
setGeneric("summaryErrors", function(x,
                                     summarize = function(x) c(mean = mean(x, na.rm = TRUE),
                                                               sd = stats::sd(x, na.rm = TRUE),
                                                               len = length(stats::na.omit(x)))) {
  standardGeneric("summaryErrors")
})

##
## 6. Setters and Getters --------------------------------------------------
##

## mtable ----
#' @exportMethod mtable
setGeneric("mtable", function(x) {
  standardGeneric("mtable")
})

#' @exportMethod mtable<-
setGeneric("mtable<-", function(x, value) {
  standardGeneric("mtable<-")
})

## annotations ----
#' @exportMethod annotations
setGeneric("annotations", function(x) {
  standardGeneric("annotations")
})

#' @exportMethod annotations<-
setGeneric("annotations<-", function(x, value) {
  standardGeneric("annotations<-")
})

## phenotype ----
#' @exportMethod phenotype
setGeneric("phenotype", function(x) {
  standardGeneric("phenotype")
})

#' @exportMethod phenotype<-
setGeneric("phenotype<-", function(x, value) {
  standardGeneric("phenotype<-")
})

## annotations ----
#' @exportMethod info
setGeneric("info", function(x) {
  standardGeneric("info")
})

#' @exportMethod info<-
setGeneric("info<-", function(x, value) {
  standardGeneric("info<-")
})

## features ----
#' @exportMethod features
setGeneric("features", function(x) {
  standardGeneric("features")
})

#' @exportMethod features<-
setGeneric("features<-", function(x, value) {
  standardGeneric("features<-")
})

## samples ----
#' @exportMethod samples
setGeneric("samples", function(x) {
  standardGeneric("samples")
})

#' @exportMethod samples<-
setGeneric("samples<-", function(x, value) {
  standardGeneric("samples<-")
})

## response ----
#' @exportMethod response
setGeneric("response", function(x) {
  standardGeneric("response")
})

#' @exportMethod response<-
setGeneric("response<-", function(x, value) {
  standardGeneric("response<-")
})

## performances ----
#' @exportMethod performances
setGeneric("performances", function(x) {
  standardGeneric("performances")
})

#' @exportMethod performances<-
setGeneric("performances<-", function(x, value) {
  standardGeneric("performances<-")
})

## clasifiers ----
#' @exportMethod classifiers
setGeneric("classifiers", function(x) {
  standardGeneric("classifiers")
})

#' @exportMethod classifiers<-
setGeneric("classifiers<-", function(x, value) {
  standardGeneric("classifiers<-")
})

## featimp ----
#' @exportMethod featimp
setGeneric("featimp", function(x) {
  standardGeneric("featimp")
})

#' @exportMethod featimp<-
setGeneric("featimp<-", function(x, value) {
  standardGeneric("featimp<-")
})

## partitions ----
#' @exportMethod partitions
setGeneric("partitions", function(x) {
  standardGeneric("partitions")
})

#' @exportMethod partitions<-
setGeneric("partitions<-", function(x, value) {
  standardGeneric("partitions<-")
})

## predictions ----
#' @exportMethod predictions
setGeneric("predictions", function(x) {
  standardGeneric("predictions")
})

#' @exportMethod predictions<-
setGeneric("predictions<-", function(x, value) {
  standardGeneric("predictions<-")
})

## errors ----
#' @exportMethod errors
setGeneric("errors", function(x) {
  standardGeneric("errors")
})

#' @exportMethod errors<-
setGeneric("errors<-", function(x, value) {
  standardGeneric("errors<-")
})

# classifierNames ----
#' @exportMethod classifierNames
setGeneric("classifierNames", function(x) {
  standardGeneric("classifierNames")
})

#' @exportMethod classifierNames<-
setGeneric("classifierNames<-", function(x, value) {
  standardGeneric("classifierNames<-")
})

# partitionNames ----
#' @exportMethod partitionNames
setGeneric("partitionNames", function(x) {
  standardGeneric("partitionNames")
})

#' @exportMethod partitionNames<-
setGeneric("partitionNames<-", function(x, value) {
  standardGeneric("partitionNames<-")
})
