################################################################################
## Machine Learning Genome-Wide Association Study (MLGWAS) Methods
################################################################################

##
## 1. Constructor --------------------------------------------------------------
##

MLGWAS <- function(performances,
                   featimp,
                   errors,
                   classifiers,
                   partitions,
                   predictions,
                   call) {

  methods::new("MLGWAS",
               performances = performances,
               classifiers = classifiers,
               featimp = featimp,
               partitions = partitions,
               predictions = predictions,
               errors = errors,
               call = call)
}

##
## 2. Main Methods -------------------------------------------------------------
##

## trainClassifier ----
#' @rdname trainClassifier
setMethod("trainClassifier", "MDT", function(mdt,
                                             id,
                                             mtable_vars,
                                             phen_vars,
                                             fill,
                                             partitions,
                                             preprocess,
                                             method,
                                             permute,
                                             validation,
                                             parallel,
                                             verbose,
                                             .export,
                                             ...) {


  # Add ClassifierID if missing
  if (missing(id) || is.null(id) || is.na(id) || length(id) == 0) {
    if (permute) id <- base::paste0(method, ".permuted")
    else id <- method
  }

  ## Check arguments types
  id <- tryCatch(as.character(id)[1],
                 error = function(e) {
                   warning(e)
                   stop("Invalid 'id' argument: must be convertible to a character")
                 })

  mtable_vars <- tryCatch(as.character(mtable_vars),
                          error = function(e) {
                            warning(e)
                            stop("Invalid 'mtable_vars' argument: must be convertible to a character")
                          })

  phen_vars <- tryCatch(as.character(phen_vars),
                        error = function(e) {
                          warning(e)
                          stop("Invalid 'phen_vars' argument: must be convertible to a character")
                        })

  partitions <- tryCatch(as.list(partitions),
                         error = function(e) {
                           warning(e)
                           stop("Invalid 'partitions' argument: must be convertible to a list")
                         })

  method <- tryCatch(as.character(method)[1],
                     error = function(e) {
                       warning(e)
                       stop("Invalid 'method' argument: must be convertible to a character")
                     })

  permute <- tryCatch(as.logical(permute)[1],
                      error = function(e) {
                        warning(e)
                        stop("Invalid 'permute' argument: must be convertible to a logical")
                      })

  validation <- tryCatch(as.logical(validation)[1],
                         error = function(e) {
                           warning(e)
                           stop("Invalid 'validation' argument: must be convertible to a logical")
                         })

  parallel <- tryCatch(as.logical(parallel)[1],
                       error = function(e) {
                         warning(e)
                         stop("Invalid 'parallel' argument: must be convertible to a logical")
                       })

  verbose <- tryCatch(as.logical(verbose)[1],
                      error = function(e) {
                        warning(e)
                        stop("Invalid 'verbose' argument: must be convertible to a logical")
                      })

  # Check if response is present
  if (nrow(MachineLearningGWAS::phenotype(mdt)) == 0) {
    stop("Invalid 'mdt' argument: no phenotype data. Please call importPhenotype.")
  }

  ## Converts character response to factor
  if (class(MachineLearningGWAS::response(mdt)) == "character") {
    warning("response(mdt) was 'character'; converting to 'factor'")
    MachineLearningGWAS::response(mdt) <- as.factor(MachineLearningGWAS::response(mdt))
  }

  if (! method %in% names(caret::getModelInfo())) {
    stop("Invalid 'method' argument: must be in names(getModelInfo())")
  }

  ## Make partitions
  if (is.null(names(partitions))) {
    names(partitions) <- base::make.names(1:length(partitions))
  }
  ## Make x and y matrix and vector for training
  x <- .makeXDf(mdt, mtable_vars, phen_vars, fill)
  y <- MachineLearningGWAS::response(mdt)

  if (base::nrow(x) != length(y)) {
    stop("For some reason, x (data) had a different number of samples than y (response).")
  }

  ## Train classifiers using bootstraped data
  if (verbose) cat("Training classifiers\n")
  if (parallel) { doer <- foreach::`%dopar%`
  } else { doer <- foreach::`%do%` }
  dp <- NULL # CRAN Compatibility


  # partitions <- list(a = c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 14, 16, 17, 18, 20, 22, 24, 27, 28, 29, 31, 33, 35, 38, 39, 40), b = c(2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 30, 31, 33, 36, 37, 38, 39))
  # partitions <- list(fail = c(1), success = c(2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 30, 31, 33, 36, 37, 38, 39))

  results <- doer(foreach::foreach(dp = partitions,
                                   .export = .export,
                                   .verbose = verbose,
                                   .packages = "MachineLearningGWAS"),
                  .trainClassifierSupport(dp = dp,
                                          x = x,
                                          y = y,
                                          permute = permute,
                                          preprocess = preprocess,
                                          method = method,
                                          validation = validation,
                                          mtable_vars = mtable_vars,
                                          phen_vars = phen_vars,
                                          fill = fill, ...))
  names(results) <- names(partitions)

  ## Handle failed parititons
  failed <- base::unlist(base::lapply(results, function(x) is.null(x[[1]])))
  failed[is.null(failed)] <- TRUE
  failed[is.na(failed)] <- TRUE

  if (any(failed)) {

    results <- results[! failed]

    if (length(results) == 0) {
      stop("Training failed for all.")
    }

    warning(base::paste("The following Partitions failed:",
                        base::paste(names(partitions)[failed], collapse = ", ")))

    partitions <- partitions[! failed]

  }

  if (verbose) base::cat("\nFormating results\n")

  results <- .reverseList(results)

  ## Define IDs
  classifier_id <- as.character(id)
  partition_id <- as.character(names(partitions))

  ## Format performances into data.table
  performances <- data.table::rbindlist(results$perf)
  performances <- cbind(ClassifierID = classifier_id,
                        PartitionID = partition_id,
                        performances)
  data.table::setkeyv(performances, c("ClassifierID", "PartitionID"))

  ## Format featimp: transform list of featimp to cast data.table
  featimp <- lapply(seq_along(results$featimp), function(i) {
    data.table::data.table(PartitionID = partition_id[i], results$featimp[[i]])
  })
  featimp <- data.table::rbindlist(featimp)
  featimp <- cbind(ClassifierID = classifier_id, featimp)
  data.table::setkeyv(featimp, c("ClassifierID", "PartitionID", "FeatureID"))

  ## Format errors
  errors <- base::lapply(seq_along(partitions[! failed]), function(i) {

    if (validation) {
      in_testset <- ! 1:methods::selectMethod("nrow", "MDT")(mdt) %in% partitions[[i]]
    } else {
      in_testset <- rep(TRUE, methods::selectMethod("nrow", "MDT")(mdt))
    }

    if (is.factor(response(mdt))) {
      error <- rep(NA_integer_, methods::selectMethod("nrow", "MDT")(mdt))
      error[in_testset] <- response(mdt)[in_testset] != results$pred[[i]]
    } else {
      error <- rep(NA_real_, methods::selectMethod("nrow", "MDT")(mdt))
      error[in_testset] <- (results$pred[[i]] - response(mdt)[in_testset])**2
    }

    error
  })
  errors <- data.table::as.data.table(errors)
  colnames(errors) <- partition_id
  errors$SampleID <- MachineLearningGWAS::samples(mdt)
  errors <- data.table::melt.data.table(errors,
                                        measure.vars = partition_id,
                                        id.vars = "SampleID",
                                        variable.name = "PartitionID",
                                        value.name = "Error")
  errors$ClassifierID <- classifier_id
  errors$SampleID <- as.character(errors$SampleID)
  errors$PartitionID <- as.character(errors$PartitionID)
  errors$Error <- as.numeric(errors$Error)
  errors <- errors[, c("ClassifierID", "SampleID", "PartitionID", "Error"), with = FALSE]
  data.table::setkeyv(errors, c("ClassifierID", "SampleID", "PartitionID"))

  ## Format classifiers into list of lists
  classifiers <- list(results$clas)
  names(classifiers) <- classifier_id

  ## Format partitions
  partitions <- list(partitions)
  names(partitions) <- classifier_id

  ## Format predictions
  predictions <- list(results$pred)
  names(predictions) <- classifier_id

  ## Format call
  call <- list(match.call())
  names(call) <- classifier_id

  ## Return MLGWAS object
  MLGWAS(performances = performances,
         featimp = featimp,
         errors = errors,
         classifiers = classifiers,
         partitions = partitions,
         predictions = predictions,
         call = call)
})

## trainClassifiers ----
#' @rdname trainClassifiers
setMethod("trainClassifiers", c(mdt = "MDT", params = "data.frame"), function(mdt,
                                                                              params,
                                                                              rbind,
                                                                              verbose,
                                                                              seed,
                                                                              .export) {

  if (is.null(base::colnames(params))) {
    stop("Invalid 'params' argument: must have colnames.")
  }

  # Make sure ClassifierID are present and unique
  if (! is.null(params$id)) {
    params$id <- base::make.names(params$id, unique = TRUE)

  } else {
    params$id <- base::paste0("Row", 1:nrow(params))
  }

  if (! is.null(seed)) {
    seed <- tryCatch(as.integer(seed[1]), error = function(e) {
      stop("Invalid 'seed' argument: must be integer(1).")
    })
  }

  ## Convert params data.frame to list
  params_list <- base::lapply(1:base::nrow(params), function(i) {

    # Convert row to list
    p <- base::lapply(params[i, ], function(x) x)

    # Unlist lists in data.frame
    for (j in names(p)) {

      q <- p[[j]]

      # value is list, function or data.frame
      if (is.list(q) | is.function(q[[1]]) | is.data.frame(q[[1]])) {
        p[[j]] <- q[[1]]

        # value if factor
      } else if ( is.factor(q) ) {
        p[[j]] <- as.character(q)

      } else if (is.na(q)) {
        p[[j]] <- NULL

        # value is not supported
      } else if ( ! ( is.list(q) | is.function(q) | is.atomic(q) ) | isS4(q) ) {
        if ( is.na(q[[1]]) ) p[[j]] <- NULL
      }

      # value is NULL, omit
      if ( is.null(q) ) {
        p[[j]] <- NULL
      }
    }
    p
  })
  names(params_list) <- base::lapply(params_list, "[[", "id")

  ## Run trainClassifier for every combination
  mylapply <- .chooseLapply(verbose)

  results <- mylapply(params_list, function(p) {
    p$mdt <- mdt
    if (! is.null(seed)) set.seed(seed)
    tryCatch(do.call(MachineLearningGWAS::trainClassifier, p),
             error = function(e) {
               warning(e)
               return(NULL)
             })
  })

  ## Remove failed runs
  success <- base::sapply(seq_along(results), function(i) {
    if (! class(results[[i]]) == "MLGWAS") {
      warning(.pn(paste("Row", i, "failed.")))
      FALSE
    } else {
      TRUE
    }
  })
  results <- results[success]

  ## Bind results if necessary and return it
  if (rbind) {
    results <- Reduce(MachineLearningGWAS::bindML, results)
  }

  results
})

##
## 3. Plots --------------------------------------------------------------------
##

## plotPerf ----
#' @rdname plotMLGWAS
setMethod("plotPerf", "MLGWAS", function(x,
                                         measure,
                                         group) {

  ## Retrieve data.
  df <- as.data.frame(MachineLearningGWAS::performances(x))

  ## Check 'measure' arguments
  if (missing(measure)) {
    measure <- base::colnames(df)[! base::colnames(df) %in% c("ClassifierID", "PartitionID")][1]
  }

  measure <- measure[1]
  if (! measure %in% base::colnames(df)[3:length(df)]) {
    stop(base::paste("Invalid 'measure' argument: choose from:",
                     base::paste(colnames(df)[3:length(df)], collapse = ", ")))
  }

  ## Format data.frame for ggplot
  df <- stats::na.omit(df[, c("ClassifierID", measure)])
  base::colnames(df)[base::colnames(df) == measure] <- "Value"

  df$ClassifierID <- as.character(df$ClassifierID)
  df$Value <- as.numeric(df$Value)

  ## Add groups
  if (! missing(group)) {
    groupdf <- base::Reduce(rbind, base::lapply(seq_along(group), function(i) {
      data.frame(ClassifierID = base::make.names(group[[i]]), Group = names(group)[i])
    }))
    df <- merge(df, groupdf, by = "ClassifierID", all.x = TRUE, all.y = FALSE)
    df <- df[base::order(df$Group), ]
    df$Group <- base::factor(df$Group, levels = base::unique(df$Group), ordered = TRUE)
    df$ClassifierID <- base::factor(df$ClassifierID, levels = base::unique(df$ClassifierID), ordered = TRUE)

    gp <- ggplot2::ggplot(data = df, ggplot2::aes_string(x = "ClassifierID", y = "Value")) +
      ggplot2::geom_boxplot(ggplot2::aes_string(fill = "Group"), outlier.colour = NA)
  } else {
    df <- df[order(df$ClassifierID), ]
    df$ClassifierID <- factor(df$ClassifierID, ordered = TRUE)
    gp <- ggplot2::ggplot(data = df, ggplot2::aes_string(x = "ClassifierID", y = "Value")) +
      ggplot2::geom_boxplot(outlier.colour = NA, fill = NA)
  }

  gp <- gp + ggplot2::geom_jitter(width = 0.2, height = 0.0) +
    ggplot2::ylab(measure) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  gp
})

## manhattanPlot ----
#' @rdname manhattanPlot
setMethod("manhattanPlot", c("MDT", "MLGWAS"), function(x,
                                                        value,
                                                        group,
                                                        pos,
                                                        f,
                                                        col) {


  group <- group[1]
  pos <- pos[1]

  if (any(! c(group, pos, value) %in% base::colnames(annotations(x)))) {
    stop("Invalid 'group' or 'pos' argument: does not match any column in annotations(x)")
  }

  dfmdt <- MachineLearningGWAS::annotations(x)[, c("FeatureID", group, pos), with = FALSE]
  base::colnames(dfmdt) <- c("FeatureID", "Group", "Position")

  if (! as.integer(dfmdt$Position)) {
    stop("Invalid 'pos' argument: must correspond to a numeric column in annotations(x)")
  }

  if (methods::selectMethod("nrow", "MLGWAS")(value) > 1) {
    id <- MachineLearningGWAS::classifierNames(value)[1]
    warning(paste0("More than one 'ClassifierID' found. Using ", id, ".\nPlease use `[` to subset MLGWAS object."))
    value <- value[id]
  }

  # Create mean values
  dfml <- base::tapply(MachineLearningGWAS::featimp(value)$Importance,
                       MachineLearningGWAS::featimp(value)$FeatureID,
                       mean, na.rm = TRUE)
  dfml <- data.table::data.table(FeatureID = names(dfml), Value = as.numeric(dfml))
  df <- .mergeDT(dfmdt, dfml, by = "FeatureID", all.x = TRUE, all.y = TRUE)
  base::rm(dfmdt, dfml)

  ## Fill in missing groups
  df$Group <- as.character(df$Group)
  df$Group[is.na(df$Group)] <- "99"

  # Fill in missing values
  df$Value[is.na(df$Value) | is.nan(df$Value) | is.infinite(df$Value)] <- 0

  if (! missing(f)) {
    df$Value <- f(df$Value)
  }

  .manhattanPlotSupport(df, col) +
    ggplot2::ylab("Variable Importance in Machine Learning (%)")
})

# plotErrors ----
#' @rdname plotMLGWAS
setMethod("plotErrors", c("MLGWAS"), function(x) {

  df <- MachineLearningGWAS::errors(x)
  df$ClassifierID <- as.factor(df$ClassifierID)
  df$SampleID <- as.factor(df$SampleID)
  df$Error <- as.numeric(df$Error)

  df <- stats::na.omit(df)

  gp <- ggplot2::ggplot(df, ggplot2::aes_string(x = "SampleID", y = "Error", color = "ClassifierID")) +
    ggplot2::geom_boxplot() +
    ggplot2::coord_flip() +
    ggplot2::ylab("Error per Sample")


  maxi <- max(df$Error, na.rm = TRUE)
  if (maxi) {
    gp <- gp + ggplot2::ylim(c(0, 1))
  } else {
    gp <- gp + ggplot2::ylim(c(0, maxi))
  }

  gp

})

# plotFeatimp ----
#' @rdname plotMLGWAS
setMethod("plotFeatimp", c("MLGWAS"), function(x, sort_by_classifiers, features, max_features) {

  ## Retrieve data
  df <- MachineLearningGWAS::featimp(x)
  if (! missing(features)) {
    if (! all(features %in% df$FeatureID)) {
      stop("Invalid 'features' argument: must correspond to FeatureID in featimp(x).")
    }
    df <- df[df$FeatureID %in% features, ]
  }

  ## Sort in order of mean variable importance
  if (missing(sort_by_classifiers)) {
    sort_by_classifiers <- classifierNames(x)
  }
  fi <- base::sort(base::tapply(unlist(df[ClassifierID %in% sort_by_classifiers, "Importance", with = FALSE]),
                                unlist(df[ClassifierID %in% sort_by_classifiers, "FeatureID", with = FALSE]),
                                stats::median, na.rm = TRUE),
                   decreasing = TRUE)
  fi <- names(fi)[1:min(length(fi), max_features)]

  ## Only retains high scoring features
  df <- df[FeatureID %in% fi, ]
  df$FeatureID <- base::factor(df$FeatureID, levels = fi, ordered = TRUE)

  ggplot2::ggplot(df, ggplot2::aes_string(x = "FeatureID",
                                          y = "Importance",
                                          color = "ClassifierID")) +
    ggplot2::geom_boxplot() +
    # ggplot2::geom_boxplot(outlier.colour = NA, fill = NA) +
    # ggplot2::geom_point(position = ggplot2::position_jitterdodge(dodge.width = 0.9)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

})

##
## 4. Summaries ----------------------------------------------------------------
##

## summaryPerf ----
#' @rdname summaryMLGWAS
setMethod("summaryPerf", "MLGWAS", function(x, summarize) {

  perf <- MachineLearningGWAS::performances(x)
  coln <- ! base::colnames(perf) %in% c("ClassifierID", "PartitionID")
  data.table::setkeyv(perf, "ClassifierID")

  sum_test <- summarize(1:3)
  names_summarize <- names(sum_test)
  if (is.null(names_summarize)) {
    names_summarize <- base::make.names(1:length(sum_test))
  }

  results <- base::lapply(MachineLearningGWAS::classifierNames(x), function(clasid) {
    subperf <- perf[clasid, coln, with = FALSE]
    r <- base::lapply(subperf, summarize)
    r <- data.table::as.data.table(r)
    data.table::data.table(ClassifierID = clasid, Summary = names_summarize, r)
  })

  data.table::rbindlist(results)
})

# summaryFeatimp ----
#' @rdname summaryMLGWAS
setMethod("summaryFeatimp", "MLGWAS", function(x, summarize) {

  featimp <- MachineLearningGWAS::featimp(x)
  data.table::setkeyv(featimp, "ClassifierID")

  sum_test <- summarize(1:3)
  names_summarize <- names(sum_test)
  if (is.null(names_summarize)) {
    names_summarize <- base::make.names(1:length(sum_test))
  }

  results <- base::lapply(MachineLearningGWAS::classifierNames(x), function(clasid) {
    subfeatimp <- featimp[clasid]
    r <- base::tapply(subfeatimp$Importance, subfeatimp$FeatureID,
                      summarize, simplify = FALSE)
    r <- data.table::data.table(ClassifierID = clasid,
                                FeatureID = names(r),
                                matrix(unlist(r),
                                       length(r),
                                       length(r[[1]]),
                                       byrow = TRUE))
    colnames(r) <- c("ClassifierID", "FeatureID", names_summarize)
    r
  })

  data.table::rbindlist(results)
})

# summaryErrors ----
#' @rdname summaryMLGWAS
setMethod("summaryErrors", "MLGWAS", function(x, summarize) {

  err <- MachineLearningGWAS::errors(x)
  data.table::setkeyv(err, "ClassifierID")

  sum_test <- summarize(1:3)
  names_summarize <- names(sum_test)
  if (is.null(names_summarize)) {
    names_summarize <- base::make.names(1:length(sum_test))
  }

  results <- base::lapply(MachineLearningGWAS::classifierNames(x), function(clasid) {
    suberr <- err[clasid]
    r <- base::tapply(suberr$Error, suberr$SampleID, summarize, simplify = FALSE)
    r <- data.table::data.table(ClassifierID = clasid,
                                SampleID = names(r),
                                base::matrix(unlist(r),
                                             length(r),
                                             length(r[[1]]),
                                             byrow = TRUE))
    base::colnames(r) <- c("ClassifierID", "SampleID", names_summarize)
    r
  })

  data.table::rbindlist(results)
})

##
## 5. Combine and Subset -------------------------------------------------------
##

# []
setMethod("[", c(x = "MLGWAS", i = "character"), function(x, i) {

  if (any(! i %in% MachineLearningGWAS::classifierNames(x))) {
    stop("Invalid 'i' argument: does not correspond to classifierNames(x).")
  }

  perf <- x@performances
  featimp <- x@featimp
  err <- x@errors
  clas <- x@classifiers
  part <- x@partitions
  pred <- x@predictions
  call <- x@call

  data.table::setkey(perf, "ClassifierID")
  perf <- perf[i, ]

  data.table::setkey(featimp, "ClassifierID")
  featimp <- featimp[i, ]
  data.table::setkeyv(featimp, c("ClassifierID", "PartitionID", "FeatureID"))

  data.table::setkey(err, "ClassifierID")
  err <- err[i, ]
  data.table::setkeyv(err, c("ClassifierID", "PartitionID", "SampleID"))

  clas <- clas[i]
  part <- part[i]
  pred <- pred[i]
  call <- call[i]

  MLGWAS(performances = perf,
         featimp = featimp,
         errors = err,
         classifiers = clas,
         partitions = part,
         predictions = pred,
         call = call)
})

## bindML ----
#' @exportMethod bindML
setMethod("bindML", "MLGWAS", function(...) {
  objects <- list(...)

  ## Retrieve times and names
  times <- base::sapply(objects, function(o) base::nrow(o))
  names <- base::sapply(objects, function(o) MachineLearningGWAS::classifierNames(o))

  ## Check for duplicated ClassifierID
  if (length(base::unique(unlist(names))) != length(base::unlist(names))) {
    stop(.pn("Duplicated ClassifierIDs. Can not procede !",
             "Please change 'ClassifierID' using classifierNames(mlgwas) <- new_name"))
  }

  ## Bind slots
  perf <- data.table::as.data.table(do.call(rbind, lapply(objects, performances)))
  featimp <- data.table::as.data.table(do.call(rbind, lapply(objects, featimp)))
  err <- data.table::as.data.table(do.call(rbind, lapply(objects, errors)))
  clas <- base::do.call(c, base::lapply(objects, classifiers))
  part <- base::do.call(c, base::lapply(objects, partitions))
  pred <- base::do.call(c, base::lapply(objects, predictions))
  call <- base::do.call(c, base::lapply(objects, methods::slot, "call"))

  MLGWAS(performances = perf,
         featimp = featimp,
         errors = err,
         classifiers = clas,
         partitions = part,
         predictions = pred,
         call = call)
})

##
## 5. General Methods ----------------------------------------------------------
##

## show ----
setMethod("show", "MLGWAS", function(object) {

  cat("--- MLGWAS object ---\n")
  .printMLGWAS(base::nrow(object), "classifiers", MachineLearningGWAS::classifierNames(object))
  .printMLGWAS(base::ncol(object), "partitions", MachineLearningGWAS::partitionNames(object))
  .printMLGWAS(length(MachineLearningGWAS::features(object)), "features", MachineLearningGWAS::features(object))
  .printMLGWAS(length(MachineLearningGWAS::samples(object)), "samples", MachineLearningGWAS::samples(object))

})

## length ----
setMethod("length", "MLGWAS", function(x) {
  methods::selectMethod("ncol", "MLGWAS")(x)
})

## dim ----
setMethod("dim", "MLGWAS", function(x) {
  c(methods::selectMethod("nrow", "MLGWAS")(x),
    methods::selectMethod("ncol", "MLGWAS")(x))
})

## nrow ----
setMethod("nrow", "MLGWAS", function(x) {
  length(base::unique(names(x@classifiers)))
})

## ncol ----
setMethod("ncol", "MLGWAS", function(x) {
  length(base::unique(MachineLearningGWAS::featimp(x)$PartitionID))
})

##
## 6. Setters and Getters ------------------------------------------------------
##

## performances ----
setMethod("performances", "MLGWAS", function(x) {
  x@performances
})

setReplaceMethod("performances", "MLGWAS", function(x, value) {
  MLGWAS(performances = value,
         featimp = x@featimp,
         errors = x@errors,
         classifiers = x@classifiers,
         partitions = x@partitions,
         predictions = x@predictions,
         call = x@call)
})

## featimp ----
setMethod("featimp", "MLGWAS", function(x) {
  x@featimp
})


setReplaceMethod("featimp", "MLGWAS", function(x, value) {
  MLGWAS(performances = x@performances,
         featimp = value,
         errors = x@errors,
         classifiers = x@classifiers,
         partitions = x@partitions,
         predictions = x@predictions,
         call = x@call)
})

## errors ----
setMethod("errors", "MLGWAS", function(x) {
  x@errors
})

setReplaceMethod("errors", "MLGWAS", function(x, value) {
  MLGWAS(performances = x@performances,
         featimp = x@featimp,
         errors = errors,
         classifiers = x@classifiers,
         partitions = x@partitions,
         predictions = x@predictions,
         call = x@call)
})

## classifiers ----
setMethod("classifiers", "MLGWAS", function(x) {
  x@classifiers
})


setReplaceMethod("classifiers", "MLGWAS", function(x, value) {
  MLGWAS(performances = x@performances,
         featimp = x@featimp,
         errors = x@errors,
         classifiers = value,
         partitions = x@partitions,
         predictions = x@predictions,
         call = x@call)
})

## partitions ----
setMethod("partitions", "MLGWAS", function(x) {
  x@partitions
})

setReplaceMethod("partitions", "MLGWAS", function(x, value) {
  MLGWAS(performances = x@performances,
         featimp = x@featimp,
         errors = x@errors,
         classifiers = x@classifiers,
         partitions = value,
         predictions = x@predictions,
         call = x@call)
})

## predictions ----
setMethod("predictions", "MLGWAS", function(x) {
  x@predictions
})

setReplaceMethod("predictions", "MLGWAS", function(x, value) {
  MLGWAS(performances = x@performances,
         featimp = x@featimp,
         errors = x@errors,
         classifiers = x@classifiers,
         partitions = x@partitions,
         predictions = value,
         call = x@call)
})

## classifierNames ----
setMethod("classifierNames", "MLGWAS", function(x) {
  names(x@classifiers)
})

setReplaceMethod("classifierNames", "MLGWAS", function(x, value) {

  value <- as.character(value)

  performances <- x@performances
  featimp <- x@featimp
  err <- x@errors
  classifiers <- x@classifiers
  partitions <- x@partitions
  predictions <- x@predictions
  call <- x@call

  performances$ClassifierID <- .changeLevels(performances$ClassifierID, value)
  featimp$ClassifierID <- .changeLevels(featimp$ClassifierID, value)
  err$ClassifierID <- .changeLevels(err$ClassifierID, value)

  names(classifiers) <- .changeLevels(names(classifiers), value)
  names(partitions) <- .changeLevels(names(partitions), value)
  names(predictions) <- .changeLevels(names(predictions), value)
  names(call) <- .changeLevels(names(call), value)

  MLGWAS(performances = performances,
         featimp = featimp,
         errors = err,
         classifiers = classifiers,
         partitions = partitions,
         predictions = predictions,
         call = call)
})

## partitionNames ----
setMethod("partitionNames", "MLGWAS", function(x) {
  base::unique(x@featimp$PartitionID)
})

setReplaceMethod("partitionNames", "MLGWAS", function(x, value) {

  value <- as.character(value)

  performances <- x@performances
  featimp <- x@featimp
  errors <- x@errors
  classifiers <- x@classifiers
  partitions <- x@partitions
  predictions <- x@predictions

  performances$PartitionID <- .changeLevels(performances$PartitionID, value)
  featimp$PartitionID <- .changeLevels(featimp$PartitionID, value)
  errors$PartitionID <- .changeLevels(errors$PartitionID, value)

  for (i in seq_along(classifiers)) {
    names(classifiers)[[i]] <- .changeLevels(names(classifiers)[[i]], value)
  }
  for (i in seq_along(partitions)) {
    names(partitions[[i]]) <- .changeLevels(names(partitions[[i]]), value)
  }
  for (i in seq_along(predictions)) {
    names(predictions[[i]]) <- .changeLevels(names(predictions[[i]]), value)
  }

  MLGWAS(performances = performances,
         featimp = featimp,
         errors = errors,
         classifiers = classifiers,
         partitions = partitions,
         predictions = predictions,
         call = x@call)
})

## features ----
setMethod("features", "MLGWAS", function(x) {
  base::unique(as.character(x@featimp$FeatureID))
})

setReplaceMethod("features", "MLGWAS", function(x, value) {
  value <- as.character(value)
  featimp <- x@featimp
  featimp$FeatureID <- .changeLevels(featimp$FeatureID, value)
  MLGWAS(performances = x@performances,
         featimp = featimp,
         errors = x@errors,
         classifiers = x@classifiers,
         partitions = x@partitions,
         predictions = x@predictions,
         call = x@call)
})

## samples ----
setMethod("samples", "MLGWAS", function(x) {
  base::unique(as.character(x@errors$SampleID))
})

setReplaceMethod("samples", "MLGWAS", function(x, value) {
  value <- as.character(value)
  x@errors$SampleID <- .changeLevels(MachineLearningGWAS::samples(x), value)
  MLGWAS(performances = x@performances,
         featimp = x@featimp,
         errors = errors,
         classifiers = x@classifiers,
         partitions = x@partitions,
         predictions = x@predictions,
         call = x@call)
})
