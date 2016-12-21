################################################################################
## All Utilities
################################################################################

##
## 1. Exported Utilities -------------------------------------------------------
##

# chrAsCharacter ----
#' Chromsome integer values as character
#'
#' Chromsome \code{integer} values as \code{character}
#'
#' Chromosome codes follow MAP files codeing as in PLINK.
#' Autosomes are coded 1 through 22.
#' Others use the following code:
#' \itemize{
#' \item  X chromosome: X or 23
#' \item  Y chromosome: Y or 24
#' \item  Pseudo-autosomal region of X: XY or 25
#' \item  Mitochondrial: MT or 26
#' }
#'
#' @param x \code{integer} \code{vector} of chromosomes.
#'
#' @return \code{character} \code{vector} of chromosomes.
#'
#' @export
#'
chrAsCharacter <- function(x) {

  if (! is.integer(x)) {
    x <- tryCatch(as.integer(x),
                  error = function(e) {
                    stop("Invalid 'x' argument: must be a integer vector.")
                  },
                  warning = function(w) {
                    stop("Invalid 'x' argument: must be a integer vector.")
                  })
  }
  y <- as.character(x)
  y[x == 23L] <- "X"
  y[x == 24L] <- "Y"
  y[x == 25L] <- "XY"
  y[x == 26L] <- "MT"

  y
}

# chrAsNumeric ----
#' Chromsome character values as integer
#'
#' Chromsome \code{character} values as \code{integer}
#'
#' Chromosome codes follow MAP files codeing as in PLINK.
#' Autosomes are coded 1 through 22.
#' Others use the following code:
#' \itemize{
#' \item  X chromosome: X or 23
#' \item  Y chromosome: Y or 24
#' \item  Pseudo-autosomal region of X: XY or 25
#' \item  Mitochondrial: MT or 26
#' }
#'
#' @param x \code{character} \code{vector} of chromosomes.
#'
#' @return \code{integer} \code{vector} of chromosomes.
#'
#' @export
#'
chrAsNumeric <- function(x) {

  if (! is.character(x)) {
    x <- tryCatch(as.character(x),
                  error = function(e) {
                    stop("Invalid 'x' argument: must be a character vector.")
                  },
                  warning = function(w) {
                    stop("Invalid 'x' argument: must be a character vector.")
                  })
  }

  x <- base::gsub("^X$", replacement = "23", x, ignore.case = TRUE)
  x <- base::gsub("^Y$", replacement = "24", x, ignore.case = TRUE)
  x <- base::gsub("^XY$", replacement = "25", x, ignore.case = TRUE)
  x <- base::gsub("^M$", replacement = "26", x, ignore.case = TRUE)
  x <- base::gsub("^MT$", replacement = "26", x, ignore.case = TRUE)

  x <- tryCatch(as.integer(x),
                error = function(e) {
                  stop("Invalid 'x' argument: contained values character values that were not 'X', 'Y', 'XY', 'M' or 'MT'.")
                },
                warning = function(w) {
                  stop("Invalid 'x' argument: contained values character values that were not 'X', 'Y', 'XY', 'M' or 'MT'.")
                })
  x
}

# readVcfDir ----
#' @title Reads directory of VCF files
#'
#' @description Real all Variant Call Format (VCF) files contained in a
#' directory using the VariantAnnotation Bioconductor package.
#'
#' @inheritParams VariantAnnotation::readVcf
#'
#' @param dir Directory containing VCF files.
#'
#' @param grep_out Regular expression for file names to be filter out.
#'
#' @param param_function Function to be evaluated and as param argument.
#' Only takes as argument filename and must return a
#' \code{\link[VariantAnnotation]{ScanVcfParam}} object.
#'
#' @param n_files If specified, will only import the first
#' \code{n_files} in the directory.
#'
#' @param verbose \code{logical} specifying if function should be run in verbose mode.
#'
#' @param ... Further arguments to pass on to
#' \code{\link[VariantAnnotation]{readVcf}}
#'
#' @return A list of \code{\link[VariantAnnotation]{VCF}} objects.
#'
#' @seealso \code{\link[VariantAnnotation]{readVcf}}
#'
#' @details Uses \code{\link{list.files}} to identify files ending with ".vcf"'.
#'
#' @export
#'
readVcfDir <- function(dir,
                       param_function = function(fn) VariantAnnotation::ScanVcfParam(),
                       grep_out,
                       n_files,
                       verbose = TRUE,
                       ...) {

  mylapply <- .chooseLapply(verbose)

  file_list <- base::list.files(path = dir, pattern = "\\.vcf$") # List VCF files

  if (! missing(grep_out)) {
    file_list <- file_list[!grepl(grep_out, file_list)] # Grep out files
  }

  if (! missing(n_files)) {
    file_list <- file_list[1:n_files]
  }

  if (length(file_list) == 0) {
    stop("No VCF files were found.")
  }

  vcfs <- mylapply(file_list, function(fn) { # Read all VCF files
    param <- param_function(fn)
    VariantAnnotation::readVcf(file = paste(dir, fn, sep = ""),
                               param = param,
                               ...)
  })
  names(vcfs) <- base::unlist(lapply(file_list,
                                     function(x) strsplit(x, ".vcf")[[1]]))
  vcfs
}

##
## 2. Internal Utilities -------------------------------------------------------
##

# changeLevels ----
# Converts a vector (x) into a factor,
# changes its levels (value) and
# reconverts it back to the original vector class.
.changeLevels <- function(x, value) {
  class <- class(x)
  x <- base::as.factor(x)
  base::levels(x) <- value
  methods::as(x, class)
}

# chooseLapply ----
.chooseLapply <- function(verbose = TRUE, ...) {
  if (verbose) { pbapply::pblapply
  } else if (! verbose) { base::lapply
  } else { stop("Invalid 'verbose' argument.") }
}

# fixGT ----
# Makes GT column consistant throughout VCF files
# By consistent, we mean that numeric values in GT should refer to the same
# ALT nucleotide(s)
.fixGT <- function(featureid, bad_alt_df) {

  ## Retrieve data
  x <- bad_alt_df[featureid]
  alt <- x$ALT # alternative allele
  gt <- x$GT # genotype

  ## Format data
  strsplit_alt <- base::strsplit(alt, ",")
  gt1 <- as.integer(base::substr(gt, 1L, 1L))
  sep <- base::substr(gt, 2L, 2L)
  gt2 <- as.integer(base::substr(gt, 3L, 3L))

  ## Map GT to Nucleotides (NA means REF, A,C,G,T means ALT)
  genotype <- base::vapply(seq_len(length(alt)), function(j) {
    if (gt1[j] == 0L) { alt1 <- NA_character_
    } else { alt1 <- strsplit_alt[[j]][gt1[j]] }
    if (gt2[j] == 0L) { alt2 <- NA_character_
    } else { alt2 <- strsplit_alt[[j]][gt2[j]] }
    c(alt1, alt2)
  }, character(2))

  ## Create new ALT based on frequency of observed ALT
  new_alt <- names(base::sort(base::table(base::as.vector(genotype)),
                              decreasing = TRUE))

  ## Map new ALT to genotype
  genotype <- base::matrix(as.character(base::match(genotype, new_alt)),
                           nrow = 2, ncol = ncol(genotype))
  genotype[is.na(genotype)] <- "0"

  ## Create new ALT and GT fields
  x$ALT <- base::paste(new_alt, collapse = ",")
  x$GT <- base::unlist(base::lapply(seq_len(ncol(genotype)), function(j) {
    if (sep[j] == "|") {
      base::paste(genotype[, j, drop = TRUE], collapse = sep[j])
    } else {
      base::paste(base::sort(genotype[, j, drop = TRUE], decreasing = TRUE),
                  collapse = sep[j])
    }
  }))

  x
}

# gtToNum ----
# Converts GT column of VCF file to numeric values as in PLINK
.gtToNum <- function(gt) {

  # Convert to character vector
  gt <- as.character(unlist(gt))

  # Take into account phased genotypes
  gt <- base::gsub("|", "/", gt, fixed = TRUE)

  # Reorder genotypes (1/0 becomes 0/1)
  notna <- !is.na(gt)
  gtmat <- base::matrix(c(base::substr(gt[notna], 1L, 1L),
                          base::substr(gt[notna], 3L, 3L)),
                        2, sum(notna), TRUE)
  gt[notna] <- base::apply(gtmat, 2, function(x) base::paste(sort(x), collapse="/"))

  # Most common genotypes ("0/0", "0/1", "1/1") get 0L, 1L and 2L
  common <- c("0/0", "0/1", "1/1")
  # common <- common[!is.na(match(common, gt))] # check if they are in the data

  # Weird ones (e.g. 2/3) get values in function of frequency
  weird <- names(base::sort(base::table(gt[! gt %in% common]), decreasing = TRUE))

  # Convert gt character to integers
  lev <- c(common, weird)
  gt[!is.na(gt)] <- seq(0, length(lev)-1)[base::match(gt[!is.na(gt)], lev)]
  gt <- as.integer(gt)
}

# makeXMatrix ----
# Creates a data.frame combining features in mtable and phenotype
.makeXDf <- function(x, mtable_vars = "all", phen_vars = "none", fill = NA) {

  if (mtable_vars[1] == "none" & phen_vars[1] == "none") {
    base::stop("Invalid combination of 'mtable_vars' and 'phen_vars'arguments: both can not be set to 'none'.")
  }

  ## Subset mtable
  mtable <- mtable(x)
  if (mtable_vars[1] == "none") {
    mtable_vars <- base::character(0)
  } else if (mtable_vars[1] == "all") {
    mtable_vars <- MachineLearningGWAS::features(x)
  } else if ((any(! mtable_vars %in% MachineLearningGWAS::features(x)))) {
    base::stop("Invalid 'mtable_vars' argument: does not match values in features(x).")
  }
  if (length(mtable_vars) != 0) {
    mtable <- mtable[mtable$FeatureID %in% base::unique(stats::na.omit(mtable_vars)), ]
    mtable <- data.table::dcast.data.table(mtable,
                                           SampleID ~ FeatureID,
                                           value.var = "VALUE")
  } else {
    mtable <- data.table::data.table(SampleID = MachineLearningGWAS::samples(x))
  }

  ## Subset phenotype
  phenotype <- MachineLearningGWAS::phenotype(x)

  if (phen_vars[1] == "none") {
    phen_vars <- base::character(0)
  } else if (phen_vars[1] == "all") {
    phen_vars <- base::colnames(phenotype)[base::colnames(phenotype) != "RESPONSE"]
  } else if ((any(! phen_vars %in% base::colnames(phenotype)))) {
    stop("Invalid 'phen_vars' argument: does not match colnames of phenotype(x).")
  } else if (any(phen_vars == "RESPONSE")) {
    stop("Invalid 'phen_vars' argument: can not use RESPONSE.")
  }
  if (length(phen_vars) != 0) {
    phenotype <- phenotype[, base::unique(stats::na.omit(c("SampleID", phen_vars))),
                           with = FALSE]
  } else {
    phenotype <- data.table::data.table(SampleID = MachineLearningGWAS::samples(x))
  }

  ## Combine mtable and phenotype
  df <- merge(mtable, phenotype, all = TRUE, by = "SampleID")
  df[, "SampleID":=NULL]

  ## Replace missing values
  fill <- fill[1]
  if (! is.null(fill) && ! is.na(fill)) {
    df[is.na(df)] <- fill
  }

  # Convert to data.frame and check if not empty
  df <- base::as.data.frame(df)

  if (base::nrow(df) == 0 | base::ncol(df) == 0) {
    stop("Predictor matrix is empty.")
  }

  df
}

# manhattanPlotSupport ----
# Support function for manhattanPlot methods
# df has three columns: Group, Position and Value
# Returns a ggplot2 Manhattan Plot
.manhattanPlotSupport <- function(df, col) {

  # Order data
  df$Group <- base::tryCatch(as.numeric(df$Group), error = function(e) df$Group)
  df$Position <- base::tryCatch(as.numeric(df$Position),
                                error = function(e) df$Position)
  df <- df[order(df$Group, df$Position), ]

  ## Make sure values are correct data.type
  df$Group <- base::as.factor(df$Group)
  df$Position <- 1L:nrow(df)
  df$Value <- as.numeric(df$Value)

  ## Remove missing data
  df <- stats::na.omit(df)

  ## Determine where ticks should be set in ggplot
  ticks <- base::tapply(df$Position, df$Group, function(x) round(base::mean(x)))

  ## Define colours
  col <- rep(col, length.out = length(levels(df$Group)))
  names(col) <- base::levels(df$Group)
  col <- base::unlist(base::lapply(base::levels(df$Group), function(x) {
    rep(col[as.character(x)], sum(df$Group == x))
  }))

  ## Plot using ggplot
  df <- as.data.frame(df)

  ggplot2::ggplot(df, ggplot2::aes_string(x = "Position",
                                          y = "Value",
                                          color = "Group")) +
    ggplot2::geom_point(colour = col) +
    ggplot2::scale_x_continuous(breaks = ticks, labels = names(ticks)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   axis.line = ggplot2::element_line(colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = .1,
                                                              color = "black"))
}

# mergeDT ----
# Fast merge for data.tables. Automatically converts to
# data.table, and removes duplicated rows and columns.
.mergeDT <- function(x, y, by, all.x = TRUE,
                     all.y = FALSE, sort = TRUE, allow.cartesian = TRUE) {

  ## Convert to data.table
  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)

  ## Get only unique rows
  data.table::setkeyv(x, NULL)
  data.table::setkeyv(y, NULL)
  x <- base::unique(x)
  y <- base::unique(y)

  ## Set key for fast merge
  data.table::setkeyv(x, by)
  data.table::setkeyv(y, by)

  ## Remove duplicated columns in x
  if (any(colnames(x)[colnames(x) != by] %in% colnames(y))) {
    x[, colnames(x)[colnames(x) != by & colnames(x) %in% colnames(y)]:=NULL]
  }

  ## Merge
  merge(x, y,
        by = by,
        all.x = all.x, all.y = all.y,
        sort = sort, allow.cartesian = allow.cartesian)
}

# parseGenomicLocations ----
# Converts GRanges object to data.table because data.tables are faster
.parseGenomicLocations <- function(gloc) {
  # Handle REF and ALT alleles: make list of characters for ALT
  gloc$REF <- as.character(gloc$REF)
  gloc$ALT <- S4Vectors::unstrsplit(IRanges::CharacterList(gloc$ALT), sep=",")

  # Convert to data.table and format chromsome names and positions into integers
  gloc <- data.table::as.data.table(gloc)
  gloc[, c("end", "width"):=NULL]
  base::colnames(gloc)[1:3] <- c("CHROM", "POS", "STRAND")

  gloc$CHROM <- base::gsub("^X$", replacement = "23", gloc$CHROM, ignore.case = TRUE)
  gloc$CHROM <- base::gsub("^Y$", replacement = "24", gloc$CHROM, ignore.case = TRUE)
  gloc$CHROM <- base::gsub("^XY$", replacement = "25", gloc$CHROM, ignore.case = TRUE)
  gloc$CHROM <- base::gsub("^M$", replacement = "26", gloc$CHROM, ignore.case = TRUE)
  gloc$CHROM <- base::gsub("^MT$", replacement = "26", gloc$CHROM, ignore.case = TRUE)
  gloc$CHROM <- as.integer(gloc$CHROM)
  gloc$CHROM[is.na(gloc$CHROM)] <- 99L

  gloc$POS <- as.integer(gloc$POS)

  gloc
}

# pn ----
.pn <- function(...) {
  base::paste(..., sep = "\n")
}

# printMLGWAS ----
# Support function printting MLGWAS objects to console
.printMLGWAS <- function(number_of_items, item_name, items, n = 5L) {
  base::cat(base::paste0(number_of_items, " ", item_name, ": ",
                         base::paste(utils::head(items, n = n), collapse = ", ")))
  if (number_of_items > n) base::cat(", ...\n")
  else base::cat("\n")
}

# reverseList ----
# For a list with the following structure
# [[1]]$a [[1]]$b [[2]]$a [[2]]$b ...,
# reverses index structure for it to have the following structure
# $a[[1]] $a[[2]] $b[[1]] $b[[2]].
.reverseList <- function(x) {
  in_names <- names(x[[1]])
  flat <- base::lapply(in_names, function(i) base::lapply(x, `[[`, i))
  names(flat) <- in_names
  flat
}

# trainClassifierSupport ----
# Arguments are the same in trainClassifier
.trainClassifierSupport <- function(dp, x, y, permute,
                                    preprocess, method, validation,
                                    mtable_vars, phen_vars, fill, ...) {

  ## Value to be returned if anything fails
  failure <- list(clas = NULL, perf = NULL, featimp = NULL, pred = NULL)

  ## Permute MDT
  if (permute) y <- base::sample(y)

  ## Partition MDT and y into test and train
  train_x <- x[dp, ]
  train_y <- y[dp]

  if (validation) {
    test_x <- x[-dp, ]
    test_y <- y[-dp]
  } else {
    test_x <- x
    test_y <- y
  }

  if (base::nrow(train_x) == 0L |
      base::length(train_y) == 0L |
      base::nrow(test_x) == 0L |
      base::length(test_y) == 0L) {
    warning("Partitioning of led to data without samples.")
    return(failure)
  }

  ## Create and apply preprocessing function
  pp <- tryCatch(preprocess(train_x, train_y),
  error = function(e) {
    warning(paste("Preprocessing failed with error:", e))
    NULL
  })
  if (is.null(pp)) return(failure)

  train_x <- tryCatch(pp(train_x),
           error = function(e) {
             warning(paste("Preprocessing could not be applied to training set; failed with error:", e))
             NULL
           })
  if (is.null(train_x)) return(failure)


  test_x <- tryCatch(pp(test_x),
           error = function(e) {
             warning(paste("Preprocessing could not be applied to test set; failed with error:", e))
             NULL
           })
  if (is.null(test_x)) return(failure)

  ## Fit classifier
  clas <- base::tryCatch(caret::train(x = train_x, y = train_y, method = method, ...),
                         error = function(e) {
                           warning(paste("Training failed with error:", e))
                           NULL
                         })
  if (is.null(clas)) {
    return(failure)
  }


  ## Compute variable importance
  featimp <- base::tryCatch(caret::varImp(clas)$importance,
                            error = function(e) {
                              warning(paste("Variable importance computation failed with error:", e))
                              NULL
                            })
  if (is.null(featimp)) return(failure)

  featimp <- data.table::data.table(FeatureID = as.character(rownames(featimp)),
                                    Importance = featimp[, 1])

  ## Estimate performance
  test_y_hat <- tryCatch(stats::predict(clas, test_x),
                         error = function(e) {
                           warning(paste("Predicting on test set failed with error:", e))
                           NULL
                         })
  if (is.null(test_y_hat)) return(failure)

  if (is.factor(y)) {

    # if (all(dp == c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 14, 16, 17, 18, 20, 22, 24, 27, 28, 29, 31, 33, 35, 38, 39, 40))) browser()

    acc <- sum(test_y_hat == test_y) / length(test_y)

    kappa <- tryCatch(psych::cohen.kappa(data.frame(test_y_hat, test_y))$kappa,
                       error = function(e) {
                         warning(e)
                         NA_real_
                       })


    sen <- tryCatch(caret::sensitivity(test_y_hat, test_y),
                error = function(e) {
                  warning(e)
                  NA_real_
                })

    spe <- tryCatch(caret::specificity(test_y_hat, test_y),
                error = function(e) {
                  warning(e)
                  NA_real_
                })

    ppv <- tryCatch(caret::posPredValue(test_y_hat, test_y),
                error = function(e) {
                  warning(e)
                  NA_real_
                })

    npv <- tryCatch(caret::negPredValue(test_y_hat, test_y),
                error = function(e) {
                  warning(e)
                  NA_real_
                })

    f1s <- tryCatch((2 * ppv * sen) / (ppv + sen),
                error = function(e) {
                  warning(e)
                  NA_real_
                })


    # Replace NaN by NA
    acc[is.nan(acc)] <- NA_real_
    kappa[is.nan(kappa)] <- NA_real_

    sen[is.nan(sen)] <- NA_real_
    spe[is.nan(spe)] <- NA_real_

    ppv[is.nan(ppv)] <- NA_real_
    npv[is.nan(npv)] <- NA_real_

    f1s[is.nan(f1s)] <- NA_real_

    perf <- data.table::data.table(Accuracy = acc,
                                   Kappa = kappa,
                                   Sensitivity = sen,
                                   Specificity = spe,
                                   PPV = ppv,
                                   NPV = npv,
                                   F1Score = f1s)
  } else {

    ssres <- sum( (test_y - test_y_hat)**2 ) # sum of squared residuals
    sstot <- sum( (test_y - mean(test_y))**2 ) # total sum of squares
    rmse <- sqrt(sstot / length(test_y))
    rsquared = 1 - ssres/sstot

    # Replace NaN by NA
    rmse[is.nan(rmse)] <- NA_real_
    rsquared[is.nan(rsquared)] <- NA_real_

    perf <- data.table::data.table(RMSE = rmse, Rsquared = rsquared)
  }

  list(clas = clas,
       perf = perf,
       featimp = featimp,
       pred = test_y_hat)
}
