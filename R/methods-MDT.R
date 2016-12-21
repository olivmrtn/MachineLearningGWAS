## Mutation Data Table (MDT) Methods -------------------------------------------

##
## 1. Constructor --------------------------------------------------------------
##
MDT <- function(mtable,
                annotations = data.table::data.table(),
                phenotype = data.table::data.table(),
                info = data.table::data.table()) {


  mtable$FeatureID <- as.character(mtable$FeatureID)
  mtable$SampleID <- as.character(mtable$SampleID)
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))

  features <- base::unique(mtable$FeatureID)
  samples <- base::unique(mtable$SampleID)

  if (nrow(annotations) != 0) {
    annotations$FeatureID <- as.character(annotations$FeatureID)
    data.table::setkeyv(annotations, "FeatureID")
    annotations <- annotations[annotations$FeatureID %in% features, ]
  }
  if (nrow(phenotype) != 0) {
    phenotype$SampleID <- as.character(phenotype$SampleID)
    data.table::setkeyv(phenotype, "SampleID")
    phenotype <- phenotype[phenotype$SampleID %in% samples, ]
  }
  if (nrow(info) != 0) {
    info$FeatureID <- as.character(info$FeatureID)
    info$SampleID <- as.character(info$SampleID)
    data.table::setkeyv(info, c("FeatureID", "SampleID"))
    info <- info[info$FeatureID %in% features & info$SampleID %in% samples, ]
  }

  methods::new("MDT",
               mtable = mtable,
               annotations = annotations,
               phenotype = phenotype,
               info = info,
               features = features,
               samples = samples)
}

##
## 2. Main Methods -------------------------------------------------------------
##

# vcfsToMDT ----
#' @rdname vcfsToMDT
setMethod("vcfsToMDT", "list", function(vcfs,
                                        fill,
                                        verbose) {

  mylapply <- .chooseLapply(verbose)

  ## Check arguments
  classes <- base::sapply(vcfs, function(vcf) {
    ! class(vcf) %in% c("CollapsedVCF", "ExtendedVCF")
  })
  if (any(classes)) {
    stop(.pn("Invalid 'vcfs' argument: must be a list of VCF objects.",
             "The following objects were not more of class 'CollapsedVCF' or 'ExtendedVCF':",
             base::paste(names(vcfs)[classes], collapse = ", ")))
  }
  rm(classes)

  ## Format VCF fields into data.tables
  cat("Converting VCF into data.table\n")
  # Iterate through VCF files in vcfs list
  formated_vcfs <- mylapply(seq_along(vcfs), function(i) {

    # Subset VCF file from VCF list.
    # May contain one or more individual
    vcf <- vcfs[[i]]

    # Define SampleIDs
    sampleid <- as.character(VariantAnnotation::samples(VariantAnnotation::header(vcf)))
    ns <- length(sampleid)

    # Extract Genomic Locations (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
    gloc <- .parseGenomicLocations(SummarizedExperiment::rowRanges(vcf))

    if (base::ncol(gloc) == 0L) {
      stop(base::paste("No genomic locations in file", names(vcfs)[i], "with sample", sampleid))
    }

    # Extract INFO field
    info <- data.table::as.data.table(VariantAnnotation::info(vcf))
    if (base::ncol(info) > 0L) {
      base::colnames(info) <- base::paste0("INFO_", base::colnames(info))
    }

    # Column bind Genomic Locations and INFO
    if (base::ncol(info) > 0L) {
      df <- cbind(gloc, info)
    } else {
      df <- gloc
    }
    base::rm(gloc, info)

    # Replicate rows in function of the number of samples
    df <- data.table::data.table(SampleID = rep(sampleid, each = base::nrow(df)),
                                 data.table::rbindlist(replicate(ns, df, simplify = FALSE)))

    # Itereate through individuals in VCF for GENO field
    geno <- data.table::rbindlist(base::lapply(seq_len(dim(vcf)[2]), function(j) {
      data.table::as.data.table(base::as.list(VariantAnnotation::geno(vcf[, j])))
    }))

    if (base::ncol(geno) == 0L) {
      stop(base::paste("No genotype in file", names(vcfs)[i], "with sample", sampleid))
    }

    cbind(df, geno)
  })


  ## Bind formated VCFs together to form info slot of MDT object
  info <- data.table::rbindlist(formated_vcfs) # bind all infos
  colnames(info) <- base::make.names(colnames(info), unique = TRUE) # make names unique

  ## Create annotations data.table
  annotations <- base::unique(info[, c("CHROM", "POS", "STRAND"), with = FALSE])
  annotations <- annotations[order(annotations$CHROM, annotations$POS), ]
  annotations$FeatureID <- as.character(1L:nrow(annotations))

  ## Merge annotations and info for FeatureIDs
  data.table::setkeyv(annotations, c("CHROM", "POS", "STRAND"))
  data.table::setkeyv(info, c("CHROM", "POS", "STRAND"))
  info <- merge(info, annotations, by = c("CHROM", "POS", "STRAND"), all.x = TRUE)
  rm(formated_vcfs)

  ## Set keys for annotations and info
  data.table::setkeyv(annotations, c("FeatureID"))
  data.table::setkeyv(info, c("FeatureID", "SampleID"))

  ## Fix ALT and GT columns in info so they are consistent
  if (verbose) base::cat("\nMaking ALT and GT columns consistent throughout data.\n")
  df <- base::unique(info[, c("FeatureID", "SampleID", "REF", "ALT", "GT"), with = FALSE])
  data.table::setkeyv(df, c("FeatureID"))

  # Determine which ALT are not consistent
  v <- base::tapply(df$ALT, df$FeatureID, function(alt) any(alt != alt[1]))
  badalt <- as.character(names(v)[v])
  base::rm(v)

  # Fix ALT and GT
  if (length(badalt) > 0) {
    df[badalt] <- data.table::rbindlist(mylapply(badalt, fixGT, df[badalt]))
  }
  data.table::setkeyv(df, c("FeatureID", "SampleID"))
  rm(badalt)

  # Reorder GT (1/0 becomes 0/1; 1|0 stays 1|0)
  reorder <- base::substr(df$GT, 2L, 2L) == "/"
  gtmat <- base::matrix(c(substr(df$GT[reorder], 1L, 1L),
                          substr(df$GT[reorder], 3L, 3L)),
                        2, sum(reorder), TRUE)
  df$GT[reorder] <- base::apply(gtmat, 2, function(x) base::paste(base::sort(x), collapse = "/"))
  base::rm(reorder, gtmat)

  ## REF and ALT are given to annotatitons; GT to info
  info[, c("CHROM", "POS", "STRAND", "REF", "ALT", "GT"):=NULL]
  info <- merge(info,
                df[, c("FeatureID", "SampleID", "GT"), with = FALSE],
                by = c("FeatureID", "SampleID"))

  annotations <- merge(annotations,
                       unique(df[, c("FeatureID", "REF", "ALT"), with = FALSE]),
                       by = "FeatureID")
  rm(df)

  ## Check if Features were duplicated and remove them with a warning
  dup <- duplicated(info[, c("FeatureID", "SampleID", "GT"), with = FALSE])
  if (any(dup)) {
    base::cat("\nCertain samples had duplicate genotypes (GT).\nThose were removed but you may want to check your VCF files.\n")
    base::print(info[dup, c("SampleID", "FeatureID", "CHROM", "POS", "GT"), with = FALSE])
    info <- info[!dup, ]
  }
  base::rm(dup)

  ## Create mtable
  if (verbose) base::cat("\nCreating mutation data.table.")
  mtable <- data.table::CJ(FeatureID = unique(info$FeatureID),
                           SampleID = unique(info$SampleID))
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))

  mtable <- merge(mtable, info[, c("FeatureID", "SampleID", "GT"), with = FALSE],
                  by = c("FeatureID", "SampleID"),
                  all.x = TRUE)
  base::colnames(mtable) <- c("FeatureID", "SampleID", "VALUE") # rename columns

  mtable$VALUE <- .gtToNum(mtable$VALUE)
  mtable$VALUE[is.na(mtable$VALUE)] <- fill # fill in misssing values
  data.table::setkeyv(info, c("FeatureID", "SampleID")) # set key

  ## Return Mutation Data Table Object
  MDT(mtable = mtable, annotations = annotations, info = info)
})

# importPhenotype ----
#' @rdname importPhenotype
setMethod("importPhenotype", "MDT", function(x,
                                             file,
                                             response_type, ...) {

  ## Check if does not already have phenotype data
  if (length(phenotype(x)) != 0) {
    warning("MDT object already had phenotype data. Overwriting.")
  }

  ## Import phenotype data
  phenotype <- data.table::fread(file, ...)

  if (nrow(phenotype) == 0) {
    stop("Imported file had no rows.")
  }

  base::colnames(phenotype)[base::colnames(phenotype) != "SampleID"] <-
    base::toupper(base::colnames(phenotype)[base::colnames(phenotype) != "SampleID"])

  if (! any(c("SampleID", "RESPONSE") %in% base::colnames(phenotype))) {
    stop("File must contain at least two column names: 'SampleID' and 'RESPONSE'")
  }

  phenotype$SampleID <- as.character(phenotype$SampleID)
  data.table::setkeyv(phenotype, "SampleID")

  ## Match SampleIDs
  phenotype <- phenotype[samples(x)]

  ## Recode response
  if (response_type == "factor") {
    phenotype$RESPONSE <- as.factor(phenotype$RESPONSE)

  } else if (response_type == "numeric") {
    phenotype$RESPONSE <- as.numeric(phenotype$RESPONSE)

  } else if (response_type == "ingeter") {
    phenotype$RESPONSE <- as.integer(phenotype$RESPONSE)

  } else if (response_type == "character") {
    phenotype$RESPONSE <- as.character(phenotype$RESPONSE)

  } else {
    stop("Invalid 'response_type' argument.")
  }

  ## Return MDT object
  MDT(mtable = x@mtable,
      annotations = x@annotations,
      phenotype = phenotype,
      info = x@info)
})

## aggregateMDT ----
#' @rdname aggregateMDT
setMethod("aggregateMDT", "MDT", function(x,
                                          group,
                                          fun.aggregate, ...) {


  ## Format group argument.
  ## Must be a data.frame and have two FeatureID columns: 'OLD' and 'NEW'
  group <- data.table::as.data.table(group)
  data.table::setkeyv(group, NULL)
  group <- base::unique(stats::na.omit(group))

  if (base::ncol(group) != 2) {
    stop(.pn("Invalid 'group' argument: must be a data.frame with two columns.",
             "The first is the OLD FeatureID, the second is the NEW FeatureID."))
  }

  base::colnames(group) <- c("OLD", "NEW")
  group$OLD <- as.character(group$OLD)
  group$NEW <- as.character(group$NEW)
  data.table::setkeyv(group, "OLD")

  if (! all(group$OLD %in% MachineLearningGWAS::features(x))) {
    warning("Invalid 'group' argument: values in the first column of 'group' do not match features(x).")
  }

  ## Merge group with mtable
  mtable <- merge(group, MachineLearningGWAS::mtable(x),
                  by.x = "OLD", by.y = "FeatureID",
                  all.x = TRUE, allow.cartesian = TRUE)
  mtable$OLD <- NULL
  base::colnames(mtable) <- c("FeatureID", "SampleID", "VALUE")
  data.table::setkeyv(mtable, NULL)
  mtable <- base::unique(mtable)
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))
  mtable <- base::with(mtable, mtable[, fun.aggregate(VALUE), by = .(FeatureID, SampleID)])
  base::colnames(mtable) <- c("FeatureID", "SampleID", "VALUE")
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))

  ## Get rid of annotations
  featid <- base::unique(as.character(mtable$FeatureID))
  sampid <- base::unique(as.character(mtable$SampleID))
  annotations <- data.table::data.table(FeatureID = featid)
  info <- data.table::data.table(base::expand.grid(FeatureID = featid, SampleID = sampid))
  data.table::setkeyv(annotations, "FeatureID")
  data.table::setkeyv(info, c("FeatureID", "SampleID"))

  MDT(mtable = mtable,
      annotations = annotations,
      phenotype = x@phenotype,
      info = info)
})

## mdtToGRanges ----
#' @rdname mdtToGRanges
setMethod("mdtToGRanges", "MDT", function(x) {

  if (any(! c("CHROM", "POS", "STRAND", "REF", "ALT") %in% colnames(annotations(x)))) {
    stop("Missing CHROM, POS, STRAND, REF, and/or ALT columns in annotations(x).")
  }

  df <- MachineLearningGWAS::annotations(x)[, c("FeatureID", "CHROM", "POS", "STRAND", "REF", "ALT"), with = FALSE]
  df <- base::unique(df)

  alt <- base::lapply(seq_len(nrow(df)), function(i) {
    data.table::data.table(FeatureID = df$FeatureID[i],
                           ALT = strsplit(df$ALT[i], ",", TRUE)[[1]])
  })
  alt <- data.table::rbindlist(alt)
  df[, "ALT":=NULL]
  data.table::setkey(alt, "FeatureID")
  df <- merge(df, alt, all = TRUE, by = "FeatureID")
  base::rm(alt)

  df <- df[base::order(df$CHROM, df$POS, df$STRAND, df$REF, df$ALT), ]

  gloc <- GenomicRanges::GRanges(S4Vectors::Rle(df$CHROM),
                                 IRanges::IRanges(df$POS, df$POS),
                                 S4Vectors::Rle(df$STRAND))
  gloc$FeatureID <- df$FeatureID
  gloc$REF <- df$REF
  gloc$ALT <- df$ALT

  gloc
})

## selectMDT ----
#' @rdname selectMDT
setMethod("selectMDT", c(x = "MDT"), function(x, keys, columns, keytype, na.rm) {

  ## Check keytype
  keytype <- keytype[1]
  if (! keytype %in% c('annotations', 'phenotype', 'info')) {
    stop("Invalid 'keytype' argument: select between 'annotations', 'phenotype' or 'info'.")
  }

  ## Convert keys to data.table if necessary
  if (! data.table::is.data.table(keys)) {

    if (keytype == "info") {
      stop("Invalid combination of 'keys' and 'keytype': 'keys' must be data.table when keytype 'info'")
    }

    keys <- base::tryCatch(data.table::as.data.table(as.character(keys)),
                           error = function(e) {
                             stop("Invalid 'keys' argument: must be convertible to data.table.")
                           })
    if (keytype == "annotations") {
      base::colnames(keys) <- "FeatureID"

    } else if (keytype == "phenotype") {
      base::colnames(keys) <- "SampleID"
    }
  }

  ## Select annotation data.table
  if (keytype == "annotations") {

    db <- MachineLearningGWAS::annotations(x)

    ## Check if keys are FeatureID
    if (! "FeatureID" %in% base::colnames(keys)) {
      stop("Invalid 'keys' argument: must contain 'FeatureID' column.")
    }

    ## Remove extra columns
    if (any(base::colnames(keys) != "FeatureID")) {
      keys[, base::colnames(keys)[base::colnames(keys) != "FeatureID"]:=NULL]
    }

  } else if (keytype == "phenotype") {

    db <- MachineLearningGWAS::phenotype(x)

    if (! "SampleID" %in% base::colnames(keys)) {
      stop("Invalid 'keys' argument: must contain 'SampleID' column.")
    }

    ## Remove extra columns
    if (any(base::colnames(keys) != "SampleID")) {
      keys[, base::colnames(keys)[base::colnames(keys) != "SampleID"]:=NULL]
    }

  } else if (keytype == "info") {

    db <- MachineLearningGWAS::info(x)

    if (all(! c("FeatureID", "SampleID") %in% base::colnames(keys))) {
      stop("Invalid 'keys' argument: must contain 'FeatureID' and/or 'SampleID' column.")
    }

    ## Remove extra columns
    if (any(! base::colnames(keys) %in% c("FeatureID", "SampleID"))) {
      keys[, base::colnames(keys)[! base::colnames(keys) %in% c("FeatureID", "SampleID")]:=NULL]
    }
  }

  ## Convert keys to character and set key
  keys <- keys[, base::lapply(.SD, as.character)]
  data.table::setkeyv(keys, base::colnames(keys))

  ## Subset annotation db according to columns
  if (any(! columns %in% base::colnames(db))) {
    stop("Invalid 'columns' argument: must correspond to colname in 'keytype'.")
  }

  ## Combine keys and annotation db
  results <- merge(keys, db[, c(base::colnames(keys), columns), with = FALSE],
                   all.x = TRUE, by = base::colnames(keys))

  if (na.rm) {
    results <- stats::na.omit(results)
  }

  ## Remove non unique rows
  data.table::setkeyv(results, NULL)
  results <- base::unique(results)
  data.table::setkeyv(results, colnames(keys))

  results

})

# fillMDT ----
#' @rdname fillMDT
setMethod("fillMDT", "MDT", function(x, fill) {

  fullmtable <- data.table::CJ(FeatureID = features(x),
                               SampleID = samples(x))
  data.table::setkeyv(fullmtable, c("FeatureID", "SampleID"))
  results <- merge(fullmtable, MachineLearningGWAS::mtable(x),
                   all.x = TRUE, by = c("FeatureID", "SampleID"))
  results$VALUE[is.na(results$VALUE)] <- fill
  MDT(mtable = results,
      annotations = x@annotations,
      phenotype = x@phenotype,
      info = x@info)
})

##
## 3. Plots  -------------------------------------------------------------------
##

## plotMDT ----
#' @rdname plotMDT
setMethod("plotMDT", "MDT", function(x,
                                     type,
                                     preprocess,
                                     space,
                                     dims,
                                     fill,
                                     ...) {


  ## Check arguments
  eout <- "Invalid 'type' argument: must be 'scatter' or 'heatmap'."
  type <- tryCatch(type[1], error = function(e) {
    stop(eout)
  })
  if (! type %in% c("scatter", "heatmap")) {
    stop(eout)
  }

  if (! is.null(preprocess)) {
    if (! any(preprocess %in% c("scale", "center", "pca"))) {
      stop("Invalid 'preprocess' argument: at least one entry does not correspond to 'scale', 'center' or 'pca'.")
    }
  }

  eout <- "Invalid 'space' argument: must be 'features' or 'samples'"
  space <- tryCatch(space[1], error = function(e) {
    stop(eout)
  })
  if (! space %in% c("features", "samples")) {
    stop(eout)
  }

  dims <- tryCatch(dims[1:2], error = function(e) {
    stop("Invalid features argument: must be two numerics or two characters")
  })

  ## Retrieve data
  mat <- MachineLearningGWAS::asMatrixMDT(x, fill)

  ## Check if numeric
  if (! is.numeric(mat)) {
    mat <- tryCatch(apply(mat, 2, as.numeric), error = function(e) {
      stop("Could not convert MDT object to numeric matrix. Check output of asMatrixMDT(x).")
    })
  }

  ## Select space
  if (space == "features") {
    mat <- base::t(mat)
  }

  ## Preprocess data
  if (! is.null(preprocess)) {

    if (all(c("scale", "center") %in% preprocess)) {
      mat <- base::scale(mat, center = TRUE, scale = TRUE)

    } else if ("scale" %in% preprocess) {
      mat <- base::scale(mat, center = FALSE, scale = TRUE)

    } else if ("center" %in% preprocess) {
      mat <- base::scale(mat, center = TRUE, scale = FALSE)
    }

    if ("pca" %in% preprocess) {
      pca <- stats::prcomp(mat, center = FALSE, scale. = FALSE)
      mat <- pca$x
    }

  }

  if (type == "heatmap") {

    if (length(response(x)) > 0 &&
        (is.factor(response(x)) | is.character(response(x))) &&
        space != "features") {

      cols <- grDevices::rainbow(length(unique(response(x))))
      gplots::heatmap.2(mat,
                        trace = "none",
                        RowSideColors = cols[response(x)], ...)

    } else {
      gplots::heatmap.2(mat,
                        col = "rainbow",
                        trace = "none",
                        ...)
    }

  } else if (type == "scatter") {

    df <- base::as.data.frame(mat)

    ## Name of corresponding label
    if (all(is.numeric(dims))) { lab <- base::colnames(df)[dims]
    } else { lab <- as.character(dims) }

    base::colnames(df) <- base::make.names(base::colnames(df), unique = TRUE, allow_ = FALSE)

    ## Name of column to select
    if (all(is.numeric(dims))) {
      sel <- base::colnames(df)[dims]
    } else {
      sel <- as.character(dims)
    }

    ## Create base plot
    p <- ggplot2::ggplot(df, ggplot2::aes_string(sel[1], sel[2]))

    if (length(MachineLearningGWAS::response(x)) == 0 | space == "features") {
      p <- p + ggplot2::geom_point(...)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes_string(color = base::eval("response(x)")), ...)
    }

    p <- p + ggplot2::xlab(lab[1]) + ggplot2::ylab(lab[2])

    p

  }
})


## manhattanPlot ----
#' @rdname manhattanPlot
setMethod("manhattanPlot", c("MDT"), function(x,
                                              value,
                                              group,
                                              pos,
                                              f,
                                              col) {

  group <- group[1]
  pos <- pos[1]
  value <- value[1]

  if (any(! c(group, pos, value) %in% base::colnames(annotations(x)))) {
    stop("Invalid 'group', 'pos' or 'value' argument: does not match any column in annotations(x)")
  }

  df <- MachineLearningGWAS::annotations(x)[, c("FeatureID", group, pos, value), with = FALSE]
  colnames(df) <- c("FeatureID", "Group", "Position", "Value")

  if (! is.numeric(df$Position)) {
    stop("Invalid 'pos' argument: must correspond to a numeric column in annotations(x)")
  }

  if (! is.numeric(df$Value)) {
    stop("Invalid 'column' argument: must correspond to a numeric column in annotations(x)")
  }

  df$Value <- f(df$Value)

  .manhattanPlotSupport(df, col)
})


##
## 4. Subset and Merge  --------------------------------------------------------
##

# filterMDT ----
#' @rdname filterMDT
setMethod("filterMDT", c("MDT"), function(x, condition, with) {

  if (missing(condition)) condition <- TRUE

  ## Select correct annotation data.table
  with <- with[1]
  if ("annotations" %in% with) {
    db <- MachineLearningGWAS::annotations(x)
    data.table::setkeyv(db, "FeatureID")
  } else if ("phenotype" %in% with) {
    db <- MachineLearningGWAS::phenotype(x)
    data.table::setkeyv(db, "SampleID")
  } else if ("info" %in% with) {
    db <- MachineLearningGWAS::info(x)
    data.table::setkeyv(db, c("FeatureID", "SampleID"))
  } else {
    stop("Invalid 'with' argument: choose between 'annotations', 'phenotype', 'info'.")
  }

  if (length(db) == 0) {
    stop("Invalid 'with' argument: empty data.table.")
  }

  ## Evaluate condition
  e <- substitute(condition)
  r <- base::eval(e, db, base::parent.frame())

  if (! is.logical(r)) {
    stop("Invalid 'condition' argument: must be logical.")
  }

  r <- r & ! is.na(r)

  ## Subset MDT
  subdb <- base::unique(db[r, key(db), with = FALSE])
  data.table::setkeyv(subdb, key(db))
  mtable <- merge(subdb, mtable(x),
                  all.x = TRUE, all.y = FALSE,
                  allow.cartesian = TRUE)

  ## Check keys that had no annotation entries and keep them
  db2 <- base::unique(db[, data.table::key(db), with = FALSE])
  db2$ISTHERE <- TRUE
  db2 <- merge(db2, mtable(x), all = TRUE)
  db2 <- db2[is.na(db2$ISTHERE), ]
  db2[, "ISTHERE":=NULL]
  mtable <- base::rbind.data.frame(mtable, db2)
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))
  base::rm(db, db2, r, e)

  if (base::nrow(mtable) == 0) {
    stop("Filtering 'condition' resulted in empty MDT.")
  }

  ## Remove irrelevant annotations
  if (length(annotations(x)) > 0 & ! is.null(subdb$FeatureID)) {
    annotations <- MachineLearningGWAS::annotations(x)[base::unique(stats::na.omit(as.character(subdb$FeatureID)))]
  } else {
    annotations <- MachineLearningGWAS::annotations(x)
  }
  if (length(MachineLearningGWAS::phenotype(x)) > 0 & ! is.null(subdb$SampleID)) {
    phenotype <- MachineLearningGWAS::phenotype(x)[base::unique(stats::na.omit(as.character(subdb$SampleID)))]
  } else {
    phenotype <- MachineLearningGWAS::phenotype(x)
  }
  if (length(info(x)) > 0) {
    info <- merge(subdb, MachineLearningGWAS::info(x),
                  all.x = TRUE, all.y = FALSE,
                  by = data.table::key(subdb),
                  allow.cartesian = TRUE)
  } else {
    info <- MachineLearningGWAS::info(x)
  }

  MDT(mtable = mtable,
      annotations = annotations,
      phenotype = phenotype,
      info = info)
})

# merge ----
#' @rdname annotateMDT
setMethod("merge", c(x = "MDT", y = "data.frame"), function(x,
                                                            y,
                                                            by = "FeatureID",
                                                            all = FALSE,
                                                            all.x = TRUE,
                                                            all.y = FALSE,
                                                            sort = TRUE,
                                                            allow.cartesian = TRUE) {
  annotations <-  .mergeDT(x = MachineLearningGWAS::annotations(x), y = y,
                           by = by,
                           all.x = all.x, all.y = all.y,
                           sort = sort,
                           allow.cartesian = allow.cartesian)
  MDT(mtable = x@mtable,
      annotations = annotations,
      phenotype = x@phenotype,
      info = x@info)
})

setMethod("[", c(x = "MDT", i = "character", j = "missing"), function(x,
                                                                      i,
                                                                      j,
                                                                      ...,
                                                                      drop = TRUE) {

  ## Retrieve data
  mtable <- MachineLearningGWAS::mtable(x)
  phenotype <- MachineLearningGWAS::phenotype(x)
  info <- MachineLearningGWAS::info(x)

  ## Check if i arg are actual SampleIDs
  i <- unique(stats::na.omit(i))
  if (any(! i %in% MachineLearningGWAS::samples(x))) {
    stop("Invalid 'i' argument: does not correspond to any SampleID")
  }

  # Subset
  mtable <- mtable[mtable$SampleID %in% i, ]

  if (length(phenotype) != 0) {
    phenotype <- phenotype[phenotype$SampleID %in% i, ]
  } else {
    phenotype <- data.table::data.table()
  }

  if (length(info) != 0) {
    info <- info[info$SampleID %in% i, ]
  } else {
    info <- data.table::data.table()
  }

  ## Return new object
  MDT(mtable = mtable,
      annotations = x@annotations,
      phenotype = phenotype,
      info = info)
})

setMethod("[", c(x = "MDT", i = "missing", j = "character"), function(x,
                                                                      i,
                                                                      j,
                                                                      ...,
                                                                      drop = TRUE) {

  ## Retrieve data
  mtable <- MachineLearningGWAS::mtable(x)
  annotations <- MachineLearningGWAS::annotations(x)
  info <- MachineLearningGWAS::info(x)

  ## Check if j arg are actual FeatureIDs
  j <- base::unique(stats::na.omit(j))
  if (any(! j %in% MachineLearningGWAS::features(x))) {
    stop("Invalid 'j' argument: does not correspond to any FeatureID")
  }

  ## Subset for FeatureID
  mtable <- mtable[mtable$FeatureID %in% j, ]

  if (length(annotations) != 0) {
    annotations <- annotations[annotations$FeatureID %in% j, ]
  } else {
    annotations <- data.table::data.table()
  }

  if (length(info) != 0) {
    info <- info[info$FeatureID %in% j, ]
  } else {
    info <- data.table::data.table()
  }

  ## Return new object
  MDT(mtable = mtable,
      annotations = annotations,
      phenotype = x@phenotype,
      info = info)
})

setMethod("[", c(x = "MDT", i = "character", j = "character"), function(x,
                                                                        i,
                                                                        j,
                                                                        ...,
                                                                        drop = TRUE) {
  x <- x[, j]
  x <- x[i, ]
  x
})

##
## 5. General Methods ----------------------------------------------------------
##

## as.matrix ----
#' @rdname asMatrixMDT
setMethod("asMatrixMDT", "MDT", function(x, fill, fun.aggregate, ...) {

  if (missing(fun.aggregate)) {
    mat <- data.table::dcast.data.table(MachineLearningGWAS::mtable(x),
                                        formula = SampleID ~ FeatureID,
                                        value.var = "VALUE",
                                        fill = fill, ...)
  } else {
    mat <- data.table::dcast.data.table(MachineLearningGWAS::mtable(x),
                                        formula = SampleID ~ FeatureID,
                                        fun.aggregate = fun.aggregate,
                                        value.var = "VALUE",
                                        fill = fill, ...)
  }
  rownames <- mat$SampleID
  mat[, "SampleID":=NULL]
  mat <- base::as.matrix(mat)
  base::rownames(mat) <- rownames
  mat
})

## show ----
setMethod("show", "MDT", function(object) {

  mtable <- MachineLearningGWAS::mtable(object)
  annotations <- MachineLearningGWAS::annotations(object)
  info <- MachineLearningGWAS::info(object)
  response <- MachineLearningGWAS::response(object)
  phenotype <- MachineLearningGWAS::phenotype(object)

  cat("--- MDT object ---\n")

  cat(paste("mtable(x):     ",
            methods::selectMethod("nrow", "MDT")(object),
            "individuals and",
            methods::selectMethod("ncol", "MDT")(object),
            "features\n"))

  if (length(annotations) != 0) {

    base::cat(base::paste("annotations(x):", base::nrow(annotations),
                          "rows in", base::ncol(annotations)-1, "columns\n"))

  }

  if (length(info) != 0) {

    base::cat(base::paste("info(x):       ", base::nrow(info),
                          "rows in", base::ncol(info)-2, "columns\n"))

  }

  if (length(phenotype) != 0) {

    if (is.character(response) | base::is.factor(response)) {
      t <- base::table(response)
      n <- base::paste0(names(t), sep = ": ")
      v <- base::apply(base::rbind(n, t), 2, paste, collapse = "")
      v <- base::paste0(v, collapse = " | ")
    } else {
      v <- base::paste("with a range from", base::paste(range(response), collapse = " to "))
    }

    base::cat("response(x):   ", class(response), v, "\n")

    base::cat(base::paste("phenotype(x):  ", base::nrow(phenotype),
                          "rows in", base::ncol(phenotype), "columns\n"))
  }
})


## nrow ----
setMethod("nrow", "MDT", function(x) {
  length(MachineLearningGWAS::samples(x))
})

## ncol ----
setMethod("ncol", "MDT", function(x) {
  length(MachineLearningGWAS::features(x))
})

## length ----
setMethod("length", "MDT", function(x) {
  methods::selectMethod("nrow", "MDT")(x)
})

## dim ----
setMethod("dim", "MDT", function(x) {
  c(methods::selectMethod("nrow", "MDT")(x),
    methods::selectMethod("ncol", "MDT")(x))
})

##
## 6. Setters and getters ------------------------------------------------------
##

# mtable ----
setMethod("mtable", "MDT", function(x) {
  x@mtable
})

setReplaceMethod("mtable", "MDT", function(x, value) {
  value <- data.table::as.data.table(value)
  data.table::setkeyv(value, c("FeatureID", "SampleID"))
  MDT(mtable = value,
      annotations = x@annotations,
      phenotype = x@phenotype,
      info = x@info)
})

# annotations ----
setMethod("annotations", "MDT", function(x) {
  x@annotations
})

setReplaceMethod("annotations", "MDT", function(x, value) {
  value <- data.table::as.data.table(value)
  data.table::setkeyv(value, "FeatureID")
  MDT(mtable = x@mtable,
      annotations = value,
      phenotype = x@phenotype,
      info = x@info)
})

# phenotype ----
setMethod("phenotype", "MDT", function(x) {
  x@phenotype
})

setReplaceMethod("phenotype", "MDT", function(x, value) {

  if (length(x@phenotype) == 0) {
    stop("No phenotype data: please use 'importPhenotype'\n")
  }

  value <- data.table::as.data.table(value)
  MDT(mtable = x@mtable,
      annotations = x@annotations,
      phenotype = value,
      info = x@info)
})

# info ----
setMethod("info", "MDT", function(x) {
  x@info
})

setReplaceMethod("info", "MDT", function(x, value) {
  value <- data.table::as.data.table(value)
  data.table::setkeyv(value, c("FeatureID", "SampleID"))
  MDT(mtable = x@mtable,
      annotations = x@annotations,
      phenotype = x@phenotype,
      info = value)
})

# response ----
setMethod("response", "MDT", function(x) {

  if (length(x@phenotype) == 0) {
    NULL
  } else {
    r <- x@phenotype$RESPONSE
    names(r) <- x@phenotype$SampleID
    r
  }
})

setReplaceMethod("response", "MDT", function(x, value) {

  phenotype <- phenotype(x)

  if (length(phenotype) == 0) {
    stop("No phenotype data: please use 'importPhenotype'\n")
  }
  phenotype$RESPONSE <- value
  MDT(mtable = x@mtable,
      annotations = x@annotations,
      phenotype = phenotype,
      info = x@info)
})

# features ----
setMethod("features", "MDT", function(x) {
  x@features
})

setReplaceMethod("features", "MDT", function(x, value) {

  if (length(value) != length(features(x))) {
    stop("Length of values is not equal to lenght of features. Use `[` to remove features.")
  }

  value <- as.character(value)

  ## Retrieve data
  mtable <- MachineLearningGWAS::mtable(x)
  annotations <- MachineLearningGWAS::annotations(x)
  info <- MachineLearningGWAS::info(x)

  ## Change values
  mtable$FeatureID <- .changeLevels(mtable$FeatureID, value)
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))

  if (length(annotations) != 0L) {
    annotations$FeatureID <- .changeLevels(annotations$FeatureID, value)
    data.table::setkeyv(annotations, "FeatureID")
  } else {
    annotations <- data.table::data.table()
  }

  if (length(info) != 0L) {
    info$FeatureID <- .changeLevels(info$FeatureID, value)
    data.table::setkeyv(info, c("FeatureID", "SampleID"))
  } else {
    info <- data.table::data.table()
  }

  ## Return new object
  MDT(mtable = mtable,
      annotations = annotations,
      phenotype = x@phenotype,
      info = info)
})

# samples ----
setMethod("samples", "MDT", function(x) {
  x@samples
})

setReplaceMethod("samples", "MDT", function(x, value) {

  if (length(value) != length(samples(x))) {
    stop("Length of values is not equal to lenght of samples Use `[` to remove samples")
  }

  value <- as.character(value)

  ## Retrieve data (no annotations)
  mtable <- MachineLearningGWAS::mtable(x)
  phenotype <- MachineLearningGWAS::phenotype(x)
  info <- MachineLearningGWAS::info(x)

  ## Change values
  mtable$SampleID <- .changeLevels(mtable$SampleID, value)
  data.table::setkeyv(mtable, c("FeatureID", "SampleID"))

  if (length(phenotype) != 0L) {
    phenotype$SampleID <- .changeLevels(phenotype$SampleID, value)
    data.table::setkeyv(phenotype, "SampleID")
  } else {
    phenotype <- data.table::data.table()
  }

  if (length(info) != 0L) {
    info$SampleID <- .changeLevels(info$SampleID, value)
    data.table::setkeyv(info, c("FeatureID", "SampleID"))
  } else {
    info <- data.table::data.table()
  }

  ## Return new object
  MDT(mtable = mtable,
      annotations = x@annotations,
      phenotype = phenotype,
      info = info)
})
