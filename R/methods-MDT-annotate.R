################################################################################
## annotate* MDT Methods
################################################################################


# annotateMAF ----
#' @rdname annotateMDT
setMethod("annotateMAF", "MDT", function(x) {

  ## Retrieve data
  geno <- base::unique(MachineLearningGWAS::mtable(x))

  ## Check if mtable(x)$VALUE is integer
  geno$VALUE <- base::tryCatch(geno$VALUE,
                               error = function(e) {
                                 warning(e)
                                 stop("mtable(x)$VALUE must be convertible to integer")
                               }, warning = function(w) {
                                 warning(w)
                                 stop("mtable(x)$VALUE must be convertible to integer")
                               })


  ## Assign non 0, 1, 2 to NA
  na <- geno$VALUE > 2L
  geno[na, ] <- NA
  if (any(na)) {
    warning("Some values in mtable(x)$VALUE were bigger than 2. These have been set to NA.")
  }

  ## Convert to matrix
  geno <- data.table::dcast.data.table(geno, SampleID ~ FeatureID, value.var = "VALUE")
  sampid <- geno$SampleID
  geno[, "SampleID":=NULL]
  geno <- base::as.matrix(geno)
  base::rownames(geno) <- sampid
  base::rm(sampid)

  ## Count genotypes per SNP
  n0 <- base::apply(geno == 0L, 2, sum, na.rm=T)
  n1 <- base::apply(geno == 1L, 2, sum, na.rm=T)
  n2 <- base::apply(geno == 2L, 2, sum, na.rm=T)
  n <- n0 + n1 + n2

  ## Calculate allele frequencies
  p <- ((2 * n0) + n1) / (2 * n)
  q <- 1 - p
  maf <- base::pmin(p, q)
  # mgf <- apply(cbind(n0, n1, n2), 1, min) / n

  ## Merge results
  annotations(x) <- .mergeDT(MachineLearningGWAS::annotations(x),
                             data.table::data.table(FeatureID = names(maf), MAF = as.numeric(maf)),
                             by = "FeatureID")
  x
})

## annotateStat ----
#' @rdname annotateMDT
setMethod("annotateStat", "MDT", function(x,
                                          function_list,
                                          fill,
                                          verbose) {


  ## Choose lapply function to use
  mylapply <- .chooseLapply(verbose)

  ## Get values
  mat <- MachineLearningGWAS::asMatrixMDT(x, fill = fill)
  response <- MachineLearningGWAS::response(x)

  ## Check arguments
  if (length(response) == 0) {
    stop(.pn("Invalid 'x' argument: no phenotype data.",
             "Please call importPhenotype."))
  }

  # Names that values will take in annotations slot
  cname <- base::make.names(names(function_list))

  ## DataFrame to contain statistical values
  results <- data.table::data.table(matrix(0,
                                           nrow = base::ncol(mat),
                                           ncol = length(function_list) + 1))
  base::colnames(results) <- c("FeatureID", cname)
  results$FeatureID <- MachineLearningGWAS::features(x)

  ## Compute statistical values
  ## If computation fails for one table, return NA_real_
  for (i in seq_along(function_list)) {
    if (verbose) base::cat(base::paste("Estimating", cname[i], "\n"))

    results[, cname[i]] <- base::unlist(mylapply(seq_len(base::ncol(mat)), function(j) {
      x <- mat[, j]
      base::tryCatch(function_list[[i]](x, response), error = function(e) {
        warning(e)
        return(NA_real_)
      })
    }))
    if (verbose) cat("\n")
  }

  ## Format results and merge to annotation slot of MDT
  if (verbose) base::cat("Merging results with annotations(x)\n")
  results <- .mergeDT(MachineLearningGWAS::annotations(x), results, by = "FeatureID")
  data.table::setkey(results, "FeatureID")
  MachineLearningGWAS::annotations(x) <- results
  x
})

## annotateLocation ----
#' @rdname annotateMDT
setMethod("annotateLocation", "MDT", function(x,
                                              txdb,
                                              region,
                                              columns,
                                              verbose) {

  ## Check arguments
  if (all(! class(txdb) %in% c("TxDb", "GRangesList"))) {
    stop(.pn("Invalid 'txdb' argument: must be of class 'TxDb' or 'GRangesList'.",
             "See 'subject' argument of 'locateVariants' function in 'VariantAnnotation' package."))
  }

  ## Retrieve data
  gloc <- MachineLearningGWAS::mdtToGRanges(x)
  GenomeInfoDb::genome(gloc) <- GenomeInfoDb::genome(txdb)
  GenomeInfoDb::seqlevelsStyle(gloc) <- GenomeInfoDb::seqlevelsStyle(txdb)

  ## Locating variants in and around genes
  if (verbose) cat("Locating variants\n")
  results <- VariantAnnotation::locateVariants(query = gloc,
                                               subject = txdb,
                                               region = region)

  ## Convert to data.table
  results <- data.table::as.data.table(stats::na.omit(SummarizedExperiment::mcols(results)))

  ## Change name of GENEID to ENTREZID
  geneid <- base::colnames(results) == "GENEID"
  base::colnames(results)[geneid] <- "ENTREZID"

  ## Remove irrelevant columns from results
  cs <- base::unique(c(base::colnames(results)[base::colnames(results) %in% columns],
                       "QUERYID"))
  results <- results[, cs, with = FALSE]

  ## Merge with annotations(x)
  if (verbose) base::cat("Merging results with annotations(x)\n")
  results$FeatureID <- gloc$FeatureID[results$QUERYID]
  results$QUERYID <- NULL
  data.table::setkey(results, "FeatureID")
  results <- .mergeDT(MachineLearningGWAS::annotations(x), results, by = "FeatureID")
  data.table::setkey(results, "FeatureID")
  MachineLearningGWAS::annotations(x) <- results

  x
})

## annotateBioconductor ----
#' @rdname annotateMDT
setMethod("annotateBioconductor", "MDT", function(x,
                                                  annotations_db,
                                                  keys_colname,
                                                  columns,
                                                  keytype,
                                                  verbose) {

  ## Check arguments
  keys_colname <- keys_colname[1]
  if (! keys_colname %in% base::colnames(MachineLearningGWAS::annotations(x))) {
    stop("Invalid 'keys_colname' argument: column not in annotations(x)")
  }

  ## Annotate variants
  keys <- stats::na.omit(base::unique(MachineLearningGWAS::annotations(x)[[keys_colname]]))

  if (verbose) cat("Annotating variants\n")
  results <- AnnotationDbi::select(annotations_db,
                                   keys = keys,
                                   columns = columns,
                                   keytype = keytype)

  if (verbose) cat("Merging results with annotations(x)\n")
  results <- .mergeDT(MachineLearningGWAS::annotations(x), results, by = keytype)
  data.table::setkey(results, "FeatureID")
  MachineLearningGWAS::annotations(x) <- results
  x
})

# annotateCoding ----
#' @rdname annotateMDT
setMethod("annotateCoding", "MDT", function(x,
                                            seqSource,
                                            txdb,
                                            columns,
                                            verbose) {


  ## Check arguments
  if (all(! class(txdb) %in% c("TxDb"))) {
    stop(.pn("Invalid 'txdb' argument: must be of class 'TxDb'.",
             "See 'subject' argument of 'predictCoding' function in 'VariantAnnotation' package."))
  }

  if (all(! class(seqSource) %in% c("BSgenome", "FaFile"))) {
    stop("Invalid 'seqSource' argument: must be of class 'BSgenome' or 'FaFile'.",
         "See 'seqSource' argument of 'predictCoding' function in 'VariantAnnotation' package.")
  }

  if (! "ALT" %in% colnames(annotations(x))) {
    stop("Missing 'ALT' (alternative allele) column in annotations(x).")
  }

  ## Retrieve data
  gloc <- MachineLearningGWAS::mdtToGRanges(x)
  GenomeInfoDb::genome(gloc) <- GenomeInfoDb::genome(txdb)
  GenomeInfoDb::seqlevelsStyle(gloc) <- GenomeInfoDb::seqlevelsStyle(txdb)
  GenomeInfoDb::seqlevelsStyle(seqSource) <- GenomeInfoDb::seqlevelsStyle(txdb)

  ## Variant Annotations
  if (verbose) cat("Predicting coding\n")
  coding <- VariantAnnotation::predictCoding(query = gloc,
                                             subject = txdb,
                                             seqSource = seqSource,
                                             varAllele = Biostrings::DNAStringSet(gloc$ALT))

  ## Keep only columns as specified by 'columns'
  columns <- c(columns, "FeatureID", "ALT")
  results <- data.table::as.data.table(SummarizedExperiment::mcols(coding))
  results[, base::colnames(results)[! base::colnames(results) %in% columns]:=NULL]
  results <- base::unique(results)
  base::colnames(results)[base::colnames(results) == "ALT"] <- "CONSEQUENCE_ALT"

  ## Merge with annotations
  if (verbose) cat("Merging results with annotations(x)\n")
  results <- .mergeDT(MachineLearningGWAS::annotations(x), results,
                      by = c("FeatureID"), all.x = TRUE)
  data.table::setkeyv(results, "FeatureID")
  MachineLearningGWAS::annotations(x) <- results

  x
})

# annotateRSID ----
#' @rdname annotateMDT
setMethod("annotateRSID", "MDT", function(x, snp_loc, verbose) {

  if (verbose) cat("Retrieving SNPs rs Identification Numbers\n")

  ## Retrieve genomic locations, remove SampleID and make unique
  gloc <- MachineLearningGWAS::mdtToGRanges(x)
  GenomeInfoDb::genome(gloc) <- GenomeInfoDb::genome(snp_loc)
  GenomeInfoDb::seqlevelsStyle(gloc) <- GenomeInfoDb::seqlevelsStyle(snp_loc)

  ## Find SNPs rs IDs by overlap between genomic regions
  snp <- BSgenome::snpsByOverlaps(snp_loc, gloc, type = "any")

  ## Identify overlapping genomic alignments
  hits <- GenomicRanges::findOverlaps(gloc, snp)
  subject_hits <- S4Vectors::subjectHits(hits)
  query_hits <- S4Vectors::queryHits(hits)

  ## Add RSID to annotations
  rsids <- SummarizedExperiment::mcols(snp)[subject_hits, "RefSNP_id"]
  rsids <- data.table::data.table(FeatureID = gloc$FeatureID[query_hits],
                                  RSID = rsids)

  ## Merge results
  if (verbose) cat("Merging results\n")
  results <- .mergeDT(MachineLearningGWAS::annotations(x), rsids, by = "FeatureID")
  data.table::setkeyv(results, "FeatureID")
  MachineLearningGWAS::annotations(x) <- results
  x
})
