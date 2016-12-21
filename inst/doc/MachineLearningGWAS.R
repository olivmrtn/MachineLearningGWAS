## ---- setup = TRUE, include = FALSE--------------------------------------
knitr::opts_chunk$set(comment = NA, fig.width = 7, fig.height = 7)

## ----eval = FALSE--------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("olivmrtn/MachineLearningGWAS")

## ----message = FALSE, warning = FALSE------------------------------------
library(MachineLearningGWAS)

## ------------------------------------------------------------------------
citation("MachineLearningGWAS")

## ------------------------------------------------------------------------
genotype_file <- system.file("extdata", "snps.vcf", 
                             package = "MachineLearningGWAS")
clinical_file <- system.file("extdata", "clinical.txt", 
                             package = "MachineLearningGWAS")

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
vcf <- VariantAnnotation::readVcf(genotype_file)
mdt <- vcfsToMDT(list(vcf))

## ------------------------------------------------------------------------
mdt <- importPhenotype(mdt, clinical_file, response_type = "factor", sep = "\t")

## ------------------------------------------------------------------------
mdt

## ------------------------------------------------------------------------
head(MachineLearningGWAS::samples(mdt))
head(MachineLearningGWAS::features(mdt))

## ------------------------------------------------------------------------
head(mtable(mdt))

## ------------------------------------------------------------------------
head(annotations(mdt))

## ------------------------------------------------------------------------
head(MachineLearningGWAS::info(mdt))

## ------------------------------------------------------------------------
head(phenotype(mdt))

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
## Annotate Fisher's test and HW equilibrium p-values
mdt <- annotateStat(mdt, function_list = list(FISHER = function(x, y) {    
  stats::fisher.test(table(factor(x, c(0L, 1L, 2L)), 
                           factor(y, c(0L, 1L))))$p.value },
  HW = function(x, y) {    
    HardyWeinberg::HWExact(table(factor(x, c(0L, 1L, 2L))), 
                           verbose = FALSE)$pval
  }), 
  fill = NA, verbose = TRUE)

# Annotate MAF
mdt <- annotateMAF(mdt)

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
mdt <- annotateLocation(mdt, TxDb.Hsapiens.UCSC.hg19.knownGene)

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
library(org.Hs.eg.db)
mdt <- annotateBioconductor(mdt, 
                            annotations_db = org.Hs.eg.db, 
                            keys_colname = "ENTREZID", 
                            columns = c("GENENAME", "GO"))

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
mdt <- annotateCoding(mdt, 
                      BSgenome.Hsapiens.UCSC.hg19,
                      TxDb.Hsapiens.UCSC.hg19.knownGene)

mdt <- annotateRSID(mdt, SNPlocs.Hsapiens.dbSNP144.GRCh37)

## ------------------------------------------------------------------------
head(annotations(mdt))

## ------------------------------------------------------------------------
df <- selectMDT(mdt, MachineLearningGWAS::features(mdt), "ENTREZID")
head(df)

## ------------------------------------------------------------------------
ag_mdt <- aggregateMDT(mdt, df, fun.aggregate = mean)
ag_mdt
head(mtable(ag_mdt))

## ------------------------------------------------------------------------
mdt[c("Sample1", "Sample2"), c("1", "2")]

## ------------------------------------------------------------------------
mdt <- filterMDT(mdt, QUAL > 50, "info")
mdt

## ------------------------------------------------------------------------
fillMDT(mdt, fill = 0L)

## ------------------------------------------------------------------------
plotMDT(mdt, 
        type = "scatter", 
        preprocess = c("scale", "center", "pca"), 
        space = c("samples"), 
        dims = c(1L, 2L))

## ----width = 10, height = 5----------------------------------------------
try(plotMDT(filterMDT(mdt, INFO_AS == 1, "info"),
        type = "heatmap", 
        preprocess = c("scale", "center")))

## ------------------------------------------------------------------------
manhattanPlot(mdt, "FISHER")

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
allfeat <- MachineLearningGWAS::features(mdt)
associated <- unique(MachineLearningGWAS::info(mdt)[INFO_AS == 1, FeatureID])
random <- sample(allfeat[! allfeat %in% associated], 15)
selected <- c(associated, random)
mdt <- mdt[, selected]

clas <- trainClassifier(mdt, id = "clas", phen_vars = "all", method = 'rf')

## ------------------------------------------------------------------------
clas

## ------------------------------------------------------------------------
performances(clas)
summaryPerf(clas)

## ------------------------------------------------------------------------
head(featimp(clas))
head(summaryFeatimp(clas))

## ------------------------------------------------------------------------
head(na.omit(errors(clas)))
head(summaryErrors(clas))

## ---- results = 'hide', message = FALSE, warning = FALSE-----------------
perm <- trainClassifier(mdt, id = "perm",
                        method = 'rf', phen_vars = "all", 
                        permute = TRUE)
fits <- bindML(clas, perm)

## ------------------------------------------------------------------------
plotPerf(fits, "Accuracy")

## ------------------------------------------------------------------------
plotFeatimp(fits)

## ----eval = FALSE--------------------------------------------------------
#  ml <- trainClassifier(mdt, id = "boot",
#                        partition = caret::createResample(response(mdt),
#                                                          times = 10),
#                        trControl = trainControl("none"),
#                        tuneGrid = data.frame(.mtry = 5L),
#                        validation = FALSE)

## ----eval = FALSE--------------------------------------------------------
#  fisherSelector <- function(x, y) {
#    pv <- apply(x, 2, function(col) {
#      fisher.test(table(factor(col, 0:2), y))$p.value
#    })
#    pv <- sort(pv)
#    selected <- names(pv)[1:5]
#    function(x) x[, selected]
#  }

## ----eval = FALSE--------------------------------------------------------
#  params <- expand.grid(method = "rf",
#                        permuted = c(TRUE, FALSE),
#                        preprocess = c(fisherSelector))
#  params$id <- c("perm", "real")
#  fits2 <- trainClassifiers(mdt, params)

## ---- echo = FALSE-------------------------------------------------------
sessionInfo()

