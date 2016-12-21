## This script simulates data for a case/control GWAS study.
## It is used for illustration purposes in the MachineLearningGWAS vignette.
## By no means should it be considered as a realistic simulation.
## It does not take into account SNP interactions or linkage dysequlibrium

simulateData <- function() {

  ## Libraries ----
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

  ## Functions ----
  simRandomSNP <- function(n) {
    q <- runif(1, 0.01, 0.5)
    p <- 1 - q
    x <- factor(sample(0L:2L, n, TRUE, c(p**2, 2*p*q, q**2)), c(0L:2L))
    levels(x) <- c("0/0", "0/1", "1/1")
    x
  }

  simAssociatedSNP <- function(response) {

    case <- response == 1L
    cont <- response == 0L

    q1 <- 0
    q0 <- 0
    while (! all(c(q0, q1) >= 0.01 & c(q0, q1) <= 0.50)) {
      q <- runif(1, 0.01, 0.5)
      e <- runif(1, 0.05, 0.20)
      q1 <- q + e
      q0 <- q - e
    }

    p1 <- 1 - q1
    p0 <- 1 - q0

    x <- numeric(length(response))
    x[case] <- sample(0L:2L, sum(case), TRUE, c(p1**2, 2*p1*q1, q1**2))
    x[cont] <- sample(0L:2L, sum(cont), TRUE, c(p0**2, 2*p0*q0, q0**2))
    x <- factor(x, c(0L:2L))
    levels(x) <- c("0/0", "0/1", "1/1")
    x
  }

  ## Parameters ----
  n <- 200 # number of patients/samples
  p_cases <- 0.5 # proportion of cases
  rand <- 995 # number of SNPs/features unrelated to response
  ascd <- 5 # number of SNPs/features related to response
  genome <- BSgenome.Hsapiens.UCSC.hg19
  snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37

  ## Select random SNPs from
  subsnps <- snpsBySeqname(snps, seqnames(snps)[19:22L])

  vcf <- subsnps[sort(sample(1L:length(subsnps), rand+ascd))]
  vcf$nuc <- IUPAC_CODE_MAP[vcf$alleles_as_ambig]
  vcf$REF <- substr(vcf$nuc, 1L, 1L)
  vcf$ALT <- substr(vcf$nuc, 2L, 2L)
  vcf <- as.data.table(vcf)
  vcf <- vcf[, .(CHROM = gsub("[^0-9]", "", as.character(seqnames)), POS = pos, ID = RefSNP_id, REF, ALT)]

  ## Generate cases and controls ----
  response <- factor(sample(0L:1L, n, TRUE, c(1-p_cases, p_cases)), c(0L, 1L))

  ## Simulate clinical variables
  age <- numeric(n)
  age[response == 1L] <- round(rnorm(sum(response == 1L), mean = 70, sd = 5))
  age[response == 0L] <- round(rnorm(sum(response == 0L), mean = 60, sd = 5))
  sex <- sample(0L:1L, n, TRUE)
  clinical <- data.table(SampleID = paste0("Sample", 1L:n), RESPONSE = response,
                         AGE = age, SEX = sex)

  ## Simulate genotypes ----
  ## Plays on MAF frequency to get association
  geno <- data.table(t(cbind(sapply(seq_len(ascd), function(i) simAssociatedSNP(response)),
                             sapply(seq_len(rand), function(i) simRandomSNP(n)))))
  colnames(geno) <- paste0("Sample", 1L:n)
  shuffle <- sample(1L:nrow(geno))
  geno <- geno[shuffle]
  associated <- integer(ascd+rand)
  associated[match(1:ascd, shuffle)] <- 1L

  ## Add different fields to VCF ----
  vcf$QUAL <- rexp(n, 0.001)
  vcf$FILTER <- "."
  vcf$INFO <- paste0(paste0("NS=", n), ";",
                     paste0("DP=", round(rexp(n, 0.01))), ";",
                     paste0("AS=", associated))
  vcf$FORMAT <- "GT"
  vcf <- cbind(vcf, geno)
  colnames(vcf)[1] <- "#CHROM"

  ## Write results
  vcf_header <- '##fileformat=VCFv4.1
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=AS,Number=1,Type=Integer,Description="Is SNP associated with response?">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  writeLines(vcf_header, "snps.vcf")
  write.table(vcf, "snps.vcf", append = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(clinical, "clinical.txt", row.names = FALSE, quote = FALSE, sep = "\t")

}
