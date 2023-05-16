snpFilePrep <- function(org) {
  file <-
    list.files(
      path = paste0(snpPath),
      pattern = paste0(org),
      full.names = TRUE
    )
  SNPfile <- read.csv(file)
  SNPfile <- SNPfile[, c("trait", "org", "chr", "bp")] # I don't have pvalues for Sorghum so not select pval yet
  # ordering--must have chromosomes and SNPs in numerical order
  SNPfile$chr <- as.integer(SNPfile$chr)
  SNPfile <-
    SNPfile[order(SNPfile[, "trait"], SNPfile[, "org"], SNPfile[, "chr"], SNPfile[, "bp"]), ]
  SNPfile$SNPend <- NA
  SNPfile$clpsRanges <- NA
  SNPfile$start <- NA
  SNPfile$end <- NA
  return(SNPfile)
}

collapsing <- function(snpTable, begin_loci = 1, done = F) {
  # not for sure why this has to be on repeat, I just know this won't work w/o it
  
  repeat {
    start_pos <- snpTable[(begin_loci), "bp"]
    for (j in 1:nrow(snpTable)) {
      current_pos <- snpTable[j, "bp"]
      next_pos <- snpTable[(j + 1), "bp"]
      
      # if statement to stop when last snp is reached
      
      if (j >= nrow(snpTable)) {
        done <- TRUE
        snpTable[(begin_loci), "clpsRanges"] <-
          as.numeric(current_pos) - snpTable[(begin_loci), "bp"]
        snpTable[(begin_loci), "start"] <-
          snpTable[(begin_loci), "bp"] - snpTable[(begin_loci), "range"]
        snpTable[(begin_loci), "end"] <-
          snpTable[(begin_loci), "bp"] + snpTable[(begin_loci), "clpsRanges"] + snpTable[(begin_loci), "range"]
        break
        
        # if not last snp, check to see if next snp is father than 2*range
        # begin_loci keeps track of the beginning of the current loci
        # if snps are < 2*range apart, it skips on and their clpsRanges/start/end remain as NAs to search for later
        # if snps are > 2*range apart, it records the clpsRanges/start/end info for the snp that begins the loci
        
      } else if ((as.numeric(next_pos) - as.numeric(current_pos)) > 2 *
                 snpTable[(begin_loci), "range"]) {
        snpTable[(begin_loci), "SNPend"] <- as.numeric(current_pos)
        snpTable[(begin_loci), "clpsRanges"] <-
          as.numeric(current_pos) - snpTable[(begin_loci), "bp"]
        snpTable[(begin_loci), "start"] <-
          snpTable[(begin_loci), "bp"] - snpTable[(begin_loci), "range"]
        snpTable[(begin_loci), "end"] <-
          snpTable[(begin_loci), "bp"] + snpTable[(begin_loci), "clpsRanges"] + snpTable[(begin_loci), "range"]
        
        # the beginning of the next loci is made to be the next_pos (==j+1) since it was > 2* range away from current snp
        
        begin_loci <- j + 1
      }
    }
    if (done == T) {
      break
    }
  }
  return(snpTable)
}

get_orthogroups<-function(x, tab_strs=og_strs){
  og<-gsub("^(OG\\d+).*","\\1", tab_strs[str_which(tab_strs, x)])
  if (!identical(og, character(0))) {
    return(data.frame(Orthogroup = og, ID = x))
  }
}