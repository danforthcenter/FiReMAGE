#Orthofinder local

################################################################################
### 0. Setting up ##############################################################
################################################################################

rm(list=ls())

packages <-
  c(
    #"docopt",
    "plyr",
    "dplyr",
    "tidyr",
    "data.table",
    "doParallel",
    "reader",
    "readr",
    "iterators",
    "ggplot2",
    "scales",
    "rbenchmark",
    "stringr"
  )
for (x in packages) {
  if (!require(x, character.only = TRUE, quietly = TRUE)) {
    install.packages(x, dep = TRUE)
    #if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

#library(docopt)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(reader)
library(readr)
library(iterators)
library(ggplot2)
library(scales)
library(rbenchmark)
library(stringr)

#set working dir
setwd("YOURPATH/FiReMAGE_v1/OrthoFinder_v/")

source("./scripts/functions_OrthoFinder_version.R")

# unique dir for traits, current code won't work if multiple files for an org exist

snpPath <- "./data/seedIonome_snps/"

# metaTable has chr and range info for ea org and allows users to select which orgs are in the comparison

metaTable <-
  fread(file = "./data/metaTables/V1_5s_metaTable.csv",
        sep = ",",
        stringsAsFactors = FALSE)

# ortholog file from OrthoFinder run

orthologs <-
  fread(
    file = "./data/OrthoFinder_orthologs/Orthogroups_5s.tsv",
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )
# collapsing orthogroups into strings speed up searches

og_strs<-apply(orthologs,1,paste0,collapse=",")

# permutations, 10 is usually good for testing

numPermutations <- 1000

## setting up cores for parallel processes, cl and cores always need to be the same number
# this is about how many I can use and still have other processes running in the background

nCore = 4
registerDoParallel(cores = nCore)


## doParallel can be funky depending on Windows/Linux/MacOS. 
## I use Windows locally and Linux on our server and doparallel works on both.
## If it wont work for you it's best to run nCore <- 1 so it won't mess w running in parallel

print(paste0("nCore working like it should?: ", getDoParWorkers() == nCore))

#name of your output dir

output <- "./results/"
dir.create(paste0(output))
dir.create(paste0(output, "graphs/"))
dir.create(paste0(output, "permutation_files/"))
dir.create(paste0(output, "permutation_files/snps/"))
dir.create(paste0(output, "permutation_files/gene_hits/"))
dir.create(paste0(output, "candidate_lists/")) 

## making sure metaTable is in alphabetical order for sanity's sake

metaTable <- metaTable[order(metaTable$orgs), ]

## gene coordinates for all species, used to find linked genes in the snps

# NA's from all_gene_coords is trying to change mitchondria (M) and chloroplast (C) chr to integers

all_gene_coords <-
  foreach(
    o = as.character(metaTable$orgs),
    .packages = packages.loaded(),
    .combine = rbind
  ) %do% {
    coords <-
      fread(
        file = paste0("./data/phytozome/", o, "_gene_coords.csv"),
        sep = ",",
        header = T,
        stringsAsFactors = F
      )
    coords[, `Chromosome Name` := gsub("Chr|Chr_", "", `Chromosome Name`),]
    coords$`Chromosome Name` <- as.integer(coords$`Chromosome Name`)
    coords$org <- o
    return(coords)
  }

rm(orthologs)

################################################################################
### 1.1 Actual data: collapsing snps into loci
################################################################################

## even if there are other species in the snp dir, should only read those listed in the metaTable

snps <-
  foreach(
    q = as.character(metaTable$orgs[!metaTable$input_as_loci]),
    .packages = packages.loaded(),
    .combine = rbind
  ) %do% {
    SNPfile <- snpFilePrep(q)
    SNPfile$range <- metaTable$range[metaTable$orgs == q]
    return(SNPfile)
  }

## collapsing SNPs: looks to see if 2 SNPs have a distance between them of more than 2*the LD (the "range").
## If not: removes the SNPs subsequent to the first one until it finds one that is > 2*the LD.

snps_split <- split(snps, list(snps$trait, snps$org, snps$chr)) # makes list of snps for each trait/org/chr combo

CollapsedSNPs <-
  foreach(s = snps_split,
          .packages = packages.loaded(),
          .combine = rbind) %do% {
            # check for empty tables, errors halt execution
            
            if (!nrow(s) == 0) {
              collapsed_s <- collapsing(s)
              return(collapsed_s)
              
            }
          }

# appending datasets input as loci

for(m in metaTable$orgs[metaTable$input_as_loci]){
  file <-
    list.files(
      path = paste0(snpPath),
      pattern = paste0(m),
      full.names = TRUE
    )
  SNPfile <- read.csv(file)
  SNPfile<-SNPfile[,c("trait","org","chr","bp","start","end")]
  ###ordering--must have chromosomes and SNPs in numerical order
  SNPfile$chr<-as.integer(SNPfile$chr)
  SNPfile<-SNPfile[order(SNPfile[,"trait"], SNPfile[,"org"], SNPfile[,"chr"], SNPfile[, "end"]),]
  SNPfile$SNPend<-NA
  SNPfile$clpsRanges<-SNPfile$end - SNPfile$start
  SNPfile$range <- metaTable$range[metaTable$orgs == m]
  CollapsedSNPs<-rbind(CollapsedSNPs,SNPfile)
}

# SNPs within 2*range of the beginning of loci are left with clpsRanges at NA, so we remove them to leave only loci

snps_sub <- CollapsedSNPs[!is.na(CollapsedSNPs$clpsRanges), ]

# The first loci will get it's start set at bp - range, which will be negative, so set those back to 1

snps_sub$start[snps_sub$start < 1] <- 1

# get rid of SNPend and clpsranges columns and reordering columns, info still in CollapsedSNPs for sanity checks if needed

snps_sub <- data.frame(snps_sub[, c("org", "chr", "bp", "trait")], apply(snps_sub[, c("start", "end")], 2, ceiling))

# making unique loci names

snps_sub$loci <- paste0(snps_sub$org,
                        "_",
                        snps_sub$trait,
                        "_",
                        snps_sub$chr,
                        "_",
                        snps_sub$bp)

# counting collapsed snps for each org/trait combo

leafcollapsedSNPs <- data.frame(unclass(table(snps_sub$trait, snps_sub$org)))
leafcollapsedSNPs <- leafcollapsedSNPs[, order(colnames(leafcollapsedSNPs))]
leafcollapsedSNPs$trait <- rownames(leafcollapsedSNPs)
rownames(leafcollapsedSNPs) <- c()

write.csv(
  leafcollapsedSNPs,
  file = paste0(output, "summary_collapsedSNPs.csv"),
  row.names = F
)

################################################################################
### 1.2 Filtering linked genes for orthologs
################################################################################

setDT(snps_sub)
setkey(snps_sub, org, chr, start, end)

snp_gene_hitTable <-
  foverlaps(
    all_gene_coords,
    snps_sub,
    by.x = c("org", "Chromosome Name", "Gene Start (bp)", "Gene End (bp)"),
    by.y = c("org", "chr", "start", "end"),
    type = "any",
    nomatch = NULL
  )
write.csv(
  snp_gene_hitTable,
  file = paste0(output, "snp_gene_hitTable.csv"),
  row.names = F
)

OrthoMerge <-
  foreach(
          i = as.character(unique(snp_gene_hitTable$`Gene Name`)),
          .packages = packages.loaded(),
          .combine = rbind
        ) %dopar% {
          og<-orthologs$Orthogroup[as.integer(unlist(Map(grep, paste(i), orthologs)))]
          if (!identical(og, character(0))) {
            return(cbind(snp_gene_hitTable[snp_gene_hitTable$`Gene Name` == i,],
                         data.frame(Orthogroup = og)))
          }
        }

## quick sanity check to make sure orthogroups merged with their genes correctly, I'll leave it here for troubleshooting

# backend<- foreach(i=1:nrow(OrthoMerge), .combine = c) %dopar% {
#   OrthoMerge$Orthogroup[i]==orthologs$Orthogroup[as.integer(unlist(Map(grep, paste(OrthoMerge$`Gene Name`[i]), orthologs)))]
# }
 
##

# present is what I named the col for counting the number of species present in a ortho group
# present is used later to rank genes

speciesCount<-as.data.frame(table(OrthoMerge$Orthogroup, OrthoMerge$org, OrthoMerge$trait))
speciesCount$Freq<-as.logical(speciesCount$Freq)
present<-aggregate(speciesCount$Freq, by= list(speciesCount$Var1, speciesCount$Var3), FUN = sum)
colnames(present)<-c("Orthogroup","trait","present")

OrthoMerge<-merge(OrthoMerge, present, by=c("Orthogroup","trait"))

# focusing on conserved genes, but this could be changed for users if they want stuff between 2 orgs

OrthoMerge <- OrthoMerge[OrthoMerge$present > 2,]

genecount <-
  foreach(s = 1:nrow(snps_sub), .packages = packages.loaded()) %dopar% {

    checking_genes <- length(as.character(unique(OrthoMerge$`Gene Name`[OrthoMerge$loci==snps_sub$loci[s]])))

    return(checking_genes)
  }
snps_sub$genecount <- unlist(genecount)

write.csv(snps_sub, file=paste0(output, "snps_subset.csv"), row.names=F)
write.table(OrthoMerge, file = paste0(output,"Orthogroup_hits.csv"), col.names = T, row.names = F, sep = ",")

rm(snps,snp_gene_hitTable)
gc()

################################################################################
### 2.1 Creating random permutations
################################################################################

## reading in chr coordinates for orgs in metaTable
## check to make sure the file has been updated to accomodate new orgs

chrLengths <-
  read.table(
    file = "./data/org_chromosome_coords.csv",
    header = TRUE,
    sep = ",",
    na.strings = NA,
    stringsAsFactors = FALSE
  )

chrLengths <- as.list(chrLengths)
chrLengths <- lapply(chrLengths, function(x) x[!is.na(x)])

##each row in leafcollapsedSNPs is a trait

AllPermuts <-
  foreach(
    row = iter(leafcollapsedSNPs, by = 'row'),
    .combine = 'rbind',
    .packages = packages.loaded()
  ) %dopar% {
    print(paste("Starting permutations for", row$trait))
    
    # tells the permutations how many chr and snps per org to replicate for this trait
    orgDetails <- data.frame(
      org = metaTable$orgs,
      nChrs = metaTable$nChrs,
      nSNPs = c(as.numeric(row[1, c(1:nrow(metaTable))]))
    )
    permuteDataset <- with(orgDetails, {
      return(data.frame(rbindlist(lapply(1:numPermutations, function(x) {
          # initiate table
          
          thisPermute <-
            data.table(
              permutation = numeric(0),
              org = character(0),
              chr = numeric(0),
              bp = numeric(0),
              stringsAsFactors = FALSE
            )
          for (j in 1:nrow(orgDetails)) {
            # j == organism iterator
            # this is meant to handle traits that don't exist in all org datasets
            # for example, rice is missing cobalt, but I want to know what the results are for the other 4 orgs
            
            if (nSNPs[j] == 0) {
              thisChr <- NA
              
            } else{
              # randomly scatters snps among chrs (but doesn't select bp yet)
              
              thisChr <- sort(sample(nChrs[j], nSNPs[j], replace = TRUE))
              
              # randomly assigns a bp depending on the chr selected above
              # accounts for each chr's length when assigning the largest bp value
              
              thisPermute <-
                with(chrLengths, {
                  rbindlist(list(
                    thisPermute,
                    data.table(
                      permutation = x,
                      org = org[j],
                      chr = thisChr,
                      bp = unlist(lapply(1:nChrs[j], function(y)
                            sample(get(as.character(org[j]))[y], 
                                   length(which(thisChr == y)), 
                                   replace = FALSE))),
                      stringsAsFactors = FALSE
                    )
                  ))
                })
            }
          }
          return(thisPermute)
        })
      )))
    })  #end permutation (with() loop)
    permuteDataset$trait <- row$trait
    
    #permutation snps also get unique loci identifiers
    
    permuteDataset$loci <- paste0(permuteDataset$org,
                                  "_",
                                  permuteDataset$trait,
                                  "_",
                                  permuteDataset$chr,
                                  "_",
                                  c(1:nrow(permuteDataset))
                                )
    return(permuteDataset)
  }
rm(chrLengths)

## getting ranges from real data

RandomLociRanges <- CollapsedSNPs[!is.na(CollapsedSNPs$clpsRanges), ]

## same split as before with the additional permutation column

AllPermuts_split <-
  split(AllPermuts,
        list(AllPermuts$org, AllPermuts$trait, AllPermuts$permutation))

## randomly shuffles ranges from each org/trait combo in each permutation

AllPermuts <-
  foreach(A = AllPermuts_split,
          .combine = rbind,
          .packages = packages.loaded()) %do% {
            c <- sample(RandomLociRanges$clpsRanges[RandomLociRanges$trait == A$trait[1] &
                                                      RandomLociRanges$org == A$org[1]],
                        size = nrow(A),
                        replace = F)
            A$start <-
              A$bp - (ceiling(.5 * c) + metaTable$range[metaTable$orgs == A$org[1]])
            
            A$end <-
              A$bp + (ceiling(.5 * c) + metaTable$range[metaTable$orgs == A$org[1]])
            
            return(A)
          }

## if the range randomly happens to hang over a chr beginning we need to reset it back to 1
## if the range hangs over a chr end it won't matter since there won't be any genes on that chr past the end point

AllPermuts$start[AllPermuts$start < 0] <- 1

## permutations need to be in order for foverlaps, still records permutations so this is fine

AllPermuts <-
  AllPermuts[order(AllPermuts[, "org"], AllPermuts[, "chr"], AllPermuts[, "bp"]), ]

# at 1000 permutations, this file gets big, so best to process in chunks
# if your datasets are much larger and would like to chunk it up more to save memory,
# this is where I would edit

if(numPermutations > 100) {
  for (i in seq(100, numPermutations, by = 100)) {
    permut_chunk <- AllPermuts[AllPermuts$permutation %in% c((i - 99):i), ]
    
    write.table(
      permut_chunk,
      file = paste0(
        output,
        "permutation_files/snps/",
        i - 99,
        "_",
        i,
        "_permutations.csv"
      ),
      sep = ",",
      col.names = TRUE,
      row.names = FALSE
    )
  }
  rm(permut_chunk)
} else{
  write.table(
    AllPermuts,
    file = paste0(output, "permutation_files/snps/",min(AllPermuts$permutation),"_",max(AllPermuts$permutation),"_permutations.csv"),
    sep = ",",
    col.names = T,
    row.names = F
  )
}
rm(AllPermuts)

################################################################################
### 2.2 Filtering linked genes for orthologs in permutation data
################################################################################

backend <-
  foreach(
    e = list.files(
      paste0(output, "permutation_files/snps"),
      pattern = "permutations.csv",
      full.names = T
    ),
    .packages = packages.loaded()
  ) %dopar% {
    p <- fread(e, sep = ",", header = T, stringsAsFactors = F)
    setDT(p)
    setkey(p, org, chr, start, end)
    
    permutation_gene_hitTable <-
      foverlaps(
        all_gene_coords,
        p,
        by.x = c("org", "Chromosome Name", "Gene Start (bp)", "Gene End (bp)"),
        by.y = c("org", "chr", "start", "end"),
        type = "any",
        nomatch = NULL
      )
    
    ## permutation_gene_hitTable was mostly useful for me when making graphs
    ##it's a large file so if not planning on making graphs and short on mem it can be commented out
    
    write.csv(
      permutation_gene_hitTable,
      file = paste0(
        output,
        "permutation_files/gene_hits/",
        min(p$permutation),
        "_",
        max(p$permutation),
        "_gene_hitTable.csv"
      ),
      row.names = F
    )
    gc()
  }

rm(RandomLociRanges)

print("Random dataset ...")

## the random dataset permutation ortholog file comparison is the most time consuming part of FiReMAGE

print(paste("Beginning permutation ortholog search:", Sys.time()))

for(e in list.files(
  paste0(output, "permutation_files/gene_hits"),
  pattern = "_hitTable.csv",
  full.names = T
)) {
  
  gene_hitTable <- fread(e, sep = ",", stringsAsFactors = F, header = T)
  
  Trait_split <-
    split(
      gene_hitTable,
      gene_hitTable$trait)
  
  ## genecounting like we did for Actual doesn't scale for permutations
  ## so I break it out by trait and save the tables at the end to read in for making gene/loci summaries
  ## stuff also gets funky when i get to a third layer of foreach loops, so this very outer one is a regular for loop
  
  for (e in Trait_split) {
    print(paste0(e$trait[1], " permutation merging"))
    
    PermutationMerge <-
      
      foreach(
        i = as.character(unique(e$`Gene Name`)),
        .packages = packages.loaded(),
        .combine = rbind
      ) %dopar% {
        og<-orthologs$Orthogroup[as.integer(unlist(Map(grep, paste(i), orthologs)))]
        if (!identical(og, character(0))) {
          return(cbind(e[e$`Gene Name` == i, ], data.frame(Orthogroup = og)))
        }
      }
    PermutationMerge$Orthogroup <- as.character(PermutationMerge$Orthogroup) 
    # present is what I named the col for counting the number of species present in a ortho group
    # present is used later to rank genes
    
    speciesCount <- as.data.frame(table(PermutationMerge$Orthogroup, PermutationMerge$org, PermutationMerge$permutation))
    speciesCount$Freq <- as.logical(speciesCount$Freq)
    present <- aggregate(speciesCount$Freq, by = list(speciesCount$Var1, speciesCount$Var3), FUN = sum)
    colnames(present) <- c("Orthogroup", "permutation","present")
    present$permutation <- as.integer(as.character(present$permutation))
    present$Orthogroup <- as.character(present$Orthogroup)
   
    PermutationMerge <- merge(PermutationMerge, present, by = c("Orthogroup","permutation"))
    
    # focusing on conserved genes, but this could be changed for users if they want stuff between 2 orgs
    
    PermutationMerge <- PermutationMerge[PermutationMerge$present > 2, ]
    
    # I write this out to use later for making graphs
    
    write.csv(
      PermutationMerge,
      file = paste0(
        output,
        "/permutation_files/",
        min(e$permutation),
        "_",
        max(e$permutation),
        "_",
        e$trait[1],
        "_PermutationMerge.csv"
      ),
      row.names = F
    )
    
    rm(PermutationMerge)
  }
}

print(paste("End permutations", Sys.time()))

################################################################################
### 3. Using permutation distributions to assess the actual loci
################################################################################

print("Summarizing actual data")
 
## counts number of genes & loci for each trait/org/# of species present in ortho group

Final_summary<-aggregate(
  cbind(`Gene Name`, loci) ~ trait + present + org,
  data = OrthoMerge,
  FUN = function(x)
    length(unique(x))
)
setnames(Final_summary, old = c("Gene Name", "loci"), new = c("GeneCount", "LociCount"))
 
Final_summary$dataset <- "Actual"

## to save memory, this does the same thing as the real but I read in each trait file separately
 
Permutation_summary <-
  foreach(
    e = list.files(
      paste0(output, "permutation_files/"),
      pattern = "PermutationMerge.csv",
      full.names = T
    ),
    .packages = packages.loaded(),
    .combine = rbind
  ) %dopar% {
    print(paste0("Summarizing ", e))
    f <- fread(e, sep = ",", stringsAsFactors = FALSE)
    if (!nrow(f) == 0) {
      fsummary <- aggregate(
        cbind(`Gene Name`, loci) ~ trait + present + org + permutation,
        data = f,
        FUN = function(x)
          length(unique(x))
      )
      setnames(fsummary, old = c("Gene Name", "loci"), new = c("GeneCount", "LociCount"))
      
      return(fsummary)
    }
  }
Permutation_summary$dataset <- "Permutations"

## only permutations that have a hit are represented in Permutation_summary at this point
## and it isn't accounting for all the permutations with no hits
## the following split and foreach loop will find the n permutations missing, and
## add n observations with value "0", this way quantiles can be calculated correctly

Permutation_summary_split <-
   split(
    Permutation_summary,
    list(
      Permutation_summary$org,
      Permutation_summary$trait,
      Permutation_summary$present
    )
  )
 
Permutation_summary <-
  foreach(p = Permutation_summary_split,
          .packages = packages.loaded(),
          .combine = rbind) %do% {
            if ((!nrow(p) == 0) & (nrow(p) < numPermutations)){
              d <-
                data.frame(
                  org = rep(p$org[1], c(numPermutations - nrow(p))),
                  trait = rep(p$trait[1], c(numPermutations - nrow(p))),
                  present = rep(p$present[1], c(numPermutations - nrow(p))),
                  dataset = rep("Permutations", c(numPermutations - nrow(p))),
                  GeneCount = rep(0, c(numPermutations - nrow(p))),
                  LociCount = rep(0, c(numPermutations - nrow(p)))
                )
              p$permutation<-NULL
              # sanity check  
              #print(sum(nrow(d),nrow(p)))
              return(rbind(p, d))
            }else{
              if(nrow(p)==numPermutations){
                p$permutation<-NULL
                return(p)
                }
            }
          }

AllSummaries <- rbind(Final_summary, Permutation_summary)
AllSummaries$GeneCount <- as.numeric(AllSummaries$GeneCount)
AllSummaries$LociCount <- as.numeric(AllSummaries$LociCount)

write.csv(AllSummaries, file = paste0(output, "AllSummaries.csv"), row.names = F)

## calculating some quantiles for the permutations and setting actual quantiles to NA

graphingDF <-
  ddply(
    AllSummaries,
    .(org, trait, dataset, present),
    summarize,
    GenesMean = mean(GeneCount),
    LociMean = mean(LociCount),
    Gene95 = quantile(GeneCount, type=3, .95),
    Loci95 = quantile(LociCount, type=3, .95),
    Gene99 = quantile(GeneCount, type=3, .99),
    Loci99 = quantile(LociCount, type=3, .99),
    Gene05 = quantile(GeneCount, type=3, .05),
    Loci05 = quantile(LociCount, type=3, .05),
    Gene01 = quantile(GeneCount, type=3, .01),
    Loci01 = quantile(LociCount, type=3, .01)
  )

graphingDF$Gene95[graphingDF$dataset == "Actual"] <- NA
graphingDF$Loci95[graphingDF$dataset == "Actual"] <- NA
graphingDF$Gene99[graphingDF$dataset == "Actual"] <- NA
graphingDF$Loci99[graphingDF$dataset == "Actual"] <- NA
graphingDF$Loci01[graphingDF$dataset == "Actual"] <- NA
graphingDF$Gene01[graphingDF$dataset == "Actual"] <- NA
graphingDF$Gene05[graphingDF$dataset == "Actual"] <- NA
graphingDF$Loci05[graphingDF$dataset == "Actual"] <- NA
 
## currently, I make graphs for each set of ortholog group species representation (ie. 3/5, 4/5, 5/5)
## so for now I make 3 different graphs for genes with CI 95 and genes with CI 99
## if you have any idea how to summarize the everything into one graph I'd be interested to know

graphingDF$trait<-as.factor(graphingDF$trait)
graphing_split <- split(graphingDF, graphingDF$present)

backend <-
  foreach(data = graphing_split, .packages = packages.loaded()) %do% {
    G95 <-
      ggplot(data,
             aes(
               x = GenesMean,
               y = trait,
               group = dataset,
               color = dataset
             )) +
      geom_errorbarh(aes(
        xmax = Gene95,
        xmin = Gene05,
        height = 0.75
      )) +
      geom_point(aes(
        color = dataset,
        size = dataset
      )) +
      scale_size_manual(values = c(4, 3)) +
      labs(x = "Overlapped Genes Returned",
           title = paste0("Ortholgs in groups with ",data$present[1],"/",nrow(metaTable)," species representation")) +
      theme_bw(base_size = 18) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#7fcdbb", "#2c7fb8")) +
      facet_grid( ~org, scales = "free")

    G99 <-
      ggplot(data,
             aes(
               x = GenesMean,
               y = trait,
               group = dataset,
               color = dataset
             )) +
      geom_errorbarh(aes(
        xmax = Gene99,
        xmin = Gene01,
        height = 0.75
      )) +
      geom_point(aes(
        color = dataset,
        size = dataset
      )) +
      scale_size_manual(values = c(4, 3)) +
      labs(x = "Overlapped Genes Returned",
           title = paste0("Ortholgs in groups with ",data$present[1],"/",nrow(metaTable)," species representation")) +
      theme_bw(base_size = 18) +
      scale_color_manual(values = c("#7fcdbb", "#2c7fb8")) +
      facet_grid( ~org, scales = "free")
 
    ggsave(
      paste0(
        output,
        "graphs/Species_combo_",
        data$present[1],
        "_Genes95.png"
      ),
      G95,
      height = 9,
      width = 16
    )

    ggsave(
      paste0(
        output,
        "graphs/Species_combo_",
        data$present[1],
        "_Genes99.png"
      ),
      G99,
      height = 9,
      width = 16
    )
   }

## making the candidateLists for each organism

candidateList <-
  foreach(
    o = as.character(metaTable$orgs),
    .packages = packages.loaded(),
    .combine = rbind
  ) %do% {
    sub_OrthoMerge <-
      subset(OrthoMerge, org==o)

    # getting all the genes with a loci hit, and getting the gene/loci counts added on

    subList <-
      distinct(merge(sub_OrthoMerge,
                     snps_sub[, c("loci", "genecount")],
                     by = "loci"))

    # bc we think that a loci should be corresponding to only one true causal gene (usually),
    # we penalize genes coming from loci with high genecounts

    subList$gFDR <- (1 - (1 / (subList$genecount)))

    # use the info from the summaries above to calculate FDR from the permutations
    # it's long lines, but it's avg # of genes from permutations / avg # of genes from real data for each org/trait/# of species present in ortho group

    pFDR_backend <-
      foreach(
        row = iter(subList, by = 'row'),
        .packages = packages.loaded(),
        .combine = rbind
      ) %do% {
        pFDRvalue <-
          (graphingDF$GenesMean[graphingDF$org == row$org &
                                  graphingDF$trait == row$trait &
                                  graphingDF$present == row$present &
                                  graphingDF$dataset == "Permutations"] /
             graphingDF$GenesMean[graphingDF$org == row$org &
                                    graphingDF$trait == row$trait &
                                    graphingDF$present == row$present &
                                    graphingDF$dataset == "Actual"])

        if (identical(pFDRvalue, numeric(0))) {
          pFDRvalue <- 0
        }
        return(data.frame(pFDR = pFDRvalue))
      }

    subList$pFDR <- pFDR_backend$pFDR
    return(subList)
  }

## ranking genes
candidateList$rank <-
  ((1 - candidateList$gFDR) * (1 - candidateList$pFDR) * (candidateList$present / nrow(metaTable)))

## we thought gFDR was too harsh on some of our known ionomics genes so we took it down a bit
## the ranges used are pretty wide for ea species so it's likely to hit multiple genes, pFDR seemed more important
## might do something similar to the total # of species represented in an ortho group

candidateList$rank_new <-
  ((1 - (candidateList$gFDR * 0.2)) * (1 - (candidateList$pFDR)) * (candidateList$present / nrow(metaTable)))

## Dividing the candidate list up by trait

candidateList_split <- split(candidateList, candidateList$trait)
for (t in candidateList_split) {
  individualcandidate <- arrange(t, desc(rank_new))
  
  write.csv(
    individualcandidate,
    file = paste0(output, "candidate_lists/", t$trait[1], "_candidateList.csv"),
    row.names = FALSE
  )
}
