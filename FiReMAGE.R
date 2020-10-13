#!/usr/bin/Rscript

"FiReMAGE script
Usage: CGASconcise.R -s <path> -m <path> -f <path> [-p <int> -c <int> -o <path>]
Options:
-h, --help                      Show this screen.
-s <path>, --snps <path>        Path to the directory containing the SNP files for the comparison
-m <path>, --metaTable <path>   Path to read in the metadata table of organisms in the comparison
-f <path>, --orthologFiles <path> Path to read in ortholog files for the comparison
-o <path>, --output <path>      Files will be written to this path [default: './FiReMAGE_output']
-p <int>, --permutations <int>  Number of iterations [default: 1000].
-c <int>, --cores <int>         Number of cores to use for parallel loop [default: 1]
" -> doc

library(docopt)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(reader)
library(readr)
library(iterators)
library(ggplot2)
library(ionomicsUtils)
library(scales)

library(rbenchmark)
source("./functions.R")

###docopt options, for use in terminal
#opts <- docopt(doc)
#str(opts)
#snpPath<-opts[["--snps"]]
#orthoFilePath<-opts[["--orthologFiles"]]
#metaTable<-fread(opts[["--metaTable"]], sep=",", stringsAsFactors = FALSE)
#output<-opts[["--output"]]
#numPermutations <- opts[["--permutations"]]
#nCore <- opts[["--cores"]]

###editable options, for use in RStudio GUI
snpPath<-"./data/el_snps/current/ionomics_6species/"
orthoFilePath<-"./data/phytozome/current/"
metaTable<-fread(file = "./data/metaTables/ArabSoySorgCornMetaTable.csv", sep = ",", stringsAsFactors = FALSE)
output<-"./data/FiReMAGEtesting/seedIonome_Test/"
numPermutations<-100
nCore<-4

#######beginning pipeline

###setting up cores for parallel processes, cl and cores always need to be the same number

registerDoParallel(cores = nCore)

dir.create(paste0(output))
dir.create(paste0(output,"graphs/"))
dir.create(paste0(output,"permutation_files/"))
dir.create(paste0(output,"priority_lists/"))

metaTable <- metaTable[order(metaTable$orgs),]

####snp files

snps<-foreach(q=as.character(metaTable$orgs), .packages=packages.loaded(), .combine=rbind) %do% {
  SNPfile<-snpFilePrep(q)
  SNPfile$range<-metaTable$range[metaTable$orgs==q]
  return(SNPfile)
}
all_gene_coords<-foreach(o=as.character(metaTable$orgs), .packages=packages.loaded(), .combine = rbind) %do% {
  coords<-fread(file=paste0("./data/phytozome/current/",o,"_gene_coords.csv"),sep = ",",header = T,stringsAsFactors = F)
  coords[,`Chromosome Name` := gsub("Chr|Chr_","",`Chromosome Name`),]
  coords$`Chromosome Name`<-as.integer(coords$`Chromosome Name`)
  coords$org<-o
  return(coords)
}
snps$trait<-gsub("[0-9]*","",snps$trait)
traits<-as.character(unique(snps$trait))

###collapsing SNPs: looks to see if 2 SNPs have a distance between them of more than 2*the LD (the "range"). If not: removes the 
###SNPs subsequent to the first one until it finds one that is > 2*the LD. 
snps$SNPend<-NA
snps$clpsRanges<-NA
snps$start<-NA
snps$end<-NA
rows_to_remove<-NA

CollapsedSNPs<-collapsing(snps)

rows_to_remove<-rows_to_remove[!is.na(rows_to_remove)]
snps_sub<-CollapsedSNPs[-rows_to_remove,]
snps_sub$start[snps_sub$start < 1]<-1
#gets rid of SNPend and clpsranges
snps_sub<-data.frame(snps_sub[,c("org", "chr", "bp", "trait")], apply(snps_sub[,c("start", "end")], 2, ceiling))
snps_sub$loci<-paste0(snps_sub$org,"_",snps_sub$trait,"_",snps_sub$chr,"_",c(1:nrow(snps_sub)))

#seeing how many collapsed snps each organism has/trait
leafcollapsedSNPs<-data.frame(unclass(table(snps_sub$trait, snps_sub$org)))
leafcollapsedSNPs<-leafcollapsedSNPs[,order(colnames(leafcollapsedSNPs))]
leafcollapsedSNPs$trait<-rownames(leafcollapsedSNPs)
rownames(leafcollapsedSNPs)<-c()
write.csv(leafcollapsedSNPs, file=paste0(output,"leafcollapsedSNPsTESTING.csv"), row.names = F)

setDT(snps_sub)
setkey(snps_sub,org,chr,start,end)
snp_gene_hitTable<-foverlaps(all_gene_coords,snps_sub, by.x=c("org","Chromosome Name","Gene Start (bp)","Gene End (bp)"),
                             by.y = c("org","chr","start","end"), type = "any", nomatch = NULL)
write.csv(snp_gene_hitTable, file = paste0(output,"snp_gene_hitTable.csv"), row.names = F)
####ortholog file prep
genefile<-GeneFilePrep(orthoFilePath)
genefile[is.na(genefile)]<-0
genefile<-subset(genefile, select = c(colnames(genefile)[-grep("Org",colnames(genefile))]))

chunk<-10000
n<-nrow(genefile)
r<-rep(1:ceiling(n/chunk),each=chunk)[1:n]
d<-split(genefile,r)
rm(genefile)
snp_split<-split(subset(snp_gene_hitTable, select = c(org,loci,trait, `Gene Name`)), snp_gene_hitTable$trait)
print("Starting overlap...")
FinalMerge<-foreach(geneTable=d, .packages = packages.loaded(), .combine=rbind) %:%
  foreach(snps=snp_split, .packages = packages.loaded(), .combine = rbind) %dopar% {
    #print(snps$trait[1])
    #geneTable<-geneTable
    for(o in as.character(unique(snp_gene_hitTable$org))) {
      if(nrow(subset(snps, org==o, select=c(loci, `Gene Name`)))==0){
        #print(paste0("No snps for ",o))
        geneTable<-cbind(geneTable, data.frame(loci=NA))
        setnames(geneTable, old = c("loci"), new = c(paste0(o,"_loci")))
      }else{
        geneTable<-merge(geneTable,subset(snps, org==o, select=c(loci, `Gene Name`)), by.x = paste0(o,"ID"), by.y = "Gene Name", all.x = T, allow.cartesian = T)
        setnames(geneTable, old = c("loci"), new = c(paste0(o,"_loci")))
      }
    }
    geneTable$present<-(nrow(metaTable)-rowSums(is.na(geneTable)))
    geneTable<-geneTable[geneTable$present > 2,]
    geneTable$trait<-snps$trait[1]
    return(geneTable)
  }

###gene per loci counting
print("Starting Actual Gene Counting")
genecount<-foreach(s=1:nrow(snps_sub), .packages = packages.loaded()) %dopar% {
  lociHit<-subset(snp_gene_hitTable, loci==snps_sub$loci[s])
  if(nrow(lociHit)==0){
    return(0)
  }else{
    checking<-FinalMerge[,c(paste0(snps_sub$org[s],"ID"),paste0(snps_sub$org[s],"_loci")),with=FALSE]
    colnames(checking)<-c("ID","loci")
    checking<-subset(checking, loci==snps_sub$loci[s])
    return(length(lociHit$`Gene Name`[which(lociHit$`Gene Name` %in% checking$ID)]))
  }
}
snps_sub$genecount<-unlist(genecount)

write.csv(snps_sub, file=paste0(output,"snps_sub.csv"), row.names = F)
write.csv(FinalMerge, file=paste0(output,"FinalMerge.csv"), row.names = F)
rm(CollapsedSNPs, snps, rows_to_remove)

#####Random_datasets.R permutations#####

####reading in chromosome coordinates for orgs in metaTable, check to make sure the file has been updated to accomodate new orgs

chrLengths<-read.table(file = "./data/org_chromosome_coords.csv", header = TRUE, sep = ",", na.strings = NA, stringsAsFactors = FALSE)
chrLengths<-as.list(chrLengths)
chrLengths<-lapply(chrLengths, function(x) x[!is.na(x)])
AllPermuts <- foreach(row=iter(leafcollapsedSNPs,by='row'), .combine='rbind', .packages=packages.loaded()) %dopar% {
  print(paste("Starting permutations for",row$trait))
  orgDetails <- data.frame(org=metaTable$orgs,nChrs=metaTable$nChrs,
                           nSNPs=c(as.numeric(row[1,c(1:nrow(metaTable))])))
  permuteDataset <- with(orgDetails,{
    return(data.frame(rbindlist(lapply(1:numPermutations,function(x){
      thisPermute <- data.table(permutation=numeric(0),org=character(0),chr=numeric(0),bp=numeric(0),stringsAsFactors = FALSE)
      for(j in 1:nrow(orgDetails)){
        if(nSNPs[j]==0){     ####this is meant to handle traits that don't exist in all datasets. Useful for ionomic datasets
          thisChr<-NA
        }else{
          thisChr <- sort(sample(nChrs[j],nSNPs[j],replace=TRUE))
          thisPermute <- with(chrLengths, {rbindlist(list(thisPermute,data.table(permutation=x,org=org[j],chr=thisChr,
                                                                                 bp=unlist(lapply(1:nChrs[j],function(y) sample(get(as.character(org[j]))[y],length(which(thisChr==y)),replace = FALSE))),stringsAsFactors = FALSE)))})
        }
      }
      return(thisPermute)
    }))))
  })  #end permutation (with() loop)
  permuteDataset$trait <- row$trait
  permuteDataset$loci<-paste0(permuteDataset$org,"_",permuteDataset$trait,"_",permuteDataset$chr,"_",c(1:nrow(permuteDataset)))
  return(permuteDataset)
}#end foreach phenotype
rm(chrLengths)

#####Random_searchNEW.R#####

Actual_snp<-foreach(q=as.character(metaTable$orgs), .packages=packages.loaded(), .combine=rbind) %do% {
  SNPfile<-snpFilePrep(q)
  SNPfile$range<-metaTable$range[metaTable$orgs==q]
  return(SNPfile)
}
Actual_snp$new_ranges<-Actual_snp[,"range"]
Actual_snp$loci<-paste0(Actual_snp$org,"_",Actual_snp$trait,"_",Actual_snp$chr,"_",c(1:nrow(Actual_snp)))

RandomLociRanges<-collapsing(Actual_snp,permutations=T) 
for(i in metaTable$orgs){
  #set.seed(1) #making random selections from the actual ranges that are reproducible
  AllPermuts$start[AllPermuts$org==i] <- AllPermuts$bp[AllPermuts$org==i]-ceiling(.5*sample(RandomLociRanges$new_ranges[RandomLociRanges$org==i],size = nrow(AllPermuts[AllPermuts$org==i,]),replace = T))
  #set.seed(1)
  AllPermuts$end[AllPermuts$org==i] <- AllPermuts$bp[AllPermuts$org==i]+ceiling(.5*sample(RandomLociRanges$new_ranges[RandomLociRanges$org==i],size = nrow(AllPermuts[AllPermuts$org==i,]),replace = T))
}
AllPermuts$start[AllPermuts$start < 0 ] <- 1
AllPermuts<-AllPermuts[order(AllPermuts[,"org"], AllPermuts[,"chr"], AllPermuts[,"bp"]),]
write.table(AllPermuts,file = paste0(output,"AllPermuts.csv"),sep=",",col.names = TRUE,row.names = FALSE)

setDT(AllPermuts)
setkey(AllPermuts,org,chr,start,end)
permutation_gene_hitTable<-foverlaps(all_gene_coords,AllPermuts, by.x=c("org","Chromosome Name","Gene Start (bp)","Gene End (bp)"),
                                     by.y = c("org","chr","start","end"), type = "any", nomatch = NULL)
write.csv(permutation_gene_hitTable,file = paste0(output,"permutation_gene_hitTable.csv"),row.names = F)

#AllPermuts_split<-split(AllPermuts, list(AllPermuts$permutation, AllPermuts$trait))
#AllPermuts_split<-split(subset(permutation_gene_hitTable, select = c(permutation,org,loci,trait, `Gene Name`)), 
                        #list(permutation_gene_hitTable$permutation, permutation_gene_hitTable$trait))
rm(Actual_snp,RandomLociRanges)
print("Random dataset ...")
####editing this chunk~~~~~~~~
#tracking<-1
print(paste("Begin permutations",Sys.time()))
Trait_split<-split(subset(permutation_gene_hitTable, select=c(permutation,org,loci,trait, `Gene Name`)), permutation_gene_hitTable$trait)
permutation_genecounting<-NULL
for(e in Trait_split){
  print(paste0(e$trait[1]," permutation merging"))
  byPermutation_split<-split(e, e$permutation)
  PermutationMerge<-foreach(geneTable=d, .packages = packages.loaded(), .combine = rbind) %:%
    foreach(snps=byPermutation_split, .packages = packages.loaded(), .combine = rbind) %dopar% {
      for(o in as.character(unique(snp_gene_hitTable$org))) {
        if(nrow(subset(snps, org==o, select=c(loci, `Gene Name`)))==0){
          #print(paste0("No snps for ",o))
          #geneTable<-cbind(geneTable, data.frame(loci=NA))
          geneTable[, loci := NA]
          setnames(geneTable, old = c("loci"), new = c(paste0(o,"_loci")))
        }else{
          geneTable<-merge(geneTable,subset(snps, org==o, select=c(loci, `Gene Name`)), by.x = paste0(o,"ID"), by.y = "Gene Name", all.x = T, allow.cartesian = T)
          setnames(geneTable, old = c("loci"), new = c(paste0(o,"_loci")))
        }
      }
      geneTable[, present := (nrow(metaTable)-rowSums(is.na(geneTable)))]
      #geneTable$present<-(nrow(metaTable)-rowSums(is.na(geneTable)))
      geneTable<-geneTable[geneTable$present > 3,]
      #geneTable$trait<-snps$trait[1]
      #geneTable$permutation<-snps$permutation[1]
      geneTable[,permutation := snps$permutation[1]]
      return(geneTable)
    }
  PermutationMerge$trait<-e$trait[1]
  write.csv(PermutationMerge, file = paste0(output,"permutation_files/",e$trait[1],"_PermutationMerge.csv"), row.names = F)
  ###calculate genes/loci stuff here
  
  genecountSummary<-foreach(snps=byPermutation_split, .packages = packages.loaded(), .combine = rbind) %dopar% {
    uniqueHits<-distinct(subset(snps, select=c(permutation,org,loci,trait)))
    MergeSubset<-subset(PermutationMerge, permutation==snps$permutation[1])
    if(nrow(MergeSubset)==0){
      print(paste0("No snp hits for permutation ", snps$permutation[1]," in ", snps$trait[1]))
      genecount<-NULL
      for(org in as.character(unique(uniqueHits$org))){
        genecount<-rbind(genecount,data.table(org=org,trait=uniqueHits$trait[1],dataset="permutation",GeneCount=0,
                                             LociCount=leafcollapsedSNPs[leafcollapsedSNPs$trait==uniqueHits$trait[1],paste(org)],check=snps$permutation[1]))
      }
      return(genecount)
    }else{
      for(s in 1:nrow(uniqueHits)){
        checking<-distinct(subset(MergeSubset, select = c(paste0(uniqueHits$org[s],"ID"),paste0(uniqueHits$org[s],"_loci"))))
        uniqueHits[s, genes := length(which(checking[,2]==uniqueHits$loci[s]))]
      }
      genecount<-NULL
      for(org in as.character(unique(uniqueHits$org))){
        summary<-as.data.table(table(uniqueHits$genes[uniqueHits$org==org]))
        summary[V1==0,N := N[V1==0]+leafcollapsedSNPs[leafcollapsedSNPs$trait==uniqueHits$trait[1],paste(org)]-sum(summary$N)]
        attachSummary<-data.table(org=org,trait=uniqueHits$trait[1],dataset="permutation",GeneCount=as.integer(as.character(summary$V1)),
                                  LociCount=summary$N, check=snps$permutation[1])
        genecount<-rbind(genecount,attachSummary)
      }
      return(genecount)
    }
  }
  permutation_genecounting<-rbind(permutation_genecounting,genecountSummary)
  ####end calculate genes stuff
  rm(PermutationMerge,byPermutation_split)
  }

print(paste("End permutations",Sys.time()))
rm(d)
print("Summarizing actual data")
Final_summary<-summaryTable(FinalMerge)   ##in old code called OrgTraitTableList
Final_summary$dataset<-"Actual"

Permutation_summary<-foreach(e=list.files(paste0(output,"permutation_files/"), pattern = "PermutationMerge.csv", full.names = T),
                             .packages = packages.loaded(), .combine = rbind) %dopar%{
  print(paste0("Summarizing ",e))
  f<-fread(e, sep = ",", stringsAsFactors = FALSE)
  if(!nrow(f)==0){
    fsummary<-summaryTable(f, random=T)
    return(fsummary)}
}

AllSummaries<-rbind(Final_summary,Permutation_summary)
AllSummaries$GeneCount<-as.numeric(AllSummaries$GeneCount)
AllSummaries$LociCount<-as.numeric(AllSummaries$LociCount)
###removing Cobalt for issues w/corn
#AllSummaries<-subset(AllSummaries, !trait=="Co")
write.csv(AllSummaries, file = paste0(output,"AllSummaries.csv"), row.names = F)

graphingDF<-ddply(AllSummaries, .(org, trait, dataset, present), summarize, GenesMean=mean(GeneCount), LociMean=mean(LociCount), 
                  Gene95=quantile(GeneCount, .95), Loci95=quantile(LociCount, .95), Gene99=quantile(GeneCount, .99), 
                  Loci99=quantile(LociCount, .99), Gene05=quantile(GeneCount, .05), Loci05=quantile(LociCount, 0.5),
                  Gene01=quantile(GeneCount, .01), Loci01=quantile(LociCount, .01))
graphingDF$Gene95[graphingDF$dataset=="Actual"]<-NA
graphingDF$Loci95[graphingDF$dataset=="Actual"]<-NA
graphingDF$Gene99[graphingDF$dataset=="Actual"]<-NA
graphingDF$Loci99[graphingDF$dataset=="Actual"]<-NA
graphingDF$Loci01[graphingDF$dataset=="Actual"]<-NA
graphingDF$Gene01[graphingDF$dataset=="Actual"]<-NA
graphingDF$Gene05[graphingDF$dataset=="Actual"]<-NA
graphingDF$Loci05[graphingDF$dataset=="Actual"]<-NA

graphing_split<-split(graphingDF, graphingDF$present)

backend<-foreach(data=graphing_split, .packages = packages.loaded()) %do% {
  
  #G95<-ggplot(data, aes(x=GenesMean, y=trait, group=dataset, color=dataset))+
    #geom_errorbarh(aes(xmax=Gene95, xmin=Gene05, height=0.75)) +
    #geom_point(aes(color=dataset, size=dataset)) + scale_size_manual(values=c(4,3)) +
    #labs(x = "Overlapped Genes Returned", title="Ortholgs in groups with 4/5 species representation") + 
    #theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5))+
    #scale_color_manual(values=c("#7fcdbb","#2c7fb8")) + facet_grid(~org, scales = "free")
  G99<-ggplot(data, aes(x=GenesMean, y=trait, group=dataset, color=dataset))+
    geom_errorbarh(aes(xmax=Gene99, xmin=Gene01, height=0.75)) +
    geom_point(aes(color=dataset, size=dataset)) + scale_size_manual(values=c(4,3)) +
    labs(x = "Overlapped Genes Returned") + theme_bw(base_size = 18) + 
    scale_color_manual(values=c("#7fcdbb","#2c7fb8"))# + facet_grid(~org, scales = "free")
  #L95<-ggplot(data, aes(x=LociMean, y=trait, group=dataset, color=dataset))+
    #geom_errorbarh(aes(xmax=Loci95, xmin=Loci05, height=0.75)) +
    #geom_point(aes(color=dataset, size=dataset)) + scale_size_manual(values=c(4,3)) +
    #labs(x = "Overlapped Genes Returned") + theme_bw(base_size = 18) + 
   # scale_color_manual(values=c("#7fcdbb","#2c7fb8")) + facet_grid(~org, scales = "free")
  #L99<-ggplot(data, aes(x=LociMean, y=trait, group=dataset, color=dataset))+
    #geom_errorbarh(aes(xmax=Loci99, xmin=Loci01, height=0.75)) +
    #geom_point(aes(color=dataset, size=dataset)) + scale_size_manual(values=c(4,3)) +
    #labs(x = "Overlapped Genes Returned") + theme_bw(base_size = 18) + 
    #scale_color_manual(values=c("#7fcdbb","#2c7fb8")) + facet_grid(~org, scales = "free")
  
 # ggsave(paste0(output,"graphs/A.thaliana_only_combo_",data$present[1],"_Genes95.png"),G95,height = 9,width=16)
  ggsave(paste0(output,"graphs/A.thaliana_only_99th_combo_",data$present[1],"_Genes99.png"),G99,height = 9,width=16)
  #ggsave(paste0(output,"graphs/A.thaliana_only_combo_",data$present[1],"_Loci95.png"),L95,height = 9,width=16)
 # ggsave(paste0(output,"graphs/A.thaliana_only_combo_",data$present[1],"_Loci99.png"),L99,height = 9,width=16)
  
}

priorityList<-foreach(o=as.character(metaTable$orgs), .packages = packages.loaded(), .combine = rbind) %do% {
  sub_FinalMerge<-FinalMerge[,grep(paste0(o,"ID|",o,"_loci|trait|present"),colnames(FinalMerge)),with=F]
  setnames(sub_FinalMerge, old=c(paste0(o,"ID"),paste0(o,"_loci")), new = c("ID","loci"))
  subList<-distinct(merge(sub_FinalMerge,snps_sub[,c("loci","genecount")],by = "loci"))
  subList$trait<-gsub("[0-9]*","",subList$trait)
  subList$gFDR<-(1-(1/(subList$genecount)))
  pFDR_backend<-foreach(row=iter(subList,by='row'), .packages = packages.loaded(), .combine = rbind) %do% {
    pFDRvalue<-(graphingDF$GenesMean[graphingDF$org==gsub("_.*","",row$loci) & graphingDF$trait==row$trait & graphingDF$present==row$present & 
                                 graphingDF$dataset=="permutations"] /
      graphingDF$GenesMean[graphingDF$org==gsub("_.*","",row$loci) & graphingDF$trait==row$trait & graphingDF$present==row$present & 
                             graphingDF$dataset=="Actual"])
    if(identical(pFDRvalue,numeric(0))){pFDRvalue<-0}
    return(data.frame(pFDR=pFDRvalue))
  }
  subList$pFDR<-pFDR_backend$pFDR
  return(subList)
}
priorityList<-priorityList[!priorityList$trait=="Co",]
###priorityList editing
#for(i in unique(Bdups$ID)){
#  subtab<-Bdups[Bdups$ID==i,]
#  if(subtab$pFDR[subtab$present==max(subtab$present)]>subtab$pFDR[subtab$present==min(subtab$present)]){print(subtab)}
#  
#}

##############################

priorityList$rank <- ((1 - priorityList$gFDR) * (1 - priorityList$pFDR) * (priorityList$present / nrow(metaTable)))
priorityList$rank_new <- ((1 - (priorityList$gFDR*0.2)) * (1 - (priorityList$pFDR)) * (priorityList$present / nrow(metaTable)))

### Dividing the priority list up by organism
#for(o in as.character(metaTable$orgs)){
#  individualPriority<-subset(priorityList, gsub("_.*","",priorityList$loci)==o)
#  individualPriority<-arrange(individualPriority, desc(rank))
#  write.csv(individualPriority, file = paste0(output,o,"priorityList.csv"), row.names = FALSE)
#}

### Dividing the priority list up by trait
priorityList_split<-split(priorityList, priorityList$trait)
for(t in priorityList_split){
  individualPriority<-arrange(t, desc(rank_new))
  uniqueIndividualPriority<-foreach(i=as.character(unique(individualPriority$ID)), .packages = packages.loaded(), .combine = rbind) %dopar% {
    uniqueSub<-individualPriority[individualPriority$ID==i,]
    if(nrow(uniqueSub)>1){
      return(uniqueSub[which.max(uniqueSub$present),])
    }else{
      return(uniqueSub)}
    }
  write.csv(uniqueIndividualPriority, file = paste0(output,"priority_lists/",t$trait[1],"_priorityList.csv"), row.names = FALSE)
}