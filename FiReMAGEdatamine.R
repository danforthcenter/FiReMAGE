###retrieving data from CGAS trials
library(readr)
library(plyr)
library(dplyr)
library(data.table)
library(doParallel)
library(reader)
library(car)
library(ggplot2)

inputdir<-"./data/FiReMAGEtesting/FM_new/"

A.thaliana_known <- read_csv("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs_Liberal/A.thaliana_knownIonomicsGenesWOrthologs_Liberal.csv")
G.max_known <- read_csv("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs_Liberal/G.max_knownIonomicsGenesWOrthologs_Liberal.csv")
S.bicolor_known <- read_csv("./knownIonomicsGenesWOrthologs_Liberal/S.bicolor_knownIonomicsGenesWOrthologs_Liberal.csv")
Z.mays_known <- read_csv("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs_Liberal/Z.mays_knownIonomicsGenesWOrthologs_Liberal.csv")
O.sativa_known <- read_csv("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs_Liberal/O.sativa_knownIonomicsGenesWOrthologs_Liberal.csv")

FinalMerge<-read.csv(file = paste0(inputdir,"FinalMerge.csv"),sep = ",",header = T, stringsAsFactors = F)
snp_gene_hitTable<-read.csv(file = paste0(inputdir,"snp_gene_hitTable"),sep = ",",header = T, stringsAsFactors = F)
orgs<-c("A.thaliana","G.max","O.sativa","S.bicolor","Z.mays")

knownList<-foreach(i=orgs, .packages = packages.loaded(), .combine = rbind) %do% {
  orgSet<-get(paste0(i,"_known"))
  setnames(orgSet, old=c(colnames(orgSet)[grep("ortholog",colnames(A.thaliana_known),invert=T)]),
           new=c(colnames(A.thaliana_known)[grep("ortholog",colnames(A.thaliana_known),invert = T)]))
  orgSet<-orgSet[,colnames(orgSet)[grep("ortholog",colnames(orgSet), invert = TRUE)]]
  return(orgSet)
}
knownList<-knownList[,-c(7,9)]
knownList$Elements<-gsub(" ","",knownList$Elements)

priorityList_files<-list.files(paste0(inputdir,"/priority_lists/"), full.names = TRUE)
priorityList<-foreach(i=priorityList_files, .packages = packages.loaded(), .combine = rbind) %do% {
  l<-read.csv(i, sep = ",", header = T)
  return(l)
}

adjacent_elements<-vector(mode="list",length = 17)
names(adjacent_elements)<-c("B","Ca","Cd","Co","Cu","Fe","K","Mg","Mn","Mo","Na","Ni","P","Rb","S","Se","Zn")
adjacent_elements[["B"]]<-c("B")
adjacent_elements[["Ca"]]<-c("Ca","Sr")
adjacent_elements[["Cd"]]<-c("Cd","Zn")
adjacent_elements[["Co"]]<-c("Co")
adjacent_elements[["Cu"]]<-c("Cu","Zn")
adjacent_elements[["Fe"]]<-c("Fe","Mn","Zn")
adjacent_elements[["K"]]<-c("K","Rb","Na")
adjacent_elements[["Mg"]]<-c("Mg")
adjacent_elements[["Mn"]]<-c("Mn","Fe")
adjacent_elements[["Mo"]]<-c("Mo")
adjacent_elements[["Na"]]<-c("Na","K","Rb")
adjacent_elements[["Ni"]]<-c("Ni")
adjacent_elements[["P"]]<-c("P")
adjacent_elements[["Rb"]]<-c("Rb","K")
adjacent_elements[["S"]]<-c("S")
adjacent_elements[["Se"]]<-c("Se","S")
adjacent_elements[["Sr"]]<-c("Ca")
adjacent_elements[["Zn"]]<-c("Zn","Cd","Fe","Cu")
###not element specific
overlaps<-merge(knownList,priorityList,by.x = "GeneID",by.y = "ID")
#write.csv(overlaps, file = "./FiReMAGE_KIG_overlapped_genes.csv", row.names = F)
lengthTable<-foreach(i=orgs, .packages = packages.loaded(), .combine = rbind) %do% {
  ll<-subset(overlaps, Species==i)
  l<-length(unique(ll$GeneID))
  return(data.frame(org=i,numGenes=l,totalGenes=nrow(knownList[knownList$Species==i,])))
}
lengthTable$fraction<-lengthTable$numGenes/lengthTable$totalGenes

ggplot(lengthTable, aes(x=org, y=fraction)) +
  geom_bar(stat = 'identity', position = 'stack') + 
  xlab("") + ylab("Fraction of KIG Genes") + theme_bw(base_size = 20)


###element specific
overlaps_trait_specific<-foreach(t=traits, .packages = packages.loaded(), .combine = rbind) %do% {
  known_sub<-knownList[grep(paste0(t),knownList$Elements),]
  priorityList_sub<-priorityList[priorityList$trait==t,]
  overlap<-merge(known_sub,priorityList_sub,by.x = "GeneID", by.y = "ID",)
  return(overlap)
}
#write.csv(overlaps_trait_specific, file = "./FiReMAGE_KIG_overlapped_genes_byElem.csv", row.names = F)
lengthTable_specific<-foreach(i=orgs, .packages = packages.loaded(), .combine = rbind) %do% {
  ll<-subset(overlaps_trait_specific, Species==i)
  l<-length(unique(ll$GeneID))
  return(data.frame(org=i,numGenes=l,totalGenes=nrow(knownList[knownList$Species==i,])))
}
lengthTable_specific$fraction<-lengthTable_specific$numGenes/lengthTable_specific$totalGenes

ggplot(lengthTable_specific, aes(x=org, y=fraction)) +
  geom_bar(stat = 'identity', position = 'stack') + 
  xlab("") + ylab("Fraction of KIG Genes") + theme_bw(base_size = 20)

lengthTable_specific$diff<- c(0.023,0.045,0.01,0.012,0.025)
####checking snp distributions over elements in our datasets
snpfiles<-list.files(path = "./data/el_snps/current/", full.names = T,recursive = F,pattern = "snps")
distr<-foreach(f=snpfiles, .packages = packages.loaded(), .combine=rbind) %do% {
  snpFile<-read.csv(f)
  snpdist<-as.data.frame(table(snpFile$trait))
  colnames(snpdist)<-c("trait","snps")
  snpdist$org<-snpFile$org[1]
  return(snpdist)
}
ggplot(distr[!(distr$org=="O.sativa"|distr$trait=="Co59"),], aes(x=trait, y=snps, fill=org))+facet_grid(~ trait, space = "free_x",scales = "free_x")+
  geom_bar(stat='identity', position='dodge')


###checking KIG conservation in all the species

testingSet<-overlaps_trait_specific[overlaps_trait_specific$Species=="A.thaliana" & !overlaps_trait_specific$trait=="Co" & overlaps_trait_specific$genecount > 1,]

outee<-foreach(l=as.character(testingSet$loci), .packages=packages.loaded()) %do% {
  underlyingGenes<-FinalMerge[which(FinalMerge$A.thaliana_loci==l),]
  length(unique(underlyingGenes$A.thalianaID))
  orthologousGenes<-foreach(g=as.character(unique(underlyingGenes$A.thalianaID)), .packages= packages.loaded(), .combine = rbind) %do% {
    orthos<-unlist(sapply(underlyingGenes[underlyingGenes$A.thalianaID==g,colnames(underlyingGenes)[grep("ID",colnames(underlyingGenes))]],
                                    FUN=unique), use.names = F)
    geneCoords<-all_gene_coords[which(all_gene_coords$`Gene Name` %in% orthos),]
    geneCoords$A.thalianaID<-g
    geneCoords<-geneCoords[-c(geneCoords$`Gene Name`==g),]
    return(geneCoords)
  }
  regions<-foreach(o=as.character(unique(geneCoords$org)), .packages = packages.loaded(), .combine = rbind) %do% {
    orgsub<-subset(orthologousGenes, org==o)
    orgsub$groups<-"a"
    orgsub<-orgsub[order(`Chromosome Name`,`Gene Start (bp)`),]
    for(i in 1:(nrow(orgsub)-1)){
      if(orgsub$`Gene Start (bp)`[i+1] - orgsub$`Gene End (bp)`[i] < metaTable$range[metaTable$orgs==o] &
         orgsub$`Chromosome Name`[i]==orgsub$`Chromosome Name`[i+1]){
        print("In range")
        orgsub$groups[i+1]<-orgsub$groups[i]
      }else{
        print("Next group")
        orgsub$groups[i+1]<-letters[which(letters==orgsub$groups[i])+1]
      }
    }
    for(group in as.character(unique(orgsub$groups))){
      if(length(which(orgsub$groups==group)) > 1 & abs(min(orgsub$`Gene End (bp)`[orgsub$groups==group]) - 
            max(orgsub$`Gene Start (bp)`[orgsub$groups==group])) > as.integer(2*1e+06)){
        print(paste0(group," is larger than 2x it's range: ",metaTable$range[metaTable$orgs==o]))
      }
    }
    return(orgsub)
  }
  return(regions)
}
names(outee)<-as.character(testingSet$loci)

    
####look at how many KIG genes have orthologs in all 4/5 species
lengthTable<-foreach(o=as.character(metaTable$orgs), .packages = packages.loaded(), .combine = rbind) %do%{
  known<-get(paste0(o,"_known"))
  known_subset<-known[,c(1,2,4,9,which(colnames(known) 
                                       %in% c("O.sativa orthologs","Z.mays orthologs","G.max orthologs",
                                              "S.bicolor orthologs","A.thaliana orthologs")))]
  complete_known_subset<-known_subset[complete.cases(known_subset),]
  not_complete_known_subset<-known_subset[!complete.cases(known_subset),]
  df<-data.frame(org=o,p=nrow(complete_known_subset)/nrow(known),
                 coverlap=length(which(complete_known_subset$GeneID %in% priorityList$ID)), 
                 complete_cases=length(which(complete_known_subset$GeneID %in% priorityList$ID))/nrow(known),
                 ncoverlap=length(which(not_complete_known_subset$GeneID %in% priorityList$ID)),
                 incomplete_cases=length(which(not_complete_known_subset$GeneID %in% priorityList$ID))/nrow(known))
  
}
lengthmelt<-melt(lengthTable[,c("org","complete_cases","incomplete_cases")])
ggplot(lengthmelt, aes(x=org, y=value, fill=variable)) +
  geom_bar(stat = 'identity', position = position_dodge()) + 
  scale_fill_manual(values = c("cadetblue","seagreen"))+
  xlab("") + ylab("Fraction of KIG Genes") + theme_bw(base_size = 20)

ggplot(lengthTable, aes(x=org, y=overlap_fraction)) + geom_bar(stat='identity', position='stack')+
  xlab("") + ylab("Fraction of KIG Genes") + theme_bw(base_size = 20)

###number of genes under snp peaks slide
###get the single snps dataset returned right after the snp_gene_hitTable is made
snp_gene_hitTable<-fread(paste0(inputdir,"snp_gene_hitTable.csv"), sep = ",", stringsAsFactors = F)
snps_sub<-fread(paste0(inputdir, "snps_sub.csv"), sep = ",", stringsAsFactors = F)
perm_snps<-fread(paste0(inputdir, "AllPermuts.csv"), sep = ",", stringsAsFactors = F)
permutation_gene_hitTable<-fread(paste0(inputdir, "permutation_gene_hitTable.csv"), sep = ",",stringsAsFactors = F)


#genecount for single dataset actual snps

genecount<-foreach(s=1:nrow(snps_sub), .packages = packages.loaded()) %dopar% {
  lociHit<-subset(snp_gene_hitTable, loci==snps_sub$loci[s])
  if(nrow(lociHit)==0){
    return(0)
  }else{
    return(nrow(lociHit))
  }
}
single_snp_sub$genecount<-unlist(genecount)

#genecount for single dataset permutation snps
genecount<-foreach(s=1:nrow(perm_snps), .packages = packages.loaded()) %dopar% {
  lociHit<-subset(permutation_gene_hitTable, loci==perm_snps$loci[s])
  if(nrow(lociHit)==0){
    return(0)
  }else{
    return(nrow(lociHit))
  }
}
single_perm_snps$genecount<-unlist(genecount)

KIGcount<-foreach(s=1:nrow(single_snps_sub), .packages = packages.loaded()) %dopar% {
  lociHit<-subset(snp_gene_hitTable, loci==single_snps_sub$loci[s])
  if(nrow(lociHit)==0){
    return(0)
  }else{
    return(length(which(lociHit$`Gene Name` %in% knownList$GeneID)))
  }
}
single_snps_sub$KIGcount<-unlist(KIGcount)

#for the usual snps_sub output from FiReMAGE.R script

KIGcount<-foreach(s=1:nrow(snps_sub), .packages = packages.loaded()) %dopar% {
  lociHit<-subset(snp_gene_hitTable, loci==snps_sub$loci[s])
  if(nrow(lociHit)==0){
    return(0)
  }else{
    checking<-FinalMerge[,c(paste0(snps_sub$org[s],"ID"),paste0(snps_sub$org[s],"_loci")),with=FALSE]
    colnames(checking)<-c("ID","loci")
    checking<-subset(checking, loci==snps_sub$loci[s])
    return(length(which(lociHit$`Gene Name`[which(lociHit$`Gene Name` %in% checking$ID)] %in% knownList$GeneID)))
  }
}
snps_sub$KIGcount<-unlist(KIGcount)
snps_sub$ratio<- (snps_sub$KIGcount / snps_sub$genecount)
single_snps_sub$ratio <- (single_snps_sub$KIGcount / single_snps_sub$genecount)
snps_sub[is.na(snps_sub)]<-0
single_snps_sub[is.na(single_snps_sub)]<-0
mean(single_snps_sub$ratio)
mean(snps_sub$ratio)
ttest<-t.test(single_snps_sub$ratio, snps_sub$ratio)
ttest$statistic

nozeros_single<-single_snps_sub[single_snps_sub$genecount > 0,]
nozeros_snp_sub<-snps_sub[snps_sub$genecount > 0,]
t.test(nozeros_single$ratio, nozeros_snp_sub$ratio)
#percentage of snps not covered by KIG genes
length(which(nozeros_snp_sub$KIGcount == 0))/nrow(nozeros_snp_sub)
length(which(nozeros_single$KIGcount==0))/nrow(nozeros_snp_sub)

###KIG under GWAS peaks
#specific to trait

snp_hitting_specific_KIG<-foreach(e=traits, .packages = packages.loaded(), .combine = rbind) %do% {
  known_sub<-knownList[grep(paste0("^",e,"$|^",e,",|,",e,",|,",e,"$"),knownList$Elements),]
  snp_hit_sub<-subset(snp_gene_hitTable, trait==e)
  snp_specific_KIG<-known_sub[which(known_sub$GeneID %in% snp_hit_sub$`Gene Name`),]
  graphing<-as.data.frame(table(snp_specific_KIG$Species))
  if(nrow(graphing)==0){
    graphing<-data.frame(Var1=c("A.thaliana","G.max","O.sativa","S.bicolor","Z.mays"),`GWAS overlap`=c(0,0,0,0,0))
  }
  colnames(graphing)[2]<-"GWAS overlap"
  graphing<-merge(graphing,as.data.frame(table(known_sub$Species)))
  colnames(graphing)[c(1,3)]<-c("Species","Total")
  graphing$trait<-e
  return(graphing)
}
aggdata<-aggregate(snp_hitting_specific_KIG$`GWAS overlap`, by=list(snp_hitting_specific_KIG$Species), FUN =sum)
#test<-aggregate(snp_hitting_specific_KIG$Total, by=list(snp_hitting_specific_KIG$Species), FUN=sum)
#aggdata<-cbind(aggdata,test$x)
aggdata$total<-c(136,268,141,135,152)
colnames(aggdata)<-c("Species","snphit","total")
aggdata$proportion<-aggdata$snphit/aggdata$total


perm_hitting_specific_KIG<-foreach(e=traits, .packages = packages.loaded(), .combine = rbind) %:%
  foreach(n=1:100, .packages = packages.loaded(), .combine = rbind) %dopar% {
    
  known_sub<-knownList[grep(paste0("^",e,"$|^",e,",|,",e,",|,",e,"$"),knownList$Elements),]
  snp_hit_sub<-subset(permutation_gene_hitTable, trait==e & permutation==n)
  snp_specific_KIG<-known_sub[which(known_sub$GeneID %in% snp_hit_sub$`Gene Name`),]
  graphing<-as.data.frame(table(snp_specific_KIG$Species))
  if(nrow(graphing)==0){
    graphing<-data.frame(Var1=c("A.thaliana","G.max","O.sativa","S.bicolor","Z.mays"),`GWAS overlap`=c(0,0,0,0,0))
  }
  colnames(graphing)[2]<-"GWAS overlap"
  graphing<-merge(graphing,as.data.frame(table(known_sub$Species)))
  colnames(graphing)[c(1,3)]<-c("Species","Total")
  graphing$trait<-e
  graphing$permutation<-n
  return(graphing)
}
permaggdata<-aggregate(perm_hitting_specific_KIG$`GWAS overlap`, by=list(perm_hitting_specific_KIG$Species,perm_hitting_specific_KIG$permutation), FUN =sum)
meanpermaggdata<-aggregate(permaggdata$x, by=list(permaggdata$Group.1), FUN=mean)
quant95<-aggregate(permaggdata$x, by=list(permaggdata$Group.1), FUN=quantile, .95)
quant5<-aggregate(permaggdata$x, by=list(permaggdata$Group.1), FUN=quantile, .05)
#test<-aggregate(perm_hitting_specific_KIG$Total, by=list(perm_hitting_specific_KIG$Species), FUN=sum)
#aggdata<-cbind(aggdata,test$x)


meanpermaggdata<-cbind(meanpermaggdata,quant95$x,quant5$x)
meanpermaggdata$total<-c(136,268,141,135,152)
colnames(meanpermaggdata)<-c("Species","snphit","95%","5%","total")
meanpermaggdata$proportion<-meanpermaggdata$snphit/meanpermaggdata$total
meanpermaggdata$proportion95<-meanpermaggdata$`95%`/meanpermaggdata$total
meanpermaggdata$proportion5<-meanpermaggdata$`5%`/meanpermaggdata$total
meltdata<-melt(meanpermaggdata)
aggdata<-cbind(aggdata,data.frame(`95%`=c(0,0,0,0,0),`5%`=c(0,0,0,0,0),proportion95=c(0,0,0,0,0),proportion5=c(0,0,0,0,0)))
colnames(aggdata)[c(5,6)]<-c("95%","5%")
aggdata$dataset<-"Actual"
meanpermaggdata$dataset<-"permutation"
aggdata<-rbind(aggdata,meanpermaggdata)
ggplot(aggdata, aes(x=Species, y=proportion, group=dataset, fill=dataset))+
  geom_bar(aes(x=Species, y=proportion),stat = "identity",position=position_dodge())+
  geom_errorbar(aes(x=Species, ymin=proportion5, ymax=proportion95), position = position_dodge())+
  xlab("")+ylab("Proportion")+theme_bw(base_size = 18)+ggtitle("KIG list genes under GWAS peaks")+
  scale_fill_manual(values=c("#7fcdbb","#2c7fb8"))
#####
##the reverse from above. Proportion of snps not covered by KIG list?
#specific to trait
inputdir<-"./data/FiReMAGEtesting/FM_4/"
snp_gene_hitTable<-fread(paste0(inputdir,"snp_gene_hitTable.csv"), sep = ",", stringsAsFactors = F)
snps_sub<-fread(paste0(inputdir, "snps_sub.csv"), sep = ",", stringsAsFactors = F)
perm_snps<-fread(paste0(inputdir, "AllPermuts.csv"), sep = ",", stringsAsFactors = F)
permutation_gene_hitTable<-fread(paste0(inputdir, "permutation_gene_hitTable.csv"), sep = ",",stringsAsFactors = F)

snp_hitting_specific_KIG<-foreach(e=traits, .packages = packages.loaded(), .combine = rbind) %do% {
  known_sub<-knownList[grep(paste0("^",e,"$|^",e,",|,",e,",|,",e,"$"),knownList$Elements),]
  snp_hit_sub<-subset(snp_gene_hitTable, trait==e)
  snp_specific_KIG<-snp_hit_sub[which(snp_hit_sub$`Gene Name` %in% known_sub$GeneID),]
  unique_loci<-distinct(subset(snp_specific_KIG, select = c("org","loci")))
  graphing<-as.data.frame(table(unique_loci$org))
  if(nrow(graphing)==0){
    graphing<-data.frame(Var1=c("A.thaliana","G.max","O.sativa","S.bicolor","Z.mays"),`GWAS overlap`=c(0,0,0,0,0))
  }
  colnames(graphing)[2]<-"GWAS overlap"
  graphing$trait<-e
  return(graphing)
}
aggdata<-aggregate(snp_hitting_specific_KIG$`GWAS overlap`, by=list(snp_hitting_specific_KIG$Var1), FUN =sum)
total<-as.data.frame(table(snps_sub$org))
colnames(total)[2]<-"total"
aggdata<-merge(aggdata,total, by.x = "Group.1", by.y = "Var1")
colnames(aggdata)<-c("Species","snphit","total")
aggdata$proportion<-aggdata$snphit/aggdata$total


perm_hitting_specific_KIG<-foreach(e=traits, .packages = packages.loaded(), .combine = rbind) %:%
  foreach(n=1:100, .packages = packages.loaded(), .combine = rbind) %dopar% {
    
    known_sub<-knownList[grep(paste0("^",e,"$|^",e,",|,",e,",|,",e,"$"),knownList$Elements),]
    snp_hit_sub<-subset(permutation_gene_hitTable, trait==e & permutation==n)
    snp_specific_KIG<-snp_hit_sub[which(snp_hit_sub$`Gene Name` %in% known_sub$GeneID),]
    unique_loci<-distinct(subset(snp_specific_KIG, select = c("org","loci")))
    graphing<-as.data.frame(table(unique_loci$org))
    if(nrow(graphing)==0){
      graphing<-data.frame(Var1=c("A.thaliana","G.max","O.sativa","S.bicolor","Z.mays"),`GWAS overlap`=c(0,0,0,0,0))
    }
    colnames(graphing)[2]<-"GWAS overlap"
    graphing$trait<-e
    graphing$permutation<-n
    return(graphing)
  }
permaggdata<-aggregate(perm_hitting_specific_KIG$`GWAS overlap`, by=list(perm_hitting_specific_KIG$Var1,perm_hitting_specific_KIG$permutation), FUN =sum)
meanpermaggdata<-aggregate(permaggdata$x, by=list(permaggdata$Group.1), FUN=mean)
quant95<-aggregate(permaggdata$x, by=list(permaggdata$Group.1), FUN=quantile, .95)
quant5<-aggregate(permaggdata$x, by=list(permaggdata$Group.1), FUN=quantile, .05)
#test<-aggregate(perm_hitting_specific_KIG$Total, by=list(perm_hitting_specific_KIG$Species), FUN=sum)
#aggdata<-cbind(aggdata,test$x)
meanpermaggdata<-cbind(meanpermaggdata,quant95$x,quant5$x)
meanpermaggdata<-merge(meanpermaggdata,total, by.x = "Group.1", by.y = "Var1")
colnames(meanpermaggdata)<-c("Species","snphit","95%","5%","total")

meanpermaggdata$proportion<-meanpermaggdata$snphit/meanpermaggdata$total
meanpermaggdata$proportion95<-meanpermaggdata$`95%`/meanpermaggdata$total
meanpermaggdata$proportion5<-meanpermaggdata$`5%`/meanpermaggdata$total
aggdata<-cbind(aggdata,data.frame(`95%`=0,`5%`=0,proportion95=0,proportion5=0))
colnames(aggdata)[c(5,6)]<-c("95%","5%")
aggdata$dataset<-"Actual"
meanpermaggdata$dataset<-"permutation"
aggdata<-rbind(aggdata,meanpermaggdata)
ggplot(aggdata, aes(x=Species, y=proportion, group=dataset, fill=dataset))+
  geom_bar(aes(x=Species, y=proportion),stat = "identity",position=position_dodge())+
  geom_errorbar(aes(x=Species, ymin=proportion5, ymax=proportion95), position = position_dodge())+
  scale_y_continuous(labels = percent)+
  xlab("")+ylab("GWAS intervals w/ KIG list hit")+theme_bw(base_size = 18)+ggtitle("")+
  scale_fill_manual(values=c("#7fcdbb","#2c7fb8"))


#non-specific to trait
snp_hitting_KIG<-knownList[which(knownList$GeneID %in% snp_gene_hitTable$`Gene Name`),]
graphing<-as.data.frame(table(snp_hitting_KIG$Species))
colnames(graphing)[2]<-"GWAS overlap"
graphing<-merge(graphing,as.data.frame(table(knownList$Species)))
colnames(graphing)[c(1,3)]<-c("Species","Total")
graphing$proportion<-graphing$`GWAS overlap`/graphing$Total
graphing
ggplot(graphing, aes(x=Species, y=proportion))+
  geom_bar(stat = "identity",fill="steelblue")+
  xlab("")+ylab("Proportion")+theme_bw(base_size = 18)+ggtitle("KIG list genes under GWAS peaks")

permutation_hitting_KIG<-knownList[which(knownList$GeneID %in% permutation_gene_hitTable$`Gene Name`),]
graphing<-as.data.frame(table(permutation_hitting_KIG$Species))
