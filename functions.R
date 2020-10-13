GeneFilePrep<-function(filePath){
  Gene1names<-c("Gene1id","Gene1org")
  Gene2names<-c("Gene2id","Gene2org")
  orgList<-foreach(org=1:(nrow(metaTable)), .packages=packages.loaded()) %:%
    foreach(i=which(!c(1:nrow(metaTable))==org)) %do% {
      print(paste0(metaTable$orgs[org],"_",metaTable$orgs[i]))
      file<-list.files(path = filePath, pattern = paste0(metaTable$orgs[org],"_",metaTable$orgs[i],"|",metaTable$orgs[i],"_",metaTable$orgs[org]), 
                       full.names = TRUE)[1]
      comparison<-fread(file, sep = ",", stringsAsFactors = FALSE, 
                        select = c("Homolog > Gene . Primary Identifier","Homolog > Gene > Organism Name",
                                   "Homolog > Ortholog _ Gene . Primary Identifier","Homolog > Ortholog _ Gene > Organism Name"))
      comparison[, `Homolog > Gene > Organism Name` := gsub(" ","",`Homolog > Gene > Organism Name`),]
      comparison[, `Homolog > Ortholog _ Gene > Organism Name` := gsub(" ","",`Homolog > Ortholog _ Gene > Organism Name`),]
      setnames(comparison,old = c("Homolog > Gene . Primary Identifier","Homolog > Gene > Organism Name",
                                  "Homolog > Ortholog _ Gene . Primary Identifier","Homolog > Ortholog _ Gene > Organism Name"),
               new = c(paste0(comparison$`Homolog > Gene > Organism Name`[1],"ID"),paste0(comparison$`Homolog > Gene > Organism Name`[1],"Org"),
                       paste0(comparison$`Homolog > Ortholog _ Gene > Organism Name`[1],"ID"),
                       paste0(comparison$`Homolog > Ortholog _ Gene > Organism Name`[1],"Org")))
      
      return(comparison)
    }
  merging<-foreach(l=orgList, .packages = packages.loaded()) %do% {
    test<-Reduce(function(x,y) merge(x=x, y=y, allow.cartesian=TRUE, all=TRUE), l)
    return(test)
  }
  IntersectsInferring<-distinct(bind_rows(merging))
  IntersectsInferring<-IntersectsInferring[rowSums(is.na(IntersectsInferring[,grep("ID",colnames(IntersectsInferring)),with=FALSE])) < 3,]
  print("Finished ortholog table merging")
  return(IntersectsInferring)
}
snpFilePrep<-function(org){
  file<-list.files(path = paste0(snpPath), pattern = paste0(org), full.names = TRUE)
  #SNPfile<-fread(file, sep = ",", stringsAsFactors = TRUE)
  SNPfile<-read.csv(file)
  ###ordering--must have chromosomes and SNPs in numerical order
  SNPfile$chr<-as.integer(SNPfile$chr)
  SNPfile<-SNPfile[order(SNPfile[,"trait"], SNPfile[,"org"], SNPfile[,"chr"], SNPfile[,"bp"]),]
  return(SNPfile)
}
collapsing<-function(snpTable,counter=1,done=F,permutations=F){
  if(permutations){
    for (j in 1:nrow(snpTable)){
      current_pos<-snpTable[j,"bp"]
      next_pos<-snpTable[(j+1),"bp"]
      next_range<-snpTable[(j+1),"range"]
      this_org<-snpTable[j,"org"]
      next_org<-snpTable[(j+1),"org"]
      this_chr<-snpTable[j,"chr"]
      next_chr<-snpTable[(j+1),"chr"]
      this_trait<-as.character(snpTable[j,"trait"])
      next_trait<-as.character(snpTable[(j+1),"trait"])
      
      if(j >= nrow(snpTable)){
        print("Bam! End of randomizing ranges.")
        done <- T
        break
        
      }else if (this_org != next_org){
        #print("org switch")
        #print(j)
        
      }else if (this_org != next_org || this_trait != next_trait || this_chr != next_chr){#for org/chr/el changes
        #print("change")
        #print(j)
        
      }else if ((as.numeric(next_pos)-as.numeric(current_pos)) < as.numeric(next_range)){
        snpTable[(j+1),"new_ranges"]<-as.numeric(next_pos)-as.numeric(current_pos) 
        
      }else if ((as.numeric(next_pos)-as.numeric(current_pos)) > as.numeric(next_range)){
        snpTable[(j+1),"new_ranges"]<-as.numeric(next_range)
        
      }
      if (done==T){
        break
      }
    }
  }else{
    repeat{
      start_pos <-snpTable[(counter),"bp"]
      for (j in 1:nrow(snpTable)){
        current_pos<-snpTable[j, "bp"]
        next_pos <- snpTable[(j+1),"bp"]
        this_org<-as.character(snpTable[j,"org"])
        next_org<-as.character(snpTable[(j+1),"org"])
        this_trait<-as.character(snpTable[j,"trait"])
        next_trait<-as.character(snpTable[(j+1),"trait"])
        this_chr<-snpTable[j,"chr"]
        next_chr<-snpTable[(j+1),"chr"]
        if(j >= nrow(snpTable)){
          print("Bam! End of Search")
          done <- TRUE
          snpTable[(counter), "clpsRanges"]<-as.numeric(current_pos)-snpTable[(counter), "bp"]
          snpTable[(counter), "start"] <- snpTable[(counter), "bp"]-snpTable[(counter), "range"]
          snpTable[(counter), "end"] <- snpTable[(counter), "bp"]+snpTable[(counter), "clpsRanges"]+snpTable[(counter), "range"]
          break
        }else if (this_org != next_org || this_trait != next_trait || this_chr != next_chr){#for chr/org/el changes
          snpTable[(counter), "SNPend"]<-as.numeric(current_pos)
          snpTable[(counter), "clpsRanges"]<-as.numeric(current_pos)-snpTable[(counter), "bp"]
          snpTable[(counter), "start"] <- snpTable[(counter), "bp"]-snpTable[(counter), "range"]
          snpTable[(counter), "end"] <- snpTable[(counter), "bp"]+snpTable[(counter), "clpsRanges"]+snpTable[(counter), "range"] 
          counter<-j+1
        }else if ((as.numeric(next_pos)-as.numeric(current_pos)) < 2*snpTable[(counter), "range"]){
          rows_to_remove<<-c(rows_to_remove, (j+1)) 
        }else if ((as.numeric(next_pos)-as.numeric(current_pos)) > 2*snpTable[(counter), "range"]){
          snpTable[(counter), "SNPend"]<-as.numeric(current_pos) 
          snpTable[(counter), "clpsRanges"]<-as.numeric(current_pos)-snpTable[(counter), "bp"]
          snpTable[(counter), "start"] <- snpTable[(counter), "bp"]-snpTable[(counter), "range"]
          snpTable[(counter), "end"] <- snpTable[(counter), "bp"]+snpTable[(counter), "clpsRanges"]+snpTable[(counter), "range"] 
          counter<- j+1
        }  
      }
      if(done == T){
        break
      }
    }
  }
  return(snpTable)
}

getOverlapIndxs <- function(org1,org2,snpTable,orthoTable){
  if(within){
    org1Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org1,],by.x=c("Gene1chr","Gene1start","Gene1end"),by.y=c("chr","start","end"),
                                         mult="first",type="within",which=TRUE))))
    org2Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org2,],by.x=c("Gene2chr","Gene2start","Gene2end"),by.y=c("chr","start","end"),
                                         mult="first",type="within",which=TRUE))))
  }else{
    org1Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org1,],by.x=c("Gene1chr","Gene1start","Gene1end"),by.y=c("chr","start","end"),
                                         mult="first",type="any",which=TRUE))))
    org2Indxs <- which(!(is.na(foverlaps(orthoTable,snpTable[org==org2,],by.x=c("Gene2chr","Gene2start","Gene2end"),by.y=c("chr","start","end"),
                                         mult="first",type="any",which=TRUE))))
  }
    intersect(org1Indxs,org2Indxs)
}
orthologousOverlapGenes<-function(snpTable, orthoTableReturn, random=F){
  setDT(snpTable)
  setkey(snpTable,chr,start,end)
  if(random==F){
    for(o in as.character(unique(snpTable$org))){
      #print("...")
      orglist<-c(paste0(o,"GeneChr"),paste0(o,"GeneStart"),paste0(o,"GeneEnd"))
      orthoTableReturn<-foverlaps(orthoTableReturn,subset(snpTable, org==o),by.x = c(orglist), by.y = c("chr","start","end"), type = "any")
      orthoTableReturn$loci[is.na(orthoTableReturn$loci)]<-""
      orthoTableReturn<- orthoTableReturn %>% unite(newloci, paste0(o,"_loci"), loci, sep = " ")
      orthoTableReturn$newloci<-gsub("\\s+",";", gsub("^\\s+|\\s+$", "",orthoTableReturn$newloci))
      orthoTableReturn[,paste0(o,"_loci")]<-orthoTableReturn$newloci
      orthoTableReturn<-subset(orthoTableReturn, select = -c(org,bp,trait,start,end,newloci))
    }
    #print(paste("Begin foverlaps: ",Sys.time()))
    #orthoTableReturn[orthoTableReturn==""]<-NA
    #orthoTableReturn<-replace_with_na_all(orthoTableReturn, condition = ~.x=="")
    #orthoTableReturn$trait<-snpTable$trait[1]
    #print(paste("End foverlaps: ",Sys.time()))
    #print(paste("Begin foverlaps 2: ",Sys.time()))
    #orthoTableReturn<-orthoTableReturn[rowSums(is.na(orthoTableReturn)) != length(grep("_loci", colnames(orthoTableReturn))),]
    orthoTableReturn$present<-(length(unique(snpTable$org))-rowSums(orthoTableReturn==""))
    orthoTableReturn<-orthoTableReturn[orthoTableReturn$present > 1,]
    #print(paste("End foverlaps 2: ",Sys.time()))
    
    return(orthoTableReturn)
  }else{
    for(o in as.character(unique(snpTable$org))){
      #print("...")
      orglist<-c(paste0(o,"GeneChr"),paste0(o,"GeneStart"),paste0(o,"GeneEnd"))
      orthoTableReturn<-foverlaps(orthoTableReturn,subset(snpTable, org==o),by.x = c(orglist), by.y = c("chr","start","end"), type = "any")
      orthoTableReturn$loci[is.na(orthoTableReturn$loci)]<-""
      orthoTableReturn<- orthoTableReturn %>% unite(newloci, paste0(o,"_loci"), loci, sep = " ")
      orthoTableReturn$newloci<-gsub("\\s+",";", gsub("^\\s+|\\s+$", "",orthoTableReturn$newloci))
      orthoTableReturn[,paste0(o,"_loci")]<-orthoTableReturn$newloci
      orthoTableReturn<-subset(orthoTableReturn, select = -c(permutation,trait,org,bp,start,end,newloci))
    }
    orthoTableReturn<-subset(orthoTableReturn, select=c(colnames(orthoTableReturn)[grep("ID|loci", colnames(orthoTableReturn))]))
    #orthoTableReturn[orthoTableReturn==""]<-NA
    #orthoTableReturn$trait<-snpTable$trait[1]
    #orthoTableReturn<-orthoTableReturn[rowSums(is.na(orthoTableReturn)) != length(grep("_loci", colnames(orthoTableReturn))),]
    orthoTableReturn$present<-(length(unique(snpTable$org))-rowSums(orthoTableReturn==""))
    orthoTableReturn<-orthoTableReturn[orthoTableReturn$present > 1,]
    snpTable<-snpTable[order(snpTable[,"org"], snpTable[,"chr"], snpTable[,"bp"]),]
    geneholder<-genecounting(snpTable,orthoTableReturn)
    snpTable$genecount<-geneholder$V1
    attachToPermSummary<-NULL
    for(org in as.character(unique(snpTable$org))){
      genecountSummary<-as.data.frame(table(snpTable$genecount[snpTable$org==org]))
      attach<-data.frame(org=org, trait=snpTable$trait[1], present=NA,
                                      GeneCount=as.integer(as.character(genecountSummary$Var1)), LociCount=genecountSummary$Freq)
      attachToPermSummary<-rbind(attachToPermSummary, attach)
    }
    #making a summary table since the permutation tables are too big
    PermutationSummary<-summaryTable(orthoTableReturn)
    PermutationSummary<-rbind(PermutationSummary, attachToPermSummary)
    PermutationSummary$dataset<-"permutations"
    return(PermutationSummary)
  }
}
genecounting<-function(snpTable,orthoTable){
  genecountMaster<-data.frame(V1=integer())
  for(o in as.character(unique(snpTable$org))){
    #print(paste0("working on ",o))
    checking<-orthoTable[,c(paste0(o,"ID"),paste0(o,"_loci")),with=FALSE]
    colnames(checking)<-c("ID","loci")
    checking<-distinct(checking[!is.na(checking$loci),])
    snp_byOrg<-subset(snpTable, org==o)
    countingGenes<-foreach(i=1:nrow(snp_byOrg), .packages = packages.loaded(), .combine = rbind) %do%{
      #print(paste0(o, i))
      genecount<-length(which(grepl(paste0("^",snp_byOrg$loci[i],"$"), checking$loci)))
    }
    genecountMaster<-rbind(genecountMaster,countingGenes)
  }
  return(genecountMaster)
}

summaryTable<-function(orgTable, random=F){
  OrgTraitTable<-foreach(q=as.character(metaTable$orgs), .packages = packages.loaded(), .combine = rbind) %do% {
    if(random==T){
      OrgSub<-distinct(droplevels(subset(orgTable, select=c(paste0(q,"ID"),"trait",paste0(q,"_loci"),"present","permutation"))))
      setnames(OrgSub, c("ID","trait","loci","present","permutation"))
      OrgSub<-OrgSub[!is.na(OrgSub$loci),]
      Orgsplit<-split(OrgSub, list(OrgSub$trait,OrgSub$present,OrgSub$permutation))
      OrgTrait<-foreach(e=Orgsplit, .packages=packages.loaded(), .combine = rbind) %do% {
        return(data.frame(org=q,trait=as.character(e$trait[1]), present=e$present[1], dataset="permutations",
                          GeneCount=as.numeric(length(unique(e$ID))),LociCount=as.numeric(length(unique(e$loci)))))
      }
      OrgTrait<-na.omit(OrgTrait)
      return(OrgTrait)
    }else{
      OrgSub<-distinct(droplevels(subset(orgTable, select=c(paste0(q,"ID"),"trait",paste0(q,"_loci"),"present"))))
      setnames(OrgSub, c("ID","trait","loci","present"))
      OrgSub<-OrgSub[!is.na(OrgSub$loci),]
      Orgsplit<-split(OrgSub, list(OrgSub$trait,OrgSub$present))
      OrgTrait<-foreach(e=Orgsplit, .packages=packages.loaded(), .combine = rbind) %do% {
        return(data.frame(org=q,trait=as.character(e$trait[1]), present=e$present[1], GeneCount=as.numeric(length(unique(e$ID))), 
                          LociCount=as.numeric(length(unique(e$loci)))))
      }
      OrgTrait<-na.omit(OrgTrait)
      return(OrgTrait)
    }
  }
  return(OrgTraitTable)
}
Frequencycreation<-function(split_list, random=F){
  combined_list<-foreach(p=split_list, .packages = packages.loaded(), .combine = rbind) %do% {
    LociSum<-sum(p$LociCount)
    if(random==T){return(ddply(p, .(org, trait, dataset, GeneCount,permutation), summarize, Freq=signif((LociCount/LociSum), 2)))
    }else{return(ddply(p, .(org, trait, dataset, GeneCount), summarize, Freq=signif((LociCount/LociSum), 2)))}
  }
}