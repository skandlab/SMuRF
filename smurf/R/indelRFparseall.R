#' indelRF-Parse
#'
#' Indel prediction model
#' Step 1 Parse
#'
#'  
#' @examples
#' 
#' 
#' @export
indelRFparseall = function(a){
  
  print("indelRF Parsing")

  #Filtering Indels from VRanges object
  mutect<-as.data.frame(a[[1]][!isSNV(a[[1]], singleAltOnly=FALSE)],row.names=NULL)
  #mutect<-mutect[mutect$sampleNames == samples(header(a[[2]]))[1],]
  mutect<-mutect[,c("seqnames","start","end", "ref", "alt", "BaseQRankSum","ClippingRankSum","FS","HCNT","MQ","MQRankSum","NLOD","ReadPosRankSum","TLOD")]
  mutect<-mutect[grep("GL*", mutect$seqnames, invert=TRUE),]
  colnames(mutect)<-c("X.CHROM", "POS", "END_POS_Mutect2", "REF_Mutect2", "ALT_Mutect2", "m2_BaseQRankSum","m2_ClippingRankSum","m2_FS","m2_HCNT","m2_MQ","m2_MQRankSum","m2_NLOD","m2_ReadPosRankSum","m2_TLOD")
  
  freebayes<-as.data.frame(a[[3]][!isSNV(a[[3]], singleAltOnly=FALSE)],row.names=NULL)
  freebayes<-freebayes[,c("seqnames","start","end", "ref", "alt", "AB","ABP","EPP","EPPR","GTI","LEN","MQM","MQMR","NS","NUMALT","ODDS","PAIRED","PAIREDR","RPP","SAP","SRP")]
  freebayes<-freebayes[grep("GL*", freebayes$seqnames, invert=TRUE),]
  #freebayes_GQT<-freebayes[freebayes$sampleNames == samples(header(a[[2]]))[1],] #[1] for tumor
  #freebayes_GQT<-freebayes_GQT[,c("seqnames","start","end", "ref", "alt","QUAL","GQ")]
  #colnames(freebayes_GQT)<-c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes","QUAL","GQ-T")
  #freebayes_GQN<-freebayes[freebayes$sampleNames == samples(header(a[[2]]))[2],]
  #freebayes_GQN<-freebayes_GQN[,c("sampleNames","GQ")]
  #colnames(freebayes_GQN)<-c("sampleNames","GQ-N")
  #freebayes<-cbind(freebayes_GQT, freebayes_GQN)
  #freebayes_GQN<-freebayes[freebayes$sampleNames == samples(header(a[[2]]))[2],] #[2] for normal
  #freebayes_GQN<-freebayes_GQN[,c("seqnames","start","end", "ref", "alt","QUAL","GQ")]
  #colnames(freebayes_GQN)<-c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes","GQ-N")
  #freebayes <- freebayes_GQN
  #freebayes<-freebayes[,c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes")]
  colnames(freebayes)<-c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes","f_AB","f_ABP","f_EPP","f_EPPR","f_GTI","f_LEN","f_MQM","f_MQMR","f_NS","f_NUMALT","f_ODDS","f_PAIRED","f_PAIREDR","f_RPP","f_SAP","f_SRP")
  
  varscan<-as.data.frame(a[[5]][!isSNV(a[[5]], singleAltOnly=FALSE)],row.names=NULL)
  varscan<-varscan[,c("seqnames","start","end", "ref", "alt", "SSC","GPV","SS","SPV")]
  varscan<-varscan[grep("GL*", varscan$seqnames, invert=TRUE),]
  colnames(varscan)<-c("X.CHROM", "POS", "END_POS_Varscan", "REF_Varscan", "ALT_Varscan", "vs_SSC","vs_GPV","vs_SS","vs_SPV")
  
  vardict<-as.data.frame(a[[7]][!isSNV(a[[7]], singleAltOnly=FALSE)],row.names=NULL)
  vardict<-vardict[,c("seqnames","start","end", "ref", "alt", "SSF","MSI","SOR")]
  vardict<-vardict[grep("GL*", vardict$seqnames, invert=TRUE),]
  colnames(vardict)<-c("X.CHROM", "POS", "END_POS_Vardict", "REF_Vardict", "ALT_Vardict", "vd_SSF","vd_MSI","vd_SOR")
  
  
  #### Coding annotation ####
  
  #Filtering Indels from VRanges object
  mutect_cds<-as.data.frame(a[[9]][!isSNV(a[[9]], singleAltOnly=FALSE)],row.names=NULL)
  mutect_cds<-mutect_cds[,c("seqnames","start","end", "ref", "alt", "BaseQRankSum","ClippingRankSum","FS","HCNT","MQ","MQRankSum","NLOD","ReadPosRankSum","TLOD")]
  mutect_cds<-mutect_cds[grep("GL*", mutect_cds$seqnames, invert=TRUE),]
  colnames(mutect_cds)<-c("X.CHROM", "POS", "END_POS_Mutect2", "REF_Mutect2", "ALT_Mutect2", "m2_BaseQRankSum","m2_ClippingRankSum","m2_FS","m2_HCNT","m2_MQ","m2_MQRankSum","m2_NLOD","m2_ReadPosRankSum","m2_TLOD")

  freebayes_cds<-as.data.frame(a[[11]][!isSNV(a[[11]], singleAltOnly=FALSE)],row.names=NULL)
  freebayes_cds<-freebayes_cds[,c("seqnames","start","end", "ref", "alt", "AB","ABP","EPP","EPPR","GTI","LEN","MQM","MQMR","NS","NUMALT","ODDS","PAIRED","PAIREDR","RPP","SAP","SRP")]
  freebayes_cds<-freebayes_cds[grep("GL*", freebayes_cds$seqnames, invert=TRUE),]
  colnames(freebayes_cds)<-c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes","f_AB","f_ABP","f_EPP","f_EPPR","f_GTI","f_LEN","f_MQM","f_MQMR","f_NS","f_NUMALT","f_ODDS","f_PAIRED","f_PAIREDR","f_RPP","f_SAP","f_SRP")

  varscan_cds<-as.data.frame(a[[13]][!isSNV(a[[13]], singleAltOnly=FALSE)],row.names=NULL)
  varscan_cds<-varscan_cds[,c("seqnames","start","end", "ref", "alt", "SSC","GPV","SS","SPV")]
  varscan_cds<-varscan_cds[grep("GL*", varscan_cds$seqnames, invert=TRUE),]
  colnames(varscan_cds)<-c("X.CHROM", "POS", "END_POS_Varscan", "REF_Varscan", "ALT_Varscan", "vs_SSC","vs_GPV","vs_SS","vs_SPV")

  vardict_cds<-as.data.frame(a[[15]][!isSNV(a[[15]], singleAltOnly=FALSE)],row.names=NULL)
  vardict_cds<-vardict_cds[,c("seqnames","start","end", "ref", "alt", "SSF","MSI","SOR")]
  vardict_cds<-vardict_cds[grep("GL*", vardict_cds$seqnames, invert=TRUE),]
  colnames(vardict_cds)<-c("X.CHROM", "POS", "END_POS_Vardict", "REF_Vardict", "ALT_Vardict", "vd_SSF","vd_MSI","vd_SOR")

  #assign column and then merge with original
  mutect_cds$Region <- "CDS"
  mutect_cds <- mutect_cds[,c(1,2,dim(mutect_cds)[2])]
  mutect_cds <- unique(mutect_cds)
  mutect<-merge(mutect, mutect_cds, by=c("X.CHROM", "POS"), all.x=TRUE, sort=TRUE)
  #mutect[which(is.na(mutect[,dim(mutect)[2]])), dim(mutect)[2]] <- "Non-coding"

  freebayes_cds$Region <- "CDS"
  freebayes_cds <- freebayes_cds[,c(1,2,dim(freebayes_cds)[2])]
  freebayes_cds <- unique(freebayes_cds)
  freebayes<-merge(freebayes, freebayes_cds, by=c("X.CHROM", "POS"), all.x=TRUE, sort=TRUE)
  #freebayes[which(is.na(freebayes[,dim(freebayes)[2]])), dim(freebayes)[2]] <- "Non-coding"

  varscan_cds$Region <- "CDS"
  varscan_cds <- varscan_cds[,c(1,2,dim(varscan_cds)[2])]
  varscan_cds <- unique(varscan_cds)
  varscan<-merge(varscan, varscan_cds, by=c("X.CHROM", "POS"), all.x=TRUE, sort=TRUE)
  #varscan[which(is.na(varscan[,dim(varscan)[2]])), dim(varscan)[2]] <- "Non-coding"

  vardict_cds$Region <- "CDS"
  vardict_cds <- vardict_cds[,c(1,2,dim(vardict_cds)[2])]
  vardict_cds <- unique(vardict_cds)
  vardict<-merge(vardict, vardict_cds, by=c("X.CHROM", "POS"), all.x=TRUE, sort=TRUE)
  #vardict[which(is.na(vardict[,dim(vardict)[2]])), dim(vardict)[2]] <- "Non-coding"

  
  #### Coordinates filtering ####
  
  
  #Filtering Indels from collapsed vcf
  mutect_Indel<-a[[2]][!isSNV(a[[2]], singleAltOnly=FALSE)]
  freebayes_Indel<-a[[4]][!isSNV(a[[4]], singleAltOnly=FALSE)]
  varscan_Indel<-a[[6]][!isSNV(a[[6]], singleAltOnly=FALSE)]
  vardict_Indel<-a[[8]][!isSNV(a[[8]], singleAltOnly=FALSE)]

  ## Filtering Passed calls from raw file
  mutect_clean_Indel<-rowRanges(mutect_Indel[mutect_Indel@fixed$FILTER=="PASS",])[,1]
  mutect_coordinates<-as.data.frame(unlist(mutect_clean_Indel),row.names=NULL)
  mutect_coordinates<-mutect_coordinates[1:2]
  mutect_coordinates$FILTER_Mutect2<-"PASS"
  colnames(mutect_coordinates)[1:2]<-c("X.CHROM", "POS")
  mutect_coordinates<-mutect_coordinates[grep("GL*", mutect_coordinates$X.CHROM, invert = TRUE),]
  
  freebayes_clean_Indel<-rowRanges(freebayes_Indel[freebayes_Indel@fixed$FILTER=="PASS",])[,1]
  freebayes_coordinates<-as.data.frame(unlist(freebayes_clean_Indel),row.names=NULL)
  freebayes_coordinates<-freebayes_coordinates[1:2]
  freebayes_coordinates$FILTER_Freebayes<-"PASS"
  colnames(freebayes_coordinates)[1:2]<-c("X.CHROM", "POS")
  freebayes_coordinates<-freebayes_coordinates[grep("GL*", freebayes_coordinates$X.CHROM, invert = TRUE),]
  
  #we are including the SpvFreq as PASSED calls
  varscan_clean_Indel<-rowRanges(varscan_Indel[varscan_Indel@fixed$FILTER=="PASS" | varscan_Indel@fixed$FILTER=="SpvFreq",])[,1]
  varscan_coordinates<-as.data.frame(unlist(varscan_clean_Indel),row.names=NULL)
  varscan_coordinates<-varscan_coordinates[1:2]
  varscan_coordinates$FILTER_Varscan<-"PASS"
  colnames(varscan_coordinates)[1:2]<-c("X.CHROM", "POS")
  varscan_coordinates<-varscan_coordinates[grep("GL*", varscan_coordinates$X.CHROM, invert = TRUE),]
  
  vardict_clean_Indel<-rowRanges(vardict_Indel[vardict_Indel@fixed$FILTER=="PASS",])[,1]
  vardict_coordinates<-as.data.frame(unlist(vardict_clean_Indel),row.names=NULL)
  vardict_coordinates<-vardict_coordinates[1:2]
  vardict_coordinates$FILTER_Vardict<-"PASS"
  colnames(vardict_coordinates)[1:2]<-c("X.CHROM", "POS")
  vardict_coordinates<-vardict_coordinates[grep("GL*", vardict_coordinates$X.CHROM, invert = TRUE),]

  #All Pass coordinates(i.e., coordinates passed by at least one caller)
  f1<-merge(mutect_coordinates, freebayes_coordinates, by=c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  f2<-merge(varscan_coordinates, vardict_coordinates, by=c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  coordinates<-merge(f1, f2, by=c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  coordinates<-unique(coordinates)

  #Pulling out calls from VRanges object for each caller using coordinates
  m2<-dim(mutect)[2]
  mutect<-merge(mutect, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  mutect<-unique(mutect)
  M_2<-which( colnames(mutect)=="FILTER_Mutect2" )
  mutect<-mutect[,c(1:m2,M_2)]
  
  fb<-dim(freebayes)[2]
  freebayes<-merge(freebayes, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  freebayes<-unique(freebayes)
  F_B<-which( colnames(freebayes)=="FILTER_Freebayes" )
  freebayes<-freebayes[,c(1:fb,F_B)]
  
  vs<-dim(varscan)[2]
  varscan<-merge(varscan, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  varscan<-unique(varscan)
  V_S<-which( colnames(varscan)=="FILTER_Varscan" )
  varscan<-varscan[,c(1:vs,V_S)]
  
  vd<-dim(vardict)[2]
  vardict<-merge(vardict, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  vardict<-unique(vardict)
  V_D<-which( colnames(vardict)=="FILTER_Vardict" )
  vardict<-vardict[,c(1:vd,V_D)]
  
  #Merging all calls from each caller 
  f1<-merge(mutect, freebayes, by = c("X.CHROM", "POS", "Region"), all=TRUE, sort=TRUE)
  f2<-merge(varscan, vardict,  by = c("X.CHROM", "POS", "Region"), all=TRUE, sort=TRUE)
  
  final<-merge(f1,f2,  by = c("X.CHROM", "POS", "Region"), all=TRUE, sort=TRUE)
  final[which(is.na(final[,c("Region")])), c("Region")] <- "Non-coding"
  
  raw <- final
  
  
  
  #Formatting
  final$REF<-paste(final$REF_Mutect2,final$REF_Freebayes,final$REF_Vardict,final$REF_Varscan, sep ="/")
  final$ALT<-paste(final$ALT_Mutect2,final$ALT_Freebayes,final$ALT_Vardict,final$ALT_Varscan, sep ="/")
  f<-subset(final,select=c("END_POS_Mutect2", "END_POS_Freebayes", "END_POS_Vardict", "END_POS_Varscan"))
  final$END_POS<-apply(f, MARGIN=1, function(x) max(x,na.rm=TRUE))
  final<-final[,c("X.CHROM", "POS", "END_POS", "REF", "ALT", "Region", "FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan",
                  "m2_BaseQRankSum","m2_ClippingRankSum","m2_FS","m2_HCNT","m2_MQ","m2_MQRankSum","m2_NLOD","m2_ReadPosRankSum","m2_TLOD",
                  "f_AB","f_ABP","f_EPP","f_EPPR","f_GTI","f_LEN","f_MQM","f_MQMR","f_NS","f_NUMALT","f_ODDS","f_PAIRED","f_PAIREDR","f_RPP","f_SAP","f_SRP",
                  "vs_SSC","vs_GPV","vs_SS","vs_SPV","vd_SSF","vd_MSI","vd_SOR")]
  colnames(final)<-c("X.CHROM", "START_POS_REF", "END_POS_REF", "REF_MFVdVs", "ALT_MFVdVs", "REGION", "FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan",
                     "m2_BaseQRankSum","m2_ClippingRankSum","m2_FS","m2_HCNT","m2_MQ","m2_MQRankSum","m2_NLOD","m2_ReadPosRankSum","m2_TLOD",
                     "f_AB","f_ABP","f_EPP","f_EPPR","f_GTI","f_LEN","f_MQM","f_MQMR","f_NS","f_NUMALT","f_ODDS","f_PAIRED","f_PAIREDR","f_RPP","f_SAP","f_SRP",
                     "vs_SSC","vs_GPV","vs_SS","vs_SPV","vd_SSF","vd_MSI","vd_SOR")
  
  
  for(i in 11:dim(final)[2])
  {
    final[,i]<-as.numeric(final[,i])
  }
  
  #Removing Duplicates
  final<-unique(final)
  
  for(i in c("FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan"))
  {
    final[which(is.na(final[,i])), i] <- FALSE
    final[which(final[,i]!="PASS"), i] <- "FALSE"
    final[which(final[,i]=="PASS"), i] <- "TRUE"
    final[,i] <- as.logical(final[,i])
  }

  final$vd_SOR[which(is.infinite(final$vd_SOR)==TRUE)] <- 0
  
  
  final<-as.matrix(final)
  z <- list(final)
  names(z)<- c("allfeatures_indel")
  return(z)
  
  #write.table(final , file = "indel-parse-allfeatures.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  #x = read.table('indel-parse-allfeatures.txt',header=T,sep = "\t")
  #write.csv(x, 'indel-parse-allfeatures.csv',row.names = F)
  
  
}
