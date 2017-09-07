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
indelRFparse = function(a){
  
  print("Parsing INDELs")
  
  #Filtering Indels from VRanges object
  mutect<-as.data.frame(a[[1]][!isSNV(a[[1]], singleAltOnly=FALSE)],row.names=NULL)
  mutect<-mutect[,c("seqnames","start","end", "ref", "alt", "MQ","MQRankSum","NLOD","TLOD")]
  mutect<-mutect[grep("GL*", mutect$seqnames, invert=TRUE),]
  colnames(mutect)<-c("X.CHROM", "POS", "END_POS_Mutect2", "REF_Mutect2", "ALT_Mutect2","m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD")
  
  freebayes<-as.data.frame(a[[3]][!isSNV(a[[3]], singleAltOnly=FALSE)],row.names=NULL)
  freebayes<-freebayes[,c("seqnames","start","end", "ref", "alt", "LEN")]
  freebayes<-freebayes[grep("GL*", freebayes$seqnames, invert=TRUE),]
  colnames(freebayes)<-c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes","f_LEN")
  
  varscan<-as.data.frame(a[[5]][!isSNV(a[[5]], singleAltOnly=FALSE)],row.names=NULL)
  varscan<-varscan[,c("seqnames","start","end", "ref", "alt", "SSC","SPV")]
  varscan<-varscan[grep("GL*", varscan$seqnames, invert=TRUE),]
  colnames(varscan)<-c("X.CHROM", "POS", "END_POS_Varscan", "REF_Varscan", "ALT_Varscan", "vs_SSC","vs_SPV")
  
  vardict<-as.data.frame(a[[7]][!isSNV(a[[7]], singleAltOnly=FALSE)],row.names=NULL)
  vardict<-vardict[,c("seqnames","start","end", "ref", "alt", "SSF","MSI","SOR")]
  vardict<-vardict[grep("GL*", vardict$seqnames, invert=TRUE),]
  colnames(vardict)<-c("X.CHROM", "POS", "END_POS_Vardict", "REF_Vardict", "ALT_Vardict", "vd_SSF","vd_MSI","vd_SOR")
  
  
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
  #varscan_clean_Indel<-rowRanges(varscan_Indel[varscan_Indel@fixed$FILTER=="PASS",])[,1]
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
  
  #Caller names
  mutect_coordinates$m_caller <- "Mutect2"
  freebayes_coordinates$f_caller <- "Freebayes"
  varscan_coordinates$vs_caller <- "Varscan"
  vardict_coordinates$vd_caller <- "Vardict"

  #All Pass coordinates(i.e., coordinates passed by at least one caller)
  f1<-merge(mutect_coordinates, freebayes_coordinates, by=c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  f2<-merge(varscan_coordinates, vardict_coordinates, by=c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  coordinates<-merge(f1, f2, by=c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  coordinates<-unique(coordinates)
  
  #### Filtering for PASSED INDELS ####
  
  #Pulling out calls from VRanges object for each caller using coordinates
  m2<-dim(mutect)[2]
  mutect<-merge(mutect, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  mutect<-unique(mutect)
  M_2<-which( colnames(mutect)=="FILTER_Mutect2" )
  M_call<-which( colnames(mutect)=="m_caller" )
  mutect<-mutect[,c(1:m2,M_2,M_call)]
  
  fb<-dim(freebayes)[2]
  freebayes<-merge(freebayes, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  freebayes<-unique(freebayes)
  F_B<-which( colnames(freebayes)=="FILTER_Freebayes" )
  F_call<-which( colnames(freebayes)=="f_caller" )
  freebayes<-freebayes[,c(1:fb,F_B,F_call)]
  
  vs<-dim(varscan)[2]
  varscan<-merge(varscan, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  varscan<-unique(varscan)
  V_S<-which( colnames(varscan)=="FILTER_Varscan" )
  VS_call<-which( colnames(varscan)=="vs_caller" )
  varscan<-varscan[,c(1:vs,V_S,VS_call)]
  
  vd<-dim(vardict)[2]
  vardict<-merge(vardict, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  vardict<-unique(vardict)
  V_D<-which( colnames(vardict)=="FILTER_Vardict" )
  VD_call<-which( colnames(vardict)=="vd_caller" )
  vardict<-vardict[,c(1:vd,V_D,VD_call)]
  
  #Merging all calls from each caller 
  f1<-merge(mutect, freebayes, by = c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  f2<-merge(varscan, vardict,  by = c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  
  final<-merge(f1,f2,  by = c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  
  #Compiling REF
  for (i in 1:nrow(final)){
    if (!is.na(final$REF_Mutect2[i])){
      final$REF[i] <- final$REF_Mutect2[i]
    } else if (!is.na(final$REF_Freebayes[i])){
      final$REF[i] <- final$REF_Freebayes[i]
    } else if (!is.na(final$REF_Vardict[i])){
      final$REF[i] <- final$REF_Vardict[i]
    } else if (!is.na(final$REF_Varscan[i])){
      final$REF[i] <- final$REF_Varscan[i]
    } 
  }
  
  #Compiling ALT
  for (i in 1:nrow(final)){
    if (!is.na(final$ALT_Mutect2[i])){
      final$ALT[i] <- final$ALT_Mutect2[i]
    } else if (!is.na(final$ALT_Freebayes[i])){
      final$ALT[i] <- final$ALT_Freebayes[i]
    } else if (!is.na(final$ALT_Vardict[i])){
      final$ALT[i] <- final$ALT_Vardict[i]
    } else if (!is.na(final$ALT_Varscan[i])){
      final$ALT[i] <- final$ALT_Varscan[i]
    } 
  }
  
  
  final$CALLS<-paste(final$m_caller,final$f_caller,final$vs_caller,final$vd_caller, sep ="/")
  
  #Formatting
  final$REF2<-paste(final$REF_Mutect2,final$REF_Freebayes,final$REF_Vardict,final$REF_Varscan, sep ="/")
  final$ALT2<-paste(final$ALT_Mutect2,final$ALT_Freebayes,final$ALT_Vardict,final$ALT_Varscan, sep ="/")
  f<-subset(final,select=c("END_POS_Mutect2", "END_POS_Freebayes", "END_POS_Vardict", "END_POS_Varscan"))
  final$END_POS<-apply(f, MARGIN=1, function(x) max(x,na.rm=TRUE))
  final<-final[,c("X.CHROM", "POS", "END_POS", "REF", "ALT", "REF2", "ALT2", 
                  "CALLS",
                  "FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan",
                  "m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD",
                  "f_LEN",
                  "vs_SSC","vs_SPV","vd_SSF","vd_MSI","vd_SOR")]
  colnames(final)<-c("X.CHROM", "START_POS_REF", "END_POS_REF", "REF", "ALT", "REF_MFVdVs", "ALT_MFVdVs", 
                     "CALLS",
                     "FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan",
                     "m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD",
                     "f_LEN",
                     "vs_SSC","vs_SPV","vd_SSF","vd_MSI","vd_SOR")
  
  #Removing Duplicates
  for(i in c("m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD",
             "f_LEN",
             "vs_SSC","vs_SPV","vd_SSF","vd_MSI", "vd_SOR"))
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
  
  raw <- final
  
  indel_parse <- final


  print("Predicting INDELs")
  
  df <- final[,c("X.CHROM","START_POS_REF","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                 "m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD","f_LEN","vs_SSC","vs_SPV","vd_SSF","vd_MSI","vd_SOR")]
  df <- as.h2o(df)

  smurfdir <- find.package("smurf")
  smurfmodeldir <- paste0(smurfdir, "/data/indel-model-combined-grid")
  indel_model <- h2o.loadModel(path = smurfmodeldir)
  
  #indel_model <- h2o.loadModel(path = "D:/Users/Tyler/Dropbox/Sscripts/smurf/smurf1.2/smurf/data/indel-model-combined-grid")
  
  predicted <- h2o.predict(object = indel_model, newdata = df)
  p <- as.data.frame(predicted)
  table <- final
  
  indel_parse <- cbind(indel_parse, p)
  
  results<- indel_parse[which(indel_parse$predict==TRUE),]
  
  
  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
    
    
    table <- results
    
    names(table)[names(table) == 'TRUE.'] <- 'SMuRF_score'
    
    # Predicted list of mutations
    names(table)[names(table) == 'X.CHROM'] <- 'Chr'
    indel_predict <- unique(table[, !names(table) %in% c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                                         "m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD","f_LEN","vs_SSC","vs_SPV","vd_SSF","vd_MSI","vd_SOR",
                                                         "predict","FALSE.")])
    
  
  
  # Generate stats
   stats <- matrix(,nrow = 9, ncol = 1)
   stats <- as.data.frame(stats)
   colnames(stats) <- c("Passed_Calls")
   rownames(stats) <- c("Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "All4", "SMuRF_INDEL")
  
   counts <- apply(table[, c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
  
   stats$Passed_Calls[1] <- length(which((table$FILTER_Mutect2==TRUE)))
   stats$Passed_Calls[2] <- length(which((table$FILTER_Freebayes==TRUE)))
   stats$Passed_Calls[3] <- length(which((table$FILTER_Vardict==TRUE)))
   stats$Passed_Calls[4] <- length(which((table$FILTER_Varscan==TRUE)))
   
   stats$Passed_Calls[5] <- length(which(counts>=1))
   stats$Passed_Calls[6] <- length(which(counts>=2))
   stats$Passed_Calls[7] <- length(which(counts>=3))
   stats$Passed_Calls[8] <- length(which(counts>=4))
  
   stats$Passed_Calls[9] <- nrow(table)
  
   
   
   stats<-as.matrix(stats)
   predict<-as.matrix(indel_predict)
   parse<-as.matrix(indel_parse)
   raw<-as.matrix(raw)
   z<- list(stats, predict, parse, raw)
   names(z)<- c("stats_indel", "predicted_indel", "parse_indel", "raw_indel")
   
   
   return(z)
   
   } else{
    print("Error: There are no predicted indel calls in this sample.")
    
     parse<-as.matrix(indel_parse)
     raw<-as.matrix(raw)
     z<- list(parse, raw)
     names(z)<- c("parse_indel", "raw_indel")
     return(z)
   } 
  
  
}
