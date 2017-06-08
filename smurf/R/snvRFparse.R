#' snvRF-Parse
#'
#' snv prediction model
#' Step 1 Parse
#'
#'  
#' @examples
#' 
#' 
#' @export
snvRFparse = function(a){
  
  #a<-x
  print("Parsing SNVs")
  
  #Filtering SNVs from VRanges object along with features
  mutect<-as.data.frame(a[[1]][isSNV(a[[1]], singleAltOnly=FALSE)],row.names=NULL)
  #mutect<-mutect[mutect$sampleNames == samples(header(a[[2]]))[1],]
  mutect<-mutect[,c("seqnames","start","end", "ref", "alt", "MQ","MQRankSum","TLOD")]
  mutect<-mutect[grep("GL*", mutect$seqnames, invert=TRUE),]
  colnames(mutect)<-c("X.CHROM", "POS", "END_POS_Mutect2", "REF_Mutect2", "ALT_Mutect2", "m2_MQ","m2_MQRankSum","m2_TLOD")
  
  freebayes<-as.data.frame(a[[3]][isSNV(a[[3]], singleAltOnly=FALSE)],row.names=NULL)
  freebayes<-freebayes[,c("seqnames","start","end", "ref", "alt","MQM","MQMR","ODDS")]
  freebayes<-freebayes[grep("GL*", freebayes$seqnames, invert=TRUE),]
  colnames(freebayes)<-c("X.CHROM", "POS", "END_POS_Freebayes", "REF_Freebayes", "ALT_Freebayes","f_MQM","f_MQMR","f_ODDS")
  
  varscan<-as.data.frame(a[[5]][isSNV(a[[5]], singleAltOnly=FALSE)],row.names=NULL)
  varscan<-varscan[,c("seqnames","start","end", "ref", "alt", "SSC","SPV")]
  varscan<-varscan[grep("GL*", varscan$seqnames, invert=TRUE),]
  colnames(varscan)<-c("X.CHROM", "POS", "END_POS_Varscan", "REF_Varscan", "ALT_Varscan", "vs_SSC","vs_SPV")
  
  vardict<-as.data.frame(a[[7]][isSNV(a[[7]], singleAltOnly=FALSE)],row.names=NULL)
  #vardict<-vardict[vardict$sampleNames == samples(header(a[[2]]))[1],] #choose only the tumor sample values (1 or 2) eg.icgc_cll-T
  vardict<-vardict[,c("seqnames","start","end", "ref", "alt", "SSF","SOR")]
  vardict<-vardict[grep("GL*", vardict$seqnames, invert=TRUE),]
  colnames(vardict)<-c("X.CHROM", "POS", "END_POS_Vardict", "REF_Vardict", "ALT_Vardict", "vd_SSF","vd_SOR")
  
  #Filtering SNV's from collapsed vcf
  mutect_SNV<-a[[2]][isSNV(a[[2]], singleAltOnly=FALSE)]
  freebayes_SNV<-a[[4]][isSNV(a[[4]], singleAltOnly=FALSE)]
  varscan_SNV<-a[[6]][isSNV(a[[6]], singleAltOnly=FALSE)]
  vardict_SNV<-a[[8]][isSNV(a[[8]], singleAltOnly=FALSE)]
  
  ## Filtering Passed calls from raw file
  mutect_clean_SNV<-rowRanges(mutect_SNV[mutect_SNV@fixed$FILTER=="PASS",])[,1]
  mutect_coordinates<-as.data.frame(unlist(mutect_clean_SNV),row.names=NULL)
  mutect_coordinates<-mutect_coordinates[1:2]
  mutect_coordinates$FILTER_Mutect2<-"PASS"
  colnames(mutect_coordinates)[1:2]<-c("X.CHROM", "POS")
  mutect_coordinates<-mutect_coordinates[grep("GL*", mutect_coordinates$X.CHROM, invert = TRUE),]
  
  freebayes_clean_SNV<-rowRanges(freebayes_SNV[freebayes_SNV@fixed$FILTER=="PASS",])[,1]
  freebayes_coordinates<-as.data.frame(unlist(freebayes_clean_SNV),row.names=NULL)
  freebayes_coordinates<-freebayes_coordinates[1:2]
  freebayes_coordinates$FILTER_Freebayes<-"PASS"
  colnames(freebayes_coordinates)[1:2]<-c("X.CHROM", "POS")
  freebayes_coordinates<-freebayes_coordinates[grep("GL*", freebayes_coordinates$X.CHROM, invert = TRUE),]
  
  varscan_clean_SNV<-rowRanges(varscan_SNV[varscan_SNV@fixed$FILTER=="PASS",])[,1]
  varscan_coordinates<-as.data.frame(unlist(varscan_clean_SNV),row.names=NULL)
  varscan_coordinates<-varscan_coordinates[1:2]
  varscan_coordinates$FILTER_Varscan<-"PASS"
  colnames(varscan_coordinates)[1:2]<-c("X.CHROM", "POS")
  varscan_coordinates<-varscan_coordinates[grep("GL*", varscan_coordinates$X.CHROM, invert = TRUE),]
  
  vardict_clean_SNV<-rowRanges(vardict_SNV[vardict_SNV@fixed$FILTER=="PASS",])[,1]
  vardict_coordinates<-as.data.frame(unlist(vardict_clean_SNV),row.names=NULL)
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
  mutect<-merge(mutect, coordinates, by=c("X.CHROM", "POS"), all.y=TRUE, sort=TRUE)
  mutect<-unique(mutect)
  mutect<-mutect[,1:(dim(mutect)[2]-3)]

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
  f1<-merge(mutect, freebayes, by = c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  f2<-merge(varscan, vardict,  by = c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  
  final<-merge(f1,f2,  by = c("X.CHROM", "POS"), all=TRUE, sort=TRUE)
  
  raw <- final
  
  #Formatting
  final$REF<-paste(final$REF_Mutect2,final$REF_Freebayes,final$REF_Vardict,final$REF_Varscan, sep ="/")
  final$ALT<-paste(final$ALT_Mutect2,final$ALT_Freebayes,final$ALT_Vardict,final$ALT_Varscan, sep ="/")
  f<-subset(final,select=c("END_POS_Mutect2", "END_POS_Freebayes", "END_POS_Vardict", "END_POS_Varscan"))
  final$END_POS<-apply(f, MARGIN=1, function(x) max(x,na.rm=TRUE))
  final<-final[,c("X.CHROM", "POS", "END_POS", "REF", "ALT", "FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan",
                  "m2_MQ","m2_MQRankSum","m2_TLOD",
                  "f_MQM","f_MQMR","f_ODDS",
                  "vs_SSC","vs_SPV","vd_SSF","vd_SOR")]
  colnames(final)<-c("X.CHROM", "START_POS_REF", "END_POS_REF", "REF_MFVdVs", "ALT_MFVdVs", "FILTER_Mutect2", "FILTER_Freebayes", "FILTER_Vardict", "FILTER_Varscan",
                     "m2_MQ","m2_MQRankSum","m2_TLOD",
                     "f_MQM","f_MQMR","f_ODDS",
                     "vs_SSC","vs_SPV","vd_SSF","vd_SOR")
  
  for(i in c("m2_MQ","m2_MQRankSum","m2_TLOD","f_MQM","f_MQMR","f_ODDS","vs_SSC","vs_SPV","vd_SSF","vd_SOR"))
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
  
  #summary(final)
  #snv_parse1<-final
  
  ##value cutoff for DP


  #final[,13] <- sapply(final[,13],function(x) {ifelse(is.na(x),NA,x >= 15) } )
  
  
  snv_parse <- final
  
  print("Predicting SNVs")

  df <- final[,c("X.CHROM","START_POS_REF","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                 "m2_MQ","m2_MQRankSum","m2_TLOD","f_MQM","f_MQMR","f_ODDS","vs_SSC","vs_SPV","vd_SSF","vd_SOR")]
  write.csv(df, 'snv_parse.csv',row.names = F)
  df <- h2o.importFile(path = normalizePath("snv_parse.csv"),header=T)
  
  smurfdir <- find.package("smurf1.0")
  smurfmodeldir <- paste0(smurfdir, "/data/snv-model-combined-grid")
  snv_model <- h2o.loadModel(path = smurfmodeldir)
  
  #snv_model <- h2o.loadModel(path = "/home/huangwt/R/x86_64-pc-linux-gnu-library/3.3/smbio/data/snv_model-combined-0221")
  #snv_model <- h2o.loadModel(path = "D:/Users/Tyler/Dropbox/Scripts/Real-v3/Results-snv-v3/snv_model-combined-0310")
  #snv_model <- h2o.loadModel(path = "C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.0/data/snv-model-combined-grid")
  
  predicted <- h2o.predict(object = snv_model, newdata = df)
  suppressMessages(file.remove("snv_parse.csv"))
  p <- as.data.frame(predicted)
  table <- final
  
  allpredictedcalls <- cbind(table[,c("X.CHROM","START_POS_REF","END_POS_REF","REF_MFVdVs","ALT_MFVdVs")],p)
  results<- allpredictedcalls[which(allpredictedcalls$predict==TRUE), c(1,2)]
  
  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
  table<-cbind(table, allpredictedcalls$TRUE.)
  truth_RF <- results
  
  # Merge the predicted indels to the sample
  truth_RF <- cbind(truth_RF, 1)
  truth_RF[,3] <- as.logical(truth_RF[,3])
  colnames(truth_RF) <- c("X.CHROM", "START_POS_REF", "TRUTH_RF")
  table1 <- merge(table, truth_RF, by = c("X.CHROM", "START_POS_REF"))
  names(table1)[names(table1) == 'allpredictedcalls$TRUE.'] <- 'SMuRF_score'
  
 
  
  # Generate stats
  stats <- matrix(,nrow = 9, ncol = 1)
  stats <- as.data.frame(stats)
  colnames(stats) <- c("Passed_Calls")
  rownames(stats) <- c("Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "All4", "SMuRF_SNV")
  
  counts <- apply(table[, c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
  
  stats$Passed_Calls[1] <- length(which((table$FILTER_Mutect2==TRUE)))
  stats$Passed_Calls[2] <- length(which((table$FILTER_Freebayes==TRUE)))
  stats$Passed_Calls[3] <- length(which((table$FILTER_Vardict==TRUE)))
  stats$Passed_Calls[4] <- length(which((table$FILTER_Varscan==TRUE)))
  
  
  stats$Passed_Calls[5] <- length(which(counts>=1))
  stats$Passed_Calls[6] <- length(which(counts>=2))
  stats$Passed_Calls[7] <- length(which(counts>=3))
  stats$Passed_Calls[8] <- length(which(counts>=4))
  
  
  stats$Passed_Calls[9] <- nrow(table1)
  
  # Predicted list of mutations
  names(table1)[names(table1) == 'X.CHROM'] <- 'Chr'
  snv_predict <- unique(table1[,c("Chr","START_POS_REF", "END_POS_REF","REF_MFVdVs","ALT_MFVdVs", "SMuRF_score")])
  
  #smbio <- as.list(indel-parse, indel-stats, indel-predict)
  stats<-as.matrix(stats)
  #parse<-as.matrix(table)
  predict<-as.matrix(snv_predict)
  #parse1<-as.matrix(snv_parse1)
  parse<-as.matrix(snv_parse)
  raw<-as.matrix(raw)
  y <- list(stats, predict, parse, raw)
  names(y)<- c("stats_snv", "predicted_snv", "parse_snv","raw_snv")
  return(y)
  
  } else{
    print("Error: There are no predicted SNV calls in this sample.")
    
    parse<-as.matrix(snv_parse)
    raw<-as.matrix(raw)
    y<- list(parse, raw)
    names(y)<- c("parse_snv", "raw_snv")
    return(y)
  } 
}
