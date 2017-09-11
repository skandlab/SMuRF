#' snvRF-Parse
#'
#' snv prediction model
#' Step 2 Predict
#'
#'  
#' @examples
#' 
#' 
#' @export
snvRFpredict = function(b){
  
  final <- b[[1]]
  
  print("Predicting SNVs")
  

  df <- final[,c("X.CHROM","START_POS_REF","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                 "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","f_MQM","f_MQMR","vs_SSC","vs_SPV","vd_SSF")]
  df <- as.h2o(df)

  smurfdir <- find.package("smurf")
  smurfmodeldir <- paste0(smurfdir, "/data/snv-model-combined-grid")
  snv_model <- h2o.loadModel(path = smurfmodeldir)
  
  #snv_model <- h2o.loadModel(path = "D:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/snv-model-combined-grid")
  
  predicted <- h2o.predict(object = snv_model, newdata = df)
  p <- as.data.frame(predicted)

  snv_parse <- cbind(final, p)
  
  results<- snv_parse[which(snv_parse$predict==TRUE),]
  
  
  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
    
    table <- results
    
    names(table)[names(table) == 'TRUE.'] <- 'SMuRF_score'
    names(table)[names(table) == 'X.CHROM'] <- 'Chr'
    
    snv_predict <- unique(table[, !names(table) %in% c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                                       "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","f_MQM","f_MQMR","vs_SSC","vs_SPV","vd_SSF",
                                                       "predict","FALSE.")])
    
  # Generate stats
  stats <- matrix(,nrow = 9, ncol = 1)
  stats <- as.data.frame(stats)
  colnames(stats) <- c("Passed_Calls")
  rownames(stats) <- c("Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "All4", "SMuRF_SNV")
  
  counts <- apply(snv_parse[, c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
  
  stats$Passed_Calls[1] <- length(which((snv_parse$FILTER_Mutect2==TRUE)))
  stats$Passed_Calls[2] <- length(which((snv_parse$FILTER_Freebayes==TRUE)))
  stats$Passed_Calls[3] <- length(which((snv_parse$FILTER_Vardict==TRUE)))
  stats$Passed_Calls[4] <- length(which((snv_parse$FILTER_Varscan==TRUE)))
  
  
  stats$Passed_Calls[5] <- length(which(counts>=1))
  stats$Passed_Calls[6] <- length(which(counts>=2))
  stats$Passed_Calls[7] <- length(which(counts>=3))
  stats$Passed_Calls[8] <- length(which(counts>=4))
  
  
  stats$Passed_Calls[9] <- nrow(snv_predict)
  

  stats<-as.matrix(stats)
  predict<-as.matrix(snv_predict)
  parse<-as.matrix(snv_parse)
  raw<-as.matrix(final)
  y <- list(stats, predict, parse, raw)
  names(y)<- c("stats_snv", "predicted_snv", "parse_snv","raw_snv")
  
  return(y)
  
  } else{
    
    print("Error: There are no predicted SNV calls in this sample.")
    
    parse<-as.matrix(snv_parse)
    raw<-as.matrix(final)
    y<- list(parse, raw)
    names(y)<- c("parse_snv", "raw_snv")
    
    return(y)
  } 
  
  
  
}
