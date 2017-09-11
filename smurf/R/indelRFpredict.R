#' indelRF-Parse
#'
#' Indel prediction model
#' Step 2 Predict
#'
#'  
#' @examples
#' 
#' 
#' @export
indelRFpredict = function(c){
  
  final <- c[[1]]
  
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

  indel_parse <- cbind(final, p)
  
  results<- indel_parse[which(indel_parse$predict==TRUE),]
  
  
  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
    
    
    table <- results
    
    names(table)[names(table) == 'TRUE.'] <- 'SMuRF_score'
    names(table)[names(table) == 'X.CHROM'] <- 'Chr'
    
    indel_predict <- unique(table[, !names(table) %in% c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                                                         "m2_MQ","m2_MQRankSum","m2_NLOD","m2_TLOD","f_LEN","vs_SSC","vs_SPV","vd_SSF","vd_MSI","vd_SOR",
                                                         "predict","FALSE.")])
    
  
  
  # Generate stats
   stats <- matrix(,nrow = 9, ncol = 1)
   stats <- as.data.frame(stats)
   colnames(stats) <- c("Passed_Calls")
   rownames(stats) <- c("Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "All4", "SMuRF_INDEL")
  
   counts <- apply(indel_parse[, c("FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
  
   stats$Passed_Calls[1] <- length(which((indel_parse$FILTER_Mutect2==TRUE)))
   stats$Passed_Calls[2] <- length(which((indel_parse$FILTER_Freebayes==TRUE)))
   stats$Passed_Calls[3] <- length(which((indel_parse$FILTER_Vardict==TRUE)))
   stats$Passed_Calls[4] <- length(which((indel_parse$FILTER_Varscan==TRUE)))
   
   stats$Passed_Calls[5] <- length(which(counts>=1))
   stats$Passed_Calls[6] <- length(which(counts>=2))
   stats$Passed_Calls[7] <- length(which(counts>=3))
   stats$Passed_Calls[8] <- length(which(counts>=4))
  
   stats$Passed_Calls[9] <- nrow(indel_predict)
  
   
   
   stats<-as.matrix(stats)
   predict<-as.matrix(indel_predict)
   parse<-as.matrix(indel_parse)
   raw<-as.matrix(final)
   z<- list(stats, predict, parse, raw)
   names(z)<- c("stats_indel", "predicted_indel", "parse_indel", "raw_indel")
   
   
   return(z)
   
   } else{
    print("Error: There are no predicted indel calls in this sample.")
    
     parse<-as.matrix(indel_parse)
     raw<-as.matrix(final)
     z<- list(parse, raw)
     names(z)<- c("parse_indel", "raw_indel")
     return(z)
   } 
  
  
}
