#' snvRF-Predict
#'
#' snv prediction model
#' Step 2 Predict
#'
#'  
#' @examples
#' 
#' 
#' @export
snvRFpredict = function(parsevcf){
  
  final <- parsevcf[[1]]
  #final <- parse_snv
  
  print("Predicting SNVs")
  
  # df <- final[,c("X.CHROM","START_POS_REF","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
  #                "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","f_MQM","f_MQMR","vs_SSC","vs_SPV","vd_SSF")]
  
  df <- as.h2o(final)

  smurfdir <- find.package("smurf")
  smurfmodeldir <- paste0(smurfdir, "/data/snv-model-combined-grid")
  snv_model <- h2o.loadModel(path = smurfmodeldir)
  
  #h2o.find_threshold_by_max_metric(h2o.performance(snv_model), "f1") #0.01860056
  #snv_model <- h2o.loadModel(path = "C:/Users/Tyler/Dropbox/Scripts/smurf/smurf1.2/smurf/data/snv-model-combined-grid")
  
  predicted <- h2o.predict(object = snv_model, newdata = df)
  p <- as.data.frame(predicted)

  snv_parse <- cbind(final, p)
  
  #set our threshold for SMuRF passing score based on 0.9 Recall
  # snv_parse$predict_adjusted <- FALSE
  # snv_parse[which(snv_parse$TRUE.>=0.005),"predict_adjusted"] <- TRUE
  # results<- snv_parse[which(snv_parse$TRUE.>=0.005),] 
  
  results<- snv_parse[which(snv_parse$predict==TRUE),]
  
  # snv_parse <- unique(snv_parse[, !names(snv_parse) %in% c("m2_refDepth", "m2_altDepth", "m2_AF",
  #                                                  "f_refDepth", "f_altDepth",
  #                                                  "vs_totalDepth", "vs_AF",
  #                                                  "vd_refDepth", "vd_altDepth", "vd_AF",
  #                                                  "REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict",
  #                                                  "ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict")])
  
  
  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
    
    table <- results
    
    names(table)[names(table) == 'TRUE.'] <- 'SMuRF_score'
    #names(table)[names(table) == 'X.CHROM'] <- 'Chr'
    
    # snv_predict <- unique(table[, !names(table) %in% c("m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","m2_BaseQRankSum","m2_FS",
    #                                                    "f_MQM","f_ODDS","vs_SSC","vd_SSF",
    #                                                    "predict","FALSE.")])
    
    snv_predict <- table[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                            "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                            "Alt_Allele_Freq",
                            "N_refDepth","N_altDepth","T_refDepth","T_altDepth",
                            "SMuRF_score")]
    
    
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
  #predict<-as.matrix(snv_predict)
  parse<-as.matrix(snv_parse)
  #raw<-as.matrix(final)
  #y <- list(stats, snv_predict, parse)
  #names(y)<- c("stats_snv", "predicted_snv", "parse_snv")
  
  return(list(stats_snv=stats, predicted_snv=snv_predict, parse_snv=parse))
  
  } else{
    
    print("Warning: There are no predicted SNV calls in this sample.")
    
    parse<-as.matrix(snv_parse)
    #raw<-as.matrix(final)
    # y<- list(parse)
    # names(y)<- c("parse_snv")
    
    return(list(parse_snv=parse))
    #return()
  } 
  
  
  
}
