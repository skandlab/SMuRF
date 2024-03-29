#' snvRF-Predict
#'
#' snv prediction model
#' Step 2 Predict
#'
#'@param parsevcf List object containing 1.snv-parse and 2.indel-parse matrices.
#'@param snv.cutoff Threshold for true predictions by SMuRF
#'@param fixed.snv.cutoff debug for when threshold is zero. 
#'  
#'@examples
#' 
#' 
#'@export
snvRFpredict = function(parsevcf, snv.cutoff, fixed.snv.cutoff=F){
  
  final <- parsevcf[[1]]
  
  print("Predicting SNVs")
  
  df <- as.h2o(final)

  smurfdir <- find.package("smurf")
  
  if(exists('snv.cutoff')==F) {snv.cutoff = 'default'}
  
  # smurfmodeldir <- paste0(smurfdir, "/data/smurf-snv-model-v6") #SMuRFv1.6
  # snv_model <- h2o.loadModel(path = smurfmodeldir)
  #smurfmodeldir <- paste0(smurfdir, "/data/snv-v7.zip") #SMuRFv2.0
  smurfmodeldir <- paste0(smurfdir, "/data/snv-v3-version3.zip") #SMuRFv3.0
  snv_model <- h2o.import_mojo(smurfmodeldir)
  
  #Define cutoffs
  if (snv.cutoff == 'default') {
    # cutoff = h2o.find_threshold_by_max_metric(h2o.performance(snv_model), "f1") #0.4627947 v7-2
    # cutoff = h2o.find_threshold_by_max_metric(h2o.performance(snv_model), "f1") #0.4587242
    #cutoff = 0.0822 #High sensitivity Recall >0.99 v7-2, version 2
     cutoff = 0.351791676448298	  
    # cutoff = 0.12760258 #High sensitivity Recall >0.99 
  } else if (snv.cutoff != 'default') {
    cutoff = snv.cutoff
  } else if (snv.cutoff == 0) {
    cutoff = 0
    fixed.snv.cutoff = T
  } 
  
  
  predicted <- h2o.predict(object = snv_model, newdata = df)
  p <- as.data.frame(predicted)

  snv_parse <- cbind(final, p)
  
  results <- snv_parse[which(snv_parse$TRUE.>cutoff),]
  
  if (fixed.snv.cutoff == T) {
    results <- snv_parse
  }
  
  
  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
    
    table <- results
    
    names(table)[names(table) == 'TRUE.'] <- 'SMuRF_score'

    snv_predict <- table[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                            "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                            "Sample_Name",
                            "Alt_Allele_Freq",
                            "N_refDepth","N_altDepth","T_refDepth","T_altDepth",
                            "SMuRF_score")]
    
    
  # Generate stats
  stats <- matrix(,nrow = 11, ncol = 1)
  stats <- as.data.frame(stats)
  colnames(stats) <- c("Passed_Calls")
  rownames(stats) <- c("Strelka2", "Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "Atleast4", "All5", "SMuRF_SNV")
  
  counts <- apply(snv_parse[, c("FILTER_Strelka2","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
  
  stats$Passed_Calls[1] <- length(which((snv_parse$FILTER_Strelka2==TRUE)))
  stats$Passed_Calls[2] <- length(which((snv_parse$FILTER_Mutect2==TRUE)))
  stats$Passed_Calls[3] <- length(which((snv_parse$FILTER_Freebayes==TRUE)))
  stats$Passed_Calls[4] <- length(which((snv_parse$FILTER_Vardict==TRUE)))
  stats$Passed_Calls[5] <- length(which((snv_parse$FILTER_Varscan==TRUE)))
  
  stats$Passed_Calls[6] <- length(which(counts>=1))
  stats$Passed_Calls[7] <- length(which(counts>=2))
  stats$Passed_Calls[8] <- length(which(counts>=3))
  stats$Passed_Calls[9] <- length(which(counts>=4))
  stats$Passed_Calls[10] <- length(which(counts>=5))
  
  stats$Passed_Calls[11] <- nrow(snv_predict)
  
  stats<-as.matrix(stats)
  parse<-as.matrix(snv_parse)

  return(list(stats_snv=stats, predicted_snv=snv_predict, parse_snv=parse))
  
  } else{
    
    print("Warning: There are no predicted SNV calls in this sample. Re-examine myresults$smurf_snv$parse_snv and set a lower snv.cutoff.")
    
    # Generate stats
    stats <- matrix(,nrow = 11, ncol = 1)
    stats <- as.data.frame(stats)
    colnames(stats) <- c("Passed_Calls")
    rownames(stats) <- c("Strelka2", "Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "Atleast4", "All5", "SMuRF_SNV")
    
    counts <- apply(snv_parse[, c("FILTER_Strelka2","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
    
    stats$Passed_Calls[1] <- length(which((snv_parse$FILTER_Strelka2==TRUE)))
    stats$Passed_Calls[2] <- length(which((snv_parse$FILTER_Mutect2==TRUE)))
    stats$Passed_Calls[3] <- length(which((snv_parse$FILTER_Freebayes==TRUE)))
    stats$Passed_Calls[4] <- length(which((snv_parse$FILTER_Vardict==TRUE)))
    stats$Passed_Calls[5] <- length(which((snv_parse$FILTER_Varscan==TRUE)))
    
    stats$Passed_Calls[6] <- length(which(counts>=1))
    stats$Passed_Calls[7] <- length(which(counts>=2))
    stats$Passed_Calls[8] <- length(which(counts>=3))
    stats$Passed_Calls[9] <- length(which(counts>=4))
    stats$Passed_Calls[10] <- length(which(counts>=5))
    
    stats$Passed_Calls[11] <- 0
    
    stats<-as.matrix(stats)
    parse<-as.matrix(snv_parse)

    return(list(stats_snv=stats, parse_snv=parse))
  }
  
}
