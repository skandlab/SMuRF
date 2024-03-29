#' indelRF-Predict
#'
#' Indel prediction model
#' Step 2 Predict
#'
#'@param parsevcf List object containing 1.snv-parse and 2.indel-parse matrices.
#'@param indel.cutoff Threshold for true predictions by SMuRF
#'@param fixed.indel.cutoff debug for when threshold is zero. 
#'
#' @examples
#' 
#' 
#' @export
indelRFpredict = function(parsevcf, indel.cutoff, fixed.indel.cutoff=F){
  
  final <- parsevcf[[2]]
  
  print("Predicting INDELs")
  
  df <- as.h2o(final)

  smurfdir <- find.package("smurf")
  
  if(exists('indel.cutoff')==F) {indel.cutoff = 'default'}
  
  # smurfmodeldir <- paste0(smurfdir, "/data/smurf-indel-model-v6") #SMuRFv1.6
  # indel_model <- h2o.loadModel(path = smurfmodeldir)
#  smurfmodeldir <- paste0(smurfdir, "/data/indel-v7.zip") #SMuRFv2.0
  smurfmodeldir <- paste0(smurfdir, "/data/indel-v3-version7.zip") #SMuRFv3.0.0
  indel_model <- h2o.import_mojo(smurfmodeldir)
  
  #Define cutoffs
  if (indel.cutoff == 'default') {
    # cutoff = h2o.find_threshold_by_max_metric(h2o.performance(indel_model), "f1") #0.4220476
    # cutoff = h2o.find_threshold_by_max_metric(h2o.performance(indel_model), "f1") #0.4794521
    #cutoff = 0.2272 #High sensitivity Recall >0.90 v7-2, version 2.0
     cutoff =  0.206274634032831 	  
    # cutoff = 0.28175877 #High sensitivity Recall >0.90
    # cutoff = 0.1875 #High sensitivity Recall >0.95
  } else if (indel.cutoff != 'default') {
    cutoff = indel.cutoff
  } else if (indel.cutoff == 0) {
    cutoff = 0
    fixed.indel.cutoff = T
  } 
  
  
  predicted <- h2o.predict(object = indel_model, newdata = df)
  p <- as.data.frame(predicted)

  indel_parse <- cbind(final, p)
  
  results <- indel_parse[which(indel_parse$TRUE.>cutoff),]
  
  if (fixed.indel.cutoff == T) {
    results <- indel_parse
  }
  

  if (dim(results)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
    
    
    table <- results
    
    names(table)[names(table) == 'TRUE.'] <- 'SMuRF_score'

    indel_predict <- table[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                            "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                            "Sample_Name",
                            "Alt_Allele_Freq",
                            "N_refDepth","N_altDepth","T_refDepth","T_altDepth",
                            "SMuRF_score")]
    
    
    # Generate stats
    stats <- matrix(,nrow = 11, ncol = 1)
    stats <- as.data.frame(stats)
    colnames(stats) <- c("Passed_Calls")
    rownames(stats) <- c("Strelka2", "Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "Atleast4", "All5", "SMuRF_INDEL")
    
    counts <- apply(indel_parse[, c("FILTER_Strelka2","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
    
    stats$Passed_Calls[1] <- length(which((indel_parse$FILTER_Strelka2==TRUE)))
    stats$Passed_Calls[2] <- length(which((indel_parse$FILTER_Mutect2==TRUE)))
    stats$Passed_Calls[3] <- length(which((indel_parse$FILTER_Freebayes==TRUE)))
    stats$Passed_Calls[4] <- length(which((indel_parse$FILTER_Vardict==TRUE)))
    stats$Passed_Calls[5] <- length(which((indel_parse$FILTER_Varscan==TRUE)))
    
    stats$Passed_Calls[6] <- length(which(counts>=1))
    stats$Passed_Calls[7] <- length(which(counts>=2))
    stats$Passed_Calls[8] <- length(which(counts>=3))
    stats$Passed_Calls[9] <- length(which(counts>=4))
    stats$Passed_Calls[10] <- length(which(counts>=5))
    
    stats$Passed_Calls[11] <- nrow(indel_predict)
    
    stats<-as.matrix(stats)
    parse<-as.matrix(indel_parse)

    return(list(stats_indel=stats, predicted_indel=indel_predict, parse_indel=parse))
   
   } else{
     
     print("Warning: There are no predicted indel calls in this sample. Re-examine myresults$smurf_indel$parse_indel and set a lower indel.cutoff.")
    
     # Generate stats
     stats <- matrix(,nrow = 11, ncol = 1)
     stats <- as.data.frame(stats)
     colnames(stats) <- c("Passed_Calls")
     rownames(stats) <- c("Strelka2", "Mutect2", "FreeBayes", "VarDict", "VarScan", "Atleast1", "Atleast2", "Atleast3", "Atleast4", "All5", "SMuRF_INDEL")
     
     counts <- apply(indel_parse[, c("FILTER_Strelka2","FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan")], 1, function(x) length(which(x=="TRUE")))
     
     stats$Passed_Calls[1] <- length(which((indel_parse$FILTER_Strelka2==TRUE)))
     stats$Passed_Calls[2] <- length(which((indel_parse$FILTER_Mutect2==TRUE)))
     stats$Passed_Calls[3] <- length(which((indel_parse$FILTER_Freebayes==TRUE)))
     stats$Passed_Calls[4] <- length(which((indel_parse$FILTER_Vardict==TRUE)))
     stats$Passed_Calls[5] <- length(which((indel_parse$FILTER_Varscan==TRUE)))
     
     stats$Passed_Calls[6] <- length(which(counts>=1))
     stats$Passed_Calls[7] <- length(which(counts>=2))
     stats$Passed_Calls[8] <- length(which(counts>=3))
     stats$Passed_Calls[9] <- length(which(counts>=4))
     stats$Passed_Calls[10] <- length(which(counts>=5))
     
     stats$Passed_Calls[11] <- 0
     
     stats<-as.matrix(stats)
     parse<-as.matrix(indel_parse)
     
     return(list(stats_indel=stats, parse_indel=parse))
   } 
  
  
}
