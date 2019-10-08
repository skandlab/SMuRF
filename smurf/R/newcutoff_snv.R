#' newcutoff_snv
#'
#' SMuRFv1.6 (R-3.3.1, R-3.5.1)
#' 
#'
#'  
#' @examples
#' 
#' 
#' @export
newcutoff_snv = function(parsevcf, snv.cutoff){
  
  
  if (getRversion()<3.5) {
      if (snv.cutoff == 'default') {
        print('Using default snv.cutoff')
        cutoff = 0.24 #High sensitivity Recall >0.95
      } else if (snv.cutoff != 'default' & snv.cutoff !=0) {
        print("Assigning new snv.cutoff")
        cutoff = snv.cutoff
      } else if (snv.cutoff == 0) {
        cutoff = 0
        # fixed.snv.cutoff = T
      } 
    }
    else {
        if (snv.cutoff == 'default') {
          print('Using default snv.cutoff')
          cutoff = 0.254 #High sensitivity Recall >0.95
        } else if (snv.cutoff != 'default' & snv.cutoff !=0) {
          print("Assigning new snv.cutoff")
          cutoff = snv.cutoff
        } else if (snv.cutoff == 0) {
          cutoff = 0
          # fixed.snv.cutoff = T
        } 
    }
    
    suppressWarnings(suppressMessages(library(dplyr)))
    
    # final <- parsevcf[[1]]
    snv_parse <- parsevcf[[1]]
    
    
    if (cutoff == 0) {
      snv_predict = snv_parse
    } else {
      snv_predict = dplyr::filter(snv_parse, TRUE.>cutoff)
    }
    
    if (dim(snv_predict)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
      
    names(snv_predict)[names(snv_predict) == 'TRUE.'] <- 'SMuRF_score'
    snv_predict <- snv_predict[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                      "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","Sample_Name",
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
  parse<-as.matrix(snv_parse)

  
  return(list(stats_snv=stats, predicted_snv=snv_predict, parse_snv=parse))
  
  } else{
    
    print("Warning: There are no predicted SNV calls in this sample. Try to set a lower snv.cutoff.")
    
    parse<-as.matrix(snv_parse)

    return(list(parse_snv=parse))
  } 
  
  
  
}
