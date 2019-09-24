#' newcutoff_indel
#'
#' SMuRFv1.6 (R-3.3.1, R-3.5.1)
#' 
#'
#'  
#' @examples
#' 
#' 
#' @export
newcutoff_indel = function(parsevcf, indel.cutoff){
  
  
  if (getRversion()<3.5) {
      if (indel.cutoff == 'default') {
        cutoff = 0.352 #High sensitivity Recall >0.80
      } else if (indel.cutoff != 'default' & indel.cutoff !=0) {
        print("Assigning new indel.cutoff")
        cutoff = indel.cutoff
      } else if (indel.cutoff == 0) {
        cutoff = 0
        fixed.indel.cutoff = T
      } 
  }
    else {
        if (indel.cutoff == 'default') {
          print('Using default indel.cutoff')
          cutoff = 0.352 #High sensitivity Recall >0.80
        } else if (indel.cutoff != 'default'& indel.cutoff !=0) {
          print("Assigning new indel.cutoff")
          cutoff = indel.cutoff
        } else if (indel.cutoff == 0) {
          cutoff = 0
          fixed.indel.cutoff = T
        } 
      }

    suppressWarnings(suppressMessages(library(dplyr)))
    
    # final <- parsevcf[[2]]
    indel_parse <- parsevcf[[2]]
    
    indel_predict = dplyr::filter(indel_parse, TRUE.>cutoff)
    
    if (fixed.indel.cutoff == T) {
      indel_predict = indel_parse
    }
    
    if (dim(indel_predict)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
      
    names(indel_predict)[names(indel_predict) == 'TRUE.'] <- 'SMuRF_score'
    indel_predict <- indel_predict[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                          "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","Sample_Name",
                                          "Alt_Allele_Freq",
                                          "N_refDepth","N_altDepth","T_refDepth","T_altDepth",
                                          "SMuRF_score")]

    
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
   parse<-as.matrix(indel_parse)

   
   return(list(stats_indel=stats, predicted_indel=indel_predict, parse_indel=parse))
   
   } else{
    print("Warning: There are no predicted indel calls in this sample. Try to set a lower indel.cutoff.")
    
     parse<-as.matrix(indel_parse)
     
     return(list(parse_indel=parse))
   } 
  
  
}
