#' newcutoffs
#'
#' 
#'
#' @param parse.dir Directory to parse files 
#' @param output.dir your output directory
#' @param snv.cutoff 
#' @param indel.cutoff
#'   
#' @examples
#' 
#' 
#' @export
#reassign cutoffs only using existing parse files

newcutoffs = function(parse.dir=NULL, snv.cutoff='default', indel.cutoff='default'){
  
  
  if (!is.null(snv.cutoff) & snv.cutoff == 'default') {
    print('Using default snv.cutoff')
    # snv.cutoff = 0.4627947 #default v7.2
    snv.cutoff = 0.0822 #High sensitivity Recall >0.99 v7-2
    # snv.cutoff = 0.4587242 #default v7
  } else if (snv.cutoff != 'default' & snv.cutoff !=0) {
    print("Assigning new snv.cutoff")
    snv.cutoff = snv.cutoff
  } else if (snv.cutoff == 0) {
    snv.cutoff = 0
  } 
  
  #Define cutoffs
  if (!is.null(indel.cutoff) & indel.cutoff == 'default') {
    # indel.cutoff = 0.4220476 #default v7.2
    indel.cutoff = 0.2272 #High sensitivity Recall >0.90 v7-2
    # indel.cutoff = 0.28175877 #default v7
  } else if (indel.cutoff != 'default') {
    indel.cutoff = indel.cutoff
  } else if (indel.cutoff == 0) {
    indel.cutoff = 0
  } 
  
  suppressWarnings(suppressMessages(library(dplyr)))
  
  if (!is.null(parse.dir)) {
    if(dir.exists(parse.dir)==F){
      stop('parse.dir path does not exist. Assigning newcutoff requires path to parse.txt files')
    } else {
      print('Checking parse directory')
      snv.parse.dir <- Sys.glob(paste0(parse.dir,'/',"snv-parse*"))
      indel.parse.dir <- Sys.glob(paste0(parse.dir,'/',"indel-parse*"))
    }
  }

      
  #SNV parse    
      if(length(snv.parse.dir)==1) {
        
        print(paste0(snv.parse.dir, ' found'))
        parse_snv = read.table(snv.parse.dir, header = T, stringsAsFactors = F, na.strings='.')
        
        snv_parse <- parse_snv
        
        if (snv.cutoff == 0) {
          snv_predict = snv_parse
        } else {
          snv_predict = dplyr::filter(snv_parse, TRUE.>snv.cutoff)
        }
        
        if (dim(snv_predict)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
          
          names(snv_predict)[names(snv_predict) == 'TRUE.'] <- 'SMuRF_score'
          snv_predict <- snv_predict[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                        "FILTER_Strelka2",    
                                        "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","Sample_Name",
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
          # parse_snv<-as.matrix(snv_parse)
          # predicted=snv_predict
          
          smurf_snv = list(stats_snv=stats,predicted_snv=snv_predict)
          
        } else{
          
          print("Warning: There are no predicted SNV calls in this sample. Try to set a lower snv.cutoff.")
          
          # parse_snv<-as.matrix(snv_parse)
          stats=NULL
          snv_predict=NULL
          
          smurf_snv = list(stats_snv=stats,predicted_snv=snv_predict)
          
        } 
        
        } else {
           print('snv-parse.txt not found. Skipping snv.')
        }
      
        
        
   #indel parse     
      if(length(indel.parse.dir)==1) {
        
        print(paste0(indel.parse.dir, ' found'))
        parse_indel = read.table(indel.parse.dir, header = T, stringsAsFactors = F, na.strings='.')
        
        indel_parse <- parse_indel
        
        if (indel.cutoff == 0) {
          indel_predict = indel_parse
        } else {
          indel_predict = dplyr::filter(indel_parse, TRUE.>indel.cutoff)
        }
        
        if (dim(indel_predict)[1] != 0) { #encountering zero predictions will exit code, output contains parse and raw file only. 
          
          names(indel_predict)[names(indel_predict) == 'TRUE.'] <- 'SMuRF_score'
          
          indel_predict <- indel_predict[,c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs",
                                    "FILTER_Strelka2",
                                    "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","Sample_Name",
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
          # parse_indel<-as.matrix(indel_parse)
          # predicted_indel=indel_predict
          
          smurf_indel = list(stats_indel=stats,predicted_indel=indel_predict)
          
          
        } else{
          print("Warning: There are no predicted indel calls in this sample. Try to set a lower indel.cutoff.")
          
          # parse<-as.matrix(indel_parse)
          stats=NULL
          indel_predict=NULL
          smurf_indel = list(stats_indel=stats,predicted_indel=indel_predict)
          
        } 
        
      } else {
        print('indel-parse.txt not found. Skipping indel.')
      }
  
  return(list(smurf_snv=smurf_snv,smurf_indel=smurf_indel))

  
  
}
