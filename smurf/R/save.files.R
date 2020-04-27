#' save.files
#'
#' 
#'
#' @param myresults List object for SMuRF output 
#' @param output.dir your output directory
#'   
#' @examples
#' 
#' 
#' @export
save.files = function (myresults, output.dir) {
  
  a<- myresults$smurf_indel$stats_indel
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/indel-stats.txt"), sep = "\t", quote = FALSE, row.names = TRUE, na = ".")
  }
  
  a<- myresults$smurf_indel$predicted_indel
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/indel-predicted.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  a<- myresults$smurf_indel$parse_indel
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/indel-parse.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  a<- myresults$smurf_snv$stats_snv
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/snv-stats.txt"), sep = "\t", quote = FALSE, row.names = TRUE, na = ".")
  }
  
  a<- myresults$smurf_snv$predicted_snv
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/snv-predicted.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  a<- myresults$smurf_snv$parse_snv
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/snv-parse.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  a<- myresults$smurf_snv_annotation$annotated
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/snv-annotated.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  a<- myresults$smurf_indel_annotation$annotated
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/indel-annotated.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  # a<- myresults$smurf_snv_roi_annotation$annotated
  # if(!is.null(a)){
  #   write.table(a , file = paste0(output.dir, "/snv-roi.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  # }
  # 
  # a<- myresults$smurf_indel_roi_annotation$annotated
  # if(!is.null(a)){
  #   write.table(a , file = paste0(output.dir, "/indel-roi.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  # }
  
  a<- myresults$parsevcf_featureselection$snv
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/snv-parse.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }

  a<- myresults$parsevcf_featureselection$indel
  if(!is.null(a)){
    write.table(a , file = paste0(output.dir, "/indel-parse.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  }
  
  #show time taken for run
  a<- myresults$time.taken
  if(!is.null(a)){
    write(a, file = paste0(output.dir, "/time.txt"))
  }
  
}


