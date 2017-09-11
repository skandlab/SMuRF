#' parsevcf
#'
#' Initial file parse and preparation before filtering
#'
#' @param filedir Sets the working directory containing the vcf.gz files from 
#' callers MuTect2, Freebayes, VarDict and VarScan.  
#' include Allele frequency here
#'   
#' @examples
#' 
#' 
#' @export
parsevcf = function(x){
  
  print("Parsing step")
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  
  #find header names
  T_m2 <- suppressWarnings(scanVcfHeader(x[[1]])@samples[grep(substrRight(scanVcfHeader(x[[1]])@samples, 2), pattern="-T")])
  T_f <- suppressWarnings(scanVcfHeader(x[[2]])@samples[grep(substrRight(scanVcfHeader(x[[2]])@samples, 2), pattern="-T")])
  T_vs <- suppressWarnings(scanVcfHeader(x[[3]])@samples[grep(substrRight(scanVcfHeader(x[[3]])@samples, 2), pattern="-T")])
  T_vd <- suppressWarnings(scanVcfHeader(x[[4]])@samples[grep(substrRight(scanVcfHeader(x[[4]])@samples, 2), pattern="-T")])
  
  # T_svp_m<-suppressWarnings(ScanVcfParam(info=c("ANN"),samples=T_m2))    #AF=AO/DP(totalDepth)
  # T_svp_f<-suppressWarnings(ScanVcfParam(info=c("ANN"),samples=T_f))    #AF=AO/DP(totalDepth)
  # T_svp_vs<-suppressWarnings(ScanVcfParam(info=c("ANN"),samples=T_vs))    #AF=AO/DP(totalDepth)
  # T_svp_vd<-suppressWarnings(ScanVcfParam(info=c("ANN"),samples=T_vd))    #AF=AO/DP(totalDepth)
  # 
  # #ScanVCF with required parameters
  # svp_m<-suppressWarnings(ScanVcfParam(info=c("MQ","MQRankSum","TLOD","NLOD"),geno=NA))   
  # svp_f<-suppressWarnings(ScanVcfParam(info=c("MQM","MQMR","LEN"),geno=NA))             
  # svp_vs<-suppressWarnings(ScanVcfParam(info=c("SSC","SPV"),geno=NA))       
  # svp_vd<-suppressWarnings(ScanVcfParam(info=c("SSF","SOR","MSI"),geno=NA))
  
  #ScanVCF with required parameters
  svp_m<-suppressWarnings(ScanVcfParam(info=c("MQ","MQRankSum","TLOD","NLOD"),samples=T_m2))   
  svp_f<-suppressWarnings(ScanVcfParam(info=c("MQM","MQMR","LEN"),samples=T_f))             
  svp_vs<-suppressWarnings(ScanVcfParam(info=c("SSC","SPV"),samples=T_vs))       
  svp_vd<-suppressWarnings(ScanVcfParam(info=c("SSF","SOR","MSI"),samples=T_vd))
  
  
  N_m2 <- suppressWarnings(scanVcfHeader(x[[1]])@samples[grep(substrRight(scanVcfHeader(x[[1]])@samples, 2), pattern="-N")])
  N_f <- suppressWarnings(scanVcfHeader(x[[2]])@samples[grep(substrRight(scanVcfHeader(x[[2]])@samples, 2), pattern="-N")])
  N_vs <- suppressWarnings(scanVcfHeader(x[[3]])@samples[grep(substrRight(scanVcfHeader(x[[3]])@samples, 2), pattern="-N")])
  N_vd <- suppressWarnings(scanVcfHeader(x[[4]])@samples[grep(substrRight(scanVcfHeader(x[[4]])@samples, 2), pattern="-N")])
  
  #ScanVCF with required parameters
  
  N_svp_m<-suppressWarnings(ScanVcfParam(info=NA,samples=N_m2))    #AF=AO/DP(totalDepth)
  N_svp_f<-suppressWarnings(ScanVcfParam(info=NA,samples=N_f))    #AF=AO/DP(totalDepth)
  N_svp_vs<-suppressWarnings(ScanVcfParam(info=NA,samples=N_vs))    #AF=AO/DP(totalDepth)
  N_svp_vd<-suppressWarnings(ScanVcfParam(info=NA,samples=N_vd))    #AF=AO/DP(totalDepth)
  
  #Read Vcf of each caller and Convert collapsed vcf to VRanges object
  obj1<- suppressWarnings(readVcf(x[[1]], "hg19", svp_m))
  obj2<-suppressWarnings(as(obj1, "VRanges"))
  obj3<- suppressWarnings(readVcf(x[[2]], "hg19",svp_f))
  obj4<-suppressWarnings(as(obj3, "VRanges"))
  obj5<- suppressWarnings(readVcf(x[[3]], "hg19", svp_vs))
  obj6<-suppressWarnings(as(obj5, "VRanges"))
  obj7<- suppressWarnings(readVcf(x[[4]], "hg19", svp_vd))
  obj8<-suppressWarnings(as(obj7, "VRanges"))
  
  #Read Vcf of each caller and Convert collapsed vcf to VRanges object
  obj1n<- suppressWarnings(readVcf(x[[1]], "hg19", N_svp_m))
  obj2n<-suppressWarnings(as(obj1n, "VRanges"))
  obj3n<- suppressWarnings(readVcf(x[[2]], "hg19",N_svp_f))
  obj4n<-suppressWarnings(as(obj3n, "VRanges"))
  obj5n<- suppressWarnings(readVcf(x[[3]], "hg19", N_svp_vs))
  obj6n<-suppressWarnings(as(obj5n, "VRanges"))
  obj7n<- suppressWarnings(readVcf(x[[4]], "hg19", N_svp_vd))
  obj8n<-suppressWarnings(as(obj7n, "VRanges"))
  
  a<-list(obj2, obj1, obj4, obj3, obj6, obj5, obj8, obj7, obj2n, obj4n, obj6n, obj8n)
  

  # a<-list("mutect_raw"=obj2, "mutect_vcf"=obj1, "freebayes_raw"=obj4, "freebayes_vcf"=obj3,
  #          "varscan_raw"=obj6, "varscan_vcf"=obj5, "vardict_raw"=obj8, "vardict_vcf"=obj7)
  return(a)
   
   
   
}


