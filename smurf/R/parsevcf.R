#' parsevcf
#'
#' Initial file parse and preparation before filtering
#'
#' @param filedir Sets the working directory containing the vcf.gz files from 
#' callers MuTect2, Freebayes, VarDict and VarScan.  
#' 
#'   
#' @examples
#' 
#' 
#' @export
parsevcf = function(x){

  #ScanVCF with required parameters
  svp_m<-ScanVcfParam(info=c("MQ","MQRankSum","TLOD","NLOD"),geno=NA)
  svp_f<-ScanVcfParam(info=c("MQM","MQMR","LEN"),geno=NA)
  svp_vs<-ScanVcfParam(info=c("SSC","SPV"), geno=NA)
  svp_vd<-ScanVcfParam(info=c("SSF","SOR","MSI"), geno=NA)
  
  #Read Vcf of each caller and Convert collapsed vcf to VRanges object
  obj1<- suppressWarnings(readVcf(x[[1]], "hg19", svp_m))
  obj2<-suppressWarnings(as(obj1, "VRanges"))
  obj3<- suppressWarnings(readVcf(x[[2]], "hg19",svp_f))
  obj4<-suppressWarnings(as(obj3, "VRanges"))
  obj5<- suppressWarnings(readVcf(x[[3]], "hg19", svp_vs))
  obj6<-suppressWarnings(as(obj5, "VRanges"))
  obj7<- suppressWarnings(readVcf(x[[4]], "hg19", svp_vd))
  obj8<-suppressWarnings(as(obj7, "VRanges"))
   
   a<-list("mutect_raw"=obj2, "mutect_vcf"=obj1, "freebayes_raw"=obj4, "freebayes_vcf"=obj3,
           "varscan_raw"=obj6, "varscan_vcf"=obj5, "vardict_raw"=obj8, "vardict_vcf"=obj7)
   return(a)
}


