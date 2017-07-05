#' parse
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
parsevcfall = function(x){

  
  #ScanVCF with required parameters
  svp_m<-ScanVcfParam(info=c("BaseQRankSum","ClippingRankSum","FS","HCNT","MQ","MQRankSum","NLOD","ReadPosRankSum","TLOD"),geno=NA)
  svp_f<-ScanVcfParam(info=c("AB","ABP","EPP","EPPR","GTI","LEN","MQM","MQMR","NS","NUMALT","ODDS","PAIRED","PAIREDR","RPP","SAP","SRP"),geno=NA)
  svp_vs<-ScanVcfParam(info=c("SSC","GPV","SS","SPV"), geno=NA)
  svp_vd<-ScanVcfParam(info=c("SSF","MSI","SOR"), geno=NA)
  
  #Read Vcf of each caller and Convert collapsed vcf to VRanges object
  obj1<- readVcf(x[[1]], "hg19", svp_m)
  obj2<-as(obj1, "VRanges")
  obj3<- readVcf(x[[2]], "hg19",svp_f)
  obj4<-as(obj3, "VRanges")
  obj5<- readVcf(x[[3]], "hg19", svp_vs)
  obj6<-as(obj5, "VRanges")
  obj7<- readVcf(x[[4]], "hg19", svp_vd)
  obj8<-as(obj7, "VRanges")
   
  

   obj1cds<- readVcf(x[[5]], "hg19", svp_m)
   obj2cds<-as(obj1cds, "VRanges")
   obj3cds<- readVcf(x[[6]], "hg19",svp_f)
   obj4cds<-as(obj3cds, "VRanges")
   obj5cds<- readVcf(x[[7]], "hg19", svp_vs)
   obj6cds<-as(obj5cds, "VRanges")
   obj7cds<- readVcf(x[[8]], "hg19", svp_vd)
   obj8cds<-as(obj7cds, "VRanges")
   
  # a<-list("mutect_raw"=obj2, "mutect_vcf"=obj1, "freebayes_raw"=obj4, "freebayes_vcf"=obj3,
  #         "varscan_raw"=obj6, "varscan_vcf"=obj5, "vardict_raw"=obj8, "vardict_vcf"=obj7
  # )
  
   a<-list("mutect_raw"=obj2, "mutect_vcf"=obj1, "freebayes_raw"=obj4, "freebayes_vcf"=obj3,
           "varscan_raw"=obj6, "varscan_vcf"=obj5, "vardict_raw"=obj8, "vardict_vcf"=obj7,
           "mutect_raw_cds"=obj2cds, "mutect_vcf_cds"=obj1cds, "freebayes_raw_cds"=obj4cds, "freebayes_vcf_cds"=obj3cds,
           "varscan_raw_cds"=obj6cds, "varscan_vcf_cds"=obj5cds, "vardict_raw_cds"=obj8cds, "vardict_vcf_cds"=obj7cds
           )
   
   return(a)
   
   
}

