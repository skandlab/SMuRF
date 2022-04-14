#' parsevcf
#'
#' Step 1 Initial parse and filtering
#'
#' @param x List object containing the four vcf.gz files from 
#' callers Strelka2, MuTect2, Freebayes, VarDict and VarScan.  
#' 
#' @param t.label Default='-T'. Identify your tumour label based on your default vcf tumour sample.
#'   
#' @examples
#' 
#' 
#' @export
parsevcf_allfeaturesall = function(x, tbi, t.label=NULL){
  
  print("Parsing step")
  
  # if(exists('t.label')==F) {t.label = NULL}
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  # substrLeft <- function(x, n){
  #   substr(x, 1, nchar(x)-n)
  # }
  
  print("reading vcfs")
  start.time=Sys.time()

  #ScanVCF with required parameters
  svp_m<-ScanVcfParam(info=c("FS","MQ","MQRankSum","NLOD","ReadPosRankSum","TLOD"), samples=suppressWarnings(scanVcfHeader(x[[1]])@samples), geno=c("AD", "AF", "DP"))
  svp_f<-ScanVcfParam(info=c("GTI","LEN","MQM","MQMR","ODDS"),samples=suppressWarnings(scanVcfHeader(x[[2]])@samples), geno=c("RO","DP"))
  svp_vs<-ScanVcfParam(info=c("SSC","GPV","SS","SPV"), samples=suppressWarnings(scanVcfHeader(x[[3]])@samples), geno=c("AD","FREQ", "DP","RD"))
  svp_vd<-ScanVcfParam(info=c("SSF","MSI","SOR"), samples=suppressWarnings(scanVcfHeader(x[[4]])@samples), geno=c("AD", "AF", "DP"))
  svp_s = ScanVcfParam(info=c("QSS","MQ","SomaticEVS","ReadPosRankSum"), samples=suppressWarnings(scanVcfHeader(x[[5]])@samples), geno=c("AD", "AF", "DP"))
  
  
  ## read only chromosomes 1-22, X, Y and M
  print("reading mutect2")
  vcf_m2<- suppressWarnings(readVcf(tbi[[1]], genome=seqinfo(scanVcfHeader(x[[1]])), svp_m))
  print("reading freebayes")
  vcf_f<- suppressWarnings(readVcf(tbi[[2]], genome=seqinfo(scanVcfHeader(x[[2]])), svp_f))
  print("reading varscan")
  vcf_vs<- suppressWarnings(readVcf(tbi[[3]], genome=seqinfo(scanVcfHeader(x[[3]])), svp_vs))
  print("reading vardict")
  vcf_vd<- suppressWarnings(readVcf(tbi[[4]], genome=seqinfo(scanVcfHeader(x[[4]])), svp_vd))
  print("reading strelka2")
  vcf_s2<- suppressWarnings(readVcf(tbi[[5]], genome=seqinfo(scanVcfHeader(x[[5]])), svp_s))
  
  
  # remove duplicates from vcfs
  vcf_m2=vcf_m2[!duplicated(rowRanges(vcf_m2))]
  vcf_f=vcf_f[!duplicated(rowRanges(vcf_f))]
  vcf_vs=vcf_vs[!duplicated(rowRanges(vcf_vs))]
  vcf_vd=vcf_vd[!duplicated(rowRanges(vcf_vd))]
  vcf_s2=vcf_s2[!duplicated(rowRanges(vcf_s2))]
  end.time=Sys.time() 
  
  print(end.time-start.time)
  
  # get vcf headers
  H_m2=header(vcf_m2)
  H_f=header(vcf_f)  
  H_vs=header(vcf_vs) 
  H_vd=header(vcf_vd)
  H_s2=header(vcf_s2)
  
  # sample name
  if (is.null(t.label)) {
    sampleid.t = H_m2@header@listData$PEDIGREE$Derived
    sampleid.n = H_m2@header@listData$PEDIGREE$Original
  } else {
    sampleid.t <-H_m2@samples[grep((H_m2@samples), pattern=t.label)]
    sampleid.n <-H_m2@samples[grep((H_m2@samples), pattern=t.label, invert=T)]
  }
  
  if(length(sampleid.t)!=1 & length(sampleid.n)!=1) {
    write(paste0("t.label='",t.label,"'"),stdout())
    stop('t.label for tumor sample is not unique, duplicated or missing')
  }
  
  print("extracting calls passed by at least 1 caller")
  start.time=Sys.time()
  # get rows that are SNVs
  snv_m2<-isSNV(vcf_m2, singleAltOnly=FALSE)
  snv_f<-isSNV(vcf_f, singleAltOnly=FALSE)
  snv_vs<-isSNV(vcf_vs, singleAltOnly=FALSE)
  snv_vd<-isSNV(vcf_vd, singleAltOnly=FALSE)  
  snv_s2<-isSNV(vcf_s2, singleAltOnly=FALSE)
  

    # extract passed calls from each caller
    pass_m2<- vcf_m2@fixed$FILTER=="PASS" | vcf_m2@fixed$FILTER=="MinAF"
    pass_f<- vcf_f@fixed$FILTER=="PASS" 
    pass_vs<- vcf_vs@fixed$FILTER=="PASS"
    # pass_vs<- vcf_vs@fixed$FILTER=="PASS" | vcf_vs@fixed$FILTER=="SpvFreq" | vcf_vs@fixed$FILTER=="REJECT;SpvFreq" #SpvFreq not found in GATK4
    pass_vd<- vcf_vd@fixed$FILTER=="PASS"
    pass_s2<- vcf_s2@fixed$FILTER=="PASS" | vcf_s2@fixed$FILTER=="MinAF"

  
  # get passed snv calls from all callers
  snv_pass=c(rowRanges(vcf_m2[pass_m2&snv_m2]), 
             rowRanges(vcf_f[pass_f&snv_f]),  
             rowRanges(vcf_vs[pass_vs&snv_vs]), 
             rowRanges(vcf_vd[pass_vd&snv_vd]),
             rowRanges(vcf_s2[pass_s2&snv_s2]),
             ignore.mcols=T)
  snv_pass= unique(snv_pass)
  
  ###### parsing SNVs
  # get the row numbers of calls that are passed by at least 1 caller in each vcf
  snv.m2.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_m2), type="equal"))
  snv.f.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_f), type="equal"))
  snv.vs.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vs), type="equal"))
  snv.vd.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vd), type="equal"))
  snv.s2.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_s2), type="equal"))
  
  end.time=Sys.time() 
  print(end.time-start.time)
  
  print("extracting meta data from VRanges")
  start.time=Sys.time() 
  
  if (length(snv.m2.index)!=0) {
  m2.na=F
    
  # convert vcf to VRanges then to Granges, keep metadata columns
  vr_m2<- as(vcf_m2[snv.m2.index], "VRanges")
  vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
  gr_m2=GRanges(vr_m2[[sampleid.t]])
  mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))

  gr_m2 <- unique(gr_m2)
  gr_m2$FILTER=vcf_m2[snv.m2.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from MuTect2.")
    m2.na=T
  }
  
  
  if (length(snv.f.index)!=0) {
  f.na=F  
  vr_f<- suppressWarnings(as(vcf_f[snv.f.index], "VRanges"))
  vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
  gr_f=GRanges(vr_f[[sampleid.t]])
  mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                            T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                            N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
  
  gr_f <- unique(gr_f)
  gr_f$FILTER=vcf_f[snv.f.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from FreeBayes.")
    f.na=T
  }
  
  
  if (length(snv.vs.index)!=0) {
  vs.na=F  
  vr_vs<- suppressWarnings(as(vcf_vs[snv.vs.index], "VRanges"))
  vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
  gr_vs=GRanges(vr_vs[[sampleid.t]])
  mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vs <- unique(gr_vs)
  gr_vs$FILTER=vcf_vs[snv.vs.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from VarScan.")
    vs.na=T
  }
  
  if (length(snv.vd.index)!=0) {
  vd.na=F  
  vr_vd<- suppressWarnings(as(vcf_vd[snv.vd.index], "VRanges"))
  vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
  gr_vd=GRanges(vr_vd[[sampleid.t]])
  mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vd <- unique(gr_vd)
  gr_vd$FILTER=vcf_vd[snv.vd.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from VarDict.")
    vd.na=T
  }
  
  
  
  if (length(snv.s2.index)!=0) {
    s2.na=F
    
    # convert vcf to VRanges then to Granges, keep metadata columns
    vr_s2<- as(vcf_s2[snv.s2.index], "VRanges")
    vr_s2=GenomicRanges::split(vr_s2, vr_s2@sampleNames)
    gr_s2=GRanges(vr_s2[[sampleid.t]])
    mcols(gr_s2)=cbind(mcols(gr_s2), data.frame(REF=ref(vr_s2[[sampleid.t]]), ALT=alt(vr_s2[[sampleid.t]]), T_totalDepth=totalDepth(vr_s2[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_s2[[sampleid.t]]), T_altDepth=altDepth(vr_s2[[sampleid.t]]),N_totalDepth=totalDepth(vr_s2[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_s2[[sampleid.n]]), N_altDepth=altDepth(vr_s2[[sampleid.n]]), stringsAsFactors=F))
    
    gr_s2 <- unique(gr_s2)
    gr_s2$FILTER=vcf_s2[snv.s2.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed SNVs from Strelka2.")
    s2.na=T
  }
  
  
  end.time=Sys.time() 
  print(end.time-start.time)
  
  print("merge 5 vcfs and meta data")
  start.time=Sys.time()
  
  ## merge 5 vcfs and meta data
  m2.cols = c(names(vcf_m2@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  f.cols = c(names(vcf_f@info@listData),'RO','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  vs.cols = c(names(vcf_vs@info@listData),'FREQ','RD','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  vd.cols = c(names(vcf_vd@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  s2.cols = c(names(vcf_s2@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  
  names_m2= paste(m2.cols, "_Mutect2", sep="")
  names_f= paste(f.cols, "_Freebayes", sep="")
  names_vs= paste(vs.cols, "_Varscan", sep="")
  names_vd= paste(vd.cols, "_Vardict", sep="")
  names_s2= paste(s2.cols, "_Strelka2", sep="")
  
  # names_m2= paste(names(mcols(gr_m2)[,-1]), "_Mutect2", sep="")
  # names_f= paste(names(mcols(gr_f)[,-1]), "_Freebayes", sep="")
  # names_vs= paste(names(mcols(gr_vs)[,-1]), "_Varscan", sep="")
  # names_vd= paste(names(mcols(gr_vd)[,-1]), "_Vardict", sep="")
  # names_s2= paste(names(mcols(gr_s2)[,-1]), "_Strelka2", sep="")
  
  meta_data=data.frame(snv_pass)[,1:3]
  meta_data[c(names_m2,names_f,names_vs,names_vd,names_s2)]=NA
  
  #do not merge when index=0, caller.na=T
  if (m2.na==F) {
    meta_data[Biostrings::match(gr_m2, snv_pass), names_m2]=data.frame(mcols(gr_m2)[,m2.cols])
  }
  if (f.na==F) {
    meta_data[Biostrings::match(gr_f, snv_pass), names_f]=data.frame(mcols(gr_f)[,f.cols])
  }
  if (vs.na==F) {
    meta_data[Biostrings::match(gr_vs, snv_pass), names_vs]=data.frame(mcols(gr_vs)[,vs.cols])
  }
  if (vd.na==F) {
    meta_data[Biostrings::match(gr_vd, snv_pass), names_vd]=data.frame(mcols(gr_vd)[,vd.cols])
  }
  if (s2.na==F) {
    meta_data[Biostrings::match(gr_s2, snv_pass), names_s2]=data.frame(mcols(gr_s2)[,s2.cols])
  }
  
  # #do not merge when the 1st row is all NAs + last column PASSED/REJECT
  # if ((rowSums(is.na(data.frame(mcols(gr_m2)[1,1:length(mcols(gr_m2))-1])))!=ncol(data.frame(mcols(gr_m2)[1,1:length(mcols(gr_m2))-1])))==TRUE) {
  #   meta_data[Biostrings::match(gr_m2, snv_pass), names_m2]=data.frame(mcols(gr_m2)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_f)[1,1:length(mcols(gr_f))-1])))!=ncol(data.frame(mcols(gr_f)[1,1:length(mcols(gr_f))-1])))==TRUE) {
  #   meta_data[Biostrings::match(gr_f, snv_pass), names_f]=data.frame(mcols(gr_f)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_vs)[1,1:length(mcols(gr_vs))-1])))!=ncol(data.frame(mcols(gr_vs)[1,1:length(mcols(gr_vs))-1])))==TRUE) {
  #   meta_data[Biostrings::match(gr_vs, snv_pass), names_vs]=data.frame(mcols(gr_vs)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_vd)[1,1:length(mcols(gr_vd))-1])))!=ncol(data.frame(mcols(gr_vd)[1,1:length(mcols(gr_vd))-1])))==TRUE) {
  #   meta_data[Biostrings::match(gr_vd, snv_pass), names_vd]=data.frame(mcols(gr_vd)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_s2)[1,1:length(mcols(gr_s2))-1])))!=ncol(data.frame(mcols(gr_s2)[1,1:length(mcols(gr_s2))-1])))==TRUE) {
  #   meta_data[Biostrings::match(gr_s2, snv_pass), names_s2]=data.frame(mcols(gr_s2)[,-1])
  # }

  end.time=Sys.time()
  print(end.time-start.time)
  
  
  print("formating")
  start.time=Sys.time()
  
  #pick out cases with some indels
  meta_indel_cases <- subset(meta_data, 
                               nchar(REF_Mutect2)!=1 |
                               nchar(REF_Freebayes)!=1 |
                               nchar(REF_Varscan)!=1 |
                               nchar(REF_Vardict)!=1 |
                               nchar(REF_Strelka2)!=1 |
                               nchar(ALT_Mutect2)!=1 |
                               nchar(ALT_Freebayes)!=1 |
                               nchar(ALT_Varscan)!=1 |
                               nchar(ALT_Vardict)!=1 |
                               nchar(ALT_Strelka2)!=1)
  
  # extract reference allele
  ref=meta_data[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict", "REF_Strelka2")] 
  ref.ind=which(!is.na(ref), arr.ind=T) # find non-NA alleles
  ref.ind=ref.ind[order(ref.ind[,1]),] #order array indices by row number
  ref.ind=ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
  meta_data$REF=ref[ref.ind]
  
  # extract alternate allele
  alt=meta_data[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict", "ALT_Strelka2")]
  suppressWarnings(alt$ALT_Mutect2[(!is.na(alt$ALT_Mutect2) & nchar(alt$ALT_Mutect2)!=1)] <- substrRight(alt$ALT_Mutect2[(!is.na(alt$ALT_Mutect2) & nchar(alt$ALT_Mutect2)!=1)],1))
  suppressWarnings(alt$ALT_Freebayes[(!is.na(alt$ALT_Freebayes) & nchar(alt$ALT_Freebayes)!=1)] <- substrRight(alt$ALT_Freebayes[(!is.na(alt$ALT_Freebayes) & nchar(alt$ALT_Freebayes)!=1)],1))
  suppressWarnings(alt$ALT_Varscan[(!is.na(alt$ALT_Varscan) & nchar(alt$ALT_Varscan)!=1)] <- substrRight(alt$ALT_Varscan[(!is.na(alt$ALT_Varscan) & nchar(alt$ALT_Varscan)!=1)],1))
  suppressWarnings(alt$ALT_Vardict[(!is.na(alt$ALT_Vardict) & nchar(alt$ALT_Vardict)!=1)] <- substrRight(alt$ALT_Vardict[(!is.na(alt$ALT_Vardict) & nchar(alt$ALT_Vardict)!=1)],1))
  suppressWarnings(alt$ALT_Strelka2[(!is.na(alt$ALT_Strelka2) & nchar(alt$ALT_Strelka2)!=1)] <- substrRight(alt$ALT_Strelka2[(!is.na(alt$ALT_Strelka2) & nchar(alt$ALT_Strelka2)!=1)],1))
  # alt$ALT_Freebayes[nchar(alt$ALT_Freebayes)!=1] <- NA
  # alt$ALT_Varscan[nchar(alt$ALT_Varscan)!=1] <- NA
  # alt$ALT_Vardict[nchar(alt$ALT_Vardict)!=1] <- NA
  alt.ind=which(!is.na(alt), arr.ind=T) #find all non-NA alleles
  alt.ind=alt.ind[order(alt.ind[,1]),]  #order array indices by row number
  alt.ind=alt.ind[!duplicated(alt.ind[,1]),]# take the first non-NA allele of each row
  meta_data$ALT=alt[alt.ind]  
  
  # check for list objects in columns
  for (i in 1:ncol(meta_data)) {
    if(class(meta_data[,i])=='list'){
      #meta_data[,i][sapply(meta_data[,i], is.null)] <- NA
      meta_data[,i] = unlist(meta_data[,i])
    }
  }
  
  # calculate alt allele frequency
  meta_data$T_altDepth_Freebayes <- meta_data$T_totalDepth_Freebayes-meta_data$RO_Freebayes
  meta_data$T_refDepth_Freebayes <- meta_data$RO_Freebayes
  meta_data$AF_Freebayes <- meta_data$T_altDepth_Freebayes/meta_data$T_totalDepth_Freebayes
  # meta_data$AF_Freebayes[meta_data$AF_Freebayes == 0] <- NA
  
  meta_data$T_refDepth_Varscan <- meta_data$RD_Varscan
  meta_data$T_altDepth_Varscan <- meta_data$T_totalDepth_Varscan-meta_data$RD_Varscan
  
  meta_data$AF_Strelka2[sapply(meta_data$AF_Strelka2, is.null)] <- NA
  meta_data$AF_Strelka2 <- unlist(meta_data$AF_Strelka2)
  meta_data$T_altDepth_Strelka2 <- round(meta_data$AF_Strelka2*meta_data$T_totalDepth_Strelka2)
  meta_data$T_refDepth_Strelka2 <- meta_data$T_totalDepth_Strelka2-meta_data$T_altDepth_Strelka2
  

  #catch vcf caller mapping errors
  if(all(is.na(meta_data$T_altDepth_Mutect2))==T | 
     all(is.na(meta_data$T_refDepth_Mutect2))==T |
     all(is.na(meta_data$N_altDepth_Mutect2))==T | 
     all(is.na(meta_data$N_refDepth_Mutect2))==T) {
    meta_data$T_refDepth_Mutect2 <- NA
    meta_data$N_refDepth_Mutect2 <- NA
    meta_data$T_altDepth_Mutect2 <- NA
    meta_data$N_altDepth_Mutect2 <- NA
  }
  if(all(is.na(meta_data$T_altDepth_Freebayes))==T | 
     all(is.na(meta_data$T_refDepth_Freebayes))==T 
     # all(is.na(meta_data$N_altDepth_Freebayes))==T | 
     # all(is.na(meta_data$N_refDepth_Freebayes))==T
     ) {
    meta_data$T_refDepth_Freebayes <- NA
    # meta_data$N_refDepth_Freebayes <- NA
    meta_data$T_altDepth_Freebayes <- NA
    # meta_data$N_altDepth_Freebayes <- NA
  }
  if(all(is.na(meta_data$T_altDepth_Varscan))==T | 
     all(is.na(meta_data$T_refDepth_Varscan))==T 
     # all(is.na(meta_data$N_altDepth_Varscan))==T | 
     # all(is.na(meta_data$N_refDepth_Varscan))==T
     ) {
    meta_data$T_refDepth_Varscan <- NA
    # meta_data$N_refDepth_Varscan <- NA
    meta_data$T_altDepth_Varscan <- NA
    # meta_data$N_altDepth_Varscan <- NA
  }
  if(all(is.na(meta_data$T_altDepth_Vardict))==T | 
     all(is.na(meta_data$T_refDepth_Vardict))==T |
     all(is.na(meta_data$N_altDepth_Vardict))==T | 
     all(is.na(meta_data$N_refDepth_Vardict))==T) {
    meta_data$T_refDepth_Vardict <- NA
    meta_data$N_refDepth_Vardict <- NA
    meta_data$T_altDepth_Vardict <- NA
    meta_data$N_altDepth_Vardict <- NA
  }
  if(all(is.na(meta_data$T_altDepth_Strelka2))==T | 
     all(is.na(meta_data$T_refDepth_Strelka2))==T 
     # all(is.na(meta_data$N_altDepth_Strelka2))==T | 
     # all(is.na(meta_data$N_refDepth_Strelka2))==T
  ) {
    meta_data$T_refDepth_Strelka2 <- NA
    # meta_data$N_refDepth_Strelka2 <- NA
    meta_data$T_altDepth_Strelka2 <- NA
    # meta_data$N_altDepth_Strelka2 <- NA
  }
  
  
  # calculate mean depth
  meta_data$N_refDepth <- round(rowMeans(meta_data[, c("N_refDepth_Mutect2", "N_refDepth_Freebayes", "N_refDepth_Varscan", "N_refDepth_Vardict","N_refDepth_Strelka2")],na.rm = TRUE))
  meta_data$N_altDepth <- round(rowMeans(meta_data[, c("N_altDepth_Mutect2", "N_altDepth_Freebayes", "N_altDepth_Varscan", "N_altDepth_Vardict","N_altDepth_Strelka2")],na.rm = TRUE))
  meta_data$T_refDepth <- round(rowMeans(meta_data[, c("T_refDepth_Mutect2", "T_refDepth_Freebayes", "T_refDepth_Varscan", "T_refDepth_Vardict","T_refDepth_Strelka2")],na.rm = TRUE))
  meta_data$T_altDepth <- round(rowMeans(meta_data[, c("T_altDepth_Mutect2", "T_altDepth_Freebayes", "T_altDepth_Varscan", "T_altDepth_Vardict","T_refDepth_Strelka2")],na.rm = TRUE))
  
  meta_data$N_refDepth[is.nan(meta_data$N_refDepth)] <- 0
  meta_data$N_altDepth[is.nan(meta_data$N_altDepth)] <- 0
  meta_data$T_refDepth[is.nan(meta_data$T_refDepth)] <- 0
  meta_data$T_altDepth[is.nan(meta_data$T_altDepth)] <- 0
  
  #re-calculate mean AF
  meta_data$AF = round(meta_data$T_altDepth/(meta_data$T_altDepth+meta_data$T_refDepth), digits = 3)
  meta_data$AF[is.na(meta_data$AF)] <- 0
  
  if(class(meta_data$AF)=='list') {
    meta_data$AF = as.numeric(meta_data$AF)
  }
  
  # make filters logical
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 != "PASS" & meta_data$FILTER_Mutect2 != "MinAF"] <- FALSE
  meta_data$FILTER_Mutect2[is.na(meta_data$FILTER_Mutect2)] <- FALSE
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "PASS"] <- TRUE
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "MinAF"] <- TRUE
  meta_data$FILTER_Mutect2 <- as.logical(meta_data$FILTER_Mutect2)
  
  meta_data$FILTER_Freebayes[meta_data$FILTER_Freebayes != "PASS"] <- FALSE
  meta_data$FILTER_Freebayes[is.na(meta_data$FILTER_Freebayes)] <- FALSE
  meta_data$FILTER_Freebayes[meta_data$FILTER_Freebayes == "PASS"] <- TRUE
  meta_data$FILTER_Freebayes <- as.logical(meta_data$FILTER_Freebayes)
  
  meta_data$FILTER_Vardict[meta_data$FILTER_Vardict != "PASS"] <- FALSE
  meta_data$FILTER_Vardict[is.na(meta_data$FILTER_Vardict)] <- FALSE
  meta_data$FILTER_Vardict[meta_data$FILTER_Vardict == "PASS"] <- TRUE
  meta_data$FILTER_Vardict <- as.logical(meta_data$FILTER_Vardict)
  
  # meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "SpvFreq"] <- "PASS"
  # meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "REJECT;SpvFreq"] <- "PASS"
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan != "PASS"] <- FALSE
  meta_data$FILTER_Varscan[is.na(meta_data$FILTER_Varscan)] <- FALSE
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "PASS"] <- TRUE
  meta_data$FILTER_Varscan <- as.logical(meta_data$FILTER_Varscan)
  
  meta_data$FILTER_Strelka2[meta_data$FILTER_Strelka2 != "PASS" & meta_data$FILTER_Strelka2 != "MinAF"] <- FALSE
  meta_data$FILTER_Strelka2[is.na(meta_data$FILTER_Strelka2)] <- FALSE
  meta_data$FILTER_Strelka2[meta_data$FILTER_Strelka2 == "PASS"] <- TRUE
  meta_data$FILTER_Strelka2[meta_data$FILTER_Strelka2 == "MinAF"] <- TRUE
  meta_data$FILTER_Strelka2 <- as.logical(meta_data$FILTER_Strelka2)
  
  # get all 5 ref/alternates
  meta_data$REF_MFVdVs<-paste(meta_data$REF_Mutect2,meta_data$REF_Freebayes,meta_data$REF_Vardict,meta_data$REF_Varscan,meta_data$REF_Strelka2, sep ="/")
  meta_data$ALT_MFVdVs<-paste(meta_data$ALT_Mutect2,meta_data$ALT_Freebayes,meta_data$ALT_Vardict,meta_data$ALT_Varscan,meta_data$ALT_Strelka2, sep ="/")
  end.time=Sys.time()
  print(end.time-start.time)
  
  # sample name
  # meta_data$Sample_Name <- substrLeft(sampleid.t, 2)
  # meta_data$Sample_Name <- gsub(t.label,'',sampleid.t)
  meta_data$Sample_Name <- sampleid.t
  
  
  #filtering low read depth (in case of targeted sequencing or poor coverage) 
  meta_data = meta_data[((meta_data$N_refDepth!=0 & meta_data$N_altDepth!=0) | (meta_data$T_refDepth!=0 & meta_data$T_altDepth!=0)),]
  
  #insert missing feature(s)
  # meta_data$BaseQRankSum_Mutect2 <- NA
  meta_data$relcov <- (meta_data$T_refDepth+meta_data$T_altDepth+meta_data$N_refDepth+meta_data$N_altDepth)/median(meta_data$T_refDepth+meta_data$T_altDepth+meta_data$N_refDepth+meta_data$N_altDepth)
  
  #numerical
  meta_data$seqnames = gsub('chr','',meta_data$seqnames)
  
  #add "m2_BaseQRankSum" and "m2_ClippingRankSum" for long-term support
  meta_data$BaseQRankSum_Mutect2 = NA
  meta_data$ClippingRankSum_Mutect2 = NA
  
  #check for critical missing columns
  if(length(grep('MQ_Mutect2', colnames(meta_data)))!=1){meta_data$MQ_Mutect2 = NA} 
  if(length(grep('MQRankSum_Mutect2', colnames(meta_data)))!=1){meta_data$MQRankSum_Mutect2 = NA} 
  if(length(grep('NLOD_Mutect2', colnames(meta_data)))!=1){meta_data$NLOD_Mutect2 = NA} 
  if(length(grep('TLOD_Mutect2', colnames(meta_data)))!=1){meta_data$TLOD_Mutect2 = NA} 
  if(length(grep('FS_Mutect2', colnames(meta_data)))!=1){meta_data$FS_Mutect2 = NA} 
  if(length(grep('ReadPosRankSum_Mutect2', colnames(meta_data)))!=1){meta_data$ReadPosRankSum_Mutect2 = NA} 
  
  if(length(grep('MQM_Freebayes', colnames(meta_data)))!=1){meta_data$MQM_Freebayes = NA} 
  if(length(grep('MQMR_Freebayes', colnames(meta_data)))!=1){meta_data$MQMR_Freebayes = NA} 
  if(length(grep('GTI_Freebayes', colnames(meta_data)))!=1){meta_data$GTI_Freebayes = NA} 
  if(length(grep('LEN_Freebayes', colnames(meta_data)))!=1){meta_data$LEN_Freebayes = NA} 
  if(length(grep('ODDS_Freebayes', colnames(meta_data)))!=1){meta_data$ODDS_Freebayes = NA} 
  
  if(length(grep('SSC_Varscan', colnames(meta_data)))!=1){meta_data$SSC_Varscan = NA} 
  if(length(grep('SPV_Varscan', colnames(meta_data)))!=1){meta_data$SPV_Varscan = NA} 
  if(length(grep('GPV_Varscan', colnames(meta_data)))!=1){meta_data$GPV_Varscan = NA} 
  if(length(grep('SS_Varscan', colnames(meta_data)))!=1){meta_data$SS_Varscan = NA} 
  
  if(length(grep('SSF_Vardict', colnames(meta_data)))!=1){meta_data$SSF_Vardict = NA} 
  if(length(grep('MSI_Vardict', colnames(meta_data)))!=1){meta_data$MSI_Vardict = NA} 
  if(length(grep('SOR_Vardict', colnames(meta_data)))!=1){meta_data$SOR_Vardict = NA} 
  
  if(length(grep('QSS_Strelka2', colnames(meta_data)))!=1){meta_data$QSS_Strelka2 = NA} 
  if(length(grep('MSI_Strelka2', colnames(meta_data)))!=1){meta_data$MSI_Strelka2 = NA} 
  if(length(grep('SomaticEVS_Strelka2', colnames(meta_data)))!=1){meta_data$SomaticEVS_Strelka2 = NA} 
  if(length(grep('ReadPosRankSum_Strelka2', colnames(meta_data)))!=1){meta_data$ReadPosRankSum_Strelka2 = NA} 
  
  meta_data[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
               "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
               "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
               "SSF_Vardict","MSI_Vardict","SOR_Vardict",
               "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
               "AF",
               "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")] <- lapply (meta_data[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                                                                                                      "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                                                                                                      "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                                                                                                      "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                                                                                                      "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                                                                                                      "AF",
                                                                                                      "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")], as.numeric)
  # meta_data$AF = round(meta_data$AF,digits=3)
  
  #keep important columns and rename columns
  parse_snv <- meta_data[,c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                            "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                            "MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                            "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                            "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                            "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                            "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                            "AF",
                            "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")]
  colnames(parse_snv) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                           "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                           "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","m2_FS","m2_ReadPosRankSum","m2_BaseQRankSum","m2_ClippingRankSum",
                           "f_MQM","f_MQMR","f_GTI","f_LEN","f_ODDS",
                           "vs_SSC","vs_SPV","vs_GPV","vs_SS",
                           "vd_SSF","vd_MSI","vd_SOR",
                           "s2_QSS","s2_MQ","s2_SomaticEVS","s2_ReadPosRankSum",
                           "Alt_Allele_Freq",
                           "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
  
  #SOR inf values
  parse_snv$vd_SOR[which(is.infinite(parse_snv$vd_SOR)==TRUE)] <-  suppressWarnings(max(parse_snv$vd_SOR[which(is.infinite(parse_snv$vd_SOR)==F)],na.rm=T)+1) #added for v2
  
  #sort table
  parse_snv = parse_snv[order(parse_snv$START_POS_REF),]
  parse_snv = parse_snv[order(parse_snv$Chr),]
  
  parse_snv = parse_snv[!duplicated(parse_snv),]
  
  ### continue from here to format the output
  
  print("Parsing INDEL")
  
  #### Parsing indels ####
  # get passed indel calls from all callers
  indel_pass=c(rowRanges(vcf_m2[pass_m2&!snv_m2]), 
               rowRanges(vcf_f[pass_f&!snv_f]),  
               rowRanges(vcf_vs[pass_vs&!snv_vs]), 
               rowRanges(vcf_vd[pass_vd&!snv_vd]),
               rowRanges(vcf_s2[pass_s2&!snv_s2]),
               ignore.mcols=T)
  indel_pass= unique(indel_pass)
  
  # get the row numbers of calls that are passed by at least 1 caller in each vcf
  indel.m2.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_m2), type="equal"))
  indel.f.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_f), type="equal"))
  indel.vs.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_vs), type="equal"))
  indel.vd.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_vd), type="equal"))
  indel.s2.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_s2), type="equal"))
  
  print("extracting meta data from VRanges for indels")
  start.time=Sys.time() 
  
  if (length(indel.m2.index)!=0) {
  m2.na=F  
  # convert vcf to VRanges then to Granges, keep metadata columns
  vr_m2<- suppressWarnings(as(vcf_m2[indel.m2.index], "VRanges"))
  vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
  gr_m2=GRanges(vr_m2[[sampleid.t]])
  mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))
  
  gr_m2 <- unique(gr_m2)
  gr_m2$FILTER=vcf_m2[indel.m2.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from MuTect2.")
    m2.na=T
  }
  
  if (length(indel.f.index)!=0) {
  f.na=F  
  vr_f<- suppressWarnings(as(vcf_f[indel.f.index], "VRanges"))
  vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
  gr_f=GRanges(vr_f[[sampleid.t]])
  mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                            T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                            N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
  
  gr_f <- unique(gr_f)
  gr_f$FILTER=vcf_f[indel.f.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from FreeBayes.")
    f.na=T
  }
  
  if (length(indel.vs.index)!=0) {
  vs.na=F  
  vr_vs<- suppressWarnings(as(vcf_vs[indel.vs.index], "VRanges"))
  vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
  gr_vs=GRanges(vr_vs[[sampleid.t]])
  mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vs <- unique(gr_vs)
  gr_vs$FILTER=vcf_vs[indel.vs.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from VarScan.")
    vs.na=T
  }
  
  if (length(indel.vd.index)!=0) {
  vd.na=F  
  vr_vd<- suppressWarnings(as(vcf_vd[indel.vd.index], "VRanges"))
  vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
  gr_vd=GRanges(vr_vd[[sampleid.t]])
  mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vd <- unique(gr_vd)
  gr_vd$FILTER=vcf_vd[indel.vd.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from VarDict.")
    vd.na=T
  }
  
  if (length(indel.s2.index)!=0) {
    s2.na=F
    
    # convert vcf to VRanges then to Granges, keep metadata columns
    vr_s2<- suppressWarnings(as(vcf_s2[indel.s2.index], "VRanges"))
    vr_s2=GenomicRanges::split(vr_s2, vr_s2@sampleNames)
    gr_s2=GRanges(vr_s2[[sampleid.t]])
    mcols(gr_s2)=cbind(mcols(gr_s2), data.frame(REF=ref(vr_s2[[sampleid.t]]), ALT=alt(vr_s2[[sampleid.t]]), T_totalDepth=totalDepth(vr_s2[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_s2[[sampleid.t]]), T_altDepth=altDepth(vr_s2[[sampleid.t]]),N_totalDepth=totalDepth(vr_s2[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_s2[[sampleid.n]]), N_altDepth=altDepth(vr_s2[[sampleid.n]]), stringsAsFactors=F))
    
    gr_s2 <- unique(gr_s2)
    gr_s2$FILTER=vcf_s2[indel.s2.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed INDELs from Strelka2.")
    s2.na=T
  }
  
  
  end.time=Sys.time() 
  print(end.time-start.time)
  
  print("merge 5 vcfs and meta data")
  start.time=Sys.time()
  
  ## merge 5 vcfs and meta data
  m2.cols = c(names(vcf_m2@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  f.cols = c(names(vcf_f@info@listData),'RO','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  vs.cols = c(names(vcf_vs@info@listData),'FREQ','RD','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  vd.cols = c(names(vcf_vd@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  s2.cols = c(names(vcf_s2@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  
  names_m2= paste(m2.cols, "_Mutect2", sep="")
  names_f= paste(f.cols, "_Freebayes", sep="")
  names_vs= paste(vs.cols, "_Varscan", sep="")
  names_vd= paste(vd.cols, "_Vardict", sep="")
  names_s2= paste(s2.cols, "_Strelka2", sep="")
  
  # names_m2= paste(names(mcols(gr_m2)[,-1]), "_Mutect2", sep="")
  # names_f= paste(names(mcols(gr_f)[,-1]), "_Freebayes", sep="")
  # names_vs= paste(names(mcols(gr_vs)[,-1]), "_Varscan", sep="")
  # names_vd= paste(names(mcols(gr_vd)[,-1]), "_Vardict", sep="")
  # names_s2= paste(names(mcols(gr_s2)[,-1]), "_Strelka2", sep="")
  
  meta_indel=data.frame(indel_pass)[,1:3]
  meta_indel[c(names_m2,names_f,names_vs,names_vd,names_s2)]=NA
  
  #do not merge when index=0, caller.na=T
  if (m2.na==F) {
    meta_indel[Biostrings::match(gr_m2, indel_pass), names_m2]=data.frame(mcols(gr_m2)[,m2.cols])
  }
  if (f.na==F) {
    meta_indel[Biostrings::match(gr_f, indel_pass), names_f]=data.frame(mcols(gr_f)[,f.cols])
  }
  if (vs.na==F) {
    meta_indel[Biostrings::match(gr_vs, indel_pass), names_vs]=data.frame(mcols(gr_vs)[,vs.cols])
  }
  if (vd.na==F) {
    meta_indel[Biostrings::match(gr_vd, indel_pass), names_vd]=data.frame(mcols(gr_vd)[,vd.cols])
  }
  if (s2.na==F) {
    meta_indel[Biostrings::match(gr_s2, indel_pass), names_s2]=data.frame(mcols(gr_s2)[,s2.cols])
  }
  
  # #do not merge when the 1st row is all NAs + last column PASSED/REJECT
  # if ((rowSums(is.na(data.frame(mcols(gr_m2)[1,1:length(mcols(gr_m2))-1])))!=ncol(data.frame(mcols(gr_m2)[1,1:length(mcols(gr_m2))-1])))==TRUE) {
  #   meta_indel[Biostrings::match(gr_m2, indel_pass), names_m2]=data.frame(mcols(gr_m2)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_f)[1,1:length(mcols(gr_f))-1])))!=ncol(data.frame(mcols(gr_f)[1,1:length(mcols(gr_f))-1])))==TRUE) {
  #   meta_indel[Biostrings::match(gr_f, indel_pass), names_f]=data.frame(mcols(gr_f)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_vs)[1,1:length(mcols(gr_vs))-1])))!=ncol(data.frame(mcols(gr_vs)[1,1:length(mcols(gr_vs))-1])))==TRUE) {
  #   meta_indel[Biostrings::match(gr_vs, indel_pass), names_vs]=data.frame(mcols(gr_vs)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_vd)[1,1:length(mcols(gr_vd))-1])))!=ncol(data.frame(mcols(gr_vd)[1,1:length(mcols(gr_vd))-1])))==TRUE) {
  #   meta_indel[Biostrings::match(gr_vd, indel_pass), names_vd]=data.frame(mcols(gr_vd)[,-1])
  # }
  # if ((rowSums(is.na(data.frame(mcols(gr_s2)[1,1:length(mcols(gr_s2))-1])))!=ncol(data.frame(mcols(gr_s2)[1,1:length(mcols(gr_s2))-1])))==TRUE) {
  #   meta_indel[Biostrings::match(gr_s2, indel_pass), names_s2]=data.frame(mcols(gr_s2)[,-1])
  # }


  end.time=Sys.time()
  print(end.time-start.time)
  
  #add in indel cases from snv matrix
  meta_indel <- rbind(meta_indel,meta_indel_cases)
  meta_indel <- unique(meta_indel) #remove all duplicate snv/indel cases
  
  print("formating")
  start.time=Sys.time()
  # extract reference allele
  # ref=meta_indel[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict")] 
  # ref.ind=which(!is.na(ref), arr.ind=T) # find non-NA alleles
  # ref.ind=ref.ind[order(ref.ind[,1]),] #order array indices by row number
  # ref.ind=ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
  # meta_indel$REF=ref[ref.ind]
  ref <- meta_indel[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict","REF_Strelka2")] 
  ref.nchar <- sapply(ref, nchar)
  ref.nchar <- as.data.frame(ref.nchar)
  ref.nchar$MAX <- apply(ref.nchar,MARGIN = 1,function(x) max(x,na.rm=TRUE))
  ref.max <- ref.nchar[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict","REF_Strelka2")]
  ref.max[,c("REF_Mutect2")] <- ref.nchar[,c("REF_Mutect2")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Freebayes")] <- ref.max[,c("REF_Freebayes")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Varscan")] <- ref.max[,c("REF_Varscan")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Vardict")] <- ref.max[,c("REF_Vardict")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Strelka2")] <- ref.nchar[,c("REF_Strelka2")] - ref.nchar[,c("MAX")]
  ref.ind <- which(!is.na(ref.max)&(ref.max)==0, arr.ind=T) # find non-NA alleles
  ref.ind <- ref.ind[order(ref.ind[,1]),] #order array indices by row number
  ref.ind <- ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
  meta_indel$REF <- ref[ref.ind]
  
  # extract alternate allele
  # alt=meta_indel[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict")]
  # alt.ind=which(!is.na(alt), arr.ind=T) #find all non-NA alleles
  # alt.ind=alt.ind[order(alt.ind[,1]),]  #order array indices by row number
  # alt.ind=alt.ind[!duplicated(alt.ind[,1]),]# take the first non-NA allele of each row
  # meta_indel$ALT=alt[alt.ind]
  alt <- meta_indel[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict","ALT_Strelka2")] 
  alt.nchar <- sapply(alt, nchar)
  alt.nchar <- as.data.frame(alt.nchar)
  alt.nchar$MAX <- apply(alt.nchar,MARGIN = 1,function(x) max(x,na.rm=TRUE))
  alt.max <- alt.nchar[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict","ALT_Strelka2")]
  alt.max[,c("ALT_Mutect2")] <- alt.nchar[,c("ALT_Mutect2")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Freebayes")] <- alt.max[,c("ALT_Freebayes")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Varscan")] <- alt.max[,c("ALT_Varscan")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Vardict")] <- alt.max[,c("ALT_Vardict")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Strelka2")] <- alt.nchar[,c("ALT_Strelka2")] - alt.nchar[,c("MAX")]
  alt.ind <- which(!is.na(alt.max)&(alt.max)==0, arr.ind=T) # find non-NA alleles
  alt.ind <- alt.ind[order(alt.ind[,1]),] #order array indices by row number
  alt.ind <- alt.ind[!duplicated(alt.ind[,1]),] #get first non-NA value for each row
  meta_indel$ALT <- alt[alt.ind]
  
  # check for list objects in columns
  for (i in 1:ncol(meta_indel)) {
    if(class(meta_indel[,i])=='list'){
      meta_indel[,i] = unlist(meta_indel[,i])
    }
  }
  
  # calculate alt allele frequency
  # meta_indel$T_refDepth_Varscan <- meta_indel$T_totalDepth_Varscan-meta_indel$T_altDepth_Varscan
  # meta_indel$N_refDepth_Varscan <- meta_indel$N_totalDepth_Varscan-meta_indel$N_altDepth_Varscan
  meta_indel$T_altDepth_Freebayes <- meta_indel$T_totalDepth_Freebayes-meta_indel$RO_Freebayes
  meta_indel$T_refDepth_Freebayes <- meta_indel$RO_Freebayes
  meta_indel$AF_Freebayes <- meta_indel$T_altDepth_Freebayes/meta_indel$T_totalDepth_Freebayes
  
  meta_indel$T_refDepth_Varscan <- meta_indel$RD_Varscan
  meta_indel$T_altDepth_Varscan <- meta_indel$T_totalDepth_Varscan-meta_indel$RD_Varscan
  
  meta_indel$AF_Strelka2[sapply(meta_indel$AF_Strelka2, is.null)] <- NA
  meta_indel$AF_Strelka2 <- unlist(meta_indel$AF_Strelka2)
  meta_indel$T_altDepth_Strelka2 <- round(meta_indel$AF_Strelka2*meta_indel$T_totalDepth_Strelka2)
  meta_indel$T_refDepth_Strelka2 <- meta_indel$T_totalDepth_Strelka2-meta_indel$T_altDepth_Strelka2
  

  #catch vcf caller mapping errors
  if(all(is.na(meta_indel$T_altDepth_Mutect2))==T | 
     all(is.na(meta_indel$T_refDepth_Mutect2))==T |
     all(is.na(meta_indel$N_altDepth_Mutect2))==T | 
     all(is.na(meta_indel$N_refDepth_Mutect2))==T) {
    meta_indel$T_refDepth_Mutect2 <- NA
    meta_indel$N_refDepth_Mutect2 <- NA
    meta_indel$T_altDepth_Mutect2 <- NA
    meta_indel$N_altDepth_Mutect2 <- NA
  }
  if(all(is.na(meta_indel$T_altDepth_Freebayes))==T | 
     all(is.na(meta_indel$T_refDepth_Freebayes))==T 
     # all(is.na(meta_indel$N_altDepth_Freebayes))==T | 
     # all(is.na(meta_indel$N_refDepth_Freebayes))==T
  ) {
    meta_indel$T_refDepth_Freebayes <- NA
    # meta_indel$N_refDepth_Freebayes <- NA
    meta_indel$T_altDepth_Freebayes <- NA
    # meta_indel$N_altDepth_Freebayes <- NA
  }
  if(all(is.na(meta_indel$T_altDepth_Varscan))==T | 
     all(is.na(meta_indel$T_refDepth_Varscan))==T 
     # all(is.na(meta_indel$N_altDepth_Varscan))==T | 
     # all(is.na(meta_indel$N_refDepth_Varscan))==T
  ) {
    meta_indel$T_refDepth_Varscan <- NA
    # meta_indel$N_refDepth_Varscan <- NA
    meta_indel$T_altDepth_Varscan <- NA
    # meta_indel$N_altDepth_Varscan <- NA
  }
  if(all(is.na(meta_indel$T_altDepth_Vardict))==T | 
     all(is.na(meta_indel$T_refDepth_Vardict))==T |
     all(is.na(meta_indel$N_altDepth_Vardict))==T | 
     all(is.na(meta_indel$N_refDepth_Vardict))==T) {
    meta_indel$T_refDepth_Vardict <- NA
    meta_indel$N_refDepth_Vardict <- NA
    meta_indel$T_altDepth_Vardict <- NA
    meta_indel$N_altDepth_Vardict <- NA
  }
  if(all(is.na(meta_indel$T_altDepth_Strelka2))==T | 
     all(is.na(meta_indel$T_refDepth_Strelka2))==T 
     # all(is.na(meta_indel$N_altDepth_Strelka2))==T | 
     # all(is.na(meta_indel$N_refDepth_Strelka2))==T
  ) {
    meta_indel$T_refDepth_Strelka2 <- NA
    # meta_indel$N_refDepth_Strelka2 <- NA
    meta_indel$T_altDepth_Strelka2 <- NA
    # meta_indel$N_altDepth_Strelka2 <- NA
  }
  
  # calculate mean depth
  meta_indel$N_refDepth <- round(rowMeans(meta_indel[, c("N_refDepth_Mutect2", "N_refDepth_Freebayes", "N_refDepth_Varscan", "N_refDepth_Vardict","N_refDepth_Strelka2")],na.rm = TRUE))
  meta_indel$N_altDepth <- round(rowMeans(meta_indel[, c("N_altDepth_Mutect2", "N_altDepth_Freebayes", "N_altDepth_Varscan", "N_altDepth_Vardict","N_altDepth_Strelka2")],na.rm = TRUE))
  meta_indel$T_refDepth <- round(rowMeans(meta_indel[, c("T_refDepth_Mutect2", "T_refDepth_Freebayes", "T_refDepth_Varscan", "T_refDepth_Vardict","T_refDepth_Strelka2")],na.rm = TRUE))
  meta_indel$T_altDepth <- round(rowMeans(meta_indel[, c("T_altDepth_Mutect2", "T_altDepth_Freebayes", "T_altDepth_Varscan", "T_altDepth_Vardict","T_altDepth_Strelka2")],na.rm = TRUE))

  meta_indel$N_refDepth[is.nan(meta_indel$N_refDepth)] <- 0
  meta_indel$N_altDepth[is.nan(meta_indel$N_altDepth)] <- 0
  meta_indel$T_refDepth[is.nan(meta_indel$T_refDepth)] <- 0
  meta_indel$T_altDepth[is.nan(meta_indel$T_altDepth)] <- 0
  
  #re-calculate mean AF
  meta_indel$AF = round(meta_indel$T_altDepth/(meta_indel$T_altDepth+meta_indel$T_refDepth), digits = 3)
  meta_indel$AF[is.na(meta_indel$AF)] <- 0
  
  if(class(meta_indel$AF)=='list') {
    meta_indel$AF = as.numeric(meta_indel$AF)
  }
  
  # make filters logical
  meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 != "PASS" & meta_indel$FILTER_Mutect2 != "MinAF"] <- FALSE
  meta_indel$FILTER_Mutect2[is.na(meta_indel$FILTER_Mutect2)] <- FALSE
  meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 == "PASS"] <- TRUE
  meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 == "MinAF"] <- TRUE
  meta_indel$FILTER_Mutect2 <- as.logical(meta_indel$FILTER_Mutect2)
  
  meta_indel$FILTER_Freebayes[meta_indel$FILTER_Freebayes != "PASS"] <- FALSE
  meta_indel$FILTER_Freebayes[is.na(meta_indel$FILTER_Freebayes)] <- FALSE
  meta_indel$FILTER_Freebayes[meta_indel$FILTER_Freebayes == "PASS"] <- TRUE
  meta_indel$FILTER_Freebayes <- as.logical(meta_indel$FILTER_Freebayes)
  
  meta_indel$FILTER_Vardict[meta_indel$FILTER_Vardict != "PASS"] <- FALSE
  meta_indel$FILTER_Vardict[is.na(meta_indel$FILTER_Vardict)] <- FALSE
  meta_indel$FILTER_Vardict[meta_indel$FILTER_Vardict == "PASS"] <- TRUE
  meta_indel$FILTER_Vardict <- as.logical(meta_indel$FILTER_Vardict)
  
  # meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan == "SpvFreq"] <- "PASS"
  # meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan == "REJECT;SpvFreq"] <- "PASS"
  meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan != "PASS"] <- FALSE
  meta_indel$FILTER_Varscan[is.na(meta_indel$FILTER_Varscan)] <- FALSE
  meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan == "PASS"] <- TRUE
  meta_indel$FILTER_Varscan <- as.logical(meta_indel$FILTER_Varscan)
  
  meta_indel$FILTER_Strelka2[meta_indel$FILTER_Strelka2 != "PASS" & meta_indel$FILTER_Strelka2 != "MinAF"] <- FALSE
  meta_indel$FILTER_Strelka2[is.na(meta_indel$FILTER_Strelka2)] <- FALSE
  meta_indel$FILTER_Strelka2[meta_indel$FILTER_Strelka2 == "PASS"] <- TRUE
  meta_indel$FILTER_Strelka2[meta_indel$FILTER_Strelka2 == "MinAF"] <- TRUE
  meta_indel$FILTER_Strelka2 <- as.logical(meta_indel$FILTER_Strelka2)
  
  # get all 4 ref/alternates
  meta_indel$REF_MFVdVs<-paste(meta_indel$REF_Mutect2,meta_indel$REF_Freebayes,meta_indel$REF_Vardict,meta_indel$REF_Varscan,meta_indel$REF_Strelka2, sep ="/")
  meta_indel$ALT_MFVdVs	<-paste(meta_indel$ALT_Mutect2,meta_indel$ALT_Freebayes,meta_indel$ALT_Vardict,meta_indel$ALT_Varscan,meta_indel$ALT_Strelka2, sep ="/")
  
  # sample name
  # meta_indel$Sample_Name <- substrLeft(sampleid.t, 2)
  # meta_indel$Sample_Name <- gsub(t.label,'',sampleid.t)
  meta_indel$Sample_Name <- sampleid.t
  
  #filtering low read depth (in case of targeted sequencing or poor coverage) 
  meta_indel = meta_indel[((meta_indel$N_refDepth!=0 & meta_indel$N_altDepth!=0) | (meta_indel$T_refDepth!=0 & meta_indel$T_altDepth!=0)),]
  
  #insert missing feature
  # meta_indel$BaseQRankSum_Mutect2 <- NA
  meta_indel$relcov <- (meta_indel$T_refDepth+meta_indel$T_altDepth+meta_indel$N_refDepth+meta_indel$N_altDepth)/median(meta_indel$T_refDepth+meta_indel$T_altDepth+meta_indel$N_refDepth+meta_indel$N_altDepth)
  
  #add "m2_BaseQRankSum" and "m2_ClippingRankSum" for long-term support
  meta_indel$BaseQRankSum_Mutect2 = NA
  meta_indel$ClippingRankSum_Mutect2 = NA
  
  #numerical
  meta_indel$seqnames = gsub('chr','',meta_indel$seqnames)
  
  #check for critical missing columns
  if(length(grep('MQ_Mutect2', colnames(meta_indel)))!=1){meta_indel$MQ_Mutect2 = NA} 
  if(length(grep('MQRankSum_Mutect2', colnames(meta_indel)))!=1){meta_indel$MQRankSum_Mutect2 = NA} 
  if(length(grep('NLOD_Mutect2', colnames(meta_indel)))!=1){meta_indel$NLOD_Mutect2 = NA} 
  if(length(grep('TLOD_Mutect2', colnames(meta_indel)))!=1){meta_indel$TLOD_Mutect2 = NA} 
  if(length(grep('FS_Mutect2', colnames(meta_indel)))!=1){meta_indel$FS_Mutect2 = NA} 
  if(length(grep('ReadPosRankSum_Mutect2', colnames(meta_indel)))!=1){meta_indel$ReadPosRankSum_Mutect2 = NA} 
  
  if(length(grep('MQM_Freebayes', colnames(meta_indel)))!=1){meta_indel$MQM_Freebayes = NA} 
  if(length(grep('MQMR_Freebayes', colnames(meta_indel)))!=1){meta_indel$MQMR_Freebayes = NA} 
  if(length(grep('GTI_Freebayes', colnames(meta_indel)))!=1){meta_indel$GTI_Freebayes = NA} 
  if(length(grep('LEN_Freebayes', colnames(meta_indel)))!=1){meta_indel$LEN_Freebayes = NA} 
  if(length(grep('ODDS_Freebayes', colnames(meta_indel)))!=1){meta_indel$ODDS_Freebayes = NA} 
  
  if(length(grep('SSC_Varscan', colnames(meta_indel)))!=1){meta_indel$SSC_Varscan = NA} 
  if(length(grep('SPV_Varscan', colnames(meta_indel)))!=1){meta_indel$SPV_Varscan = NA} 
  if(length(grep('GPV_Varscan', colnames(meta_indel)))!=1){meta_indel$GPV_Varscan = NA} 
  if(length(grep('SS_Varscan', colnames(meta_indel)))!=1){meta_indel$SS_Varscan = NA} 
  
  if(length(grep('SSF_Vardict', colnames(meta_indel)))!=1){meta_indel$SSF_Vardict = NA} 
  if(length(grep('MSI_Vardict', colnames(meta_indel)))!=1){meta_indel$MSI_Vardict = NA} 
  if(length(grep('SOR_Vardict', colnames(meta_indel)))!=1){meta_indel$SOR_Vardict = NA} 
  
  if(length(grep('QSS_Strelka2', colnames(meta_indel)))!=1){meta_indel$QSS_Strelka2 = NA} 
  if(length(grep('MSI_Strelka2', colnames(meta_indel)))!=1){meta_indel$MSI_Strelka2 = NA} 
  if(length(grep('SomaticEVS_Strelka2', colnames(meta_indel)))!=1){meta_indel$SomaticEVS_Strelka2 = NA} 
  if(length(grep('ReadPosRankSum_Strelka2', colnames(meta_indel)))!=1){meta_indel$ReadPosRankSum_Strelka2 = NA} 
  
  meta_indel[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                "AF",
                "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")] <- lapply (meta_indel[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                                                                                                        "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                                                                                                        "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                                                                                                        "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                                                                                                        "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                                                                                                        "AF",
                                                                                                        "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")], as.numeric)
  
  parse_indel <- meta_indel[,c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                               "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                               "MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                               "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                               "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                               "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                               "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                               "AF",
                               "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")]
  
  colnames(parse_indel) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                             "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                             "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","m2_FS","m2_ReadPosRankSum","m2_BaseQRankSum","m2_ClippingRankSum",
                             "f_MQM","f_MQMR","f_GTI","f_LEN","f_ODDS",
                             "vs_SSC","vs_SPV","vs_GPV","vs_SS",
                             "vd_SSF","vd_MSI","vd_SOR",
                             "s2_QSS","s2_MQ","s2_SomaticEVS","s2_ReadPosRankSum",
                             "Alt_Allele_Freq",
                             "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
  
  #SOR inf values
  parse_indel$vd_SOR[which(is.infinite(parse_indel$vd_SOR)==TRUE)] <-  suppressWarnings(max(parse_indel$vd_SOR[which(is.infinite(parse_indel$vd_SOR)==F)],na.rm=T)+1) #added for v2
  
  #sort table
  parse_indel = parse_indel[order(parse_indel$START_POS_REF),]
  parse_indel = parse_indel[order(parse_indel$Chr),]
  
  parse_indel = parse_indel[!duplicated(parse_indel),]
  
  return(list(snv=parse_snv, indel=parse_indel))
}



