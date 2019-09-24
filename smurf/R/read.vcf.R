#' parsevcf
#'
#' Step 1 Initial parse and filtering
#'
#'@param x List object containing the four vcf.gz files from 
#'callers MuTect2, Freebayes, VarDict and VarScan.
#'@param roi set to TRUE to enable this function. Alternatively read.vcf is a VCF parser in development.
#'@param roi.dir Directory to BED file of genomic coordinates.
#'@param t.label Tumor sample name label     
#' @examples
#' 
#' 
#' @export
read.vcf = function(x, roi=F, roi.dir=NULL, t.label){
  
  print("Parsing step")
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  substrLeft <- function(x, n){
    substr(x, 1, nchar(x)-n)
  }
  
  print("reading vcfs")
  start.time=Sys.time()

  #ScanVCF with required parameters
  svp_m<-ScanVcfParam(info=c("FS","MQ","MQRankSum","NLOD","ReadPosRankSum","TLOD"), samples=scanVcfHeader(x[[1]])@samples, geno=c("AD", "AF", "DP"))
  svp_f<-ScanVcfParam(info=c("GTI","LEN","MQM","MQMR","ODDS"),samples=scanVcfHeader(x[[2]])@samples, geno=c("RO","DP"))
  svp_vs<-ScanVcfParam(info=c("SSC","GPV","SS","SPV"), samples=scanVcfHeader(x[[3]])@samples, geno=c("AD","FREQ", "DP"))
  svp_vd<-ScanVcfParam(info=c("SSF","MSI","SOR"), samples=scanVcfHeader(x[[4]])@samples, geno=c("AD", "AF", "DP"))
  
  
  ## read only chromosomes 1-22, X, Y and M
  print("reading mutect2")
  vcf_m2<- readVcf(x[[1]], "hg19", svp_m)
  print("reading freebayes")
  vcf_f<- readVcf(x[[2]], "hg19", svp_f)
  print("reading varscan")
  vcf_vs<- readVcf(x[[3]], "hg19", svp_vs)
  print("reading vardict")
  vcf_vd<- readVcf(x[[4]], "hg19", svp_vd)
  
  
  # remove duplicates from vcfs
  vcf_m2=vcf_m2[!duplicated(rowRanges(vcf_m2))]
  vcf_f=vcf_f[!duplicated(rowRanges(vcf_f))]
  vcf_vs=vcf_vs[!duplicated(rowRanges(vcf_vs))]
  vcf_vd=vcf_vd[!duplicated(rowRanges(vcf_vd))]
  end.time=Sys.time() 
  
  print(end.time-start.time)
  
  # get vcf headers
  H_m2=header(vcf_m2)
  H_f=header(vcf_f)  
  H_vs=header(vcf_vs) 
  H_vd=header(vcf_vd)
  
  # sample name
  # sampleid.t <-H_m2@samples[grep(substrRight(H_m2@samples, 2), pattern="-T")]
  # sampleid.n <-H_m2@samples[grep(substrRight(H_m2@samples, 2), pattern="-T", invert=T)]
  sampleid.t <-H_m2@samples[grep((H_m2@samples), pattern=t.label)]
  sampleid.n <-H_m2@samples[grep((H_m2@samples), pattern=t.label, invert=T)]
  
  if(length(sampleid.t)!=1 & length(sampleid.n)!=1) {
    stop('t.label for tumor sample is duplicated or missing')
  }
  
  # print("extracting calls passed by at least 1 caller")
  start.time=Sys.time()
  # get rows that are SNVs
  snv_m2<-isSNV(vcf_m2, singleAltOnly=FALSE)
  snv_f<-isSNV(vcf_f, singleAltOnly=FALSE)
  snv_vs<-isSNV(vcf_vs, singleAltOnly=FALSE)
  snv_vd<-isSNV(vcf_vd, singleAltOnly=FALSE)  
  
  if (roi==F) { 
  
  # extract passed calls from each caller
  pass_m2<- vcf_m2@fixed$FILTER=="PASS"
  pass_f<- vcf_f@fixed$FILTER=="PASS" 
  pass_vs<- vcf_vs@fixed$FILTER=="PASS" | vcf_vs@fixed$FILTER=="SpvFreq" | vcf_vs@fixed$FILTER=="REJECT;SpvFreq"
  pass_vd<- vcf_vd@fixed$FILTER=="PASS"
  
  } else if (roi==T) {
    
    # if (is.null(roi.dir)) {
    #   
    #   print('Subsetting genomic coordinate')
    #   
    #   pass_m2<- as.numeric(seqnames(vcf_m2@rowRanges))==chr & start(vcf_m2@rowRanges)==start & end(vcf_m2@rowRanges)==end
    #   pass_f<- as.numeric(seqnames(vcf_f@rowRanges))==chr & start(vcf_f@rowRanges)==start & end(vcf_f@rowRanges)==end
    #   pass_vs<- as.numeric(seqnames(vcf_vs@rowRanges))==chr & start(vcf_vs@rowRanges)==start & end(vcf_vs@rowRanges)==end
    #   pass_vd<- as.numeric(seqnames(vcf_vd@rowRanges))==chr & start(vcf_vd@rowRanges)==start & end(vcf_vd@rowRanges)==end
    # 
    # } else {
      
    print('Subsetting ROI region')
    
    roi.df = read.delim(roi.dir, header = F, colClasses = c('factor','integer','integer'))
    colnames(roi.df) = c('chrom','start','end')
    roi.df$chrom = gsub('chr', '', roi.df$chrom)
    roi.gr <- with(roi.df, GRanges(chrom, IRanges(start, end)))
    roi.gr = GenomicRanges::reduce(roi.gr)
    
    ovl_m2 = suppressWarnings(IRanges::findOverlaps(vcf_m2@rowRanges, roi.gr))
    ovl_f = suppressWarnings(IRanges::findOverlaps(vcf_f@rowRanges, roi.gr))
    ovl_vs = suppressWarnings(IRanges::findOverlaps(vcf_vs@rowRanges, roi.gr))
    ovl_vd = suppressWarnings(IRanges::findOverlaps(vcf_vd@rowRanges, roi.gr))
    
    # pass_m2 = suppressWarnings(vcf_m2@rowRanges == vcf_m2@rowRanges[queryHits(ovl_m2),])
    # pass_f = suppressWarnings(vcf_f@rowRanges == vcf_f@rowRanges[queryHits(ovl_f),])
    # pass_vs = suppressWarnings(vcf_vs@rowRanges == vcf_vs@rowRanges[queryHits(ovl_vs),])
    # pass_vd = suppressWarnings(vcf_vd@rowRanges == vcf_vd@rowRanges[queryHits(ovl_vd),])
    
    # pass_m2 = names(vcf_m2@rowRanges) %in% names(vcf_m2@rowRanges[queryHits(ovl_m2),])
    # pass_f = names(vcf_f@rowRanges) %in% names(vcf_f@rowRanges[queryHits(ovl_f),])
    # pass_vs = names(vcf_vs@rowRanges) %in% names(vcf_vs@rowRanges[queryHits(ovl_vs),])
    # pass_vd = names(vcf_vd@rowRanges) %in% names(vcf_vd@rowRanges[queryHits(ovl_vd),])
    
    pass_m2 = rep(F, length(vcf_m2@rowRanges))
    pass_f = rep(F, length(vcf_f@rowRanges))
    pass_vs = rep(F, length(vcf_vs@rowRanges))
    pass_vd = rep(F, length(vcf_vd@rowRanges))
    
    pass_m2[queryHits(ovl_m2)] = T 
    pass_f[queryHits(ovl_f)] = T
    pass_vs[queryHits(ovl_vs)] = T 
    pass_vd[queryHits(ovl_vd)] = T

  # }

  }
  
  # get passed snv calls from all callers
  snv_pass=c(rowRanges(vcf_m2[pass_m2&snv_m2]), rowRanges(vcf_f[pass_f&snv_f]),  rowRanges(vcf_vs[pass_vs&snv_vs]), rowRanges(vcf_vd[pass_vd&snv_vd]),ignore.mcols=T)
  snv_pass= unique(snv_pass)
  
  if(length(snv_pass)!=0) {
    
  ###### parsing SNVs
  # get the row numbers of calls that are passed by at least 1 caller in each vcf
  snv.m2.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_m2), type="equal"))
  snv.f.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_f), type="equal"))
  snv.vs.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vs), type="equal"))
  snv.vd.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vd), type="equal"))
  
  end.time=Sys.time() 
  print(end.time-start.time)
  
  print("extracting meta data from VRanges")
  start.time=Sys.time() 
  
  if (length(snv.m2.index)!=0) {
  # convert vcf to VRanges then to Granges, keep metadata columns
  vr_m2<- as(vcf_m2[snv.m2.index], "VRanges")
  #mcols(vr_m2)=vr_m2[,c("MQ","MQRankSum","TLOD","NLOD","AF")]
  
  #vr_m2=split(vr_m2, vr_m2@sampleNames)
  #vr_m2=GenomicRanges::split(vr_m2, levels(vr_m2@sampleNames))
  vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
  

  gr_m2=GRanges(vr_m2[[sampleid.t]])

  mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))

  gr_m2 <- unique(gr_m2)

  gr_m2$FILTER=vcf_m2[snv.m2.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from MuTect2.")
    vr_m2<- as(vcf_m2, "VRanges")
    vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
    gr_m2=GRanges(vr_m2[[sampleid.t]])
    mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                              T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                              N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_m2 <- unique(gr_m2)
    gr_m2$FILTER=vcf_m2@fixed$FILTER
    
    gr_m2 <- gr_m2[1,]
    gr_m2@elementMetadata@listData$QUAL <- NA
    gr_m2@elementMetadata@listData$MQ <- NA
    gr_m2@elementMetadata@listData$MQRankSum <- NA
    gr_m2@elementMetadata@listData$NLOD <- NA
    gr_m2@elementMetadata@listData$TLOD <- NA
    # gr_m2@elementMetadata@listData$ClippingRankSum <- NA
    gr_m2@elementMetadata@listData$ReadPosRankSum <- NA
    gr_m2@elementMetadata@listData$FS <- NA
    gr_m2@elementMetadata@listData$AF <- NA
    

  }
  
  
  if (length(snv.f.index)!=0) {
    
  vr_f<- suppressWarnings(as(vcf_f[snv.f.index], "VRanges"))
  #mcols(vr_f)=vr_f[,c("MQM","MQMR")]
  vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
  gr_f=GRanges(vr_f[[sampleid.t]])
  mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                            T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                            N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
  
  gr_f <- unique(gr_f)
  gr_f$FILTER=vcf_f[snv.f.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from FreeBayes.")
    vr_f<- suppressWarnings(as(vcf_f, "VRanges"))
    vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
    gr_f=GRanges(vr_f[[sampleid.t]])
    mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                              T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                              N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_f <- unique(gr_f)
    gr_f$FILTER=vcf_f@fixed$FILTER
    
    gr_f <- gr_f[1,]
    gr_f@elementMetadata@listData$QUAL <- NA
    gr_f@elementMetadata@listData$LEN <- NA
    gr_f@elementMetadata@listData$MQM <- NA
    gr_f@elementMetadata@listData$MQMR <- NA
    gr_f@elementMetadata@listData$ODDS <- NA
    gr_f@elementMetadata@listData$GTI <- NA
    gr_f@elementMetadata@listData$RO <- NA
    
  }
  
  
  if (length(snv.vs.index)!=0) {
    
  vr_vs<- suppressWarnings(as(vcf_vs[snv.vs.index], "VRanges"))
  #mcols(vr_vs)=vr_vs[,c("SSC","SPV", "FREQ")]
  vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
  gr_vs=GRanges(vr_vs[[sampleid.t]])
  mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vs <- unique(gr_vs)
  gr_vs$FILTER=vcf_vs[snv.vs.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from VarScan.")
    vr_vs<- suppressWarnings(as(vcf_vs, "VRanges"))
    vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
    gr_vs=GRanges(vr_vs[[sampleid.t]])
    mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                                T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                                N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_vs <- unique(gr_vs)
    gr_vs$FILTER=vcf_vs@fixed$FILTER
    
    gr_vs <- gr_vs[1,]
    gr_vs@elementMetadata@listData$QUAL <- NA
    gr_vs@elementMetadata@listData$SSC <- NA
    gr_vs@elementMetadata@listData$SPV <- NA
    gr_vs@elementMetadata@listData$GPV <- NA
    gr_vs@elementMetadata@listData$SS <- NA
    gr_vs@elementMetadata@listData$FREQ <- NA
    
  }
  
  if (length(snv.vd.index)!=0) {
    
  vr_vd<- suppressWarnings(as(vcf_vd[snv.vd.index], "VRanges"))
  #mcols(vr_vd)=vr_vd[,c("SSF","SOR","MSI", "AF")]
  vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
  gr_vd=GRanges(vr_vd[[sampleid.t]])
  mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vd <- unique(gr_vd)
  gr_vd$FILTER=vcf_vd[snv.vd.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed SNVs from VarDict.")
    vr_vd<- suppressWarnings(as(vcf_vd, "VRanges"))
    vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
    gr_vd=GRanges(vr_vd[[sampleid.t]])
    mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                                T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                                N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_vd <- unique(gr_vd)
    gr_vd$FILTER=vcf_vd@fixed$FILTER
    
    gr_vd <- gr_vd[1,]
    gr_vd@elementMetadata@listData$QUAL <- NA
    gr_vd@elementMetadata@listData$SSF <- NA
    gr_vd@elementMetadata@listData$MSI <- NA
    gr_vd@elementMetadata@listData$SOR <- NA
    gr_vd@elementMetadata@listData$AF <- NA
    

  }
  
  end.time=Sys.time() 
  print(end.time-start.time)
  
  print("merge 4 vcfs and meta data")
  start.time=Sys.time()
  
  ## merge 4 vcfs and meta data
  names_m2= paste(names(mcols(gr_m2)[,-1]), "_Mutect2", sep="")
  names_f= paste(names(mcols(gr_f)[,-1]), "_Freebayes", sep="")
  names_vs= paste(names(mcols(gr_vs)[,-1]), "_Varscan", sep="")
  names_vd= paste(names(mcols(gr_vd)[,-1]), "_Vardict", sep="")
  
  meta_data=data.frame(snv_pass)[,1:3]
  meta_data[c(names_m2,names_f,names_vs,names_vd)]=NA
  #do not merge when the 1st row is all NAs + last column PASSED/REJECT
  if ((rowSums(is.na(data.frame(mcols(gr_m2)[1,-1])))!=ncol(data.frame(mcols(gr_m2)[1,-1]))-1)==TRUE) {
    meta_data[Biostrings::match(gr_m2, snv_pass), names_m2]=data.frame(mcols(gr_m2)[,-1])
  }
  if ((rowSums(is.na(data.frame(mcols(gr_f)[1,-1])))!=ncol(data.frame(mcols(gr_f)[1,-1]))-1)==TRUE) {
    meta_data[Biostrings::match(gr_f, snv_pass), names_f]=data.frame(mcols(gr_f)[,-1])
  }
  if ((rowSums(is.na(data.frame(mcols(gr_vs)[1,-1])))!=ncol(data.frame(mcols(gr_vs)[1,-1]))-1)==TRUE) {
    meta_data[Biostrings::match(gr_vs, snv_pass), names_vs]=data.frame(mcols(gr_vs)[,-1])
  }
  if ((rowSums(is.na(data.frame(mcols(gr_vd)[1,-1])))!=ncol(data.frame(mcols(gr_vd)[1,-1]))-1)==TRUE) {
    meta_data[Biostrings::match(gr_vd, snv_pass), names_vd]=data.frame(mcols(gr_vd)[,-1])
  }
  end.time=Sys.time()
  print(end.time-start.time)
  
  
  print("formating")
  start.time=Sys.time()
  
  #pick out cases with some indels
  meta_indel_cases <- subset(meta_data, nchar(REF_Mutect2)!=1 |
                               nchar(REF_Freebayes)!=1 |
                               nchar(REF_Varscan)!=1 |
                               nchar(REF_Vardict)!=1 |
                               nchar(ALT_Mutect2)!=1 |
                               nchar(ALT_Freebayes)!=1 |
                               nchar(ALT_Varscan)!=1 |
                               nchar(ALT_Vardict)!=1)
  
  # extract reference allele
  ref=meta_data[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict")] 
  ref.ind=which(!is.na(ref), arr.ind=T) # find non-NA alleles
  ref.ind=ref.ind[order(ref.ind[,1]),] #order array indices by row number
  ref.ind=ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
  meta_data$REF=ref[ref.ind]
  
  # extract alternate allele
  alt=meta_data[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict")]
  suppressWarnings(alt$ALT_Mutect2[(!is.na(alt$ALT_Mutect2) & nchar(alt$ALT_Mutect2)!=1)] <- substrRight(alt$ALT_Mutect2[(!is.na(alt$ALT_Mutect2) & nchar(alt$ALT_Mutect2)!=1)],1))
  suppressWarnings(alt$ALT_Freebayes[(!is.na(alt$ALT_Freebayes) & nchar(alt$ALT_Freebayes)!=1)] <- substrRight(alt$ALT_Freebayes[(!is.na(alt$ALT_Freebayes) & nchar(alt$ALT_Freebayes)!=1)],1))
  suppressWarnings(alt$ALT_Varscan[(!is.na(alt$ALT_Varscan) & nchar(alt$ALT_Varscan)!=1)] <- substrRight(alt$ALT_Varscan[(!is.na(alt$ALT_Varscan) & nchar(alt$ALT_Varscan)!=1)],1))
  suppressWarnings(alt$ALT_Vardict[(!is.na(alt$ALT_Vardict) & nchar(alt$ALT_Vardict)!=1)] <- substrRight(alt$ALT_Vardict[(!is.na(alt$ALT_Vardict) & nchar(alt$ALT_Vardict)!=1)],1))
  # alt$ALT_Freebayes[nchar(alt$ALT_Freebayes)!=1] <- NA
  # alt$ALT_Varscan[nchar(alt$ALT_Varscan)!=1] <- NA
  # alt$ALT_Vardict[nchar(alt$ALT_Vardict)!=1] <- NA
  alt.ind=which(!is.na(alt), arr.ind=T) #find all non-NA alleles
  alt.ind=alt.ind[order(alt.ind[,1]),]  #order array indices by row number
  alt.ind=alt.ind[!duplicated(alt.ind[,1]),]# take the first non-NA allele of each row
  meta_data$ALT=alt[alt.ind]  
  
  # calculate alt allele frequency
  meta_data$T_refDepth_Varscan <- meta_data$T_totalDepth_Varscan-meta_data$T_altDepth_Varscan
  meta_data$N_refDepth_Varscan <- meta_data$N_totalDepth_Varscan-meta_data$N_altDepth_Varscan
  meta_data$T_altDepth_Freebayes <- meta_data$T_totalDepth_Freebayes-meta_data$RO_Freebayes
  meta_data$T_refDepth_Freebayes <- meta_data$RO_Freebayes
  meta_data$AF_Freebayes <- meta_data$T_altDepth_Freebayes/meta_data$T_totalDepth_Freebayes
  # meta_data$AF_Freebayes[meta_data$AF_Freebayes == 0] <- NA
  
  af=meta_data[,c("AF_Mutect2", "FREQ_Varscan", "AF_Vardict", "AF_Freebayes")]
  af.ind=which(!is.na(af), arr.ind=T) #find all non-NA allele freq
  af.ind=af.ind[order(af.ind[,1]),]  #order array indices by row number
  af.ind=af.ind[!duplicated(af.ind[,1]),]# take the first non-NA allele freq of each row
  meta_data$AF=af[af.ind]
  meta_data$AF[is.na(meta_data$AF)] <- 0
  
  # meta_data$AF <- rowMeans(meta_data[, c("AF_Mutect2", "AF_Freebayes", "FREQ_Varscan", "AF_Vardict")],na.rm = TRUE)
  
  # calculate mean depth
  meta_data$N_refDepth <- round(rowMeans(meta_data[, c("N_refDepth_Mutect2", "N_refDepth_Freebayes", "N_refDepth_Varscan", "N_refDepth_Vardict")],na.rm = TRUE))
  meta_data$N_altDepth <- round(rowMeans(meta_data[, c("N_altDepth_Mutect2", "N_altDepth_Freebayes", "N_altDepth_Varscan", "N_altDepth_Vardict")],na.rm = TRUE))
  meta_data$T_refDepth <- round(rowMeans(meta_data[, c("T_refDepth_Mutect2", "T_refDepth_Freebayes", "T_refDepth_Varscan", "T_refDepth_Vardict")],na.rm = TRUE))
  meta_data$T_altDepth <- round(rowMeans(meta_data[, c("T_altDepth_Mutect2", "T_altDepth_Freebayes", "T_altDepth_Varscan", "T_altDepth_Vardict")],na.rm = TRUE))
  
  meta_data$N_refDepth[is.nan(meta_data$N_refDepth)] <- 0
  meta_data$N_altDepth[is.nan(meta_data$N_altDepth)] <- 0
  meta_data$T_refDepth[is.nan(meta_data$T_refDepth)] <- 0
  meta_data$T_altDepth[is.nan(meta_data$T_altDepth)] <- 0
  
  # make filters logical
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 != "PASS"] <- FALSE
  meta_data$FILTER_Mutect2[is.na(meta_data$FILTER_Mutect2)] <- FALSE
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "PASS"] <- TRUE
  meta_data$FILTER_Mutect2 <- as.logical(meta_data$FILTER_Mutect2)
  
  meta_data$FILTER_Freebayes[meta_data$FILTER_Freebayes != "PASS"] <- FALSE
  meta_data$FILTER_Freebayes[is.na(meta_data$FILTER_Freebayes)] <- FALSE
  meta_data$FILTER_Freebayes[meta_data$FILTER_Freebayes == "PASS"] <- TRUE
  meta_data$FILTER_Freebayes <- as.logical(meta_data$FILTER_Freebayes)
  
  meta_data$FILTER_Vardict[meta_data$FILTER_Vardict != "PASS"] <- FALSE
  meta_data$FILTER_Vardict[is.na(meta_data$FILTER_Vardict)] <- FALSE
  meta_data$FILTER_Vardict[meta_data$FILTER_Vardict == "PASS"] <- TRUE
  meta_data$FILTER_Vardict <- as.logical(meta_data$FILTER_Vardict)
  
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "SpvFreq"] <- "PASS"
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "REJECT;SpvFreq"] <- "PASS"
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan != "PASS"] <- FALSE
  meta_data$FILTER_Varscan[is.na(meta_data$FILTER_Varscan)] <- FALSE
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "PASS"] <- TRUE
  meta_data$FILTER_Varscan <- as.logical(meta_data$FILTER_Varscan)
  
  # get all 4 ref/alternates
  meta_data$REF_MFVdVs<-paste(meta_data$REF_Mutect2,meta_data$REF_Freebayes,meta_data$REF_Vardict,meta_data$REF_Varscan, sep ="/")
  meta_data$ALT_MFVdVs	<-paste(meta_data$ALT_Mutect2,meta_data$ALT_Freebayes,meta_data$ALT_Vardict,meta_data$ALT_Varscan, sep ="/")
  end.time=Sys.time()
  print(end.time-start.time)
  
  # sample name
  # meta_data$Sample_Name <- substrLeft(sampleid.t, 2)
  meta_data$Sample_Name <- gsub(t.label,'',sampleid.t)
  
  #insert missing feature(s)
  # meta_data$BaseQRankSum_Mutect2 <- NA
  meta_data$relcov <- (meta_data$T_refDepth+meta_data$T_altDepth+meta_data$N_refDepth+meta_data$N_altDepth)/median(meta_data$T_refDepth+meta_data$T_altDepth+meta_data$N_refDepth+meta_data$N_altDepth)
  
  #numerical
  
  #add "m2_BaseQRankSum" and "m2_ClippingRankSum" for long-term support
  meta_data$BaseQRankSum_Mutect2 = NA
  meta_data$ClippingRankSum_Mutect2 = NA
  
  meta_data[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
               "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
               "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
               "SSF_Vardict","MSI_Vardict","SOR_Vardict",
               "relcov")] <- lapply (meta_data[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                                                                                   "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                                                                                   "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                                                                                   "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                                                  "relcov")], as.numeric)
  
  #keep important columns and rename columns
  parse_snv <- meta_data[,c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                            "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                            "MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                            "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                            "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                            "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                            "AF",
                            "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")]
  colnames(parse_snv) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                           "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                           "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","m2_FS","m2_ReadPosRankSum","m2_BaseQRankSum","m2_ClippingRankSum",
                           "f_MQM","f_MQMR","f_GTI","f_LEN","f_ODDS",
                           "vs_SSC","vs_SPV","vs_GPV","vs_SS",
                           "vd_SSF","vd_MSI","vd_SOR",
                           "Alt_Allele_Freq",
                           "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
  
  } else { 
    parse_snv=NULL
    print('There are no SNVs in the region.')}
  
  ### continue from here to format the output
  
  print("Parsing INDEL")
  
  #### Parsing indels ####
  # get passed indel calls from all callers
  indel_pass=c(rowRanges(vcf_m2[pass_m2&!snv_m2]), rowRanges(vcf_f[pass_f&!snv_f]),  rowRanges(vcf_vs[pass_vs&!snv_vs]), rowRanges(vcf_vd[pass_vd&!snv_vd]),ignore.mcols=T)
  indel_pass= unique(indel_pass)
  
  if(length(indel_pass)!=0) {
  
  # get the row numbers of calls that are passed by at least 1 caller in each vcf
  indel.m2.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_m2), type="equal"))
  indel.f.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_f), type="equal"))
  indel.vs.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_vs), type="equal"))
  indel.vd.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_vd), type="equal"))
  
  print("extracting meta data from VRanges for indels")
  start.time=Sys.time() 
  
  if (length(indel.m2.index)!=0) {
    
  # convert vcf to VRanges then to Granges, keep metadata columns
  vr_m2<- suppressWarnings(as(vcf_m2[indel.m2.index], "VRanges"))
  #mcols(vr_m2)=vr_m2[,c("MQ","MQRankSum","TLOD","NLOD","AF")]
  vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
  gr_m2=GRanges(vr_m2[[sampleid.t]])
  mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))
  
  gr_m2 <- unique(gr_m2)
  gr_m2$FILTER=vcf_m2[indel.m2.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from MuTect2.")
    vr_m2<- as(vcf_m2, "VRanges")
    vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
    gr_m2=GRanges(vr_m2[[sampleid.t]])
    mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                                T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                                N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_m2 <- unique(gr_m2)
    gr_m2$FILTER=vcf_m2@fixed$FILTER
    
    gr_m2 <- gr_m2[1,]
    gr_m2@elementMetadata@listData$QUAL <- NA
    gr_m2@elementMetadata@listData$MQ <- NA
    gr_m2@elementMetadata@listData$MQRankSum <- NA
    gr_m2@elementMetadata@listData$NLOD <- NA
    gr_m2@elementMetadata@listData$TLOD <- NA
    # gr_m2@elementMetadata@listData$ClippingRankSum <- NA
    gr_m2@elementMetadata@listData$ReadPosRankSum <- NA
    gr_m2@elementMetadata@listData$FS <- NA
    gr_m2@elementMetadata@listData$AF <- NA
    
    
  }
  
  if (length(indel.f.index)!=0) {
    
  vr_f<- suppressWarnings(as(vcf_f[indel.f.index], "VRanges"))
  #mcols(vr_f)=vr_f[,c("MQM","MQMR")]
  vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
  gr_f=GRanges(vr_f[[sampleid.t]])
  mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                            T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                            N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
  
  gr_f <- unique(gr_f)
  gr_f$FILTER=vcf_f[indel.f.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from FreeBayes.")
    vr_f<- suppressWarnings(as(vcf_f, "VRanges"))
    vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
    gr_f=GRanges(vr_f[[sampleid.t]])
    mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                              T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                              N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_f <- unique(gr_f)
    gr_f$FILTER=vcf_f@fixed$FILTER
    
    gr_f <- gr_f[1,]
    gr_f@elementMetadata@listData$QUAL <- NA
    gr_f@elementMetadata@listData$LEN <- NA
    gr_f@elementMetadata@listData$MQM <- NA
    gr_f@elementMetadata@listData$MQMR <- NA
    gr_f@elementMetadata@listData$ODDS <- NA
    gr_f@elementMetadata@listData$GTI <- NA
    gr_f@elementMetadata@listData$RO <- NA
    
  }
  
  if (length(indel.vs.index)!=0) {
    
  vr_vs<- suppressWarnings(as(vcf_vs[indel.vs.index], "VRanges"))
  #mcols(vr_vs)=vr_vs[,c("SSC","SPV", "FREQ")]
  vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
  gr_vs=GRanges(vr_vs[[sampleid.t]])
  mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vs <- unique(gr_vs)
  gr_vs$FILTER=vcf_vs[indel.vs.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from VarScan.")
    vr_vs<- suppressWarnings(as(vcf_vs, "VRanges"))
    vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
    gr_vs=GRanges(vr_vs[[sampleid.t]])
    mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                                T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                                N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_vs <- unique(gr_vs)
    gr_vs$FILTER=vcf_vs@fixed$FILTER
    
    gr_vs <- gr_vs[1,]
    gr_vs@elementMetadata@listData$QUAL <- NA
    gr_vs@elementMetadata@listData$SSC <- NA
    gr_vs@elementMetadata@listData$SPV <- NA
    gr_vs@elementMetadata@listData$GPV <- NA
    gr_vs@elementMetadata@listData$SS <- NA
    gr_vs@elementMetadata@listData$FREQ <- NA
    
  }
  
  if (length(indel.vd.index)!=0) {
    
  vr_vd<- suppressWarnings(as(vcf_vd[indel.vd.index], "VRanges"))
  #mcols(vr_vd)=vr_vd[,c("SSF","SOR","MSI", "AF")]
  vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
  gr_vd=GRanges(vr_vd[[sampleid.t]])
  mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
  
  gr_vd <- unique(gr_vd)
  gr_vd$FILTER=vcf_vd[indel.vd.index]@fixed$FILTER
  
  } else {
    print("Warning: There are no passed INDELs from VarDict.")
    vr_vd<- suppressWarnings(as(vcf_vd, "VRanges"))
    vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
    gr_vd=GRanges(vr_vd[[sampleid.t]])
    mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=NA, ALT=NA, T_totalDepth=NA, 
                                                T_refDepth=NA, T_altDepth=NA,N_totalDepth=NA, 
                                                N_refDepth=NA, N_altDepth=NA, stringsAsFactors=F))
    gr_vd <- unique(gr_vd)
    gr_vd$FILTER=vcf_vd@fixed$FILTER
    
    gr_vd <- gr_vd[1,]
    gr_vd@elementMetadata@listData$QUAL <- NA
    gr_vd@elementMetadata@listData$SSF <- NA
    gr_vd@elementMetadata@listData$MSI <- NA
    gr_vd@elementMetadata@listData$SOR <- NA
    gr_vd@elementMetadata@listData$AF <- NA
    
    
  }
  
  end.time=Sys.time() 
  print(end.time-start.time)
  
  print("merge 4 vcfs and meta data")
  start.time=Sys.time()
  
  ## merge 4 vcfs and meta data
  names_m2= paste(names(mcols(gr_m2)[,-1]), "_Mutect2", sep="")
  names_f= paste(names(mcols(gr_f)[,-1]), "_Freebayes", sep="")
  names_vs= paste(names(mcols(gr_vs)[,-1]), "_Varscan", sep="")
  names_vd= paste(names(mcols(gr_vd)[,-1]), "_Vardict", sep="")
  
  meta_indel=data.frame(indel_pass)[,1:3]
  meta_indel[c(names_m2,names_f,names_vs,names_vd)]=NA
  #do not merge when the 1st row is all NAs + last column PASSED/REJECT
  if ((rowSums(is.na(data.frame(mcols(gr_m2)[1,-1])))!=ncol(data.frame(mcols(gr_m2)[1,-1]))-1)==TRUE) {
    meta_indel[Biostrings::match(gr_m2, indel_pass), names_m2]=data.frame(mcols(gr_m2)[,-1])
  }
  if ((rowSums(is.na(data.frame(mcols(gr_f)[1,-1])))!=ncol(data.frame(mcols(gr_f)[1,-1]))-1)==TRUE) {
    meta_indel[Biostrings::match(gr_f, indel_pass), names_f]=data.frame(mcols(gr_f)[,-1])
  }
  if ((rowSums(is.na(data.frame(mcols(gr_vs)[1,-1])))!=ncol(data.frame(mcols(gr_vs)[1,-1]))-1)==TRUE) {
    meta_indel[Biostrings::match(gr_vs, indel_pass), names_vs]=data.frame(mcols(gr_vs)[,-1])
  }
  if ((rowSums(is.na(data.frame(mcols(gr_vd)[1,-1])))!=ncol(data.frame(mcols(gr_vd)[1,-1]))-1)==TRUE) {
    meta_indel[Biostrings::match(gr_vd, indel_pass), names_vd]=data.frame(mcols(gr_vd)[,-1])
  }
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
  ref <- meta_indel[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict")] 
  ref.nchar <- sapply(ref, nchar)
  ref.nchar <- as.data.frame(ref.nchar)
  ref.nchar$MAX <- apply(ref.nchar,MARGIN = 1,function(x) max(x,na.rm=TRUE))
  ref.max <- ref.nchar[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict")]
  ref.max[,c("REF_Mutect2")] <- ref.nchar[,c("REF_Mutect2")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Freebayes")] <- ref.max[,c("REF_Freebayes")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Varscan")] <- ref.max[,c("REF_Varscan")] - ref.nchar[,c("MAX")]
  ref.max[,c("REF_Vardict")] <- ref.max[,c("REF_Vardict")] - ref.nchar[,c("MAX")]
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
  alt <- meta_indel[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict")] 
  alt.nchar <- sapply(alt, nchar)
  alt.nchar <- as.data.frame(alt.nchar)
  alt.nchar$MAX <- apply(alt.nchar,MARGIN = 1,function(x) max(x,na.rm=TRUE))
  alt.max <- alt.nchar[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict")]
  alt.max[,c("ALT_Mutect2")] <- alt.nchar[,c("ALT_Mutect2")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Freebayes")] <- alt.max[,c("ALT_Freebayes")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Varscan")] <- alt.max[,c("ALT_Varscan")] - alt.nchar[,c("MAX")]
  alt.max[,c("ALT_Vardict")] <- alt.max[,c("ALT_Vardict")] - alt.nchar[,c("MAX")]
  alt.ind <- which(!is.na(alt.max)&(alt.max)==0, arr.ind=T) # find non-NA alleles
  alt.ind <- alt.ind[order(alt.ind[,1]),] #order array indices by row number
  alt.ind <- alt.ind[!duplicated(alt.ind[,1]),] #get first non-NA value for each row
  meta_indel$ALT <- alt[alt.ind]
  
  # calculate alt allele frequency
  meta_indel$T_refDepth_Varscan <- meta_indel$T_totalDepth_Varscan-meta_indel$T_altDepth_Varscan
  meta_indel$N_refDepth_Varscan <- meta_indel$N_totalDepth_Varscan-meta_indel$N_altDepth_Varscan
  meta_indel$T_altDepth_Freebayes <- meta_indel$T_totalDepth_Freebayes-meta_indel$RO_Freebayes
  meta_indel$T_refDepth_Freebayes <- meta_indel$RO_Freebayes
  meta_indel$AF_Freebayes <- meta_indel$T_altDepth_Freebayes/meta_indel$T_totalDepth_Freebayes
  
  af=meta_indel[,c("AF_Mutect2", "FREQ_Varscan", "AF_Vardict", "AF_Freebayes")]
  af.ind=which(!is.na(af), arr.ind=T) #find all non-NA allele freq
  af.ind=af.ind[order(af.ind[,1]),]  #order array indices by row number
  af.ind=af.ind[!duplicated(af.ind[,1]),]# take the first non-NA allele freq of each row
  meta_indel$AF=af[af.ind]  
  meta_indel$AF[is.na(meta_indel$AF)] <- 0
  
  # calculate mean depth
  meta_indel$N_refDepth <- round(rowMeans(meta_indel[, c("N_refDepth_Mutect2", "N_refDepth_Freebayes", "N_refDepth_Varscan", "N_refDepth_Vardict")],na.rm = TRUE))
  meta_indel$N_altDepth <- round(rowMeans(meta_indel[, c("N_altDepth_Mutect2", "N_altDepth_Freebayes", "N_altDepth_Varscan", "N_altDepth_Vardict")],na.rm = TRUE))
  meta_indel$T_refDepth <- round(rowMeans(meta_indel[, c("T_refDepth_Mutect2", "T_refDepth_Freebayes", "T_refDepth_Varscan", "T_refDepth_Vardict")],na.rm = TRUE))
  meta_indel$T_altDepth <- round(rowMeans(meta_indel[, c("T_altDepth_Mutect2", "T_altDepth_Freebayes", "T_altDepth_Varscan", "T_altDepth_Vardict")],na.rm = TRUE))

  meta_indel$N_refDepth[is.nan(meta_indel$N_refDepth)] <- 0
  meta_indel$N_altDepth[is.nan(meta_indel$N_altDepth)] <- 0
  meta_indel$T_refDepth[is.nan(meta_indel$T_refDepth)] <- 0
  meta_indel$T_altDepth[is.nan(meta_indel$T_altDepth)] <- 0
  
  # make filters logical
  meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 != "PASS"] <- FALSE
  meta_indel$FILTER_Mutect2[is.na(meta_indel$FILTER_Mutect2)] <- FALSE
  meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 == "PASS"] <- TRUE
  meta_indel$FILTER_Mutect2 <- as.logical(meta_indel$FILTER_Mutect2)
  
  meta_indel$FILTER_Freebayes[meta_indel$FILTER_Freebayes != "PASS"] <- FALSE
  meta_indel$FILTER_Freebayes[is.na(meta_indel$FILTER_Freebayes)] <- FALSE
  meta_indel$FILTER_Freebayes[meta_indel$FILTER_Freebayes == "PASS"] <- TRUE
  meta_indel$FILTER_Freebayes <- as.logical(meta_indel$FILTER_Freebayes)
  
  meta_indel$FILTER_Vardict[meta_indel$FILTER_Vardict != "PASS"] <- FALSE
  meta_indel$FILTER_Vardict[is.na(meta_indel$FILTER_Vardict)] <- FALSE
  meta_indel$FILTER_Vardict[meta_indel$FILTER_Vardict == "PASS"] <- TRUE
  meta_indel$FILTER_Vardict <- as.logical(meta_indel$FILTER_Vardict)
  
  meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan == "SpvFreq"] <- "PASS"
  meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan == "REJECT;SpvFreq"] <- "PASS"
  meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan != "PASS"] <- FALSE
  meta_indel$FILTER_Varscan[is.na(meta_indel$FILTER_Varscan)] <- FALSE
  meta_indel$FILTER_Varscan[meta_indel$FILTER_Varscan == "PASS"] <- TRUE
  meta_indel$FILTER_Varscan <- as.logical(meta_indel$FILTER_Varscan)
  
  # get all 4 ref/alternates
  meta_indel$REF_MFVdVs<-paste(meta_indel$REF_Mutect2,meta_indel$REF_Freebayes,meta_indel$REF_Vardict,meta_indel$REF_Varscan, sep ="/")
  meta_indel$ALT_MFVdVs	<-paste(meta_indel$ALT_Mutect2,meta_indel$ALT_Freebayes,meta_indel$ALT_Vardict,meta_indel$ALT_Varscan, sep ="/")
  
  # sample name
  # meta_indel$Sample_Name <- substrLeft(sampleid.t, 2)
  meta_indel$Sample_Name <- gsub(t.label,'',sampleid.t)
  
  #insert missing feature
  # meta_indel$BaseQRankSum_Mutect2 <- NA
  meta_indel$relcov <- (meta_indel$T_refDepth+meta_indel$T_altDepth+meta_indel$N_refDepth+meta_indel$N_altDepth)/median(meta_indel$T_refDepth+meta_indel$T_altDepth+meta_indel$N_refDepth+meta_indel$N_altDepth)
  
  #add "m2_BaseQRankSum" and "m2_ClippingRankSum" for long-term support
  meta_indel$BaseQRankSum_Mutect2 = NA
  meta_indel$ClippingRankSum_Mutect2 = NA
  
  #numerical
  meta_indel[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                "relcov")] <- lapply (meta_indel[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                                                                                     "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                                                                                     "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                                                                                     "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                                                    "relcov")], as.numeric)
  
  parse_indel <- meta_indel[,c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                               "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                               "MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                               "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                               "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                               "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                               "AF",
                               "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")]
  colnames(parse_indel) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                             "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                             "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","m2_FS","m2_ReadPosRankSum","m2_BaseQRankSum","m2_ClippingRankSum",
                             "f_MQM","f_MQMR","f_GTI","f_LEN","f_ODDS",
                             "vs_SSC","vs_SPV","vs_GPV","vs_SS",
                             "vd_SSF","vd_MSI","vd_SOR",
                             "Alt_Allele_Freq",
                             "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
  
  } else { 
    parse_indel=NULL
    print('There are no indels in the region.')}
  
  return(list(snv=parse_snv, indel=parse_indel))
}



