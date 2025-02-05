#' parsevcf
#'
#' Step 1 Initial parse and filtering
#'
#' @param x List object containing the four vcf.gz files from 
#' callers Strelka2, MuTect2, Freebayes, VarDict and varscan.  
#' 
#' @param t.label Default='-T'. Identify your tumour label based on your default vcf tumour sample.
#'   
#' @examples
#' 
#' 
#' @export
parsevcf_allfeaturesall = function(x, tbi, t.label=NULL){
    print('Parsing step...')
    print("Extracting genome coordinate of the calls passed by at least 1 caller")
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }

    start.time=Sys.time()
    
    svp_m<-ScanVcfParam(info=c("FS"), samples=suppressWarnings(scanVcfHeader(x[[1]])@samples), geno=c(""), fixed = c('REF', 'ALT', 'FILTER'))
    svp_f<-ScanVcfParam(info=c("GTI"),samples=suppressWarnings(scanVcfHeader(x[[2]])@samples), geno=c(""), fixed = c('REF', 'ALT', 'FILTER'))
    svp_vs<-ScanVcfParam(info=c("SSC"), samples=suppressWarnings(scanVcfHeader(x[[3]])@samples), geno=c(""), fixed = c('REF', 'ALT', 'FILTER'))
    svp_vd<-ScanVcfParam(info=c("SSF"), samples=suppressWarnings(scanVcfHeader(x[[4]])@samples), geno=c(""), fixed = c('REF', 'ALT', 'FILTER'))
    svp_s = ScanVcfParam(info=c("QSS"), samples=suppressWarnings(scanVcfHeader(x[[5]])@samples), geno=c(""), fixed = c('REF', 'ALT', 'FILTER'))

    
    ## read only chromosomes 1-22, X, Y and M
    vcf_m2<- suppressWarnings(readVcf(tbi[[1]], genome=seqinfo(scanVcfHeader(x[[1]])), svp_m))
    vcf_f<- suppressWarnings(readVcf(tbi[[2]], genome=seqinfo(scanVcfHeader(x[[2]])), svp_f))
    vcf_vs<- suppressWarnings(readVcf(tbi[[3]], genome=seqinfo(scanVcfHeader(x[[3]])), svp_vs))
    vcf_vd<- suppressWarnings(readVcf(tbi[[4]], genome=seqinfo(scanVcfHeader(x[[4]])), svp_vd))
    vcf_s2<- suppressWarnings(readVcf(tbi[[5]], genome=seqinfo(scanVcfHeader(x[[5]])), svp_s))
    rm(svp_f, svp_m, svp_vs, svp_vd, svp_s) 
    
    # remove duplicates from vcfs
    vcf_m2=vcf_m2[!duplicated(rowRanges(vcf_m2))]
    vcf_f=vcf_f[!duplicated(rowRanges(vcf_f))]
    vcf_vs=vcf_vs[!duplicated(rowRanges(vcf_vs))]
    vcf_vd=vcf_vd[!duplicated(rowRanges(vcf_vd))]
    vcf_s2=vcf_s2[!duplicated(rowRanges(vcf_s2))]
    
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
    
    # get rows that are SNVs
    snv_m2<-isSNV(vcf_m2, singleAltOnly=FALSE)
    snv_f<-isSNV(vcf_f, singleAltOnly=FALSE)
    snv_vs<-isSNV(vcf_vs, singleAltOnly=FALSE)
    snv_vd<-isSNV(vcf_vd, singleAltOnly=FALSE)  
    snv_s2<-isSNV(vcf_s2, singleAltOnly=FALSE)
    
    
    # extract passed calls from each caller
    pass_m2<- vcf_m2@fixed$FILTER=="PASS"# | vcf_m2@fixed$FILTER=="MinAF"
    pass_f<- vcf_f@fixed$FILTER=="PASS" 
    pass_vs<- vcf_vs@fixed$FILTER=="PASS"
    # pass_vs<- vcf_vs@fixed$FILTER=="PASS" | vcf_vs@fixed$FILTER=="SpvFreq" | vcf_vs@fixed$FILTER=="REJECT;SpvFreq" #SpvFreq not found in GATK4
    pass_vd<- vcf_vd@fixed$FILTER=="PASS"
    pass_s2<- vcf_s2@fixed$FILTER=="PASS" #| vcf_s2@fixed$FILTER=="MinAF"
    
    
    # get passed snv calls from all callers
    snv_pass=c(rowRanges(vcf_m2[pass_m2&snv_m2]), 
               rowRanges(vcf_f[pass_f&snv_f]),  
               rowRanges(vcf_vs[pass_vs&snv_vs]), 
               rowRanges(vcf_vd[pass_vd&snv_vd]),
               rowRanges(vcf_s2[pass_s2&snv_s2]),
               ignore.mcols=T)
    snv_pass= unique(snv_pass)
    
    # get passed indel calls from all callers
    indel_pass=c(rowRanges(vcf_m2[pass_m2&!snv_m2]), 
                 rowRanges(vcf_f[pass_f&!snv_f]),  
                 rowRanges(vcf_vs[pass_vs&!snv_vs]), 
                 rowRanges(vcf_vd[pass_vd&!snv_vd]),
                 rowRanges(vcf_s2[pass_s2&!snv_s2]),
                 ignore.mcols=T)
    indel_pass= unique(indel_pass)
    rm(pass_f, pass_vs, pass_vd, pass_s2, pass_m2, snv_f, snv_m2, snv_s2, snv_vd, snv_vs, vcf_f, vcf_m2, vcf_s2, vcf_vd, vcf_vs) 
    end.time = Sys.time()
    print(end.time - start.time)
    
    print("Reading vcf and extracting calls passed by at least 1 caller")

    start.time=Sys.time() 
    
    print('Pasing and extracting data from Mutect2')
    #scan vcf param
    svp_m<-ScanVcfParam(info=c("FS","MQ","MQRankSum","NLOD","ReadPosRankSum","TLOD"), samples=suppressWarnings(scanVcfHeader(x[[1]])@samples), geno=c("AD", "AF", "DP"))
    #read the vcf
    vcf_m2<- suppressWarnings(readVcf(tbi[[1]], genome=seqinfo(scanVcfHeader(x[[1]])), svp_m))
    rm(svp_m)
    #remove duplicated
    vcf_m2=vcf_m2[!duplicated(rowRanges(vcf_m2))]
    #get the column names
    m2.cols = c(names(vcf_m2@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
    names_m2= paste(m2.cols, "_Mutect2", sep="")
    
    #get the position of calls that passed at least 1 callers.
    snv.m2.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_m2), type="equal"))
    indel.m2.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_m2), type="equal"))
    
    #mutect2 SNV
    if (length(snv.m2.index)!=0) {
      m2.na.snv=F
      
      # convert vcf to VRanges then to Granges, keep metadata columns
      vr_m2<- as(vcf_m2[snv.m2.index], "VRanges")
      vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
      gr_m2=GRanges(vr_m2[[sampleid.t]])
      mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_m2) 
      gr_m2 <- unique(gr_m2)
      gr_m2$FILTER=vcf_m2[snv.m2.index]@fixed$FILTER
      gr_m2_snv = gr_m2
      rm(snv.m2.index, gr_m2)
    } else {
      print("Warning: There are no passed SNVs from MuTect2.")
      m2.na.snv=T
    }
    
    #mutect2 indel
    if (length(indel.m2.index)!=0) {
      m2.na=F  
      # convert vcf to VRanges then to Granges, keep metadata columns
      vr_m2<- suppressWarnings(as(vcf_m2[indel.m2.index], "VRanges"))
      vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
      gr_m2=GRanges(vr_m2[[sampleid.t]])
      mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_m2)
      gr_m2 <- unique(gr_m2)
      gr_m2$FILTER=vcf_m2[indel.m2.index]@fixed$FILTER
      rm(indel.m2.index) 
    } else {
      print("Warning: There are no passed INDELs from MuTect2.")
      m2.na=T
    }
    rm(vcf_m2)
    
    print('Pasing and extracting data from FreeBayes')
    #scan vcf param
    svp_f<-ScanVcfParam(info=c("GTI","LEN","MQM","MQMR","ODDS"),samples=suppressWarnings(scanVcfHeader(x[[2]])@samples), geno=c("RO","DP"))
    #read the vcf
    vcf_f<- suppressWarnings(readVcf(tbi[[2]], genome=seqinfo(scanVcfHeader(x[[2]])), svp_f))
    rm(svp_f)
    #remove duplicated
    vcf_f=vcf_f[!duplicated(rowRanges(vcf_f))]
    #get the column names
    f.cols = c(names(vcf_f@info@listData),'RO','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
    names_f= paste(f.cols, "_Freebayes", sep="")
    #get the position of calls that passed at least 1 callers.
    snv.f.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_f), type="equal"))
    indel.f.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_f), type="equal"))
    
    if (length(snv.f.index)!=0) {
      f.na.snv=F  
      vr_f<- suppressWarnings(as(vcf_f[snv.f.index], "VRanges"))
      vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
      gr_f=GRanges(vr_f[[sampleid.t]])
      mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_f) 
      gr_f <- unique(gr_f)
      gr_f$FILTER=vcf_f[snv.f.index]@fixed$FILTER
      gr_f_snv = gr_f		
      rm(snv.f.index, gr_f)  
    } else {
      print("Warning: There are no passed SNVs from FreeBayes.")
      f.na.snv=T
    }
    
    
    
    #Freebayes indel
    if (length(indel.f.index)!=0) {
      f.na=F  
      vr_f<- suppressWarnings(as(vcf_f[indel.f.index], "VRanges"))
      vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
      gr_f=GRanges(vr_f[[sampleid.t]])
      mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_f)
      gr_f <- unique(gr_f)
      gr_f$FILTER=vcf_f[indel.f.index]@fixed$FILTER
      rm(indel.f.index) 
    } else {
      print("Warning: There are no passed INDELs from FreeBayes.")
      f.na=T
    }
    
    rm(vcf_f)
    
    print('Pasing and extracting data from VarScan')
    #scan vcf param
    svp_vs<-ScanVcfParam(info=c("SSC","GPV","SS","SPV"), samples=suppressWarnings(scanVcfHeader(x[[3]])@samples), geno=c("AD","FREQ", "DP","RD"))
    #read the vcf
    vcf_vs<- suppressWarnings(readVcf(tbi[[3]], genome=seqinfo(scanVcfHeader(x[[3]])), svp_vs))
    rm(svp_vs)
    #remove duplicated
    vcf_vs=vcf_vs[!duplicated(rowRanges(vcf_vs))]
    #get the column names
    vs.cols = c(names(vcf_vs@info@listData),'FREQ','RD','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
    names_vs= paste(vs.cols, "_Varscan", sep="")
    #get the position of calls that passed at least 1 callers.
    snv.vs.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vs), type="equal"))
    indel.vs.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_vs), type="equal"))
    
    #For varscan SNV
    if (length(snv.vs.index)!=0) {
      vs.na.snv=F  
      vr_vs<- suppressWarnings(as(vcf_vs[snv.vs.index], "VRanges"))
      vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
      gr_vs=GRanges(vr_vs[[sampleid.t]])
      mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
      
      rm(vr_vs)
      gr_vs <- unique(gr_vs)
      gr_vs$FILTER=vcf_vs[snv.vs.index]@fixed$FILTER
      gr_vs_snv = gr_vs
      rm(snv.vs.index, gr_vs)
    } else {
      print("Warning: There are no passed SNVs from VarScan.")
      vs.na.snv=T
    }
    
    #For Varscan Indel
    if (length(indel.vs.index)!=0) {
      vs.na=F  
      vr_vs<- suppressWarnings(as(vcf_vs[indel.vs.index], "VRanges"))
      vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
      gr_vs=GRanges(vr_vs[[sampleid.t]])
      mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_vs)
      gr_vs <- unique(gr_vs)
      gr_vs$FILTER=vcf_vs[indel.vs.index]@fixed$FILTER
      rm(indel.vs.index)
    } else {
      print("Warning: There are no passed INDELs from VarScan.")
      vs.na=T
    }
    
    rm(vcf_vs)
    
    print("Pasing and extracting data from VarDict")
    
    #scan vcf param
    svp_vd<-ScanVcfParam(info=c("SSF","MSI","SOR"), samples=suppressWarnings(scanVcfHeader(x[[4]])@samples), geno=c("AD", "AF", "DP"))
    #read the vcf
    vcf_vd<- suppressWarnings(readVcf(tbi[[4]], genome=seqinfo(scanVcfHeader(x[[4]])), svp_vd))
    rm(svp_vd)
    #remove duplicated
    vcf_vd=vcf_vd[!duplicated(rowRanges(vcf_vd))]
    #get the column names
    vd.cols = c(names(vcf_vd@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
    names_vd= paste(vd.cols, "_Vardict", sep="")
    #get the position of calls that passed at least 1 callers.
    snv.vd.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vd), type="equal"))
    indel.vd.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_vd), type="equal"))
    
    
    if (length(snv.vd.index)!=0) {
      vd.na.snv = F  
      vr_vd<- suppressWarnings(as(vcf_vd[snv.vd.index], "VRanges"))
      vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
      gr_vd=GRanges(vr_vd[[sampleid.t]])
      mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_vd)
      gr_vd <- unique(gr_vd)
      gr_vd$FILTER=vcf_vd[snv.vd.index]@fixed$FILTER
      gr_vd_snv = gr_vd	 
      rm(snv.vd.index, gr_vd)
    } else {
      print("Warning: There are no passed SNVs from VarDict.")
      vd.na.snv=T
    }
    
    #Vardict indel
    if (length(indel.vd.index)!=0) {
      vd.na=F  
      vr_vd<- suppressWarnings(as(vcf_vd[indel.vd.index], "VRanges"))
      vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
      gr_vd=GRanges(vr_vd[[sampleid.t]])
      mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_vd)
      gr_vd <- unique(gr_vd)
      gr_vd$FILTER=vcf_vd[indel.vd.index]@fixed$FILTER
      rm(indel.vd.index)
    } else {
      print("Warning: There are no passed INDELs from VarDict.")
      vd.na=T
    }
    
    rm(vcf_vd)
    
    print('Pasing and extracting data from Strelka2')
    #scan vcf param
    svp_s = ScanVcfParam(info=c("QSS","MQ","SomaticEVS","ReadPosRankSum"), samples=suppressWarnings(scanVcfHeader(x[[5]])@samples), geno=c("AD", "AF", "DP"))
    #read the vcf
    vcf_s2<- suppressWarnings(readVcf(tbi[[5]], genome=seqinfo(scanVcfHeader(x[[5]])), svp_s)) 
    rm(svp_s)
    #remove duplicated
    vcf_s2=vcf_s2[!duplicated(rowRanges(vcf_s2))]
    #get the column names
    s2.cols = c(names(vcf_s2@info@listData),'AF','REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
    names_s2= paste(s2.cols, "_Strelka2", sep="") 
    #get the position of calls that passed at least 1 callers.
    snv.s2.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_s2), type="equal"))
    indel.s2.index <- subjectHits(findOverlaps(indel_pass, rowRanges(vcf_s2), type="equal"))
    #For SNV
    if (length(snv.s2.index)!=0) {
      s2.na.snv=F
      
      # convert vcf to VRanges then to Granges, keep metadata columns
      vr_s2<- as(vcf_s2[snv.s2.index], "VRanges")
      vr_s2=GenomicRanges::split(vr_s2, vr_s2@sampleNames)
      gr_s2=GRanges(vr_s2[[sampleid.t]])
      mcols(gr_s2)=cbind(mcols(gr_s2), data.frame(REF=ref(vr_s2[[sampleid.t]]), ALT=alt(vr_s2[[sampleid.t]]), T_totalDepth=totalDepth(vr_s2[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_s2[[sampleid.t]]), T_altDepth=altDepth(vr_s2[[sampleid.t]]),N_totalDepth=totalDepth(vr_s2[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_s2[[sampleid.n]]), N_altDepth=altDepth(vr_s2[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_s2)
      gr_s2 <- unique(gr_s2)
      gr_s2$FILTER=vcf_s2[snv.s2.index]@fixed$FILTER
      gr_s2_snv = gr_s2		
      rm(snv.s2.index, gr_s2)
    } else {
      print("Warning: There are no passed SNVs from Strelka2.")
      s2.na.snv=T
    }
    
    #For indel
    if (length(indel.s2.index)!=0) {
      s2.na=F
      
      # convert vcf to VRanges then to Granges, keep metadata columns
      vr_s2<- suppressWarnings(as(vcf_s2[indel.s2.index], "VRanges"))
      vr_s2=GenomicRanges::split(vr_s2, vr_s2@sampleNames)
      gr_s2=GRanges(vr_s2[[sampleid.t]])
      mcols(gr_s2)=cbind(mcols(gr_s2), data.frame(REF=ref(vr_s2[[sampleid.t]]), ALT=alt(vr_s2[[sampleid.t]]), T_totalDepth=totalDepth(vr_s2[[sampleid.t]]), 
                                                  T_refDepth=refDepth(vr_s2[[sampleid.t]]), T_altDepth=altDepth(vr_s2[[sampleid.t]]),N_totalDepth=totalDepth(vr_s2[[sampleid.n]]), 
                                                  N_refDepth=refDepth(vr_s2[[sampleid.n]]), N_altDepth=altDepth(vr_s2[[sampleid.n]]), stringsAsFactors=F))
      rm(vr_s2)
      gr_s2 <- unique(gr_s2)		
      gr_s2$FILTER=vcf_s2[indel.s2.index]@fixed$FILTER
      rm(indel.s2.index)
    } else {
      print("Warning: There are no passed INDELs from Strelka2.")
      s2.na=T
    }
    rm(vcf_s2)
    
    end.time=Sys.time() 
    print(end.time-start.time)
    
    
    print("Merge 5 vcfs and meta data for SNV")
    start.time=Sys.time()
    
    meta_data=data.frame(snv_pass)[,1:3]
    
    #change meta_data from data.frame to data.table to deal with large dataset.
    #print('change data fortmat to data.table')
    setDT(meta_data)
    
    #meta_data[,c(names_m2,names_f,names_vs,names_vd,names_s2):=NA] 
    #meta_data[c(names_m2,names_f,names_vs,names_vd,names_s2)]=NA
    
    #do not merge when index=0, caller.na=T
    if (m2.na.snv==F) {
      #print(m2.na)
      meta_data[Biostrings::match(gr_m2_snv, snv_pass), names_m2]=data.frame(mcols(gr_m2_snv)[,m2.cols])
      rm(gr_m2_snv)
      #meta_data[Biostrings::match(gr_m2, snv_pass), names_m2:=data.frame(mcols(gr_m2)[,m2.cols])]
      #print(colnames(mcols(gr_m2)))
      #print(m2.cols)
      #print(head(data.frame(mcols(gr_m2)[,m2.cols])))
    }else{
      meta_data[,names_m2] = NA

    }

    if (f.na.snv==F) {
      meta_data[Biostrings::match(gr_f_snv, snv_pass), names_f]=data.frame(mcols(gr_f_snv)[,f.cols])
      rm(gr_f_snv)
      #meta_data[Biostrings::match(gr_f, snv_pass), names_f:=data.frame(mcols(gr_f)[,f.cols])]
    }else{
      meta_data[, names_f] = NA
    }

    if (vs.na.snv==F) {
      meta_data[Biostrings::match(gr_vs_snv, snv_pass), names_vs]=data.frame(mcols(gr_vs_snv)[,vs.cols])
      rm(gr_vs_snv)
      #meta_data[Biostrings::match(gr_vs, snv_pass), names_vs:=data.frame(mcols(gr_vs)[,vs.cols])]
    }else{
      meta_data[, names_vs] = NA
    }

    if (vd.na.snv==F) {
      #meta_data[Biostrings::match(gr_vd, snv_pass), names_vd:=data.frame(mcols(gr_vd)[,vd.cols])]
      meta_data[Biostrings::match(gr_vd_snv, snv_pass), names_vd]=data.frame(mcols(gr_vd_snv)[,vd.cols])
      rm(gr_vd_snv)
    }else{
      meta_data[, names_vd] = NA 
    }

    if (s2.na.snv==F) {
      meta_data[Biostrings::match(gr_s2_snv, snv_pass), names_s2]=data.frame(mcols(gr_s2_snv)[,s2.cols])
      # meta_data[Biostrings::match(gr_s2, snv_pass), names_s2:=data.frame(mcols(gr_s2)[,s2.cols])]
      rm(gr_s2_snv)
    }else{
      meta_data[, names_s2] = NA
    }
    #print(meta_data) 
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
    #change the format of meta_data from data.table back to data.frame for the follwing calculation.
    # meta_data = as.data.frame(meta_data)
    rm(snv_pass)
    end.time=Sys.time()
    print(end.time-start.time)
    
    print("Formating for SNV")
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
    ref = as.data.frame(ref)
    ref.ind=which(!is.na(ref), arr.ind=T) # find non-NA alleles
    ref.ind=ref.ind[order(ref.ind[,1]),] #order array indices by row number
    ref.ind=ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
    meta_data$REF=ref[ref.ind]
    rm(ref, ref.ind)
    #print('ref extraction done')
    #print(class(meta_data))
    
    # extract alternate allele
    alt=meta_data[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict", "ALT_Strelka2")]
    alt = as.data.frame(alt)
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
    rm(alt, alt.ind)
    #print('alt extraction done')  
    #print(class(meta_data))
    
    # check for list objects in columns
    #meta_data = as.data.frame(meta_data)
    # for (i in 1:ncol(meta_data)) {
    #   if(class(meta_data[[i]])=='list'){
    #     #meta_data[,i][sapply(meta_data[,i], is.null)] <- NA
    #     meta_data[,..i] = unlist(meta_data[,..i])
    #   }
    # }
    #fwrite(meta_data, 'meta_data_fortest.csv')
    #meta_data[, names(meta_data) := lapply(.SD, unlist)]
    for (i in colnames(meta_data) ){
      #print(class(meta_data[[i]]))
      
      if(class(meta_data[[i]])=='AsIs'){
        #meta_data[,i][sapply(meta_data[,i], is.null)] <- NA
        j1 <- which(lengths(meta_data[[i]]) == 0) 
        set(meta_data, i = j1, j = i, value = list(NA)) 
        meta_data[[i]] = unlist(meta_data[[i]])
      }
    }
    
    # calculate alt allele frequency
    # print('calculate alt allele frequency')
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
    
    #print('catch vcf caller mapping errors')
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
    #print('calculate mean depth')
    meta_data$N_refDepth <- round(rowMeans(meta_data[, c("N_refDepth_Mutect2", "N_refDepth_Freebayes", "N_refDepth_Varscan", "N_refDepth_Vardict","N_refDepth_Strelka2")],na.rm = TRUE))
    meta_data$N_altDepth <- round(rowMeans(meta_data[, c("N_altDepth_Mutect2", "N_altDepth_Freebayes", "N_altDepth_Varscan", "N_altDepth_Vardict","N_altDepth_Strelka2")],na.rm = TRUE))
    meta_data$T_refDepth <- round(rowMeans(meta_data[, c("T_refDepth_Mutect2", "T_refDepth_Freebayes", "T_refDepth_Varscan", "T_refDepth_Vardict","T_refDepth_Strelka2")],na.rm = TRUE))
    meta_data$T_altDepth <- round(rowMeans(meta_data[, c("T_altDepth_Mutect2", "T_altDepth_Freebayes", "T_altDepth_Varscan", "T_altDepth_Vardict","T_altDepth_Strelka2")],na.rm = TRUE))
    
    meta_data$N_refDepth[is.nan(meta_data$N_refDepth)] <- 0
    meta_data$N_altDepth[is.nan(meta_data$N_altDepth)] <- 0
    meta_data$T_refDepth[is.nan(meta_data$T_refDepth)] <- 0
    meta_data$T_altDepth[is.nan(meta_data$T_altDepth)] <- 0
    
    #re-calculate mean AF
    #print('re-calculate mean AF')
    meta_data$AF = round(meta_data$T_altDepth/(meta_data$T_altDepth+meta_data$T_refDepth), digits = 3)
    meta_data$AF[is.na(meta_data$AF)] <- 0
    
    if(class(meta_data$AF) %in% c('list', 'AsIs')) {
      meta_data$AF = as.numeric(meta_data$AF)
    }
    
    # make filters logical
    #print('make filters logical')
    meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 != "PASS"] <- FALSE
    meta_data$FILTER_Mutect2[is.na(meta_data$FILTER_Mutect2)] <- FALSE
    meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "PASS"] <- TRUE
    # meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "MinAF"] <- TRUE
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
    
    meta_data$FILTER_Strelka2[meta_data$FILTER_Strelka2 != "PASS"] <- FALSE
    meta_data$FILTER_Strelka2[is.na(meta_data$FILTER_Strelka2)] <- FALSE
    meta_data$FILTER_Strelka2[meta_data$FILTER_Strelka2 == "PASS"] <- TRUE
    # meta_data$FILTER_Strelka2[meta_data$FILTER_Strelka2 == "MinAF"] <- TRUE
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
    #setDT(meta_data)
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
    
    cols_to_convert <- c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                         "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                         "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                         "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                         "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                         "AF",
                         "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
    
    #  meta_data[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]
    #  meta_data[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
    #               "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
    #               "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
    #               "SSF_Vardict","MSI_Vardict","SOR_Vardict",
    #               "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
    #               "AF",
    #               "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")] <- lapply (meta_data[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
    #                                                                                                      "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
    #                                                                                                      "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
    #                                                                                                      "SSF_Vardict","MSI_Vardict","SOR_Vardict",
    #                                                                                                      "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
    #                                                                                                      "AF",
    #                                                                                                      "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")], as.numeric)   
    #
    #
    for (i in colnames(meta_data) ){
      #print(class(meta_data[[i]]))

      if(class(meta_data[[i]]) %in% c('AsIs', 'list')){
        #meta_data[,i][sapply(meta_data[,i], is.null)] <- NA
        j1 <- which(lengths(meta_data[[i]]) == 0)
        set(meta_data, i = j1, j = i, value = list(NA))
        meta_data[[i]] = unlist(meta_data[[i]])
      }
    }
 
   
    ind <- match(cols_to_convert , names(meta_data))
    

    for (i in seq_along(ind)) {
     # print(i)
     # print(meta_data[[ind[i]]])    
      set(meta_data, NULL, ind[i], as.numeric(meta_data[[ind[i]]]))
    }  
    rm(ind) 
    # meta_data$AF = round(meta_data$AF,digits=3)
    
    #keep important columns and rename columns
    cols_to_select <- c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                        "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                        "MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                        "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                        "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                        "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                        "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                        "AF",
                        "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
    
    meta_data <- setDT(unclass(meta_data)[cols_to_select])[]
    parse_snv <- meta_data
    rm(meta_data) 
    
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
    parse_snv <- setkey(parse_snv, Chr, START_POS_REF)
    
    # Remove duplicated rows
    parse_snv <- unique(parse_snv)
    
    ### continue from here to format the output
    print("Merge 5 vcfs and meta data for INDEL")
    start.time=Sys.time()
    
    meta_indel=data.frame(indel_pass)[,1:3]
    setDT(meta_indel)
    #meta_indel[c(names_m2,names_f,names_vs,names_vd,names_s2)]=NA
    
    #do not merge when index=0, caller.na=T
    if (m2.na==F) {
      meta_indel[Biostrings::match(gr_m2, indel_pass), names_m2]=data.frame(mcols(gr_m2)[,m2.cols])
      rm(gr_m2)
    }else{
      meta_indel[,names_m2] = NA
    }

    if (f.na==F) {
      meta_indel[Biostrings::match(gr_f, indel_pass), names_f]=data.frame(mcols(gr_f)[,f.cols])
      rm(gr_f)  
    }else{
      meta_indel[, names_f] = NA	   
    }

    if (vs.na==F) {
      meta_indel[Biostrings::match(gr_vs, indel_pass), names_vs]=data.frame(mcols(gr_vs)[,vs.cols])
      rm(gr_vs)  
    }else{
      meta_indel[, names_vs] = NA
    }

    if (vd.na==F) {
      meta_indel[Biostrings::match(gr_vd, indel_pass), names_vd]=data.frame(mcols(gr_vd)[,vd.cols])
      rm(gr_vd)
    }else{
      meta_indel[, names_vd] = NA
    }

    if (s2.na==F) {
      meta_indel[Biostrings::match(gr_s2, indel_pass), names_s2]=data.frame(mcols(gr_s2)[,s2.cols])
      rm(gr_s2)
    }else{
     meta_indel[, names_s2] = NA
    }
    
    rm(indel_pass) 
    
    end.time=Sys.time()
    print(end.time-start.time)
    
    for (i in colnames(meta_indel) ){
      #print(class(meta_indel[[i]]))
      
      if(class(meta_indel[[i]]) %in% c('AsIs', 'list')){
        #meta_indel[,i][sapply(meta_indel[,i], is.null)] <- NA
        j1 <- which(lengths(meta_indel[[i]]) == 0) 
        set(meta_indel, i = j1, j = i, value = list(NA)) 
        meta_indel[[i]] = unlist(meta_indel[[i]])
      }
    }
    for (i in colnames(meta_indel_cases) ){
      
      if(class(meta_indel_cases[[i]]) %in% c('AsIs', 'list')){
        #meta_indel[,i][sapply(meta_indel[,i], is.null)] <- NA
        j1 <- which(lengths(meta_indel_cases[[i]]) == 0) 
        set(meta_indel_cases, i = j1, j = i, value = list(NA)) 
        meta_indel_cases[[i]] = unlist(meta_indel_cases[[i]])
      }
    }
     
    #add in indel cases from snv matrix
   
    meta_indel <- dplyr::bind_rows(meta_indel,meta_indel_cases)
    rm(meta_indel_cases)  
 
    meta_indel <- unique(meta_indel) #remove all duplicate snv/indel cases
    
    print("Formating for INDEL")
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
    ref = as.data.frame(ref)
    meta_indel$REF <- ref[ref.ind]
    rm(ref, ref.ind, ref.max, ref.nchar) 
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
    alt = as.data.frame(alt)
    meta_indel$ALT <- alt[alt.ind]
    rm(alt, alt.ind, alt.max, alt.nchar)
    # check for list objects in columns
    for (i in colnames(meta_indel) ){
      #print(class(meta_indel[[i]]))
      
      if(class(meta_indel[[i]]) %in% c('AsIs', 'list')){
        #meta_indel[,i][sapply(meta_indel[,i], is.null)] <- NA
        j1 <- which(lengths(meta_indel[[i]]) == 0) 
        set(meta_indel, i = j1, j = i, value = list(NA)) 
        meta_indel[[i]] = unlist(meta_indel[[i]])
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
    
    if(class(meta_indel$AF) %in% c('list', 'AsIs')) {
      meta_indel$AF = as.numeric(meta_indel$AF)
    }
    
    # make filters logical
    meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 != "PASS"] <- FALSE
    meta_indel$FILTER_Mutect2[is.na(meta_indel$FILTER_Mutect2)] <- FALSE
    meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 == "PASS"] <- TRUE
    #  meta_indel$FILTER_Mutect2[meta_indel$FILTER_Mutect2 == "MinAF"] <- TRUE
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
    
    meta_indel$FILTER_Strelka2[meta_indel$FILTER_Strelka2 != "PASS"] <- FALSE
    meta_indel$FILTER_Strelka2[is.na(meta_indel$FILTER_Strelka2)] <- FALSE
    meta_indel$FILTER_Strelka2[meta_indel$FILTER_Strelka2 == "PASS"] <- TRUE
    # meta_indel$FILTER_Strelka2[meta_indel$FILTER_Strelka2 == "MinAF"] <- TRUE
    meta_indel$FILTER_Strelka2 <- as.logical(meta_indel$FILTER_Strelka2)
    
    # get all 4 ref/alternates
    meta_indel$REF_MFVdVs<-paste(meta_indel$REF_Mutect2,meta_indel$REF_Freebayes,meta_indel$REF_Vardict,meta_indel$REF_Varscan,meta_indel$REF_Strelka2, sep ="/")
    meta_indel$ALT_MFVdVs	<-paste(meta_indel$ALT_Mutect2,meta_indel$ALT_Freebayes,meta_indel$ALT_Vardict,meta_indel$ALT_Varscan,meta_indel$ALT_Strelka2, sep ="/")
    end.time=Sys.time()
    print(end.time-start.time)
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
    
    
    cols_to_convert <- c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                         "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                         "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                         "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                         "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                         "AF",
                         "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
    
    #meta_indel[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]
    #  meta_indel[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
    #               "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
    #               "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
    #               "SSF_Vardict","MSI_Vardict","SOR_Vardict",
    #               "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
    #               "AF",
    #               "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")] <- lapply (meta_indel[,c("MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
    #                                                                                                       "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
    #                                                                                                       "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
    #                                                                                                       "SSF_Vardict","MSI_Vardict","SOR_Vardict",
    #                                                                                                       "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
    #                                                                                                       "AF",
    #                                                                                                       "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")], as.numeric)
    #
    for (i in colnames(meta_indel) ){
      #print(class(meta_data[[i]]))

      if(class(meta_indel[[i]]) %in% c('AsIs', 'list')){
        #meta_data[,i][sapply(meta_data[,i], is.null)] <- NA
        j1 <- which(lengths(meta_indel[[i]]) == 0)
        set(meta_indel, i = j1, j = i, value = list(NA))
        meta_indel[[i]] = unlist(meta_indel[[i]])
      }
    }

    ind <- match(cols_to_convert , names(meta_indel))
    
    for (i in seq_along(ind)) {
     # print(i)
     # print(meta_indel[[ind[i]]])    
      set(meta_indel, NULL, ind[i], as.numeric(meta_indel[[ind[i]]]))
    }
    rm(ind)
    cols_to_select <- c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                        "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                        "MQ_Mutect2","MQRankSum_Mutect2","TLOD_Mutect2","NLOD_Mutect2","FS_Mutect2","ReadPosRankSum_Mutect2","BaseQRankSum_Mutect2","ClippingRankSum_Mutect2",
                        "MQM_Freebayes","MQMR_Freebayes","GTI_Freebayes","LEN_Freebayes","ODDS_Freebayes",
                        "SSC_Varscan","SPV_Varscan","GPV_Varscan","SS_Varscan",
                        "SSF_Vardict","MSI_Vardict","SOR_Vardict",
                        "QSS_Strelka2","MQ_Strelka2","SomaticEVS_Strelka2","ReadPosRankSum_Strelka2",
                        "AF",
                        "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
    
    
    meta_indel <- setDT(unclass(meta_indel)[cols_to_select])[]
    parse_indel <- meta_indel
    
    colnames(parse_indel) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                               "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan","FILTER_Strelka2",
                               "m2_MQ","m2_MQRankSum","m2_TLOD","m2_NLOD","m2_FS","m2_ReadPosRankSum","m2_BaseQRankSum","m2_ClippingRankSum",
                               "f_MQM","f_MQMR","f_GTI","f_LEN","f_ODDS",
                               "vs_SSC","vs_SPV","vs_GPV","vs_SS",
                               "vd_SSF","vd_MSI","vd_SOR",
                               "s2_QSS","s2_MQ","s2_SomaticEVS","s2_ReadPosRankSum",
                               "Alt_Allele_Freq",
                               "N_refDepth","N_altDepth","T_refDepth","T_altDepth","relcov")
    rm(meta_indel) 
    #SOR inf values
    parse_indel$vd_SOR[which(is.infinite(parse_indel$vd_SOR)==TRUE)] <-  suppressWarnings(max(parse_indel$vd_SOR[which(is.infinite(parse_indel$vd_SOR)==F)],na.rm=T)+1) #added for v2
    
    #sort table
    parse_indel <- setkey(parse_indel, Chr, START_POS_REF)
    
    # Remove duplicated rows
    parse_indel <- unique(parse_indel)
    
    
    parse_snv = as.data.frame(parse_snv) 
    parse_indel = as.data.frame(parse_indel)
    #print(parse_snv)

    #print('this is indel')
    #print(parse_indel)
    return(list(snv=parse_snv, indel=parse_indel))

}
