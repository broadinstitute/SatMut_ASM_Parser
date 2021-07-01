library(GGally)
library(tidyverse)
library(stringdist)
library(cowplot)
library(ggrepel)

############## ANALYSIS SETUP ###########


# This amplicon should be the PCR products inclusive of PCR primer sequences. The boundary of ORF are marked by '[' before the start codon and ']' after the stop codon. This has to be the EXACTLY the 'reference file' used by AnalyzeSaturationMutagenesis to produce the .codonCounts and .variantCounts.
# This amplicon should be the PCR products inclusive of PCR primer sequences. The boundary of ORF are marked by '[' before the start codon and ']' after the stop codon. This has to be the EXACTLY the 'reference file' used by AnalyzeSaturationMutagenesis to produce the .codonCounts and .variantCounts.

#pMT025-TP53_XY304_305_Amplicon

ORF_Amplicon="ATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCATTCTCAAGCCTCAGACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAGGTGTCGTGAGGCTAGCGCCACC[ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCGCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTAG]ACGCGTTAAGTCGACAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATG"

##dir for .variantCounts
#Make your .variantCounts file names in this format, and place them in a folder for the parser to access.
#1_CJNC2.1.Sample01_KRASG12C.variantCounts
#where 1_CJNC2.1. is lane.flowcell.lane.
#where Sample01 is generic sample number.

dir_variantCounts="~/Desktop/SatMut_ASM_Parser-1.9/Training_data/variantCounts/"

dir_cdnCounts="~/Desktop/SatMut_ASM_Parser-1.9/Training_data/codonCounts/"

sampleAnnot = read_csv("~/Desktop/SatMut_ASM_Parser-1.9/Training_data/P53_Lane13_A4C_SampleAnnot.csv")
#Sample mapping file: two columns named 'Sample' and 'Experiment'
#Example
# Sample	Experiment
# Sample01	pMT025_SMARCB1sat_plasmidPool_PCR
# Sample02	pMT025_SMARCB1sat_plasmidPool_straight

codonDesigned <- read_csv("~/Desktop/SatMut_ASM_Parser-1.9/Training_data/P53_pMT025_codonList.csv")
#planned codon changes: one column named 'key'
#Example
# key
# 1|AAA
# 1|AAT
# 1|ACT
#formatted as aa_position|planned_codon. You can comma-delimit multiple codon changes when a variant was designed to have such changes.


screenNM<-"P53_pMT025_scrn_train"

gene <- 'TP53' # It is part of input file name, so has to be exact.

clonalSample<-NULL #"Sample13" #specify clonal sample number if there is one. Otherwise set clonalSample<-NULL

pDNASample<-c('Sample14') ##specify pDNA library sample number - you should alway carry one. If not, use an ETP sample

refSamples<-c('Sample13','Sample14','Sample13') ## # It needs 3 elements. The reference samples are those that were not selected. e.g. early time point (ETP) samples. If you don't have 3 reference samples, repeat your reference sample name(s) to make a 3-element array.


lowCountCutForRef=2 # counts equal or below this will be filtered out.  0 allows all species

lowCountCutForTreatment=2 # counts equal or below this will be filtered out.  0 allows all species



######################### end of analysis setup

codons.per.pos=22 #if If the 'codonDesigned' is nor populated, this will use to calculate top abund variants as the 'intended'
pos.off.set<-0 ##Set it 0 for mutagenesis starts at codon #1, otherwise set accordingly

samples.for.rank=c(pDNASample,refSamples) # pick 3 samples of pDNA, ETP. These samples may be used to trim the data by removing millions of low-count variants. If you don't have 3 samples, you may repeat one trice.

contrl.pos=NULL

cat("Done with analysis_setup")

############## EEND OF ANALYSIS SETUP ###########

############## FUNCTIONS #################
codonTable<-tibble(
  "CODON"=c('GCT','GCC','GCA','GCG','TGT','TGC','GAT','GAC','GAA','GAG','TTT','TTC','GGT','GGC','GGA','GGG','CAT','CAC','ATT','ATC','ATA','AAA','AAG','TTA','TTG','CTT','CTC','CTA','CTG','ATG','AAT','AAC','CCT','CCC','CCA','CCG','CAA','CAG','CGT','CGC','CGA','CGG','AGA','AGG','TCT','TCC','TCA','TCG','AGT','AGC','TAA','TAG','TGA','ACT','ACC','ACA','ACG','GTT','GTC','GTA','GTG','TGG','TAT','TAC'),
  "AA"=c('A','A','A','A','C','C','D','D','E','E','F','F','G','G','G','G','H','H','I','I','I','K','K','L','L','L','L','L','L','M','N','N','P','P','P','P','Q','Q','R','R','R','R','R','R','S','S','S','S','S','S','X','X','X','T','T','T','T','V','V','V','V','W','Y','Y'))

FromVariantCountToCodon_return9cols_callWithMatrix2<-function(input_str, wt_mtx=wt_lt, off.set){
  ###takes in 806:G>A, 807:A>G, 808:T>C and produces GAT.227.AGC/10
  str_in<-as.character(str_split(input_str,"/")[[1]][1])
  iteration<-as.numeric(str_split(input_str,"/")[[1]][2])
  if(iteration %% (floor(N/1000)*100) == 0){
    cat(paste("Done: ", round((iteration/N) * 100,0),'%\n'))
  } 
  
  #cat(paste0("FromVariant_return9cols: ", str_in,"\n"))
  #set wt nt matrix with 5 additional positions to catch the end of ORF
  wt_lt<-as.matrix(wt_mtx)
  orf_size<-dim(wt_lt)[1]-5
  
  #cat("input orfsize=", orf_size, ", mat_dim =", dim(wt_lt), "\n")
  #str_in=c("128:->A, 128:->T, 128:->G")
  num.calls<-str_count(str_in,">")
  calls<-c()
  if(str_detect(str_in,",")){
    calls<-str_split(str_in, ",")
  }else{
    calls<-c(str_in)
  }
  cdn_sum<-list()
  cdn_type<-list()
  cdnPOS<-list()
  which.nth<-list()
  positions<-c()
  POS_list_to=list()
  POS_list_from=list()
  all<-c()
  key<-c()
  ps<-c()
  froms<-c()
  tos<-c()
  p_cv<-c()
  for(n in 1:num.calls){
    call_single<-NULL
    ntPOS<-NULL
    cdnP<-NULL
    call_single<-calls[[1]][n]
    ntPOS=as.numeric(substr(call_single,1,str_locate(call_single,":")[,1]-1))-off.set
    cdnP=as.numeric(((ntPOS-0.1) %/% 3) +1)
    positions<-c(positions,cdnP)
    POS_list_to[[cdnP]]=c('*','*','*')
    POS_list_from[[cdnP]]=c('*','*','*')
  }
  uniq_POS=unique(positions)
  
  for(n in 1:num.calls){
    cdnPOS<-NULL
    call_single<-NULL
    ntPOS<-NULL
    change_from<-NULL
    change_to<-NULL
    nth_nt_in_cdn<-NULL
    
    call_single<-calls[[1]][n]
    # groupType=ifelse(str_detect(call_single,">-"), "del", 
    #                  ifelse(str_detect(str_in,"->"),"ins","sub"))
    ntPOS=as.numeric(substr(call_single,1,str_locate(call_single,":")[,1]-1))-off.set
    
    change_from<-substr(call_single,nchar(call_single)-2,nchar(call_single)-2)
    change_to<-substr(call_single,nchar(call_single),nchar(call_single))
    nth_nt_in_cdn=ifelse(ntPOS %% 3 ==0, 3,ntPOS%% 3)
    cdnPOS=as.numeric(((ntPOS-0.1) %/% 3) +1)
    if(cdnPOS>dim(wt_lt)[1]){next}
    
    if(change_from != "-" & change_from != wt_lt[cdnPOS,nth_nt_in_cdn]){
      cat("Change_from is not reference base:", call_single, "reference=",wt_lt[cdnPOS,nth_nt_in_cdn],"\n")
      stopifnot(wt_lt[cdnPOS,nth_nt_in_cdn] == 0)
      #       if(wt_lt[cdnPOS,nth_nt_in_cdn] !=0){
      #   exit 
      # }
    }
    
    if(POS_list_to[[cdnPOS]][nth_nt_in_cdn]!="*"){
      POS_list_to[[cdnPOS]][nth_nt_in_cdn]=paste0(POS_list_to[[cdnPOS]][nth_nt_in_cdn],change_to)
    } else{
      POS_list_to[[cdnPOS]][nth_nt_in_cdn]=change_to
    }
    if(POS_list_from[[cdnPOS]][nth_nt_in_cdn]!="*"){
      POS_list_from[[cdnPOS]][nth_nt_in_cdn]=paste0(POS_list_from[[cdnPOS]][nth_nt_in_cdn],change_from)
    } else{
      POS_list_from[[cdnPOS]][nth_nt_in_cdn]=change_from
    }
  }
  for(p in uniq_POS){
    if(p>orf_size || p<0){
      if(is.null(key)){
        all=paste0("OUTBOUND|",p,"|OUTBOUND")
        key=paste0(p,"|OUTBOUND")
        ps=p
        froms="OUTBOUND"
        tos="OUTBOUND"
      }else{
        
        all=paste0(all,",",paste0("OUTBOUND|",p,"|OUTBOUND"))
        key=paste0(key,",",paste0(p,"|OUTBOUND"))
        ps=paste0(ps,",",p)
        froms=paste0(froms,",","OUTBOUND")
        tos=paste0(tos,",","OUTBOUND")
      }
      next
    }
    
    if(POS_list_from[[p]][1]=="*" ||
       is.na(POS_list_from[[p]][1]) ||
       is.null(POS_list_from[[p]][1])){
      POS_list_from[[p]][1]=wt_lt[p,1]
    }
    if(POS_list_from[[p]][2]=="*"|| 
       is.na(POS_list_from[[p]][2])|| 
       is.null(POS_list_from[[p]][2])){
      POS_list_from[[p]][2]=wt_lt[p,2]
    }
    if(POS_list_from[[p]][3]=="*"|| 
       is.na(POS_list_from[[p]][3])|| 
       is.null(POS_list_from[[p]][3])){
      POS_list_from[[p]][3]=wt_lt[p,3]
    }
    
    
    if(POS_list_from[[p]][1] == strrep("-",nchar(POS_list_from[[p]][1]))){
      POS_list_from[[p]][1]=paste0(POS_list_from[[p]][1],wt_lt[p,1])
      POS_list_to[[p]][1]=paste0(POS_list_to[[p]][1],wt_lt[p,1])
    }
    
    
    if(POS_list_from[[p]][2] == strrep("-",nchar(POS_list_from[[p]][2]))){
      POS_list_from[[p]][2]=paste0(POS_list_from[[p]][2],wt_lt[p,2])
      POS_list_to[[p]][2]=paste0(POS_list_to[[p]][2],wt_lt[p,2])
    }
    
    # if(!is.na(str_locate(POS_list_from[[p]][2],"-")[1])){
    #   POS_list_from[[p]][2]=paste0(POS_list_from[[p]][2],wt_lt[p,2])
    #   POS_list_to[[p]][2]=paste0(POS_list_to[[p]][2],wt_lt[p,2])
    # }
    
    if(POS_list_from[[p]][3] == strrep("-",nchar(POS_list_from[[p]][3]))){
      POS_list_from[[p]][3]=paste0(POS_list_from[[p]][3],wt_lt[p,3])
      POS_list_to[[p]][3]=paste0(POS_list_to[[p]][3],wt_lt[p,3])
    }
    
    if(POS_list_to[[p]][1]=="*" ||
       is.na(POS_list_to[[p]][1])||
       is.null(POS_list_to[[p]][1])){
      POS_list_to[[p]][1]=wt_lt[p,1]
    }
    if(POS_list_to[[p]][2]=="*" ||
       is.na(POS_list_to[[p]][2])||
       is.null(POS_list_to[[p]][2])){
      POS_list_to[[p]][2]=wt_lt[p,2]
    }
    if(POS_list_to[[p]][3]=="*" ||
       is.na(POS_list_to[[p]][3])||
       is.null(POS_list_to[[p]][3])){
      POS_list_to[[p]][3]=wt_lt[p,3]
    }
    if(is.null(key)){
      key=paste0(p,"|",paste(POS_list_to[[p]][1],POS_list_to[[p]][2],POS_list_to[[p]][3],sep=""))
      all=paste0(
        paste(POS_list_from[[p]][1],POS_list_from[[p]][2],POS_list_from[[p]][3],sep=""),"|",p,"|",paste(POS_list_to[[p]][1],POS_list_to[[p]][2],POS_list_to[[p]][3],sep=""))
      ps=p
      froms<-paste0(POS_list_from[[p]][1],POS_list_from[[p]][2],POS_list_from[[p]][3])
      tos<-paste0(POS_list_to[[p]][1],POS_list_to[[p]][2],POS_list_to[[p]][3])
    }else{
      key=paste0(key,",",paste0(
        p,"|",paste(POS_list_to[[p]][1],POS_list_to[[p]][2],POS_list_to[[p]][3],sep="")))
      all=paste0(all,",",paste0(
        paste(POS_list_from[[p]][1],POS_list_from[[p]][2],POS_list_from[[p]][3],sep=""),"|",p,"|",paste(POS_list_to[[p]][1],POS_list_to[[p]][2],POS_list_to[[p]][3],sep="")))
      ps<-paste(ps,p,sep=",")
      froms<-paste(froms,paste0(POS_list_from[[p]][1],POS_list_from[[p]][2],POS_list_from[[p]][3]),sep=",")
      tos<-paste(tos,paste0(POS_list_to[[p]][1],POS_list_to[[p]][2],POS_list_to[[p]][3]),sep=",")
    }
  }
  nt_ins=str_count(froms,"-")
  nt_del=str_count(tos,"-")
  ps_mn=mean(uniq_POS)
  fs<-0
  if((nt_ins+nt_del) %% 3 !=0){ fs<- 1}
  
  return(c(str_in, all, key, ps, froms, tos, nt_ins, nt_del, ps_mn, fs))
  
}
# 
# FromVariantCountToCodon_return9cols_callWithMatrix2("input_str"="2419:->G, 2419:->T, 2419:->T, 2419:A>T, 2432:G>C/10", "wt_mtx"=wt_lt, "off.set"=off.set)
# 
# 
# FromVariantCountToCodon_return9cols_callWithMatrix2("input_str"="2419:->G, 2419:->T, 2432:G>C/10", "wt_mtx"=wt_lt, "off.set"=off.set)

FromRawToPercent_cbindRawCounts <- function(rawFile, codonTb, colStart, colEnd) {
  ###take rawfile(POS/CODON/COUNT) to get fractionfiles, w or wo wt, cbind to the raw columns
  codonTb<-codonTb
  wf0<-rawFile
  colStart<-colStart
  colEnd<-colEnd
  
  wf<-merge(codonTb, wf0, by.x='CODON', by.y='CODON',all=TRUE)
  wf$AA3<-NULL
  wf$Aalong<-NULL
  wf <- wf[order(wf$POS, wf$CODON),]
  wf$X<-NULL
  
  print(str(wf))
  samples<-names(wf[-c(1,2,3)])
  samples_raw<-str_c(samples,'_raw',sep='')
  samples_frctn<-str_c(samples,'_frctn',sep='')
  
  ###########
  #iris %>% mutate_at(.cols=vars(-Species), .funs=funs(mysum = sum(.))) %>% head()
  fraction<-function (x){x/sum(x)}
  
  wf0<-wf
  names(wf0)[-c(1,2,3)]<-samples_raw
  wf1<-wf
  names(wf1)[-c(1,2,3)]<-samples_frctn
  wf2<-cbind(wf0,wf1[,-c(1,2,3)])
  
  wf2$sampleForRank<-wf2[[samples_raw[1]]]
  wf_inPerct_w_wt <-
    wf2 %>%
    dplyr::group_by (POS) %>%
    dplyr::mutate(Rank=rank(-sampleForRank)) %>%
    dplyr::mutate_at(samples_frctn, funs(fraction)) 
  
  wf00<-wf2
  wf00$sampleForRank<-wf00[[samples_raw[1]]]
  
  sub_wf0<-
    wf00 %>%
    dplyr::group_by (POS) %>%
    dplyr::mutate(Rank=rank(-sampleForRank)) %>%
    dplyr::filter(Rank ==1) %>%
    dplyr::select(CODON, AA, POS)
  
  names(sub_wf0)[1:3]<-c('Wt_codon','Wt_aa', 'POS')
  names(wf_inPerct_w_wt)[1:3]<-c('Vt_codon', 'Vt_aa', 'POS')
  
  wf_inPerct_w_wt<-merge(sub_wf0, wf_inPerct_w_wt, by='POS', all=TRUE)
  wf_inPerct_wo_wt<-wf_inPerct_w_wt[wf_inPerct_w_wt$Rank!=1,]
  wf_inPerct_w_wt$sampleForRank<-NULL 
  wf_inPerct_w_wt$X<-NULL
  wf_inPerct_w_wt$Rank<-NULL
  wf_inPerct_wo_wt$sampleForRank<-NULL 
  wf_inPerct_wo_wt$X<-NULL
  wf_inPerct_wo_wt$Rank<-NULL
  
  # wf_inPerct<-wf_inPerct[-1] ##remove X column
  return(list("wf_fraction_wo_wt"=wf_inPerct_wo_wt, "wf_fraction_w_wt"=wf_inPerct_w_wt))
}
#####

GetPrimaryKey_theKey_theDeltaNt <- function(keys_itr,keyList=IntendedKey,wtCdns=wtCdns){
  ###from key it returns primarykeys(all intended pasted), theKey(best deltaNT one), and the deltaNT of theKey
  keys_itr=as.character(keys_itr)
  str_in<-as.character(str_split(keys_itr, "/")[[1]][1])
  iteration<-as.numeric(str_split(keys_itr, "/")[[1]][2])
  
  if(iteration %% (floor(M/1000)*100) == 0){
    cat(paste("Done: ", round((iteration/M) * 100,0),'%\n'))
  } 
  #cat(paste0("GetPrimaryKey_theKey_theDeltaNt: ", str_in,"\n"))
  keyList=as.vector(keyList)
  wtCdns=as.vector(wtCdns)
  num.calls<-str_count(str_in,",")+1
  calls<-c()
  primaryKeys<-c()
  delta<-NA
  theKey<-NA
  primaryKey<-NA
  outstr<-NULL
  if(str_in %in% keyList){
    ##cat(str_in)
    delta<-0
    primaryKey<-str_in
    theKey<-str_in
    if(str_detect(str_in,",")){  #for doubles
      calls<-as.vector(str_split(str_in, ","))
      for(i in 1:num.calls){
        k=calls[[1]][i]
        pos<-NULL
        pos<-as.numeric(str_split(k, "\\|")[[1]][1])
        vtCdn<-as.character(str_split(k, "\\|")[[1]][2])
        delta<-delta+stringdist(wtCdns[pos],vtCdn,method = "hamming")
      }
      #      return(list(primarykey,theKey,delta))
      
    }else{
      pos<-as.numeric(str_split(str_in, "\\|")[[1]][1])
      vtCdn<-as.character(str_split(str_in, "\\|")[[1]][2])
      delta<-stringdist(wtCdns[pos],vtCdn,method = "hamming")
      #      return(list(primarykey,theKey,delta))
    }
    #return(list(primaryKey=primaryKey,theKey=theKey,maxDelta=delta))
    outstr<-paste(primaryKey,theKey,delta,sep=";")
    return(outstr)
  }
  
  calls<-c()
  primaryKeys<-c()
  delta<-NA
  theKey<-NA
  primaryKey<-NA
  outstr<-NULL
  
  if(str_detect(str_in,",")){
    calls<-as.vector(str_split(str_in, ","))
  }else{
    calls<-as.vector(c(str_in))
  }
  
  deltas<-c()
  delta<-0
  ##if multiple intended codon present in same mol, use the highest deltant as the anchor.
  for(i in 1:num.calls){
    if(calls[[1]][i] %in% keyList){
      ky=calls[[1]][i]
      pos<-NULL
      pos<-as.numeric(str_split(ky, "\\|")[[1]][1])
      vtCdn<-as.character(str_split(ky, "\\|")[[1]][2])
      deltas<-c(deltas,stringdist(wtCdns[pos],vtCdn,method = "hamming"))
      primaryKeys<-c(primaryKeys,ky)
      if(is.na(primaryKey)){
        primaryKey<-ky
      }else{
        primaryKey<-paste(primaryKey,ky,sep=',')
      }
    }else{
      next
    }
  }
  
  #cat("primaryKey=")
  #cat(primaryKey)
  # cat("\n")
  
  if(is.na(primaryKey)==TRUE){
    outstr<-"NA;NA;NA"
    return(outstr)
    #return(list(primaryKey=NA,theKey=NA,maxDelta=NA))
  }else{
    maxDelta=max(deltas,na.rm=TRUE)
    theKey<-primaryKeys[which.max(deltas)]
    outstr<-paste(primaryKey,theKey,maxDelta,sep=";")
    return(outstr)
    #return(list(primaryKey=primaryKey,theKey=theKey,maxDelta=max(deltas,na.rm=TRUE)))
  }
}



func_hist_from_raw <- function (fileIn, column, col='black', txhi=0, breaks=50){
  
  workFile1 <- fileIn;  
  column <- column;
  col1 <- col
  breaks <- breaks
  hight <-txhi
  library(ggplot2)
  
  workFile1$G0SumNorm<-log2(workFile1[[column]]/sum(workFile1[[column]])*dim(workFile1)[1]*1000+1)
  workFile2<-workFile1[workFile1$G0SumNorm>0,] ###count out 0 for avg...
  avg<-round(mean(workFile2[[column]]),1)
  
  g<-ggplot(workFile1, aes(G0SumNorm))+xlab(paste('log2Norm','(',column, '); avg ', avg,'reads', by=''))
  g<-g+theme(axis.text=element_text(size=8),
             axis.title=element_text(size=10,face="bold"))
  g<-g+geom_histogram(breaks=breaks, colour ="black",fill=col1) 
  g<-g+geom_vline(xintercept=mean(workFile2$G0SumNorm)+
                    sd(workFile2$G0SumNorm),
                  col='red', lty=6, lwd=0.5)
  
  g<-g+geom_text(label=round(mean(workFile2$G0SumNorm)+
                               sd(workFile2$G0SumNorm),1), 
                 x=mean(workFile2$G0SumNorm)+sd(workFile2$G0SumNorm), 
                 y=hight+100,  size=5,colour='red',angle=90)
  ###
  g<-g+geom_vline(xintercept=mean(workFile2$G0SumNorm)-
                    sd(workFile2$G0SumNorm),
                  col='red', lty=6, lwd=0.5)
  
  g<-g+geom_text(label=round(mean(workFile2$G0SumNorm)-
                               sd(workFile2$G0SumNorm),1), 
                 x=mean(workFile2$G0SumNorm)-sd(workFile2$G0SumNorm), 
                 y=hight+100,size=5,colour='red',angle=90)
  ####
  g<-g+geom_vline(xintercept=mean(workFile2$G0SumNorm)-
                    sd(workFile2$G0SumNorm*2),
                  col='blue', lty=6, lwd=0.5)
  
  g<-g+geom_text(label=round(mean(workFile2$G0SumNorm)-
                               sd(workFile2$G0SumNorm)*2,1), 
                 x=mean(workFile2$G0SumNorm)-sd(workFile2$G0SumNorm)*2, 
                 y=hight+100, size=5,colour='red',angle=90)
  ###
  g<-g+geom_vline(xintercept=mean(workFile2$G0SumNorm)+
                    sd(workFile2$G0SumNorm*2),
                  col='blue', lty=6, lwd=0.5)
  
  g<-g+geom_text(label=round(mean(workFile2$G0SumNorm)+
                               sd(workFile2$G0SumNorm)*2,1), 
                 x=mean(workFile2$G0SumNorm)+sd(workFile2$G0SumNorm)*2, 
                 y=hight+100,  size=5,colour='red',angle=90)
  ###
  
  g<-g+geom_vline(xintercept=mean(workFile2$G0SumNorm),
                  col='green', lty=6, lwd=0.5)
  g<-g+geom_text(label=round(mean(workFile2$G0SumNorm),1), 
                 x=mean(workFile2$G0SumNorm), 
                 y=hight+100, size=5,colour='red',angle=90)
  return (list("avgReads"=avg, "hist"=g));
}



fracFunNoBlueNoLabel_colArg_generic_auc <- function (fileIn, columns, filterColumn='', subset='', col='black', cexin=1, xlim=c(0,1), ylim=c(0,1),alpha=1){
  inputFile <- fileIn;
  columnList <- columns;
  filterCol <- filterColumn;
  filter <- subset;
  col1 <- col;
  pch1 <- cexin;
  xlm <-xlim;
  ylm <-ylim;
  i <-alpha;
  ##before looking at summation, check if replicates are similar
  list <- paste(columnList, collapse="/");
  if(length(columnList)==2){
    plot(inputFile[[columnList[1]]] ~ inputFile[[columnList[2]]], main = paste("Replicates reproducible?\n(", list, ")"));
  }
  if(filterCol != "" && filter != ""){
    subset_inputFile <- inputFile[inputFile[[filterCol]] == filter,];
  }
  else {
    subset_inputFile <- inputFile;
  }
  subset_inputFile$tmp_sum<-subset_inputFile[[columnList[1]]];
  #str(subset_inputFile);
  if(length(columnList)>1){
    for (i in 2:length(columnList)){
      col <- columnList[i];
      subset_inputFile$tmp_sum <- subset_inputFile$tmp_sum + subset_inputFile[[col]];
    }
  }
  #str(subset_inputFile);
  sorted_subset_inputFile <- subset_inputFile[order(-subset_inputFile$tmp_sum),]
  sorted_subset_inputFile$readFraction <- sorted_subset_inputFile$tmp_sum/sum(sorted_subset_inputFile$tmp_sum)
  sorted_subset_inputFile$cumsum <- cumsum(sorted_subset_inputFile$tmp_sum)
  sorted_subset_inputFile$cumsumFraction <- sorted_subset_inputFile$cumsum/sum(sorted_subset_inputFile$tmp_sum)
  sorted_subset_inputFile$cloneFraction <- rank(-sorted_subset_inputFile$tmp_sum, ties.method='first')/length(sorted_subset_inputFile$tmp_sum)
  
  y2max <- max(sorted_subset_inputFile$cumsumFraction) 
  require(MESS)
  library(ggplot2)
  require(ggplot2)
  
  a<-round(auc(sorted_subset_inputFile$cloneFraction,sorted_subset_inputFile$cumsumFraction,
               from=xlm[1], to=xlm[2],type = 'spline'),4)
  sampleNM<-paste(columnList, 'auc=', a, sep=' ')
  g=ggplot(data=data.frame(sorted_subset_inputFile), aes(x=cloneFraction,y=cumsumFraction))+
    xlab(paste('CloneFraction','(',sampleNM, ')', by=''))+
    geom_point(colour=col1, pch=1, cex=ifelse(sorted_subset_inputFile$cumsumFraction==1, 5*pch1, pch1))
  g<-g+theme(axis.text=element_text(size=8),
             axis.title=element_text(size=10,face="bold"))
  #g=g+geom_text(label=paste(columnList, 'auc=', a, sep=' '),cex=3, x=0.5, y=0.2);
  return (list("auc"=a,"gph"=g))
}

renameSample2Experiment<-function (dataFile, annotFile){
  df<-data.frame(dataFile)
  annotF<-data.frame(annotFile)
  nms=names(df)
  for(i in 1:dim(df)[2]){
    if(str_sub(nms[i],1,6)=='Sample'){
      nms=gsub(str_sub(nms[i],1,8),annotF[sampleAnnot$Sample==str_sub(nms[i],1,8),]$Experiment, nms)
    }
  }
  names(df)<-nms
  return(df)
}



############## END OF FUNCTIONS #################


















if(!dir.exists(paste0(dir_variantCounts,"outbox_",screenNM,"/"))){
  dir.create(file.path(paste0(dir_variantCounts,"outbox_",screenNM,"/")))
}
dir_wrt_variantCounts=paste0(dir_variantCounts,"outbox_",screenNM,"/")

if(!dir.exists(paste0(dir_cdnCounts,"outbox_",screenNM,"/"))){
  dir.create(file.path(paste0(dir_cdnCounts,"outbox_",screenNM,"/")))
}
dir_wrt_cdnCounts=paste0(dir_cdnCounts,"outbox_",screenNM,"/")


variantCountsfiles =
  list.files(path = dir_variantCounts, pattern = '.variantCounts$', all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


outBox_VrtCt <- paste0(dir_variantCounts,"outBox_",lowCountCutForRef,"refSampleCtCutOff_",Sys.Date(),"/")
if(!dir.exists(path=outBox_VrtCt)){
  print ("new folder need to make")
  dir.create(path=outBox_VrtCt,mode="0777", showWarnings = TRUE)
}


vtCtTblfolder <- paste0(outBox_VrtCt,"variantCountsTables_",
                  lowCountCutForRef,"RefCtCutOff_", lowCountCutForTreatment,"others/")

if(!dir.exists(path=vtCtTblfolder)){
  print ("new folder need to make")
  dir.create(path=vtCtTblfolder,mode="0777", showWarnings = TRUE)
}

#for(i in 1:3){
for(i in 1:length(variantCountsfiles)){
  if(!str_detect(variantCountsfiles[i],'.variantCounts$')){next}
  sampl=str_sub(variantCountsfiles[i],str_locate(variantCountsfiles[i],'Sample')[1],str_locate(variantCountsfiles[i],'Sample')[1]+7)
  tmp_df<-NULL
  path0=paste0(dir_variantCounts,variantCountsfiles[i])
  tmp_df<-read.table(file=path0,fill=TRUE,sep='\t', header=FALSE,col.names = paste0("V",seq_len(9)))
  tmp_df$X<-NULL
  tmp_df$X1<-NULL
  
  names(tmp_df)
  tmp_df<-tmp_df%>%
    mutate(WtLen=V1*V3) %>%
    select(V1,V2,WtLen,V4,V5,V6,V7,V8,V9)
  if((sampl %in% refSamples) | (sampl %in% pDNASample)){
    tmp_df<-tmp_df%>%
      filter(as.numeric(V1)>lowCountCutForRef)
  }else{
    tmp_df<-tmp_df%>%
      filter(as.numeric(V1)>lowCountCutForTreatment)
  }
  
  print(paste0("Done reading in :",variantCountsfiles[i]))
  print(dim(tmp_df))
  
  write.csv(tmp_df,paste0(vtCtTblfolder,variantCountsfiles[i],".tbl"))
}

###Done parsing .varriantCoounts into a 9-cols. table
## ready to consolidate down to samples.

variantCountsTblfiles =
  list.files(path = vtCtTblfolder, pattern = NULL, all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


vtCtTblConslFolder <- paste0(outBox_VrtCt,"variantCountsTables_consolidated_",lowCountCutForRef,"RefCtCutOff_", lowCountCutForTreatment,"others/")

if(!dir.exists(path=vtCtTblConslFolder)){
  print ("new folder need to make")
  dir.create(path=vtCtTblConslFolder,mode="0777", showWarnings = TRUE)
}

samplesInData<-unique(str_sub(variantCountsTblfiles,
                              str_locate(variantCountsTblfiles,'Sample')[1],str_locate(variantCountsTblfiles,'Sample')[1]+7))


for(j in (1:length(samplesInData))){
  df_consl_tmp<-data.frame()
  for(k in (1:length(variantCountsTblfiles))){
    if(!str_detect(variantCountsTblfiles[k],samplesInData[j])){next}
    df_in_tmp<-read_csv(paste0(vtCtTblfolder,variantCountsTblfiles[k]))%>%
      select(V1,V2,WtLen,V4,V5,V6,V7,V8,V9)
    summary(df_in_tmp)
    print(paste0("Done with: ",variantCountsTblfiles[k]))
    print(dim(df_in_tmp))
    
    if(dim(df_consl_tmp)[1]<1){
      df_consl_tmp=df_in_tmp
    }else{
      df_consl_tmp=rbind(df_consl_tmp,df_in_tmp)
    }
  }
  
  df_consl_tmp_aggd<-df_consl_tmp%>%
    select(V1,V2,WtLen,V4,V5,V6,V7,V8,V9)%>%
    group_by(V4,V5,V6,V7,V8,V9)%>%
    summarise(V1=sum(V1,na.rm = TRUE),
              V2=sum(V2,na.rm = TRUE),
              WtLen=sum(WtLen,na.rm = TRUE))
  
  write.csv(df_consl_tmp_aggd,paste0(vtCtTblConslFolder,samplesInData[j],".consolidatedTbl"))
  print(paste0("Done with: ",samplesInData[j]))
  print(dim(df_consl_tmp_aggd))
}

####Done with consolidating multiiple lanes of data.  
### Ready to do merging into a single .csv

consol_sample_files =
  list.files(path = vtCtTblConslFolder, pattern = NULL, all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

# outBox_VrtCt <- paste0(dir_variantCounts,"outBox_",lowCountCutForRef,"refSampleCtCutOff_",Sys.Date(),"/")
# if(!dir.exists(path=outBox_VrtCt)){
#   print ("new folder need to make")
#   dir.create(path=outBox_VrtCt,mode="0777", showWarnings = TRUE)
# }


df_in1<-data.frame()
for(i in(1:length(consol_sample_files))){
  tmp_df<-NULL
  smp<-str_sub(consol_sample_files[i],1,8)
  tmp_df<-read.csv(paste0(vtCtTblConslFolder,consol_sample_files[i] ))%>%
    select(V1,V2,WtLen,V4,V5,V6,V7,V8,V9)
  
  #remove incoming data rows that have NA in any cell.
  row.has.na.tmp <- apply(tmp_df[,1:6], 1, function(x){any(is.na(x))})
  tmp_df<-tmp_df[!row.has.na.tmp,]
  tmp_df<-as_tibble(tmp_df)
  tmp_df$V1<-as.numeric(as.character(tmp_df$V1))
  tmp_df$V2<-as.numeric(as.character(tmp_df$V2))
  tmp_nonNA_dim<-dim(tmp_df)[1]
  
  names(tmp_df)<-paste0(smp,".",c("ct", "ref_ct", "len","changes","variant","cdn.changes.by.Ted", "Vt.Cdn.By.Ted","Vt.aa.By.Ted", "StdNom.Vt.aa.By.Ted"))
  names(tmp_df)[c(4:9)]<-c("changes","variant","cdn.changes.by.Ted","Vt.Cdn.By.Ted","Vt.aa.By.Ted","StdNom.Vt.aa.By.Ted")
  if(dim(df_in1)[1]<1){
    df_in1<-tmp_df
  }else{
    df_in1<-merge(df_in1,tmp_df, by.x=c("changes","variant","cdn.changes.by.Ted","Vt.Cdn.By.Ted","Vt.aa.By.Ted","StdNom.Vt.aa.By.Ted"), by.y=c("changes","variant","cdn.changes.by.Ted","Vt.Cdn.By.Ted","Vt.aa.By.Ted","StdNom.Vt.aa.By.Ted"), all=TRUE)
  }
  cat(paste0("done with:",smp," nrow=", dim(df_in1)[1],"\n"))
}


df_in<-df_in1 %>%
  select(variant,changes,everything())

df_in$X<-NULL
names(df_in)


write.csv(df_in, file=paste0(dir_wrt_variantCounts,"All_variantCounts_aggd.csv"))

##########33##########################





#similarly to .variantCounts
##dir for .codonCounts
#make your .codonCounts file names in this format, and place them in a folder for the parser to access.
#1_CJNC2.1.Sample01_KRASG12C.codonCounts
#where 1_CJNC2.1. is lane.flowcell.lane.
#where Sample01 is generic sample number.



codonCountsfiles =
  list.files(path = dir_cdnCounts, pattern =".codonCounts$", all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


samplesInCdnData<-unique(str_sub(codonCountsfiles,
                                 str_locate(codonCountsfiles,'Sample')[1],str_locate(codonCountsfiles,'Sample')[1]+7))

cdnCount_outBox <- paste0(dir_cdnCounts,"codonCounts_outBox/")
if(!dir.exists(path=cdnCount_outBox)){
  print ("new folder need to make")
  dir.create(path=cdnCount_outBox,mode="0777", showWarnings = TRUE)
}

for(j in 1:length(samplesInCdnData)){
  tmp_tbl_aggd<-data.frame()
  for (k in 1:length(codonCountsfiles)){
    if(!str_detect(codonCountsfiles[k],'.codonCounts$')){next}
    if(!str_detect(codonCountsfiles[k],samplesInCdnData[j])){next}
    tmp_tbl=read.delim(paste0(dir_cdnCounts,codonCountsfiles[k]))[1:64]
    if(dim(tmp_tbl_aggd)[1]<1){
      tmp_tbl_aggd=tmp_tbl
    }else{
      tmp_tbl_aggd=tmp_tbl_aggd+tmp_tbl
    }
  }
  write_csv(tmp_tbl_aggd,paste0(cdnCount_outBox,samplesInCdnData[j],".cdnConsolidated"))
}

codonCountsConsolfiles =
  list.files(path = cdnCount_outBox, pattern =".cdnConsolidated$", all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

cdn_df<-data.frame()
for(i in (1:length(codonCountsConsolfiles))){
  smp=str_sub(codonCountsConsolfiles[i],1,8)
  cdntf_in=read_csv(paste0(cdnCount_outBox,codonCountsConsolfiles[i]))
  cdntf_in$X<-NULL
  cdntf_in$X1<-NULL
  cdntf_in$POS=row.names(cdntf_in)
  tmp.m=melt(cdntf_in,value.name = smp, variable.name = 'AA_codon')
  
  if(dim(cdn_df)[1]<1){
    cdn_df=tmp.m
  } else{
    cdn_df<-merge(cdn_df,tmp.m,by=c('POS','AA_codon'), all=TRUE)
  }
}

##rename the col.names

if(exists("sampleAnnot")){
  NewNamescdn<-c()
  for (x in 1:dim(cdn_df)[2]){
    NewNamescdn<-c(NewNamescdn,
                   ifelse(!is.element(names(cdn_df)[x], sampleAnnot$Sample),names(cdn_df)[x], as.character(sampleAnnot[sampleAnnot$Sample==names(cdn_df)[x],]$Experiment)))
  }
  names(cdn_df)=NewNamescdn
}

write.csv(cdn_df,paste0(cdnCount_outBox,"codonCount_theSingleFile.csv"))




######################## check if orf template is right
wf_0cdn<-read.csv(paste0(cdnCount_outBox,"codonCount_theSingleFile.csv"))
wf_0cdn$X<-NULL

wtCdn_detected<-wf_0cdn%>%
  dplyr::group_by(POS) %>%
  dplyr::mutate(rk=rank(-get(pDNASample))) %>%
  filter(rk==1)%>%
  select(POS,AA_codon) %>%
  arrange(as.numeric(POS))

TemplateORFSeq_bypDNASeq=paste(wtCdn_detected$AA_codon, collapse="")


###function
list.string.diff<-function(a,b,ignore.case=TRUE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case==TRUE)
  {
    a<-toupper(a)
    b<-toupper(b)
  }
  seq.a<-unlist(strsplit(a,split=""))
  seq.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(seq.a,seq.b)
  only.diff<-diff.d[,diff.d[1,]!=diff.d[2,]]
  pos<-which(diff.d[1,]!=diff.d[2,])
  only.diff<-rbind(pos,only.diff)
  return(only.diff)
}
#############

##Cleanup input data
ORF_Amplicon<-gsub("\\t","",gsub("\\s","",gsub("\\n","",toupper(ORF_Amplicon))))
ORF_Amplicon

if(exists("codonDesigned")){
  IntendedKey <- as.vector(as.character(codonDesigned$key))
} else{IntendedKey <-NULL}


#cat (ORF_Amplicon)

left_flnk=strsplit(ORF_Amplicon,"\\[")[[1]][1]
orfseq=strsplit(strsplit(ORF_Amplicon,"\\[")[[1]][2],"\\]")[[1]][1]
right_flnk=strsplit(strsplit(ORF_Amplicon,"\\[")[[1]][2],"\\]")[[1]][2]
off.set=nchar(left_flnk)

orf_size=nchar(orfseq)/3
nts <- as.vector(strsplit(orfseq,""))


wt_list<-matrix(0, nrow=orf_size, ncol=3,
                byrow=TRUE)
wtCdns<-c()

counter=1
for(r in 1:orf_size){
  cdn<-""
  for(c in 1:3){
    wt_list[r,c]=str_sub(orfseq,counter,counter)
    cdn<-paste0(cdn,str_sub(orfseq,counter,counter))
    counter=counter+1
  }
  wtCdns<-c(wtCdns,cdn)
}

###
orfseq=orfseq
orf_size=nchar(orfseq)/3
nts <- as.vector(strsplit(orfseq,""))
wt_lt<-matrix(0, nrow=orf_size+5, ncol=3,
              byrow=TRUE)

counter=1
for(r in 1:orf_size){
  for(c in 1:3){
    wt_lt[r,c]=str_sub(orfseq,counter,counter)
    counter=counter+1
  }
}

dim(wt_lt)
#### wt_lt is now a matrix holding single nt. will be used to call parser
####
#1_HKWYHDSXX.1.Sample01_PPM1D.variantCounts
date_str<-format(Sys.time(), "%Y%m%d")

cat("Done with Chunks prior to read_in_variantCounts\n")



if(TemplateORFSeq_bypDNASeq==orfseq){
  cat("\nTemplateORFSeq_bypDNASeq MATCHES: GOOD NEWS!\n")
} else{  cat("\nTemplateORFSeq_bypDNASeq DO NOT MATCH ****\nTemplateORFSeq_bypDNASeq:\n")
  cat(TemplateORFSeq_bypDNASeq)
  cat("\nORFseq assumption:\n")
  cat(orfseq)
  cat("\n")
  nt_diff<-as.data.frame(list.string.diff("a"=orfseq, "b"=TemplateORFSeq_bypDNASeq))
  cat("List of discrepancies: seq.a=assumed; seq.b=actual\n")
  print(nt_diff)
}
#########################




if(!exists("df_in")){
  df_in<-read.csv(file=paste0(dir_wrt_variantCounts,"All_variantCounts_aggd.csv"))
}
df_in$X<-NULL

#names(df_in)

df_in_wGroupType<-df_in%>%
  mutate(groupType=ifelse(str_detect(variant,">-") & str_detect(variant,"->"), "ins_del", ifelse(str_detect(variant,">-"),"del", ifelse(str_detect(variant,"->"),"ins","sub"))))%>%
  mutate(nt_changes=str_count(variant,",")+1)

sampleNMs<-sampleAnnot[!str_detect(str_to_upper(sampleAnnot$Experiment),'HOP')]$Sample
smp_ls<-paste0(sampleNMs,".ct")
smp_ls_ref<-paste0(sampleNMs,".ref_ct")
smp_ls_frctn<-paste0(sampleNMs,".frctn")

##############
df_in$X<-NULL
#names(df_in)

top.nth=1000000 # down to top million 

row.has.low.refct.all <- apply(df_in[,seq(8,dim(df_in)[2]-1,by=3)], 1, function(x){all(is.na(x) | x<10000)})

sum(row.has.low.refct.all)

#names(row.has.low.refct.all)

samples.for.rank.ct=str_c(samples.for.rank, ".ct", sep = '')
samples.for.rank.ref=str_c(samples.for.rank, ".ref_ct", sep = '')



df_in_ex_low_ct<-df_in[!row.has.low.refct.all,]

nms<-names(df_in_ex_low_ct)
df_in_annotated_frctn<-df_in_ex_low_ct
ct.col.indx<-which(str_detect(names(df_in_annotated_frctn),'\\.ct$'))
for(i in ct.col.indx){
  df_in_annotated_frctn[[paste0(nms[i],"_frctn")]]=
    df_in_annotated_frctn[[nms[i]]]/df_in_annotated_frctn[[nms[i+1]]]
}

###for correction
nms<-names(df_in_annotated_frctn)
df_in_annotated_frctn_corrctd<-df_in_annotated_frctn
fractionNMs<-names(df_in_annotated_frctn_corrctd)[which(str_detect(names(df_in_annotated_frctn_corrctd),'frctn$'))]
correctionNMs<-gsub('frctn$', 'frctn_corrctd', fractionNMs)

if(!is.null(clonalSample)){
  clonalCol<-which(nms==paste0(clonalSample,".ct_frctn"))
  for(i in (1:length(correctionNMs))){
    df_in_annotated_frctn_corrctd[[correctionNMs[i]]]<- 
      ifelse(is.na(df_in_annotated_frctn_corrctd[[names(df_in_annotated_frctn_corrctd)[clonalCol]]]) |
               (df_in_annotated_frctn_corrctd[[names(df_in_annotated_frctn_corrctd)[clonalCol]]] == 'NA'),
             df_in_annotated_frctn_corrctd[[fractionNMs[i]]],    ifelse(as.numeric(df_in_annotated_frctn_corrctd[[fractionNMs[i]]]-df_in_annotated_frctn_corrctd[[names(df_in_annotated_frctn_corrctd)[clonalCol]]])>=0,
                                                                        as.numeric(df_in_annotated_frctn_corrctd[[fractionNMs[i]]]-df_in_annotated_frctn_corrctd[[names(df_in_annotated_frctn_corrctd)[clonalCol]]]),
                                                                        df_in_annotated_frctn_corrctd[[fractionNMs[i]]]))
  }
}else{
  for(i in (1:length(correctionNMs))){
    df_in_annotated_frctn_corrctd[[correctionNMs[i]]]<-df_in_annotated_frctn_corrctd[[fractionNMs[i]]]
  }
}

######Done with correction

#names(df_in_annotated_frctn_corrctd)
dim(df_in_annotated_frctn_corrctd)
pDNASample.ct<-paste0(pDNASample, ".ct")
pDNASample.corrctd<-paste0(pDNASample, ".ct_frctn_corrctd")
clonalSample.frctn<-NULL
if(!is.null(clonalSample)){
  clonalSample.frctn<-paste0(clonalSample, ".ct_frctn")
}

pDNASample.frctn<-paste0(pDNASample, ".ct_frctn")

df_in_annotated_frctn_corrctd<-df_in_annotated_frctn_corrctd%>% 
  mutate(groupType=ifelse(str_detect(variant,">-") & str_detect(variant,"->"), "ins_del", ifelse(str_detect(variant,">-"),"del", ifelse(str_detect(variant,"->"),"ins","sub"))))%>%
  mutate(nt_change=str_count(variant,">"))%>%
  mutate(ntPOS=as.numeric(substr(variant,1,str_locate(variant,":")[,1]-1))-off.set)%>%
  filter(ntPOS>0 & ntPOS<=orf_size*3+5)%>%
  mutate(nth_nt_in_cdn=ifelse(ntPOS %% 3 ==0, 3,ntPOS %% 3))%>%
  mutate(POS=as.numeric(((ntPOS-0.1) %/% 3) +1))%>% 
  mutate(rank.by.pDNA=as.numeric(rank(-get(pDNASample.frctn),na.last = TRUE, ties.method = 'random'))) %>%  #rank by pDNA library sequencing
  mutate(rank.by.corrected.pDNA=as.numeric(rank(-get(pDNASample.corrctd),na.last = TRUE, ties.method = 'random')))%>%
  select(variant,POS,rank.by.pDNA,rank.by.corrected.pDNA,groupType,nt_change,ntPOS,nth_nt_in_cdn,everything()) %>%
  arrange(POS)

#names(df_in_annotated_frctn_corrctd)#678097
table(df_in_annotated_frctn_corrctd$groupType)

# 
#     del     ins ins_del     sub 
#  244963   78216   45845  421587 

cat("Done with: plotting abundance \n\n Next up: 10-column annotation\n\n")

#names(df_in_annotated_frctn_corrctd)
df9<-df_in_annotated_frctn_corrctd

df9_1<-data.frame(cbind("variant"=as.character(df9[,1]),"iteration"=1:dim(df9)[1])) %>%
  mutate(variant_itr=str_c(variant,iteration,sep="/")) %>%
  select(variant_itr)

N=dim(df9_1)[1]
df_variant2details<-data.frame(t(as.data.frame(apply(df9_1, 1, FUN=function(x){FromVariantCountToCodon_return9cols_callWithMatrix2(input_str = x, wt_mtx =wt_lt, off.set = off.set)}))))

names(df_variant2details)<-c("variant","variant_newName","keys", "cdn.positions", "change.from", "change.to", "nt_ins", "nt_del", "lesion.center", "frame.shift")

dim(df9)
dim(df_variant2details)
############
write.csv(df_variant2details, file=paste0(dir_wrt_variantCounts,"df_variant2details.csv"))

df_variant2details$variant<-as.character(df_variant2details$variant)
df9$variant<-as.character(df9$variant)

df10<-merge(df_variant2details,df9,by="variant",all=TRUE)

cat("Done with: filterig data and adding 9-col annotations and df_variant2details.csv is saved\n")

# below is for cases that the designed codon information is not available
group.names<-c()
if (!exists("codonDesigned")){
  group.names<-c("Intended_byAbndCut", "Unintended_byAbndCut")
} else{
  group.names<-c("Intended", "Unintended")
}

IntendedKey<-NULL

if (!exists("codonDesigned")){
  df20<-NULL
  
  names(df10)
  
  cdn.per.pos=codons.per.pos
  filtered.for.planned.vt<-NULL
  
  filtered.for.planned.vt<-df10 %>%
    filter(groupType=='sub') %>%
    filter(cdn.changes.by.Ted==1) %>%
    filter(str_detect(variant_newName,",")==FALSE) %>% ##single codon changes
    filter(!is.na(get(pDNASample.ct))) %>%
    dplyr::group_by(POS)%>%
    dplyr::mutate(rank_wi_pos=as.numeric(rank(-as.numeric(get(pDNASample.ct)), na.last = TRUE, ties.method ='random'))) %>%
    mutate(made.rank.cut=ifelse(rank_wi_pos < cdn.per.pos+3, 'madeCut','belowCut'))%>%
    arrange(POS)
  
  summary(filtered.for.planned.vt)
  table(filtered.for.planned.vt$rank_wi_pos)
  table(filtered.for.planned.vt$POS)
  
  
  table(as.character(filtered.for.planned.vt$change.to))
  
  filtered.for.planned.vt_m<-merge(filtered.for.planned.vt, codonTable, by.x='change.to', by.y="CODON", all=FALSE)
  
  #####################Done with codonDesign restruction
  
  
  filtered.as.planned.vt<-filtered.for.planned.vt%>%
    filter(made.rank.cut == "madeCut")
  
  
  write.csv(filtered.as.planned.vt, file=paste0(dir_wrt_variantCounts,"DataFileToDeduceCodonDesigns.csv"))
  
  IntendedKey<-as.vector(as.character(filtered.as.planned.vt$keys))
  
  
  cat("using abundance of detected variants to determine the said variant is planned or not\n")
} else {
  IntendedKey<-as.vector(as.character(codonDesigned$key))
  cat("using design file for Intended/Unintended classification\n")}

cat(paste0("Intended Variants: ", length(IntendedKey)))

cat("\nDone with: library design info\n Next up: process  cdonCounts for ORFCall version analysis, including ORF template sequence validation, clustering, AUC curves etc\n")


######################
library(reshape2)
###When need to rerun from saved file ..all_aggrgated, you will need restruct this file wf_0cdn
#names(wf_scaledToFraction_cdn_woWt)

######Very important file, should write out.

###sum up codon counts in .cdn files
wf_renamedCdn<-wf_0cdn
smpNMs<-names(wf_renamedCdn)[-c(1,2)]
wf_0cdn_colSum <- wf_renamedCdn %>%
  select(-POS,-AA_codon) %>%
  summarise_at(smpNMs,sum)
cdnCountSum<-data.frame(t(data.frame(wf_0cdn_colSum)))
names(cdnCountSum)<-c("cdnSum")
cdnCountSum$sample <- rownames(cdnCountSum)
cdnCountSum <- cdnCountSum %>%
  select("sample"=sample, "cdnSum"=cdnSum)

sampleNMs<-sampleAnnot$Experiment
levelOrder<-sampleNMs
cdnCountSum$sample<-factor(cdnCountSum$sample, levels=levelOrder)

write.csv(cdnCountSum, file=paste0(dir_wrt_cdnCounts,"ORFcallOfAnalyzeSatMut_cdnCountSumPerSample.csv"))

summary(cdnCountSum)

G1<-ggplot(cdnCountSum)+
  geom_col(mapping=aes(x=as.factor(sample), y=cdnSum))+
  coord_flip()
ggsave(filename=paste0(dir_wrt_cdnCounts,"ORFcallOfAnalyzeSatMut_cdnCountSumPerSample.pdf"))

G1

sum(cdnCountSum$cdnSum)

sum(cdnCountSum$cdnSum)*3/150


####SCALE CDN DATA
wf_tmpCdn<-wf_renamedCdn %>%
  select("POS"=POS, "CODON"=AA_codon, everything())

names(wf_tmpCdn)


ls_scaled<-FromRawToPercent_cbindRawCounts("rawFile"=wf_tmpCdn,
                                           "codonTb"=codonTable, 
                                           "colStart"=3,
                                           "colEnd"=dim(wf_tmpCdn)[2])


wf_scaledToFraction_cdn_woWt<-ls_scaled$wf_fraction_wo_wt  
wf_scaledToFraction_cdn_wWt<-ls_scaled$wf_fraction_w_wt #### for bookkeeping

wf_scaledToFraction_cdn_woWt_group_deltant<-wf_scaledToFraction_cdn_woWt %>%
  mutate(group=ifelse(is.element(str_c(POS,Vt_codon,sep="|"),IntendedKey),paste0(group.names[1],"_VT"), paste0(group.names[2],"_VT"))) %>%
  mutate(delta_nt=stringdist(Vt_codon,Wt_codon, method="hamming")) %>%
  mutate(group=ifelse(Wt_aa==Vt_aa, gsub("VT","WT",group),group)) %>%
  mutate(Vt_aa=ifelse(group ==paste0(group.names[1],"_WT"),"B", ifelse(group ==paste0(group.names[2],"_WT"),"b",as.character(Vt_aa)))) %>%
  mutate(Vt_aa_3G=ifelse(group ==paste0(group.names[1],"_WT"),"B", ifelse(group ==paste0(group.names[2],"_WT"),"b",ifelse(as.character(Vt_aa)=='X','X', 'others')))) %>%
  select(POS,Wt_codon,Wt_aa,Vt_codon,Vt_aa, group,delta_nt,Vt_aa_3G,everything())

summary(wf_scaledToFraction_cdn_woWt_group_deltant)

wf_scaledToFraction_cdn_wWt_group_deltant<-wf_scaledToFraction_cdn_wWt %>%
  mutate(group=ifelse(is.element(str_c(POS,Vt_codon,sep="|"),IntendedKey),group.names[1], group.names[2])) %>%
  mutate(delta_nt=stringdist(Vt_codon,Wt_codon, method="hamming")) %>%
  mutate(group=ifelse(Wt_aa==Vt_aa, gsub("VT","WT",group),group)) %>%
  mutate(Vt_aa=ifelse(group ==paste0(group.names[1],"_WT"),"B", ifelse(group ==paste0(group.names[2],"_WT"),"b",as.character(Vt_aa)))) %>%
  mutate(Vt_aa_3G=ifelse(group ==paste0(group.names[1],"_WT"),"B", ifelse(group ==paste0(group.names[2],"_WT"),"b",ifelse(as.character(Vt_aa)=='X','X', 'others')))) %>%
  select(POS,Wt_codon,Wt_aa,Vt_codon,Vt_aa, group,delta_nt,Vt_aa_3G,everything())


#names(wf_scaledToFraction_cdn_woWt_group_deltant)

nms2callpase<-names(wf_scaledToFraction_cdn_woWt_group_deltant)[9:dim(wf_scaledToFraction_cdn_woWt_group_deltant)[2]]
wf_collapsedToAA_FromScaledCdnFraction_andRaw_woWt <- wf_scaledToFraction_cdn_woWt_group_deltant %>%
  dplyr::group_by(POS,Wt_codon,Wt_aa,Vt_aa,group,Vt_aa_3G) %>%
  #summarise_at(.vars=nms2callpase,.funs=sum)
  summarise_at(.vars=nms2callpase, list(~sum(.)))

write.csv(wf_scaledToFraction_cdn_woWt_group_deltant,file=paste0(dir_wrt_cdnCounts,"File0_ORFcallOfAnalyzeSatMut_scaledToFraction_cdn_woWt_group_deltant.csv"))

write.csv(wf_scaledToFraction_cdn_wWt_group_deltant,file=paste0(dir_wrt_cdnCounts,"ORFcallOfAnalyzeSatMut_scaledToFraction_cdn_wWt_group_deltant.csv"))

write.csv(wf_0cdn_colSum,file=paste0(dir_wrt_cdnCounts,"ORFcallOfAnalyzeSatMut_perSampleDetectedCdnSum.csv"))


write.csv(wf_collapsedToAA_FromScaledCdnFraction_andRaw_woWt,file=paste0(dir_wrt_cdnCounts,"ORFcallOfAnalyzeSatMut_collapsedToAA_fromScaledCdnFractionAndRaw_woWt.csv"))

#########clustering
forClust<-NULL
tmp4clust<-NULL
forClust<-wf_scaledToFraction_cdn_woWt_group_deltant
tmp4clust<-forClust[which(str_detect(names(forClust),'_raw') & !str_detect(str_to_upper(str_sub(names(forClust),1,3)),'HOP'))]
if(dim(tmp4clust)[1]!=0){
  tmp<-scale(tmp4clust, center = TRUE)
  row.has.na<- apply(tmp, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  tmp<-data.frame(tmp[!row.has.na,])
  clustTmp12<- as.data.frame(t(tmp))
  d <- dist(clustTmp12, method='euclidean')
  fit <- hclust(d, method='ward.D')
  pdf(file=paste0(dir_wrt_cdnCounts,"cluster_cdnCt_raw.pdf"), width=6, height=6) 
  old.par=par(mfrow=c(1,1))
  plot(fit, cex=0.7)
  par(old.par)
  dev.off()
}
forClust<-NULL
tmp4clust<-NULL
forClust<-wf_scaledToFraction_cdn_woWt_group_deltant%>%
  filter(group %in% c(paste0(group.names[1],"_VT"), paste0(group.names[1],"_WT")))
tmp4clust<-forClust[which(str_detect(names(forClust),'_raw') & !str_detect(str_to_upper(str_sub(names(forClust),1,3)),'HOP'))]
if(dim(tmp4clust)[1]!=0){
  tmp<-scale(tmp4clust, center = TRUE)
  row.has.na<- apply(tmp, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  tmp<-data.frame(tmp[!row.has.na,])
  clustTmp12<- as.data.frame(t(tmp))
  d <- dist(clustTmp12, method='euclidean')
  fit <- hclust(d, method='ward.D')
  pdf(file=paste0(dir_wrt_cdnCounts,"cluster_cdnCt_Intended_raw.pdf"), width=6, height=6) 
  old.par=par(mfrow=c(1,1))
  plot(fit, cex=0.7)
  par(old.par)
  dev.off()
}


###frctn
#names(wf_scaledToFraction_cdn_woWt_group_deltant)
forClust<-NULL
tmp4clust<-NULL
tmp4clust<-wf_scaledToFraction_cdn_woWt_group_deltant[which(str_detect(names(wf_scaledToFraction_cdn_woWt_group_deltant),'_frctn') & !str_detect(str_to_upper(str_sub(names(wf_scaledToFraction_cdn_woWt_group_deltant),1,3)),'HOP'))]
tmp<-scale(tmp4clust, center = TRUE)
row.has.na<- apply(tmp, 1, function(x){any(is.na(x))})
sum(row.has.na)
tmp<-data.frame(tmp[!row.has.na,])
clustTmp12<- as.data.frame(t(tmp))
d <- dist(clustTmp12, method='euclidean')
fit <- hclust(d, method='ward.D')
pdf(file=paste0(dir_wrt_cdnCounts,"cluster_cdnCt_fraction.pdf"), width=6, height=6) 
old.par=par(mfrow=c(1,1))
plot(fit, cex=0.7)
par(old.par)
dev.off()


forClust<-NULL
tmp4clust<-NULL
forClust<-wf_scaledToFraction_cdn_woWt_group_deltant%>%
  filter(group %in% c(paste0(group.names[1],"_VT"), paste0(group.names[1],"_WT")))
tmp4clust<-forClust[which(str_detect(names(forClust),'_frctn') & !str_detect(str_to_upper(str_sub(names(forClust),1,3)),'HOP'))]

if(dim(tmp4clust)[1]!=0){
  tmp<-scale(tmp4clust, center = TRUE)
  row.has.na<- apply(tmp, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  tmp<-data.frame(tmp[!row.has.na,])
  clustTmp12<- as.data.frame(t(tmp))
  d <- dist(clustTmp12, method='euclidean')
  fit <- hclust(d, method='ward.D')
  
  pdf(file=paste0(dir_wrt_cdnCounts,"cluster_cdnCt_Intended_fraction.pdf"), width=6, height=6) 
  old.par=par(mfrow=c(1,1))
  plot(fit, cex=0.7)
  par(old.par)
  dev.off()
}
forClust<-NULL
tmp4clust<-NULL
########Distribution using .cdnCounts.

expNM_ct<-paste0(sampleAnnot[sampleAnnot$Sample==pDNASample,]$Experiment,"_raw")
expNM_frctn<-paste0(sampleAnnot[sampleAnnot$Sample==pDNASample,]$Experiment,"_frctn")

table(wf_scaledToFraction_cdn_woWt_group_deltant$group)

intended_cdnCounts_forAUC<-wf_scaledToFraction_cdn_woWt_group_deltant%>%
  filter(group %in% c(paste0(group.names[1],"_VT"), paste0(group.names[1],"_WT"))) %>%
  filter(expNM_ct>0)

#names(intended_cdnCounts_forAUC)

hg2<-NULL
hg1<-NULL

list2<-fracFunNoBlueNoLabel_colArg_generic_auc(fileIn=intended_cdnCounts_forAUC,columns=expNM_frctn,col='black',cex=0.5,alpha=1)
a<-list2$auc
hg2<-list2$gph
hg2=hg2+geom_abline(aes(intercept=0,slope=1),
                    colour='red', linetype=2, alpha=0.2, size=1)+
  theme(axis.title=element_text(size=10,face="bold"))

hg2


list1<-func_hist_from_raw(fileIn=intended_cdnCounts_forAUC,column=expNM_frctn, txhi=0,breaks=seq(0,20, by=0.1)) 
avg<-list1$avgReads
hg1<-list1$hist+
  theme(axis.title=element_text(size=10,face="bold"))
hg1

ggsave(filename = paste0(dir_wrt_cdnCounts,paste0(screenNM,"_histogram_intendedVt_frctn.pdf")),hg1)
ggsave(filename = paste0(dir_wrt_cdnCounts,paste0(screenNM,"_AUCcurve_intendedVt_frctn.pdf")),hg2)

detected_intCdns_cdnCounts<-paste(intended_cdnCounts_forAUC$POS,intended_cdnCounts_forAUC$Vt_codon, sep="|")

missing_cdns_cdnCounts<-IntendedKey[which(!(IntendedKey %in% detected_intCdns_cdnCounts))]

write(missing_cdns_cdnCounts, file=paste0(dir_wrt_cdnCounts,paste0(screenNM,"_AnalyzeSatMut_missingCodons_cdnCounts.txt")))

table(intended_cdnCounts_forAUC$POS)


expNM_ct<-NULL
expNM_frctn<-NULL

cat("Done with chunk4_cdnCounts\n")

####Below we are back to .variantCounts data
##df10 is a datafile annotated with 9 columns 
df10_1<-df10%>%
  mutate(group.plan=ifelse(keys %in% IntendedKey, group.names[1], group.names[2])) %>%
  mutate(delta_nt=stringdist(change.from,change.to, method = "hamming"))

table(df10_1$group.plan)

gplot<-ggplot(df10_1)+
  geom_point(mapping = aes(x=rank.by.pDNA, y=get(pDNASample.ct),col=group.plan), cex=1, pch=19, alpha=0.2)+
  scale_y_log10()+
  ggtitle("Annotated group.plan")+
  coord_cartesian(xlim=c(0,orf_size*63*2))
ggsave(gplot, filename = paste0(dir_wrt_variantCounts,"abundance_byGroup.pdf"))
gplot
gplot<-NULL

table(df10_1$group.plan,df10_1$groupType)

cdn.tb2c=codonTable[,1:2]

df11<-merge(cdn.tb2c,df10_1, by.x="CODON",by.y="change.from",all=TRUE)
names(df11)[1:2]=c("change.from","wt_aa")

table(df11$group.plan,df11$groupType)

df12<-merge(cdn.tb2c,df11, by.x="CODON",by.y="change.to",all=TRUE)
names(df12)[1:2]=c("change.to","vt_aa")

table(df12$group.plan,df12$groupType)


#names(df12)
#dim(df12)
df_variant2details_annotated_frctn_corrctd <- df12%>%
  mutate(vt_aa = ifelse(!is.na(vt_aa) & wt_aa==vt_aa, "B", as.character(vt_aa)))%>%
  select(lesion.center, variant_newName,group.plan,delta_nt, everything()) %>%
  arrange(lesion.center)

####################

write.csv(df_variant2details_annotated_frctn_corrctd, file=paste0(dir_wrt_variantCounts,"All_variantCounts_merged_annotatedWithDetails_frctn_corrctd.csv"))

cat("Done with: ORFcall version, saved df_variant2details_annotated_frctn_corrctd\n\nNext up: add key (species with planned vatiants) columns\n")

##################

#names(df_variant2details_annotated_frctn_corrctd)
df_variant2details_annotated_frctn_corrctd$X<-NULL
dim(df_variant2details_annotated_frctn_corrctd)
df0<-df_variant2details_annotated_frctn_corrctd %>%
  filter(change.from!="OUTBOUND" & !is.na(keys))
###tag PrimaryKey infomation

df0_1<-data.frame(cbind("keys"=as.character(df0[,10]),"iteration"=1:dim(df0)[1])) %>%
  mutate(keys_itr=str_c(keys,iteration,sep="/")) %>%
  select(keys_itr)

M=dim(df0_1)[1]

df0$prim.key.maxDelta<-
  apply(data.frame(df0_1), 1, FUN=function (x){ GetPrimaryKey_theKey_theDeltaNt("keys_itr"=x,"keyList"=IntendedKey,"wtCdns"=wtCdns)})

df0<-df0 %>%
  select(prim.key.maxDelta,everything())
#names(df0)

df<-df0

df$all.int.keys<-apply(data.frame(df[,which(names(df)=="prim.key.maxDelta")]),1,FUN=function(x) { str_split(as.character(x),"\\;")[[1]][1]})
#names(df)

df$the.key<-apply(data.frame(df[,which(names(df)=="prim.key.maxDelta")]),1,FUN=function(x) { str_split(as.character(x),"\\;")[[1]][2]})

df$the.key.deltaNT<-apply(data.frame(df[,which(names(df)=="prim.key.maxDelta")]),1,FUN=function(x) { str_split(as.character(x),"\\;")[[1]][3]})

df<-df %>%
  select(all.int.keys,the.key,the.key.deltaNT,prim.key.maxDelta,everything())

df[,1:4][df[,1:4]=='NA']<-NA


cat("Done with: adding all.int.keys,the.key,the.key.deltaNT \n\nNext up: add the.key.pos and the.key.vt.codon\n")

df_forKey<-df
P=dim(df_forKey)[1]
df_forKey_0 <- df_forKey %>%
  mutate(itr=1:P) %>%
  #mutate(the.key_itr=str_c(the.key,itr,sep='/')) %>%
  mutate(the.key_itr=ifelse(is.na(the.key), str_c("",itr,sep='/'), str_c(the.key,itr,sep='/'))) %>%
  select(the.key_itr)

df_forKey <- cbind(df_forKey,
                   "the.key.pos"=apply(df_forKey_0, 1, FUN = function(x) {
                     iteration=as.numeric(str_split(as.character(x),"/")[[1]][2])
                     the.key=as.character(str_split(as.character(x),"/")[[1]][1])
                     if((iteration %% (floor(P/1000)*100)) == 0){
                       cat(paste("Done: ", round((iteration/P) * 100,0),'%\n'))
                     } 
                     return(ifelse(str_detect(the.key,"\\|"),as.numeric(paste(t(data.frame(str_split(str_split(the.key,'\\,')[[1]],"\\|")))[,1],collapse = ",")),NA))}),
                   
                   "the.key.vt_codon"=apply(df_forKey_0, 1, FUN = function(x) {
                     iteration=as.numeric(str_split(as.character(x),"/")[[1]][2])
                     the.key=as.character(str_split(as.character(x),"/")[[1]][1])
                     if((iteration %% (floor(P/1000)*100)) == 0){
                       cat(paste("Done: ", round((iteration/P) * 100,0),'%\n'))
                     } 
                     return(ifelse(str_detect(the.key,"\\|"),paste(t(data.frame(str_split(str_split(the.key,'\\,')[[1]],"\\|")))[,2],collapse = ","),NA))}))

cat("\nDone with adding the.key.pos and the.key.vt_codon\n\nNext up: add the.key.wt_codon\n")

Q=dim(df_forKey)[1]
df_forKey_00 <- df_forKey %>%
  mutate(itr=1:Q) %>%
  #mutate(the.key_itr=str_c(the.key,itr,sep='/')) %>%
  mutate(the.key.pos_itr=ifelse(is.na(the.key.pos), str_c("",itr,sep='/'), str_c(as.character(the.key.pos),itr,sep="/"))) %>%
  select(the.key.pos_itr)

df_forKey <- cbind(df_forKey,        
                   "the.key.wt_codon"=apply(df_forKey_00,1,FUN = function(x) { 
                     iteration=as.numeric(str_split(as.character(x),"/")[[1]][2])
                     the.key.pos=str_split(as.character(x),"/")[[1]][1]
                     if((iteration %% (floor(Q/1000)*100)) == 0){
                       cat(paste("Done: ", round((iteration/Q) * 100,0),'%\n'))
                     } 
                     return(ifelse(!is.na(the.key.pos),paste(wtCdns[as.numeric(str_split(the.key.pos,'\\,')[[1]])],collapse = ","),NA))}))

#####
cat("\nDone with adding the.key.wt_codon\n")

#names(df_forKey)

sum(table(df_forKey$the.key.wt_codon))

sum(table(df_forKey$the.key.vt_codon))

m.df_forKey_vt<-merge(x=df_forKey, y=as.data.frame(codonTable[,1:2]), by.x="the.key.vt_codon", by.y="CODON", all=TRUE) %>%
  select(the.key.pos,the.key.vt_codon,"the.key.vt_aa"=AA, everything())

#names(m.df_forKey_vt)

table(m.df_forKey_vt$vt_aa)
table(m.df_forKey_vt$the.key.vt_aa)

table(m.df_forKey_vt$the.key.vt_aa,m.df_forKey_vt$vt_aa)

m.df_forKey_vt_wt<-merge(m.df_forKey_vt, as.data.frame(codonTable[,1:2]), by.x="the.key.wt_codon", by.y="CODON", all=TRUE) %>%
  select(the.key.pos,the.key.wt_codon,"the.key.wt_aa"=AA, the.key.vt_codon, the.key.vt_aa, everything())

table(m.df_forKey_vt_wt$the.key.vt_aa,m.df_forKey_vt_wt$vt_aa)


df_fullSet_31LeadCols<-m.df_forKey_vt_wt %>%
  #mutate(the.key.deltaNT=stringdist(the.key.wt_codon,the.key.vt_codon,method = "hamming"))%>%
  select(the.key.pos,the.key.wt_codon,the.key.wt_aa, the.key.vt_codon, the.key.vt_aa,
         the.key.deltaNT,  everything())

table(is.na(m.df_forKey_vt_wt$the.key))

#names(df_fullSet_31LeadCols)

sum(table(df_fullSet_31LeadCols$the.key.vt_aa,df_fullSet_31LeadCols$vt_aa))

table(df_fullSet_31LeadCols$the.key.vt_aa,df_fullSet_31LeadCols$vt_aa)

df_fullSet_35LeadCols<-df_fullSet_31LeadCols

write.csv(df_fullSet_35LeadCols, file=paste0(dir_wrt_variantCounts,"df_fullSet_35LeadCols.csv"))

cat("Done with: adding the.key columns, df_fullSet_35LeadCols.csv is saved\n")

##############rename
df_fullSet_35LeadCols_renamed<-df_fullSet_35LeadCols
nms1=names(df_fullSet_35LeadCols_renamed)
for(i in 1:dim(df_fullSet_35LeadCols_renamed)[2]){
  if(str_sub(nms1[i],1,6)=='Sample'){
    nms1=gsub(str_sub(nms1[i],1,8), sampleAnnot[sampleAnnot$Sample==str_sub(nms1[i],1,8),]$Experiment, nms1)
  }
}
names(df_fullSet_35LeadCols_renamed)<-nms1
cat("Done with: 35-col annotations. Ready to create File1 and File2 for handoff\n")

table(df_fullSet_35LeadCols_renamed$groupType)

df_toHandoff <- df_fullSet_35LeadCols_renamed

#names(df_toHandoff)

df_toHandoff<-df_toHandoff[,-which(str_detect(names(df_toHandoff),"_corrctd")==TRUE)]
#df_toHandoff<-df_toHandoff[,-which(str_detect(names(df_toHandoff),".len")==TRUE)]

names(df_toHandoff)

df_toHandoff_sub_indel_21leadCols<-df_toHandoff[,c(1:7,11,12,18,19,21:23,27,29:dim(df_toHandoff)[2])]

df_toHandoff_intended_11leadCols<-df_toHandoff[,c(1:6,11,12,27,35:dim(df_toHandoff)[2])] %>%
  mutate(variant.by.aa=str_c(the.key.wt_aa, the.key.pos, the.key.vt_aa, sep=''))%>%
  filter(group.plan==group.names[1]) %>%
  mutate(mutationType=ifelse(as.character(the.key.vt_aa) == 'B' | (as.character(the.key.wt_codon) != as.character(the.key.vt_codon) & as.character(the.key.wt_aa) == as.character(the.key.vt_aa) & group.plan==group.names[1]) , 'silent', ifelse(the.key.vt_aa == 'X', 'nonsense', 'missense'))) %>%
  select("variant.by.nt.at.codon.level"=variant_newName,
         variant.by.aa,
         "POS"=the.key.pos,
         "Wt_codon"=the.key.wt_codon,
         "Vt_codon"=the.key.vt_codon,
         "Wt_aa"=the.key.wt_aa,
         "Vt_aa"=the.key.vt_aa,
         "delta_nt"=the.key.deltaNT,
         "is.the.whole.mol.planned"=group.plan,
         "del_ins_or_sub"=groupType,
         mutationType,everything()
  )%>%
  mutate_at(.vars="POS",list(~as.numeric(.)))%>% ##conforming with new dplyr funs(name = f(.)) to list(name = ~f(.))
  arrange(POS)

df_toHandoff_intended_11leadCols<-renameSample2Experiment(df_toHandoff_intended_11leadCols,sampleAnnot)

table(df_toHandoff_intended_11leadCols$groupType)

table(df_toHandoff_intended_11leadCols$POS)

table(df_toHandoff_intended_11leadCols$is.the.whole.mol.planned)

table(df_toHandoff_intended_11leadCols$mutationType)

table(df_toHandoff_intended_11leadCols$Vt_aa)

df_toHandoff_sub_indel_23leadCols<-renameSample2Experiment(df_toHandoff_sub_indel_21leadCols,sampleAnnot) %>%
  mutate(Vt.aa.By.Ted=gsub("I:->","II>",Vt.aa.By.Ted)) %>%
  mutate(Vt.aa.By.Ted=gsub(":-",">-",Vt.aa.By.Ted)) %>%
  mutate(Vt.aa.By.Ted=gsub("II>","I:->",Vt.aa.By.Ted))%>%
  mutate(frame.shift=str_detect(Vt.aa.By.Ted,"FS"))%>%
  mutate(mutationType=ifelse(as.character(the.key.vt_aa) == 'B' | (as.character(the.key.wt_codon) != as.character(the.key.vt_codon) & as.character(the.key.wt_aa) == as.character(the.key.vt_aa) & group.plan==group.names[1]) , 'silent', ifelse(the.key.vt_aa == 'X', 'nonsense', 'missense'))) %>%
  mutate(pos.of.cdnChange.if.whole.mol.is.intended =ifelse(group.plan==group.names[1],the.key.pos,''))%>%
  select(
    #"variant.by.nt.at.codon.level"=variant_newName, 
    "variant.description.by.unit.of.wtCdn.cdnPos"=variant_newName, 
    "variant.description.by.unit.of.wt.nt_Ted"=variant,                                         
    "variant.description.by.unit.of.cdnPos"=keys,
    "num.whole.mol.nt.changes_Ted"=changes,                                                  
    "num.whole.mol.cdn.changes_upto.a.stop_Ted"=cdn.changes.by.Ted,                             
    "list.whole.mol.all.cdn.changes_upto.a.stop_Ted"=Vt.Cdn.By.Ted,                             
    "list.whole.mol.all.aa.changes_upto.a.stop_Ted"=Vt.aa.By.Ted, 
    StdNom.Vt.aa.By.Ted,
    "is.the.whole.mol.planned"=group.plan, 
    pos.of.cdnChange.if.whole.mol.is.intended,
    "del_ins_or_sub"=groupType,
    mutationType,
    "is.there.frame.shift"=frame.shift,
    "subList.all.intended.cdnChanges.ifAny.detected.in.the.mol"= all.int.keys, 
    "pick.the.highest.intendedCdnDeltaNtChange.cdnPos"=the.key.pos,
    "pick.the.highest.intendedCdnDeltaNtChange.wtCdn"=the.key.wt_codon,
    "pick.the.highest.intendedCdnDeltaNtChange.wtAA"=the.key.wt_aa,
    "pick.the.highest.intendedCdnDeltaNtChange.vtCdn"=the.key.vt_codon,                         
    "pick.the.highest.intendedCdnDeltaNtChange.vtAA"=the.key.vt_aa,                                
    "pick.the.highest.intendedCdnDeltaNtChange.deltaNT"=the.key.deltaNT,                           
    "num.insertion.nt"=nt_ins,                                                   
    "num.deletion.nt"=nt_del,                                                         
    "the.first.mutation.event.occurs.at.ntPos"=ntPOS,                                              
    "the.first.mutation.event.occurs.at.nthNT.in.a.cdn"=nth_nt_in_cdn,                           
    everything()
  ) %>%
  mutate_at(.vars="pick.the.highest.intendedCdnDeltaNtChange.cdnPos",list(~as.numeric(.)))%>% 
  arrange(pick.the.highest.intendedCdnDeltaNtChange.cdnPos)


table(df_toHandoff_sub_indel_23leadCols$is.the.whole.mol.planned,df_toHandoff_sub_indel_23leadCols$mutationType)


generic_samples<-sampleAnnot$Sample
experiments<-sampleAnnot$Experiment

experiments[which(as.vector(generic_samples)==as.character(pDNASample))]

table(df_toHandoff_sub_indel_23leadCols$is.the.whole.mol.planned,df_toHandoff_sub_indel_23leadCols[[paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct')]])

table(df_toHandoff_sub_indel_23leadCols[df_toHandoff_sub_indel_23leadCols$is.the.whole.mol.planned %in% c('Intended','Intended_byAbndCut'),][[paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct')]])

###########Using raw count to look at library distribution. Utimately, we should using .fraction for the real distribution.
raw.ct.of.intended<-as.data.frame(table(df_toHandoff_sub_indel_23leadCols[df_toHandoff_sub_indel_23leadCols$is.the.whole.mol.planned %in% c('Intended','Intended_byAbndCut'),][[paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct')]]))
names(raw.ct.of.intended)<-c("Count","Freq")

intended.vt.detected=sum(raw.ct.of.intended$Freq)
intended.vt.attempted=length(IntendedKey)

perct=round(intended.vt.detected/intended.vt.attempted,3)

int.rawCt.pDNA.hist=ggplot(raw.ct.of.intended)+
  geom_col(mapping=aes(x=as.numeric(Count), y=as.numeric(Freq)))+
  ggtitle(paste0("pDNA library: detected = ",intended.vt.detected,"(",perct, ") , out of attempted = ",intended.vt.attempted))

ggsave(plot = int.rawCt.pDNA.hist, filename =paste0(dir_wrt_variantCounts,"IntendedVt_rawCountHist_pDNA.pdf"))
#####

if(!is.null(twistSample)){
  raw.ct.of.intended.twist<-as.data.frame(table(df_toHandoff_sub_indel_23leadCols[df_toHandoff_sub_indel_23leadCols$is.the.whole.mol.planned=='Intended',][[paste0(experiments[which(as.vector(generic_samples)==as.character(twistSample))],'.ct')]]))
  names(raw.ct.of.intended.twist)<-c("Count","Freq")
  
  intended.vt.detected.twist=sum(raw.ct.of.intended.twist$Freq)
  perct.twist=round(intended.vt.detected.twist/intended.vt.attempted,3)
  
  int.rawCt.twist.hist=ggplot(raw.ct.of.intended.twist)+
    geom_col(mapping=aes(x=as.numeric(Count), y=as.numeric(Freq)))+
    ggtitle(paste0("Twist library: detected = ",intended.vt.detected.twist,"(",perct.twist, ") , out of attempted = ",intended.vt.attempted))
  
  ggsave(plot = int.rawCt.twist.hist, filename =paste0(dir_wrt_variantCounts,"IntendedVt_rawCountHist_TwistPool.pdf"))
}

#########Lib distribution using pDNA fraction

df_int_for_distr.pdna<-df_toHandoff_intended_11leadCols %>%
  filter(!is.na(get(paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct_frctn'))))


detected_int_variantCounts<-paste(df_int_for_distr.pdna$POS,df_int_for_distr.pdna$Vt_codon, sep="|")

missing_vt_variantCounts.pdna<-IntendedKey[which(!(IntendedKey %in% detected_int_variantCounts))]

write(missing_vt_variantCounts.pdna, file=paste0(dir_wrt_variantCounts,paste0(screenNM,"_AnalyzeSatMut_missingCodons_pDNA_variantCounts.txt")))

num.missing.vt.pdna<-length(missing_vt_variantCounts.pdna)

num.missing.vt.pdna

pcnt.miss.pdna<-round(length(missing_vt_variantCounts.pdna)/length(IntendedKey),3)

names(df_int_for_distr.pdna)

summary(df_int_for_distr.pdna[[paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct_frctn')]])

hgg2<-NULL
hgg1<-NULL
list2<-NULL
list1<-NULL

list2<-fracFunNoBlueNoLabel_colArg_generic_auc(fileIn=df_int_for_distr.pdna,columns=c(paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct_frctn')),col='black',cex=0.5,alpha=1)
a<-list2$auc
hgg2<-list2$gph
hgg2=hgg2+geom_abline(aes(intercept=0,slope=1),
                      colour='red', linetype=2, alpha=0.2, size=1)+
  theme(axis.title=element_text(size=10,face="bold"))+
  ggtitle(paste0("pDNA library: missing = ",num.missing.vt.pdna,"(",pcnt.miss.pdna, ") , out of attempted = ",length(IntendedKey)))

hgg2

list1<-func_hist_from_raw(fileIn=df_int_for_distr.pdna,column=c(paste0(experiments[which(as.vector(generic_samples)==as.character(pDNASample))],'.ct_frctn')), txhi=0,breaks=seq(0,20, by=0.1)) 
avg<-list1$avgReads
hgg1<-list1$hist+
  theme(axis.title=element_text(size=10,face="bold"))+
  ggtitle(paste0("pDNA library: missing = ",num.missing.vt.pdna,"(",pcnt.miss.pdna, ") , out of attempted = ",length(IntendedKey)))
hgg1

ggsave(filename = paste0(dir_wrt_variantCounts,paste0(screenNM,"_histogram_intendedVt_pDNA_ct_frctn.pdf")),hgg1)
ggsave(filename = paste0(dir_wrt_variantCounts,paste0(screenNM,"_AUCcurve_intendedVt_pDNA_ct_frctn.pdf")),hgg2)


#######




write.csv(df_toHandoff_intended_11leadCols, file=paste0(dir_wrt_variantCounts, paste0("File1_",screenNM,"_AnalyzeSatMut_theIntended_11leadCols.csv")))

write.csv(df_toHandoff_sub_indel_23leadCols, 
          file=paste0(dir_wrt_variantCounts, paste0("File2_",screenNM,"_AnalyzeSatMut_fullSet_23leadCols.csv")))

cat("\nDone with: File1 and File2 \n")

#######File1 replicates

File1<-read_csv(paste0(dir_wrt_variantCounts,paste0("File1_",screenNM,"_AnalyzeSatMut_theIntended_11leadCols.csv")))
File1$X1<-NULL
names(File1)


df_clust<-File1[which(
  str_detect(names(File1),'\\.ct$'))] 


names(df_clust)

library(corrplot)



row.has.na<- apply(df_clust, 1, function(x){any(is.na(x))})
sum(row.has.na)
df_clust<-data.frame(df_clust[!row.has.na,])

  png(paste0(dir_wrt_variantCounts,"File1_replicate_correlation_hClust_sortByReplicateName.png"),width = 11, height = 9, units = 'in', res=500)


oldpar<-par(mfrow=c(1,1), ps=8)
cor1<-cor(df_clust, method='pearson')
corrplot(cor1,method="number", title = "By replicate names", 
         order='alphabet',hclust.method = 'ward.D2',
         type='upper',
         mar=c(0,0,4,0), number.cex = 1, 
         number.digits=2)
par(oldpar)


  dev.off()

  png(paste0(dir_wrt_variantCounts,"File1_replicate_correlation_hClust_sortByPearson.png"),width = 11, height = 9, units = 'in', res=500)


oldpar<-par(mfrow=c(1,1), ps=8)
corrplot(cor1,method="number", title = "By Pearson", 
         order='hclust',hclust.method = 'ward.D2',
         type='upper',
         mar=c(0,0,4,0), number.cex = 1, 
         number.digits=3)

par(oldpar)

  dev.off()

###################################More correlation graphs
  tmp<-scale(df_clust, center = TRUE)
  row.has.na<- apply(tmp, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  tmp<-data.frame(tmp[!row.has.na,])
  clustTmp12<- as.data.frame(t(tmp))
  d <- dist(clustTmp12, method='euclidean')
  fit <- hclust(d, method='ward.D')
  pdf(file=paste0(dir_wrt_variantCounts,"cluste_File1_rawCt.pdf"), width=6, height=6) 
  old.par=par(mfrow=c(1,1))
  plot(fit, cex=0.7)
  par(old.par)
  dev.off()
##################


df_clust_scaled<-scale(df_clust, center = TRUE)
mat<-as.matrix(df_clust_scaled)

getUpperTri <- function(cor.mat){
  cor.mat[lower.tri(cor.mat)] <- NA
  return(cor.mat)
}

getBothTri <- function(cor.mat){
  # cor.mat[lower.tri(cor.mat)] <- NA
  return(cor.mat)
}

reorderCormat <- function(cor.mat){
  dist.mat <- as.dist((1-cor.mat)/2)
  hc <- hclust(dist.mat)
  cor.mat <-cor.mat[hc$order,hc$order]
}




#cor.df <- reshape2::melt(getUpperTri(reorderCormat(cor(mat))),na.rm=TRUE,value.name="correlation",varnames=c("well_1","well_2"))


cor.df <- reshape2::melt(getBothTri(reorderCormat(cor(mat))),na.rm=TRUE,value.name="correlation",varnames=c("well_1","well_2"))
names(cor.df)


hmapRepl<-
  ggplot(cor.df, aes(well_2, well_1, fill=correlation))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limit=c(-1,1), space="Lab", name="Pearson\nCorrelation")+
  geom_text(aes(well_2, well_1,label=round(correlation,3)),cex=5, angle=45)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=1, size=7, hjust=1),
        axis.text.y=element_text(size=7))+
  coord_fixed()+
  labs(x="",y="")



  ggsave(filename=paste0(dir_wrt_variantCounts,"File1_replicateHeatMap.pdf"), plot=hmapRepl)


hmapRepl







