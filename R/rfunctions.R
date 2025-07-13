#' Runs the PTMoreR Shiny web application.
#' @export
PTMoreR_app <- function(appDir=system.file('PTMoreRapp', package='PTMoreR'),host="127.0.0.1", port =8989) {
  shiny::runApp(appDir=appDir,host=host,port=port)
}

#' Pre-alignment.
#' @export
prealign<-function(data,datatype="Normal",central.amino.acid="STY",
                   label.of.modification="#",width=7,
                   species.fasta.file="10116.fasta"){
  library(Biostrings)
  library(stringi)
  library(stringr)
  library(tidyr)
  dataread<-data
  if(ncol(dataread)==1){
    datareadx<-dataread
  }else{
    datareadx<-dataread[,2,drop=F]
  }
  origidatatypex<-datatype
  Peptides<-vector()
  if(origidatatypex=="MaxQuant"){
    for(i in 1:nrow(datareadx)){
      pep1<-datareadx[[1]][i]
      Peptidesi1<-strsplit(gsub("[^0-9.]", ";", pep1),";")[[1]]
      Peptidesi2<-Peptidesi1[Peptidesi1!=""]
      for(ii in 1:length(Peptidesi2)){
        if(as.numeric(Peptidesi2[ii])>=0.75){
          pep1<-gsub(paste0("\\(",Peptidesi2[ii],"\\)"),label.of.modification,pep1)#"#"
        }else{
          pep1<-gsub(paste0("\\(",Peptidesi2[ii],"\\)"),"",pep1)
        }
      }
      Peptides[i]<-pep1
    }
    dataconvert<-data.frame(AnnotatedPeps=Peptides,stringsAsFactors = FALSE)
  }else if(origidatatypex=="Spectronaut"){
    phosphoindex<-grep("\\[Phospho \\(STY\\)\\]",datareadx[[1]], perl = TRUE)
    uploaddata1<-datareadx[phosphoindex,]
    Peptidesx<-gsub("_","",uploaddata1, perl = TRUE)
    Peptidesx<-gsub("\\[Phospho \\(STY\\)\\]",label.of.modification,Peptidesx, perl = TRUE)#"#"
    Peptidesx2<-str_replace_all(Peptidesx,"\\[.*?\\]","")
    dataconvert<-data.frame(AnnotatedPeps=Peptidesx2,stringsAsFactors = FALSE)
  }else{
    dataconvert<-datareadx
  }
  dataconvertindex<-grep(label.of.modification,dataconvert[[1]])
  uploaddata2<-dataconvert[dataconvertindex,,drop=FALSE]
  if(ncol(uploaddata2)==1){
    uploaddata1<-datareaddq<-uploaddata2
  }else{
    uploaddata1<-datareaddq<-uploaddata2[,2,drop=F]
  }
  centralres1<-strsplit(label.of.modification,";|")[[1]]
  centralres<-centralres1[centralres1!=""]
  centralres2<-paste(centralres,collapse = "|")
  uploaddata1$Stripped.pep<-gsub(paste0("_|",centralres2),"",datareaddq[[1]], perl = TRUE)
  EGindex<-lapply(datareaddq[[1]],function(x){
    xx4<-gregexpr(centralres[1],x)[[1]]
    xx3<-gregexpr(centralres2,x)[[1]]
    xx5<-unlist(lapply(xx4,function(x) which(x==xx3)))
    xx1<-1:length(xx3)
    xx6<-as.numeric(xx3)-xx1
    xx2<-paste(xx6[xx5],collapse = ";")
    xx2
  })
  EGindex1<-lapply(datareaddq[[1]],function(x){
    xx3<-gregexpr(centralres2,x)[[1]]
    xx1<-1:length(xx3)
    xx6<-as.numeric(xx3)-xx1
    xx2<-paste(xx6,collapse = ";")
    xx2
  })
  centeranjisuan<-lapply(datareaddq[[1]],function(x){
    pepi<-strsplit(gsub(centralres2,"",x),"")[[1]]
    xx4<-gregexpr(centralres[1],x)[[1]]
    xx3<-gregexpr(centralres2,x)[[1]]
    xx5<-unlist(lapply(xx4,function(x) which(x==xx3)))
    xx1<-1:length(xx3)
    xx6<-as.numeric(xx3)-xx1
    xx2<-paste(pepi[xx6[xx5]],collapse = ";")
    xx2
  })
  uploaddata1$Pep.main.index<-unlist(EGindex)
  uploaddata1$Pep.all.index<-unlist(EGindex1)
  uploaddata1$Center.amino.acid<-unlist(centeranjisuan)
  colnames(uploaddata1)<-c("Pep.upload","Stripped.pep","Pep.main.index","Pep.all.index","Center.amino.acid")
  datafasta<-readAAStringSet(species.fasta.file)
  n_data_fasta<-length(datafasta@ranges@NAMES)
  pro_seqdf<-as.data.frame(datafasta)
  pro_seqdfnames<-unlist(lapply(rownames(pro_seqdf),function(x) strsplit(x,"\\|")[[1]][2]))
  if(sum(is.na(pro_seqdfnames))>n_data_fasta/2){
    pro_seqdfnames<-unlist(lapply(rownames(pro_seqdf),function(x) strsplit(x,"\\ ")[[1]][1]))
  }
  danlength<-width
  seqseqall<-vector()
  proidall<-vector()
  proidindexall<-vector()
  PRO.CombinedID<-vector()
  for(i in 1:nrow(uploaddata1)){
    seqindex1<-grep(uploaddata1$Stripped.pep[i],pro_seqdf$x, perl = TRUE)
    seqindex3<-as.numeric(strsplit(uploaddata1$Pep.main.index[i],";")[[1]])
    seqseqall1<-vector()
    proidindexall1<-vector()
    procomb<-vector()
    if(length(seqindex1)>0 & length(seqindex3)>0){
      for(k in 1:length(seqindex1)){
        seqindex2<-stri_locate_first(pattern = uploaddata1$Stripped.pep[i], pro_seqdf$x[seqindex1[k]], fixed = TRUE)[[1]][1]#[,1]
        seqnchar<-nchar(pro_seqdf$x[seqindex1[k]])
        indexjian<-unlist(lapply(seqindex2, function(x){x+seqindex3-1}))
        seqseq<-vector()
        for(j in 1:length(indexjian)){
          indexjian1<-indexjian[j]-danlength
          indexjian2<-indexjian[j]+danlength
          if(indexjian1<=0){
            xhx1<-paste(rep("_",abs(indexjian1)+1),collapse ="")
            xhx2<-stri_sub(pro_seqdf$x[seqindex1[k]],from = 0,to=indexjian2)
            xhx3<-paste0(xhx1,xhx2)
          }
          else if(indexjian2>seqnchar){
            xhx1<-paste(rep("_",(indexjian2-seqnchar)),collapse="")
            xhx2<-stri_sub(pro_seqdf$x[seqindex1[k]],from = indexjian1,to=seqnchar)
            xhx3<-paste0(xhx2,xhx1)
          }
          else{
            xhx3<-stri_sub(pro_seqdf$x[seqindex1[k]],from = indexjian1,to=indexjian2)
          }
          seqseq[j]<-xhx3
        }
        seqseqall1[k]<-paste(seqseq,collapse = ";")
        proidindexall1[k]<-paste(indexjian,collapse = ";")
        procomb[k]<-paste(paste0(pro_seqdfnames[seqindex1[k]],"_",strsplit(uploaddata1$Center.amino.acid[i],";")[[1]],indexjian),collapse = ";")
      }
      seqseqall[i]<-paste(seqseqall1,collapse = "::")
      proidall[i]<-paste(pro_seqdfnames[seqindex1],collapse = "::")
      proidindexall[i]<-paste(proidindexall1,collapse = "::")
      PRO.CombinedID[i]<-paste(procomb,collapse = "::")
    }else{
      seqseqall[i]<-"No Match"
      proidall[i]<-"No Match"
      proidindexall[i]<-"No Match"
      PRO.CombinedID[i]<-"No Match"
    }
  }
  uploaddata1$Seqwindows<-seqseqall
  uploaddata1$PRO.from.Database<-proidall
  uploaddata1$PROindex.from.Database<-proidindexall
  uploaddata1$PRO.CombinedID<-PRO.CombinedID
  datareaddq<-unique(uploaddata1)
  datareaddq1<-separate_rows(datareaddq,6:9,sep = "::")
  datareaddq2<-unique(separate_rows(datareaddq1,3:9,sep = ";"))
  colnames(datareaddq2)<-c("Uploaded Peptides","Peptide Skeleton",
                          "Main Modification Site in Uploaded Peptides",
                          "All Modification Site in Uploaded Peptides",
                          "Center amino acid","Aligned Standard Peptides",
                          "Protein ID from Database",
                          "Modified Amino Acid Site in Protein Sequence",
                          "Combined Protein Identifier")
  datareaddq3<-datareaddq2[,c(1,7,5,8,9,6,2,3,4)]
  datareaddq3
}
