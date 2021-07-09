
  setwd("~/GitHub/TargetSupercontigs")

library(tidyverse)
library(ggrepel)

rhodo_mm<-as_tibble(read.table("test_data/blastn_markerminer_to_Rgrier.txt",header = T))
rhodo_mm

<<<<<<< Updated upstream
pdf("ggplot_segments_Rgrier.pdf",width=10)
for(i in unique(rhodo_mm$qseqid)){
  
  temp<-rhodo_mm %>% filter(qseqid==i)
  
  introns <-  temp %>% group_by(qseqid,sseqid)%>% 
    mutate(p_sstart=ifelse(sstart<send,sstart,send),
           p_send=ifelse(p_sstart==sstart,send,sstart),
           orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
    arrange(p_sstart,.by_group = TRUE) %>%
    summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
              intron=(shift_sstart-shift_send),
              qpos=c(NA,qend),shift_qstart=c(qstart,NA),
              orientation=c(orientation,NA),
              .groups = "keep")%>%
    filter(intron>1000,intron<1e6)%>%
    na.exclude()
  
  p<-ggplot(temp)+
    geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=sseqid),
                 show.legend = FALSE)+
    facet_wrap(~sseqid,scales = "free")+
    xlim(c(0,unique(temp$qlen)))+
    theme_bw()+
    labs(x="Target position",y="Genome position")+
    ggtitle(unique(temp$qseqid))
  if(nrow(introns)>0){
    p<-p+geom_rect(data=introns,
                   aes(xmin=-Inf,xmax=Inf,ymin=shift_sstart,ymax=shift_send),
                   colour="darkgreen",alpha=0.05,size=NA)+
      geom_label_repel(data=introns,
                       aes(x=qpos,
                           y=shift_sstart-(intron/2),
                           label=intron),fill="white",
                       force_pull = 0.1,
                       size=2,
                       min.segment.length=unit(5,"mm"),
                       box.padding = 0.05,label.padding=0.07,
                       label.r=0.01,seed=1234,
                       max.overlaps=15)
  }
  
  
  print(p)
}
dev.off()
dev.off()
=======
library(tidyverse)
library(ggrepel)
>>>>>>> Stashed changes

rhodo_mm<-as_tibble(read.table("test_data/blastn_markerminer_to_Rgrier.txt",header = T))
rhodo_mm

<<<<<<< Updated upstream
min_fragment_length=100
max_intron_length=1000
### create filtering data.frame

## NOTE: need to distinguish true introns from distances between paralogs on the same chromosome
## Best done by first identifying same-chromosome paralogs, then only doing FindIntrons on the rest


FindIntrons<-function(data,min_fragment_length,max_intron_length){
  introns_out<-data.frame()
  introns_flag_out<-data.frame()
  for(i in unique(data$qseqid)){
    
    temp<-data %>% 
      filter(qseqid==i,
             length >= min_fragment_length)
    
    introns <-  temp %>% group_by(qseqid,sseqid)%>% 
      mutate(p_sstart=ifelse(sstart<send,sstart,send),
             p_send=ifelse(p_sstart==sstart,send,sstart),
             orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
      arrange(p_sstart,.by_group = TRUE) %>%
      summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
                intron_length=(shift_sstart-shift_send),
                orientation=c(orientation,NA),
                .groups = "keep")%>%
      na.exclude() %>%
      select(qseqid,sseqid,intron_length) 
    
    ## write to data.frame
    introns_out<-rbind(introns_out,introns)
    
    introns_flag<-introns %>%
      filter(intron_length>=max_intron_length)%>%
      ungroup()%>%
      select(qseqid)%>%unique()
    introns_flag_out<-rbind(introns_flag_out,introns_flag)
  }
  return(list(intron_details=introns_out,targets_with_large_introns=introns_flag_out))
}
introns_out
introns_flag_out

=======
pdf("ggplot_segments_Rgrier.pdf",width=10)
for(i in unique(rhodo_mm$qseqid)){
  
  temp<-rhodo_mm %>% filter(qseqid==i)
  
  introns <-  temp %>% group_by(qseqid,sseqid)%>% 
    mutate(p_sstart=ifelse(sstart<send,sstart,send),
           p_send=ifelse(p_sstart==sstart,send,sstart),
           orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
    arrange(p_sstart,.by_group = TRUE) %>%
    summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
              intron=(shift_sstart-shift_send),
              qpos=c(NA,qend),shift_qstart=c(qstart,NA),
              orientation=c(orientation,NA),
              .groups = "keep")%>%
    filter(intron>1000,intron<1e6)%>%
    na.exclude()
  
  p<-ggplot(temp)+
    geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=sseqid),
                 show.legend = FALSE)+
    facet_wrap(~sseqid,scales = "free")+
    xlim(c(0,unique(temp$qlen)))+
    theme_bw()+
    labs(x="Target position",y="Genome position")+
    ggtitle(unique(temp$qseqid))
  if(nrow(introns)>0){
    p<-p+geom_rect(data=introns,
                   aes(xmin=-Inf,xmax=Inf,ymin=shift_sstart,ymax=shift_send),
                   colour="darkgreen",alpha=0.05,size=NA)+
      geom_label_repel(data=introns,
                 aes(x=qpos,
                     y=shift_sstart-(intron/2),
                     label=intron),fill="white",
                 force_pull = 0.1,
                 size=2,
                 min.segment.length=unit(5,"mm"),
                 box.padding = 0.05,label.padding=0.07,
                 label.r=0.01,seed=1234,
                 max.overlaps=15)
  }

  
  print(p)
}
dev.off()
dev.off()


min_fragment_length=100
max_intron_length=1000
### create filtering data.frame

## NOTE: need to distinguish true introns from distances between paralogs on the same chromosome
## Best done by first identifying same-chromosome paralogs, then only doing FindIntrons on the rest


FindIntrons<-function(data,min_fragment_length,max_intron_length){
  introns_out<-data.frame()
  introns_flag_out<-data.frame()
  for(i in unique(data$qseqid)){
    
    temp<-data %>% 
      filter(qseqid==i,
             length >= min_fragment_length)
    
    introns <-  temp %>% group_by(qseqid,sseqid)%>% 
      mutate(p_sstart=ifelse(sstart<send,sstart,send),
             p_send=ifelse(p_sstart==sstart,send,sstart),
             orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
      arrange(p_sstart,.by_group = TRUE) %>%
      summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
                intron_length=(shift_sstart-shift_send),
                orientation=c(orientation,NA),
                .groups = "keep")%>%
      na.exclude() %>%
      select(qseqid,sseqid,intron_length) 
    
    ## write to data.frame
    introns_out<-rbind(introns_out,introns)
    
    introns_flag<-introns %>%
      filter(intron_length>=max_intron_length)%>%
      ungroup()%>%
      select(qseqid)%>%unique()
    introns_flag_out<-rbind(introns_flag_out,introns_flag)
  }
  return(list(intron_details=introns_out,targets_with_large_introns=introns_flag_out))
}
introns_out
introns_flag_out

>>>>>>> Stashed changes
## find same-chromosome paralogs
coverage_summary<-data.frame()
for(i in unique(rhodo_mm$qseqid)[120:129]){
  
  temp<-rhodo_mm %>% 
    filter(qseqid==i)%>%
    arrange(qstart)
  temp
  
  for(m in unique(temp$sseqid)){
    # print(m)
    temp2<-temp %>% filter(sseqid==m)
    # print(temp2)
    target_seq<-seq(1:unique(temp2$qlen))
    for(f in 1:nrow(temp2)){
      target_seq<-c(target_seq,seq(temp2[f,]$qstart,temp2[f,]$qend))
    }
    # print(table(target_seq)[table(target_seq)>1])
    same_chrom_overlap<-length(table(target_seq)[table(target_seq)>2])/unique(temp2$qlen)
    full_coverage<-length(table(target_seq)[table(target_seq)>1])/unique(temp2$qlen)
    unique_coverage<-length(table(target_seq)[table(target_seq)==2])/unique(temp2$qlen)
    coverage_summary<-rbind(coverage_summary,
                            data.frame(sseqid=unique(temp2$sseqid),
                                       qseqid=unique(temp2$qseqid),
                                       same_chromosome_overlap=same_chrom_overlap,
                                       full_coverage=full_coverage,
                                       unique_coverage=unique_coverage))
  }
  ## write to data.frame
}
<<<<<<< Updated upstream
=======


# #####
# 
# i=unique(rhodo_mm$qseqid)[125]
# pdf("ggplot_segments_Rgrier.pdf",width=10)
# for(i in unique(rhodo_mm$qseqid)){
#   
#   temp<-rhodo_mm %>% filter(qseqid==i)
#   
#   introns <-  temp %>% group_by(qseqid,sseqid)%>% 
#     mutate(p_sstart=ifelse(sstart<send,sstart,send),
#            p_send=ifelse(p_sstart==sstart,send,sstart))%>% 
#     arrange(p_sstart,.by_group = TRUE) %>%
#     summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),intron=(shift_sstart-shift_send),
#               qpos=c(NA,qend),shift_qstart=c(qstart,NA),
#               .groups = "keep")%>%
#     filter(intron>1000,intron<1e6)%>%
#     na.exclude()
#   
#   p<-ggplot(temp)+
#     geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=sseqid,size=pident/100),
#                  show.legend = FALSE)+
#     # geom_label(aes(x=qend,y=send,label=abs(qend-qstart)))+
#     geom_segment(data=introns,
#                  aes(x=qpos,xend=shift_qstart,
#                      y=shift_send,yend=shift_sstart),
#                  lty=2,colour="darkgreen")+
#     geom_label(data=introns,
#                aes(x=qpos,
#                    y=shift_sstart-(intron/2),
#                    label=intron),fill="white")+
#     facet_wrap(~sseqid,scales = "free")+
#     xlim(c(0,unique(temp$qlen)))+
#     theme_bw()+
#     labs(x="Target position",y="Genome position")+
#     ggtitle(unique(temp$qseqid))
#   if(nrow(temp)>1){
#     p<-p + scale_size_binned(range=c(1,2))
#   }
#   print(p)
# }
# dev.off()


##

# i=unique(rhodo_mm$qseqid)[125]
# pdf("ggplot_segments_Rgrier.pdf",width=10)
# for(i in unique(rhodo_mm$qseqid)){
#   
#   temp<-rhodo_mm %>% filter(qseqid==i)
#   
#   introns <-  temp %>% group_by(qseqid,sseqid)%>% 
#     mutate(p_sstart=ifelse(sstart<send,sstart,send),
#            p_send=ifelse(p_sstart==sstart,send,sstart))%>% 
#     arrange(p_sstart,.by_group = TRUE) %>%
#     summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),intron=(shift_sstart-shift_send),
#               qpos=c(NA,qend),shift_qstart=c(qstart,NA),
#               .groups = "keep")%>%
#     filter(intron>1000,intron<1e6)%>%
#     na.exclude()
#   
#   p<-ggplot(temp)+
#     geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=sseqid),
#                  show.legend = FALSE)+
#     geom_segment(data=introns,
#                  aes(x=qpos,xend=shift_qstart,
#                      y=shift_send,yend=shift_sstart),
#                  lty=2,colour="darkgreen")+
#     geom_label(data=introns,
#                aes(x=qpos,
#                    y=shift_sstart-(intron/2),
#                    label=intron),fill="white")+
#     facet_wrap(~sseqid,scales = "free")+
#     xlim(c(0,unique(temp$qlen)))+
#     theme_bw()+
#     labs(x="Target position",y="Genome position")+
#     ggtitle(unique(temp$qseqid))
#   print(p)
# }
# dev.off()
# dev.off()





















cinerea_prot<-as_tibble(read.table("test_data/tblastn_mega353_BYO_longest_Rhodo_match_length1000_not_in_Erica134_or_MarkerMiner_e1e-6_maxt50k_to_cinerea_k96_merged_deduped-8_minlen500_withHeader.txt",
                         header = T,stringsAsFactors = T))
cinerea_prot
cinerea_prot$qseqid<-gsub("cds_1","cds",cinerea_prot$qseqid)
cinerea_prot$sseqid<-as.character(cinerea_prot$sseqid)
cinerea_prot$qseqid<-as.character(cinerea_prot$qseqid)

cinerea_prot_good<-cinerea_prot[cinerea_prot$pident>70,]
dim(cinerea_prot_good)

cinerea_prot_good_summary<-cinerea_prot_good %>% 
  group_by(qseqid,sseqid) %>%
  summarise(count=n(),.groups = "keep",
            matchlen=sum(abs(qend-qstart+1)),
            qlen=unique(qlen),
            prop_qlen=sum(matchlen)/qlen,
            mean_pident=mean(pident),
            sd_pidend=sd(pident))
cinerea_prot_good_summary

######################################################################
pdf("test_data/Contig_segment_plots.pdf")
for(t in unique(cinerea_prot_good_summary$qseqid)){
  temp<-cinerea_prot_good[cinerea_prot_good$qseqid==t,]
  outx<-vector()
  outy<-vector()
  f<-length(unique(temp$sseqid))
  minident<-min(temp$pident)
  maxident<-max(temp$pident)
  sc<-100-minident
  plot( c(-50,max(temp$qlen)),c(0.75,f), type="n",xlab="Position",ylab="Contig of origin",
        main=paste(t,"\nmin_pident =",round(minident,0),"max_pident =",round(maxident,0),sep=" "),
        yaxt="n")
  for(i in 1:(nrow(temp))){
    x1<-c(temp[i,"qstart"],temp[i,"qend"])
    outx<-c(x1)
    ID<-which(unique(temp$sseqid) %in% temp[i,"sseqid"])
    ident<-temp[i,"pident"]
    segments(x0=outx$qstart,x1=outx$qend,y0=ID,y1=ID,lwd=4,col=rgb(0,0,0,1-((maxident-ident)/sc)))
    ID2<-unique(temp$sseqid)[ID]
    text(x=-25,y=ID,labels=ID2)
  }
}
dev.off()

## ggplot version
temp



######################################################################
pdf("test_data/Contig_segment_plots_inverted.pdf")
for(t in unique(cinerea_prot_good_summary$sseqid)){
  temp<-cinerea_prot_good[cinerea_prot_good$sseqid==t,]
  outx<-vector()
  f<-length(unique(temp$qseqid))
  minident<-min(temp$pident)
  maxident<-max(temp$pident)
  p1<-min(temp[,c("sstart","send")])
  p2<-max(temp[,c("sstart","send")])
  truelength<-abs(p2-p1)
  tpos<-p1+((p2-p1)/2)
  sc<-100-minident
  plot(c(p1,p2),
       c(0.75,f), 
       yaxt="n",
       type="n",xlab="Position",ylab="Transcript ID",
       main=paste0("contig: ",t,"\nmin_pident = ",round(minident,0),", max_pident = ",round(maxident,0),"\ntrue length = ",truelength))
  for(i in 1:(nrow(temp))){
    outx<-unlist(c(temp[i,"sstart"],temp[i,"send"]))
    ID<-which(unique(temp$qseqid) %in% temp[i,"qseqid"])
    ident<-temp[i,"pident"]
    if(!is.na((maxident-ident)/sc)){
      myalpha = (1.1-((maxident-ident)/sc))/1.1
      mycolour=rgb(0,0,0,myalpha)
    }else{
      mycolour=rgb(0,0,0,1)
    }
    segments(x0=outx[1],x1=outx[2],y0=ID,y1=ID,lwd=4,col=mycolour)
    ID2<-unique(temp$qseqid)[ID]
    text(x=tpos,y=ID-0.05,labels=ID2)
  }
}
dev.off()
dev.off()
dev.off()
>>>>>>> Stashed changes


# #####
# 
# i=unique(rhodo_mm$qseqid)[125]
# pdf("ggplot_segments_Rgrier.pdf",width=10)
# for(i in unique(rhodo_mm$qseqid)){
#   
#   temp<-rhodo_mm %>% filter(qseqid==i)
#   
#   introns <-  temp %>% group_by(qseqid,sseqid)%>% 
#     mutate(p_sstart=ifelse(sstart<send,sstart,send),
#            p_send=ifelse(p_sstart==sstart,send,sstart))%>% 
#     arrange(p_sstart,.by_group = TRUE) %>%
#     summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),intron=(shift_sstart-shift_send),
#               qpos=c(NA,qend),shift_qstart=c(qstart,NA),
#               .groups = "keep")%>%
#     filter(intron>1000,intron<1e6)%>%
#     na.exclude()
#   
#   p<-ggplot(temp)+
#     geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=sseqid,size=pident/100),
#                  show.legend = FALSE)+
#     # geom_label(aes(x=qend,y=send,label=abs(qend-qstart)))+
#     geom_segment(data=introns,
#                  aes(x=qpos,xend=shift_qstart,
#                      y=shift_send,yend=shift_sstart),
#                  lty=2,colour="darkgreen")+
#     geom_label(data=introns,
#                aes(x=qpos,
#                    y=shift_sstart-(intron/2),
#                    label=intron),fill="white")+
#     facet_wrap(~sseqid,scales = "free")+
#     xlim(c(0,unique(temp$qlen)))+
#     theme_bw()+
#     labs(x="Target position",y="Genome position")+
#     ggtitle(unique(temp$qseqid))
#   if(nrow(temp)>1){
#     p<-p + scale_size_binned(range=c(1,2))
#   }
#   print(p)
# }
# dev.off()


##

# i=unique(rhodo_mm$qseqid)[125]
# pdf("ggplot_segments_Rgrier.pdf",width=10)
# for(i in unique(rhodo_mm$qseqid)){
#   
#   temp<-rhodo_mm %>% filter(qseqid==i)
#   
#   introns <-  temp %>% group_by(qseqid,sseqid)%>% 
#     mutate(p_sstart=ifelse(sstart<send,sstart,send),
#            p_send=ifelse(p_sstart==sstart,send,sstart))%>% 
#     arrange(p_sstart,.by_group = TRUE) %>%
#     summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),intron=(shift_sstart-shift_send),
#               qpos=c(NA,qend),shift_qstart=c(qstart,NA),
#               .groups = "keep")%>%
#     filter(intron>1000,intron<1e6)%>%
#     na.exclude()
#   
#   p<-ggplot(temp)+
#     geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=sseqid),
#                  show.legend = FALSE)+
#     geom_segment(data=introns,
#                  aes(x=qpos,xend=shift_qstart,
#                      y=shift_send,yend=shift_sstart),
#                  lty=2,colour="darkgreen")+
#     geom_label(data=introns,
#                aes(x=qpos,
#                    y=shift_sstart-(intron/2),
#                    label=intron),fill="white")+
#     facet_wrap(~sseqid,scales = "free")+
#     xlim(c(0,unique(temp$qlen)))+
#     theme_bw()+
#     labs(x="Target position",y="Genome position")+
#     ggtitle(unique(temp$qseqid))
#   print(p)
# }
# dev.off()
# dev.off()





















cinerea_prot<-as_tibble(read.table("test_data/tblastn_mega353_BYO_longest_Rhodo_match_length1000_not_in_Erica134_or_MarkerMiner_e1e-6_maxt50k_to_cinerea_k96_merged_deduped-8_minlen500_withHeader.txt",
                                   header = T,stringsAsFactors = T))
cinerea_prot
cinerea_prot$qseqid<-gsub("cds_1","cds",cinerea_prot$qseqid)
cinerea_prot$sseqid<-as.character(cinerea_prot$sseqid)
cinerea_prot$qseqid<-as.character(cinerea_prot$qseqid)

cinerea_prot_good<-cinerea_prot[cinerea_prot$pident>70,]
dim(cinerea_prot_good)

cinerea_prot_good_summary<-cinerea_prot_good %>% 
  group_by(qseqid,sseqid) %>%
  summarise(count=n(),.groups = "keep",
            matchlen=sum(abs(qend-qstart+1)),
            qlen=unique(qlen),
            prop_qlen=sum(matchlen)/qlen,
            mean_pident=mean(pident),
            sd_pidend=sd(pident))
cinerea_prot_good_summary

######################################################################
pdf("test_data/Contig_segment_plots.pdf")
for(t in unique(cinerea_prot_good_summary$qseqid)){
  temp<-cinerea_prot_good[cinerea_prot_good$qseqid==t,]
  outx<-vector()
  outy<-vector()
  f<-length(unique(temp$sseqid))
  minident<-min(temp$pident)
  maxident<-max(temp$pident)
  sc<-100-minident
  plot( c(-50,max(temp$qlen)),c(0.75,f), type="n",xlab="Position",ylab="Contig of origin",
        main=paste(t,"\nmin_pident =",round(minident,0),"max_pident =",round(maxident,0),sep=" "),
        yaxt="n")
  for(i in 1:(nrow(temp))){
    x1<-c(temp[i,"qstart"],temp[i,"qend"])
    outx<-c(x1)
    ID<-which(unique(temp$sseqid) %in% temp[i,"sseqid"])
    ident<-temp[i,"pident"]
    segments(x0=outx$qstart,x1=outx$qend,y0=ID,y1=ID,lwd=4,col=rgb(0,0,0,1-((maxident-ident)/sc)))
    ID2<-unique(temp$sseqid)[ID]
    text(x=-25,y=ID,labels=ID2)
  }
}
dev.off()

## ggplot version
temp



######################################################################
pdf("test_data/Contig_segment_plots_inverted.pdf")
for(t in unique(cinerea_prot_good_summary$sseqid)){
  temp<-cinerea_prot_good[cinerea_prot_good$sseqid==t,]
  outx<-vector()
  f<-length(unique(temp$qseqid))
  minident<-min(temp$pident)
  maxident<-max(temp$pident)
  p1<-min(temp[,c("sstart","send")])
  p2<-max(temp[,c("sstart","send")])
  truelength<-abs(p2-p1)
  tpos<-p1+((p2-p1)/2)
  sc<-100-minident
  plot(c(p1,p2),
       c(0.75,f), 
       yaxt="n",
       type="n",xlab="Position",ylab="Transcript ID",
       main=paste0("contig: ",t,"\nmin_pident = ",round(minident,0),", max_pident = ",round(maxident,0),"\ntrue length = ",truelength))
  for(i in 1:(nrow(temp))){
    outx<-unlist(c(temp[i,"sstart"],temp[i,"send"]))
    ID<-which(unique(temp$qseqid) %in% temp[i,"qseqid"])
    ident<-temp[i,"pident"]
    if(!is.na((maxident-ident)/sc)){
      myalpha = (1.1-((maxident-ident)/sc))/1.1
      mycolour=rgb(0,0,0,myalpha)
    }else{
      mycolour=rgb(0,0,0,1)
    }
    segments(x0=outx[1],x1=outx[2],y0=ID,y1=ID,lwd=4,col=mycolour)
    ID2<-unique(temp$qseqid)[ID]
    text(x=tpos,y=ID-0.05,labels=ID2)
  }
}
dev.off()
dev.off()
dev.off()



@ @@