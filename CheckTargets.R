### This script will take a tabular blast result (-outfmt 6) with your target exons as the query and (draft) genome as the subject
###  and do the following:
###   1. Find potential paralogs by identifying genes with good matches from more than one contig/chromosome segments, or multiple positions within them, in the reference genome
###   2. Find potentially missing genes by identifying genes with no or few good matches in the reference genome
###   3. Find genes spanning huge introns
###  then
###   4. Create sets of simple segment plots to allow for visual inspection of blast alignments of genes failing and passing the above checks
###   5. Output files listing names of targets passing or failing each check, and a file of 'clean' targets passing all checks

### Arguments:
###  --blast_file <tab-delimited blast result, target=query, genome=subject>
###  --min_pident < % identity below which to discard blast matches prior to running checks >
###  --max_intron_length < maximum allowable length of any given intron in a gene >

###  --min_fragment_length < minimum length of a blast hit, ignore anything smaller >
###  --output_prefix < prefix to name results files >

## not yet implemented
###  --max_intron_prop < maximum allowable proportion (0-1) of a supercontig's length consisting of introns >

GetCoverageStats<-function(data){
  coverage_summary<-data.frame()
  for(i in unique(data$qseqid)){
    temp<-data %>% 
      filter(qseqid==i)%>%
      arrange(qstart) 
    
    for(m in unique(temp$sseqid)){
      temp2<-temp %>% filter(sseqid==m)
      target_seq<-seq(1:unique(temp2$qlen))
      for(f in 1:nrow(temp2)){
        target_seq<-c(target_seq,seq(temp2[f,]$qstart,temp2[f,]$qend))
      }
      mycounts<-table(target_seq)-1 #subtract 1 because we've appended to the original sequence
      paralog_coverage<-length(mycounts[mycounts>1])/unique(temp2$qlen)
      full_coverage<-length(mycounts[mycounts>0])/unique(temp2$qlen)
      unique_coverage<-length(mycounts[mycounts==1])/unique(temp2$qlen)
      coverage_summary<-rbind(coverage_summary,
                              data.frame(qseqid=unique(temp2$qseqid),
                                         sseqid=unique(temp2$sseqid),
                                         paralog_coverage=paralog_coverage,
                                         full_coverage=full_coverage,
                                         unique_coverage=unique_coverage))
    }
  }
  return(coverage_summary)
}

FindIntrons<-function(data,max_intron_length){
  introns_out<-data.frame()
  introns_flag_out<-data.frame()
  for(i in unique(data$qseqid)){
    
    temp<-data %>% 
      filter(qseqid==i)
    
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
      mutate(sum_introns=sum(intron_length))%>%
      select(qseqid)%>%unique()
    introns_flag_out<-rbind(introns_flag_out,introns_flag)
  }
  return(list(intron_details=introns_out,targets_with_large_introns=introns_flag_out))
}

PlotTargets<-function(data,output_prefix){
  pdf(paste0(output_prefix,"_plots.pdf"),width=10)
  for(i in unique(data$qseqid)){
    
    temp<-data %>% 
      filter(qseqid==i)
    
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
      ggtitle(i)
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
  
}

CheckTargets<-function(blast_file,min_pident,min_fragment_length,max_intron_length,output_prefix){
  dat<-as_tibble(read.table(blast_file,header=T))
  dat<-dat %>%
    filter(pident >= min_pident,
           length >= min_fragment_length)
  print(dat)
  PlotTargets(dat,output_prefix)
  
  cov_stats<-GetCoverageStats(dat)
  write.table(cov_stats,paste0(output_prefix,"_CoverageStats.txt"),quote = F,row.names = F,col.names = TRUE)
  
  intron_stats<-FindIntrons(data=dat,max_intron_length)
  write.table(intron_stats$intron_details,paste0(output_prefix,"_IntronStats.txt"),quote = F,row.names = F,col.names = TRUE)
  write.table(intron_stats$targets_with_large_introns,paste0(output_prefix,"_TargetsWithLargeIntrons.txt"),quote = F,row.names = F,col.names = TRUE)
}

#### DONE DEFINING FUNCTIONS
suppressMessages(suppressWarnings(require(optparse)))

p <- OptionParser(usage=" This script will take a tabular blast result (-outfmt 6) with your target exons as the query and (draft) genome as the subject\n
   ***WITH A HEADER LINE**\n
  and do the following:\n
   1. Find potential paralogs by identifying genes with good matches from more than one contig/chromosome segments in the reference genome\n
   2. Find potentially missing genes by identifying genes with no or few good matches in the reference genome\n
   3. Find genes spanning huge introns\n
  then\n
   4. Create sets of simple segment plots to allow for visual inspection of blast alignments of genes failing and passing the above checks\n
   5. Output files listing names of targets passing or failing each check, and a file of 'clean' targets passing all checks\n
  Run using Rscript, e.g.\n
  Rscript CheckTargets.R --blast_file blastn_targets_to_genome.txt --min_pident 80 --max_intron_length 10000")
# Add a positional argument
p <- add_option(p, c("--blast_file"), help="<tab-delimited blast result, target=query, genome=subject>",type="character")
p <- add_option(p, c("--min_fragment_length"), help="< minimum length of a blast hit, ignore anything smaller >",type="numeric")
p <- add_option(p, c("--min_pident"), help="% identity below which to discard blast matches prior to running checks >",type="numeric")
p <- add_option(p, c("--max_intron_length"), help="< maximum allowable length of any given intron in a gene >",type="numeric")
# p <- add_option(p, c("--max_intron_percent"), help="< maximum allowable percentage of a supercontig's length consisting of introns >",type="numeric")
p <- add_option(p, c("--output_prefix"), help="< prefix to name results files >",type="character")

# parse
args<-parse_args(p)

suppressMessages(suppressWarnings(require(tidyverse,quietly =TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(ggrepel,quietly =TRUE,warn.conflicts=FALSE)))

## RUN
CheckTargets(blast_file=args$blast_file,
             min_pident=args$min_pident,
             min_fragment_length=args$min_fragment_length,
             max_intron_length=args$max_intron_length,
             output_prefix=args$output_prefix)
             