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
###  --max_intron_percent < maximum allowable percentage of a supercontig's length consisting of introns >

GetCoverageStats<-function(data){
  coverage_summary_chromosome_aware<-data.frame()
  coverage_summary_chromosome_UNaware<-data.frame()
  for(i in unique(data$qseqid)){
    temp<-data %>% 
      filter(qseqid==i)%>%
      arrange(qstart) 
   
    # per-target+chromosome combo 
    for(m in unique(temp$sseqid)){
      temp2<-temp %>% filter(sseqid==m)
      target_seq<-seq(1:unique(temp2$qlen))
      for(f in 1:nrow(temp2)){
        target_seq<-c(target_seq,seq(temp2[f,]$qstart,temp2[f,]$qend))
      }
      mycounts<-table(target_seq)-1 #subtract 1 because we've appended to the original sequence
      paralog_percent<-length(mycounts[mycounts>1])/unique(temp2$qlen)
      mean_paralogy_rate<-mean(mycounts)
      full_percent<-length(mycounts[mycounts>0])/unique(temp2$qlen)
      unique_percent<-length(mycounts[mycounts==1])/unique(temp2$qlen)
      missing_percent<-length(mycounts[mycounts==0])/unique(temp2$qlen)
      coverage_summary_chromosome_aware<-rbind(coverage_summary_chromosome_aware,
                              data.frame(qseqid=unique(temp2$qseqid),
                                         sseqid=unique(temp2$sseqid),
                                         mean_paralogy_rate=round(mean_paralogy_rate,1),
                                         paralog_percent=paralog_percent,
                                         full_percent=full_percent,
                                         unique_percent=unique_percent,
                                         missing_percent=missing_percent))
    }
    
    # per-target only (ignoring which chromosome(s) it matches to)
    target_seq2<-seq(1:unique(temp$qlen))
    for(f in 1:nrow(temp)){
      target_seq2<-c(target_seq2,seq(temp[f,]$qstart,temp[f,]$qend))
    }
    mycounts<-table(target_seq2)-1 #subtract 1 because we've appended to the original sequence
    paralog_percent<-length(mycounts[mycounts>1])/unique(temp$qlen)
    # weighted_paralog_percent<-(length(mycounts[mycounts>1])*mean(mycounts[mycounts>1])/unique(temp$qlen))
    mean_paralogy_rate<-mean(mycounts)
    full_percent<-length(mycounts[mycounts>0])/unique(temp$qlen)
    unique_percent<-length(mycounts[mycounts==1])/unique(temp$qlen)
    missing_percent<-length(mycounts[mycounts==0])/unique(temp$qlen)
    coverage_summary_chromosome_UNaware<-rbind(coverage_summary_chromosome_UNaware,
                                             data.frame(qseqid=unique(temp$qseqid),
                                                        mean_paralogy_rate=round(xend,1),
                                                        paralog_percent=paralog_percent,
                                                        full_percent=full_percent,
                                                        unique_percent=unique_percent,
                                                        missing_percent=missing_percent))
  }
  # round and convert to percentage
  coverage_summary_chromosome_aware<-data.frame(cbind(coverage_summary_chromosome_aware[,1:3],apply(coverage_summary_chromosome_aware[,4:ncol(coverage_summary_chromosome_aware)],2,function(x) round(x,3)*100)))
  coverage_summary_chromosome_UNaware<-data.frame(cbind(coverage_summary_chromosome_UNaware[,1:2],apply(coverage_summary_chromosome_UNaware[,3:ncol(coverage_summary_chromosome_UNaware)],2,function(x) round(x,3)*100)))
  
  coverage_summary_nchroms<-coverage_summary_chromosome_aware %>%
    group_by(qseqid) %>%
    summarise(n_chroms=n(),
              .groups = "keep")
  coverage_summary_chromosome_UNaware<-left_join(coverage_summary_chromosome_UNaware,coverage_summary_nchroms,"qseqid")
  
  return(list(coverage_summary_chromosome_aware=coverage_summary_chromosome_aware,coverage_summary_chromosome_unaware=coverage_summary_chromosome_UNaware))
}

# GetCoverageStats(rhodo_mm[rhodo_mm$qseqid%in%unique(rhodo_mm$qseqid)[1:20],])

# FindIntrons<-function(data,max_intron_length){
#   introns_out<-data.frame()
#   introns_flag_out<-data.frame()
#   for(i in unique(data$qseqid)){
#     
#     temp<-data %>% filter(qseqid==i)
#     
#     introns <-  temp %>% group_by(qseqid,sseqid)%>% 
#       mutate(p_sstart=ifelse(sstart<send,sstart,send),
#              p_send=ifelse(p_sstart==sstart,send,sstart),
#              orientation=ifelse(p_sstart==sstart,"forward","reverse"))%>% 
#       arrange(p_sstart,.by_group = TRUE) %>%
#       summarise(shift_sstart=c(p_sstart,NA),shift_send=c(NA,p_send),
#                 intron_length=(shift_sstart-shift_send),
#                 orientation=c(orientation,NA),
#                 intron_position_on_target=c(qend,NA),
#                 .groups = "keep")%>%
#       na.exclude() %>%
#       select(qseqid,sseqid,intron_length,intron_position_on_target) 
#     ## write to data.frame
#     introns_out<-rbind(introns_out,introns)
#     
#     ## flag targets with an intron exceeding max_intron_length
#     introns_flag<-introns %>%
#       filter(intron_length>=max_intron_length)%>%
#       ungroup()%>%
#       mutate(sum_introns=sum(intron_length))%>%
#       select(qseqid)%>%
#       unique()
#     ## write to data.frame
#     introns_flag_out<-rbind(introns_flag_out,introns_flag)
#   }
# 
#   return(list(intron_details=introns_out,targets_with_large_introns=introns_flag_out))
# }

FindIntrons<-function(data,max_intron_length,max_intron_percent){
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
                intron_length=shift_sstart-shift_send,
                orientation=c(orientation,NA),
                intron_position_on_target=c(qend,NA),
                supercontig_length=max(p_send)-min(p_sstart),
                .groups = "keep")%>%
      na.exclude() %>%
      select(qseqid,sseqid,intron_length,intron_position_on_target,supercontig_length)
    ## write to data.frame
    introns_out<-rbind(introns_out,introns)
    
    introns_flag<-introns %>%
      summarise(longest_intron_length=max(intron_length),
                total_intron_length=sum(intron_length),
                largest_intron_percent=max(sum(intron_length)/unique(supercontig_length))*100,
                exceeds_max_intron_length=ifelse(any(intron_length>=max_intron_length),TRUE,FALSE),
                exceeds_max_intron_percent=ifelse(largest_intron_percent>=max_intron_percent,TRUE,FALSE),
                summed_supercontig_length=unique(supercontig_length),
                .groups="keep") %>%
      ungroup() %>%
      select(qseqid,sseqid,summed_supercontig_length,
             longest_intron_length,total_intron_length,exceeds_max_intron_length,
             largest_intron_percent,exceeds_max_intron_percent) %>%
      unique()
    
    introns_flag_out<-rbind(introns_flag_out,introns_flag)
  }
  return(list(intron_details=introns_out,targets_with_intron_flags=introns_flag_out))
}

PlotTargets<-function(data,output_prefix,min_display_intron,max_display_intron){
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
      filter(intron>min_display_intron,intron<max_display_intron)%>%
      na.exclude()
    
    p<-ggplot(temp)+
      geom_segment(aes(y=sstart,yend=send,x=qstart,xend=qend,colour=bitscore),
                   show.legend = T,size=2.5)+
      scale_colour_viridis_c("Sequence\nsimilarity")+
      facet_wrap(~sseqid,scales = "free")+
      xlim(c(0,unique(temp$qlen)))+
      theme_bw()+
      labs(x="Target position",y="Genome position")+
      ggtitle(i)+
      coord_flip()
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

CheckTargets<-function(blast_file,
                       min_pident,min_fragment_length,
                       max_intron_length,max_intron_percent,
                       output_prefix,
                       min_display_intron,max_display_intron){
  dat<-as_tibble(read.table(blast_file,header=T))
  dat<-dat %>%
    filter(pident >= min_pident,
           length >= min_fragment_length)
  
  cov_stats<-GetCoverageStats(dat)
  write.table(cov_stats$coverage_summary_chromosome_aware,paste0(output_prefix,"_CoverageStats_PerChromosome.txt"),quote = F,row.names = F,col.names = TRUE)
  write.table(cov_stats$coverage_summary_chromosome_unaware,paste0(output_prefix,"_CoverageStats_AcrossChromosomes.txt"),quote = F,row.names = F,col.names = TRUE)
  
  intron_stats<-FindIntrons(data=dat,max_intron_length,max_intron_percent)
  write.table(intron_stats$intron_details,paste0(output_prefix,"_IntronStats.txt"),quote = F,row.names = F,col.names = TRUE)
  write.table(intron_stats$targets_with_intron_flags,paste0(output_prefix,"_IntronFlags.txt"),quote = F,row.names = F,col.names = TRUE)
  
  PlotTargets(dat,output_prefix,min_display_intron,max_display_intron)
  
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
p <- add_option(p, c("--blast_file"), help="<Required: tab-delimited blast result, target=query, genome=subject>",type="character")
p <- add_option(p, c("--min_fragment_length"), help="<Required: minimum length of a blast hit, ignore anything shorter>",type="numeric")
p <- add_option(p, c("--min_pident"), help="<Required: minimum % identity of a blast hit, ignore anything lower>",type="numeric")
p <- add_option(p, c("--max_intron_length"), help="<Required: length threshold of the largest intron in a gene; flag any gene exceeding>",type="numeric")
p <- add_option(p, c("--max_intron_percent"), help="<Required: percentage of a supercontig's length consisting of introns; flag any gene exceeding>",type="numeric")
p <- add_option(p, c("--output_prefix"), help="<prefix to name results files; defaults to CheckTargets_results>",type="character",default = "CheckTargets_results")
p <- add_option(p, c("--min_display_intron"), help="<intron length above which to annotate on plots; default 1kb>",type="numeric",default=1000)
p <- add_option(p, c("--max_display_intron"), help="<don't annotate introns longer than this; default 1Mb>",type="numeric",default=1e6)

# parse
args<-parse_args(p)

suppressMessages(suppressWarnings(require(tidyverse,quietly=TRUE,warn.conflicts=FALSE)))
suppressMessages(suppressWarnings(require(ggrepel,quietly=TRUE,warn.conflicts=FALSE)))

## RUN
CheckTargets(blast_file=args$blast_file,
             min_pident=args$min_pident,
             min_fragment_length=args$min_fragment_length,
             max_intron_length=args$max_intron_length,
             max_intron_percent=args$max_intron_percent,
             output_prefix=args$output_prefix,
             min_display_intron=args$min_display_intron,
             max_display_intron=args$max_display_intron)

sink(paste0(args$output_prefix,".warnings.txt"))
print(warnings())
sink()



             