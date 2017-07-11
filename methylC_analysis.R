#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
#test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Sample name must be given", call.=FALSE)
}

#define plotting function for features
plot_features <- function(df){
  require(ggplot2)
  ggplot(df, aes(x=Bin)) + geom_line(aes(y=mCG), color="dodgerblue4", size=0.8) +
         geom_line(aes(y=mCHG), color="olivedrab", size=0.8) +
         geom_line(aes(y=mCHH), color="hotpink4", size=0.8) +
         theme(panel.background=element_blank(), panel.grid=element_blank(),
         axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
         axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
         legend.position="none", axis.line=element_line(color="black")) +
         ylab("Percent methylation") + xlab( "" ) +
         scale_y_continuous(limits=c(0,1), expand=c(0,0),
         breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
         geom_vline(xintercept=20, linetype="longdash", color="grey55") +
         geom_vline(xintercept=40, linetype="longdash", color="grey55")
}

#define plotting function for gene-mC correlations
plot_correlations <- function(df,mC,mC_color){
  require(ggplot2)
    ggplot(df, aes(genes, mC)) + geom_point(color=mC_color) +
            theme(panel.background=element_blank(), panel.grid=element_blank(),
            axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
            axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
            legend.position="none", axis.line=element_line(color="black")) +
            ylab("CG methylation level") + xlab("Number of genes") +
            scale_x_continuous(expand=c(0,0)) +
            stat_smooth(aes(genes,mC), color="black", method="lm", se=F)
}

#Total methylation
if(file.exists("results/total_weighted_mC.tsv")){
  print("Plotting genome-wide total methylation")
  df <- read.table("results/total_weighted_mC.tsv",header=T,sep="\t")
  plot <- ggplot(df[1:3,],aes(x=Context,y=Weighted_mC,fill=Context)) +
          geom_bar(stat = "identity") +
          theme(panel.background=element_blank(), panel.grid=element_blank(),
          axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
          axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
          legend.position="none", axis.line=element_line(color="black")) +
          ylab("Percent methylation") + xlab( "" ) +
          scale_y_continuous(limits=c(0,1), expand=c(0,0),
          breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
          scale_fill_manual("",values=c("dodgerblue4","olivedrab","hotpink4"))
  filename=paste("figures_tables/Fvesca_", args[1], "_total_methylation.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=3, useDingbats=F)
}

#Per-site methylation
filename=paste("results/", args[1], "_per_site_boxplot.tsv", sep="")
if(!file.exists(filename)){
  print("Combining per-site methylation data")
  x <- data.frame()
  context <- c('CG','CGA','CGT','CGC','CGG','CHG','CAG','CTG','CCG','CHH',
               'CAA','CAT','CAC','CTA','CTT','CTC','CCA','CCT','CCC')
  for(seq in context){
      input=paste("results/", seq, "_site_methylation.txt", sep="")
      df <- read.table(input,header=T)
      x <- rbind(x,cbind(rbind(boxplot.stats(df$mc_level)$stats),rbind(boxplot.stats(df$mc_level)$conf)))
      file.remove(input)
  }
  row.names(x) <- context
  colnames(x) <- c('ymin','lower','middle','upper','ymax','notchlower','notchupper')
  write.table(x,filename,quote=F,sep='\t')
}

if(file.exists(filename)){
  print("Making per-site methylation plot")
  x <- read.table(filename,header=T,sep="\t")
  x$order <- c(1:19)
  x$color <- c(1,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3)
  plot <- ggplot(x, aes(reorder(row.names(x),order), fill=factor(color))) +
          geom_boxplot(aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,
            notchlower=notchlower,notchupper=notchupper),stat="identity",notch=T) +
            theme(panel.background=element_blank(), panel.grid=element_blank(),
            axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
            axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
            legend.position="none", axis.line=element_line(color="black")) +
            ylab("Methylation level") + xlab("Context") +
            scale_y_continuous(limits=c(0,1.05), expand=c(0,0),
            breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
            scale_fill_manual("",values = c("dodgerblue4","olivedrab","hotpink4"))
  filename=paste("figures_tables/Fvesca_", args[1], "_per_site_boxplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=8, useDingbats=F)
}

#genome distributions
if(file.exists("results/genome_windows_data.tsv")){
  print("Making genome distribution plots")
  df <- read.table("results/genome_windows_data.tsv",header=T,sep="\t")

  #plot mCG gene correlations
  plot <- plot_correlations(df,df$mCG,"dodgerblue4") +
          scale_y_continuous(limits=c(0,1), expand=c(0,0),
          breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%"))
  filename=paste("figures_tables/Fvesca_", args[1], "_gene_mCG_correlation.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(plot)

  #plot mCHG gene correlations
  plot <- plot_correlations(df,df$mCHG,"olivedrab") +
          scale_y_continuous(limits=c(0,1), expand=c(0,0),
          breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%"))
  filename=paste("figures_tables/Fvesca_", args[1], "_gene_mCHG_correlation.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(plot)

  #plot mCHH gene correlations
  plot <- plot_correlations(df,df$mCHH,"hotpink4") +
          scale_y_continuous(limits=c(0,0.1), expand=c(0,0),
          breaks=c(0.025,0.05,0.075,0.1),labels=c("2.5%","5.0%","7.5%","10%"))
  filename=paste("figures_tables/Fvesca_", args[1], "_gene_mCHH_correlation.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(plot)

  #plot chr metaplots
  print("Making chromosome plots")
  library(stringi)
  df$chr <- stri_replace_last(df$window,'',regex="_[0-9]+.*?")
  df$window <- as.numeric(stri_replace(df$window,'',regex=".*_"))
  chr <- unique(df$chr)
  for(i in chr){
     tmp <- df[df$chr == i,]
     if (max(tmp$window) > 10){
       break_points <- c(2,4,6,8,10)*max(tmp$window)/10
       label_points <- c(2,4,6,8,10)*max(tmp$window)/20
       plot = ggplot(tmp, aes(x=window, group=1)) +
              geom_line(aes(y=mCG), color="dodgerblue4", size=0.8) +
              geom_line(aes(y=mCHG), color="olivedrab", size=0.8) +
              geom_line(aes(y=mCHH), color="hotpink4", size=0.8) +
              theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
              axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
              legend.position="none", axis.line=element_line(color="black")) +
              ylab("Percent methylation") + xlab( "Kbps" ) +
              scale_y_continuous(limits=c(0,1), expand=c(0,0),
              breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
              scale_x_continuous(expand=c(0,0), breaks=break_points, labels=label_points)
       filename=paste("figures_tables/Fvesca_", args[1], "_Chr", i,"_metaplot.pdf", sep="")
       ggsave(filename=filename, plot, height=4, width=8, useDingbats=F)
       rm(tmp,plot)
     }
  }
  rm(df)
}

#plot repeat metaplots
if(file.exists("results/repeat_metaplot.tsv")){
  print("Making repeat metaplot")
  df <- read.table("results/repeat_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","5'","3'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/Fvesca_", args[1], "_repeat_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}

#plot gene metaplots
if(file.exists("results/gene_metaplot.tsv")){
  print("Making gene metaplot")
  df <- read.table("results/gene_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/Fvesca_", args[1], "_gene_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}
