#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
require(cowplot)
pdf(NULL)

ProVsPSG     <- read_tsv('result/HumanH3N2_PSG.tsv')
textsize <- 14
dodge <- position_dodge(width=0.9)
p_Pro <- ggplot(ProVsPSG,aes(x=`Passage Number`,y=`Proportion of Pro`,color='#e2a3a3')) +
           geom_point(size=2,color='#9e9ed3') +
           geom_errorbar(aes(ymax=`Proportion of Pro`+SE, ymin=`Proportion of Pro`-SE), width=0.25, color='#9e9ed3',size=0.8) +
           geom_line(color='#9e9ed3',size=1.3) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 axis.text.x=element_text(size=textsize,face="bold",vjust=0.5),
                 legend.title=element_blank(),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none') +
           xlab(bquote(bold('Passage Number'))) +
           ylab(bquote(bold('% of Pro')))
ggsave(filename='graph/HumanH3N2_ProVsPSG.png',p_Pro,height=3.5,width=5)
