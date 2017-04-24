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


PlotYearVsPSG <- function(infilename, graphname, ColPalette, textsize, h, w){
  YearVsAA     <- read_tsv(infilename) %>%
		    filter(Year>=2007) %>%
		    mutate(AA=factor(AA,levels=c('I','M','P','L')))
  Egg_YearVsAA <- YearVsAA %>% 
		    filter(Passage=='Egg')
  Ori_YearVsAA <- YearVsAA %>% 
		    filter(Passage=='Ori')
  print (paste('After 2007: ', sum(Egg_YearVsAA$Count),' egg-passaged sequences',sep=''))
  print (paste('After 2007: ', sum(Ori_YearVsAA$Count),' unpassaged sequences',sep=''))
  p_egg <- ggplot(Egg_YearVsAA,aes(x=Year,y=Count,fill=AA)) + 
	     geom_bar(stat="identity", position="fill") +
	     scale_fill_manual(values=ColPalette) +
	     theme(axis.title=element_text(size=textsize,face="bold"),
		   axis.text=element_text(size=textsize,face="bold"),
		   axis.text.x=element_blank(),
		   axis.title.x=element_blank(),
		   legend.title=element_blank(),
		   legend.key.size = unit(0.5, 'lines'),
		   legend.text=element_text(size=textsize,face="bold"),
		   legend.position='top') +
	     xlab(bquote(bold(Year))) +
	     ylab(bquote(bold(Percentage))) +
	     scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,50,100))
  p_ori <- ggplot(Ori_YearVsAA,aes(x=Year,y=Count,fill=AA)) +
	     geom_bar(stat="identity", position="fill") +
	     scale_fill_manual(values='#8dc68d') +
	     theme(axis.title=element_text(size=textsize,face="bold"),
		   axis.text=element_text(size=textsize,face="bold"),
		   axis.text.x=element_text(size=textsize,face="bold",angle=90,vjust=0.5),
		   legend.title=element_blank(),
		   legend.text=element_text(size=textsize,face="bold"),
		   legend.position='none') +
	     xlab(bquote(bold(Year))) +
	     ylab(bquote(bold(Percentage))) +
	     scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,50,100))
  p_Comb <- grid.arrange(rbind(ggplotGrob(p_egg),ggplotGrob(p_ori),size="last"))
  ggsave(filename=graphname,p_Comb,height=h,width=w)
  }

PlotYearVsPSG('result/HumanH3N2_Pos194YearVsPSG.tsv','graph/H3N2_YearVsAA_resi194.png',
              c('#e2a3a3','#d8ac52','#9e9ed3','#8dc68d'),7,1.8,1.4)
PlotYearVsPSG('result/pdmH1N1_Pos194YearVsPSG.tsv','graph/pdmH1N1_YearVsAA_resi194.png',
              c('#e2a3a3','#8dc68d'),11,3.6,2.8)
