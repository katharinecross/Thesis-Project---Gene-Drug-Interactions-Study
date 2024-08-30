#install.packages("devtools")
#devtools::install_github("rlbarter/superheat")

library(superheat)
data1 = read.csv("C:/Users/student/OneDrive - Bryant University/Desktop/Honors Thesis/HonorThesisKatie/Cluster1/subGeneMtx01.csv", header=FALSE)
data2 = read.csv("C:/Users/student/OneDrive - Bryant University/Desktop/Honors Thesis/HonorThesisKatie/Cluster2/subGeneMtx01.csv", header=FALSE)


levels1 = c('Vorinostat','DMSO')


cluster1 =  factor(c(rep('Vorinostat',(dim(data1)[1]/2)), rep('DMSO', (dim(data1)[1]/2))),levels=levels1)
cluster2 =  factor(c(rep('Vorinostat',(dim(data2)[1]/2)), rep('DMSO', (dim(data2)[1]/2))),levels=levels1)

library("RColorBrewer")
#install.packages("gplots")
library(gplots)
heatmap.2(data.matrix(data1),col=brewer.pal(9,"RdBu"), dendrogram="none",
          Rowv=NULL,Cov=F, trace="none",
          breaks = c(-56, -30,-20, -10, -5, 0, 5, 10, 30, 56))
varNames1 = read.csv("C:/Users/student/OneDrive - Bryant University/Desktop/Honors Thesis/HonorThesisKatie/Cluster1/allSelectedGenes01.csv", header=TRUE)
geneNames1 = varNames['pr_gene_symbol']
names(data1) = unlist(geneNames1)

Ykmean1 = read.csv("C:/Users/student/OneDrive - Bryant University/Desktop/Honors Thesis/HonorThesisKatie/KnockoffAverage1.csv", header=FALSE)
Ykselected1 = Ykmean1[, varNames1[,1]]
names(data1) = unlist(geneNames1)
png("C:/Users/student/OneDrive - Bryant University/Desktop/Honors Thesis/Plots/superheatFDR01knockoff.png", height = 900, width = 900)
superheat(Ykselected1,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster1,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,
          bottom.label.size = 0.12,#legend.vspace = 0,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(10,"RdBu"),
          heat.pal.values = c(0, 0.1,0.2, 0.25, 0.3, 0.38, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
          #title='Heatmap of Robust Z-scores for Selected Gene Knockoffs',
          title.size=5,
          # legend.num.ticks = 6
          #heat.pal.values = c( 0,0.45, 0.55,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1)
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25, 40, 60)
)
dev.off()

#-----
png("C:/Users/student/OneDrive - Bryant University/Desktop/Honors Thesis/Plots/superheatFDR01Full.png", height = 900, width = 900)
superheat(data1,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster1,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,# legend.vspace = 0,
          bottom.label.size = 0.12,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(10,"RdBu"),
          heat.pal.values = c( 0,0.2, 0.35,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1),
          #title='Heatmap of Robust Z-scores for Selected Genes',
          title.size=5,
          #legend.num.ticks = 6
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25)
)
dev.off()
#----
varNames2 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/trichostatin-a/allSelectedGenes01.csv", header=TRUE)
geneNames2 = varNames2['pr_gene_symbol']
names(data2) = unlist(geneNames2)

Ykmean2 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/trichostatin-a/KnockoffAverage.csv", header=FALSE)
Ykselected2 = Ykmean2[, varNames2[,1]]
names(Ykselected2) = unlist(geneNames2)

png("/Users/tingtingzhao/Documents/data/rubicinAndOther/trichostatin-a/superheatFDR02Full.png", height = 900, width = 900)
superheat(data2,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster2,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,# legend.vspace = 0,
          bottom.label.size = 0.12,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(10,"RdBu"),
          heat.pal.values = c( 0,0.2, 0.35,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1),
          #title='Heatmap of Robust Z-scores for Selected Genes',
          title.size=5,
          #legend.num.ticks = 6
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25)
)
dev.off()

png("/Users/tingtingzhao/Documents/data/rubicinAndOther/trichostatin-a/superheatFDR01Knockoff.png", height = 900, width = 900)
superheat(Ykselected2,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster2,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,
          bottom.label.size = 0.12,#legend.vspace = 0,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(10,"RdBu"),
          heat.pal.values = c(0, 0.1,0.2, 0.25, 0.3, 0.38, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
          #title='Heatmap of Robust Z-scores for Selected Gene Knockoffs',
          title.size=5,
          # legend.num.ticks = 6
          #heat.pal.values = c( 0,0.45, 0.55,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1)
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25, 40, 60)
)
dev.off()
#-------------------
varNames3 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/wortmannin/allSelectedGenes01.csv", header=TRUE)
geneNames3 = varNames3['pr_gene_symbol']
names(data3) = unlist(geneNames3)

Ykmean3 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/wortmannin/KnockoffAverage.csv", header=FALSE)
varNames3[which(varNames3[,1]==0), 1]=1
Ykselected3 = Ykmean3[, (varNames3[,1])]
names(Ykselected3) = unlist(geneNames3)

quantile(c(as.matrix(data3)),1:9/10)
#----
png("/Users/tingtingzhao/Documents/data/rubicinAndOther/wortmannin/superheatFDR02Full.png", height = 900, width = 900)
superheat(data3,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster3,left.label.text.alignment = 'center',
          left.label.text.size = 5,
          left.label.size = 0.15,# legend.vspace = 0,
          bottom.label.size = 0.12,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(9,"RdBu"),
          # heat.pal.values = c( 0, 0.1,0.15, 0.2, 0.35, 0.4, 0.5, 0.65, 0.8, 0.9, 1),
          heat.pal.values = c(0.2, 0.4,0.45, 0.48, 0.5, 0.52, 0.55, 0.6, 0.8),
          #title='Heatmap of Robust Z-scores for Selected Genes',
          title.size=5,
          #legend.num.ticks = 6
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25)
)
dev.off()
#----
#hist(c(as.matrix(Ykselected3)))
vis_twopart<-function(data){
  n=nrow(data)
  x1=c(as.matrix(data[1:n/2,]))
  x2=c(as.matrix(data[(n/2+1):n,]))
  q1=quantile(x1,1:9/10)
  q2=quantile(x2,1:9/10)
  rbind(q1,q2)
}
vis_twopart(Ykselected3)
vis_twopart(data3)

show_col(brewer.pal(11,'RdBu'))
#----
png("/Users/tingtingzhao/Documents/data/rubicinAndOther/wortmannin/superheatFDR01Knockoff.png", height = 900, width = 900)
superheat(Ykselected3,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster3,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,
          bottom.label.size = 0.12,#legend.vspace = 0,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(9,"RdBu"),
          heat.pal.values = c(0, 0.1,0.2, 0.3, 0.5,0.7, 0.8, 0.9, 1),
          #title='Heatmap of Robust Z-scores for Selected Gene Knockoffs',
          title.size=5,
          # legend.num.ticks = 6
          #heat.pal.values = c( 0,0.45, 0.55,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1)
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25, 40, 60)
)
dev.off()
#---------------------------------------
varNames4 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/geldanamycin/allSelectedGenes01.csv", header=TRUE)
geneNames4 = varNames4['pr_gene_symbol']
names(data4) = unlist(geneNames4)

Ykmean4 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/geldanamycin/KnockoffAverage.csv", header=FALSE)
Ykselected4 = Ykmean4[, (varNames4[,1])]
names(Ykselected4) = unlist(geneNames4)

#----
show_col(brewer.pal(9,'RdBu'))
png("/Users/tingtingzhao/Documents/data/rubicinAndOther/geldanamycin/superheatFDR02Full.png", height = 900, width = 900)
superheat(data4,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster4,left.label.text.alignment = 'center',
          left.label.text.size = 5,
          left.label.size = 0.15,# legend.vspace = 0,
          bottom.label.size = 0.12,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(10,"RdBu"),
          heat.pal.values = c(0, 0.1,0.2, 0.3,0.4,0.5,0.55, 0.6, 0.75, 1),
          #title='Heatmap of Robust Z-scores for Selected Genes',
          title.size=5
          #legend.num.ticks = 6
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25)
)
dev.off()

#----
png("/Users/tingtingzhao/Documents/data/rubicinAndOther/geldanamycin/superheatFDR01Knockoff.png", height = 900, width = 900)
superheat(Ykselected4,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster4,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,
          bottom.label.size = 0.12,#legend.vspace = 0,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(9,"RdBu"),
          heat.pal.values = c(0, 0.1,0.2, 0.3, 0.5,0.6, 0.7, 0.8, 1),
          #title='Heatmap of Robust Z-scores for Selected Gene Knockoffs',
          title.size=5,
          # legend.num.ticks = 6
          #heat.pal.values = c( 0,0.45, 0.55,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1)
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25, 40, 60)
)
dev.off()
#---------------------------------------
varNames5 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/Sirolimus/allSelectedGenes01.csv", header=TRUE)
geneNames5 = varNames5['pr_gene_symbol']
names(data5) = unlist(geneNames5)

Ykmean5 = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/Sirolimus/KnockoffAverage.csv", header=FALSE)
Ykselected5 = Ykmean5[, (varNames5[,1])]
names(Ykselected5) = unlist(geneNames5)

#----
png("/Users/tingtingzhao/Documents/data/rubicinAndOther/Sirolimus/superheatFDR02Full.png", height = 900, width = 900)
superheat(data5,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster5,left.label.text.alignment = 'center',
          left.label.text.size = 5,
          left.label.size = 0.15,# legend.vspace = 0,
          bottom.label.size = 0.12,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(9,"RdBu"),
          heat.pal.values = c(0, 0.4, 0.45,0.48,0.5,0.52, 0.55, 0.6, 1),
          #title='Heatmap of Robust Z-scores for Selected Genes',
          title.size=5,
          #legend.num.ticks = 6
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25)
)
dev.off()

#----
png("/Users/tingtingzhao/Documents/data/rubicinAndOther/Sirolimus/superheatFDR01Knockoff.png", height = 900, width = 900)
superheat(Ykselected5,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.cols = FALSE,
          #bottom.label='none',
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3.5,
          # plot miles per gallon next to the rows
          # plot correlation with mpg above columns
          # increase size of left labels
          # left.label.size = 0.45
          membership.rows = cluster5,left.label.text.alignment = 'center',
          left.label.text.size = 6,
          left.label.size = 0.15,
          bottom.label.size = 0.12,#legend.vspace = 0,
          grid.hline.col = "black",
          grid.vline = FALSE,
          # heat.pal = c("#b35806", "",  "white", "#542788"),
          heat.pal = brewer.pal(9,"RdBu"),
          heat.pal.values = c(0, 0.1,0.2, 0.3, 0.45,0.55, 0.65, 0.75, 1),
          #title='Heatmap of Robust Z-scores for Selected Gene Knockoffs',
          title.size=5,
          # legend.num.ticks = 6
          #heat.pal.values = c( 0,0.45, 0.55,0.6, 0.65, 0.7, 0.75,0.8, 0.9, 1)
          #heat.pal = c('blue', 'pink', 'black', 'pink'),
          #heat.pal.values = c(-60,-25, 0, 25, 40, 60)
)
dev.off()
#--------
df = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/vorinostat/allSelectedGenes01.csv", header=TRUE)
df = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/trichostatin-a/allSelectedGenes01.csv", header=TRUE)
df = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/wortmannin/allSelectedGenes01.csv", header=TRUE)
df = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/geldanamycin/allSelectedGenes01.csv", header=TRUE)
df = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/sirolimus/allSelectedGenes01.csv", header=TRUE)

genename= df$pr_gene_symbol[1:12]
df_sub = data5[, 1:12]
df_sub$Method = cluster5
# cbPalette = c( "#FFCCFF", "#66CCFF") #"#CC79A7",,
cbPalette = c("#C5E7F5", "#F9AB27")
library(tidyverse)
df_long = df_sub%>%pivot_longer(cols=1:12)#%>%mutate(name=factor(name,level=paste0('V',1:12),labels=genename))
library(ggplot2)
library(wesanderson)
mytheme = theme_bw()+ theme(strip.text.x = element_text(size = 12))+
  theme(axis.title= element_text(size=14))+
  theme(axis.text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(df_long, aes(x=Method,y=value))+geom_boxplot()+facet_wrap(~name)
ggplot(df_long, aes(x=value,fill=Method))+geom_density(alpha=0.8)+facet_wrap(~name)+
  xlim(-10, 10)+
  ylim(0, 0.65)+
  geom_vline(aes(xintercept = 0), df_long, color='blue')+
  mytheme+
  theme(strip.background = element_blank())+
  #theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Density of robust z-scores")+
  xlab("Robust z-score values truncated at (-10, 10)")+
  theme(legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal')+
  #scale_fill_manual(values=wes_palette(n=2, name="Moonrise3"))+
  #scale_fill_brewer(palette = "Pastel1")
  scale_fill_manual(values=cbPalette)+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  theme(legend.title=element_text(size=12))+
  theme(legend.text=element_text(size=10))
ggsave("/Users/tingtingzhao/Documents/data/rubicinAndOther/sirolimus/TrucatedDensityFDR02Full.pdf", width = 16, height= 12, units = "cm")

## compute pairwise correlation
corr =  read.csv("/Users/tingtingzhao/Documents/data/rubicin/result/OverallCorr.csv", header=FALSE)
corr_df = c(unlist(c(corr)))
corr_df = data.frame(corr=corr_df)
library(scales)
library(dplyr)
summ <- corr_df %>% 
  summarize(min = min(corr), max = max(corr), 
            mean = mean(corr), median = median(corr),
            sd = sd(corr)) %>% mutate(lab = paste0("min = ", round(min,3), "\nmax = ", round(max,3),
                                                   "\nmean = ", round(mean,3),  "\nmedian = ", round(median,3),
                                                   "\nsd = ",round(sd,3)))%>% select(lab)

ggplot()+
  geom_histogram(data=corr_df, aes(x=corr, y = stat(count) / sum(count)), color='orange',fill="skyblue", alpha=0.5)+
  mytheme+geom_text(data = summ, aes(label = lab), x=-Inf, y=Inf, hjust=-0.05, vjust=1.2,size=6)+
  xlab('Pairwise correlation')+
  ylab('Relative frequency') +
  ggsave("/Users/tingtingzhao/Documents/data/rubicin/result/heatmap/PairwiseCorrHist.pdf", width = 16, height= 12, units = "cm")
#####################################################
## create error bar plots
df = read.csv("/Users/tingtingzhao/Documents/data/rubicinAndOther/sirolimus/allSelectedGenes01.csv", header=TRUE)
errorbar_df = df[1:12,]
errorbar_df$pr_gene_sym=factor(errorbar_df$pr_gene_symbol,levels = errorbar_df$pr_gene_symbol[12:1])

ggplot(errorbar_df) +
  geom_bar( aes(x=pr_gene_sym, y=mean), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar( aes(x=pr_gene_symbol, ymin=mean-std, ymax=mean+std), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  coord_flip()+theme_bw()+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.title= element_text(size=14))+
  theme(axis.text = element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab( "Mean and standard error of coefficients")+
  xlab("Top 12 selected genes")
ggsave("/Users/tingtingzhao/Documents/data/rubicinAndOther/sirolimus/Coef.pdf", width = 16, height= 12, units="cm")




