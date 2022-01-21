setwd("/Volumes/My_disk/Wang_hb/Liu_my/H3K4me3_2nd/")
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicAlignments)


########################     Menin 4-14-2021  #########################
setwd("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata")
WT_Menin = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/Menin_normalized_2_peaks.txt", sep = "\t", header = TRUE)
WT_input = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/F3_FKDL202605429-1aAligned_peaks.txt", sep = "\t", header = TRUE)

WT_input_gr=GRanges(WT_input)
WT_Menin_gr=GRanges(WT_Menin)

txdb=makeTxDbFromGFF("/Users/wenbodeng/E/Index/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf", 
                     format=c("auto"), 
                     organism="Mus musculus",circ_seqs=character())

anno_WT_Menin =annotatePeak(WT_Menin_gr,tssRegion = c(-3000,3000),
                            TxDb=txdb,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

anno_WT_input =annotatePeak(WT_input_gr,tssRegion = c(-3000,3000),
                            TxDb=txdb,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

anno_WT_Menin_gr = anno_WT_Menin@anno
anno_WT_input_gr = anno_WT_input@anno

m=findOverlaps(anno_WT_Menin_gr, anno_WT_input_gr)
anno_WT_Menin_gr_specific = anno_WT_Menin_gr[-queryHits(m),]


Ac_anno_con_df=data.frame(anno_WT_Menin_gr_specific)
Ac_anno_con_df$annotation=gsub("Intron (.*)", "Intron", Ac_anno_con_df$annotation)
Ac_anno_con_df$annotation=gsub("Promoter (.*)", "Promoter", Ac_anno_con_df$annotation)
Ac_anno_con_df$annotation=gsub("Exon (.*)", "Exon", Ac_anno_con_df$annotation)
Ac_anno_con_df$annotation=gsub("Downstream (.*)", "Downstream", Ac_anno_con_df$annotation)
table(Ac_anno_con_df$annotation)


library(ggplot2)
display.brewer.all()

library("RColorBrewer")

my_col = rev(colorRampPalette(brewer.pal(name="Dark2", n = 8))(7))  ## Good
my_col = c("#4DAF4A" ,"#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02" ,"#A6761D", "#4DAF4A")
my_col = c("#4DAF4A" ,"#D95F02" ,"#984EA3" ,"#E7298A" ,"#66A61E", "#E6AB02", "#A6761D" ,"#666666")

pie <- ggplot(Ac_anno_con_df, aes(x = factor(1),fill=factor(annotation))) +geom_bar(width=1,color="white")
pie + coord_polar("y",direction = 1)+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +
  scale_fill_manual(values =my_col)+ 
  theme(panel.grid=element_blank())  + 
  theme(axis.text = element_blank(),axis.ticks = element_blank(), panel.grid  = element_blank()
)





###########  count ###########
library(GenomicAlignments)

me3_WT_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/WT_FKDL190744256-1aAligned.sortedByCoord.out.bam")
me3_KO_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/KO_FKDL190744258-1aAligned.sortedByCoord.out.bam")
menin_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/F1_FKDL202605427-1aAligned.sortedByCoord.out.bam")

anno_WT_Menin_gr_specific$annotation=gsub("Intron (.*)", "Intron", anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific$annotation=gsub("Promoter (.*)", "Promoter", anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific$annotation=gsub("Exon (.*)", "Exon", anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific$annotation=gsub("Downstream (.*)", "Downstream", anno_WT_Menin_gr_specific$annotation)
table(anno_WT_Menin_gr_specific$annotation)

anno_WT_Menin_gr_specific = anno_WT_Menin_gr_specific[anno_WT_Menin_gr_specific$annotation %in% "Promoter",]

counts_men_WT <- countOverlaps(anno_WT_Menin_gr_specific, menin_aligns)
width_men_WT = width(anno_WT_Menin_gr_specific)  # length for each peak
rpkm_men_WT = counts_men_WT / ( width_men_WT/1000 * sum(counts_men_WT)/ 1e+06 )   # RPKM
anno_WT_Menin_gr_specific$men_reads = counts_men_WT
anno_WT_Menin_gr_specific$men_rpkm = rpkm_men_WT

counts_me3_WT <- countOverlaps(anno_WT_Menin_gr_specific, me3_WT_aligns)
width_me3_WT = width(anno_WT_Menin_gr_specific)  # length for each peak
rpkm_me3_WT = counts_me3_WT / ( width_me3_WT/1000 * sum(counts_me3_WT)/ 1e+06 )   # RPKM
anno_WT_Menin_gr_specific$WT_me3_reads = counts_me3_WT
anno_WT_Menin_gr_specific$WT_me3_rpkm = rpkm_me3_WT

counts_me3_KO <- countOverlaps(anno_WT_Menin_gr_specific, me3_KO_aligns)
width_me3_KO = width(anno_WT_Menin_gr_specific)  # length for each peak
rpkm_me3_KO = counts_me3_KO / ( width_me3_KO/1000 * sum(counts_me3_KO)/ 1e+06 )   # RPKM
anno_WT_Menin_gr_specific$KO_me3_reads = counts_me3_KO
anno_WT_Menin_gr_specific$KO_me3_rpkm = rpkm_me3_KO

anno_WT_Menin_df = as.data.frame(anno_WT_Menin_gr_specific)


p = ggplot(anno_WT_Menin_df,aes(log2(men_rpkm+1),log2(WT_me3_rpkm+1)))+ geom_point()+ylim(0,20)+xlim(0,20) 
p+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank())  

p <- ggplot(anno_WT_Menin_df, aes(log(men_reads,2), log(WT_me3_reads,2)))+ geom_point()+xlim(-1,20)+ylim(-1,20) +    # rpkm
  stat_smooth(method = lm)
p


a = anno_WT_Menin_df[,c(24:25)]
a$flag = "WT"
colnames(a) = c("me3_reads", "me3_rpkm", "flag")
b = anno_WT_Menin_df[,c(26:27)]
b$flag = "KO"
colnames(b) = c("me3_reads", "me3_rpkm", "flag")
c = rbind(a,b)
c$me3_reads_log = log2(c$me3_reads+1)
c$me3_rpkm_log = log2(c$me3_rpkm+1)

library(ggpubr)
my_comparisons <- list( c("WT", "KO")) 

p = ggviolin(c, x = "flag", y = "me3_reads_log", 
             fill = "flag",
             color =  "flag",
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#FC4E01","#FC4EA7"),
             add = "boxplot", add.params = list(fill = NULL,color = "black")
) 

p+stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 20)  






###########   smoothScatter of H3K4me3  8-2-2021   ###########

#WT_Menin_1 = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/Menin_normalized_peaks.txt", sep = "\t", header = TRUE)
WT_Menin = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/Menin_normalized_2_peaks.txt", sep = "\t", header = TRUE)

WT_input = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/F3_FKDL202605429-1aAligned_peaks.txt", sep = "\t", header = TRUE)
#WT_H3K27me3 = read.table("/Volumes/My_disk/Wang_hb/Liu_my/H3K4me3_2nd/LMY12A1_FKDL190766338_1a_without_summit_peaks.txt", sep = "\t", header = TRUE)
WT_input_gr=GRanges(WT_input)
WT_Menin_gr=GRanges(WT_Menin)
txdb=makeTxDbFromGFF("/Users/wenbodeng/E/Index/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf", 
                     format=c("auto"), 
                     organism="Mus musculus",circ_seqs=character())

anno_WT_Menin =annotatePeak(WT_Menin_gr,tssRegion = c(-3000,3000),
                            TxDb=txdb,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

anno_WT_input =annotatePeak(WT_input_gr,tssRegion = c(-3000,3000),
                            TxDb=txdb,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

#anno_WT_H3K27me3_gr = anno_WT_H3K27me3@anno
anno_WT_Menin_gr = anno_WT_Menin@anno
anno_WT_input_gr = anno_WT_input@anno

m=findOverlaps(anno_WT_Menin_gr, anno_WT_input_gr)
anno_WT_Menin_gr_specific = anno_WT_Menin_gr[-queryHits(m),]

anno_WT_Menin_gr_specific$annotation=gsub("Intron (.*)", "Intron", anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific$annotation=gsub("Promoter (.*)", "Promoter", anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific$annotation=gsub("Exon (.*)", "Exon", anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific$annotation=gsub("Downstream (.*)", "Downstream", anno_WT_Menin_gr_specific$annotation)
table(anno_WT_Menin_gr_specific$annotation)
anno_WT_Menin_gr_specific = anno_WT_Menin_gr_specific[anno_WT_Menin_gr_specific$annotation %in% "Promoter",]


library(GenomicAlignments)
menin_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/F1_FKDL202605427-1aAligned.sortedByCoord.out.bam")
counts_men_WT <- countOverlaps(anno_WT_Menin_gr_specific, menin_aligns)
width_men_WT = width(anno_WT_Menin_gr_specific)  # length for each peak
rpkm_men_WT = counts_men_WT / ( width_men_WT/1000 * sum(counts_men_WT)/ 1e+06 )   # RPKM
anno_WT_Menin_gr_specific$men_reads = counts_men_WT
anno_WT_Menin_gr_specific$men_rpkm = rpkm_men_WT


H3K4me3_WT = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/WT_FKDL190744256-1a_peaks.xls", sep = "\t", header = TRUE)
H3K4me3_WT_gr=GRanges(H3K4me3_WT)

anno_H3K4me3_WT =annotatePeak(H3K4me3_WT_gr,tssRegion = c(-3000,3000),
                            TxDb=txdb,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))

anno_H3K4me3_WT_gr = anno_H3K4me3_WT@anno

m=findOverlaps(anno_H3K4me3_WT_gr, anno_WT_input_gr)
anno_H3K4me3_WT_gr_specific = anno_H3K4me3_WT_gr[-queryHits(m),]

anno_H3K4me3_WT_gr_specific$annotation=gsub("Intron (.*)", "Intron", anno_H3K4me3_WT_gr_specific$annotation)
anno_H3K4me3_WT_gr_specific$annotation=gsub("Promoter (.*)", "Promoter", anno_H3K4me3_WT_gr_specific$annotation)
anno_H3K4me3_WT_gr_specific$annotation=gsub("Exon (.*)", "Exon", anno_H3K4me3_WT_gr_specific$annotation)
anno_H3K4me3_WT_gr_specific$annotation=gsub("Downstream (.*)", "Downstream", anno_H3K4me3_WT_gr_specific$annotation)
table(anno_H3K4me3_WT_gr_specific$annotation)

anno_H3K4me3_WT_gr_specific = anno_H3K4me3_WT_gr_specific[anno_H3K4me3_WT_gr_specific$annotation %in% "Promoter",]


library(GenomicAlignments)
me3_WT_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/WT_FKDL190744256-1aAligned.sortedByCoord.out.bam")
counts_me3_WT <- countOverlaps(anno_H3K4me3_WT_gr_specific, me3_WT_aligns)
width_me3_WT = width(anno_H3K4me3_WT_gr_specific)  # length for each peak
rpkm_me3_WT = counts_me3_WT / ( width_me3_WT/1000 * sum(counts_me3_WT)/ 1e+06 )   # RPKM
anno_H3K4me3_WT_gr_specific$WT_me3_reads = counts_me3_WT
anno_H3K4me3_WT_gr_specific$WT_me3_rpkm = rpkm_me3_WT



mm=findOverlaps(anno_WT_Menin_gr_specific, anno_H3K4me3_WT_gr_specific)
anno_WT_Menin_gr_specific_common = anno_WT_Menin_gr_specific[queryHits(mm),]
anno_H3K4me3_WT_gr_specific_common = anno_H3K4me3_WT_gr_specific[subjectHits(mm),]
anno_WT_Menin_gr_specific_common_df = data.frame(anno_WT_Menin_gr_specific_common)
anno_H3K4me3_WT_gr_specific_common_df = data.frame(anno_H3K4me3_WT_gr_specific_common)
common_peaks = cbind (anno_WT_Menin_gr_specific_common_df, anno_H3K4me3_WT_gr_specific_common_df)

row.names(common_peaks) = paste(common_peaks$seqnames,common_peaks$start,common_peaks$end,common_peaks$WT_me3_rpkm)

smoothScatter(log(common_peaks$men_reads,2),log(common_peaks$WT_me3_reads,2), 
              nbin = 512 ,
              transformation = function(x) x^.2,
              xlim = c(0,12),
              ylim = c(0,15)) #####  density plot


cor(log(common_peaks$men_reads,2),log(common_peaks$WT_me3_reads,2),method = c("pearson"))   ##  0.273329

aa = cbind(common_peaks$men_reads,common_peaks$WT_me3_reads)
colnames(aa) = c("men_reads","WT_me3_reads")
row.names(aa) = row.names(common_peaks)
aa = as.data.frame(aa)
p <- ggplot(aa, aes(log(men_reads,2), log(WT_me3_reads,2)))+ geom_point()+xlim(2,13)+ylim(5,15) +    # rpkm
  stat_smooth(method = lm)
p






########### bed file for ngsplotr   8-2-2021   ##############
wt_peak  = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/WT_FKDL190744256-1a_peaks.xls", sep = "\t", header = T)
wt_peak = GRanges(wt_peak)
library(GenomicAlignments)
me3_WT_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/WT_FKDL190744256-1aAligned.sortedByCoord.out.bam")
counts_me3_WT <- countOverlaps(wt_peak, me3_WT_aligns)
width_me3_WT = width(wt_peak)  # length for each peak
rpkm_me3_WT = counts_me3_WT / ( width_me3_WT/1000 * sum(counts_me3_WT)/ 1e+06 )   # RPKM
wt_peak$WT_me3_reads = counts_me3_WT
wt_peak$WT_me3_rpkm = rpkm_me3_WT


KO_peak  = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/KO_FKDL190744258-1a_peaks.xls", sep = "\t", header = T)
KO_peak = GRanges(KO_peak)
library(GenomicAlignments)
me3_KO_aligns <- readGAlignments("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/H3K4me3/KO_FKDL190744258-1aAligned.sortedByCoord.out.bam")
counts_me3_WT <- countOverlaps(KO_peak, me3_KO_aligns)
width_me3_WT = width(KO_peak)  # length for each peak
rpkm_me3_WT = counts_me3_WT / ( width_me3_WT/1000 * sum(counts_me3_WT)/ 1e+06 )   # RPKM
KO_peak$KO_me3_reads = counts_me3_WT
KO_peak$KO_me3_rpkm = rpkm_me3_WT


mm=findOverlaps(wt_peak, KO_peak)
wt_peak_common = wt_peak[queryHits(mm),]
KO_peak_common = KO_peak[subjectHits(mm),]
wt_peak_common_df = data.frame(wt_peak_common)
KO_peak_common_df = data.frame(KO_peak_common)
H3K4me3_common_peaks = cbind (wt_peak_common_df, KO_peak_common_df)

H3K4me3_common_peaks$fold = H3K4me3_common_peaks$KO_me3_reads / H3K4me3_common_peaks$WT_me3_reads

H3K4me3_common_peaks_2 = H3K4me3_common_peaks[order(H3K4me3_common_peaks$fold,decreasing = T),]
H3K4me3_common_peaks_2$flag = cut(H3K4me3_common_peaks_2$fold, breaks = c(-Inf,0.75,1.5,Inf), labels = c ("-1","0","1"))

H3K4me3_common_peaks_s1 = H3K4me3_common_peaks_2[H3K4me3_common_peaks_2$flag == 1,]
H3K4me3_common_peaks_s1 = H3K4me3_common_peaks_s1[order(H3K4me3_common_peaks_s1$WT_me3_reads, decreasing = T),]

H3K4me3_common_peaks_s2 = H3K4me3_common_peaks_2[H3K4me3_common_peaks_2$flag == 0,]
H3K4me3_common_peaks_s2 = H3K4me3_common_peaks_s2[order(H3K4me3_common_peaks_s2$WT_me3_reads, decreasing = T),]

H3K4me3_common_peaks_s3 = H3K4me3_common_peaks_2[H3K4me3_common_peaks_2$flag == -1,]
H3K4me3_common_peaks_s3 = H3K4me3_common_peaks_s3[order(H3K4me3_common_peaks_s3$WT_me3_reads, decreasing = T),]

H3K4me3_common_peaks_new = rbind(H3K4me3_common_peaks_s1 ,H3K4me3_common_peaks_s2, H3K4me3_common_peaks_s3)

write.table(H3K4me3_common_peaks_new[,1:3], "bed_for_ngsplotr.txt",row.names = F, col.names = F, sep = "\t", quote = F)



KO_peak = GRanges(KO_peak)

H3K4me3_common_peaks_s3_anno =annotatePeak(GRanges(H3K4me3_common_peaks_s3[,1:3]),tssRegion = c(-3000,3000),
                            TxDb=txdb,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
H3K4me3_common_peaks_s3_anno = as.data.frame(H3K4me3_common_peaks_s3_anno)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Intron (.*)", "Intron", H3K4me3_common_peaks_s3_anno$annotation)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Promoter (.*)", "Promoter", H3K4me3_common_peaks_s3_anno$annotation)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Exon (.*)", "Exon", H3K4me3_common_peaks_s3_anno$annotation)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Downstream (.*)", "Downstream", H3K4me3_common_peaks_s3_anno$annotation)

H3K4me3_common_peaks_s3_anno_2 = H3K4me3_common_peaks_s3_anno[H3K4me3_common_peaks_s3_anno$annotation %in% "Promoter",]

my_gene = as.data.frame(unique(H3K4me3_common_peaks_s3_anno_2$geneId))
write.table(my_gene, "down_gene_in_promoter.txt", col.names = F, row.names = F,quote= F)



H3K4me3_common_peaks_s3_anno =annotatePeak(GRanges(H3K4me3_common_peaks_s2[,1:3]),tssRegion = c(-3000,3000),
                                           TxDb=txdb,
                                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
H3K4me3_common_peaks_s3_anno = as.data.frame(H3K4me3_common_peaks_s3_anno)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Intron (.*)", "Intron", H3K4me3_common_peaks_s3_anno$annotation)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Promoter (.*)", "Promoter", H3K4me3_common_peaks_s3_anno$annotation)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Exon (.*)", "Exon", H3K4me3_common_peaks_s3_anno$annotation)
H3K4me3_common_peaks_s3_anno$annotation=gsub("Downstream (.*)", "Downstream", H3K4me3_common_peaks_s3_anno$annotation)

H3K4me3_common_peaks_s3_anno_2 = H3K4me3_common_peaks_s3_anno[H3K4me3_common_peaks_s3_anno$annotation %in% "Promoter",]

my_gene = as.data.frame(unique(H3K4me3_common_peaks_s3_anno_2$geneId))
write.table(my_gene, "same_gene_in_promoter.txt", col.names = F, row.names = F,quote= F)





H3K4me3_common_peaks_s3$my_name = paste(H3K4me3_common_peaks_s3[,1],H3K4me3_common_peaks_s3[,2], H3K4me3_common_peaks_s3[,3], H3K4me3_common_peaks_s3[,15],H3K4me3_common_peaks_s3[,16],H3K4me3_common_peaks_s3[,17],    sep = "_")


H3K4me3_common_peaks_small  = H3K4me3_common_peaks_s3[,c(31,8,13,14,22,27,28)]

mm = melt(H3K4me3_common_peaks_small,id.vars =c("my_name", "pileup","WT_me3_rpkm","pileup.1","KO_me3_rpkm"))
mm$value = log2(mm$value)
library(ggpubr)
my_comparisons <- list( c("WT_me3_reads", "KO_me3_reads")) 

p = ggviolin(mm, x = "variable", y = "value", 
             fill = "variable",
             color =  "variable",
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#FC4E01","#FC4EA7"),
             add = "boxplot", add.params = list(fill = NULL,color = "black")
) 

p+stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 20)  



################   peak overlap between Menin and H3K4me3  ##########

library(tidyr)
H3k4me3_peak = read.table("/Volumes/My_disk/Wang_hb/Liu_my/H3K4me3_D8/2.cleandata/mouse/H3K4me3-WT_FKDL210271653-1a_peaks_annoed.xls")
menin_peak = read.table("/Volumes/My_disk/Wang_hb/Liu_my/Menin_ChIP_2nd/2.cleandata/Menin_normalized_2_peaks_specific_annoed.txt", sep = "\t", header = T, quote = "")

H3k4me3_peak_gr = GRanges(H3k4me3_peak)
menin_peak_gr = GRanges(menin_peak)

m=findOverlaps(H3k4me3_peak_gr, menin_peak_gr)
common = 10388

library(VennDiagram)
grid.newpage();
venn.plot <- draw.pairwise.venn(               area1           = length(H3k4me3_peak_gr),
                                               area2           = length(menin_peak_gr),
                                               cross.area      = common,
                                               category        = c("H3K4me3", "Menin"),
                                               fill            = c("#f35e5a", "#17b3b7"),
                                               scaled          = TRUE,
                                               lty             = "blank",
                                               cex             = c(1,1,1),
                                               # euler.d = T,
                                               ext.text        =0,
                                               cat.cex         = TRUE,
                                               cat.pos         = c(285, 105),
                                               cat.dist        = 0.09,
                                               alpha           = 0.8
                                               #cat.just        = list(c(-1, -1), c(1, 1))
                                               #ext.pos         = 30,
                                               #ext.dist        = -0.05,
                                               #ext.length      = 0.85,
                                               #ext.line.lwd    = 2,
                                               #ext.line.lty    = "dashed"
)






#######   downregulated gene and upregulated gene  ##########
H3K4me3_WT_anno_con_df$flag = "WT"
H3K4me3_KO_anno_con_df$flag = "KO"

H3K4me3_WT_KO = rbind(H3K4me3_WT_anno_con_df, H3K4me3_KO_anno_con_df)

DEG = read.table("/Users/wenbodeng/Desktop/D8 DEGs.txt", sep = "\t", header = T, quote = "")
DEG_up = DEG[DEG$logFC>0,]
DEG_down = DEG[DEG$logFC<0,]

H3K4me3_WT_KO_up = H3K4me3_WT_KO[H3K4me3_WT_KO$geneId %in% DEG_up$mgi_symbol,]
H3K4me3_WT_KO_down = H3K4me3_WT_KO[H3K4me3_WT_KO$geneId %in% DEG_down$mgi_symbol,]


library(ggpubr)
my_comparisons <- list( c("WT", "KO") )

H3K4me3_WT_KO_up$log_reads = log(H3K4me3_WT_KO_up$reads,2)
H3K4me3_WT_KO_up$log_rpkm = log(H3K4me3_WT_KO_up$rpkm,2)

H3K4me3_WT_KO_down$log_reads = log(H3K4me3_WT_KO_down$reads,2)
H3K4me3_WT_KO_down$log_rpkm = log(H3K4me3_WT_KO_down$rpkm,2)


p = ggviolin(H3K4me3_WT_KO_up, x = "flag", y = "log_reads", 
             fill = "flag",
             ylim = c(0,20),
             color =  "flag",
             palette = c("#00AFBB", "#E7B800"),
             add = "boxplot", 
             add.params = list(fill = NULL,color = "black")
             #trim = T
             
) 

p+stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 4)  



aa_small = aa[aa$peaks_1_geneId %in% DEG_up$mgi_symbol |aa$peaks_2_geneId %in% DEG_up$mgi_symbol ,] 

aa_small_s = aa_small[,c(9,14,15,34,39,40,50)]








##############  D6 RNA-Seq   12-2-2019  #########

setwd("/Volumes/My_disk/Wang_hb/Liu_my/D6_D8_RNASeq/CleanData/RNA_Seq")

library (DESeq2)
library(Rsubread)
library(edgeR)

######   STAR    #########

files=dir("/Volumes/My_disk/Wang_hb/Liu_my/D6_D8_RNASeq/CleanData/RNA_Seq", recursive= T, patter="accepted_hits.bam$")

files = files[grep("D6",files)]
files = files[c(3,4,1,2,7,8,5,6)]
Exp=featureCounts(files, 
                  annot.ext = "/Users/wenbodeng/E/Index/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf",
                  GTF.attrType="gene_name",
                  #GTF.attrType.extra = "gene_id",
                  isGTFAnnotationFile = TRUE,
                  nthreads = 24)

########        RPKM         #########   
D=DGEList(counts=Exp$counts)
D=calcNormFactors(D)
rpkm_exp=rpkm(D,gene.length = Exp$annotation$Length)
D$counts=rpkm_exp

D=estimateCommonDisp(D)
D=estimateTagwiseDisp(D)

D$samples$group=c(1:2)
et=exactTest(D)
p_val=cbind(2^-et$table$logFC,et$table$PValue)
colnames(p_val) = c("Fold","pValue")
data_value=cbind(D$pseudo.counts,p_val)

colnames(data_value) = c("D6KO_1","D6KO_2","D6WT_1","D6WT_2","Fold","pValue")

keep1=!grepl("Mir",row.names(data_value))

all_exp_1=data_value[keep1,]

all_exp_1 = as.data.frame(all_exp_1)

all_exp_1$Gene_ID = row.names(all_exp_1)

m=grepl("^Marc[h0-9]*$",all_exp_1$Gene_ID)  

all_exp_1[m,]$Gene_ID=  paste("'",all_exp_1[m,]$Gene_ID,sep = "")  

all_exp_2 = all_exp_1

ms=grepl("^Sept[0-9]*$",all_exp_2$Gene_ID)

all_exp_2[ms,]$Gene_ID =  paste("'",all_exp_2[ms,]$Gene_ID,sep = "")  

all_exp_3 = all_exp_2


write.table (all_exp_2, "D6_WT_KO_12022019.txt", sep = "\t")


md=grepl("^Sep15",all_exp_3$Gene_ID)

all_exp_3[md,]$Gene_ID =  paste("'",all_exp_3[md,]$Gene_ID,sep = "")  


#### volcano
deg = data_value[,5:6]
colnames(deg) = c("fold","pvalue")
deg = as.data.frame(deg)
deg$tag1 = cut(deg$fold, breaks = c(-Inf,0.5,2,Inf), labels = c(1,0,2))   
deg$tag2 = cut(deg$pvalue, breaks = c(-Inf,0.05,Inf), labels = c(3,0))
deg$tag3 = as.numeric(as.character(deg$tag1)) + as.numeric(as.character(deg$tag2))
deg$flag =cut(deg$tag3, breaks = c(-Inf,3,4,Inf), labels = c(0,1,2))
#write.table(deg,"volcanole.txt", sep = "\t")

library("RColorBrewer")
col = c( "grey", "#17b3b7", "#fc44ae")

v = ggplot(deg,aes(log2(fold),-log(pvalue,10), colour = factor(flag)))+geom_point()+
  scale_color_manual(values = col)
v
r04=v+xlim(-10,10)

r05=r04+geom_vline(xintercept=c(-1,1),linetype="dotted",size=0.7, col = "grey")+geom_hline(yintercept=-log10(0.05),col="grey",linetype="dotted",size=0.7)
r05
r05+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank()) 



####  PCA 8-2-2021  ###

library (DESeq2)
library(Rsubread)
library(edgeR)
library(dplyr)
library(ggplot2)

condition <- factor(c("WT_D6","WT_D6","KO_D6","KO_D6",
                      "WT_D8","WT_D8","KO_D8","KO_D8"
))
group_list = condition
#condition <- factor(file_name)
#cds <- DESeqDataSet(Exp$counts, conditions = Exp$targets)
sample_info <- data.frame(row.names = NULL, sample_name =c("WT_D6_1","WT_D6_2","KO_D6_1","KO_D6_2",
                                                           "WT_D8_1","WT_D8_2","KO_D8_1","KO_D8_2"
), group_list=condition)

dds = DESeqDataSetFromMatrix(countData = D$counts, 
                             colData = sample_info,
                             design = ~ group_list, 
                             tidy = FALSE)

#normalized_counts <- counts(dds, normalized=TRUE)

dds <- DESeq(dds) 
res <- results(dds)
dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd,intgroup=c("group_list"))
pcaData <- plotPCA(vsd, intgroup=c("group_list"), returnData=TRUE)

rld <- rlog(dds, blind = F)
plotPCA(rld,intgroup=c("group_list"))

pcaData <- plotPCA(rld,intgroup=c("group_list"), returnData=TRUE)

df <- rbind(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")
)

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

rld <- rlog(dds, blind = FALSE)
plotPCA(rld,intgroup=c("group_list"))


pcaData <- plotPCA(rld,intgroup=c("group_list"), returnData=TRUE)

library(ggplot2)

percentVar <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
library(ggrepel)

pcaData$name = c("WT_D6_1","WT_D6_2","KO_D6_1","KO_D6_2",
                 "WT_D8_1","WT_D8_2","KO_D8_1","KO_D8_2"
)
p=ggplot(pcaData, aes(PC1, PC2, color= group_list)) +scale_shape_identity()+
  geom_point(size= 4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")
  ) + 
  geom_text_repel(
    data = pcaData,
    aes(label = pcaData$name),
    size = 4,
    box.padding = unit(0.5, "lines"), 
    point.padding = unit(0.5, "lines")
  )


p+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank()) 







