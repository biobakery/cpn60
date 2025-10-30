library("stringr")
library("dplyr")
library("ggplot2")
library("GUniFrac")
library("ape")
library("matrixStats")
library("phytools")
library("RColorBrewer")
library("reshape2")
library("scales")
library("gplots")
library("gtools")
library("vegan")
library("phyloseq")
library("data.table")
library("tidyr")
library("seqinr")
library("ggtree")

setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/full_length_tree/translated1k-3k_MSA")
dist = fread("cpn_gene_fullLength_AA_cleaned_dist.csv")
id = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/cpn_gene_20556_utaxID.csv", sep = ",", header = T)
dist$P = id[match(dist$V1, id$id), "P"]
dist$P = gsub("p:", "", dist$P)
who = c("Actinobacteria", "Bacteroidetes","Cyanobacteria","Firmicutes", "Proteobacteria")
dist_1 = dist[dist$P %in% who,]
dist_1$P = ordered(dist_1$P, who)
dist_1$V1 = gsub(";.*", "", dist_1$V1)
names(dist_1) = gsub(";.*", "", names(dist_1))
dist_1 = arrange(dist_1, P)
library(pheatmap)
col_order = c(dist_1$V1)
dist_1 = as.data.frame(dist_1)
dist_1 = dist_1[,col_order]

heatmap(as.matrix(dist_1), 
        Rowv=TRUE, 
        Colv=TRUE, 
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale="none",
        margins=c(10, 10),
        main="Heatmap of Identities")
ggsave("heatmap.pdf", width = 6, height = 6, units = "in", dpi = 300)



library(ComplexHeatmap)
library(seriation)
library(magick)
library(circlize)
col_Taxa = c("Actinobacteria" ="#e3857b" , "Bacteroidetes"="#ba9e41", "Cyanobacteria" = "#51b848","Firmicutes" ="#5bc9b1" , "Proteobacteria" ="#d689d4" ,"Other" = "#76bed6" )
ordered.sample.anno = HeatmapAnnotation( df = temp$P,
                                         name=c("P"),
                                         col = list(P = col_Taxa),
                                         gap=unit(0.2,"mm"),
                                         show_annotation_name = TRUE)
p = Heatmap(matrix_cpn, col=colorRamp2(c(50,75,100),c("blue","white", "red")),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            show_column_names=FALSE,
            show_row_names=FALSE,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 8),
            row_title="Taxa",
            row_title_side = "right",
            heatmap_legend_param = list(at = c(50,75,100),
                                        labels = c("0.5",  "0.75", "1.0"),
                                        color_bar = "continuous",
                                        legend_direction="horizontal",
                                        labels_gp = gpar(fontsize = 8),
                                        legend_width = unit(6, "cm"),
                                        title_position = "topcenter",
                                        title = "Pairwise identity"),
            top_annotation = ordered.sample.anno,
            width = unit(40, "cm"),
            rect_gp = gpar(col = "grey", lwd = 0.5))
plot(p)

pdf(plot(p), '~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Cpn60_paper_20220727/heatmap_cpn.pdf', width = 1400, height = 800, pointsize = 12,
    bg = "white")
dev.off()


###boxplot of the per-phylum average identities otus
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/full_length_tree/1k-3k_tree")
otu = fread("/net/hutlab11/srv/export/hutlab11/share_root/users/lea/data/cpn60-seq/BioProjectPRJEB43503/cpn_gene_23545_0.9/all_samples_taxonomy_closed_reference_90_newDB.tsv")
taxa = fread("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/cpn_gene_23545_utaxID.tsv")

otu = otu$taxonomy
dist = fread("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/amplicon_tree/17678Tree/cpn_gene_23545_amplicon_dist.csv")
dist = as.data.frame(dist)
dist = melt(dist)
names(dist) = c("Ref", "OTU", "Value")
dist = dist[dist$OTU %in% otu,] #100 out of 113 OTUs were aligned in the dist file here
dist$Ref_P = taxa[match(dist$Ref, taxa$id), "P"]
dist$OTU_P = taxa[match(dist$OTU, taxa$id), "P"]
dist$OTU_P = dist$OTU_P$P
dist$Ref_P = dist$Ref_P$P
dist = dist[!dist$Ref == dist$OTU,]

dist_otu = dist[, c("OTU", "Ref_P", "Value")] %>% group_by(OTU, Ref_P) %>% summarise_all(median)
dist_otu = as.data.frame(dist_otu)
dist_otu$OTU_P = dist[match(dist_otu$OTU, dist$OTU), "OTU_P"]
who = c("p:Actinobacteria", "p:Bacteroidetes", "p:Firmicutes", "p:Proteobacteria")
dist_otu = dist_otu[dist_otu$Ref_P %in% who,]

setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/nucleotideMSA/")
dist_otu = fread("dist_otu_UT.csv")
dist_Firmicutes = dist_otu[dist_otu$OTU_P =="p:Firmicutes",]
dist_Firmicutes$color = "A"
dist_Firmicutes$color[grep("p:Firmicutes", dist_Firmicutes$Ref_P)] = "B"
dist_Firmicutes$Ref_P = ordered(dist_Firmicutes$Ref_P, c(who))
p = ggplot(dist_Firmicutes, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("Average identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(40, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p



#boxplot of the 90th identities
dist_otu = dist[, c("OTU", "Ref_P", "Value")] %>% group_by(OTU, Ref_P) %>%summarise("value (the 90th percentile)" = quantile(Value, probs=0.9, na.rm=TRUE))
dist_otu = as.data.frame(dist_otu)
dist_otu$OTU_P = dist[match(dist_otu$OTU, dist$OTU), "OTU_P"]
who = c("p:Actinobacteria", "p:Bacteroidetes", "p:Firmicutes", "p:Proteobacteria")
dist_otu = dist_otu[dist_otu$Ref_P %in% who,]
write.csv(dist_otu, "dist_otu_90th.csv")

setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/nucleotideMSA/")
dist_otu = fread("dist_otu_90th.csv")
names(dist_otu) = c("V1","OTU", "Ref_P", "Value", "OTU_P")
dist_Firmicutes = dist_otu[dist_otu$OTU_P =="p:Firmicutes",]
dist_Firmicutes$color = "A"
dist_Firmicutes$color[grep("p:Firmicutes", dist_Firmicutes$Ref_P)] = "B"
dist_Firmicutes$Ref_P = ordered(dist_Firmicutes$Ref_P, c(who))
p = ggplot(dist_Firmicutes, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("90th percentile pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(40, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p

#boxplot of the pairwise percentile identities
who = c("p:Actinobacteria", "p:Bacteroidetes", "p:Firmicutes", "p:Proteobacteria")
dist_otu = dist[dist$Ref_P %in% who, ]
write.csv(dist_otu, "dist_otu_pairwise.csv")
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/nucleotideMSA/")
dist_otu = fread("dist_otu_pairwise.csv")
dist_Actinobacteria = dist_otu[dist_otu$OTU_P =="p:Actinobacteria",]
dist_Actinobacteria$color = "A"
dist_Actinobacteria$color[grep("p:Actinobacteria", dist_Actinobacteria$Ref_P)] = "B"
dist_Actinobacteria$Ref_P = ordered(dist_Actinobacteria$Ref_P, c(who))
p = ggplot(dist_Actinobacteria, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.3) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(0, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p

#statistics
library(tidyverse)
library(ggpubr)
dist_otu <- read_csv("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/nucleotideMSA/dist_otu_UT.csv", col_types = cols(
        OTU        = col_character(),
        Ref_P      = col_character(),
        Value      = col_double(),
        OTU_P      = col_character()))

df2 <- dist_otu %>% mutate(comparison = if_else(Ref_P == OTU_P, "within", "between"))
df2 %>%group_by(comparison) %>%summarize(
        mean_id = mean(Value),
        sd_id   = sd(Value),
        n       = n())
withins  <- df2$Value[df2$comparison == "within"]
betweens <- df2$Value[df2$comparison == "between"]

wilcox.test(x = withins, y = betweens, alternative = "greater")

###do the boxplot for 16s (silva database)
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/SILVA/V4_tree")
otu = fread("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/Biobakery_output/biobakery_output_16s_90/all_samples_taxonomy_closed_reference_silva.tsv")
taxa = fread("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/SILVA/silva_132_97_utaxID.csv")
otu = otu$`# OTU`

dist = fread("Silva_20305_V4_dist.csv")
dist = as.data.frame(dist)
dist = melt(dist)
names(dist) = c("Ref", "OTU", "Value")
dist = dist[dist$OTU %in% otu,] #335 out of 341 OTUs were aligned in the dist file here
dist$Ref_P = taxa[match(dist$Ref, taxa$original_names), "P"]
dist$OTU_P = taxa[match(dist$OTU, taxa$original_names), "P"]
dist$OTU_P = dist$OTU_P$P
dist$Ref_P = dist$Ref_P$P
dist = dist[!dist$Ref == dist$OTU,]

dist_otu = dist[, c("OTU", "Ref_P", "Value")] %>% group_by(OTU, Ref_P) %>%summarise("value (the 90th percentile)" = quantile(Value, probs=0.9, na.rm=TRUE))
dist_otu = as.data.frame(dist_otu)
dist_otu$OTU_P = dist[match(dist_otu$OTU, dist$OTU), "OTU_P"]
dist_otu$Ref_P = gsub(",p:", "p:", dist_otu$Ref_P)
dist_otu$OTU_P = gsub(",p:", "p:", dist_otu$OTU_P)
who = c("p:Actinobacteria", "p:Bacteroidetes", "p:Firmicutes", "p:Proteobacteria")
dist_otu = dist_otu[dist_otu$Ref_P %in% who,]
write.csv(dist_otu, "dist_otu_90th.csv")

setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/nucleotideMSA/")
dist_otu = fread("silva_dist_otu_90th.csv")
names(dist_otu) = c("V1","OTU", "Ref_P", "Value", "OTU_P")
dist_Actinobacteria = dist_otu[dist_otu$OTU_P =="p:Actinobacteria",]
dist_Actinobacteria$color = "A"
dist_Actinobacteria$color[grep("p:Actinobacteria", dist_Actinobacteria$Ref_P)] = "B"
dist_Actinobacteria$Ref_P = ordered(dist_Actinobacteria$Ref_P, c(who))
p = ggplot(dist_Actinobacteria, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("90th percentile pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(40, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p

##silva pairwise
who = c("p:Actinobacteria", "p:Bacteroidetes", "p:Firmicutes", "p:Proteobacteria")
dist$Ref_P = gsub(",p:", "p:", dist$Ref_P)
dist$OTU_P = gsub(",p:", "p:", dist$OTU_P)
dist_otu = dist[dist$Ref_P %in% who, ]
write.csv(dist_otu, "silva_dist_otu_pairwise.csv")

setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/nucleotideMSA/")
dist_otu = fread("silva_dist_otu_pairwise.csv")
dist_Actinobacteria = dist_otu[dist_otu$OTU_P =="p:Actinobacteria",]
dist_Actinobacteria$color = "A"
dist_Actinobacteria$color[grep("p:Actinobacteria", dist_Actinobacteria$Ref_P)] = "B"
dist_Actinobacteria$Ref_P = ordered(dist_Actinobacteria$Ref_P, c(who))
p = ggplot(dist_Actinobacteria, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("90th percentile pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(0, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p

dist_Bacteroidetes = dist_otu[dist_otu$OTU_P =="p:Bacteroidetes",]
dist_Bacteroidetes$color = "A"
dist_Bacteroidetes$color[grep("p:Bacteroidetes", dist_Bacteroidetes$Ref_P)] = "B"
dist_Bacteroidetes$Ref_P = ordered(dist_Bacteroidetes$Ref_P, c(who))
p = ggplot(dist_Bacteroidetes, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("90th percentile pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(0, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p

dist_Firmicutes = dist_otu[dist_otu$OTU_P =="p:Firmicutes",]
dist_Firmicutes$color = "A"
dist_Firmicutes$color[grep("p:Firmicutes", dist_Firmicutes$Ref_P)] = "B"
dist_Firmicutes$Ref_P = ordered(dist_Firmicutes$Ref_P, c(who))
p = ggplot(dist_Firmicutes, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("90th percentile pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(0, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p
dist_Proteobacteria = dist_otu[dist_otu$OTU_P =="p:Proteobacteria",]
dist_Proteobacteria$color = "A"
dist_Proteobacteria$color[grep("p:Proteobacteria", dist_Proteobacteria$Ref_P)] = "B"
dist_Proteobacteria$Ref_P = ordered(dist_Proteobacteria$Ref_P, c(who))
p = ggplot(dist_Proteobacteria, aes(x = Ref_P, y = Value, color = color)) + 
    geom_jitter(aes(color = color), position=position_jitter(0.2), size = 2, alpha = 0.5) +
    geom_boxplot(position=position_dodge())  +theme_bw()+
    xlab("Phyla") + ylab("90th percentile pairwise identities") + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0)) + 
    ylim(50, 100) + theme(legend.position = "none")
p$data$Ref_P = factor(p$data$Ref_P, ordered = TRUE, levels = who)
p