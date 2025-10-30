#Per phylum heatmap
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/full_length_tree/1k-3k_tree")
dist = fread("cpn_fullLength_cleaned_dist.csv")
dist = melt(dist)
dist = dist[!dist$V1 == dist$variable,]
id = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/full_length_tree/id_20556_utax.csv", sep = ",", header = T)
id = id %>% separate(id, sep = ",", c("D", "P", "C", "O", "F", "G", "S"))
id$P = gsub("p:", "", id$P)
dist$V1_P = id[match(dist$V1, id$Taxa), "P"]
dist$variable_P = id[match(dist$variable, id$Taxa), "P"]
dist_P = dist[,3:5] %>% group_by(V1_P, variable_P) %>% summarise_all(list(median))
dist_P = as.data.frame(dist_P)
dist_P = dist_P[!dist_P$V1_P == "Bacteria_unclassified",]
dist_P = dist_P[!dist_P$variable_P == "Bacteria_unclassified",]
who = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Abditibacteriota","Acidobacteria","Aquificae","Armatimonadetes","Ascomycota","Balneolaeota","Calditrichaeota","candidate_division_AD3","candidate_division_NC10" ,"candidate_division_WPS_2", "candidate_division_Zixibacteria","Candidatus_Aminicenantes","Candidatus_Cloacimonetes","Candidatus_Eisenbacteria", "Candidatus_Goldbacteria","Candidatus_Riflebacteria", "Candidatus_Rokubacteria","Candidatus_Sumerlaeota", "Candidatus_Tectomicrobia","Candidatus_Wallbacteria", "Chlamydiae","Chlorobi", "Chloroflexi","Chrysiogenetes", "Deferribacteres","Deinococcus_Thermus","Elusimicrobia",  "Euryarchaeota","Fibrobacteres", "Fusobacteria", "Gemmatimonadetes" ,"Ignavibacteriae",  "Lentisphaerae","Nitrospinae"   ,  "Nitrospirae","Planctomycetes","Rhodothermaeota", "Spirochaetes","Synergistetes", "Thermodesulfobacteria","Thermotogae", "Verrucomicrobia")
dist_P = dist_P %>% pivot_wider(names_from = variable_P, values_from = value,  values_fill = list(n = 0))
dist_P = as.data.frame(dist_P)
row.names(dist_P) = dist_P$V1_P
dist_P$V1_P = NULL
dist_P = dist_P[who, who]
dist_P = na.omit(dist_P)
dist_P = dist_P[, row.names(dist_P)]
#heatmap
heatmap(as.matrix(dist_P),Rowv = NA,Colv = NA) # draw a heatmap without tree
ggsave("phylum_heatmap.pdf", width = 6, height = 6, units = "in", dpi = 300)

write.csv(dist_P, "dist_phylum.csv")

##every species heatmap, align the heatmap by the order or branches in the tree
x = read.tree("cpn_fullLength_tree.newick")
id_P = id[,c("Taxa", "P")]
id_P$Color = "Other"
id_P$Color[grep("Actinobacteria", id_P$P)] = "Actinobacteria"
id_P$Color[grep("Bacteroidetes", id_P$P)] = "Bacteroidetes"
id_P$Color[grep("Firmicutes", id_P$P)] = "Firmicutes"
id_P$Color[grep("Proteobacteria", id_P$P)] = "Proteobacteria"
id_P$Color[grep("Cyanobacteria", id_P$P)] = "Cyanobacteria"
p = ggtree(x, size = 0.4) %<+% id_P + geom_tippoint(aes(color=Color), size=2, alpha = 0.5)+ geom_treescale()

dist = fread("cpn_fullLength_cleaned_dist.csv")
dist <- replace(dist, is.na(dist), 100)
tip_labels <- x$tip.label
dist = as.data.frame(dist)
row.names(dist) = dist$V1
dist$V1 = NULL
order <- match(tip_labels, rownames(dist))
dist <- dist[tip_labels, tip_labels]
matrix = as.matrix(dist)

heatmap(matrix,Rowv = NA,Colv = NA) 
heatmap(matrix, col = heat.colors(256), xlab = "", ylab = "", main = "Heatmap")
ggsave("heatmap.pdf", width = 10, height = 10, units = "in", dpi = 300)

library(ComplexHeatmap)
library(seriation)
library(magick)
library(circlize)
id_P = id_P[id_P$Taxa %in% names(dist),]
who = row.names(dist)
row.names(id_P) = id_P$Taxa
id_P = id_P[who,]
col_Taxa = c("Actinobacteria" ="#e3857b" , "Bacteroidetes"="#ba9e41", "Cyanobacteria" = "#51b848","Firmicutes" ="#5bc9b1" , "Proteobacteria" ="#d689d4" ,"Other" = "#76bed6" )

ordered.sample.anno = HeatmapAnnotation( Color = id_P$Color,
                                         name=c("Phylum"),
                                         col = list(Color = col_Taxa),
                                         gap=unit(0.2,"mm"),
                                         show_annotation_name = TRUE)

p = Heatmap(matrix, col=colorRamp2(c(50,75,100),c("blue","white", "red")),
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
                                        labels_gp = gpar(fontsize = 6),
                                        legend_width = unit(4, "cm"),
                                        title_position = "topcenter",
                                        title = "Pairwise identity"),
            top_annotation = ordered.sample.anno,
            width = unit(30, "cm"),
            rect_gp = gpar(col = "grey", lwd = 0.5))
pdf('heatmap_cpn.pdf', width = 1400, height = 800, pointsize = 12, bg = "white")
print(p)
dev.off()

