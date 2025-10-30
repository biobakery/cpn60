setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/HMP2")

##dotplot MGX cpn vs. metaphlan3.1
mgx = fread("MGX/mgx_all_melt.tsv")
check = mgx[, c("variable", "value")] %>% group_by(variable) %>% summarise_all(sum)
check = as.data.frame(check)
mgx$sum = check[match(mgx$variable, check$variable),"value"]
mgx$RA = mgx$value/ mgx$sum
mgx$variable = gsub("_Abundance-RPKs", "", mgx$variable)
mgx[,c("variable","RA")] %>% group_by(variable)%>% summarise_all(sum) # check on the sum
mgx$library = "MGX"

metaphlan = fread("MGX/HMP2_metaphlan3.1.0_metaphlan_taxonomic_profiles.tsv")
names(metaphlan) = gsub("_taxonomic_profile", "", names(metaphlan))
sample = unique(mgx$variable)
sample = c("# taxonomy", sample)
metaphlan = as.data.frame(metaphlan)
metaphlan = metaphlan[ , c(sample)]
metaphlan_G = metaphlan[grepl("\\|g__", metaphlan$`# taxonomy`),]
unclassified = metaphlan[metaphlan$`# taxonomy` == "UNKNOWN",]
metaphlan_G = rbind(metaphlan_G, unclassified)
metaphlan_G = metaphlan_G[!grepl("\\|s__", metaphlan_G$`# taxonomy`),]
metaphlan_G = metaphlan_G %>% separate('# taxonomy', sep = "\\|", c("K", "P", "C", "O", "F", "G")) %>% select(-K, -P, -C, -O, -F)
metaphlan_G = as.data.frame(t(metaphlan_G))
names(metaphlan_G) = metaphlan_G[1,]
metaphlan_G = metaphlan_G[-1,]
metaphlan_G[] <- lapply(metaphlan_G, function(x) as.numeric(as.character(x)))
names(metaphlan_G)[274] = "unclassified"
metaphlan_G$sample = row.names(metaphlan_G)
metaphlan_G = melt(metaphlan_G)
names(metaphlan_G) = c("variable", "G", "RA")
metaphlan_G[,c("variable","RA")] %>% group_by(variable)%>% summarise_all(sum)
metaphlan_G = as.data.frame(metaphlan_G)
metaphlan_G$RA = metaphlan_G$RA /100
metaphlan_G$library = "metaphlan"

dotplot_metaphlan = rbind(mgx[,c("variable", "G", "RA", "library")], metaphlan_G[,c("variable", "G", "RA", "library")])
dotplot_metaphlan$sample_G_library = paste(dotplot_metaphlan$variable, dotplot_metaphlan$G, dotplot_metaphlan$library, sep = "-")
dotplot_G = dotplot_metaphlan[, c("sample_G_library","RA")] %>% group_by(sample_G_library) %>% summarise_all(sum) 
dotplot_G = as.data.frame(dotplot_G)
dotplot_G = dotplot_G %>% separate(sample_G_library, sep = "-", c("sample", "G", "library")) 
dotplot_G$sample_G = paste(dotplot_G$sample, dotplot_G$G, sep = "-")
dotplot_G = dotplot_G[,c("sample_G", "RA", "library")] %>% pivot_wider(names_from = sample_G, values_from = RA,  values_fill = list(n = 0))
dotplot_G = as.data.frame(t(dotplot_G))
names(dotplot_G) = dotplot_G[1,]
dotplot_G = dotplot_G[-1,]
dotplot_G$MGX = as.numeric(as.character(dotplot_G$MGX))
dotplot_G$metaphlan = as.numeric(as.character(dotplot_G$metaphlan))
dotplot_G[is.na(dotplot_G)] <- 0
write.table(dotplot_G, file = "MGX/dotplot_all_G.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

p = ggplot(dotplot_G, aes(x=MGX, y=metaphlan)) + geom_point(size = 2, color = "dark green", alpha=0.5) + scale_y_sqrt() + scale_x_sqrt() + theme_classic()+
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black")+
        labs(title="every genus in every sample",
             x="MGX cpn60 taxa RA", y = "Metaphlan3.1 RA")
ggsave("MGX/dotplot_metaphlan.pdf", width =6, height = 5)
#label the ratio >2 dots
dotplot_G$ratio = dotplot_G$metaphlan / dotplot_G$MGX
is.nan.data.frame <- function(x)
        do.call(cbind, lapply(x, is.nan))
dotplot_G[is.nan(dotplot_G)] <- 0
p = ggplot(dotplot_G, aes(x=MGX, y=metaphlan)) + geom_point(size = 2, color = "dark green", alpha=0.5) + scale_y_sqrt() + scale_x_sqrt() + theme_classic()+
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black")+
        labs(title="every genus in every sample",
             x="MGX cpn60 taxa RA", y = "Metaphlan3.1 RA")+ geom_text(aes(label=ifelse(ratio>2,as.character(G),'')),size = 2, hjust=0,vjust=1)
ggsave("dotplot_metaphlan_label.pdf", width =6, height = 5)


#remove the unclassified and plot again
dotplot_G = fread("Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/MGX/dotplot_G.csv")
dotplot_G = dotplot_G %>% separate(V1, sep = "-", c("sample", "G"))
dotplot_G = dotplot_G[!dotplot_G$G == "unclassified",]
dotplot_G$MGX_newSum = check[match(dotplot_G$sample, check$sample), "MGX"]
dotplot_G$metaphlan_newSum = check[match(dotplot_G$sample, check$sample), "metaphlan"]
dotplot_G$MGX_newRA = dotplot_G$MGX / dotplot_G$MGX_newSum
dotplot_G$metaphlan_newRA = dotplot_G$metaphlan / dotplot_G$metaphlan_newSum
check = dotplot_G[, c("sample", "MGX_newRA", "metaphlan_newRA")] %>% group_by(sample) %>% summarise_all(sum)
check = as.data.frame(check)

p = ggplot(dotplot_G, aes(x=MGX_newRA, y=metaphlan_newRA)) + geom_point(size = 2, color = "dark green", alpha=0.5) + scale_y_sqrt() + scale_x_sqrt() + theme_classic()+
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black")+
        labs(title="every genus in every sample",
             x="MGX cpn60 taxa RA", y = "Metaphlan3.1 RA")
ggsave("MGX/dotplot_metaphlan_classified.pdf", width =6, height = 5)
#write.csv(dotplot_G, "Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/MGX/dotplot_G_classified.csv")


##PCoA plot
metaphlan = fread("Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/MGX/dotplot_G_classified.csv")
cpn = metaphlan[, c("sample", "G", "MGX_newRA")]
metaphlan = metaphlan[,c("sample", "G", "metaphlan_newRA")]
cpn$sample = paste(cpn$sample, "cpn", sep = "-")
metaphlan$sample = paste(metaphlan$sample, "metaphlan", sep = "-")
names(cpn) = c("sample", "G", "RA")
names(metaphlan) = c("sample", "G", "RA")
tax_filt = rbind(metaphlan, cpn)
tax_filt = tax_filt %>% pivot_wider(names_from = G, values_from = RA,  values_fill = list(n = 0))
tax_filt = as.data.frame(tax_filt)
row.names(tax_filt) = tax_filt$sample
tax_filt$sample= NULL
meta = as.data.frame(row.names(tax_filt))
names(meta) = "temp"
meta$name = meta$temp
meta = meta %>% separate(temp, sep = "-", c("match", "method"))

bray = vegdist(tax_filt, "bray")
pc = capscale(bray~1, comm = tax_filt)
cap = data.frame(pc$CA$u)
cap = merge(meta, cap, by.x = "name", by.y = "row.names")
rownames(cap) = cap[,1]
cap = cap[,2:ncol(cap)]
s = summary(pc)
col = c("cpn" = "#486eb5", "metaphlan" = "#b8781f")
p = ggplot(cap, aes(MDS1, MDS2,  group = match))+ geom_line(size = 0.4) + geom_point(aes(color = method), size = 2, alpha = 0.6)  + theme_bw(base_size = 8)   + scale_color_manual(values = col) +
        labs(color = "method", title="", size = 8,
             x = paste("PC 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
             y = paste("PC 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))
p #4x3 inch

##permanova
library(dplyr)
library(vegan)
set.seed(123)  
meta2 <- meta %>% select(name, method, match) %>% distinct(name, .keep_all = TRUE)
rownames(meta2) <- meta2$name
meta2 <- meta2[labels(bray), , drop = FALSE]  # ensure same order as bray

## PERMANOVA
perm_res <- adonis2(bray ~ method, data = meta2,
                    permutations = 9999,
                    strata = meta2$match)

R2_method <- as.numeric(perm_res$R2[1])
p_method  <- as.numeric(perm_res$`Pr(>F)`[1])
cat(sprintf("PERMANOVA: R2 = %.3f, p = %.4f\n", R2_method, p_method))

##dotplot MGX vs. MTX
mgx = fread("MGX/mgx_all_melt.tsv")
mtx = fread("MTX/mtx_all_melt.tsv")
check = mgx[, c("variable", "value")] %>% group_by(variable) %>% summarise_all(sum)
check = as.data.frame(check)
mgx$sum = check[match(mgx$variable, check$variable),"value"]
mgx$RA = mgx$value/ mgx$sum
check = mtx[, c("variable", "value")] %>% group_by(variable) %>% summarise_all(sum)
check = as.data.frame(check)
mtx$sum = check[match(mtx$variable, check$variable),"value"]
mtx$RA = mtx$value/ mtx$sum

mgx$library = "MGX"
mtx$library = "MTX"
dotplot = rbind(mgx, mtx)
dotplot$variable = gsub("_Abundance-RPKs", "", dotplot$variable)
dotplot$sample_G_library = paste(dotplot$variable, dotplot$G, dotplot$library, sep = "-")
dotplot_G = dotplot[, c("sample_G_library","RA")] %>% group_by(sample_G_library) %>% summarise_all(sum) 
dotplot_G = as.data.frame(dotplot_G)
dotplot_G = dotplot_G %>% separate(sample_G_library, sep = "-", c("sample", "G", "library")) 
dotplot_G$sample_G = paste(dotplot_G$sample, dotplot_G$G, sep = "-")
dotplot_G = dotplot_G[,c("sample_G", "RA", "library")] %>% pivot_wider(names_from = sample_G, values_from = RA,  values_fill = list(n = 0))
dotplot_G = as.data.frame(t(dotplot_G))
names(dotplot_G) = dotplot_G[1,]
dotplot_G = dotplot_G[-1,]
dotplot_G$MGX = as.numeric(as.character(dotplot_G$MGX))
dotplot_G[is.na(dotplot_G)] <- 0
#write.csv(dotplot_G, "MTX/dotplot_all_G.tsv")
dotplot_G = fread("MTX/dotplot_all_G.tsv")

p = ggplot(dotplot_G, aes(x=MGX, y=MTX)) + geom_point(size = 2, color = "dark green", alpha=0.5) + scale_y_sqrt() + scale_x_sqrt() + theme_classic()+
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black")+
        labs(title="every genus in every sample",
             x="MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA")
ggsave("dotplot_G.pdf", width =6, height = 5)

#remove the unclassified and plot again
dotplot_G$temp = dotplot_G$V1
dotplot_G = dotplot_G %>% separate(temp, sep = "-", c("sample", "G"))
dotplot_G = dotplot_G[!dotplot_G$G == "unclassified",]
check = dotplot_G[, c("sample", "MGX", "MTX")] %>% group_by(sample) %>% summarise_all(sum)
check = as.data.frame(check)
dotplot_G$MGX_newSum = check[match(dotplot_G$sample, check$sample), "MGX"]
dotplot_G$MTX_newSum = check[match(dotplot_G$sample, check$sample), "MTX"]
dotplot_G$MGX_newMGX = dotplot_G$MGX / dotplot_G$MGX_newSum
dotplot_G$MTX_newMTX = dotplot_G$MTX / dotplot_G$MTX_newSum
p = ggplot(dotplot_G, aes(x=MGX_newMGX, y=MTX_newMTX)) + geom_point(size = 2, color = "dark green", alpha=0.5) + scale_y_sqrt() + scale_x_sqrt() + theme_classic()+
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black")+
        labs(title="every genus in every sample",
             x="MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA")
ggsave("dotplot_G_classified.pdf", width =6, height = 5)
##PCoA plot
MGX = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/MTX/dotplot_all_G_classified.tsv")
MGX = MGX[, c(2:10)] %>% separate(V1, sep = "-", c("sample", "G"))
MTX = MGX[, c(1,2,4)]
MGX = MGX[,c(1:3)]
MTX$sample = paste(MTX$sample, "MTX", sep = "-")
MGX$sample = paste(MGX$sample, "MGX", sep = "-")
names(MTX) = c("sample", "G", "RA")
names(MGX) = c("sample", "G", "RA")
tax_filt = rbind(MGX, MTX)
tax_filt = tax_filt %>% pivot_wider(names_from = G, values_from = RA,  values_fill = list(n = 0))
tax_filt = as.data.frame(tax_filt)
row.names(tax_filt) = tax_filt$sample
tax_filt$sample= NULL
meta = as.data.frame(row.names(tax_filt))
names(meta) = "temp"
meta$name = meta$temp
meta = meta %>% separate(temp, sep = "-", c("match", "method"))

tax_filt$sum = rowSums(tax_filt)
tax_filt = tax_filt[!tax_filt$sum == "NaN",]
tax_filt = na.omit(tax_filt)
tax_filt$sum = NULL
bray = vegdist(tax_filt, "bray")
pc = capscale(bray~1, comm = tax_filt)
cap = data.frame(pc$CA$u)
cap = merge(meta, cap, by.x = "name", by.y = "row.names")
rownames(cap) = cap[,1]
cap = cap[,2:ncol(cap)]
s = summary(pc)
col = c("MGX" = "#486eb5", "MTX" = "#a1312d")
p = ggplot(cap, aes(MDS1, MDS2,  group = match))+ geom_line(size = 0.4) + geom_point(aes(color = method), size = 2, alpha = 0.6)  + theme_bw(base_size = 8)   + scale_color_manual(values = col) +
        labs(color = "method", title="", size = 8,
             x = paste("PC 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
             y = paste("PC 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))#+geom_text(aes(label = match),hjust=0, vjust=0, size = 2)
p #4x3 inch

#adonis
metadata = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/hmp2_metadata.csv")
metadata = metadata[, c("External ID", "data_type", "IntervalName", "Education Level")]
metadata$data_type = gsub("metagenomics", "MGX", metadata$data_type)
metadata$data_type = gsub("metatranscriptomics", "MTX", metadata$data_type)
metadata$names = paste(metadata$`External ID`, metadata$data_type, sep = "-")
metadata = metadata[metadata$names %in% row.names(tax_filt),]
metadata=metadata[order(match(metadata$names, row.names(tax_filt))),]
metadata = as.data.frame(metadata)
row.names(metadata) = metadata$names
names(metadata) = c("sample", "data_type", "collection", "education", "names")

tax_filt$temp = rowSums(tax_filt)
tax_filt = tax_filt[!tax_filt$temp == 0,]
metadata = metadata[row.names(metadata) %in% row.names(tax_filt),]
tax_filt$temp = NULL
adonis_res_pval = vector()
adonis_res_rsq = vector()

for (col in names(metadata[, c(2:4)])){
        adonis.univ = adonis(as.formula(paste("tax_filt~", col)), data = metadata[, c(2:4)], permutations = 999, method = 'bray')
        adonis_res_pval[col] = adonis.univ$aov.tab[1,]$`Pr(>F)`
        adonis_res_rsq[col] = adonis.univ$aov.tab[1,]$R2
}

univar_res_tax = rbind(adonis_res_pval, adonis_res_rsq)
univar_res_tax = as.data.frame(t(univar_res_tax))
names(univar_res_tax) = c("P-Value", "R2")
univar_res_tax$`P-Value` = as.numeric(univar_res_tax$`P-Value`)
univar_res_tax$R2 = as.numeric(univar_res_tax$R2)
univar_res_tax$p_adj = p.adjust(univar_res_tax$`P-Value`, "fdr")

univar_res_tax #Univariate tests on Bray-Curtis dissimilarity of all samples in this study (n=192). Supplementary_3_AdonisTests
            P-Value        R2  p_adj
data_type    0.001 0.03792815 0.0015
collection   0.009 0.01884873 0.0090
education    0.001 0.08043966 0.0015
write.csv(univar_res_tax, "statistics.csv")



#2023-03-29

##cleveland plot
dotplot_G = fread("MTX/dotplot_all_G_classified.tsv")
share = dotplot_G[,2:ncol(dotplot_G)]
share = share[,c("G", "MGX_newMGX", "MTX_newMTX")] %>% group_by(G) %>% summarise_all(list(mean))
share = as.data.frame(share)
share$mean = (share$MGX_newMGX + share$MTX_newMTX) /2
share = share[!share$mean == 0,]
share = arrange(share, mean)
who_mean = share$G
share = melt(share[,1:3])
names(share) = c("G", "library", "value")
col = c("MGX_newMGX" = "#486eb5", "MTX_newMTX" = "#a1312d")
share = share[!grepl("unclassified", share$G),] #selecting only the rows where the # Gene Family column contains the "|" character

share$G = ordered(share$G, c(who_mean))
p = ggplot(share, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2.5, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Genus", color = "Dataset") + scale_x_log10() + theme_bw()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean))
ggsave("MTX/cleveland_all.pdf", width =6, height = 14)

#oral and gut bugs
oralvsgut = fread("oralvsgut.csv")
oralvsgut = oralvsgut %>% separate(feature, sep = "_", c("G", "S"))
oralvsgut$G = paste("g_", oralvsgut$G, sep = "_")
share$major_site = oralvsgut[match(share$G, oralvsgut$G), "major_site"]
share$library = gsub("_newMGX","", share$library)
share$library = gsub("_newMTX","", share$library)
share_oral = share[share$major_site == "oral",]
share_gut = share[share$major_site == "gut",]
share_oral = na.omit(share_oral)
share_gut = na.omit(share_gut)
share_oral$G = ordered(share_oral$G, c(who_mean))
share_gut$G = ordered(share_gut$G, c(who_mean))
col = c("MGX" = "#486eb5", "MTX" = "#a1312d")

p = ggplot(share_oral, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2.5, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Oral genus", color = "Dataset") + scale_x_log10() + theme_bw()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean))
ggsave("MTX/cleveland_oral.pdf", width =7, height = 5)

p = ggplot(share_gut, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2.5, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Gut genus", color = "Dataset") + scale_x_log10() + theme_bw()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean))
ggsave("MTX/cleveland_gut.pdf", width =7, height = 7)



#nost non viable 25
share = dotplot_G[,2:ncol(dotplot_G)]
share = share[,c("G", "MGX_newMGX", "MTX_newMTX")] %>% group_by(G) %>% summarise_all(list(mean))
share = as.data.frame(share)
share$ratio = share$MTX_newMTX / share$MGX_newMGX
share = arrange(share, ratio)
share = share[!share$ratio == "NaN",]
share = share[!share$ratio ==0,]
who_ratio = share$G[1:11] #the smallest ratio (non-viable)
share = melt(share[,1:3])
share = share[share$G %in% who_ratio,]
names(share) = c("G", "library", "value")
share = share[!share$G=="g__Clostridiales_Family_XIII_Incertae_Sedis_unclassified",]
share$G = ordered(share$G, c(who_mean))
p = ggplot(share, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Genus", color = "Dataset") + scale_x_continuous(trans = "sqrt", labels = scales::scientific_format())+ theme_bw()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean)) 
ggsave("MTX/cleveland_low_nonZero.pdf", width =7.4, height = 3)
#nost viable 25
share = dotplot_G[,2:ncol(dotplot_G)]
share = share[,c("G", "MGX_newMGX", "MTX_newMTX")] %>% group_by(G) %>% summarise_all(list(mean))
share = as.data.frame(share)
share$ratio = share$MTX_newMTX / share$MGX_newMGX
share = arrange(share, ratio)
who_ratio = share$G[111:122] #the smallest ratio (non-viable)
share = melt(share[,1:3])
share = share[share$G %in% who_ratio,]
names(share) = c("G", "library", "value")
share$G = ordered(share$G, c(who_mean))
p = ggplot(share, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Genus", color = "Dataset") + scale_x_sqrt() +theme_bw()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean))
ggsave("MTX/cleveland_top.pdf", width =7, height = 3)

#least MGX abundant 25
share = dotplot_G[,2:ncol(dotplot_G)]
share = share[,c("G", "MGX_newMGX", "MTX_newMTX")] %>% group_by(G) %>% summarise_all(list(mean))
share = as.data.frame(share)
share = arrange(share, MGX_newMGX)
share = share[!share$G==0,]
who = share$G[1:20] #the smallest RA
share = melt(share[,1:3])
share = share[share$G %in% who,]
names(share) = c("G", "library", "value")
share$G = ordered(share$G, c(who_mean))
p = ggplot(share, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Genus", color = "Dataset") + scale_x_sqrt()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean))
ggsave("MTX/cleveland_low.pdf", width =7, height = 4)
#nost viable 25
share = dotplot_G[,2:ncol(dotplot_G)]
share = share[,c("G", "MGX_newMGX", "MTX_newMTX")] %>% group_by(G) %>% summarise_all(list(mean))
share = as.data.frame(share)
share = arrange(share, MGX_newMGX)
who = share$G[129:149] #the smallest ratio (non-viable)
share = melt(share[,1:3])
share = share[share$G %in% who,]
names(share) = c("G", "library", "value")
share$G = ordered(share$G, c(who_mean))
p = ggplot(share, aes(x = value, y = G, color = library, group = G))  +  geom_point(size = 2, alpha = 0.7) +geom_line() + scale_color_manual(values = col)+
        labs(x = "Value", y = "Genus", color = "Dataset") + scale_x_sqrt() +theme_bw()
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who_mean))
ggsave("MTX/cleveland_top.pdf", width =7, height = 4)

dotplot_G = fread("MTX/dotplot_all_G_classified.tsv")

#plot the most life and most dead bugs
dotplot_Barnesiella = dotplot_G[dotplot_G$G=="g__Barnesiella",2:10]
cor_test <- cor.test(dotplot_Barnesiella$MGX_newMGX,dotplot_Barnesiella$MTX_newMTX)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_newMTX ~ MGX_newMGX, data = dotplot_Barnesiella)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
write.csv(dotplot_Barnesiella, "MTX/dotplot_Barnesiella.csv")
ggplot(dotplot_Barnesiella, aes(x = MGX_newMGX, y = MTX_newMTX)) +geom_point(color = "dark green", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightgreen") +
        annotate("text", x = max(dotplot_Barnesiella$MGX_newMGX), y = max(dotplot_Barnesiella$MTX_newMTX), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(dotplot_Barnesiella$MGX_newMGX), y = min(dotplot_Barnesiella$MTX_newMTX),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MGX", y = "MTX", title = "Barnesiella") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
ggsave("MTX/dotplot_Barnesiella.pdf", width =3.5, height = 3.5)

dotplot_Odoribacter = dotplot_G[dotplot_G$G=="g__Odoribacter",2:10]
cor_test <- cor.test(dotplot_Odoribacter$MGX_newMGX,dotplot_Odoribacter$MTX_newMTX)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_newMTX ~ MGX_newMGX, data = dotplot_Odoribacter)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(dotplot_Odoribacter, aes(x = MGX_newMGX, y = MTX_newMTX)) +geom_point(color = "dark green", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightgreen") +
        annotate("text", x = max(dotplot_Odoribacter$MGX_newMGX), y = max(dotplot_Odoribacter$MTX_newMTX), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(dotplot_Odoribacter$MGX_newMGX), y = min(dotplot_Odoribacter$MTX_newMTX),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MGX", y = "MTX", title = "Odoribacter") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
ggsave("MTX/dotplot_Odoribacter.pdf", width =4, height = 3.5)

dotplot_Sutterella = dotplot_G[dotplot_G$G=="g__Sutterella",2:10]
cor_test <- cor.test(dotplot_Sutterella$MGX_newMGX,dotplot_Sutterella$MTX_newMTX)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_newMTX ~ MGX_newMGX, data = dotplot_Sutterella)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(dotplot_Sutterella, aes(x = MGX_newMGX, y = MTX_newMTX)) +geom_point(color = "dark green", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightgreen") +
        annotate("text", x = max(dotplot_Sutterella$MGX_newMGX), y = max(dotplot_Sutterella$MTX_newMTX), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(dotplot_Sutterella$MGX_newMGX), y = min(dotplot_Sutterella$MTX_newMTX),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MGX", y = "MTX", title = "Sutterella") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
ggsave("MTX/dotplot_Sutterella.pdf", width =4, height = 3.5)

dotplot_Bacteroides = dotplot_G[dotplot_G$G=="g__Bacteroides",2:10]
cor_test <- cor.test(dotplot_Bacteroides$MGX_newMGX,dotplot_Bacteroides$MTX_newMTX)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_newMTX ~ MGX_newMGX, data = dotplot_Bacteroides)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
write.csv(dotplot_Barnesiella, "MTX/dotplot_Bacteroides.csv")
ggplot(dotplot_Bacteroides, aes(x = MGX_newMGX, y = MTX_newMTX)) +geom_point(color = "dark green", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightgreen") +
        annotate("text", x = max(dotplot_Bacteroides$MGX_newMGX), y = max(dotplot_Bacteroides$MTX_newMTX), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(dotplot_Bacteroides$MGX_newMGX), y = min(dotplot_Bacteroides$MTX_newMTX),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MGX", y = "MTX", title = "Bacteroides") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
ggsave("MTX/dotplot_Bacteroides.pdf", width =3.5, height = 3.5)

dotplot_Dielma = dotplot_G[dotplot_G$G=="g__Dielma",2:10]
cor_test <- cor.test(dotplot_Dielma$MGX_newMGX,dotplot_Dielma$MTX_newMTX)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_newMTX ~ MGX_newMGX, data = dotplot_Dielma)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(dotplot_Dielma, aes(x = MGX_newMGX, y = MTX_newMTX)) +geom_point(color = "dark green", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightgreen") +
        annotate("text", x = max(dotplot_Dielma$MGX_newMGX), y = max(dotplot_Dielma$MTX_newMTX), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(dotplot_Dielma$MGX_newMGX), y = min(dotplot_Dielma$MTX_newMTX),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MGX", y = "MTX", title = "Dielma") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
ggsave("MTX/dotplot_Dielma.pdf", width =4, height = 3.5)

dotplot_Dialister = dotplot_G[dotplot_G$G=="g__Dialister",2:10]
cor_test <- cor.test(dotplot_Dialister$MGX_newMGX,dotplot_Dialister$MTX_newMTX)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_newMTX ~ MGX_newMGX, data = dotplot_Dialister)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
write.csv(dotplot_Dialister, "MTX/dotplot_Dialister.csv")
ggplot(dotplot_Dialister, aes(x = MGX_newMGX, y = MTX_newMTX)) +geom_point(color = "dark green", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightgreen") +
        annotate("text", x = max(dotplot_Dialister$MGX_newMGX), y = max(dotplot_Dialister$MTX_newMTX), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(dotplot_Dialister$MGX_newMGX), y = min(dotplot_Dialister$MTX_newMTX),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MGX", y = "MTX", title = "Dialister") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
ggsave("MTX/dotplot_Dialister.pdf", width =3.5, height = 3.5)


###Get the dotplot and cleveland plot in species
mgx = fread("MGX/mgx_all_melt.tsv")
mtx = fread("MTX/mtx_all_melt.tsv")
check = mgx[, c("variable", "value")] %>% group_by(variable) %>% summarise_all(sum)
check = as.data.frame(check)
mgx$sum = check[match(mgx$variable, check$variable),"value"]
mgx$RA = mgx$value/ mgx$sum
check = mtx[, c("variable", "value")] %>% group_by(variable) %>% summarise_all(sum)
check = as.data.frame(check)
mtx$sum = check[match(mtx$variable, check$variable),"value"]
mtx$RA = mtx$value/ mtx$sum

mgx$library = "MGX"
mtx$library = "MTX"
dotplot = rbind(mgx, mtx)
dotplot$variable = gsub("_Abundance-RPKs", "", dotplot$variable)
dotplot$sample_S_library = paste(dotplot$variable, dotplot$S, dotplot$library, sep = "-")
dotplot_S = dotplot[, c("sample_S_library","RA")] %>% group_by(sample_S_library) %>% summarise_all(sum) 
dotplot_S = as.data.frame(dotplot_S)
dotplot_S = dotplot_S %>% separate(sample_S_library, sep = "-", c("sample", "S", "library")) 
dotplot_S$sample_S = paste(dotplot_S$sample, dotplot_S$S, sep = "-")
dotplot_S = dotplot_S[,c("sample_S", "RA", "library")] %>% pivot_wider(names_from = sample_S, values_from = RA,  values_fill = list(n = 0))
dotplot_S = as.data.frame(t(dotplot_S))
names(dotplot_S) = dotplot_S[1,]
dotplot_S = dotplot_S[-1,]
dotplot_S$MGX = as.numeric(as.character(dotplot_S$MGX))
dotplot_S[is.na(dotplot_S)] <- 0
#write.table(dotplot_S, "MTX/dotplot_all_S.tsv", sep = "\t", row.names = TRUE, quote = FALSE) 
