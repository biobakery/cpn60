library("stringr")
library("dplyr")
library("ggplot2")
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
library(ggpmisc)
library("Bios2cor")
library("Biostrings")

##use the otu from our biobakery workflow (uparse), classified using the RDP classifier
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/")
otu = read.csv("original/otu_RDPClassifier_all.csv", sep = ",", header = F)
otu = otu[, c(1,6,8,9,11,12,14,15,17,18,20,21,23,24,26)]
names(otu) = c("OTU", "K", "K_value", "P", "P_value", "C", "C_value", "O", "O_value", "F", "F_value", "G", "G_value", "S", "S_value")
otu$K  = ifelse(otu$K_value < 0.8, "", otu$K)
otu$P  = ifelse(otu$P_value < 0.8, "", otu$P)
otu$C  = ifelse(otu$C_value < 0.8, "", otu$C)
otu$O  = ifelse(otu$O_value < 0.8, "", otu$O)
otu$`F`  = ifelse(otu$F_value < 0.8, "", otu$`F`)
otu$G  = ifelse(otu$G_value < 0.8, "", otu$G)
otu$S  = ifelse(otu$S_value < 0.8, "", otu$S)

otu$K = paste("k",otu$K, sep = ":" )
otu$P = paste("p",otu$P, sep = ":" )
otu$C = paste("c",otu$C, sep = ":" )
otu$O = paste("o",otu$O, sep = ":" )
otu$`F` = paste("f",otu$`F`, sep = ":" )
otu$G = paste("g",otu$G, sep = ":" )
otu$S = paste("s",otu$S, sep = ":" )

otu$taxonomy = paste(otu$K, otu$P, otu$C, otu$O, otu$`F`, otu$G, otu$S, sep = ",")
otu = otu %>% separate(OTU, sep = ";", c("OTU", "seqs"))

all_samples = read.csv("original/all_samples_otu_mapping_results_all.tsv", sep = "\t", header = T)
names(all_samples)[1] = c("OTU")
names(all_samples) = gsub("_cpn_R1_001", "", names(all_samples))

all_samples$taxonomy = otu[match(all_samples$OTU, otu$OTU), "taxonomy"]
all_samples$P = otu[match(all_samples$OTU, otu$OTU ), "P"]
all_samples$P  = ifelse(all_samples$taxonomy=="k:,p:,c:,o:,f:,g:,s:", "unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:,p:,c:,o:,f:,g:,s:NA", "unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:,p:,c:,o:,f:,g:NA,s:NA", "unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:", "p:unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:NA", "p:unclassified", all_samples$P)

otu_taxonomy = all_samples[, c(2:112,113)]
write.csv(otu_taxonomy, "otu_taxonomy.csv")

all_samples_P = all_samples[, c(2:112,114)]
all_samples_P = all_samples_P %>% group_by(P)%>% summarise_all(list(sum))
all_samples_P = as.data.frame(all_samples_P)
row.names(all_samples_P) = all_samples_P$P
all_samples_P$P = NULL
all_samples_P = as.data.frame(t(all_samples_P))
all_samples_P = all_samples_P / rowSums(all_samples_P)
all_samples_P$sample = row.names(all_samples_P)
all_samples_P = melt(all_samples_P)
who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Deinococcus-Thermus", "p:Nitrospirae","p:Spirochaetes","p:Verrucomicrobia","p:Gemmatimonadetes" ,"p:Fungi","p:unclassified","p:", "unclassified")               
#who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Deinococcus-Thermus", "p:Nitrospirae","p:Spirochaetes","p:Verrucomicrobia","p:Gemmatimonadetes","p:Fibrobacteres","p:Planctomycetes","p:Fungi", "p:unclassified","p:","unclassified" )

all_samples_P$variable = ordered(all_samples_P$variable, who)
col = c("p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4","p:Deinococcus_Thermus"="#c8c9cc","p:Nitrospirae"="#7a6845", "p:Spirochaetes"="#197029","p:Verrucomicrobia"="#1e26bd","p:Gemmatimonadetes"="#b8ae54" ,"p:Fungi"="#ebd936", "p:unclassified"="#76bed6","p:"="#a3a3a3","unclassified"="#a3a3a3")
#col = c("p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4","p:Deinococcus_Thermus"="#c8c9cc","p:Nitrospirae"="#7a6845", "p:Spirochaetes"="#197029","p:Verrucomicrobia"="#1e26bd","p:Gemmatimonadetes"="#b8ae54" ,"p:Fibrobacteres"="#a969b8","p:Planctomycetes"="#61fa4d","p:Fungi"="#ebd936", "p:unclassified"="#76bed6","p:"="#a3a3a3","unclassified"="#a3a3a3")

metadata = read.csv("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/meta.csv", sep = ",", header = T)
all_samples_P = all_samples_P[all_samples_P$sample %in% metadata$Sample.name,]
all_samples_P$collection = metadata[ match(all_samples_P$sample, metadata$Sample.name), "Type" ]
all_samples_P$library = metadata[ match(all_samples_P$sample, metadata$Sample.name ), "Library" ]
all_samples_P$target = metadata[ match(all_samples_P$sample, metadata$Sample.name ), "Target" ]
all_samples_P$collection = ordered(all_samples_P$collection, c("mammalian skin","stool", "keyboard","blank"))

write.csv(all_samples_P, "all_samples_P_otu.csv")
p = ggplot(all_samples_P, aes(x = sample, y = value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20)) + theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

#sort the stacks by Actinobacteria
sample_order = all_samples[all_samples$P == "p:Actinobacteria",c(-1,-113)]
sample_order = sample_order %>% group_by(P)%>% summarise_all(list(sum))
sample_order$P = NULL
sample_order = as.data.frame(t(sample_order))
sample_order = arrange(sample_order, V1)
sample_order = row.names(sample_order)
all_samples_P$sample = ordered(all_samples_P$sample,rev(sample_order) )

p = ggplot(all_samples_P, aes(x = sample, y = value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20)) + theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p$data$sample = factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))

p # export in 8x3

#p = p + theme(panel.spacing = unit(0.01, "lines"))
#library(grid)
#gt = ggplot_gtable(ggplot_build(p))
#gt$widths[8] = 20*gt$widths[8]
#gt$widths[12] = 60*gt$widths[12]
#gt$widths[16] = 20*gt$widths[16]
#gt$widths[20] = 60*gt$widths[20]
#grid.draw(gt)


###--------------------compare with our direct output cpn60——-------------------###
metadata = read.csv("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/meta.csv", sep = ",", header = T)
tax = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Clade_checker/all_samples_taxonomy_closed_reference_90_newDB_all.tsv")
tax$`# OTU` = gsub(",s:.*", "", tax$`# OTU`)
tax$taxonomy = NULL
names(tax)[1] = "taxonomy"
tax = tax %>% group_by(taxonomy) %>% summarise_all(list(sum))
tax = tax %>% separate(taxonomy, sep = ";",c("Taxa", "taxonomy"))
tax = tax %>% separate(taxonomy, sep = ",",c("kindom", "phylum", "class", "order", "family", "genus"))
tax = tax[,c(3,8:118)]
tax = tax %>% group_by(phylum) %>% summarise_all(list(sum))
row.names(tax) = tax$phylum
tax = as.data.frame(t(tax))
tax = tax[2:112,]
row.names(tax) = gsub("_cpn_R1_001", "", row.names(tax))

read_count = fread("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Clade_checker/all_samples_read_counts_90_newDB_all.tsv")
read_count$total_reads = rowSums(read_count[,3:4])
read_count$`# sample` = gsub("_cpn_R1_001", "", read_count$`# sample` )
tax$reads = read_count[match(row.names(tax), read_count$`# sample`), "total_reads"]
tax$unclassified = read_count[match(row.names(tax), read_count$`# sample`), "reads mapping to unclassifed OTU"]
tax$reads = tax$reads$total_reads
tax$unclassified = tax$unclassified$`reads mapping to unclassifed OTU`
tax[] <- lapply(tax, function(x) as.numeric(as.character(x)))
tax = tax[order(colMeans(tax), decreasing = T)]
who = names(tax)
tax$Other = tax$reads - rowSums(tax[,c(2:ncol(tax))]) 
tax = data.frame(t(tax))
head(tax)
tax$Taxa <- row.names(tax)
m = melt(tax)
who = c("reads", "p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria", "p:Deinococcus_Thermus", "p:Basidiomycota","p:Chloroflexi",
        "p:Fusobacteria","p:Ascomycota","Other","unclassified")            
col = c("reads" = "#111211","p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4" , "p:Deinococcus_Thermus"="#c8c9cc","p:Basidiomycota" = "#8b7bbd","p:Chloroflexi"="#e3cc96", "p:Fusobacteria"="#afbce3", "p:Ascomycota"="#106145", "Other"="#76bed6", "unclassified"="#a3a3a3" )

m = subset(m, !m$Taxa=="reads")
m$collection = metadata[ match(m$variable, metadata$Sample.name), "Type" ]
m$library = metadata[ match(m$variable, metadata$Sample.name ), "Library" ]
m$target = metadata[ match(m$variable, metadata$Sample.name ), "Target" ]
m$collection = ordered(m$collection, c("mammalian skin","stool", "keyboard","blank"))
m = m[m$variable %in% sample_order,]
m$variable = ordered(m$variable, rev(sample_order))

p = ggplot(m, aes(x = variable, y = value, fill = Taxa)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Number of reads") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20)) + ylim(0, 40000) + theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )

p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(sample_order))

###relative abundance
tax = tax %>% select(-Taxa)
tax = as.data.frame(t(tax))
tax$reads = NULL
tax_gen = tax/rowSums(tax)
rowSums(tax_gen)
head(tax_gen)
tax_gen = as.data.frame(t(tax_gen))
tax_gen$Taxa = row.names(tax_gen)
m = melt(tax_gen)
m = subset(m, !m$Taxa=="reads")
m$collection = metadata[ match(m$variable, metadata$Sample.name), "Type" ]
m$library = metadata[ match(m$variable, metadata$Sample.name ), "Library" ]
m$target = metadata[ match(m$variable, metadata$Sample.name ), "Target" ]
m$collection = ordered(m$collection, c("mammalian skin","stool", "keyboard","blank"))
m = m[m$variable %in% sample_order,]
m$variable = ordered(m$variable, rev(sample_order))
m = na.omit(m)
write.csv(m, "cpn_biobakery_all.csv")
p = ggplot(m, aes(x = variable, y = value, fill = Taxa)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))+ theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )

p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(sample_order))
p
#only our samples cpn60
m_1 = m[!m$collection == "mammalian skin",]
m_1 = na.omit(m_1)
m_1$variable = ordered(m_1$variable, c("S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20"))
p = ggplot(m_1, aes(x = variable, y = value, fill = Taxa)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection+variable, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))

p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = rev(who))
p
p = p + theme(panel.spacing = unit(0.01, "lines"))
library(grid)
gt = ggplot_gtable(ggplot_build(p))
gt$widths[8] = 20*gt$widths[8]
gt$widths[12] = 60*gt$widths[12]
gt$widths[16] = 20*gt$widths[16]
gt$widths[20] = 60*gt$widths[20]
grid.draw(gt) #export in 4.5 x 3


##vsearch 16S
tax = fread("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/S1-S10/all_samples_taxonomy_closed_reference_silva.tsv")
tax$taxonomy = gsub("; s__.*", "", tax$taxonomy)
tax = tax[,2:12] %>% group_by(taxonomy) %>% summarise_all(list(sum))
tax_levels = strsplit(as.character(tax$taxonomy), split = "; ", fixed = TRUE)
tax_levels = data.frame(t(sapply(tax_levels, function(x) x)))
row.names(tax_levels) = tax$taxonomy
names(tax_levels) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
tax_levels = data.frame(tax_levels)
tax_levels = subset(tax_levels, ! Domain == "k__")
tax_levels$Domain = as.factor(tax_levels$Domain)
tax_levels$Domain = droplevels(tax_levels$Domain)
tax_levels$Phylum = ifelse(tax_levels$Phylum == "p__", "p__Unclassified", as.character(tax_levels$Phylum))
tax = tax[tax$taxonomy %in% row.names(tax_levels), ]
tax = as.data.frame(tax)
tax$phylum = tax_levels[match(tax$taxonomy, row.names(tax_levels)), "Phylum"]
read_count = fread("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/S1-S10/all_samples_read_counts_silva.tsv")
read_count$sample = read_count$`# sample`
row.names(read_count) = read_count$sample
read_count$`# sample`=NULL

###Check absolute
tax_abs = as.data.frame(tax)
tax_abs$taxonomy = NULL
tax_abs = tax_abs %>% group_by(phylum) %>% summarise_all(list(sum))
tax_abs = as.data.frame(tax_abs)
row.names(tax_abs) = tax_abs$phylum
tax_abs$phylum = NULL
tax_abs = as.data.frame(t(tax_abs))
who = names(sort(colMeans(tax_abs), decreasing = TRUE))
f = tax_abs[,names(tax_abs) %in% who]
row.names(f)=gsub("3_", "", row.names(f))
read_count$total_reads = rowSums(read_count[,2:3])
read_count$sample = gsub("3_", "", read_count$sample)
f$reads = read_count[match(row.names(f), read_count$sample), "total_reads"]
f$unclassified = read_count[match(row.names(f), read_count$sample), "reads mapping to unclassifed OTU"]
f = f[order(colMeans(f), decreasing = T)]
who = c(who, "unclassified", "reads")
f = data.frame(t(f))
head(f)
f$Taxa <- row.names(f)

##Relative abundance
tax_gen = f %>% select(-Taxa)
tax_gen = as.data.frame(t(tax_gen)) %>% select(-reads)
tax_gen = tax_gen/rowSums(tax_gen)
rowSums(tax_gen)
head(tax_gen)
tax_gen = as.data.frame(t(tax_gen))
tax_gen$Taxa = row.names(tax_gen)
m = melt(tax_gen)
m$collection = metadata[ match(m$variable, metadata$Sample.name), "Type" ]
m$library = metadata[ match(m$variable, metadata$Sample.name ), "Library" ]
m$target = metadata[ match(m$variable, metadata$Sample.name ), "Target" ]
m$collection = ordered(m$collection, c("stool", "keyboard","blank"))
m$variable = ordered(m$variable, c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))
m = subset(m, !m$Taxa=="reads")
who = c("p__Actinobacteria", "p__Bacteroidetes","p__Firmicutes" ,"p__Proteobacteria","p__Cyanobacteria", "p__Fusobacteria" ,"p__Deinococcus-Thermus", "p__Epsilonbacteraeota",  "p__Acidobacteria","p__Verrucomicrobia" ,   
        "p__Planctomycetes", "p__Spirochaetes", "p__Gemmatimonadetes","p__Thaumarchaeota","unclassified" ,"reads"  )
col = c("p__Actinobacteria"="#e3857b" ,"p__Bacteroidetes" = "#a3946a" ,"p__Firmicutes" = "#5bc9b1","p__Proteobacteria"="#d689d4" , "p__Cyanobacteria" ="#f5dce8","p__Fusobacteria" ="#afbce3","p__Deinococcus-Thermus"="#c8c9cc", 
        "p__Thaumarchaeota" ="#76bed6","p__Gemmatimonadetes"= "#76bed6","p__Spirochaetes"="#76bed6","p__Planctomycetes"="#76bed6",
        "p__Verrucomicrobia"="#76bed6","p__Acidobacteria"="#76bed6","p__Epsilonbacteraeota"="#76bed6", "unclassified"="#a3a3a3")

p = ggplot(m, aes(x = variable, y = value, fill = Taxa)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection+variable, space = "free", scales = "free") + 
        ylab("Number of reads") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))

p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = rev(who))
p = p + theme(panel.spacing = unit(0.01, "lines"))
library(grid)
gt = ggplot_gtable(ggplot_build(p))
gt$widths[8] = 20*gt$widths[8]
gt$widths[12] = 60*gt$widths[12]
gt$widths[16] = 20*gt$widths[16]
gt$widths[20] = 60*gt$widths[20]
grid.draw(gt)


####sintax cpn60 S11-20
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/sintax/")
S11 = read.csv("S11.sintax", sep = "\t", header = F)
S11 = S11[,c("V1", "V2", "V4")] %>% separate(V2, c("V3", "taxa"), sep = "\\),c")
S11 = S11 %>% separate(V3, c("D", "P"), sep = "\\),")
S11 = S11 %>% separate(P, c("P", "P_value"), sep = "\\(")
S11$P  = ifelse(S11$P_value < 0.8, "p:unclassified", S11$P)
S11 = S11 %>% separate(D, c("D", "D_value"), sep = "\\(")
S11$D  = ifelse(S11$D_value < 0.8, "unknown", S11$D)
S11$P[grep("unknown", S11$D)] = "unknown"
S11_sintax = S11[,c("V1", "P")] %>% group_by(P) %>% summarise(n = n()) %>% mutate(sample = "S11")
S11_sintax = as.data.frame(S11_sintax)

#write a loop for all the samples
sample_names <- c("S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20")
sample_sintax_list <- list() # Create an empty list to store the results for each sample

# Loop through each sample and perform the operations
for (sample_name in sample_names) {
        file_name <- paste0(sample_name, ".sintax")  # Read the data for the current sample
        df <- read.csv(file_name, sep = "\t", header = FALSE)
        df <- df[, c("V1", "V2", "V4")] %>%    # Perform the data manipulation steps as shown in your updated code
                separate(V2, c("V3", "taxa"), sep = "\\),c") %>%
                separate(V3, c("D", "P"), sep = "\\),") %>%
                separate(P, c("P", "P_value"), sep = "\\(")
        
        df$P <- ifelse(df$P_value < 0.8, "p:unclassified", df$P)
        df <- df %>% separate(D, c("D", "D_value"), sep = "\\(")
        df$D <- ifelse(df$D_value < 0.8, "unknown", df$D)
        df$P[grep("unknown", df$D)] <- "unknown"
        sintax_df <- df[, c("V1", "P")] %>%
                group_by(P) %>%
                summarise(n = n()) %>%
                mutate(sample = sample_name) %>%
                as.data.frame()
        sample_sintax_list[[sample_name]] <- sintax_df    # Store the result for the current sample in the list
}

# Combine the results from all samples into a single dataframe
cpn_sintax <- do.call(rbind, sample_sintax_list)
total = cpn_sintax[, c("sample", "n")] %>% group_by(sample) %>% summarise_all(sum)
cpn_sintax$sum = total[match(cpn_sintax$sample, total$sample), "n"]
cpn_sintax$sum = cpn_sintax$sum$n
cpn_sintax$RA = cpn_sintax$n / cpn_sintax$sum
write.csv(cpn_sintax, "cpn_sintax_P_S11-20.csv")
#plot the stack plot
cpn_sintax$collection = metadata[ match(cpn_sintax$sample, metadata$Sample.name), "Type" ]
cpn_sintax$library = metadata[ match(cpn_sintax$sample, metadata$Sample.name ), "Library" ]
cpn_sintax$target = metadata[ match(cpn_sintax$sample, metadata$Sample.name ), "Target" ]
cpn_sintax$collection = ordered(cpn_sintax$collection, c("stool", "keyboard","blank"))
cpn_sintax$sample = ordered(cpn_sintax$sample, c("S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20"))

col = c("p:Actinobacteria"="#e3857b" ,"p:Bacteroidetes" = "#a3946a" ,"p:Firmicutes" = "#5bc9b1","p:Proteobacteria"="#d689d4" , "p:Cyanobacteria" ="#f5dce8","p:Ascomycota"="#106145",
        "p:Apicomplexa"="#b81a71","p:Synergistetes"="#bdaa99","p:Acidobacteria"="#e85115", "p:Lentisphaerae"="#455959","p:[Thermi]"="#8a49ab","p:Nematoda"="#ad7424","p:Deinococcus_Thermus"="#c8c9cc", 
        "p:Verrucomicrobia"="#1e26bd","p:Basidiomycota"="#8b7bbd","p:Chloroflexi"="#e3cc96","p:unclassified"="#76bed6","unknown"="#a3a3a3")
who = c("p:Actinobacteria","p:Bacteroidetes","p:Firmicutes","p:Proteobacteria","p:Cyanobacteria","p:Ascomycota","p:Apicomplexa","p:Synergistetes","p:Acidobacteria", "p:Lentisphaerae","p:[Thermi]","p:Nematoda","p:Deinococcus_Thermus", "p:Verrucomicrobia","p:Basidiomycota","p:Chloroflexi","p:unclassified","unknown")                                    

p = ggplot(cpn_sintax, aes(x = sample, y = RA, fill = P)) + geom_bar(stat="identity") +
        theme_bw(base_size = 8) + facet_grid(~collection+sample, space = "free", scales = "free") + scale_fill_manual(values = col)+
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))+ facet_grid(~collection, space = "free", scales = "free")
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p

####sintax on all samples
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/sintax/all_analysis/")
#write a loop for all the samples
sample_names <- c("S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                  "ERR6365578","ERR6365579","ERR6365580","ERR6365581","ERR6365582","ERR6365583","ERR6365584","ERR6365585","ERR6365586","ERR6365587","ERR6365588",
                  "ERR6365589","ERR6365590","ERR6365591","ERR6365592","ERR6365593","ERR6365594","ERR6365595","ERR6365596","ERR6365597","ERR6365598","ERR6365599",
                  "ERR6365600","ERR6365601","ERR6365602","ERR6365603","ERR6365604","ERR6365605","ERR6365606","ERR6365607","ERR6365608","ERR6365609","ERR6365610",
                  "ERR6365611","ERR6365612","ERR6365613","ERR6365614","ERR6365615","ERR6365616","ERR6365617","ERR6365618","ERR6365619","ERR6365620","ERR6365621",
                  "ERR6365622","ERR6365623","ERR6365624","ERR6365625","ERR6365626","ERR6365627","ERR6365628","ERR6365629","ERR6365630","ERR6365631","ERR6365632",
                  "ERR6365633","ERR6365634","ERR6365635","ERR6365636","ERR6365637","ERR6365638","ERR6365639","ERR6365640","ERR6365641","ERR6365642","ERR6365643",
                  "ERR6365644","ERR6365645","ERR6365646","ERR6365647","ERR6365648","ERR6365649","ERR6365650","ERR6365651","ERR6365652","ERR6365653","ERR6365654",
                  "ERR6365655","ERR6365656","ERR6365657","ERR6365658","ERR6365659","ERR6365660","ERR6365661","ERR6365662","ERR6365663","ERR6365664","ERR6365665",
                  "ERR6365666","ERR6365667","ERR6365668")
sample_sintax_list <- list() # Create an empty list to store the results for each sample

# Loop through each sample and perform the operations
for (sample_name in sample_names) {
        file_name <- paste0(sample_name, ".sintax")  # Read the data for the current sample
        df <- read.csv(file_name, sep = "\t", header = FALSE)
        df <- df[, c("V1", "V2", "V4")] %>%    # Perform the data manipulation steps as shown in your updated code
                separate(V2, c("V3", "taxa"), sep = "\\),c") %>%
                separate(V3, c("D", "P"), sep = "\\),") %>%
                separate(P, c("P", "P_value"), sep = "\\(")
        
        df$P <- ifelse(df$P_value < 0.8, "p:unclassified", df$P)
        df <- df %>% separate(D, c("D", "D_value"), sep = "\\(")
        df$D <- ifelse(df$D_value < 0.8, "unknown", df$D)
        df$P[grep("unknown", df$D)] <- "unknown"
        sintax_df <- df[, c("V1", "P")] %>%
                group_by(P) %>%
                summarise(n = n()) %>%
                mutate(sample = sample_name) %>%
                as.data.frame()
        sample_sintax_list[[sample_name]] <- sintax_df    # Store the result for the current sample in the list
}

# Combine the results from all samples into a single dataframe
cpn_sintax <- do.call(rbind, sample_sintax_list)
total = cpn_sintax[, c("sample", "n")] %>% group_by(sample) %>% summarise_all(sum)
cpn_sintax$sum = total[match(cpn_sintax$sample, total$sample), "n"]
cpn_sintax$sum = cpn_sintax$sum$n
cpn_sintax$RA = cpn_sintax$n / cpn_sintax$sum
#plot the stack plot
cpn_sintax$collection = metadata[ match(cpn_sintax$sample, metadata$Sample.name), "Type" ]
cpn_sintax$library = metadata[ match(cpn_sintax$sample, metadata$Sample.name ), "Library" ]
cpn_sintax$target = metadata[ match(cpn_sintax$sample, metadata$Sample.name ), "Target" ]
cpn_sintax$collection = ordered(cpn_sintax$collection, c("mammalian skin","stool", "keyboard","blank"))
cpn_sintax$sample = ordered(cpn_sintax$sample, rev(sample_order))

col = c("p:Actinobacteria"="#e3857b" ,"p:Bacteroidetes" = "#a3946a" ,"p:Firmicutes" = "#5bc9b1","p:Proteobacteria"="#d689d4" , "p:Cyanobacteria" ="#f5dce8","p:Ascomycota"="#106145",
        "p:Apicomplexa"="#b81a71","p:Synergistetes"="#bdaa99","p:Acidobacteria"="#e85115", "p:Lentisphaerae"="#455959","p:[Thermi]"="#8a49ab","p:Nematoda"="#ad7424","p:Deinococcus_Thermus"="#c8c9cc", 
        "p:Verrucomicrobia"="#1e26bd","p:Basidiomycota"="#8b7bbd","p:Chloroflexi"="#e3cc96","p:unclassified"="#76bed6","unknown"="#a3a3a3")

who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Ascomycota","p:Apicomplexa","p:Synergistetes","p:Acidobacteria", "p:Lentisphaerae","p:[Thermi]","p:Nematoda","p:Deinococcus_Thermus","p:Verrucomicrobia","p:Basidiomycota", "p:Chloroflexi",                       
 "p:Spirochaetes","p:Nitrospirae" , "p:Gemmatimonadetes","p:Fibrobacteres", "p:Chordata" , "p:Abditibacteriota" , "p:Fusobacteria" ,"p:Candidatus_Tectomicrobia", "p:Elusimicrobia",  "p:Planctomycetes","p:Chlamydiae" , "p:Candidatus_Cloacimonetes","p:Candidatus_Riflebacteria","p:unclassified" ,"unknown")
cpn_sintax = na.omit(cpn_sintax)
p = ggplot(cpn_sintax, aes(x = sample, y = RA, fill = P)) + geom_bar(stat="identity") +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + scale_fill_manual(values = col)+
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))+ facet_grid(~collection, space = "free", scales = "free")+ theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$sample= factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))
p
write.csv(cpn_sintax, "cpn_sintax_P.csv")


###Check on the original ASV output from PMID: 37419988
asv = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/cpn60_ASV_table.tsv")
asv = asv[, 2:39]
asv$Consensus.Lineage = gsub(";D_2.*", "", asv$Consensus.Lineage)
asv = asv %>% group_by(Consensus.Lineage) %>% summarise_all(sum)
asv = as.data.frame(asv)
row.names(asv) = asv$Consensus.Lineage
asv$Consensus.Lineage = NULL
asv = as.data.frame(t(asv))
asv = asv /  rowSums(asv)
asv$sample  = row.names(asv)
asv = melt(asv)
asv$sample = gsub("_B", "", asv$sample)
asv$variable = gsub("D_0__", "k:", asv$variable)
asv$variable = gsub(";D_1__", ",p:", asv$variable)
who = c("k:Bacteria,p:Actinobacteria","k:Bacteria,p:Firmicutes","k:Bacteria,p:Bacteroidetes","k:Bacteria,p:Cyanobacteria","k:Bacteria,p:Proteobacteria", "k:Bacteria,p:Acidobacteria",                          
        "k:Bacteria,p:Deinococcus-Thermus","k:Bacteria,p:CandidatusSaccharibacteria","k:Bacteria,p:Spirochaetes","k:Bacteria,p:Verrucomicrobia","k:Bacteria,p:Chloroflexi","k:Bacteria,p:Elusimicrobia",             
      "k:Bacteria,p:Fibrobacteres", "k:Bacteria,p:Gemmatimonadetes","k:Bacteria,p:Kiritimatiellaeota","k:Bacteria,p:Nitrospirae" ,"k:Bacteria,p:Planctomycetes",  "k:Bacteria" )

col = c("k:Bacteria,p:Actinobacteria"="#e3857b","k:Bacteria,p:Firmicutes"="#5bc9b1","k:Bacteria,p:Bacteroidetes"="#a3946a","k:Bacteria,p:Proteobacteria"="#d689d4","k:Bacteria,p:Cyanobacteria"="#f5dce8", "k:Bacteria,p:Acidobacteria"="#e85115",                          
        "k:Bacteria,p:Deinococcus-Thermus"="#c8c9cc","k:Bacteria,p:CandidatusSaccharibacteria"="#f21d88","k:Bacteria,p:Spirochaetes"="#197029","k:Bacteria,p:Verrucomicrobia"="#1e26bd","k:Bacteria,p:Chloroflexi"="#e3cc96","k:Bacteria,p:Elusimicrobia"="#f77111",             
        "k:Bacteria,p:Fibrobacteres"="#a969b8", "k:Bacteria,p:Gemmatimonadetes"="#b8ae54","k:Bacteria,p:Kiritimatiellaeota"="#fafa4d","k:Bacteria,p:Nitrospirae"="#7a6845" ,"k:Bacteria,p:Planctomycetes"="#61fa4d",  "k:Bacteria"="#76bed6" )
metadata = read.csv("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/meta.csv", sep = ",", header = T)
asv$sample = metadata[match(asv$sample, metadata$Source), "Sample.name"]
asv$collection = metadata[ match(asv$sample, metadata$Sample.name), "Type" ]
asv$library = metadata[ match(asv$sample, metadata$Sample.name ), "Library" ]
asv$target = metadata[ match(asv$sample, metadata$Sample.name ), "Target" ]
asv$collection = ordered(asv$collection, c("mammalian skin","stool", "keyboard","blank"))
asv$sample = ordered(asv$sample, rev(sample_order))

p = ggplot(asv, aes(x = sample, y = value, fill = variable)) + geom_bar(stat="identity") +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + scale_fill_manual(values = col)+
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))+ facet_grid(~collection, space = "free", scales = "free")+ theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p$data$sample= factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))
p

#only plot the 37 skin samples in the all_analysis sintax
skin = cpn_sintax[cpn_sintax$collection == "mammalian skin",]
skin = skin[skin$sample %in% asv$sample,]
cpn_sintax = rbind(skin, cpn_sintax[!cpn_sintax$collection =="mammalian skin",])
cpn_sintax$sample = ordered(cpn_sintax$sample, c(sample_order))
col = c("p:Actinobacteria"="#e3857b" ,"p:Bacteroidetes" = "#a3946a" ,"p:Firmicutes" = "#5bc9b1","p:Proteobacteria"="#d689d4" , "p:Cyanobacteria" ="#f5dce8","p:Ascomycota"="#106145",
        "p:Apicomplexa"="#b81a71","p:Synergistetes"="#bdaa99","p:Acidobacteria"="#e85115", "p:Lentisphaerae"="#455959","p:[Thermi]"="#8a49ab","p:Nematoda"="#ad7424","p:Deinococcus_Thermus"="#c8c9cc", 
        "p:Verrucomicrobia"="#1e26bd","p:Basidiomycota"="#8b7bbd","p:Chloroflexi"="#e3cc96","p:unclassified"="#76bed6","unknown"="#a3a3a3")

who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Ascomycota","p:Apicomplexa","p:Synergistetes","p:Acidobacteria", "p:Lentisphaerae","p:[Thermi]","p:Nematoda","p:Deinococcus_Thermus","p:Verrucomicrobia","p:Basidiomycota", "p:Chloroflexi",                       
        "p:Spirochaetes","p:Nitrospirae" , "p:Gemmatimonadetes","p:Fibrobacteres", "p:Chordata" , "p:Abditibacteriota" , "p:Fusobacteria" ,"p:Candidatus_Tectomicrobia", "p:Elusimicrobia",  "p:Planctomycetes","p:Chlamydiae" , "p:Candidatus_Cloacimonetes","p:Candidatus_Riflebacteria","p:unclassified" ,"unknown")
p = ggplot(cpn_sintax, aes(x = sample, y = RA, fill = P)) + geom_bar(stat="identity") +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + scale_fill_manual(values = col)+
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))+ facet_grid(~collection, space = "free", scales = "free")+ theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$sample= factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))
p

#only plot the 37 skin samples in the all_analysis RDP classifier (otu)
all_samples_P = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/all_samples_P_otu.csv")
who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Deinococcus-Thermus", "p:Nitrospirae","p:Spirochaetes","p:Verrucomicrobia","p:Gemmatimonadetes" ,"p:Fungi","p:unclassified","p:", "unclassified")               
col = c("p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4","p:Deinococcus_Thermus"="#c8c9cc","p:Nitrospirae"="#7a6845", "p:Spirochaetes"="#197029","p:Verrucomicrobia"="#1e26bd","p:Gemmatimonadetes"="#b8ae54" ,"p:Fungi"="#ebd936", "p:unclassified"="#76bed6","p:"="#a3a3a3","unclassified"="#a3a3a3")
all_samples_P$variable = ordered(all_samples_P$variable, who)
all_samples_P$sample = ordered(all_samples_P$sample, rev(sample_order))
all_samples_P$collection = ordered(all_samples_P$collection, c("mammalian skin","stool", "keyboard","blank"))
all_samples_P = all_samples_P[all_samples_P$sample %in%cpn_sintax$sample,]
p = ggplot(all_samples_P, aes(x = sample, y = value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20)) + theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p$data$sample = factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))
p # export in 8x3

#only plot the 37 skin samples in the total biobakery workflow
m = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/sintax/all_analysis/cpn_biobakery_all.csv")
who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria", "p:Deinococcus_Thermus", "p:Basidiomycota","p:Chloroflexi",
        "p:Fusobacteria","p:Ascomycota","Other","unclassified")            
col = c("p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4" , "p:Deinococcus_Thermus"="#c8c9cc","p:Basidiomycota" = "#8b7bbd","p:Chloroflexi"="#e3cc96", "p:Fusobacteria"="#afbce3", "p:Ascomycota"="#106145", "Other"="#76bed6", "unclassified"="#a3a3a3" )
m = m[m$variable %in% cpn_sintax$sample,]
m$collection = ordered(m$collection, c("mammalian skin", "stool", "keyboard", "blank"))
p = ggplot(m, aes(x = variable, y = value, fill = Taxa)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20))+ theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )

p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(sample_order))
p # export in 6x3


##use the ASV from the paper (PMID: 37419988), re-classified using the RDP classifier
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/")
asv = read.csv("original/asv_RDPClassifier.csv", sep = ",", header = F)
asv = asv[, c(1,6,8,9,11,12,14,15,17,18,20,21,23,24,26)]
names(asv) = c("asv", "K", "K_value", "P", "P_value", "C", "C_value", "O", "O_value", "F", "F_value", "G", "G_value", "S", "S_value")
asv$K  = ifelse(asv$K_value < 0.8, "", asv$K)
asv$P  = ifelse(asv$P_value < 0.8, "", asv$P)
asv$C  = ifelse(asv$C_value < 0.8, "", asv$C)
asv$O  = ifelse(asv$O_value < 0.8, "", asv$O)
asv$`F`  = ifelse(asv$F_value < 0.8, "", asv$`F`)
asv$G  = ifelse(asv$G_value < 0.8, "", asv$G)
asv$S  = ifelse(asv$S_value < 0.8, "", asv$S)

asv$K = paste("k",asv$K, sep = ":" )
asv$P = paste("p",asv$P, sep = ":" )
asv$C = paste("c",asv$C, sep = ":" )
asv$O = paste("o",asv$O, sep = ":" )
asv$`F` = paste("f",asv$`F`, sep = ":" )
asv$G = paste("g",asv$G, sep = ":" )
asv$S = paste("s",asv$S, sep = ":" )

asv$taxonomy = paste(asv$K, asv$P, asv$C, asv$O, asv$`F`, asv$G, asv$S, sep = ",")
asv = asv %>% separate(asv, sep = ";", c("asv", "seqs"))

all_samples = read.csv("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/cpn60_ASV_table.tsv", sep = "\t", header = T)
names(all_samples)[1] = c("asv")
names(all_samples) = gsub("_B", "", names(all_samples))
all_samples$asv = paste("ASV", row.names(all_samples))
all_samples$asv = gsub("ASV ", "ASV", all_samples$asv)
all_samples$taxonomy = asv[match(all_samples$asv, asv$asv), "taxonomy"]
all_samples$P = asv[match(all_samples$asv, asv$asv ), "P"]
all_samples$P  = ifelse(all_samples$taxonomy=="k:,p:,c:,o:,f:,g:,s:", "unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:,p:,c:,o:,f:,g:,s:NA", "unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:,p:,c:,o:,f:,g:NA,s:NA", "unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:", "p:unclassified", all_samples$P)
all_samples$P  = ifelse(all_samples$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:NA", "p:unclassified", all_samples$P)

asv_taxonomy = all_samples[, c(2:38, 41)]
write.csv(asv_taxonomy, "asv_taxonomy.csv")

all_samples_P = all_samples[, c(2:38,42)]
all_samples_P = all_samples_P %>% group_by(P)%>% summarise_all(list(sum))
all_samples_P = as.data.frame(all_samples_P)
row.names(all_samples_P) = all_samples_P$P
all_samples_P$P = NULL
all_samples_P = as.data.frame(t(all_samples_P))
all_samples_P = all_samples_P / rowSums(all_samples_P)
all_samples_P$sample = row.names(all_samples_P)
all_samples_P = melt(all_samples_P)
who = c("p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Deinococcus-Thermus", "p:Nitrospirae","p:Spirochaetes","p:Verrucomicrobia","p:Gemmatimonadetes","p:Fibrobacteres","p:Planctomycetes","p:Fungi", "p:unclassified","p:","unclassified" )
all_samples_P$variable = ordered(all_samples_P$variable, who)
col = c("p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4","p:Deinococcus_Thermus"="#c8c9cc","p:Nitrospirae"="#7a6845", "p:Spirochaetes"="#197029","p:Verrucomicrobia"="#1e26bd","p:Gemmatimonadetes"="#b8ae54" ,"p:Fibrobacteres"="#a969b8","p:Planctomycetes"="#61fa4d","p:Fungi"="#ebd936", "p:unclassified"="#76bed6","p:"="#a3a3a3","unclassified"="#a3a3a3")

metadata = read.csv("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/meta.csv", sep = ",", header = T)
all_samples_P$sample = metadata[match(all_samples_P$sample, metadata$Source), "Sample.name"]
all_samples_P = all_samples_P[all_samples_P$sample %in% metadata$Sample.name,]
all_samples_P$collection = metadata[ match(all_samples_P$sample, metadata$Sample.name), "Type" ]
all_samples_P$library = metadata[ match(all_samples_P$sample, metadata$Sample.name ), "Library" ]
all_samples_P$target = metadata[ match(all_samples_P$sample, metadata$Sample.name ), "Target" ]
all_samples_P$collection = ordered(all_samples_P$collection, c("mammalian skin","stool", "keyboard","blank"))

write.csv(all_samples_P, "all_samples_P_asv.csv")

p = ggplot(all_samples_P, aes(x = sample, y = value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20)) + theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p$data$sample = factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))

p # export in 8x3

##what if we don't filter the 0.8, instead we just take everything identified here
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/")
asv = read.csv("original/asv_RDPClassifier.csv", sep = ",", header = F)
asv = asv[, c(1,6,8,9,11,12,14,15,17,18,20,21,23,24,26)]
names(asv) = c("asv", "K", "K_value", "P", "P_value", "C", "C_value", "O", "O_value", "F", "F_value", "G", "G_value", "S", "S_value")
#asv$K  = ifelse(asv$K_value < 0.8, "", asv$K)
#asv$P  = ifelse(asv$P_value < 0.8, "", asv$P)
asv$K = paste("k",asv$K, sep = ":" )
asv$P = paste("p",asv$P, sep = ":" )
asv$taxonomy = paste(asv$K, asv$P, sep = ",")

all_samples = read.csv("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/cpn60_ASV_table.tsv", sep = "\t", header = T)
names(all_samples)[1] = c("asv")
names(all_samples) = gsub("_B", "", names(all_samples))
all_samples$asv = paste("ASV", row.names(all_samples))
all_samples$asv = gsub("ASV ", "ASV", all_samples$asv)
all_samples$P = asv[match(all_samples$asv, asv$asv ), "P"]
all_samples$taxonomy = asv[match(all_samples$asv, asv$asv), "taxonomy"]

all_samples_P = all_samples[, c(2:38,41)]
all_samples_P = all_samples_P %>% group_by(P)%>% summarise_all(list(sum))
all_samples_P = as.data.frame(all_samples_P)
row.names(all_samples_P) = all_samples_P$P
all_samples_P$P = NULL
all_samples_P = as.data.frame(t(all_samples_P))
all_samples_P = all_samples_P / rowSums(all_samples_P)
all_samples_P$sample = row.names(all_samples_P)
all_samples_P = reshape2::melt(all_samples_P)
who = c( "p:Actinobacteria","p:Firmicutes","p:Bacteroidetes","p:Cyanobacteria","p:Proteobacteria","p:Acidobacteria" ,
         "p:Deinococcus-Thermus","p:Abditibacteriota","p:Apicomplexa","p:Aquificae","p:Armatimonadetes","p:Bacillariophyta","p:Balneolaeota",            
        "p:Candidatus Cloacimonetes","p:Candidatus Tectomicrobia", "p:Chlamydiae","p:Chlorobi","p:Chloroflexi","p:Choanoflagellata","p:Deferribacteres",               
        "p:Dictyoglomi","p:Dinophyceae","p:Elusimicrobia","p:Euglenozoa","p:Euryarchaeota","p:Fibrobacteres","p:Fungi","p:Fusobacteria",           
        "p:Gemmatimonadetes","p:Haptista","p:Ignavibacteriae","p:Lentisphaerae","p:Metazoa","p:Nitrospinae","p:Nitrospirae","p:Pelagophyceae","p:Phaeophyceae",            
        "p:Planctomycetes","p:Rhodophyta","p:Rhodothermaeota","p:Spirochaetes","p:Synergistetes","p:Tenericutes","p:Thermotogae","p:Verrucomicrobia","p:Viridiplantae"  )
all_samples_P$variable = ordered(all_samples_P$variable, who)
#col = c("p:Actinobacteria"="#e3857b" ,"p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#a3946a" ,"p:Cyanobacteria" ="#f5dce8","p:Proteobacteria"="#d689d4","p:Deinococcus_Thermus"="#c8c9cc","p:Nitrospirae"="#7a6845", "p:Spirochaetes"="#197029","p:Verrucomicrobia"="#1e26bd","p:Gemmatimonadetes"="#b8ae54" ,"p:Fibrobacteres"="#a969b8","p:Planctomycetes"="#61fa4d","p:Fungi"="#ebd936", "p:unclassified"="#76bed6","p:"="#a3a3a3","unclassified"="#a3a3a3")
col = c("p:Actinobacteria"="#e3857b","p:Firmicutes"="#5bc9b1","p:Bacteroidetes"="#a3946a","p:Proteobacteria"="#d689d4","p:Cyanobacteria"="#f5dce8", "p:Acidobacteria"="#e85115",                          
        "p:Deinococcus-Thermus"="#c8c9cc","p:CandidatusSaccharibacteria"="#f21d88","p:Spirochaetes"="#197029","p:Verrucomicrobia"="#1e26bd","p:Chloroflexi"="#e3cc96","p:Elusimicrobia"="#f77111",             
        "p:Fibrobacteres"="#a969b8", "p:Gemmatimonadetes"="#b8ae54","p:Kiritimatiellaeota"="#fafa4d","p:Nitrospirae"="#7a6845" ,"p:Planctomycetes"="#61fa4d",  "k:Bacteria"="#76bed6" )

metadata = read.csv("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/meta.csv", sep = ",", header = T)
all_samples_P$sample = metadata[match(all_samples_P$sample, metadata$Source), "Sample.name"]
all_samples_P = all_samples_P[all_samples_P$sample %in% metadata$Sample.name,]
all_samples_P$collection = metadata[ match(all_samples_P$sample, metadata$Sample.name), "Type" ]
all_samples_P$library = metadata[ match(all_samples_P$sample, metadata$Sample.name ), "Library" ]
all_samples_P$target = metadata[ match(all_samples_P$sample, metadata$Sample.name ), "Target" ]
all_samples_P$collection = ordered(all_samples_P$collection, c("mammalian skin","stool", "keyboard","blank"))

p = ggplot(all_samples_P, aes(x = sample, y = value, fill = variable)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) + facet_grid(~collection, space = "free", scales = "free") + 
        ylab("Relative abundance") + theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) + 
        labs(title="") + theme(plot.title = element_text(size = 20)) + theme(axis.text.x = element_blank(), axis.ticks.x =element_blank() )
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p$data$sample = factor(p$data$sample, ordered = TRUE, levels = rev(sample_order))

p # export in 8x3


####a dotplot showing that there's no big difference between otu vs. asv
otu_taxonomy = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/otu_taxonomy.csv")
asv_taxonomy = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/asv_taxonomy.csv")
otu_taxonomy = as.data.frame(otu_taxonomy)
asv_taxonomy = as.data.frame(asv_taxonomy)
names(asv_taxonomy) = metadata[match(names(asv_taxonomy), metadata$Source), "Sample.name"]
who = names(asv_taxonomy) 
otu_taxonomy = otu_taxonomy[, c(who)]

asv_taxonomy$V1 = NULL
otu_taxonomy$V1 = NULL

asv = melt(asv_taxonomy)
asv$temp = asv$taxonomy

asv$temp  = ifelse(asv$taxonomy=="k:,p:,c:,o:,f:,g:,s:", "unclassified", asv$temp)
asv$temp  = ifelse(asv$taxonomy=="k:,p:,c:,o:,f:,g:,s:NA", "unclassified", asv$temp)
asv$temp  = ifelse(asv$taxonomy=="k:,p:,c:,o:,f:,g:NA,s:NA", "unclassified", asv$temp)
asv$temp  = ifelse(asv$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:", "k:Bacteria,p:unclassified", asv$temp)
asv$temp  = ifelse(asv$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:NA", "k:Bacteria,p:unclassified", asv$temp)
asv$temp  = ifelse(asv$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:NA,s:NA", "k:Bacteria,p:unclassified", asv$temp)
asv$temp  = ifelse(asv$taxonomy=="k:Bacteria,p:,c:,o:,f:NA,g:NA,s:NA", "k:Bacteria,p:unclassified", asv$temp)

asv$temp = gsub("NA", "", asv$temp)
asv$temp = gsub(",c:,o:,f:,g:,s:", "", asv$temp)
asv$temp = gsub(",o:,f:,g:,s:", "", asv$temp)
asv$temp = gsub(",f:,g:,s:", "", asv$temp)
asv$temp = gsub(",g:,s:", "", asv$temp)
asv$temp = gsub(",c:,o:,f:,g:,s:", "", asv$temp)

asv = asv[, c("variable", "temp", "value")] %>% group_by(variable, temp) %>% summarise_all(sum)
asv_sum = asv %>% group_by(variable) %>% summarise( total = sum(value))
asv$sum = asv_sum[match(asv$variable, asv_sum$variable), "total"]
asv$total = asv$sum$total
asv$RA = asv$value / asv$total

otu = melt(otu_taxonomy)
otu$temp = otu$taxonomy
otu$temp  = ifelse(otu$taxonomy=="k:,p:,c:,o:,f:,g:,s:", "unclassified", otu$temp)
otu$temp  = ifelse(otu$taxonomy=="k:,p:,c:,o:,f:,g:,s:NA", "unclassified", otu$temp)
otu$temp  = ifelse(otu$taxonomy=="k:,p:,c:,o:,f:,g:NA,s:NA", "unclassified", otu$temp)
otu$temp  = ifelse(otu$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:", "k:Bacteria,p:unclassified", otu$temp)
otu$temp  = ifelse(otu$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:,s:NA", "k:Bacteria,p:unclassified", otu$temp)
otu$temp  = ifelse(otu$taxonomy=="k:Bacteria,p:,c:,o:,f:,g:NA,s:NA", "k:Bacteria,p:unclassified", otu$temp)
otu$temp  = ifelse(otu$taxonomy=="k:Bacteria,p:,c:,o:,f:NA,g:NA,s:NA", "k:Bacteria,p:unclassified", otu$temp)

otu$temp = gsub("NA", "", otu$temp)
otu$temp = gsub(",c:,o:,f:,g:,s:", "", otu$temp)
otu$temp = gsub(",o:,f:,g:,s:", "", otu$temp)
otu$temp = gsub(",f:,g:,s:", "", otu$temp)
otu$temp = gsub(",g:,s:", "", otu$temp)
otu$temp = gsub(",c:,o:,f:,g:,s:", "", otu$temp)

otu = otu[, c("variable", "temp", "value")] %>% group_by(variable, temp) %>% summarise_all(sum)
otu_sum = otu %>% group_by(variable) %>% summarise( total = sum(value))
otu$sum = otu_sum[match(otu$variable, otu_sum$variable), "total"]
otu$total = otu$sum$total
otu$RA = otu$value / otu$total

asv$sort = "asv"
asv$sample_T = paste(asv$variable, asv$temp, sep = "_")
otu$sort = "otu"
otu$sample_T = paste(otu$variable, otu$temp, sep = "_")

total = otu[, c("variable", "temp", "RA", "sample_T")]
total$asv_RA = asv[match(total$sample_T, asv$sample_T), "RA"]
total$asv_RA = total$asv_RA$RA
names(total)[3] = "otu_RA"
total = na.omit(total)
p = ggplot(total, aes(x=otu_RA, y=asv_RA)) + geom_point(size = 2, color = "blue", alpha=0.5) + scale_y_sqrt() + scale_x_sqrt() + theme_classic()+
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black")+
        labs(title="terminus taxa in every sample",
             x="UPARSE OTU", y = "DADA2 ASV")
p

total = rbind(asv[, c("sample_T", "RA", "sort")], otu[, c("sample_T", "RA", "sort")])
total = total %>% pivot_wider(names_from = sort, values_from = RA,  values_fill = list(n = 0))
total[is.na(total)] <- 0
total$check = total$asv +total$otu
total = total[!total$check == 0, ]
p = ggplot(total, aes(x=otu, y=asv)) + 
        geom_point(size = 2, color = "blue", alpha=0.5) + 
        scale_y_sqrt() + 
        scale_x_sqrt() + 
        theme_classic() +
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = "black", aes(group = 1)) +  # Set group to ensure one line
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # 1:1 diagonal line
        labs(title="terminus taxa in every sample",
             x="UPARSE OTU", y = "DADA2 ASV")
p


##----------------------statistical comparison on the taxa identidied----------------------##########
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(data.table)
library(vegan)   
cpn_RDP = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/all_samples_P_otu.csv")
cpn_sintax = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/sintax/all_analysis/cpn_sintax_P.csv")
cpn_amplicon_sintax = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/RDPClassifier/sintax/all_analysis/cpn_sintax_P_amplicon.csv")

cpn_RDP = cpn_RDP[, c("variable", "sample","value", "collection", "library")]
names(cpn_RDP) = c("P", "sample", "RA", "collection", "library")

cpn_RDP$sort = "RDP"
cpn_sintax = cpn_sintax[, c("P", "sample", "RA", "collection", "library")]%>% mutate(sort = "sintax")
cpn_amplicon_sintax = cpn_amplicon_sintax[, c("P", "sample", "RA", "collection", "library")] %>% mutate(sort = "sintax_amplicon")

cpn_all = rbind(cpn_RDP, cpn_sintax, cpn_amplicon_sintax)
cpn_all$P[cpn_all$P %in% c("p:", "unknown")] <- "unclassified"

## Ensure RA is numeric
cpn_all$RA <- suppressWarnings(as.numeric(cpn_all$RA))

## 1) Collapse any duplicates within sample × method × phylum
collapsed <- cpn_all %>%
        group_by(sample, sort, P) %>%
        summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

## 2) Decide which method pairs to compare (only those present)
methods_present <- sort(unique(collapsed$sort))
method_pairs <- expand.grid(m1 = methods_present, m2 = methods_present, stringsAsFactors = FALSE) %>%
        dplyr::filter(m1 < m2) %>%
        dplyr::filter(m1 %in% c("RDP","sintax","sintax_amplicon"),
                      m2 %in% c("RDP","sintax","sintax_amplicon"))

## 3) Compute Bray–Curtis distance (and similarity) per sample × method-pair
bc_by_sample <- collapsed %>%
        group_by(sample) %>%
        group_split() %>%
        map_dfr(function(df_s) {
                wide <- df_s %>%
                        select(P, sort, RA) %>%
                        pivot_wider(
                                names_from  = sort,
                                values_from = RA,
                                values_fill = 0,           # zero-fill missing phyla
                                values_fn   = list(RA = sum)
                        ) %>%
                        arrange(P)
                
                map2_dfr(method_pairs$m1, method_pairs$m2, ~{
                        if (!(.x %in% names(wide)) || !(.y %in% names(wide))) {
                                return(tibble(sample = unique(df_s$sample), m1 = .x, m2 = .y,
                                              bray_curtis_distance = NA_real_, bray_curtis_similarity = NA_real_))
                        }
                        a <- as.numeric(wide[[.x]])
                        b <- as.numeric(wide[[.y]])
                        
                        # If your RA values are already relative abundances, no need to renormalize.
                        # (If not, uncomment the next two lines to normalize each column to sum 1)
                        # a <- if ((sa <- sum(a)) > 0) a/sa else a
                        # b <- if ((sb <- sum(b)) > 0) b/sb else b
                        
                        d <- as.numeric(vegdist(rbind(a, b), method = "bray"))
                        tibble(sample = unique(df_s$sample), m1 = .x, m2 = .y,
                               bray_curtis_distance = d,
                               bray_curtis_similarity = 1 - d)
                })
        })

## 4) Summarize across samples (median, IQR)
bc_summary <- bc_by_sample %>%
        group_by(m1, m2) %>%
        summarise(
                n_samples      = sum(!is.na(bray_curtis_distance)),
                bc_dist_median = median(bray_curtis_distance, na.rm = TRUE),
                bc_dist_IQR_l  = quantile(bray_curtis_distance, 0.25, na.rm = TRUE),
                bc_dist_IQR_u  = quantile(bray_curtis_distance, 0.75, na.rm = TRUE),
                bc_sim_median  = median(bray_curtis_similarity, na.rm = TRUE),
                bc_sim_IQR_l   = quantile(bray_curtis_similarity, 0.25, na.rm = TRUE),
                bc_sim_IQR_u   = quantile(bray_curtis_similarity, 0.75, na.rm = TRUE),
                .groups = "drop"
        )

print(bc_summary)