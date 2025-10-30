
library("data.table")
library("httr")
library("jsonlite")
library("tidyr")
library("seqinr")
library("Biostrings")


##check the uniref90 IDs
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/")
lca = fread("uniref90_lca.dat")
unirefs = fread("final_uniq.txt", header = FALSE)
mapped <- fread("uniref_seq_pangenomes_20210718.txt", sep=";", header = FALSE)
mapped = mapped %>% separate(V1, sep = "\\|", c("D", "Taxa", "uniref90", "uniref50"))
mapped = mapped %>% separate(Taxa, sep = "\\.", c("K", "P","C", "O", "F", "G", "S"))
mapped = mapped[, c("K", "P","C", "O", "F", "G", "uniref90")]
mapped = unique(mapped)
unirefs$P = mapped[match(unirefs$V1, mapped$uniref90),"P"]
unirefs$K = mapped[match(unirefs$V1, mapped$uniref90),"K"]
mapped = na.omit(mapped)
mapped$humann_lca = lca[match(mapped$uniref90, lca$V1), "V4"]

unmapped = unirefs[!unirefs$V1 %in% mapped$uniref90,]
unmapped$humann_lca = lca[match(unmapped$V1, lca$V1), "V4"]
unmapped = fread("unmapped_with_taxid.tsv")
tax = fread("uniref90-tol-lca.dat.bz2")
unmapped$S = tax[match(unmapped$humann_lca, tax$TAXID), "NAME"]
unfound = fread("unfound_all.csv")
unmapped$P = unfound[match(unmapped$V1, unfound$FAMILY), "phylum"]
unmapped$K = unfound[match(unmapped$V1, unfound$FAMILY), "Kingdom"]
names(unmapped)[1] = "uniref90"
mapped$sort = "mapped"
unmapped$sort = "unmapped"
all = rbind(mapped[,c("uniref90", "humann_lca", "K", "P", "sort")], unmapped[,c("uniref90", "humann_lca", "K", "P", "sort")])
write.table(unmapped, file = "unmapped_with_taxid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all, file = "all_cpn_uniref90_taxid.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


##Check the taxonomy of DNA sequences
#kindom level
id = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/cpn_gene_20556_utaxID.csv", sep = ",", header = T)
pie = data.frame(table(id$kindom))
pie$labels = pie$Freq / sum(pie$Freq)
is.num <- sapply(pie, is.numeric)
pie[is.num] <- lapply(pie[is.num], round, 2)
k = c("k__Archaea"="#ba829d" , "k__Bacteria" ="#7bb581", "k__Eukaryota"="#e3c65f") 
ggplot(pie, aes(x = "", y = labels, fill = Var1)) +
        geom_col(color = "black") +
        geom_label(aes(label = labels), color = c(1, 1, 1),
                   position = position_stack(vjust = 0.5),
                   show.legend = FALSE) +
        guides(fill = guide_legend(title = "Answer")) +
        scale_fill_manual(values = k) +
        coord_polar(theta = "y") 

#phylum level
pie = data.frame(table(id$P))
pie$labels = pie$Freq / sum(pie$Freq)

pie_sub = subset(pie, !pie$labels < 0.01)
row.names(pie_sub) = pie_sub$Var1
pie_sub$Var1 = NULL
pie_sub = as.data.frame(t(pie_sub))
row.names(pie) =pie$Var1
pie$Var1 = NULL
pie_sub$Others = colSums(pie)- rowSums(pie_sub)
is.num <- sapply(pie_sub, is.numeric)
pie_sub[is.num] <- lapply(pie_sub[is.num], round, 2)
pie_sub = as.data.frame(t(pie_sub))
pie_sub$kindom = id[match(row.names(pie_sub), id$P), "kindom"]
pie_sub$phylum =row.names(pie_sub)
pie_sub = pie_sub[order(pie_sub$labels, decreasing = T),]
pie_sub$phylum = ordered(pie_sub$phylum, c("p:Proteobacteria", "p:Actinobacteria", "p:Firmicutes", "p:Bacteroidetes",  "p:Cyanobacteria",  "p:Euryarchaeota", "p:Ascomycota", "Others" ))

who = pie_sub$phylum
col = c("p:Euryarchaeota"="#f7b2b0", "p:Actinobacteria"="#d1cb94", "p:Bacteroidetes" = "#f5b767","p:Cyanobacteria"="#abd690", "p:Firmicutes" = "#9abdf5", "p:Proteobacteria" = "#8abfb4",  "p:Ascomycota" = "#ba95cc", "Others" = "#613782")

p = ggplot(pie_sub, aes(x = "", y = labels, fill = pie_sub$phylum)) + geom_col(color = "black") +scale_fill_manual(values = col) + guides(fill = guide_legend(title = "phylum")) + coord_polar(theta = "y") +geom_label(aes(label = labels), position = position_stack(vjust = 0.5)) 
p$data$phylum= factor(p$data$phylum, ordered = TRUE, levels = rev(who))
p

#de-duplicated database
fasta <- readDNAStringSet(filepath = "~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Clade_checker/cpn_gene_23487DeDuplicate.fa", format="fasta")
id =as.data.frame(names(fasta))
names(id) = "id"
id$taxa = id$id
id = id %>% separate(id, c("kindom", 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ";")
#kindom level
pie = data.frame(table(id$kindom))
pie$labels = pie$Freq / sum(pie$Freq)
is.num <- sapply(pie, is.numeric)
pie[is.num] <- lapply(pie[is.num], round, 2)
k = c("k__Archaea"="#ba829d" , "k__Bacteria" ="#7bb581", "k__Eukaryota"="#e3c65f") 
ggplot(pie, aes(x = "", y = labels, fill = Var1)) +
        geom_col(color = "black") +
        geom_label(aes(label = labels), color = c(1, 1, 1), position = position_stack(vjust = 0.5), show.legend = FALSE) +
        guides(fill = guide_legend(title = "Answer")) +
        scale_fill_manual(values = k) +
        coord_polar(theta = "y") 

#phylum level
pie = data.frame(table(id$phylum))
pie$labels = pie$Freq / sum(pie$Freq)
pie_sub = subset(pie, !pie$labels < 0.01)
row.names(pie_sub) = pie_sub$Var1
pie_sub$Var1 = NULL
pie_sub = as.data.frame(t(pie_sub))
row.names(pie) =pie$Var1
pie$Var1 = NULL
pie_sub$Others = colSums(pie)- rowSums(pie_sub)
is.num <- sapply(pie_sub, is.numeric)
pie_sub[is.num] <- lapply(pie_sub[is.num], round, 2)
pie_sub = as.data.frame(t(pie_sub))
pie_sub$kindom = id[match(row.names(pie_sub), id$phylum), "kindom"]
pie_sub$phylum =row.names(pie_sub)
pie_sub = pie_sub[order(pie_sub$labels, decreasing = T),]
pie_sub$phylum = ordered(pie_sub$phylum, c("p__Proteobacteria", "p__Actinobacteria", "p__Firmicutes",     "p__Bacteroidetes",  "p__Cyanobacteria",  "p__Spirochaetes" , "p__Euryarchaeota", "p__Ascomycota", "Others" ))
who = pie_sub$phylum
col = c("p__Euryarchaeota"="#f7b2b0", "p__Actinobacteria"="#d1cb94", "p__Bacteroidetes" = "#f5b767","p__Cyanobacteria"="#abd690", "p__Firmicutes" = "#9abdf5", "p__Proteobacteria" = "#8abfb4", "p__Spirochaetes" = "#7e8691", "p__Ascomycota" = "#ba95cc", "Others" = "#613782")

p = ggplot(pie_sub, aes(x = "", y = labels, fill = pie_sub$phylum)) + geom_col(color = "black") +scale_fill_manual(values = col) + guides(fill = guide_legend(title = "phylum")) + coord_polar(theta = "y") +geom_label(aes(label = labels), position = position_stack(vjust = 0.5)) 
p$data$phylum= factor(p$data$phylum, ordered = TRUE, levels = rev(who))
p

#Silva Database
id = read.csv("~/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Clade_checker/taxonomy_7_levels.csv", sep = ",", header = F)
id = id %>% separate(V2, c("kindom", 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ";")
#kindom level
pie = data.frame(table(id$kindom))
pie$labels = pie$Freq / sum(pie$Freq)
is.num <- sapply(pie, is.numeric)
pie[is.num] <- lapply(pie[is.num], round, 2)
k = c("k__Archaea"="#ba829d" , "k__Bacteria" ="#7bb581") 
ggplot(pie, aes(x = "", y = labels, fill = Var1)) +
        geom_col(color = "black") +
        geom_label(aes(label = labels), color = c(1, 1), position = position_stack(vjust = 0.5), show.legend = FALSE) +
        guides(fill = guide_legend(title = "Answer")) +
        scale_fill_manual(values = k) +
        coord_polar(theta = "y") 
#phylum level
pie = data.frame(table(id$phylum))
pie$labels = pie$Freq / sum(pie$Freq)
pie_sub = subset(pie, !pie$labels < 0.01)
row.names(pie_sub) = pie_sub$Var1
pie_sub$Var1 = NULL
pie_sub = as.data.frame(t(pie_sub))
row.names(pie) =pie$Var1
pie$Var1 = NULL
pie_sub$Others = colSums(pie)- rowSums(pie_sub)
is.num <- sapply(pie_sub, is.numeric)
pie_sub[is.num] <- lapply(pie_sub[is.num], round, 2)
pie_sub = as.data.frame(t(pie_sub))
pie_sub$kindom = id[match(row.names(pie_sub), id$phylum), "kindom"]
pie_sub$phylum =row.names(pie_sub)
pie_sub = pie_sub[order(pie_sub$labels, decreasing = T),]
pie_sub$phylum = ordered(pie_sub$phylum, c(" p__Proteobacteria"," p__Firmicutes", " p__Bacteroidetes", " p__Actinobacteria", " p__Acidobacteria", " p__Chloroflexi", " p__Cyanobacteria",
                                           " p__Planctomycetes", " p__Patescibacteria", " p__Crenarchaeota", " p__Epsilonbacteraeota", " p__Spirochaetes", " p__Verrucomicrobia", " p__Euryarchaeota", " Others"))
who = pie_sub$phylum
col = c(" p__Proteobacteria" = "#8abfb4"," p__Firmicutes"="#9abdf5", " p__Bacteroidetes"="#f5b767", " p__Actinobacteria"="#d1cb94", " p__Acidobacteria"="#fc5203", " p__Chloroflexi"="#522917", " p__Cyanobacteria"="#173652",
        " p__Planctomycetes"="#8e4ea3", " p__Patescibacteria"="#f5ea25", " p__Crenarchaeota"="#92e092", " p__Epsilonbacteraeota"="#e0bc92", " p__Spirochaetes"="#a8681d", " p__Verrucomicrobia"="#b94bbf", " p__Euryarchaeota"="#f7b2b0", " Others"="#613782")

p = ggplot(pie_sub, aes(x = "", y = labels, fill = pie_sub$phylum)) + geom_col(color = "black") +scale_fill_manual(values = col) + guides(fill = guide_legend(title = "phylum")) + coord_polar(theta = "y") +geom_label(aes(label = labels), position = position_stack(vjust = 0.5)) 
p$data$phylum= factor(p$data$phylum, ordered = TRUE, levels = rev(who))
p

###Phylum level barplot, before vs. after
#original cpnDB_nr
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60")
cpn_7095 = fread("cpndb_7095_UtaxRenamed.csv")
names(cpn_7095) =c("NA", "Old_name", "Sequence", "New_name")
cpn_7095 = cpn_7095[-1,]
id_7095 = cpn_7095[,c("Old_name", "New_name")] %>% separate(New_name, sep = ",", c("Temp", "P", "C", "O", "F", "G", "S"))
id_7095 = id_7095 %>% separate(Temp, sep = ";", c("Taxa", "K"))
id_7095$P = gsub("p:", "", id_7095$P)
id_7095 = id_7095 %>% group_by(P) %>% summarize(n = n())
id_7095 = as.data.frame(id_7095)
id_7095 = arrange(id_7095, n)
id_7095 = as.data.frame(t(id_7095))
names(id_7095) = id_7095[1,]
id_7095 = id_7095[-1,]
id_7095 <- as.data.frame(lapply(id_7095, function(x) as.numeric(as.character(x))))
id_7095_10 = id_7095[, 72:82]
id_7095_10$Others = rowSums(id_7095) - rowSums(id_7095_10)
id_7095_10 = id_7095_10 / rowSums(id_7095_10)
id_7095_10 = as.data.frame(t(id_7095_10)) %>% mutate(Sort = "Original")
names(id_7095_10) = c("n", "Sort")
id_7095_10$P = row.names(id_7095_10)
#20556 nucleotide database
id_20556 = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/full_length_tree/id_20556_utax.csv", sep = ",", header = T)
id_20556 = id_20556[, c("Taxa", "id")] %>% separate(id, sep = ",", c("Temp", "P", "C", "O", "F", "G", "S"))
id_20556 = id_20556 %>% separate(Temp, sep = ";", c("Taxa", "K"))
id_20556$P = gsub("p:", "", id_20556$P)
id_20556 = id_20556 %>% group_by(P) %>% summarize(n = n())
id_20556 = as.data.frame(id_20556)
id_20556 = arrange(id_20556, n)
id_20556 = as.data.frame(t(id_20556))
names(id_20556) = id_20556[1,]
id_20556 = id_20556[-1,]
id_20556 = id_20556 %>% mutate_if(is.character, as.numeric)
id_20556_10 = id_20556[, 108:118]
id_20556_10$Others = rowSums(id_20556) - rowSums(id_20556_10)
id_20556_10 = id_20556_10 / rowSums(id_20556_10)
id_20556_10 = as.data.frame(t(id_20556_10)) %>% mutate(Sort = "Nucleotides")
row.names(id_20556_10) = gsub(";", "Unknown", row.names(id_20556_10))
id_20556_10$P = row.names(id_20556_10)

#uniref90 protein database
id_uniref = fread("all_cpn_uniref90_taxid.tsv")
id_uniref = replace(id_uniref, is.na(id_uniref), "Unknown")
id_uniref$P = gsub("p__", "", id_uniref$P)
id_uniref = id_uniref %>% group_by(P) %>% summarize(n = n())
id_uniref = as.data.frame(id_uniref)
id_uniref = arrange(id_uniref, n)
id_uniref = as.data.frame(t(id_uniref))
names(id_uniref) = id_uniref[1,]
id_uniref = id_uniref[-1,]
id_uniref = id_uniref %>% mutate_if(is.character, as.numeric)
id_uniref_11 = id_uniref[, 212:223]
id_uniref_11$Others = rowSums(id_uniref) - rowSums(id_uniref_11)
id_uniref_11 = id_uniref_11 / rowSums(id_uniref_11)
id_uniref_11 = as.data.frame(t(id_uniref_11)) %>% mutate(Sort = "Protein ID")
id_uniref_11$P = row.names(id_uniref_11)

id_all =rbind(id_7095_10, id_20556_10, id_uniref_11)
unique(id_all$P)
who = c("Proteobacteria",  "Firmicutes","Actinobacteria","Bacteroidetes","Cyanobacteria",
        "Ascomycota","X.Thermi.","Rhodophyta","Spirochaetes","Chordata",              
        "Crenarchaeota","Basidiomycota","Arthropoda",  "X98", "Others","Euryarchaeota", "Eukaryota_unclassified","Bacteria candidate","Unknown")
sample_order = c("Original", "Nucleotides", "Protein ID")
col = c("Proteobacteria"="#d689d4",  "Firmicutes"="#5bc9b1","Actinobacteria"="#e3857b","Bacteroidetes"="#ba9e41","Cyanobacteria"="#51b848",
        "Ascomycota"="#947289","X.Thermi."="#f49a52","Rhodophyta"="#3332e7","Spirochaetes"="#ff335d","Chordata"="#acb33d",              
        "Crenarchaeota"="#b9f3e7","Basidiomycota"="#4c74f0","Arthropoda"="#258fb4",  "X98"="#ebdf60", "Others"="#76bed6","Euryarchaeota"="#827c70", "Eukaryota_unclassified"="#827c70","Bacteria candidate"="#8f8f8f","Unknown"="#8f8f8f")
#id_all  = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Database.csv")
p = ggplot(id_all, aes(x = Sort, y = n, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Percentage") + xlab("Database")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$Sort = factor(p$data$Sort, ordered = TRUE, levels = sample_order)
p
ggsave("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Databases.pdf", width =4, height = 6)
write.csv(id_all, "Database.csv")

#####-----------updated database on 2022------------------
###Phylum level barplot, before vs. after
#original cpnDB_nr
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60")
cpndb_2023 = fread("cpn60_ref_full_seq.csv")
cpndb_2023_17279 = fread("cpndb_2023_17279.tsv")
cpndb_2023$New_name = cpndb_2023_17279[match(cpndb_2023$Sequence, cpndb_2023_17279$Sequence), "utax"]


id_2023 = cpndb_2023[,c("Name", "New_name")] %>% separate(New_name, sep = ",", c("K", "P", "C", "O", "F", "G", "S"))
id_2023$P = gsub("p:", "", id_2023$P)
id_2023$P = gsub(" Actinobacteria", "Actinobacteria", id_2023$P)
id_2023$P = gsub(" Firmicutes" , "Firmicutes", id_2023$P)
id_2023 = id_2023 %>% group_by(P) %>% summarize(n = n())
id_2023 = as.data.frame(id_2023)
id_2023 = arrange(id_2023, n)
id_2023 = as.data.frame(t(id_2023))
names(id_2023) = id_2023[1,]
id_2023 = id_2023[-1,]
id_2023 <- as.data.frame(lapply(id_2023, function(x) as.numeric(as.character(x))))
id_2023_10 = id_2023[, 86:96]
id_2023_10$Others = rowSums(id_2023) - rowSums(id_2023_10)
id_2023_10 = id_2023_10 / rowSums(id_2023_10)
id_2023_10 = as.data.frame(t(id_2023_10)) %>% mutate(Sort = "Original")
names(id_2023_10) = c("n", "Sort")
id_2023_10$P = row.names(id_2023_10)

#23545 nucleotide database
id_23545 = fread("cpn_gene_23545.tsv")
id_23545 = id_23545[, c("utax")] %>% separate(utax, sep = ",", c("Temp", "P", "C", "O", "F", "G", "S"))
id_23545 = id_23545 %>% separate(Temp, sep = ";", c("Taxa", "K"))
id_23545$P = gsub("p:", "", id_23545$P)

id_23545 = id_23545 %>% group_by(P) %>% summarize(n = n())
id_23545 = as.data.frame(id_23545)
id_23545 = arrange(id_23545, n)
id_23545 = as.data.frame(t(id_23545))
names(id_23545) = id_23545[1,]
id_23545 = id_23545[-1,]
id_23545_10 = id_23545[, 112:122]
id_23545_10 <- id_23545_10 %>%mutate_if(is.character, as.numeric)
id_23545_10$Others = 23545 - rowSums(id_23545_10)
id_23545_10 = id_23545_10 / rowSums(id_23545_10)
id_23545_10 = as.data.frame(t(id_23545_10)) %>% mutate(Sort = "Nucleotides")
row.names(id_23545_10) = gsub(";", "Unknown", row.names(id_23545_10))
id_23545_10$P = row.names(id_23545_10)

#uniref90 protein database
id_uniref = fread("all_cpn_uniref90_taxid.tsv")
id_uniref = replace(id_uniref, is.na(id_uniref), "Unknown")
id_uniref$P = gsub("p__", "", id_uniref$P)
id_uniref = id_uniref %>% group_by(P) %>% summarize(n = n())
id_uniref = as.data.frame(id_uniref)
id_uniref = arrange(id_uniref, n)
id_uniref = as.data.frame(t(id_uniref))
names(id_uniref) = id_uniref[1,]
id_uniref = id_uniref[-1,]
id_uniref = id_uniref %>% mutate_if(is.character, as.numeric)
id_uniref_11 = id_uniref[, 212:223]
id_uniref_11$Others = rowSums(id_uniref) - rowSums(id_uniref_11)
id_uniref_11 = id_uniref_11 / rowSums(id_uniref_11)
id_uniref_11 = as.data.frame(t(id_uniref_11)) %>% mutate(Sort = "Protein ID")
id_uniref_11$P = row.names(id_uniref_11)

id_all =rbind(id_2023_10, id_23545_10, id_uniref_11)
unique(id_all$P)
who = c("Proteobacteria","Firmicutes","Actinobacteria","Bacteroidetes","Cyanobacteria","Ascomycota",  
        "Planctomycetes","X94","Tenericutes" ,"Chlamydiae","Basidiomycota","Spirochaetes","Arthropoda",           
        "Chordata", "Others" ,"Euryarchaeota","Eukaryota_unclassified","Bacteria candidate","Unknown")
sample_order = c("Original", "Nucleotides", "Protein ID")
col = c("Proteobacteria"="#d689d4",  "Firmicutes"="#5bc9b1","Actinobacteria"="#e3857b","Bacteroidetes"="#ba9e41","Cyanobacteria"="#51b848","Ascomycota"="#947289",
        "Planctomycetes"="#f49a52","X94"="#3332e7","Tenericutes"="#ff335d","Chlamydiae"="#acb33d", "Basidiomycota"="#4c74f0", "Spirochaetes"="#9e1510","Arthropoda"="#258fb4",            
        "Chordata"="#b9f3e7",  "Others"="#76bed6","Euryarchaeota"="#827c70", "Eukaryota_unclassified"="#827c70","Bacteria candidate"="#8f8f8f","Unknown"="#8f8f8f")
#id_all  = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Database.csv")
p = ggplot(id_all, aes(x = Sort, y = n, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Percentage") + xlab("Database")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$Sort = factor(p$data$Sort, ordered = TRUE, levels = sample_order)
p

#remove the unknowns in protein database
id_uniref_10 = id_uniref_11[!id_uniref_11$P == "Unknown",]
id_uniref_10$new_n = id_uniref_10$n / sum(id_uniref_10$n)
id_uniref_10 = id_uniref_10[,c(2:4)]
id_uniref_10$Sort = "Protein ID (Known)"
names(id_uniref_10) = c("Sort", "P", "n")
id_all =rbind(id_2023_10, id_23545_10, id_uniref_11, id_uniref_10)
sample_order = c("Original", "Nucleotides", "Protein ID", "Protein ID (Known)")
#id_all  = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Database.csv")
p = ggplot(id_all, aes(x = Sort, y = n, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Percentage") + xlab("Database")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$Sort = factor(p$data$Sort, ordered = TRUE, levels = sample_order)
p  #export in 4x4


#####full nucleotide database vs. amplicon database vs. unmapped sequences
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/")
id_23545 = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/cpn_gene_23545.tsv")
id_23545 = id_23545[, c("utax")] %>% separate(utax, sep = ",", c("Temp", "P", "C", "O", "F", "G", "S"))
id_23545 = id_23545 %>% separate(Temp, sep = ";", c("Taxa", "K"))
id_23545$P = gsub("p:", "", id_23545$P)
id_23545 = id_23545 %>% group_by(P) %>% summarize(n = n())
id_23545 = as.data.frame(id_23545)
id_23545 = arrange(id_23545, n)
id_23545 = as.data.frame(t(id_23545))
names(id_23545) = id_23545[1,]
id_23545 = id_23545[-1,]
id_23545_4 = id_23545[, 119:122]
id_23545_4 <- id_23545_4 %>%mutate_if(is.character, as.numeric)
id_23545_4$Others = 23545 - rowSums(id_23545_4)
id_23545_4 = as.data.frame(t(id_23545_4)) %>% mutate(Sort = "Nucleotides")
row.names(id_23545_4) = gsub(";", "Unknown", row.names(id_23545_4))
id_23545_4$P = row.names(id_23545_4)

##amplicons
fasta <- readDNAStringSet(filepath = "nucleotideMSA/cpn_gene_23545.primer.fasta", format="fasta")
id_17678 =as.data.frame(names(fasta))
names(id_17678) = "id"
id_17678$taxonomy = id_17678$id_17678
id_17678 = id_17678 %>% separate(id, c("K", 'P', 'C', 'O', 'F', 'G', 'S'), sep = ",")
id_17678$P = gsub("p_", "", id_17678$P)
id_amplicon = id_17678[,c("P","K")] %>% group_by(P) %>% summarize(n = n())
id_amplicon = as.data.frame(id_amplicon)
id_amplicon[1,1] = "Others"
id_amplicon = arrange(id_amplicon, n)
id_amplicon = as.data.frame(t(id_amplicon))
names(id_amplicon) = id_amplicon[1,]
id_amplicon = id_amplicon[-1,]
id_amplicon_4 = id_amplicon[, 97:100]
id_amplicon_4 <- id_amplicon_4 %>%mutate_if(is.character, as.numeric)
id_amplicon_4$Others = 17678 - rowSums(id_amplicon_4)
id_amplicon_4 = as.data.frame(t(id_amplicon_4)) %>% mutate(Sort = "amplicons")
row.names(id_amplicon_4) = gsub(";", "Unknown", row.names(id_amplicon_4))
id_amplicon_4$P = row.names(id_amplicon_4)

##unmapped
id_23545 = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/cpn_gene_23545.tsv")
id_23545 = id_23545[, c("utax")] %>% separate(utax, sep = ",", c("Temp", "P", "C", "O", "F", "G", "S"))
id_23545 = id_23545 %>% separate(Temp, sep = ";", c("Taxa", "K"))
id_23545$P = gsub("p:", "", id_23545$P)

id_17678 = id_17678 %>% separate(K, c("Taxa", "K"), sep = ";d_")
id_17678$Taxa = gsub('"', '', id_17678$Taxa)

id_unmapped = id_23545[!id_23545$Taxa %in% id_17678$Taxa,]
id_unmapped = id_unmapped %>% group_by(P) %>% summarize(n = n())
id_unmapped = as.data.frame(id_unmapped)
id_unmapped = arrange(id_unmapped, n)
id_unmapped = as.data.frame(t(id_unmapped))
names(id_unmapped) = id_unmapped[1,]
id_unmapped = id_unmapped[-1,]
id_unmapped_4 = id_unmapped[, c("Actinobacteria", "Bacteroidetes","Firmicutes","Proteobacteria")]
id_unmapped_4 <- id_unmapped_4 %>%mutate_if(is.character, as.numeric)
id_unmapped_4$Others = 5867 - rowSums(id_unmapped_4)
id_unmapped_4 = as.data.frame(t(id_unmapped_4)) %>% mutate(Sort = "unmapped")
row.names(id_unmapped_4) = gsub(";", "Unknown", row.names(id_unmapped_4))
id_unmapped_4$P = row.names(id_unmapped_4)

id_all =rbind(id_23545_4, id_amplicon_4, id_unmapped_4)
unique(id_all$P)
who = c("Actinobacteria", "Bacteroidetes","Firmicutes","Proteobacteria", "Others")
sample_order = c("Nucleotides", "amplicons", "unmapped")
col = c("Actinobacteria"="#e3857b","Bacteroidetes"="#ba9e41", "Firmicutes"="#5bc9b1","Proteobacteria"="#d689d4", "Others"="#76bed6")
#id_all  = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/Database.csv")
p = ggplot(id_all, aes(x = Sort, y = n, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Percentage") + xlab("Database")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$Sort = factor(p$data$Sort, ordered = TRUE, levels = sample_order)
p
