##Plot the phylo tree (calculated on server)
library(ape)
library(phyloseq)
library(ggjoy)
library(dplyr)
library("ggtree")

setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/1k-3k_tree")
x = read.tree("cpn_gene_23545_cleaned_tree.newick")
id = fread("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/sintax/cpn_gene_23545_utaxID.tsv")
id_P = id[,c("id", "P")]
id_P$P = gsub("p:", "", id_P$P)
tip_labels <- x$tip.label
id_P = id_P[id_P$id %in% tip_labels,]
id_P = as.data.frame(id_P)
row.names(id_P) = id_P$id
id_P = id_P[tip_labels,]

id_P$Color = "Other"
id_P$Color[grep("Actinobacteria", id_P$P)] = "Actinobacteria"
id_P$Color[grep("Bacteroidetes", id_P$P)] = "Bacteroidetes"
id_P$Color[grep("Firmicutes", id_P$P)] = "Firmicutes"
id_P$Color[grep("Proteobacteria", id_P$P)] = "Proteobacteria"
id_P$Color[grep("Cyanobacteria", id_P$P)] = "Cyanobacteria"

ggtree(x,size = 0.4) %<+% id_P + geom_tippoint(aes(color=Color), size=2, alpha = 0.5)+ geom_treescale()
ggsave("cpn_23545_tree.pdf", width =7, height = 12)

#reorder the tree tips
who = rev(tip_labels)
y <- ape::rotateConstr(x, c(who))
ggtree(y, ladderize=F,size = 0.4) %<+% id_P + geom_tippoint(aes(color=Color), size=2, alpha = 0.5)+ geom_treescale()
ggsave("tree_test.pdf", width =7, height = 12)

#circular
ggtree(x, layout="circular", size = 0.4) %<+% id_P + geom_tippoint(aes(color=Color), size=1, alpha = 0.5)
ggsave("tree.pdf", width =10, height = 10)


#####amplicon tree
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/tree/amplicon_tree/17678Tree")
x = read.tree("cpn_23545_amplicon_cleaned_tree.newick")
tip_labels = x$tip.label

fasta<- read.fasta(file = "cpn_gene_23545_amplicon_cleaned.fasta", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
id = as.data.frame(names(fasta))
names(id) = "names"
id = id %>% separate(names, sep = ";d_", c("id", "taxonomy"))
id$id = gsub('"', '', id$id)
id$taxonomy = gsub('"', '', id$taxonomy)
id = as.data.frame(id)
row.names(id) = id$id
id = id[tip_labels,]
id = id %>% separate(taxonomy, sep = ",", c("D", "P", "C", "O", "F", "G", "S"))
id_P = id[,c("id", "P")]
id_P$P = gsub("p_", "", id_P$P)
id_P$Color = "Other"
id_P$Color[grep("Actinobacteria", id_P$P)] = "Actinobacteria"
id_P$Color[grep("Bacteroidetes", id_P$P)] = "Bacteroidetes"
id_P$Color[grep("Firmicutes", id_P$P)] = "Firmicutes"
id_P$Color[grep("Proteobacteria", id_P$P)] = "Proteobacteria"
id_P$Color[grep("Cyanobacteria", id_P$P)] = "Cyanobacteria"
ggtree(x,size = 0.4) %<+% id_P + geom_tippoint(aes(color=Color), size=2, alpha = 0.5)+ geom_treescale()
ggsave("cpn_23545_amplicon_tree.pdf", width =7, height = 12)

who = rev(tip_labels)
y <- ape::rotateConstr(x, c(who))
ggtree(y, ladderize=F,size = 0.4) %<+% id_P + geom_tippoint(aes(color=Color), size=2, alpha = 0.5)+ geom_treescale()
ggsave("cpn_23545_amplicon_tree_aligned.pdf", width =7, height = 12)



#mark those "special" taxa in the tree
m_Bacteroidetes_P = read.csv("m_Bacteroidetes_B.csv", sep = ",", header = T)
m_Bacteroidetes_nonP = read.csv("m_Bacteroidetes_nonB.csv", sep = ",", header = T)
who = unique(m_Bacteroidetes_P$variable)
who_1 = unique(m_Bacteroidetes_nonP$variable)
who = as.data.frame(who)
who_1 = as.data.frame(who_1)
who$Bacteroidetes = "lowID_within"
who_1$Bacteroidetes = "highID_between"
names(who_1) = c("who", "Bacteroidetes")
who_1$Class = m_name[match(who_1$who, m_name$id), "C"]
who_1$Phylum = m_name[match(who_1$who, m_name$id), "P"]
who_1$Class = paste(who_1$Phylum, who_1$Class, sep = ";")
who_1$Class = gsub("p__", "", who_1$Class)
who_1$Class = gsub("c__", "", who_1$Class)
who_1$Phylum = NULL
who$Class = "z"
who = rbind(who, who_1)

id_P$Bacteroidetes = who[match(id_P$id, who$who), "Bacteroidetes"]
id_P$Class = who[match(id_P$id, who$who), "Class"]

ggtree(x, size = 0.4) %<+% id_P + geom_tippoint(aes(color=P), size=2) + geom_tippoint(aes(shape = shape), size = 2, alpha = 0.5) 
ggtree(x, layout="circular", size = 0.4) %<+% id_P + geom_tippoint(aes(color=P), size=2, alpha = 0.8) + geom_tippoint(aes(shape = shape), size = 2, alpha = 0.4) + geom_text(aes(label = outliers), na.rm = TRUE, size = 3, hjust=0,vjust=1)
ggsave("test.pdf", width =7, height = 12)

####SILVA tree
x = read.tree("silva_V3V4_tree.newick")
id = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/SILVA/subset_20000/id_P_SILVA_20000.csv", sep = ",", header = T)
id_P = id[,c("V1", "P")]
id_P$P = gsub("p__", "", id_P$P)
id_P$P = gsub("Other", " Other", id_P$P)
ggtree(x) %<+% id_P + geom_tippoint(aes(color=P), size=1)
#circular
ggtree(x, layout="circular", size = 0.4) %<+% id_P + geom_tippoint(aes(color=P), size=1, alpha = 0.5)
ggtree(x, size = 0.4) %<+% id_P + geom_tippoint(aes(color=P), size=2, alpha = 0.5)+ geom_treescale()
ggsave("SILVA_V4_tree.pdf", width =7, height = 10)
