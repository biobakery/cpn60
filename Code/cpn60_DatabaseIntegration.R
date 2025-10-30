library("data.table")
library("httr")
library("jsonlite")
library("tidyr")
library("seqinr")
setwd("Documents/Lea/Harvard/MBTA_RNA/Cpn60/")

##check on the old and the new cpndb
cpndb_ori = fread("cpndb_nr_nut_seq_7095.csv")
cpndb_2023 = fread("cpn60_ref_full_seq.csv")
cpndb_7095 = fread("cpndb_7095_UtaxRenamed.csv")
cpndb_2023_ann = read.fasta(file = "cpn60_ref_full_seq.txt", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)

cpndb_ori$utax = cpndb_7095[match(cpndb_ori$Sequence, cpndb_7095$Sequence),"Name"]
cpndb_2023$utax = cpndb_ori[match(cpndb_2023$Name, cpndb_ori$Name), "utax"]
cpndb_temp = cpndb_ori[!cpndb_ori$Name %in% cpndb_2023$Name,] ##these amplicons are not included in the new cpndb. Need to mannually include
cpndb_2023 = rbind(cpndb_2023, cpndb_temp)

cpndb_2023$utax = gsub(".*\\;tax=", "",cpndb_2023$utax)
cpndb_2023_ann_id = as.data.frame(names(cpndb_2023_ann))
cpndb_2023_ann_id = cpndb_2023_ann_id %>% separate(`names(cpndb_2023_ann)`, sep = "-", c("id", "Accession", "Genus", "Species", "Strain", "temp"))
cpndb_2023_ann_id$S = paste(cpndb_2023_ann_id$Genus, cpndb_2023_ann_id$Species, sep = "_")


##merge the new updated cpndb with the 16280 sequences
#rename the 16280 sequences
cpn_gene_16280 = read.fasta(file = "cpn_gene_16280.fasta", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)

cpn_gene_16280_id = as.data.frame(names(cpn_gene_16280))
cpn_gene_16280_id = cpn_gene_16280_id %>% separate(`names(cpn_gene_16280)`, sep = "\\|", c("Gene", "Taxa", "Uniref90", "Uniref50"))
cpn_gene_16280_id$Taxa = gsub("\\.", ",", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$Taxa = gsub("k__", "d:", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$Taxa = gsub("p__", "p:", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$Taxa = gsub("c__", "c:", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$Taxa = gsub("f__", "f:", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$Taxa = gsub("g__", "g:", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$Taxa = gsub("s__", "s:", cpn_gene_16280_id$Taxa)
cpn_gene_16280_id$id = paste("Taxa", row.names(cpn_gene_16280_id), sep = "_")
cpn_gene_16280_id$Taxa = paste(cpn_gene_16280_id$id, cpn_gene_16280_id$Taxa, sep = ";")
cpn_gene_16280_id$Taxa = gsub("d:", "tax=d:", cpn_gene_16280_id$Taxa)
names(cpn_gene_16280) = cpn_gene_16280_id$Taxa
write.fasta(as.list(cpn_gene_16280), names = names(cpn_gene_16280), file.out = "cpn_gene_16280.fasta")
write.table(cpn_gene_16280_id, file = "cpn_gene_16280_id.tsv", sep = "\t", quote = FALSE, row.names = FALSE )

#rename the cpndb_2023 database
cpn_gene_16280_id$S = gsub(".*\\,s:", "",cpn_gene_16280_id$Taxa)
cpndb_2023$S = cpndb_2023_ann_id[match(cpndb_2023$Name, cpndb_2023_ann_id$id), "S"]
cpndb_2023$temp = cpn_gene_16280_id[match(cpndb_2023$S, cpn_gene_16280_id$S), "Taxa"]
cpndb_2023$temp = gsub(".*\\;tax=", "",cpndb_2023$temp)
write.table(cpndb_2023, file = "cpndb_2023_17279.tsv", sep = "\t", quote = FALSE, row.names = FALSE )

cpndb_2023_17279 = fread("cpndb_2023_17279.tsv")

cpn_gene_16280 = fread("cpn_gene_16280.csv")
names(cpn_gene_16280) = c("utax", "Sequence")
cpn_gene_16280$utax = gsub(".*\\;tax=", "",cpn_gene_16280$utax)

check = cpndb_2023_17279[cpndb_2023_17279$Sequence %in% cpn_gene_16280$Sequence,]
check$ann_16280 = cpn_gene_16280[match(check$Sequence, cpn_gene_16280$Sequence), "utax"] ##check on the annotations and see if they match

cpndb_2023_17279 = cpndb_2023_17279[!cpndb_2023_17279$Sequence %in% check$Sequence,]
cpn_gene_28496 = rbind(cpn_gene_16280, cpndb_2023_17279[, c("utax", "Sequence")])
cpn_gene_28496 = unique(cpn_gene_28496)
#write.table(cpn_gene_28496, file = "cpn_gene_28496.tsv", sep = "\t", quote = FALSE, row.names = FALSE )

cpn_gene_28496 = fread("cpn_gene_28496.tsv")
cpn_gene_28496 = unique(cpn_gene_28496)
cpn_gene_23545 = cpn_gene_28496
cpn_gene_23545$Taxa = paste("Taxa", row.names(cpn_gene_23545), sep = "_")
cpn_gene_23545$utax = paste(cpn_gene_23545$Taxa, cpn_gene_23545$utax, sep = ";")

write.table(cpn_gene_23545[, c("utax", "Sequence")], file = "cpn_gene_23545.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
write.table(cpn_gene_23545[, c("Taxa", "Sequence")], file = "cpn_gene_23545_rename.tsv", sep = "\t", quote = FALSE, row.names = FALSE )

#extract 1k to 3k sequences for MSA and phylo tree
library(Biostrings)
input_file <- "cpn_gene_23545_rename.fasta"
output_file <- "cpn_gene_23545_1k-3k.fasta"
min_length <- 1000
max_length <- 3000

sequences <- readDNAStringSet(input_file)
filtered_sequences <- sequences[width(sequences) >= min_length & width(sequences) <= max_length]
writeXStringSet(filtered_sequences, output_file)

##extract the amplicons using primer1 and primer2 and then merge+deduplicate
primer1 <- read.fasta(file = "nucleotideMSA/cpn_gene_23545.primer1.fasta", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
primer2 <- read.fasta(file = "nucleotideMSA/cpn_gene_23545.primer2.fasta", seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
id_1 = names(primer1)
id_2 = names(primer2)
id_add = id_2[!id_2 %in% id_1]
add = primer2[names(primer2) %in% id_add]
primer = c(primer1, add)
write.fasta(sequences = primer, names = names(primer), file.out = "nucleotideMSA/cpn_gene_23545.primer.fasta")

##for RDP classifier training
amplicon = fread("RDPClassifier/cpn_gene_23545_primer.tsv")
amplicon_id = amplicon[,"Name"] %>% separate(Name, sep = ";d", c("id", "taxonomy"))
amplicon_id = amplicon_id %>% separate(taxonomy, sep = ",", c("superkingdom","phylum","class", "order", "family", "geus", "species"))
write.table(amplicon_id, file = "RDPClassifier/cpn_gene_23545_primer_id.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
amplicon$Name = gsub(";.*", "",amplicon$Name)
write.table(amplicon, file = "RDPClassifier/cpn_gene_23545_primer_rename.tsv", sep = "\t", quote = FALSE, row.names = FALSE )


