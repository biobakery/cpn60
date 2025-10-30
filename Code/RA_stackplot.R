library(data.table)



##Extract 23230 uniref90s from HMP2 MGX dataset
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/HMP2")
tax_file = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/cpn_gene_20556_utaxID.csv", sep = ",", header = T)
tax_file$G = gsub("g:_", "g__", tax_file$G)
uniref90 = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/HMP2/uniref90_total_23230.csv", sep = ",", header = T)
who = c("# Gene Family", "UNMAPPED",uniref90$who_total_23230 )

meta = read.csv("hmp2_metadata.csv", sep = ",", header = T)
mtx = meta[meta$data_type == "metatranscriptomics",]
mgx = meta[meta$data_type == "metagenomics",]
mtx_sample = unique(mtx$External.ID)
mgx_sample = unique(mgx$External.ID)
#sample = mtx_sample[mtx_sample %in% mgx_sample]
sample = paste(mgx_sample, "Abundance-RPKs", sep = "_")

sample_200 = sample[1:200]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_200))
mgx_200 = arrange(test, `# Gene Family`)
mgx_200$uniref90 = mgx_200$`# Gene Family`
mgx_200$uniref90 = gsub("\\|.*", "",mgx_200$uniref90)
mgx_200 = mgx_200[mgx_200$uniref90 %in% who,]
dim(mgx_200)
write.table(mgx_200, file = "mgx_200.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_400 = sample[201:400]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_400)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_400 = arrange(test, `# Gene Family`)
mgx_400$uniref90 = mgx_400$`# Gene Family`
mgx_400$uniref90 = gsub("\\|.*", "",mgx_400$uniref90)
mgx_400 = mgx_400[mgx_400$uniref90 %in% who,]
dim(mgx_400)
write.table(mgx_400, file = "mgx_400.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_600 = sample[401:600]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_600)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_600 = arrange(test, `# Gene Family`)
mgx_600$uniref90 = mgx_600$`# Gene Family`
mgx_600$uniref90 = gsub("\\|.*", "",mgx_600$uniref90)
mgx_600 = mgx_600[mgx_600$uniref90 %in% who,]
dim(mgx_600)
write.table(mgx_600, file = "mgx_600.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_800 = sample[601:800]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_800)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_800 = arrange(test, `# Gene Family`)
mgx_800$uniref90 = mgx_800$`# Gene Family`
mgx_800$uniref90 = gsub("\\|.*", "",mgx_800$uniref90)
mgx_800 = mgx_800[mgx_800$uniref90 %in% who,]
dim(mgx_800)
write.table(mgx_800, file = "mgx_800.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_1000 = sample[801:1000]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_1000)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_1000 = arrange(test, `# Gene Family`)
mgx_1000$uniref90 = mgx_1000$`# Gene Family`
mgx_1000$uniref90 = gsub("\\|.*", "",mgx_1000$uniref90)
mgx_1000 = mgx_1000[mgx_1000$uniref90 %in% who,]
dim(mgx_1000)
write.table(mgx_1000, file = "mgx_1000.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_1200 = sample[1001:1200]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_1200)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_1200 = arrange(test, `# Gene Family`)
mgx_1200$uniref90 = mgx_1200$`# Gene Family`
mgx_1200$uniref90 = gsub("\\|.*", "",mgx_1200$uniref90)
mgx_1200 = mgx_1200[mgx_1200$uniref90 %in% who,]
dim(mgx_1200)
write.table(mgx_1200, file = "mgx_1200.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_1400 = sample[1201:1400]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_1400)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_1400 = arrange(test, `# Gene Family`)
mgx_1400$uniref90 = mgx_1400$`# Gene Family`
mgx_1400$uniref90 = gsub("\\|.*", "",mgx_1400$uniref90)
mgx_1400 = mgx_1400[mgx_1400$uniref90 %in% who,]
dim(mgx_1400)
write.table(mgx_1400, file = "mgx_1400.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample_1638 = sample[1401:1638]
test = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", sample_1638)) #the complete file is huge, ~30G; here we only have the mgx file so I'll first extract these 200
mgx_1638 = arrange(test, `# Gene Family`)
mgx_1638$uniref90 = mgx_1638$`# Gene Family`
mgx_1638$uniref90 = gsub("\\|.*", "",mgx_1638$uniref90)
mgx_1638 = mgx_1638[mgx_1638$uniref90 %in% who,]
dim(mgx_1638)
write.table(mgx_1638, file = "mgx_1638.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
mgx_all <- cbind(mgx_200, mgx_400, mgx_600, mgx_800, mgx_1000, mgx_1200, mgx_1400, mgx_1638)
write.table(mgx_all, file = "mgx_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

mgx_all = fread("mgx_all.tsv")
dim(mgx_all)
mgx_all$temp = rowSums(mgx_all[,c(2:1638)])
mgx_all = mgx_all[!mgx_all$temp ==0, ]
mgx_all = mgx_all[grepl("\\|", mgx_all$`# Gene Family`),] #selecting only the rows where the # Gene Family column contains the "|" character
mgx_all$taxa = mgx_all$`# Gene Family`
mgx_all$uniref90 = NULL
mgx_all = mgx_all %>% separate(taxa, sep = "\\|", c("uniref90", "taxa"))


##plot the stackplot
mgx = fread("MGX/mgx_all.tsv")
mgx = mgx[, -c("# Gene Family", "temp", "uniref90")]
mgx = mgx %>% group_by(taxa) %>% summarise_all(sum)
mgx = as.data.frame(mgx)
sample_order= names(sort(colMeans(mgx[2:1638]), decreasing = TRUE))
mgx = melt(mgx) #sort the order of samples before melt
mgx$G = mgx$taxa
mgx = mgx %>% separate(taxa, sep = "\\.", c("G", "S"))
mgx$P = tax_file[match(mgx$G, tax_file$G), "P"]
mgx$P[grep("unclassified", mgx$G)] = "unclassified"
unique(mgx$P) #check if there's any "NA", if yes, write csv and mannually correct the taxa

mgx = fread("MGX/mgx_all_melt.tsv")
who = unique(mgx$P)
col = c("unclassified"="#f7ae4f","p:Spirochaetes" = "#ed078d","p:Synergistetes"="#b107a4",    "p:Lentisphaerae"="#2a3150","p:Euryarchaeota"="#f7b2b0", "p:Fusobacteria"=  "#afbce3","p:Verrucomicrobia" = "#6f6285", "p:Actinobacteria"="#e3857b",  "p:Firmicutes" = "#5bc9b1","p:Bacteroidetes" = "#ba9e41", "p:Proteobacteria" ="#d689d4")
mgx$variable = gsub("_Abundance-RPKs", "", mgx$variable)
meta = meta[meta$External.ID %in% mgx$variable,]
#relative abundance
mgx$sum <- ave(mgx$value, mgx$variable, FUN=sum)
mgx$RA = mgx$value / mgx$sum

#sort the sample by Bacteroidetes
sample_order = mgx[,c("variable", "P", "RA")] %>% group_by(variable, P) %>% summarise_all(sum)
sample_order = as.data.frame(sample_order[sample_order$P == "p:Bacteroidetes",])
sample_order = arrange(sample_order, RA)
sample_order = rev(sample_order$variable)
mgx$variable = ordered(mgx$variable, c(sample_order))
p = ggplot(mgx, aes(x = variable, y = RA, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("Samples")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank() ) #+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_all_ra.pdf", width =6, height = 3)

#mgx_genus level
mgx_G = mgx[,c("G", "variable", "RA")] %>% group_by(variable, G) %>% summarise_all(sum)
mgx_G = as.data.frame(mgx_G)
mgx_G = na.omit(mgx_G)
who = mgx_G[,c("G", "RA")] %>% group_by(G) %>% summarise_all(sum)
who = as.data.frame(who)
who = arrange(who, RA)
who = rev(who$G)[1:15]
mgx_G = mgx_G[mgx_G$G %in% who, ]
mgx_G = mgx_G %>% pivot_wider(names_from = variable, values_from=RA, values_fill = list(ID = 0))
mgx_G = as.data.frame(t(mgx_G))
names(mgx_G) = mgx_G[1,]
mgx_G = mgx_G[-1,]
mgx_G[] <- lapply(mgx_G, function(x) as.numeric(as.character(x)))
mgx_G$Others = 1 -rowSums(mgx_G)
mgx_G$variable = row.names(mgx_G)
mgx_G = melt(mgx_G)
names(mgx_G) = c("variable", "G", "RA")
who = c("g__Bacteroides", "g__Faecalibacterium", "g__Prevotella", "g__Parabacteroides", "g__Roseburia",
        "g__Alistipes","g__Lachnospiraceae_unclassified","g__Eubacterium","g__Firmicutes_unclassified" ,"g__Escherichia","g__Blautia",
        "g__Dialister", "g__Clostridium" ,"g__Ruminococcaceae_unclassified","Others",  "unclassified" )
col = c("g__Bacteroides"="#84cf26", "g__Faecalibacterium"="#d07b09", "g__Prevotella"="#9ff8e5","g__Parabacteroides"="#657f26","g__Roseburia"="#5e4210",
        "g__Alistipes"="#9c4afd","g__Lachnospiraceae_unclassified"="#58c6f9","g__Eubacterium"="#e887bb","g__Firmicutes_unclassified"="#02ce88", "g__Escherichia"="#dee020",
        "g__Blautia"="#d85b75","g__Dialister"="#0d8561","g__Clostridium"="#c3802e", "g__Ruminococcaceae_unclassified"="#2918cd", "Others"="#fceac5","unclassified"="#f7ae4f" )
mgx_G$variable = ordered(mgx_G$variable, c(sample_order))
mgx_G$G = ordered(mgx_G$G, c(who))

p = ggplot(mgx_G, aes(x = variable, y = RA, fill = G)) + geom_bar(stat="identity")  + scale_fill_manual(values = col)+
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("Samples")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() )
p$data$G = factor(p$data$G, ordered = TRUE, levels = who)
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_all_ra_cpn_G.pdf", width =6, height = 3)

###comparie this with metaphlan profile
metaphlan = fread("MGX/HMP2_metaphlan3.1.0_metaphlan_taxonomic_profiles.tsv")
names(metaphlan) = gsub("_taxonomic_profile", "", names(metaphlan))
sample = unique(mgx$variable)
sample = c("# taxonomy", sample_order)
metaphlan = as.data.frame(metaphlan)
metaphlan = metaphlan[ , c(sample)]
#check on the taxonomy
metaphlan[1:40, 1:3]
metaphlan = metaphlan[c(1,6:17,20,24), ]
metaphlan =as.data.frame(metaphlan)
row.names(metaphlan) = metaphlan$`# taxonomy`
metaphlan$`# taxonomy` =  NULL
metaphlan = as.data.frame(t(metaphlan))
metaphlan = metaphlan / 100
metaphlan$Others = 1 - rowSums(metaphlan)
metaphlan = as.data.frame(t(metaphlan)) 
metaphlan[metaphlan < 0] <- 0
metaphlan$Taxa = row.names(metaphlan)

metaphlan = melt(metaphlan, id.vars = "Taxa")
metaphlan$Taxa = gsub(".*\\|", "",metaphlan$Taxa)
names(metaphlan) = c("P", "variable", "value")
metaphlan$P = gsub("UNKNOWN", "unclassified", metaphlan$P)
metaphlan$P = gsub("p__", "p:", metaphlan$P)
who = unique(metaphlan$P)
who = c("p:Bacteroidetes","p:Firmicutes","p:Actinobacteria","p:Proteobacteria","p:Ascomycota","p:Tenericutes","p:Spirochaetes","p:Euryarchaeota","p:Fusobacteria","p:Synergistetes",             
"p:Lentisphaerae","p:Verrucomicrobia","p:Eukaryota_unclassified","p:Candidatus_Melainabacteria", "Others","unclassified")
col = c("p:Bacteroidetes" = "#ba9e41","p:Firmicutes" = "#5bc9b1","p:Actinobacteria"="#e3857b", "p:Proteobacteria" ="#d689d4","p:Ascomycota"="#578f82","p:Tenericutes"="#905bc9", "p:Spirochaetes" = "#ed078d","p:Euryarchaeota"="#f7b2b0","p:Fusobacteria"="#afbce3","p:Synergistetes"="#b107a4",
        "p:Lentisphaerae"="#2a3150", "p:Verrucomicrobia" = "#6f6285",  "p:Eukaryota_unclassified"="#f7b2b0","p:Candidatus_Melainabacteria"="#e8e161","Others"="#fceac5","unclassified"="#f7ae4f")
metaphlan$variable = ordered(metaphlan$variable, c(sample_order))
metaphlan$P = ordered(metaphlan$P, c(who))

p = ggplot(metaphlan, aes(x = variable, y = value, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("metaphlan3.1")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank() )
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_all_ra_metaphlan.pdf", width =6, height = 3)

p = ggplot(mgx, aes(x = variable, y = RA, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("Samples")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank() ) #+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_all_ra.pdf", width =6, height = 3)

#genus level
metaphlan = fread("MGX/HMP2_metaphlan3.1.0_metaphlan_taxonomic_profiles.tsv")
names(metaphlan) = gsub("_taxonomic_profile", "", names(metaphlan))
sample = unique(mgx$variable)
sample = c("# taxonomy", sample)
metaphlan = as.data.frame(metaphlan)
metaphlan = metaphlan[ , c(sample)]
metaphlan_G = metaphlan[grepl("\\|g__", metaphlan$`# taxonomy`),]
unclassified = metaphlan[metaphlan$`# taxonomy` == "UNKNOWN",]
metaphlan_G = rbind(metaphlan_G, unclassified)
metaphlan_G = metaphlan_G[!grepl("\\|s__", metaphlan_G$`# taxonomy`),] #selecting only the rows where the # taxonomy column does not contain the string "|s__".
metaphlan_G = metaphlan_G %>% separate('# taxonomy', sep = "\\|", c("K", "P", "C", "O", "F", "G")) %>% select(-K, -P, -C, -O, -F)
metaphlan_G = as.data.frame(t(metaphlan_G))
names(metaphlan_G) = metaphlan_G[1,]
metaphlan_G = metaphlan_G[-1,]
metaphlan_G[] <- lapply(metaphlan_G, function(x) as.numeric(as.character(x)))
names(metaphlan_G)[274] = "unclassified"
who = names(sort(colMeans(metaphlan_G), decreasing = TRUE))[1:15]
metaphlan_G = metaphlan_G[,names(metaphlan_G) %in% who]
metaphlan_G$Others = 100 - rowSums(metaphlan_G)

metaphlan_G$variable = row.names(metaphlan_G)
metaphlan_G = melt(metaphlan_G)
names(metaphlan_G) = c("variable", "G", "RA")
who = c("g__Bacteroides","g__Faecalibacterium","g__Prevotella", "g__Alistipes","g__Roseburia","g__Parabacteroides","g__Lachnospiraceae_unclassified",
        "g__Eubacterium","g__Escherichia","g__Clostridium","g__Firmicutes_unclassified","g__Akkermansia","g__Blautia" ,"g__Flavonifractor", "Others",  "unclassified" )
col = c("g__Bacteroides"="#84cf26", "g__Faecalibacterium"="#d07b09", "g__Prevotella"="#9ff8e5","g__Alistipes"="#9c4afd","g__Roseburia"="#5e4210","g__Parabacteroides"="#657f26",
        "g__Lachnospiraceae_unclassified"="#58c6f9","g__Eubacterium"="#e887bb","g__Escherichia"="#dee020","g__Clostridium"="#c3802e", "g__Firmicutes_unclassified"="#02ce88", 
        "g__Akkermansia"="#3e344d","g__Blautia"="#d85b75","g__Flavonifractor"="#4f2e05","g__Dialister"="#0d8561","g__Ruminococcaceae_unclassified"="#2918cd", "Others"="#fceac5","unclassified"="#f7ae4f" )
metaphlan_G$variable = ordered(metaphlan_G$variable, c(sample_order))
metaphlan_G$G = ordered(metaphlan_G$G, c(who))

p = ggplot(metaphlan_G, aes(x = variable, y = RA, fill = G)) + geom_bar(stat="identity")  + scale_fill_manual(values = col)+
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("metaphlan4.2")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank() )
p$data$G = factor(p$data$G, ordered = TRUE, levels = who)
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_all_ra_metaphlan_G.pdf", width =6.2, height = 3)

#species level
metaphlan = fread("MGX/HMP2_metaphlan3.1.0_metaphlan_taxonomic_profiles.tsv")
names(metaphlan) = gsub("_taxonomic_profile", "", names(metaphlan))
metaphlan = as.data.frame(metaphlan)
metaphlan_S = metaphlan[grepl("\\|s__", metaphlan$`# taxonomy`),]
unclassified = metaphlan[metaphlan$`# taxonomy` == "UNKNOWN",]
metaphlan_S = rbind(metaphlan_S, unclassified)
metaphlan_S = metaphlan_S %>% separate('# taxonomy', sep = "\\|", c("K", "P", "C", "O", "F", "G", "S")) %>% select(-K, -P, -C, -O, -F, -G)
metaphlan_S = as.data.frame(t(metaphlan_S))
names(metaphlan_S) = metaphlan_S[1,]
metaphlan_S = metaphlan_S[-1,]
metaphlan_S[] <- lapply(metaphlan_S, function(x) as.numeric(as.character(x)))
names(metaphlan_S)[955] = "unclassified"
who = names(sort(colMeans(metaphlan_S), decreasing = TRUE))

metaphlan_S$variable = row.names(metaphlan_S)
metaphlan_S = melt(metaphlan_S)
names(metaphlan_S) = c("variable", "S", "RA")
#save the metaphlan_S for comparison

###do the same for MTX samples
setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/HMP2")
tax_file = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/cpn_gene_20556_utaxID.csv", sep = ",", header = T)
tax_file$G = gsub("g:_", "g__", tax_file$G)
uniref90 = read.csv("uniref90_total_23230.csv", sep = ",", header = T)
who = c(uniref90$who_total_23230)
meta = read.csv("hmp2_metadata.csv", sep = ",", header = T)
mtx = meta[meta$data_type == "metatranscriptomics",]
mtx_sample = unique(mtx$External.ID)

#sample = mtx_sample[mtx_sample %in% mgx_sample]
sample = paste(mtx_sample, "Abundance-RPKs", sep = "_")
sample_200 = sample[1:200]
mtx_200 = fread("MTX/genefamilies.tsv", select = c("# Gene Family", sample_200))
mtx_200 = arrange(mtx_200, `# Gene Family`)
mtx_200$uniref90 = mtx_200$`# Gene Family`
mtx_200$uniref90 = gsub("\\|.*", "",mtx_200$uniref90)
mtx_200 = mtx_200[mtx_200$uniref90 %in% who,]
dim(mtx_200)
write.table(mtx_200, file = "mtx_200.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample = paste(mtx_sample, "Abundance-RPKs", sep = "_")
sample_400 = sample[201:400]
mtx_400 = fread("MTX/genefamilies.tsv", select = c("# Gene Family", sample_400))
mtx_400 = arrange(mtx_400, `# Gene Family`)
mtx_400$uniref90 = mtx_400$`# Gene Family`
mtx_400$uniref90 = gsub("\\|.*", "",mtx_400$uniref90)
mtx_400 = mtx_400[mtx_400$uniref90 %in% who,]
dim(mtx_400)
write.table(mtx_400, file = "mtx_400.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample = paste(mtx_sample, "Abundance-RPKs", sep = "_")
sample_600 = sample[401:600]
mtx_600 = fread("MTX/genefamilies.tsv", select = c("# Gene Family", sample_600))
mtx_600 = arrange(mtx_600, `# Gene Family`)
mtx_600$uniref90 = mtx_600$`# Gene Family`
mtx_600$uniref90 = gsub("\\|.*", "",mtx_600$uniref90)
mtx_600 = mtx_600[mtx_600$uniref90 %in% who,]
dim(mtx_600)
write.table(mtx_600, file = "mtx_600.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sample = paste(mtx_sample, "Abundance-RPKs", sep = "_")
sample_835 = sample[601:835]
mtx_835 = fread("MTX/genefamilies.tsv", select = c("# Gene Family", sample_835))
mtx_835 = arrange(mtx_835, `# Gene Family`)
mtx_835$uniref90 = mtx_835$`# Gene Family`
mtx_835$uniref90 = gsub("\\|.*", "",mtx_835$uniref90)
mtx_835 = mtx_835[mtx_835$uniref90 %in% who,]
dim(mtx_835)
write.table(mtx_835, file = "mtx_835.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

mtx_all <- cbind(mtx_200, mtx_400, mtx_600, mtx_835)
write.table(mtx_all, file = "mtx_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

mtx_all = fread("MTX/mtx_all.tsv")
dim(mtx_all)
mtx_all$temp = rowSums(mtx_all[,c(2:814)])
mtx_all = mtx_all[!mtx_all$temp ==0, ]
mtx_all = mtx_all[grepl("\\|", mtx_all$`# Gene Family`),] #selecting only the rows where the # Gene Family column contains the "|" character
mtx_all$taxa = mtx_all$`# Gene Family`
mtx_all$uniref90 = NULL
mtx_all = mtx_all %>% separate(taxa, sep = "\\|", c("uniref90", "taxa"))

mtx = arrange(mtx_all, `# Gene Family`)
mtx$uniref90 = mtx$`# Gene Family`
mtx$uniref90 = gsub("\\|.*", "",mtx$uniref90)
mtx = mtx[mtx$uniref90 %in% who,]
dim(mtx)
write.table(mtx, file = "mtx_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

##plot the stackplot
mtx = fread("MTX/mtx_all.tsv")
mtx = mtx[, -c("# Gene Family", "temp", "uniref90")]
mtx = mtx %>% group_by(taxa) %>% summarise_all(sum)
mtx = as.data.frame(mtx)
mtx = melt(mtx) 
mtx$G = mtx$taxa
mtx = mtx %>% separate(taxa, sep = "\\.", c("G", "S"))
mtx$P = tax_file[match(mtx$G, tax_file$G), "P"]
mtx$P[grep("unclassified", mtx$G)] = "unclassified"
unique(mtx$P) #check if there's any "NA", if yes, write csv and mannually correct the taxa

mtx = fread("MTX/mtx_all_melt.tsv")
mtx$variable = gsub("_Abundance-RPKs", "", mtx$variable)

#relative abundance
mtx$sum <- ave(mtx$value, mtx$variable, FUN=sum)
mtx$RA = mtx$value / mtx$sum
mtx$variable = gsub("_Abundance-RPKs", "", mtx$variable)
who = c("p:Bacteroidetes","p:Firmicutes","p:Actinobacteria","p:Proteobacteria","p:Ascomycota","p:Tenericutes","p:Spirochaetes","p:Euryarchaeota","p:Fusobacteria","p:Synergistetes",             
        "p:Lentisphaerae","p:Verrucomicrobia","p:Eukaryota_unclassified","p:Candidatus_Melainabacteria", "Others","unclassified")
col = c("p:Bacteroidetes" = "#ba9e41","p:Firmicutes" = "#5bc9b1","p:Actinobacteria"="#e3857b", "p:Proteobacteria" ="#d689d4","p:Ascomycota"="#578f82","p:Tenericutes"="#905bc9", "p:Spirochaetes" = "#ed078d","p:Euryarchaeota"="#f7b2b0","p:Fusobacteria"="#afbce3","p:Synergistetes"="#b107a4",
        "p:Lentisphaerae"="#2a3150", "p:Verrucomicrobia" = "#6f6285",  "p:Eukaryota_unclassified"="#f7b2b0","p:Candidatus_Melainabacteria"="#e8e161","Others"="#fceac5","unclassified"="#f7ae4f")

#sort the sample by Bacteroidetes in MGX samples
mgx = fread("MGX/mgx_all_melt.tsv")
who = unique(mgx$P)
mgx$variable = gsub("_Abundance-RPKs", "", mgx$variable)
meta = meta[meta$External.ID %in% mgx$variable,]
mgx$sum <- ave(mgx$value, mgx$variable, FUN=sum)
mgx$RA = mgx$value / mgx$sum
sample_order = mgx[,c("variable", "P", "RA")] %>% group_by(variable, P) %>% summarise_all(sum)
sample_order = as.data.frame(sample_order[sample_order$P == "p:Bacteroidetes",])
sample_order = arrange(sample_order, RA)
sample_order = rev(sample_order$variable)

who = c("p:Bacteroidetes","p:Firmicutes","p:Actinobacteria","p:Proteobacteria","p:Ascomycota","p:Tenericutes","p:Spirochaetes","p:Euryarchaeota","p:Fusobacteria","p:Synergistetes",             
        "p:Lentisphaerae","p:Verrucomicrobia","p:Eukaryota_unclassified","p:Candidatus_Melainabacteria", "Others","unclassified")
col = c("p:Bacteroidetes" = "#ba9e41","p:Firmicutes" = "#5bc9b1","p:Actinobacteria"="#e3857b", "p:Proteobacteria" ="#d689d4","p:Ascomycota"="#578f82","p:Tenericutes"="#905bc9", "p:Spirochaetes" = "#ed078d","p:Euryarchaeota"="#f7b2b0","p:Fusobacteria"="#afbce3","p:Synergistetes"="#b107a4",
        "p:Lentisphaerae"="#2a3150", "p:Verrucomicrobia" = "#6f6285",  "p:Eukaryota_unclassified"="#f7b2b0","p:Candidatus_Melainabacteria"="#e8e161","Others"="#fceac5","unclassified"="#f7ae4f")
mtx$P = ordered(mtx$P, c(who))
mtx$variable = ordered(mtx$variable, c(sample_order))

p = ggplot(mtx, aes(x = variable, y = RA, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Read counts") + xlab("cpn60 uniref90")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank() ) #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mtx_all_ra.pdf", width =5, height = 3)

mgx = mgx[mgx$variable %in% mtx$variable,]
mgx$P = ordered(mgx$P, c(who))
mgx$variable = ordered(mgx$variable, c(sample_order))
p = ggplot(mgx, aes(x = variable, y = RA, fill = P)) + geom_bar(stat="identity") + scale_fill_manual(values = col) +
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("MGX cpn60 uniref90")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() ) #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p$data$P = factor(p$data$P, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_813_ra.pdf", width =5, height = 3)

#mtx_genus level
mtx_G = mtx[,c("G", "variable", "RA")] %>% group_by(variable, G) %>% summarise_all(sum)
mtx_G = as.data.frame(mtx_G)
mtx_G = mtx_G[!mtx_G$RA == "NaN",]
who = mtx_G[,c("G", "RA")] %>% group_by(G) %>% summarise_all(sum)
who = as.data.frame(who)
who = arrange(who, RA)
who = rev(who$G)[1:15]
mtx_G = mtx_G[mtx_G$G %in% who, ]
mtx_G = mtx_G %>% pivot_wider(names_from = variable, values_from=RA, values_fill = list(ID = 0))
mtx_G = as.data.frame(t(mtx_G))
names(mtx_G) = mtx_G[1,]
mtx_G = mtx_G[-1,]
mtx_G[] <- lapply(mtx_G, function(x) as.numeric(as.character(x)))
mtx_G$Others = 1 -rowSums(mtx_G)
mtx_G$variable = row.names(mtx_G)
mtx_G = melt(mtx_G)
names(mtx_G) = c("variable", "G", "RA")
who = c("g__Bacteroides","g__Prevotella","g__Alistipes","g__Parabacteroides","g__Faecalibacterium","g__Roseburia","g__Flavonifractor",
        "g__Akkermansia","g__Escherichia","g__Parasutterella", "g__Barnesiella","g__Clostridium", "g__Lachnoclostridium","g__Odoribacter","Others","unclassified")
mtx_G$variable = ordered(mtx_G$variable, c(sample_order))
mtx_G$G = ordered(mtx_G$G, c(who))
p = ggplot(mtx_G, aes(x = variable, y = RA, fill = G)) + geom_bar(stat="identity")  + scale_fill_manual(values = col)+
        theme_bw(base_size = 8) +  ylab("Read counts") + xlab("cpn60 uniref90")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank() )#+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p$data$G = factor(p$data$G, ordered = TRUE, levels = rev(who))
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mtx_all_ra_cpn_G.pdf", width =5, height = 3)

#mgx_genus level
mgx_G = mgx[,c("G", "variable", "RA")] %>% group_by(variable, G) %>% summarise_all(sum)
mgx_G = as.data.frame(mgx_G)
mgx_G = na.omit(mgx_G)
who = mgx_G[,c("G", "RA")] %>% group_by(G) %>% summarise_all(sum)
who = as.data.frame(who)
who = arrange(who, RA)
who = rev(who$G)[1:15]
mgx_G = mgx_G[mgx_G$G %in% who, ]
mgx_G = mgx_G %>% pivot_wider(names_from = variable, values_from=RA, values_fill = list(ID = 0))
mgx_G = as.data.frame(t(mgx_G))
names(mgx_G) = mgx_G[1,]
mgx_G = mgx_G[-1,]
mgx_G[] <- lapply(mgx_G, function(x) as.numeric(as.character(x)))
mgx_G$Others = 1 -rowSums(mgx_G)
mgx_G$variable = row.names(mgx_G)
mgx_G = melt(mgx_G)
names(mgx_G) = c("variable", "G", "RA")
mgx_G= mgx_G[mgx_G$variable %in% mtx_G$variable,]
who = c("g__Bacteroides", "g__Faecalibacterium", "g__Prevotella", "g__Parabacteroides", "g__Roseburia","g__Alistipes","g__Lachnospiraceae_unclassified","g__Eubacterium","g__Firmicutes_unclassified" ,"g__Escherichia","g__Blautia",
        "g__Dialister", "g__Clostridium" ,"g__Ruminococcaceae_unclassified","Others",  "unclassified" )
col = c("g__Bacteroides"="#84cf26", "g__Faecalibacterium"="#d07b09", "g__Prevotella"="#9ff8e5","g__Parabacteroides"="#657f26","g__Roseburia"="#5e4210",
        "g__Alistipes"="#9c4afd","g__Lachnospiraceae_unclassified"="#58c6f9","g__Eubacterium"="#e887bb","g__Firmicutes_unclassified"="#02ce88", "g__Escherichia"="#dee020","g__Blautia"="#d85b75","g__Dialister"="#0d8561","g__Clostridium"="#c3802e", "g__Ruminococcaceae_unclassified"="#2918cd", 
        "Others"="#fceac5","unclassified"="#f7ae4f" )
mgx_G$variable = ordered(mgx_G$variable, c(sample_order))
mgx_G$G = ordered(mgx_G$G, c(who))
p = ggplot(mgx_G, aes(x = variable, y = RA, fill = G)) + geom_bar(stat="identity")  + scale_fill_manual(values = col)+
        theme_bw(base_size = 8) +  ylab("Relative abundance") + xlab("Samples")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() )
p$data$G = factor(p$data$G, ordered = TRUE, levels = who)
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mgx_801_ra_cpn_G.pdf", width =5.5, height = 3)

##hmp2 MGX species level
mgx = read.csv("MGX/hmp2_melt.csv", sep = ",", header = T)
mgx$variable = gsub("_Abundance-RPKs", "", mgx$variable)
meta = meta[meta$External.ID %in% mgx$variable,]
mgx$collection = meta[match(mgx$variable, meta$External.ID), "IntervalName"]
mgx$sum <- ave(mgx$value, mgx$variable, FUN=sum)
mgx$RA = mgx$value / mgx$sum

sample_order = mgx[,c("variable", "S", "RA")] %>% group_by(variable, S) %>% summarise_all(sum)
sample_order = as.data.frame(sample_order[sample_order$S == "s__Bacteroides_vulgatus",])
sample_order = arrange(sample_order, RA)
sample_order = rev(sample_order$variable)
mgx$variable = ordered(mgx$variable, c(sample_order))
mgx_S = mgx[,c("S", "variable", "RA")] %>% group_by(variable, S) %>% summarise_all(sum)
who = mgx_S[,c("S", "RA")] %>% group_by(S) %>% summarise_all(sum)
who = as.data.frame(who)
who = arrange(who, RA)
who = rev(who$S)[1:25]
mgx_S = mgx_S[mgx_S$S %in% who, ]
mgx_S = mgx_S %>% pivot_wider(names_from = variable, values_from=RA, values_fill = list(ID = 0))
mgx_S = as.data.frame(t(mgx_S))
names(mgx_S) = mgx_S[1,]
mgx_S = mgx_S[-1,]
mgx_S[] <- lapply(mgx_S, function(x) as.numeric(as.character(x)))
mgx_S$Others = 1 -rowSums(mgx_S)
mgx_S$variable = row.names(mgx_S)
mgx_S = melt(mgx_S)
names(mgx_S) = c("variable", "S", "RA")
who = c("s__Bacteroides_vulgatus","s__Bacteroides_dorei","s__Bacteroides_uniformis","s__Bacteroides_stercoris","s__Faecalibacterium_prausnitzii","s__Bacteroides_ovatus", "s__Prevotella_copri",            
        "s__Eubacterium_rectale","s__Parabacteroides_distasonis","s__Hungatella_hathewayi", "s__Bacteroides_thetaiotaomicron", "s__Escherichia_coli","s__Parabacteroides_merdae","s__Alistipes_putredinis","s__Roseburia_intestinalis",      
        "s__Bacteroides_caccae","s__Bacteroides_fragilis","s__Veillonella_parvula","s__Flavonifractor_plautii","s__Alistipes_finegoldii","s__Alistipes_sp_CAG_268","s__Roseburia_faecis",             "s__Bacteroides_intestinalis",    
        "s__Bacteroides_xylanisolvens", "Others",  "unclassified" )
col = c("s__Bacteroides_vulgatus"="#376220","s__Bacteroides_dorei"="#218c9c","s__Bacteroides_uniformis"="#10e6a0","s__Bacteroides_stercoris"="#cddd55","s__Bacteroides_ovatus"="#9cc55d", "s__Bacteroides_thetaiotaomicron"="#6ca9d3","s__Bacteroides_caccae"="#c8fa57","s__Bacteroides_fragilis"="#cddd55","s__Bacteroides_intestinalis"="#d3d27a","s__Bacteroides_xylanisolvens"="#0bcdca", 
        "s__Faecalibacterium_prausnitzii"="#d5239d","s__Prevotella_copri"="#c473dd",            
        "s__Eubacterium_rectale"="#e6cce3","s__Parabacteroides_distasonis"="#5f208b","s__Hungatella_hathewayi"="#8ff8f7", "s__Escherichia_coli"="#d4e923","s__Parabacteroides_merdae"="#760842","s__Alistipes_putredinis"="#f6aa6f","s__Roseburia_intestinalis"="#eed325",      
        "s__Veillonella_parvula"="#424742","s__Flavonifractor_plautii"="#cd9775","s__Alistipes_finegoldii"="#f77b5e","s__Alistipes_sp_CAG_268"="#d67636","s__Roseburia_faecis"="#e8d491",                 
        "Others"="#fceac5","unclassified"="#f7ae4f")
mgx_S$variable = ordered(mgx_S$variable, c(sample_order))
mgx_S$S = ordered(mgx_S$S, c(who))

p = ggplot(mgx_S, aes(x = variable, y = RA, fill = S)) + geom_bar(stat="identity")  + scale_fill_manual(values = col)+
        theme_bw(base_size = 8) +  ylab("Read counts") + xlab("cpn60 uniref90")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() )
p$data$S = factor(p$data$S, ordered = TRUE, levels = who)
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("MGX/mgx_mgx100_ra_cpnUniref90_S.pdf", width =7, height = 8)

##hmp2 MTX species level
mtx = read.csv("hmp2_melt.csv", sep = ",", header = T)
mtx$variable = gsub("_Abundance-RPKs", "", mtx$variable)
meta = meta[meta$External.ID %in% mtx$variable,]
mtx$collection = meta[match(mtx$variable, meta$External.ID), "IntervalName"]
mtx$sum <- ave(mtx$value, mtx$variable, FUN=sum)
mtx$RA = mtx$value / mtx$sum
mtx$variable = ordered(mtx$variable, c(sample_order))
mtx = mtx[!mtx$variable == "CSM67UEW",] ##need to remove this failed sample
mtx = mtx[!mtx$variable == "CSM7KON8",] ##need to remove this failed sample
mtx_S = mtx[,c("S", "variable", "RA")] %>% group_by(variable, S) %>% summarise_all(sum)
mtx_S = as.data.frame(mtx_S)
who = mtx_S[,c("S", "RA")] %>% group_by(S) %>% summarise_all(sum)
who = as.data.frame(who)
who = arrange(who, RA)
who = rev(who$S)[1:25]
mtx_S = mtx_S[mtx_S$S %in% who, ]
mtx_S = mtx_S %>% pivot_wider(names_from = variable, values_from=RA, values_fill = list(ID = 0))
mtx_S = as.data.frame(t(mtx_S))
names(mtx_S) = mtx_S[1,]
mtx_S = mtx_S[-1,]
mtx_S[] <- lapply(mtx_S, function(x) as.numeric(as.character(x)))
mtx_S$Others = 1 -rowSums(mtx_S)
mtx_S$variable = row.names(mtx_S)
mtx_S = melt(mtx_S)
names(mtx_S) = c("variable", "S", "RA")
who = c("s__Bacteroides_vulgatus", "s__Bacteroides_uniformis","s__Bacteroides_stercoris","s__Bacteroides_dorei","s__Prevotella_copri","s__Bacteroides_ovatus","s__Parasutterella_excrementihominis"
        ,"s__Alistipes_putredinis","s__Flavonifractor_plautii"          
        ,"s__Parabacteroides_merdae","s__Alistipes_finegoldii"            
        ,"s__Bacteroides_thetaiotaomicron","s__Akkermansia_muciniphila"         
        ,"s__Parabacteroides_distasonis","s__Faecalibacterium_prausnitzii"    
        ,"s__Hungatella_hathewayi","s__Roseburia_intestinalis"          
        ,"s__Bacteroides_intestinalis","s__Escherichia_coli"                
        ,"s__Bacteroides_caccae","s__Clostridium_neonatale"           
        ,"s__Bacteroides_fragilis","s__Bacteroides_xylanisolvens"       
        ,"s__Dialister_sp_CAG_357","Others",  "unclassified" )
col = c("s__Bacteroides_vulgatus"="#376220","s__Bacteroides_dorei"="#218c9c","s__Bacteroides_uniformis"="#10e6a0","s__Bacteroides_stercoris"="#cddd55","s__Bacteroides_ovatus"="#9cc55d", "s__Bacteroides_thetaiotaomicron"="#6ca9d3","s__Bacteroides_caccae"="#c8fa57","s__Bacteroides_fragilis"="#cddd55","s__Bacteroides_intestinalis"="#d3d27a","s__Bacteroides_xylanisolvens"="#0bcdca", 
        "s__Faecalibacterium_prausnitzii"="#d5239d","s__Prevotella_copri"="#c473dd","s__Akkermansia_muciniphila"="#2822dc",          
        "s__Clostridium_neonatale"="#a45d16","s__Dialister_sp_CAG_357"="#84a5e8",           
        "s__Eubacterium_rectale"="#e6cce3","s__Parabacteroides_distasonis"="#5f208b","s__Parasutterella_excrementihominis"="#be8231","s__Hungatella_hathewayi"="#8ff8f7", "s__Escherichia_coli"="#d4e923","s__Parabacteroides_merdae"="#760842","s__Alistipes_putredinis"="#f6aa6f","s__Roseburia_intestinalis"="#eed325",      
        "s__Veillonella_parvula"="#424742","s__Flavonifractor_plautii"="#cd9775","s__Alistipes_finegoldii"="#f77b5e","s__Alistipes_sp_CAG_268"="#d67636","s__Roseburia_faecis"="#e8d491",                 
        "Others"="#fceac5","unclassified"="#f7ae4f")
mtx_S$variable = ordered(mtx_S$variable, c(sample_order))
mtx_S$S = ordered(mtx_S$S, c(who))

p = ggplot(mtx_S, aes(x = variable, y = RA, fill = S)) + geom_bar(stat="identity")  + scale_fill_manual(values = col)+
        theme_bw(base_size = 8) +  ylab("Read counts") + xlab("cpn60 uniref90")+ theme(legend.text = element_text(face = "italic")) + 
        guides(fill = guide_legend(ncol = 1, reverse=FALSE, keyheight = 0.55)) +
        labs(title="") + theme(plot.title = element_text(size = 20))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank() )
p$data$S = factor(p$data$S, ordered = TRUE, levels = who)
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = sample_order)
ggsave("mtx_mtx98_ra_cpnUniref90_S.pdf", width =7, height = 8)
