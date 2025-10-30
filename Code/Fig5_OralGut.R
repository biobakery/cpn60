library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#setwd("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/HMP2")
setwd("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2")


#species clevelend plot
dotplot_S <- fread("MTX/dotplot_all_S_classified.tsv")
dotplot_S = dotplot_S[,2:10]
dotplot_S[is.na(dotplot_S)] <- 0
dotplot_S$S <- gsub("s__", "", dotplot_S$S)

share <- dotplot_S %>%
        select(S, sample, MGX_RA, MTX_RA) %>%      # keep your sample identifier
        pivot_longer(cols = c(MGX_RA, MTX_RA),
                     names_to  = "library",
                     values_to = "value") %>%
        mutate(library = sub("_RA$", "", library))

oralvsgut <- fread("oralvsgut.csv")
oralvsgut$S <- gsub("s__", "", paste0("s_", oralvsgut$feature))

share$major_site <- oralvsgut$major_site[match(share$S, oralvsgut$feature)]
share <- share[!is.na(share$major_site), ]

abbr <- function(x) sub("([A-Za-z])[^_]*_([^_]*)", "\\1.\\2", x)
share <- share %>%mutate(S = abbr(S))
bugs_list <- fread("bugs_list.csv")[,2:6] %>%  mutate(S = gsub("s__", "", S))
bugs_list <- bugs_list %>%mutate(S = abbr(S))

gut_sp   <- share %>% filter(major_site=="gut")   %>% pull(S) %>% unique()
oral_sp  <- share %>% filter(major_site=="oral")  %>% pull(S) %>% unique()

gut_keep <- bugs_list %>%
        filter(S %in% gut_sp) %>%
        mutate(total = MGX + MTX) %>%
        filter(total != 0, MGX >= 0.01) %>%
        group_by(S) %>%
        summarize(n = n(), .groups = "drop") %>%
        filter(n >= 16) %>%
        pull(S)

oral_keep <- bugs_list %>%
        filter(S %in% oral_sp) %>%
        mutate(total = MGX + MTX) %>%
        filter(total != 0, MGX >= 0.01) %>%
        group_by(S) %>%
        summarize(n = n(), .groups = "drop")  %>%
        pull(S)

share <- share %>%
        filter((major_site=="gut"  & S %in% gut_keep) |
                       (major_site=="oral" & S %in% oral_keep))

sig_tbl <- share %>%
        group_by(major_site, S) %>%
        summarize(
                p_val = wilcox.test(
                        value[library=="MGX"],
                        value[library=="MTX"],
                        alternative = "two.sided"
                )$p.value,
                .groups="drop"
        ) %>%
        mutate(sig = p_val < 0.05)

share <- share %>% left_join(sig_tbl, by=c("major_site","S"))

# 8) Apply your who_mean ordering if you have it
who_mean <- read.csv("Fig5_who_mean.csv", stringsAsFactors=FALSE)[,1] %>%
        gsub("s__", "", .) %>% abbr() %>% unique()
share$S <- factor(share$S, levels=who_mean, ordered=TRUE)

share_means <- share %>%
        group_by(major_site, S, library) %>%
        summarise(mean_value = mean(value), .groups="drop")

# 2) Merge in the per‐species significance flag (sig_tbl from earlier)
#    sig_tbl has: major_site, S, p_val, sig
share_means <- share_means %>%
        left_join(sig_tbl, by = c("major_site","S"))

# 3) Set up plotting function for means
make_mean_cleveland <- function(df){
        ggplot(df, aes(x = mean_value, y = S, color = library, group = S, alpha = sig)) +
                geom_point(size = 3) +
                geom_line() +
                scale_color_manual(values = c("MGX" = "steelblue", "MTX" = "#c22d61")) +
                scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = FALSE) +
                scale_x_log10() +
                coord_cartesian(xlim = c(1e-7, 1)) +
                theme_bw() +
                labs(x = "Mean relative abundance", y = NULL, color = "Dataset")
}

# 4) Split by site and plot
mean_gut  <- share_means %>% filter(major_site=="gut")  %>% 
        mutate(S = factor(S, levels = who_mean, ordered=TRUE)) %>%
        make_mean_cleveland() + ggtitle("Gut")

mean_oral <- share_means %>% filter(major_site=="oral") %>% 
        mutate(S = factor(S, levels = who_mean, ordered=TRUE)) %>%
        make_mean_cleveland() + ggtitle("Oral")

# 5) Combine
(mean_gut / mean_oral) + plot_layout(heights = c(2,1)) # 5 x 12




#dotplot oral vs. gut
dotplot_S$major_site = oralvsgut[match(dotplot_S$S, oralvsgut$S), "major_site"]
dotplot_S = dotplot_S[,2:ncol(dotplot_S)]
meta = read.csv("hmp2_metadata.csv", sep = ",", header = T)
dotplot_S$health = meta[match(dotplot_S$sample, meta$`External.ID`), "diagnosis"]
dotplot_oral = dotplot_S[dotplot_S$major_site == "oral",]
dotplot_gut = dotplot_S[dotplot_S$major_site == "gut",]
dotplot_color = na.omit(dotplot_S)
dotplot_color = dotplot_color[!dotplot_color$major_site == "neither",]

filtered_dotplot_color_nonZero = filter(dotplot_color, MGX_RA !=0, MTX_RA !=0)
#filtered_dotplot_color = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/filtered_dotplot_color_nonZero.csv")

cor_results <- data.frame(
        major_site = character(),
        cor_val = numeric(),
        p_val = numeric(),
        p_adj = numeric(),
        slope = numeric(),
        intercept = numeric(),
        eqn = character(),
        stringsAsFactors = FALSE
)

# Loop through each major_site (oral and gut) to calculate statistics
for (site in unique(filtered_dotplot_color_nonZero$major_site)) {
        subset_data <- filtered_dotplot_color_nonZero %>% filter(major_site == site)
        cor_test <- cor.test(subset_data$MGX_RA, subset_data$MTX_RA, method = "pearson")
        cor_val <- round(cor_test$estimate, 2)
        p_val <- round(cor_test$p.value, 4)
        lm_fit <- lm(MTX_RA ~ MGX_RA, data = subset_data)
        slope <- round(coef(lm_fit)[2], 2)
        intercept <- round(coef(lm_fit)[1], 2)
        p_adj <- p.adjust(p_val, method = "fdr")
        eqn <- paste("y =", slope, "*x +", intercept, "; R² =", cor_val^2, ", p(FDR) =", signif(p_adj, 3))
        cor_results <- rbind(cor_results, data.frame(
                major_site = site,
                cor_val = cor_val,
                p_val = p_val,
                p_adj = p_adj,
                slope = slope,
                intercept = intercept,
                eqn = eqn,
                stringsAsFactors = FALSE
        ))
}

print(cor_results)

col = c("oral" = "#2a64b0", "gut" = "#ebb734")
filtered_dotplot_color_nonZero$V1 <- sub("-s__([A-Za-z]+)_([a-z]+)", "-\\1.\\2", filtered_dotplot_color_nonZero$V1)
filtered_dotplot_color_nonZero <- filtered_dotplot_color_nonZero %>%mutate(V1 = str_replace(V1, "([A-Za-z0-9_]+)-([A-Z])[a-z]+\\.([a-z]+)", "\\1-\\2.\\3"))
p <- ggplot(filtered_dotplot_color_nonZero, aes(x = MGX_RA, y = MTX_RA, color = major_site)) + scale_color_manual(values = col)+
        geom_point(size = 2, alpha = 0.5) +
        scale_y_sqrt() +scale_x_sqrt() + theme_light() +
        geom_smooth(aes(group = major_site), method = lm, se = TRUE, fullrange = TRUE) + 
        labs(title = "Every Species in Every Sample",
             x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") 
p <- p + annotate("text", x = 0.1, y = max(filtered_dotplot_color_nonZero$MTX_RA), 
                  label = cor_results$eqn[cor_results$major_site == "oral"], color = "#2a64b0", hjust = 0) +
        annotate("text", x = 0.1, y = max(filtered_dotplot_color_nonZero$MTX_RA) * 0.9, 
                 label = cor_results$eqn[cor_results$major_site == "gut"], color = "#ebb734", hjust = 0)

#oral_to_label <- filtered_dotplot_color_nonZero %>% dplyr::filter(major_site == "oral" & (MTX_RA > 0.1 & MTX_RA > MGX_RA ))
#p <- p + geom_text_repel(data = oral_to_label,aes(label = V1), size = 3, color = "#2a64b0", max.overlaps = 100)
p
ggsave("dotplot_S_gut&oral_nonZero.pdf", width =6, height = 5)




######2025_03-04. revised the Figure5 as suggested
dotplot_color$diagnosis = meta[match(dotplot_color$sample, meta$External.I), "diagnosis"]

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

species_summary <- dotplot_color %>%
        filter(MGX_RA > 0 | MTX_RA > 0) %>%  # Select non-zero species
        group_by(diagnosis, major_site) %>%
        summarise(
                num_species_MGX = n_distinct(S[MGX_RA > 0]),
                num_species_MTX = n_distinct(S[MTX_RA > 0]),
                .groups = "drop"
        ) %>%
        pivot_longer(cols = c(num_species_MGX, num_species_MTX), 
                     names_to = "Dataset", values_to = "num_species") %>%
        mutate(Dataset = recode(Dataset, 
                                "num_species_MGX" = "MGX", 
                                "num_species_MTX" = "MTX"))  # Rename for clarity

# Create individual plots for each subset
#species_summary = fread("species_summary.csv")
species_summary$diagnosis = ordered(species_summary$diagnosis, c("UC", "CD", "nonIBD"))

p1 <- ggplot(species_summary %>% filter(major_site == "gut", Dataset == "MGX"), 
             aes(x = diagnosis, y = num_species, fill = diagnosis)) +
        geom_bar(stat = "identity") + ylim(0, 125)+
        labs(title = "Gut Species in MGX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p2 <- ggplot(species_summary %>% filter(major_site == "gut", Dataset == "MTX"), 
             aes(x = diagnosis, y = num_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +ylim(0, 125)+
        labs(title = "Gut Species in MTX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p3 <- ggplot(species_summary %>% filter(major_site == "oral", Dataset == "MGX"), 
             aes(x = diagnosis, y = num_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +ylim(0, 125)+
        labs(title = "Oral Species in MGX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p4 <- ggplot(species_summary %>% filter(major_site == "oral", Dataset == "MTX"), 
             aes(x = diagnosis, y = num_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +ylim(0, 125)+
        labs(title = "Oral Species in MTX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
combined_plot <- (p1 | p3) / (p2 | p4) 
combined_plot #export in 10 x 3 inch, as species_summary.pdf

## Plot the average of species in each sample in each dataset, since the sample numbers are different in CD vs. UC vs. control
species_summary <- dotplot_color %>% filter(MGX_RA > 0 | MTX_RA > 0) %>%  # Select non-zero species
        group_by(diagnosis, major_site, sample) %>%  # Group by sample within each group
        summarise(
                num_species_MGX = sum(MGX_RA > 0),  # Count non-zero species in MGX for each sample
                num_species_MTX = sum(MTX_RA > 0),  # Count non-zero species in MTX for each sample
                .groups = "drop"
        ) %>%
        pivot_longer(cols = c(num_species_MGX, num_species_MTX), 
                     names_to = "Dataset", values_to = "num_species") %>%
        mutate(Dataset = recode(Dataset, 
                                "num_species_MGX" = "MGX", 
                                "num_species_MTX" = "MTX"))  # Rename for clarity

# Calculate mean and SEM (Standard Error of the Mean)
species_summary_avg <- species_summary %>%
        group_by(diagnosis, major_site, Dataset) %>%
        summarise(
                mean_species = mean(num_species),  # Mean number of species per sample
                sem_species = sd(num_species) / sqrt(n()),  # SEM calculation
                .groups = "drop"
        )

#species_summary_avg = fread("species_summary_avg.csv")
species_summary_avg$diagnosis = ordered(species_summary_avg$diagnosis, c("UC", "CD", "nonIBD"))

p1 <- ggplot(species_summary_avg %>% filter(major_site == "gut", Dataset == "MGX"), 
             aes(x = diagnosis, y = mean_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_species - sem_species, ymax = mean_species + sem_species), width = 0.2) + ylim(0,25)+
        labs(title = "Gut Species in MGX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p2 <- ggplot(species_summary_avg %>% filter(major_site == "gut", Dataset == "MTX"), 
             aes(x = diagnosis, y = mean_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_species - sem_species, ymax = mean_species + sem_species), width = 0.2) + ylim(0,25)+
        labs(title = "Gut Species in MTX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()

p3 <- ggplot(species_summary_avg %>% filter(major_site == "oral", Dataset == "MGX"), 
             aes(x = diagnosis, y = mean_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_species - sem_species, ymax = mean_species + sem_species), width = 0.2) + ylim(0,3)+
        labs(title = "Oral Species in MGX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()

p4 <-ggplot(species_summary_avg %>% filter(major_site == "oral", Dataset == "MTX"), 
             aes(x = diagnosis, y = mean_species, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_species - sem_species, ymax = mean_species + sem_species), width = 0.2) + ylim(0,3)+
        labs(title = "Oral Species in MTX", x = "Diagnosis", y = "Number of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()

combined_plot <- (p1 | p3) / (p2 | p4) 
combined_plot #export in 10 x 3 inch, as species_summary_avg.pdf
ggsave("species_summary_avg.pdf", width =10, height = 3)

# Summarize the number of species where MGX_RA > MTX_RA and MGX_RA < MTX_RA per diagnosis and major_site
library(dplyr)
species_comparison <- dotplot_color %>%
        filter(!(MGX_RA == 0 & MTX_RA == 0)) %>%  # Exclude species with zero abundance in both MGX and MTX
        group_by(diagnosis, major_site) %>%
        summarise(
                total_species = n(),  # Total species (sum of every species in every sample) in each group
                MGX_greater_MTX = sum(MGX_RA > MTX_RA) / total_species * 100,  # % of species where MGX > MTX
                MGX_less_MTX = sum(MGX_RA < MTX_RA) / total_species * 100,  # % of species where MGX < MTX
                .groups = "drop"
        ) %>%
        select(diagnosis, major_site, MGX_greater_MTX, MGX_less_MTX)  # Keep only relevant columns

# Print the cleaned summary table
print(species_comparison)
#species_comparison = fread("species_comparison.csv")
species_comparison$diagnosis = ordered(species_comparison$diagnosis, c("UC", "CD", "nonIBD"))

p1 <- ggplot(species_comparison %>% filter(major_site == "gut"), 
             aes(x = diagnosis, y = MGX_greater_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity")  + ylim(0,100)+
        labs(title = "Gut Species MGX>MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p2 <- ggplot(species_comparison %>% filter(major_site == "gut"), 
             aes(x = diagnosis, y = MGX_less_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity")  + ylim(0,100)+
        labs(title = "Gut Species MGX<MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p3 <- ggplot(species_comparison %>% filter(major_site == "oral"), 
             aes(x = diagnosis, y = MGX_greater_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity")  + ylim(0,100)+
        labs(title = "Oral Species MGX>MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p4 <- ggplot(species_comparison %>% filter(major_site == "oral"), 
             aes(x = diagnosis, y = MGX_less_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity")  + ylim(0,100)+
        labs(title = "Oral Species MGX<MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
combined_plot <- (p1 | p3) / (p2 | p4) 
combined_plot #export in 10 x 3 inch, as species_summary_avg.pdf

#in each sample
species_comparison <- dotplot_color %>%
        filter(!(MGX_RA == 0 & MTX_RA == 0)) %>%  # Remove completely absent species
        group_by(diagnosis, major_site, sample) %>%  # Group at the sample level
        summarise(
                total_species = n(),  # Total species in each sample
                MGX_greater_MTX = sum(MGX_RA > MTX_RA) / total_species * 100,  # % per sample where MGX > MTX
                MGX_less_MTX = sum(MGX_RA < MTX_RA) / total_species * 100,  # % per sample where MGX < MTX
                .groups = "drop"
        ) %>%
        group_by(diagnosis, major_site) %>%  # Aggregate at the group level
        summarise(
                mean_MGX_greater_MTX = mean(MGX_greater_MTX),  # Average % per sample
                sem_MGX_greater_MTX = sd(MGX_greater_MTX) / sqrt(n()),  # SEM for MGX > MTX
                mean_MGX_less_MTX = mean(MGX_less_MTX),  # Average % per sample
                sem_MGX_less_MTX = sd(MGX_less_MTX) / sqrt(n()),  # SEM for MGX < MTX
                .groups = "drop"
        )
#species_comparison_avg = fread("species_comparison_avg.csv")
species_comparison_avg$diagnosis = ordered(species_comparison_avg$diagnosis, c("UC", "CD", "nonIBD"))

p1 <- ggplot(species_comparison_avg %>% filter(major_site == "gut"), 
             aes(x = diagnosis, y = mean_MGX_greater_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_MGX_greater_MTX - sem_MGX_greater_MTX, ymax = mean_MGX_greater_MTX + sem_MGX_greater_MTX), width = 0.2) + 
        labs(title = "Gut Species MGX > MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p2 <- ggplot(species_comparison_avg %>% filter(major_site == "gut"), 
             aes(x = diagnosis, y = mean_MGX_less_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_MGX_less_MTX - sem_MGX_less_MTX, ymax = mean_MGX_less_MTX + sem_MGX_less_MTX), width = 0.2) + 
        labs(title = "Gut Species MGX < MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()

p3 <- ggplot(species_comparison_avg %>% filter(major_site == "oral"), 
             aes(x = diagnosis, y = mean_MGX_greater_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_MGX_greater_MTX - sem_MGX_greater_MTX, ymax = mean_MGX_greater_MTX + sem_MGX_greater_MTX), width = 0.2) + 
        labs(title = "Oral Species MGX > MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()
p4 <- ggplot(species_comparison_avg %>% filter(major_site == "oral"), 
             aes(x = diagnosis, y = mean_MGX_less_MTX, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = mean_MGX_less_MTX - sem_MGX_less_MTX, ymax = mean_MGX_less_MTX + sem_MGX_less_MTX), width = 0.2) + 
        labs(title = "Oral Species MGX < MTX", x = "Diagnosis", y = "% of Species") +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_light()+ coord_flip()

combined_plot <- (p1 | p3) / (p2 | p4) 
combined_plot #export in 10 x 3 inch, as species_comparison_avg.pdf

# Print the summary table with averages and SEM
print(species_comparison)


# non-zero dotplot of oral vs. gut species
library("ggpubr")
filtered_dotplot_color = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/filtered_dotplot_color_nonZero.csv")
filtered_dotplot_color$V1 = NULL
filtered_dotplot_color$diagnosis = filtered_dotplot_color$health
meta = fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/hmp2_metadata.csv")
library(dplyr)

# Identify species that are non-zero in at least one sample in both nonIBD and CD/UC
species_in_both_groups <- filtered_dotplot_color %>%
        group_by(S, diagnosis) %>%
        summarise(
                non_zero_count = sum(MGX_RA > 0 | MTX_RA > 0),  # Count non-zero occurrences in MGX or MTX
                .groups = "drop"
        ) %>%
        filter(diagnosis %in% c("CD", "UC", "nonIBD")) %>%  # Keep only relevant diagnoses
        group_by(S) %>%
        summarise(
                present_in_CD = any(diagnosis == "CD" & non_zero_count > 0),
                present_in_UC = any(diagnosis == "UC" & non_zero_count > 0),
                present_in_nonIBD = any(diagnosis == "nonIBD" & non_zero_count > 0),
                .groups = "drop"
        ) %>%
        filter((present_in_CD | present_in_UC) & present_in_nonIBD) %>%  # Keep species present in both case and control
        pull(S)  # Extract species names

filtered_species_table <- filtered_dotplot_color %>%
        filter(S %in% species_in_both_groups)

#remove any 0
species_in_all_groups <- filtered_dotplot_color %>% group_by(S, diagnosis) %>%
        summarise(non_zero_MGX = any(MGX_RA > 0),   non_zero_MTX = any(MTX_RA > 0),  .groups = "drop") %>%
        filter(diagnosis %in% c("CD", "UC", "nonIBD")) %>%  # Keep only relevant diagnoses
        group_by(S) %>% summarise(
                present_in_CD_MGX = any(diagnosis == "CD" & non_zero_MGX),
                present_in_CD_MTX = any(diagnosis == "CD" & non_zero_MTX),
                present_in_UC_MGX = any(diagnosis == "UC" & non_zero_MGX),
                present_in_UC_MTX = any(diagnosis == "UC" & non_zero_MTX),
                present_in_nonIBD_MGX = any(diagnosis == "nonIBD" & non_zero_MGX),
                present_in_nonIBD_MTX = any(diagnosis == "nonIBD" & non_zero_MTX),
                .groups = "drop" ) %>%
        filter(
                (present_in_CD_MGX & present_in_CD_MTX) | (present_in_UC_MGX & present_in_UC_MTX),  # Must be present in both MGX and MTX in at least one case group (CD or UC)
                present_in_nonIBD_MGX & present_in_nonIBD_MTX  # Must be present in both MGX and MTX in nonIBD
        ) %>%
        pull(S)  # Extract species names

filtered_species_table_nonZero <- filtered_dotplot_color %>%
        filter(S %in% species_in_all_groups)

#filter the filtered_dotplot_color
#filtered_dotplot_color = filtered_dotplot_color[filtered_dotplot_color$S %in% filtered_species_table$S,]
#filtered_dotplot_color_nonZero = filter(filtered_dotplot_color, MGX_RA > 0.01 )
#filtered_dotplot_color_nonZero = filtered_dotplot_color
filtered_dotplot_color_nonZero = filter(filtered_dotplot_color_nonZero, health == "CD")


# Initialize an empty data frame to store results
cor_results <- data.frame(
        major_site = character(),
        cor_val = numeric(),
        p_val = numeric(),
        p_adj = numeric(),
        slope = numeric(),
        intercept = numeric(),
        eqn = character(),
        stringsAsFactors = FALSE
)

# Loop through each major_site (oral and gut) to calculate statistics
for (site in unique(filtered_dotplot_color_nonZero$major_site)) {
        subset_data <- filtered_dotplot_color_nonZero %>% filter(major_site == site)
        cor_test <- cor.test(subset_data$MGX_RA, subset_data$MTX_RA, method = "pearson")
        cor_val <- round(cor_test$estimate, 2)
        p_val <- round(cor_test$p.value, 4)
        lm_fit <- lm(MTX_RA ~ MGX_RA, data = subset_data)
        slope <- round(coef(lm_fit)[2], 2)
        intercept <- round(coef(lm_fit)[1], 2)
        p_adj <- p.adjust(p_val, method = "fdr")
        eqn <- paste("y =", slope, "*x +", intercept, "; R² =", cor_val^2, ", p(FDR) =", signif(p_adj, 3))
        cor_results <- rbind(cor_results, data.frame(
                major_site = site,
                cor_val = cor_val,
                p_val = p_val,
                p_adj = p_adj,
                slope = slope,
                intercept = intercept,
                eqn = eqn,
                stringsAsFactors = FALSE
        ))
}

# Print results
print(cor_results)

col = c("oral" = "#2a64b0", "gut" = "#ebb734")
filtered_dotplot_color_nonZero$V1 <- sub("-s__([A-Za-z]+)_([a-z]+)", "-\\1.\\2", filtered_dotplot_color_nonZero$V1)
filtered_dotplot_color_nonZero <- filtered_dotplot_color_nonZero %>%mutate(V1 = str_replace(V1, "([A-Za-z0-9_]+)-([A-Z])[a-z]+\\.([a-z]+)", "\\1-\\2.\\3"))
p <- ggplot(filtered_dotplot_color_nonZero, aes(x = MGX_RA, y = MTX_RA, color = major_site)) + scale_color_manual(values = col)+
        geom_point(size = 2, alpha = 0.5) +
        scale_y_sqrt() +scale_x_sqrt() + theme_light() +
        geom_smooth(aes(group = major_site), method = lm, se = TRUE, fullrange = TRUE) + 
        labs(title = "Every Species in Every Sample",
             x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") 
p <- p + annotate("text", x = 0.1, y = max(filtered_dotplot_color_nonZero$MTX_RA), 
                  label = cor_results$eqn[cor_results$major_site == "oral"], color = "#2a64b0", hjust = 0) +
        annotate("text", x = 0.1, y = max(filtered_dotplot_color_nonZero$MTX_RA) * 0.9, 
                 label = cor_results$eqn[cor_results$major_site == "gut"], color = "#ebb734", hjust = 0)

#oral_to_label <- filtered_dotplot_color_nonZero %>% dplyr::filter(major_site == "oral" & (MTX_RA > 0.1 & MTX_RA > MGX_RA ))
#p <- p + geom_text_repel(data = oral_to_label,aes(label = V1), size = 3, color = "#2a64b0", max.overlaps = 100)
p
ggsave("dotplot_S_gut&oral_nonZero.pdf", width =5, height = 10)

library(ggplot2)
library(dplyr)

# Define custom colors for health conditions
col <- c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")
fill_col <- c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")  # Same colors for shading

# Function to compute correlation statistics for a given subset
compute_correlations <- function(data, site_label) {
        cor_results <- data.frame(
                health = character(),
                cor_val = numeric(),
                p_val = numeric(),
                p_adj = numeric(),
                slope = numeric(),
                intercept = numeric(),
                eqn = character(),
                stringsAsFactors = FALSE
        )
        
        for (group in unique(data$health)) {
                subset_data <- data %>% filter(health == group)
                
                if (nrow(subset_data) > 2) {  # Ensure enough points for correlation
                        cor_test <- cor.test(subset_data$MGX_RA, subset_data$MTX_RA, method = "spearman")
                        cor_val <- round(cor_test$estimate, 2)
                        p_val <- round(cor_test$p.value, 4)
                        lm_fit <- lm(MTX_RA ~ MGX_RA, data = subset_data)
                        slope <- round(coef(lm_fit)[2], 2)
                        intercept <- round(coef(lm_fit)[1], 2)
                        
                        p_adj <- p.adjust(p_val, method = "fdr")
                        eqn <- paste("y =", slope, "*x +", intercept, "; R² =", round(cor_val^2, 3), ", p(FDR) =", signif(p_adj, 3))
                        cor_results <- rbind(cor_results, data.frame(
                                health = group,
                                cor_val = cor_val,
                                p_val = p_val,
                                p_adj = p_adj,
                                slope = slope,
                                intercept = intercept,
                                eqn = eqn,
                                stringsAsFactors = FALSE
                        ))
                }
        }
        return(cor_results)
}

# Function to get species with ≥2 samples per health category
filter_common_species_all_groups <- function(df) {
        df <- as.data.frame(df)  # convert from data.table if needed
        species_to_keep <- df %>%
                group_by(S, health) %>%
                summarize(n_samples = n_distinct(sample), .groups = "drop") %>%
                filter(n_samples >= 3) %>%
                group_by(S) %>%
                summarize(n_groups = n_distinct(health), .groups = "drop") %>%
                filter(n_groups == 3) %>%   # must appear in all 3 groups
                pull(S)
        df %>% filter(S %in% species_to_keep)
}

# Create gut and oral subsets
#filtered_dotplot_color = filtered_dotplot_color %>%filter(MGX_RA != 0, MTX_RA != 0) #add this if plotting the non-zero version
gut <- filtered_dotplot_color_nonZero %>% filter(major_site == "gut") %>% filter(MGX_RA >0 & MTX_RA > 0)
oral <- filtered_dotplot_color_nonZero %>% filter(major_site == "oral")%>% filter(MGX_RA >0 & MTX_RA > 0)

gut <- filter_common_species_all_groups(gut)
oral <- filter_common_species_all_groups(oral)

# Compute correlations for gut and oral
cor_results_gut <- compute_correlations(gut, "Gut")
cor_results_oral <- compute_correlations(oral, "Oral")

# Function to create scatter plot with linear regression and shading
create_plot <- function(data, cor_results, site_label) {
        p <- ggplot(data, aes(x = MGX_RA, y = MTX_RA, color = health, fill = health)) + 
                scale_color_manual(values = col) +  # Apply custom colors
                scale_fill_manual(values = fill_col) +  # Fill colors for shading
                geom_point(size = 2, alpha = 0.5) +  # Scatter points
                scale_y_log10() + scale_x_log10() +  # Square root scaling for better visualization
                theme_light() +
                geom_smooth(aes(group = health), method = lm, se = TRUE, fullrange = TRUE, alpha = 0.3) +  # Regression with shading
                labs(title = paste("Every Species in Every Sample -", site_label),
                     x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") +
                theme(legend.title = element_blank())
        
        # Add annotations for each health condition
        for (i in 1:nrow(cor_results)) {
                p <- p + annotate("text", x = 0.1, y = max(data$MTX_RA) * (0.9 - 0.1 * i),
                                  label = cor_results$eqn[i], color = col[cor_results$health[i]], hjust = 0)
        }
        
        return(p)
}

# Generate plots for gut and oral with shaded regression lines
p_gut <- create_plot(gut, cor_results_gut, "Gut")
p_oral <- create_plot(oral, cor_results_oral, "Oral")

combined_plot = p_gut / p_oral
combined_plot
ggsave("dotplot_S_gut&oral_filterSpecies.pdf", width =5, height = 10)

#another way to compare, RvsD
gut$RvsD = gut$MTX_RA / gut$MGX_RA
gut <- gut %>%mutate(RvsD = ifelse(is.infinite(RvsD), 1e12, RvsD))

gut_RvsD <- gut %>% group_by(diagnosis) %>% summarise(total_cases = n(), pct_RvsD_greater_1 = sum(RvsD > 1) / total_cases * 100, pct_RvsD_less_1 = sum(RvsD < 1) / total_cases * 100, 
                                                          .groups = "drop" )
print(gut_RvsD)

oral$RvsD = oral$MTX_RA / oral$MGX_RA
oral <- oral %>%mutate(RvsD = ifelse(is.infinite(RvsD), 1e12, RvsD))
oral_RvsD <- oral %>% group_by(diagnosis) %>% summarise(total_cases = n(), pct_RvsD_greater_1 = sum(RvsD > 1) / total_cases * 100, pct_RvsD_less_1 = sum(RvsD < 1) / total_cases * 100, 
                                                      .groups = "drop" )
print(oral_RvsD)

##for single species in gut or oral
dotplot_Bvulgatus = filtered_dotplot_color[filtered_dotplot_color$S=="s__Bacteroides_vulgatus",]
dotplot_Fprausnitzii= filtered_dotplot_color[filtered_dotplot_color$S=="s__Faecalibacterium_prausnitzii",]
dotplot_Kpneumonia = filtered_dotplot_color[filtered_dotplot_color$S=="s__Klebsiella_pneumoniae",]
dotplot_Hparainfluenzae = filtered_dotplot_color[filtered_dotplot_color$S=="s__Haemophilus_parainfluenzae",]
dotplot_Dmossii = filtered_dotplot_color[filtered_dotplot_color$S=="s__Dysgonomonas_mossii",]
dotplot_Sthermophilus = filtered_dotplot_color[filtered_dotplot_color$S=="s__Streptococcus_thermophilus",]
dotplot_Bdentium = filtered_dotplot_color[filtered_dotplot_color$S=="s__Bifidobacterium_dentium",]
dotplot_Hpattmaniae = filtered_dotplot_color[filtered_dotplot_color$S=="s__Haemophilus_pittmaniae",]
dotplot_Pmerdae = filtered_dotplot_color[filtered_dotplot_color$S=="s__Parabacteroides_merdae",]

cor_results_Bvulgatus <- compute_correlations(dotplot_Bvulgatus, "B. vulgatus")
cor_results_Fprausnitzii <- compute_correlations(dotplot_Fprausnitzii, "F. prausnitzii")
cor_results_Kpneumonia <- compute_correlations(dotplot_Kpneumonia, "K. pneumoniae")
cor_results_Hparainfluenzae <- compute_correlations(dotplot_Hparainfluenzae, "H. parainfluenzae")
cor_results_Dmossii <- compute_correlations(dotplot_Dmossii, "D. mossii")
cor_results_Sthermophilus <- compute_correlations(dotplot_Sthermophilus, "S. thermophilus")
cor_results_Bdentium <- compute_correlations(dotplot_Bdentium, "B. dentium")
cor_results_Hpattmaniae <- compute_correlations(dotplot_Hpattmaniae, "H. pittmaniae")
cor_results_Pmerdae <- compute_correlations(dotplot_Pmerdae, "P. merdae")

# Generate scatter plots for each species
p_Bvulgatus <- create_plot(dotplot_Bvulgatus, cor_results_Bvulgatus, "B. vulgatus")
p_Fprausnitzii <- create_plot(dotplot_Fprausnitzii, cor_results_Fprausnitzii, "F. prausnitzii")
p_Kpneumonia <- create_plot(dotplot_Kpneumonia, cor_results_Kpneumonia, "K. pneumoniae")
p_Hparainfluenzae <- create_plot(dotplot_Hparainfluenzae, cor_results_Hparainfluenzae, "H. parainfluenzae")
p_Dmossii <- create_plot(dotplot_Dmossii, cor_results_Dmossii, "D. mossii")
p_Sthermophilus <- create_plot(dotplot_Sthermophilus, cor_results_Sthermophilus, "S. thermophilus")
p_Bdentium <- create_plot(dotplot_Bdentium, cor_results_Bdentium, "B. dentium")
p_Hpattmaniae <- create_plot(dotplot_Hpattmaniae, cor_results_Hpattmaniae, "H. pittmaniae")
p_Pmerdae <- create_plot(dotplot_Pmerdae, cor_results_Pmerdae, "P. merdae")


dotplot_Dgadei = filtered_dotplot_color_nonZero[filtered_dotplot_color_nonZero$S=="s__Dysgonomonas_gadei",]
cor_results_Dgadei <- compute_correlations(dotplot_Dgadei, "Dgadei")
p_Dgadei <- create_plot(dotplot_Dgadei, cor_results_Dgadei, "Dgadei")
p_Dgadei
#calculate the cor_results for all the species
library(ggplot2)
library(dplyr)

# Function to compute Spearman correlation statistics for all species
compute_correlations <- function(data, species_name) {
        cor_results <- data.frame(
                species = character(),
                health = character(),
                cor_val = numeric(),
                p_val = numeric(),
                p_adj = numeric(),
                slope = numeric(),
                intercept = numeric(),
                eqn = character(),
                stringsAsFactors = FALSE
        )
        
        for (group in unique(data$health)) {
                subset_data <- data %>% filter(health == group)
                
                if (nrow(subset_data) > 2) {  # Ensure enough data points for correlation
                        cor_test <- cor.test(subset_data$MGX_RA, subset_data$MTX_RA, method = "spearman")
                        cor_val <- round(cor_test$estimate, 2)
                        p_val <- round(cor_test$p.value, 4)
                        
                        lm_fit <- lm(MTX_RA ~ MGX_RA, data = subset_data)
                        slope <- round(coef(lm_fit)[2], 2)
                        intercept <- round(coef(lm_fit)[1], 2)
                        
                        p_adj <- p.adjust(p_val, method = "fdr")
                        
                        eqn <- paste("y =", slope, "*x +", intercept, "; R² =", round(cor_val^2, 3), ", p(FDR) =", signif(p_adj, 3))
                        
                        cor_results <- rbind(cor_results, data.frame(
                                species = species_name,  # Ensure species name is included
                                health = group,
                                cor_val = cor_val,
                                p_val = p_val,
                                p_adj = p_adj,
                                slope = slope,
                                intercept = intercept,
                                eqn = eqn,
                                stringsAsFactors = FALSE
                        ))
                }
        }
        return(cor_results)
}

# Initialize an empty table to store all species results
all_species_results <- data.frame()

# Loop through each species in filtered_dotplot_color$S
for (species in unique(filtered_dotplot_color$S)) {
        species_data <- filtered_dotplot_color %>% filter(S == species)
        
        if (nrow(species_data) > 0) {
                cor_results_species <- compute_correlations(species_data, species)  # Now species names are properly stored
                all_species_results <- rbind(all_species_results, cor_results_species)
        }
}

# Print the final combined results table with real species names
print(all_species_results)

# Save the results to a CSV file
write.csv(all_species_results, "species_correlation_results.csv", row.names = FALSE)

library(patchwork)
layout <- c(
        area(1, 1, 2, 1),  # p1 (top-left, column 1)
        area(1, 2, 2, 5),  # p_gut (top-right, columns 2-5)
        area(3, 1, 4, 1),  # p2 (bottom-left, column 1)
        area(3, 2, 4, 5)   # p_oral (bottom-right, columns 2-5)
)

combined_p <- p1 + p_gut + p2 + p_oral + plot_layout(design = layout)
print(combined_p)

############Other vairables

library(data.table)  # for fread()
library(janitor)     # for clean_names()
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

filtered_dotplot_color <-  fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/filtered_dotplot_color_nonZero.csv") %>%
        as.data.frame() %>%            # if you prefer tibbles later you could tibble::as_tibble()
        select(-V1) %>%                # drop the auto‐created V1
        mutate(diagnosis = health) %>% # duplicate health → diagnosis
        unique()                       # remove exact row duplicates

meta <- fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/hmp2_metadata.csv") %>%
        as.data.frame()

meta_subset <-  meta %>%
        select(`External ID`, 28:ncol(meta)) %>%
        clean_names() %>%               # make all names snake_case & unique
        rename(sample = external_id)    # match the join key

filtered_dotplot_color <- filtered_dotplot_color %>% left_join(meta_subset, by = "sample", suffix = c("", "_meta"))

if (anyDuplicated(names(filtered_dotplot_color)) > 0) {
        stop("There are still duplicated column names!")
}


filtered_dotplot_color$Participant = meta[match(filtered_dotplot_color$sample, meta$`External ID`), "Participant ID"]
filtered_dotplot_color$date = meta[match(filtered_dotplot_color$sample, meta$`External ID`), "date_of_receipt"]
filtered_dotplot_color$Participant_date = paste(filtered_dotplot_color$Participant, filtered_dotplot_color$date, sep = "_")
meta$Participant_date = paste(meta$`Participant ID`, meta$date_of_receipt, sep = "_")
filtered_dotplot_color$Cipro = meta[match(filtered_dotplot_color$Participant_date, meta$Participant_date), "Cipro (Ciprofloxin)"]
filtered_dotplot_color <- filtered_dotplot_color %>%select(where(~ !(all(is.na(.)) || (is.character(.) && all(. == ""))))) ##reduce number of columns by removing the all "NA" and all empty columns
filtered_dotplot_color = unique(filtered_dotplot_color)

###Maasline3
library(maaslin3)
library(tibble)
library(dplyr)
library(tidyr)

filtered_dotplot_color <-  fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/filtered_dotplot_color_nonZero.csv") %>%
        as.data.frame() %>%            # if you prefer tibbles later you could tibble::as_tibble()
        select(-V1) %>%                # drop the auto‐created V1
        mutate(diagnosis = health) %>% # duplicate health → diagnosis
        unique()                       # remove exact row duplicates

meta <- fread("/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/hmp2_metadata.csv") %>%as.data.frame()

meta_subset <-  meta %>%
        select(`External ID`,`Participant ID`, 28:ncol(meta)) %>%
        clean_names() %>%               # make all names snake_case & unique
        rename(sample = external_id)    # match the join key
meta_subset = meta_subset[meta_subset$sample %in% filtered_dotplot_color$sample,]
meta_subset <- meta_subset %>%select(where(~ !(all(is.na(.)) || (is.character(.) && all(. == ""))))) ##reduce number of columns by removing the all "NA" and all empty columns

filtered_dotplot_color$Participant = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "participant_id"]
filtered_dotplot_color$Age = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "age_at_diagnosis"]
filtered_dotplot_color$race = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "race"]
filtered_dotplot_color$general_wellbeing = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "general_wellbeing"]
filtered_dotplot_color$sex = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "sex"]

filtered_dotplot_color$Education = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "education_level"]
filtered_dotplot_color$Occupation = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "occupation"]
filtered_dotplot_color$appendectomy = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "has_the_subject_had_an_appendectomy"]
filtered_dotplot_color$tonsillectomy = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "has_the_subject_had_a_tonsillectomy"]
filtered_dotplot_color$water = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "water"]
filtered_dotplot_color$alcohol = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "alcohol_beer_brandy_spirits_hard_liquor_wine_aperitif_etc"]
filtered_dotplot_color$yogurt = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "yogurt_or_other_foods_containing_active_bacterial_cultures_kefir_sauerkraut"]
filtered_dotplot_color$dairy = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "dairy_milk_cream_ice_cream_cheese_cream_cheese"]
filtered_dotplot_color$probiotic = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "probiotic"]
filtered_dotplot_color$fruits = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "fruits_no_juice_apples_raisins_bananas_oranges_strawberries_blueberries"]
filtered_dotplot_color$vegetables = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "vegetables_salad_tomatoes_onions_greens_carrots_peppers_green_beans_etc"]
filtered_dotplot_color$beans = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "beans_tofu_soy_soy_burgers_lentils_mexican_beans_lima_beans_etc"]
filtered_dotplot_color$wheat = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "whole_grains_wheat_oats_brown_rice_rye_quinoa_wheat_bread_wheat_pasta"]
filtered_dotplot_color$starch = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "starch_white_rice_bread_pizza_potatoes_yams_cereals_pancakes_etc"]
filtered_dotplot_color$eggs = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "eggs"]
filtered_dotplot_color$processed_meat= meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "processed_meat_other_red_or_white_meat_such_as_lunch_meat_ham_salami_bologna"]
filtered_dotplot_color$red_meat = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "red_meat_beef_hamburger_pork_lamb"]
filtered_dotplot_color$shellfish = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "shellfish_shrimp_lobster_scallops_etc"]
filtered_dotplot_color$fish = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "fish_fish_nuggets_breaded_fish_fish_cakes_salmon_tuna_etc" ]
filtered_dotplot_color$sweets = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "sweets_pies_jam_chocolate_cake_cookies_etc"]
filtered_dotplot_color$tea_or_coffee = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "tea_or_coffee_no_sugar_and_no_sugar_replacement"]

filtered_dotplot_color$preemie = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "were_you_born_prematurely_more_than_3_weeks_early"]
filtered_dotplot_color$daycare = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "did_you_attend_daycare_as_a_child"]
filtered_dotplot_color$cigarette_child = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "were_you_exposed_to_cigarette_smoke_as_a_child"]
filtered_dotplot_color$farm = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "did_you_grow_up_on_a_farm"]

filtered_dotplot_color$marijuana = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "do_you_currently_smoke_marijuana"]

filtered_dotplot_color$irritable_bowel_syndrome = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "irritable_bowel_syndrome"]
filtered_dotplot_color$T1D = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "type_i_diabetes_juvenile_diabetes"]
filtered_dotplot_color$rheumatoid_arthritis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "rheumatoid_arthritis"]

filtered_dotplot_color$cancer_breast = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_breast"]
filtered_dotplot_color$cancer_cholangiocarcinoma = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_cholangiocarcinoma"]
filtered_dotplot_color$cancer_colon_or_rectum = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_colon_or_rectum"]
filtered_dotplot_color$cancer_hodgkins_lymphoma = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_hodgkins_lymphoma"]
filtered_dotplot_color$cancer_liver = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_liver"]
filtered_dotplot_color$cancer_lung = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_lung"]
filtered_dotplot_color$lymphoma = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_lymphoma_not_otherwise_specified"]
filtered_dotplot_color$cancer_ovarian = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_ovarian"]
filtered_dotplot_color$cancer_prostate = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "cancer_prostate"]
filtered_dotplot_color$celiac_sprue = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "celiac_sprue"]
filtered_dotplot_color$chronic_bronchitis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "chronic_bronchitis"]
filtered_dotplot_color$familial_mediterranean_fever = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "familial_mediterranean_fever"]
filtered_dotplot_color$graves_disease = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "graves_disease"]
filtered_dotplot_color$guillian_barre_syndrome = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "guillian_barre_syndrome"]
filtered_dotplot_color$thyroiditis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "hashimotos_autoimmune_thyroiditis"]
filtered_dotplot_color$idiopathic_pulmonary_fibrosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "idiopathic_pulmonary_fibrosis"]
filtered_dotplot_color$idiopathic_thrombocytopenia_purpura = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "idiopathic_thrombocytopenia_purpura"]
filtered_dotplot_color$alopecia_areata = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "alopecia_areata"]
filtered_dotplot_color$multiple_sclerosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "multiple_sclerosis"]
filtered_dotplot_color$myasthenia_gravis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "myasthenia_gravis"]
filtered_dotplot_color$myocarditis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "myocarditis"]
filtered_dotplot_color$pericarditis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "pericarditis"]
filtered_dotplot_color$pemphigus_vulgaris = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "pemphigus_vulgaris"]
filtered_dotplot_color$pernicious_anemia = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "pernicious_anemia"]
filtered_dotplot_color$primary_biliary_cirrhosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "primary_biliary_cirrhosis"]
filtered_dotplot_color$primary_sclerosing_cholangitis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "primary_sclerosing_cholangitis"]
filtered_dotplot_color$ankylosing_spondylitis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "ankylosing_spondylitis"]
filtered_dotplot_color$psoriasis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "psoriasis"]
filtered_dotplot_color$sarcoidosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "sarcoidosis"]
filtered_dotplot_color$scleroderma_systemic_sclerosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "scleroderma_systemic_sclerosis"]
filtered_dotplot_color$sjogrens_syndrome = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "sjogrens_syndrome"]
filtered_dotplot_color$systemic_lupus_erythematosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "systemic_lupus_erythematosis"]
filtered_dotplot_color$temporal_arteritis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "temporal_arteritis"]
filtered_dotplot_color$vitiligo = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "vitiligo"]
filtered_dotplot_color$wegeners_granulomatosis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "wegeners_granulomatosis"]
filtered_dotplot_color$asthma = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "asthma"]
filtered_dotplot_color$autoimmune_hemolytic_anemia = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "autoimmune_hemolytic_anemia"]
filtered_dotplot_color$autoimmune_hepatitis = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "autoimmune_hepatitis"]


filtered_dotplot_color$antibiotics = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "antibiotics"]
filtered_dotplot_color$surgery = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "has_the_subject_had_a_prior_abdominal_surgery_other"]
filtered_dotplot_color$chemotherapy = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "chemotherapy"]
filtered_dotplot_color$immunosuppressants = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "immunosuppressants_e_g_oral_corticosteroids"]
filtered_dotplot_color$cholecystectomy = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "has_the_subject_had_a_cholecystectomy"]
filtered_dotplot_color$colonoscopy = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "x2_in_the_past_2_weeks_have_you_undergone_a_colonoscopy_or_other_procedure"]
filtered_dotplot_color$oral_contrast = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "x3_in_the_past_2_weeks_have_you_used_an_oral_contrast"]
filtered_dotplot_color$diarrhea = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "x4_in_the_past_2_weeks_have_you_had_diarrhea"]
filtered_dotplot_color$hospitalized = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "x5_in_the_past_2_weeks_have_you_been_hospitalized"]
filtered_dotplot_color$bowel_surgery = meta_subset[match(filtered_dotplot_color$sample, meta_subset$sample), "x6_have_you_ever_had_bowel_surgery"]

filtered_dotplot_color = unique(filtered_dotplot_color)

# Feature table: one row per sample, one column per species (MTX_RA)
feature_table_mtx <- filtered_dotplot_color %>%
        select(sample, S, MTX_RA) %>%
        filter(!is.na(MTX_RA)) %>%
        group_by(sample, S) %>%
        summarise(MTX_RA = mean(MTX_RA), .groups = "drop") %>%
        pivot_wider(names_from = S, values_from = MTX_RA) %>%
        column_to_rownames("sample")

# STEP 3: Match samples in both tables
metadata <- filtered_dotplot_color[, c(3, 11:84)] %>%
        distinct(sample, .keep_all = TRUE) %>%
        column_to_rownames("sample")
metadata = as.data.frame(metadata)
metadata = unique(metadata)

common_samples <- intersect(row.names(feature_table_mtx), row.names(metadata))
metadata <- metadata %>% filter(row.names(metadata) %in% common_samples)
feature_table_mtx <- feature_table_mtx[common_samples,]




metadata[] <- lapply(metadata, function(col) {
        if (is.character(col) || is.factor(col)) {
                col <- as.character(col)
                col <- gsub("[ ,]", "_", col)        # Replace commas and spaces with underscores
                col <- gsub("_+", "_", col)          # Replace multiple underscores with a single one
                col <- gsub("^_|_$", "", col)        # Remove leading/trailing underscores
                col <- ifelse(col == "", "unknown", col)  # Replace blanks with 'unknown'
                factor(col)
        } else {
                col  # leave numeric columns as-is
        }
})

#check on the categories and remove those with only one level
for (col in colnames(metadata)) {
        if (is.character(metadata[[col]]) || is.factor(metadata[[col]])) {
                cat("\n---", col, "---\n")
                print(unique(metadata[[col]]))
        }
}
metadata <- metadata[, sapply(metadata, function(col) length(unique(col)) > 1)]


# STEP 4: change the c-variates (except age) into factors
categorical_vars <- c(        "race", "sex", "general_wellbeing", 
        "water", "alcohol", "yogurt", "dairy", "probiotic", "fruits", "vegetables", "beans", "wheat", "starch", "eggs",
        "processed_meat", "red_meat", "shellfish", "fish", "sweets", "tea_or_coffee",
        "diagnosis", "antibiotics", "chemotherapy", "immunosuppressants", "colonoscopy", "oral_contrast", "diarrhea", "hospitalized", "bowel_surgery")
metadata[categorical_vars] <- lapply(metadata[categorical_vars], factor)


# Define the reference levels for categorical variables
fixed_effects = c("Age", "race","sex",  "general_wellbeing", 
        "water", "alcohol", "yogurt", "dairy", "probiotic", "fruits", "vegetables", "beans", "wheat", "starch", "eggs",
        "processed_meat", "red_meat", "shellfish", "fish", "sweets", "tea_or_coffee",
        "diagnosis", "antibiotics", "chemotherapy", "immunosuppressants", "colonoscopy", "oral_contrast", "diarrhea", "hospitalized", "bowel_surgery")

reference_levels <- c(
        "diagnosis,nonIBD",
        "race,White",
        "general_wellbeing,Very_Well",
        "alcohol,No_I_did_not_consume_these_products_in_the_last_7_days",
        "tea_or_coffee,No_I_did_not_consume_these_products_in_the_last_7_days",
        "yogurt,No_I_did_not_consume_these_products_in_the_last_7_days",
        "dairy,No_I_did_not_consume_these_products_in_the_last_7_days",
        "probiotic,No_I_did_not_consume_these_products_in_the_last_7_days",
        "fruits,No_I_did_not_consume_these_products_in_the_last_7_days",
        "vegetables,No_I_did_not_consume_these_products_in_the_last_7_days",
        "beans,No_I_did_not_consume_these_products_in_the_last_7_days",
        "wheat,No_I_did_not_consume_these_products_in_the_last_7_days",
        "starch,No_I_did_not_consume_these_products_in_the_last_7_days",
        "eggs,No_I_did_not_consume_these_products_in_the_last_7_days",
        "processed_meat,No_I_did_not_consume_these_products_in_the_last_7_days",
        "red_meat,No_I_did_not_consume_these_products_in_the_last_7_days",
        "shellfish,No_I_did_not_consume_these_products_in_the_last_7_days",
        "fish,No_I_did_not_consume_these_products_in_the_last_7_days",
        "sweets,No_I_did_not_consume_these_products_in_the_last_7_days",
        "water,No_I_did_not_consume_these_products_in_the_last_7_days"
)
# Run MaAsLin3
set.seed(1)

fit_data <- maaslin3(
        input_data = feature_table_mtx,
        input_metadata = metadata,
        output = "maaslin3_MTX_RA",
        fixed_effects = fixed_effects,
        random_effects = c("Participant"),
        normalization = "TSS",
        transform = "LOG",
        max_significance = 0.1,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        reference = reference_levels,
        standardize = TRUE
) 

##write into a R.script and run on server
feature_table_mtx$sample = row.names(feature_table_mtx)
metadata$sample = row.names(metadata)
write.table(feature_table_mtx, file = "/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/feature_table_mtx.tsv", sep = "\t", quote = FALSE, row.names =  FALSE)
write.table(metadata, file = "/Users/leawang/Documents/Lea/Harvard/MBTA_RNA/Cpn60/HMP2/metadata_maaslin3.tsv", sep = "\t", quote = FALSE, row.names =  FALSE)

# run_maaslin3.R

# 1. Load libraries
library(Maaslin3)
library(tidyverse)
library(data.table)

# 2. Load feature table (wide format: samples as rows, species as columns)
feature_table_mtx <- fread("feature_table_mtx.tsv") %>%
        column_to_rownames("sample") %>%
        as.data.frame()
# 3. Load metadata (samples as rows, variables as columns)
metadata <- fread("metadata_maaslin3.tsv") %>%
        column_to_rownames("sample") %>%
        as.data.frame()

# 4. Convert character columns to factors and clean up
categorical_vars <- c(
        "race", "sex", "general_wellbeing", 
        "water", "alcohol", "yogurt", "dairy", "probiotic", "fruits", "vegetables", 
        "beans", "wheat", "starch", "eggs", "processed_meat", "red_meat", "shellfish", 
        "fish", "sweets", "tea_or_coffee", "diagnosis", "antibiotics", "chemotherapy", 
        "immunosuppressants", "colonoscopy", "oral_contrast", "diarrhea", 
        "hospitalized", "bowel_surgery"
)
metadata[categorical_vars] <- lapply(metadata[categorical_vars], factor)

# 5. Reference levels (update if needed)
reference_levels <- c(
        "diagnosis,nonIBD",
        "race,White",
        "general_wellbeing,Very_Well",
        "alcohol,No_I_did_not_consume_these_products_in_the_last_7_days",
        "tea_or_coffee,No_I_did_not_consume_these_products_in_the_last_7_days",
        "yogurt,No_I_did_not_consume_these_products_in_the_last_7_days",
        "dairy,No_I_did_not_consume_these_products_in_the_last_7_days",
        "probiotic,No_I_did_not_consume_these_products_in_the_last_7_days",
        "fruits,No_I_did_not_consume_these_products_in_the_last_7_days",
        "vegetables,No_I_did_not_consume_these_products_in_the_last_7_days",
        "beans,No_I_did_not_consume_these_products_in_the_last_7_days",
        "wheat,No_I_did_not_consume_these_products_in_the_last_7_days",
        "starch,No_I_did_not_consume_these_products_in_the_last_7_days",
        "eggs,No_I_did_not_consume_these_products_in_the_last_7_days",
        "processed_meat,No_I_did_not_consume_these_products_in_the_last_7_days",
        "red_meat,No_I_did_not_consume_these_products_in_the_last_7_days",
        "shellfish,No_I_did_not_consume_these_products_in_the_last_7_days",
        "fish,No_I_did_not_consume_these_products_in_the_last_7_days",
        "sweets,No_I_did_not_consume_these_products_in_the_last_7_days",
        "water,No_I_did_not_consume_these_products_in_the_last_7_days",
        "colonoscopy,No",
        "oral_contrast,No",
        "diarrhea,No",
        "hospitalized,No",
        "bowel_surgery,No")

# 6. Fixed effects
fixed_effects <- c(
        "Age", "race", "sex", "general_wellbeing", 
        "water", "alcohol", "yogurt", "dairy", "probiotic", "fruits", "vegetables", 
        "beans", "wheat", "starch", "eggs", "processed_meat", "red_meat", "shellfish", 
        "fish", "sweets", "tea_or_coffee", "diagnosis", "antibiotics", "chemotherapy", 
        "immunosuppressants", "colonoscopy", "oral_contrast", "diarrhea", 
        "hospitalized", "bowel_surgery"
)


# 8. Run MaAsLin3
set.seed(1)
fit_data <- maaslin3(
        input_data = feature_table_mtx,
        input_metadata = metadata,
        output = "maaslin3_MTX_RA",
        fixed_effects = fixed_effects,
        random_effects = c("Participant"),
        normalization = "TSS",
        transform = "LOG",
        max_significance = 0.1,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        reference = reference_levels,
        standardize = TRUE
)


####maasline3 on MTX/MGX ratios




###Plot the most variable activity species (across participant over time)
filtered_dotplot_color <- filtered_dotplot_color %>% mutate(ratio = MTX_RA / MGX_RA)

subject_species_means <- filtered_dotplot_color %>% group_by(S, Participant) %>%
        summarise(mean_ratio = mean(ratio, na.rm = TRUE), .groups = "drop")

species_variability <- subject_species_means %>%group_by(S) %>%
        summarise(between_subject_var = var(mean_ratio, na.rm = TRUE), n_subjects = n()) %>%
        filter(n_subjects >= 5) %>%
        arrange(desc(between_subject_var))

top_species <- head(species_variability$S, 50)

filtered_top <- filtered_dotplot_color %>% filter(S %in% top_species)
filtered_top$date <- as.Date(filtered_top$date)
filtered_top$diagnosis <- factor(filtered_top$diagnosis, levels = c("nonIBD", "UC", "CD"))

#overall variation
filtered_top$diagnosis <- factor(filtered_top$diagnosis, levels = c("nonIBD", "UC", "CD"))
species_pvals <- filtered_top %>%
        group_by(S) %>%
        summarise(p_value = kruskal.test(ratio ~ diagnosis)$p.value) %>%
        mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
        filter(p_adj < 0.05)

filtered_sig <- filtered_top %>%filter(S %in% species_pvals$S)

ggplot(filtered_sig, aes(x = diagnosis, y = ratio, fill = diagnosis)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1, aes(color = diagnosis)) +
        facet_wrap(~ S, scales = "free_y") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.5) +
        scale_y_log10() +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        scale_color_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none",
                strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 0, hjust = 0.5)) +
        labs( title = "Significantly Different RNA:DNA Ratios by Diagnosis",
                x = "Diagnosis",
                y = "RNA:DNA Ratio (log10 scale)")

#change over time
ggplot(filtered_top, aes(x = date, y = ratio, group = Participant, color = diagnosis, fill = diagnosis)) +
        geom_ribbon(aes(ymin = 1, ymax = ratio), alpha = 0.15, color = NA) +
                geom_line(alpha = 0.8, size = 0.7) +
                geom_point(size = 1.2, alpha = 0.6) +
                facet_wrap(~ S, scales = "free_y") +
                scale_color_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
                scale_y_log10() +
                geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.5) +
                theme_bw(base_size = 12) +
        theme(legend.position = "top",
                strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1) ) +
                labs(title = "Top Species: RNA:DNA Ratio Trajectories Over Time by Diagnosis",
                x = "Sample Collection Date",
                y = "RNA:DNA Ratio (log10 scale)",
                color = "Diagnosis",
                fill = "Diagnosis")

#remove individual line
filtered_top <- filtered_dotplot_color %>% filter(S %in% top_species)
filtered_top$date <- as.Date(filtered_top$date)
filtered_top$diagnosis <- factor(filtered_top$diagnosis, levels = c("nonIBD", "UC", "CD"))

ggplot(filtered_top, aes(x = date, y = ratio, color = diagnosis, fill = diagnosis)) +
        geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
                geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.5) +
                facet_wrap(~ S, scales = "free_y") +
        scale_color_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        scale_fill_manual(values = c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")) +
        scale_y_log10() +
        theme_bw(base_size = 12) +
        theme(legend.position = "top",
                strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs( title = "Smoothed RNA:DNA Ratio Trends Over Time by Species and Diagnosis",
                x = "Sample Collection Date",
                y = "RNA:DNA Ratio (log10 scale)",
                color = "Diagnosis",
                fill = "Diagnosis" )



##peak plots
df_wide <- filtered_dotplot_color %>% select(sample, S, MGX_RA, MTX_RA, diagnosis)
df_wide <- df_wide %>% group_by(sample, S, diagnosis) %>%
        summarize(MGX_RA = mean(MGX_RA, na.rm = TRUE), MTX_RA = mean(MTX_RA, na.rm = TRUE),.groups = "drop" )

df_long <- df_wide %>%pivot_longer( cols= c(MGX_RA, MTX_RA),names_to  = "variable",values_to = "RA" )

my_cols <- c(MGX_RA = "steelblue", MTX_RA = "pink")

make_ridge <- function(diag_name) {
        sub_wide <- df_wide %>% filter(diagnosis == diag_name)
        sub_long <- df_long %>% filter(diagnosis == diag_name)
        
        ordering <- sub_wide %>%
                group_by(S) %>%
                summarize(mean_total = mean(MGX_RA + MTX_RA, na.rm = TRUE), .groups = "drop") %>%
                arrange(desc(mean_total)) %>%
                pull(S)
        
        sub_long <- sub_long %>%
                mutate(S = factor(S, levels = rev(ordering)))
                
        ggplot(sub_long, aes(x = RA, y = S, fill = variable)) +
                geom_density_ridges(alpha = 0.5, scale = 1, size = 0.2) +
                scale_fill_manual(values = my_cols) +
                coord_cartesian(xlim = c(0, 1)) +
                labs(title = paste0("Diagnosis: ", diag_name),
                        x     = "Relative abundance (RA)",
                        y     = "Species") +
                theme_ridges(center_axis_labels = TRUE) +scale_x_sqrt()+
                theme(legend.title = element_blank(),
                        axis.text.y  = element_text(size = 6))
}

p_nonIBD <- make_ridge("nonIBD")
p_CD     <- make_ridge("CD")
p_UC     <- make_ridge("UC")

print(p_nonIBD)
print(p_CD)
print(p_UC)

#mean and SD
stats <- df_wide %>% 
        pivot_longer(c(MGX_RA, MTX_RA), names_to="measure", values_to="RA") %>%
        group_by(diagnosis, S, measure) %>%
        summarize(mean_RA   = mean(RA, na.rm=TRUE),
                sd_RA     = sd(RA,   na.rm=TRUE),
                n_nonzero = sum(RA > 0, na.rm=TRUE),
                n_total   = n(),
                .groups = "drop" )

stats_clean <- stats %>% filter(n_nonzero > 0.1)

ggplot(stats_clean, aes(x = mean_RA, y = fct_reorder(S, mean_RA), color = measure)) +
        geom_point(position = position_dodge(width=0.6), size=2) +
        geom_errorbarh(aes(xmin = mean_RA - sd_RA, xmax = mean_RA + sd_RA),
                       position = position_dodge(width=0.6),height = 0) +
        facet_wrap(~ diagnosis, scales = "free_y") +
        scale_x_continuous(limits = c(0,1)) +
        labs(x = "Mean RA ± SD", y = "Species",color = "" ) +
        theme_light() + scale_x_sqrt()+
        theme(strip.text = element_text(face="bold"),axis.text.y = element_text(size=6))

#mean and SD in statistically significant species
# 1) Compute per-species mean absolute activity difference
mean_diff <- filtered_dotplot_color %>% mutate(deltaRA = MTX_RA - MGX_RA) %>% group_by(S) %>%
        summarize(mean_abs = mean(abs(deltaRA), na.rm=TRUE), .groups="drop")

top15 <- mean_diff %>% arrange(desc(mean_abs)) %>% slice(1:25) %>% pull(S) # 2) Top 15 by mean_abs

# 3) Combine with your sig_species vector
species_to_plot <- union(top15, sig_species)

# 4) Build summary stats for MGX_RA & MTX_RA
stats_plot <- filtered_dotplot_color %>%
        filter(S %in% species_to_plot) %>%
        group_by(diagnosis, S) %>%
        summarize(MGX_mean = mean(MGX_RA, na.rm=TRUE),
                MGX_sd   = sd(  MGX_RA, na.rm=TRUE),
                MTX_mean = mean(MTX_RA, na.rm=TRUE),
                MTX_sd   = sd(  MTX_RA, na.rm=TRUE),
                .groups  = "drop" ) %>%
        # reshape to long form with mean + sd columns
        pivot_longer( cols      = c(MGX_mean, MGX_sd, MTX_mean, MTX_sd),
                names_to  = c("measure","stat"),
                names_sep = "_",
                values_to = "value") %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        rename(mean_RA = mean, sd_RA = sd) %>%
        mutate(measure = recode(measure, MGX = "MGX_RA", MTX = "MTX_RA"),
                is_sig  = S %in% sig_species )

# 5) Determine the global ordering of species
species_order <- mean_diff %>%
        filter(S %in% species_to_plot) %>%
        arrange(desc(mean_abs)) %>% pull(S)

# 6) Turn S into a factor with that order
stats_plot <- stats_plot %>%
        mutate(S = factor(S, levels = rev(species_order)))

ggplot(stats_plot, aes( x     = mean_RA, y     = S, color = measure, alpha = is_sig)) +
        geom_point(position = position_dodge(width = 0.6), size = 2) +
        geom_errorbarh(aes(xmin = mean_RA - sd_RA,xmax = mean_RA + sd_RA ),
        height   = 0,
        position = position_dodge(width = 0.6)) +
        scale_color_manual( values = c(MGX_RA = "steelblue", MTX_RA = "pink") ) +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
        facet_wrap(~ diagnosis, scales = "free_y") +
        labs(x     = "Mean RA ± SD",
                y     = "Species",
                color = "",
                alpha = "Significant" ) +
        theme_light() +scale_x_sqrt()+
        theme(strip.text  = element_text(face = "bold"),
                axis.text.y = element_text(size = 6))

####show mean and sd RA plus spearman correlation rho, only sig species in spearman correlation
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(broom)
library(purrr)
library(tidytext)

# Step 1: Spearman correlation (same as before)
cor_res <- filtered_dotplot_color %>%
        group_by(diagnosis, S) %>%
        summarize( test = list(cor.test(MGX_RA, MTX_RA, method = "spearman", exact = FALSE)),  .groups = "drop" ) %>%
        mutate(rho     = map_dbl(test, ~ .x$estimate), p_value = map_dbl(test, ~ .x$p.value) ) %>%
        group_by(diagnosis) %>%
        mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
        ungroup() %>%
        filter(p_adj < 0.05)

# Step 2: Compute mean ± SD, and ordering by overall mean RA
stats_rho <- filtered_dotplot_color %>%
        semi_join(cor_res, by = c("diagnosis", "S")) %>%
        group_by(diagnosis, S) %>%
        summarize( MGX_mean = mean(MGX_RA, na.rm = TRUE),
                MGX_sd   = sd(MGX_RA,   na.rm = TRUE),
                MTX_mean = mean(MTX_RA, na.rm = TRUE),
                MTX_sd   = sd(MTX_RA,   na.rm = TRUE),
                total_RA = mean(MGX_RA + MTX_RA, na.rm = TRUE),
                .groups = "drop") %>%
        pivot_longer(cols = c(MGX_mean, MGX_sd, MTX_mean, MTX_sd),
                     names_to = c("Measure", "stat"),
                     names_sep = "_",
                     values_to = "value") %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        mutate(hape = ifelse(Measure == "MGX", 16, 17))

# Step 3: Merge with rho values and reorder species by total RA
stats_rho <- stats_rho %>%
        left_join(cor_res %>% select(diagnosis, S, rho), by = c("diagnosis", "S")) %>%
        mutate(S_ordered = reorder_within(S, total_RA, diagnosis))

# Step 4: Plot
ggplot(stats_rho, aes( x= mean, y = S_ordered,shape = Measure,color = rho)) +
        geom_point(position = position_dodge(width = 0.6), size = 2) +
        geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), height = 0, position = position_dodge(width = 0.6)) +
        scale_shape_manual(values = c(MGX = 16, MTX = 17)) +
        scale_color_gradient(low = "steelblue", high = "firebrick") +
        tidytext::scale_y_reordered() +
        scale_x_sqrt() +
        facet_wrap(~ diagnosis, scales = "free_y") +
        labs(x = "√(Mean RA ± SD)",y = "Species (sorted by RA)",color = "Spearman ρ",shape = "RA Type") +
        theme_light() +
        theme(strip.text  = element_text(face = "bold"), axis.text.y = element_text(size = 6) )

###show mean and sd RA plus activity,  only sig species in spearman correlation
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

eps <- 1e-6

# Step 1: Compute mean ± SD and log2(MTX/MGX) for all species × diagnosis
stats_activity <- filtered_dotplot_color %>%
        mutate(activity = log2((MTX_RA + eps) / (MGX_RA + eps))) %>%
        group_by(diagnosis, S) %>%
        summarize( MGX_mean = mean(MGX_RA, na.rm = TRUE),
                MGX_sd   = sd(MGX_RA,   na.rm = TRUE),
                MTX_mean = mean(MTX_RA, na.rm = TRUE),
                MTX_sd   = sd(MTX_RA,   na.rm = TRUE),
                activity = mean(activity, na.rm = TRUE),
                total_RA = mean(MGX_RA + MTX_RA, na.rm = TRUE),
                .groups  = "drop" )

# Step 2: Get reversed species order by nonIBD total RA (most abundant = top)
species_order <- stats_activity %>%filter(diagnosis == "nonIBD") %>% arrange(total_RA) %>%  pull(S)

# Step 3: Clean species labels for display
species_labels <- species_order %>%
        unique() %>%
        setNames(nm = ., object = gsub("_", ".", gsub("^s__", "", .)))

# Step 4: Prepare for plotting
plot_df <- stats_activity %>%
        pivot_longer(cols = c(MGX_mean, MGX_sd, MTX_mean, MTX_sd),
                     names_to = c("Measure", "stat"),
                     names_sep = "_",
                     values_to = "value") %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        mutate(Shape = ifelse(Measure == "MGX", 16, 17),
                S     = factor(S, levels = species_order))

ggplot(plot_df, aes( x     = mean, y     = S, shape = Measure,color = activity)) +
        geom_point(position = position_dodge(width = 0.6), size = 2) +
        geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                       height = 0, position = position_dodge(width = 0.6)) +
        scale_shape_manual(values = c(MGX = 16, MTX = 17)) +
        scale_color_gradient2( low = "steelblue", mid = "gray90",  high = "firebrick",  midpoint = 0, name = "log₂(MTX/MGX)") +
        scale_x_sqrt() +
        scale_y_discrete(labels = species_labels, name = "Species (sorted by RA in nonIBD)" ) +
        facet_wrap(~ diagnosis, scales = "free_y") +
        labs(x = "√(Mean RA ± SD)",shape = "RA Type" ) +
        theme_light(base_size = 10) +
        theme(strip.text  = element_text(face = "bold"),axis.text.y = element_text(face = "italic", size = 6) )



####antibiotics
No <- filtered_dotplot_color %>% filter(antibiotics == "No") %>% filter(MGX_RA >0 & MTX_RA > 0)
Yes <- filtered_dotplot_color %>% filter(antibiotics == "Yes")%>% filter(MGX_RA >0 & MTX_RA > 0)

No <- filter_common_species_all_groups(No)
Yes <- filter_common_species_all_groups(Yes)

cor_results_No <- compute_correlations(No, "No")
cor_results_Yes <- compute_correlations(Yes, "Yes")

col <- c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")
fill_col <- c("nonIBD" = "#45a8d9", "UC" = "#c49831", "CD" = "#963924")  # Same colors for shading
create_plot <- function(data, cor_results, site_label) {
        p <- ggplot(data, aes(x = MGX_RA, y = MTX_RA, color = health, fill = health)) + 
                scale_color_manual(values = col) +   scale_fill_manual(values = fill_col) +  # Fill colors for shading
                geom_point(size = 2, alpha = 0.5) +  
                theme_light() +
                geom_smooth(aes(group = health), method = lm, se = TRUE, fullrange = TRUE, alpha = 0.3) +  # Regression with shading
                labs(title = paste("Every Species in Every Sample -", site_label),
                     x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") +
                theme(legend.title = element_blank())
                for (i in 1:nrow(cor_results)) {
                p <- p + annotate("text", x = 0.1, y = max(data$MTX_RA) * (0.9 - 0.1 * i),
                                  label = cor_results$eqn[i], color = col[cor_results$health[i]], hjust = 0)
        }
        
        return(p)
}

p_No <- create_plot(No, cor_results_No, "No")
p_Yes <- create_plot(Yes, cor_results_Yes, "Yes")

library(dplyr)

#line plot and boxplot
filtered <- filtered_dotplot_color %>% mutate(date = as.Date(date), antibiotics_flag = ifelse(antibiotics == "Yes", "Yes", "No"))

first_abx <- filtered %>% filter(antibiotics_flag == "Yes") %>%  group_by(Participant) %>% summarise(first_abx_date = min(date))  # first time antibiotics used

filtered_abx_status <- filtered %>%  left_join(first_abx, by = "Participant") %>% mutate(abx_timing = case_when(
                is.na(first_abx_date) ~ "no_abx",                  # never had antibiotics
                date < first_abx_date ~ "pre",
                date >= first_abx_date ~ "post")) %>%
        filter(abx_timing %in% c("pre", "post"))  # exclude 'no_abx'

ratio_change <- filtered_abx_status %>%
        group_by(S, Participant, abx_timing) %>%
        summarise(mean_ratio = mean(ratio, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = abx_timing, values_from = mean_ratio) %>%
        filter(!is.na(pre) & !is.na(post)) %>%
        mutate(direction = case_when(
                post > pre ~ "increase",
                post < pre ~ "decrease",
                TRUE ~ "no_change"))

species_trends <- ratio_change %>%
        group_by(S) %>%
        summarise( n = n(),
                n_increase = sum(direction == "increase"),
                n_decrease = sum(direction == "decrease"),
                n_consistent = max(n_increase, n_decrease),
                percent_consistent = n_consistent / n
        ) %>%
        filter(n >= 3, percent_consistent >= 0.7) %>%  # at least 3 participants and ≥70% consistent change
        arrange(desc(percent_consistent))

print(species_trends)

top_species_consistent <- species_trends$S
paired_data <- ratio_change %>% filter(S %in% top_species_consistent)
paired_long <- paired_data %>%
        pivot_longer(cols = c(pre, post), names_to = "timing", values_to = "ratio") %>%
        mutate(timing = factor(timing, levels = c("pre", "post")))

ggplot(paired_long, aes(x = timing, y = ratio, group = Participant)) +
        geom_boxplot(aes(group = timing), fill = "gray90", color = "black", outlier.shape = NA) +
        geom_line(aes(color = timing), alpha = 0.5, linewidth = 0.8) +
        geom_point(aes(color = timing), size = 1.5) +
        facet_wrap(~ S, scales = "free_y") +
        scale_y_log10() +
        scale_color_manual(values = c("pre" = "steelblue", "post" = "firebrick")) +
        theme_bw(base_size = 12) +
        theme( strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 0, hjust = 0.5),
                legend.position = "none" ) +
        labs(title = "Per-Individual RNA:DNA Ratio Before vs. After Antibiotics",
                x = "Antibiotic Timing",
                y = "RNA:DNA Ratio (log10 scale)" )

# oberall line plot
filtered_top <- filtered_top %>%
        mutate( antibiotics_flag = ifelse(antibiotics == "Yes", "Yes", "No"),
                antibiotics_flag = factor(antibiotics_flag, levels = c("No", "Yes")))

ggplot(filtered_top, aes(x = date, y = ratio, color = antibiotics_flag, fill = antibiotics_flag)) +
        geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
                geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.5) +
                facet_wrap(~ S, scales = "free_y") +
                scale_color_manual(values = c("No" = "steelblue", "Yes" = "firebrick")) +
        scale_fill_manual(values = c("No" = "steelblue", "Yes" = "firebrick")) +
                scale_y_log10() +
        theme_bw(base_size = 12) +
        theme(legend.position = "top",
                strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1) ) +
        labs( title = "Smoothed RNA:DNA Ratio Trends Over Time by Antibiotic Use",
                x = "Sample Collection Date",
                y = "RNA:DNA Ratio (log10 scale)",
                color = "Antibiotics",
                fill = "Antibiotics")


#blood
Little <- filtered_dotplot_color %>% filter(blood == "None" | blood == "Trace (a little blood)") %>% filter(MGX_RA >0 & MTX_RA > 0)
Frank <- filtered_dotplot_color %>% filter(blood == "Occasionally frank (ocassionally a lot of blood)" | blood == "Usually frank (usually a lot of blood)")%>% filter(MGX_RA >0 & MTX_RA > 0)

filter_common_species_all_groups <- function(df) {
        df <- as.data.frame(df)  # convert from data.table if needed
        species_to_keep <- df %>%
                group_by(S, health) %>%
                summarize(n_samples = n_distinct(sample), .groups = "drop") %>%
                filter(n_samples >= 2) %>%
                group_by(S) %>%
                summarize(n_groups = n_distinct(health), .groups = "drop") %>%
                filter(n_groups == 2) %>%   # must appear in 2 groups, since Little health control here
                pull(S)
        df %>% filter(S %in% species_to_keep)
}

Little <- filter_common_species_all_groups(Little)
Frank <- filter_common_species_all_groups(Frank)

cor_results_Little <- compute_correlations(Little, "Little")
cor_results_Frank <- compute_correlations(Frank, "Frank")

create_plot <- function(data, cor_results, site_label) {
        p <- ggplot(data, aes(x = MGX_RA, y = MTX_RA, color = health, fill = health)) + 
                scale_color_manual(values = col) +   scale_fill_manual(values = fill_col) +  # Fill colors for shading
                geom_point(size = 2, alpha = 0.5) +  
                scale_y_sqrt() + scale_x_sqrt() +  # Square root scaling for better visualization
                theme_light() +
                geom_smooth(aes(group = health), method = lm, se = TRUE, fullrange = TRUE, alpha = 0.3) +  # Regression with shading
                labs(title = paste("Every Species in Every Sample -", site_label),
                     x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") +
                theme(legend.title = element_blank())
        for (i in 1:nrow(cor_results)) {
                p <- p + annotate("text", x = 0.1, y = max(data$MTX_RA) * (0.9 - 0.1 * i),
                                  label = cor_results$eqn[i], color = col[cor_results$health[i]], hjust = 0)
        }
        return(p)
}

p_Little <- create_plot(Little, cor_results_Little, "Little")
p_Frank <- create_plot(Frank, cor_results_Frank, "Frank")

#blood_yes or no
no <- filtered_dotplot_color %>% filter(blood == "None") %>% filter(MGX_RA >0 & MTX_RA > 0)
yes <- filtered_dotplot_color %>% filter( blood == "Trace (a little blood)" | blood == "Occasionally yes (ocassionally a lot of blood)" | blood == "Usually yes (usually a lot of blood)")%>% filter(MGX_RA >0 & MTX_RA > 0)
cor_results_no <- compute_correlations(no, "no")
cor_results_yes <- compute_correlations(yes, "yes")
p_no <- create_plot(no, cor_results_no, "no")
p_yes <- create_plot(yes, cor_results_yes, "yes")

#Immunosuppressants
no <- filtered_dotplot_color %>% filter(`Immunosuppressants (e.g. oral corticosteroids)` == "No") %>% filter(MGX_RA >0 & MTX_RA > 0)
yes <- filtered_dotplot_color %>% filter( `Immunosuppressants (e.g. oral corticosteroids)`== "Yes")%>% filter(MGX_RA >0 & MTX_RA > 0)
cor_results_no <- compute_correlations(no, "no")
cor_results_yes <- compute_correlations(yes, "yes")
p_no <- create_plot(no, cor_results_no, "no")
p_yes <- create_plot(yes, cor_results_yes, "yes")

#Diarrhea
no <- filtered_dotplot_color %>% filter(`4) In the past 2 weeks, have you had diarrhea?` == "No") %>% filter(MGX_RA >0 & MTX_RA > 0)
yes <- filtered_dotplot_color %>% filter( `4) In the past 2 weeks, have you had diarrhea?`== "Yes")%>% filter(MGX_RA >0 & MTX_RA > 0)
cor_results_no <- compute_correlations(no, "no")
cor_results_yes <- compute_correlations(yes, "yes")
p_no <- create_plot(no, cor_results_no, "no")
p_yes <- create_plot(yes, cor_results_yes, "yes")

#Arthralgias
No <- filtered_dotplot_color %>% filter(Arthralgias.x == "No") %>% filter(MGX_RA >0 & MTX_RA > 0)
Yes <- filtered_dotplot_color %>% filter(Arthralgias.x == "Yes")%>% filter(MGX_RA >0 & MTX_RA > 0)

No <- filter_common_species_all_groups(No)
Yes <- filter_common_species_all_groups(Yes) ##only got UC and CD, no healthy controls

cor_results_No <- compute_correlations(No, "No")
cor_results_Yes <- compute_correlations(Yes, "Yes")
create_plot <- function(data, cor_results, site_label) {
        p <- ggplot(data, aes(x = MGX_RA, y = MTX_RA, color = health, fill = health)) + 
                scale_color_manual(values = col) +   scale_fill_manual(values = fill_col) +  # Fill colors for shading
                geom_point(size = 2, alpha = 0.5) +  
                scale_y_sqrt() + scale_x_sqrt() +  # Square root scaling for better visualization
                theme_light() +
                geom_smooth(aes(group = health), method = lm, se = TRUE, fullrange = TRUE, alpha = 0.3) +  # Regression with shading
                labs(title = paste("Every Species in Every Sample -", site_label),
                     x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") +
                theme(legend.title = element_blank())
        for (i in 1:nrow(cor_results)) {
                p <- p + annotate("text", x = 0.1, y = max(data$MTX_RA) * (0.9 - 0.1 * i),
                                  label = cor_results$eqn[i], color = col[cor_results$health[i]], hjust = 0)
        }
        return(p)
}
p_No <- create_plot(No, cor_results_No, "No")
p_Yes <- create_plot(Yes, cor_results_Yes, "Yes")

#soft drinks
No <- filtered_dotplot_color %>% filter(`Soft drinks, tea or coffee with sugar (corn syrup, maple syrup, cane sugar, etc)` == "No, I did not consume these products in the last 7 days") %>% filter(MGX_RA >0 & MTX_RA > 0)
Yes <- filtered_dotplot_color %>% filter(!`Soft drinks, tea or coffee with sugar (corn syrup, maple syrup, cane sugar, etc)` == "No, I did not consume these products in the last 7 days" | `Soft drinks, tea or coffee with sugar (corn syrup, maple syrup, cane sugar, etc)` == "")%>% filter(MGX_RA >0 & MTX_RA > 0)
filter_common_species_all_groups <- function(df) {
        df <- as.data.frame(df)  # convert from data.table if needed
        species_to_keep <- df %>%
                group_by(S, health) %>%
                summarize(n_samples = n_distinct(sample), .groups = "drop") %>%
                filter(n_samples >= 2) %>%
                group_by(S) %>%
                summarize(n_groups = n_distinct(health), .groups = "drop") %>%
                filter(n_groups == 2) %>%   # must appear in 2 groups, since Little health control here
                pull(S)
        df %>% filter(S %in% species_to_keep)
}

No <- filter_common_species_all_groups(No)
Yes <- filter_common_species_all_groups(Yes) ##only got UC and CD, no healthy controls

cor_results_No <- compute_correlations(No, "No")
cor_results_Yes <- compute_correlations(Yes, "Yes")
create_plot <- function(data, cor_results, site_label) {
        p <- ggplot(data, aes(x = MGX_RA, y = MTX_RA, color = health, fill = health)) + 
                scale_color_manual(values = col) +   scale_fill_manual(values = fill_col) +  # Fill colors for shading
                geom_point(size = 2, alpha = 0.5) +  
                scale_y_sqrt() + scale_x_sqrt() +  # Square root scaling for better visualization
                theme_light() +
                geom_smooth(aes(group = health), method = lm, se = TRUE, fullrange = TRUE, alpha = 0.3) +  # Regression with shading
                labs(title = paste("Every Species in Every Sample -", site_label),
                     x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") +
                theme(legend.title = element_blank())
        for (i in 1:nrow(cor_results)) {
                p <- p + annotate("text", x = 0.1, y = max(data$MTX_RA) * (0.9 - 0.1 * i),
                                  label = cor_results$eqn[i], color = col[cor_results$health[i]], hjust = 0)
        }
        return(p)
}

p_No <- create_plot(No, cor_results_No, "No")
p_Yes <- create_plot(Yes, cor_results_Yes, "Yes")

#diet drinks
No <- filtered_dotplot_color %>% filter(`Diet soft drinks, tea or coffee with sugar (Stevia, Equal, Splenda etc)` == "No, I did not consume these products in the last 7 days") %>% filter(MGX_RA >0 & MTX_RA > 0)
Yes <- filtered_dotplot_color %>% filter(!`Diet soft drinks, tea or coffee with sugar (Stevia, Equal, Splenda etc)` == "No, I did not consume these products in the last 7 days" | `Diet soft drinks, tea or coffee with sugar (Stevia, Equal, Splenda etc)` == "")%>% filter(MGX_RA >0 & MTX_RA > 0)
filter_common_species_all_groups <- function(df) {
        df <- as.data.frame(df)  # convert from data.table if needed
        species_to_keep <- df %>%
                group_by(S, health) %>%
                summarize(n_samples = n_distinct(sample), .groups = "drop") %>%
                filter(n_samples >= 2) %>%
                group_by(S) %>%
                summarize(n_groups = n_distinct(health), .groups = "drop") %>%
                filter(n_groups == 2) %>%   # must appear in 2 groups, since Little health control here
                pull(S)
        df %>% filter(S %in% species_to_keep)
}

No <- filter_common_species_all_groups(No)
Yes <- filter_common_species_all_groups(Yes) ##only got UC and CD, no healthy controls

cor_results_No <- compute_correlations(No, "No")
cor_results_Yes <- compute_correlations(Yes, "Yes")

p_No <- create_plot(No, cor_results_No, "No")
p_Yes <- create_plot(Yes, cor_results_Yes, "Yes")

#Alcohol
no <- filtered_dotplot_color %>% filter(`Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)` == "No, I did not consume these products in the last 7 days" ) %>% filter(MGX_RA >0 & MTX_RA > 0)
yes <- filtered_dotplot_color %>% filter(`Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)` == "Within the past 4 to 7 days" |`Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)` == "Within the past 2 to 3 days" | `Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)` == "Yesterday, 1 to 2 times" | `Alcohol (beer, brandy, spirits, hard liquor, wine, aperitif, etc.)` =="Yesterday, 3 or more times")%>% filter(MGX_RA >0 & MTX_RA > 0)

no <- filter_common_species_all_groups(no)
yes <- filter_common_species_all_groups(yes)

cor_results_no <- compute_correlations(no, "no")
cor_results_yes <- compute_correlations(yes, "yes")

p_no <- create_plot(no, cor_results_no, "no")
p_yes <- create_plot(yes, cor_results_yes, "yes")

#yogurt
no <- filtered_dotplot_color %>% filter(`Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)` == "No, I did not consume these products in the last 7 days" ) %>% filter(MGX_RA >0 & MTX_RA > 0)
yes <- filtered_dotplot_color %>% filter(`Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)` == "Within the past 4 to 7 days" |`Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)` == "Within the past 2 to 3 days" | `Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)` == "Yesterday, 1 to 2 times" | `Yogurt or other foods containing active bacterial cultures (kefir, sauerkraut)` =="Yesterday, 3 or more times")%>% filter(MGX_RA >0 & MTX_RA > 0)

no <- filter_common_species_all_groups(no)
yes <- filter_common_species_all_groups(yes)

cor_results_no <- compute_correlations(no, "no")
cor_results_yes <- compute_correlations(yes, "yes")

p_no <- create_plot(no, cor_results_no, "no")
p_yes <- create_plot(yes, cor_results_yes, "yes")

#chemotherapy
No <- filtered_dotplot_color %>% filter(Chemotherapy== "No") %>% filter(MGX_RA >0 & MTX_RA > 0)
Yes <- filtered_dotplot_color %>% filter(Chemotherapy == "Yes")%>% filter(MGX_RA >0 & MTX_RA > 0)
filter_common_species_all_groups <- function(df) {
        df <- as.data.frame(df)  # convert from data.table if needed
        species_to_keep <- df %>%
                group_by(S, health) %>%
                summarize(n_samples = n_distinct(sample), .groups = "drop") %>%
                filter(n_samples >= 3) %>%
                group_by(S) %>%
                summarize(n_groups = n_distinct(health), .groups = "drop") %>%
                filter(n_groups == 2) %>%   # must appear in 2 groups, since Little health control here
                pull(S)
        df %>% filter(S %in% species_to_keep)
}

No <- filter_common_species_all_groups(No)
Yes <- filter_common_species_all_groups(Yes) ##only got UC and CD, no healthy controls

cor_results_No <- compute_correlations(No, "No")
cor_results_Yes <- compute_correlations(Yes, "Yes")

p_No <- create_plot(No, cor_results_No, "No")
p_Yes <- create_plot(Yes, cor_results_Yes, "Yes")

#chemotherapy--by disease
nonIBD <- filtered_dotplot_color %>% filter(health == "nonIBD") %>% filter(MGX_RA >0 & MTX_RA > 0)
UC <- filtered_dotplot_color %>% filter(health == "UC")%>% filter(MGX_RA >0 & MTX_RA > 0)
CD <- filtered_dotplot_color %>% filter(health == "CD")%>% filter(MGX_RA >0 & MTX_RA > 0)

compute_correlations <- function(data, species_name) {
        cor_results <- data.frame(
                species = character(),
                Chemotherapy = character(),
                cor_val = numeric(),
                p_val = numeric(),
                p_adj = numeric(),
                slope = numeric(),
                intercept = numeric(),
                eqn = character(),
                stringsAsFactors = FALSE
        )
        for (group in unique(data$Chemotherapy)) {
                subset_data <- data %>% filter(Chemotherapy == group)
                if (nrow(subset_data) > 2) {  # Ensure enough data points for correlation
                        cor_test <- cor.test(subset_data$MGX_RA, subset_data$MTX_RA, method = "spearman")
                        cor_val <- round(cor_test$estimate, 2)
                        p_val <- round(cor_test$p.value, 4)
                        lm_fit <- lm(MTX_RA ~ MGX_RA, data = subset_data)
                        slope <- round(coef(lm_fit)[2], 2)
                        intercept <- round(coef(lm_fit)[1], 2)
                        p_adj <- p.adjust(p_val, method = "fdr")
                        eqn <- paste("y =", slope, "*x +", intercept, "; R² =", round(cor_val^2, 3), ", p(FDR) =", signif(p_adj, 3))
                        cor_results <- rbind(cor_results, data.frame(
                                species = species_name,
                                Chemotherapy = group,
                                cor_val = cor_val,
                                p_val = p_val,
                                p_adj = p_adj,
                                slope = slope,
                                intercept = intercept,
                                eqn = eqn,
                                stringsAsFactors = FALSE
                        ))
                }
        }
        return(cor_results)
}

cor_results_nonIBD <- compute_correlations(nonIBD, "nonIBD")
cor_results_CD <- compute_correlations(CD, "CD")
cor_results_UC <- compute_correlations(UC, "UC")

col <- c("No" = "steelblue", "Yes" = "darkorange")
fill_col <- c("No" = "steelblue", "Yes" = "darkorange")  # Same colors for shading
create_plot <- function(data, cor_results, site_label) {
        p <- ggplot(data, aes(x = MGX_RA, y = MTX_RA, color = Chemotherapy, fill = Chemotherapy)) + 
                scale_color_manual(values = col) +   scale_fill_manual(values = fill_col) +  # Fill colors for shading
                geom_point(size = 2, alpha = 0.3) +  
                scale_y_sqrt() + scale_x_sqrt() +  # Square root scaling for better visualization
                theme_light() +
                geom_smooth(aes(group = Chemotherapy), method = lm, se = TRUE, fullrange = TRUE, alpha = 0.3) +  # Regression with shading
                labs(title = paste("Every Species in Every Sample -", site_label),
                     x = "MGX cpn60 taxa RA", y = "MTX cpn60 taxa RA") +
                theme(legend.title = element_blank())
        for (i in 1:nrow(cor_results)) {
                p <- p + annotate("text", x = 0.1, y = max(data$MTX_RA) * (0.9 - 0.1 * i),
                                  label = cor_results$eqn[i], color = col[cor_results$health[i]], hjust = 0)
        }
        return(p)
}

p_nonIBD <- create_plot(nonIBD, cor_results_nonIBD, "nonIBD")
p_CD <- create_plot(CD, cor_results_CD, "CD")
p_UC <- create_plot(UC, cor_results_CD, "UC")


####check all_gene MTX vs. cpn_gene MTX in oral bugs
irep = fread("/Users/leawang/Dropbox (Harvard University)/hutlab/Lea/Manuscript/Cpn60_manuscript/data/Supple_5_iRep.csv")
irep$sample_S = paste(irep$sample_1, irep$S, sep = "-")
irep$major_site = filtered_dotplot_color[match(irep$sample_S, filtered_dotplot_color$V1),"major_site"]
irep$health = filtered_dotplot_color[match(irep$sample_1, filtered_dotplot_color$sample), "health"]

# Perform Spearman correlation test
cor_test <- cor.test(irep$MTX_cpn_RA, irep$MTX_humann_RA, method = "spearman")
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)

lm_fit <- lm(MTX_humann_RA ~ MTX_cpn_RA, data = irep)
slope <- round(coef(lm_fit)[2], 2)
intercept <- round(coef(lm_fit)[1], 2)
eqn <- paste("y =", slope, "*x +", intercept)

x_max <- max(irep$MTX_cpn_RA, na.rm = TRUE)
y_max <- max(irep$MTX_humann_RA, na.rm = TRUE)
y_min <- min(irep$MTX_humann_RA, na.rm = TRUE)

p <- ggplot(irep, aes(x = MTX_cpn_RA, y = MTX_humann_RA)) +
        geom_point(color = "dark green", alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, fill = "lightgreen", alpha = 0.5) +
        annotate("text", x = x_max * 0.9, y = y_max * 0.9, 
                 label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1, size = 5) +
        annotate("text", x = x_max * 0.9, y = y_min * 1.1, 
                 label = eqn, hjust = 1, vjust = 0, size = 5) +
        labs(x = "MTX cpn RA", y = "MTX humann RA", title = "Comparison of MTX_cpn_RA and MTX_humann_RA") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_sqrt() + 
        scale_x_sqrt()
ggsave("check_1.pdf", plot = p, width = 5, height = 5)

#MTX/MGX ratio in cpn vs. humann
cor_test <- cor.test(irep$ratio, irep$humann_MTX_MGX_ratio, method = "spearman")
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(humann_MTX_MGX_ratio ~ ratio, data = irep)
slope <- round(coef(lm_fit)[2], 2)
intercept <- round(coef(lm_fit)[1], 2)
eqn <- paste("y =", slope, "*x +", intercept)

x_max <- max(irep$ratio, na.rm = TRUE)
y_max <- max(irep$humann_MTX_MGX_ratio, na.rm = TRUE)
y_min <- min(irep$humann_MTX_MGX_ratio, na.rm = TRUE)
irep_1 = irep[!irep$humann_MTX_MGX_ratio == "NA",]
irep_1 = irep_1[!irep_1$humann_MTX_MGX_ratio == 0,]
irep_1 = irep_1[!irep_1$ratio == 0,]

p <- ggplot(irep_1, aes(x = ratio, y = humann_MTX_MGX_ratio)) +
        geom_point(color = "dark green", alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, fill = "lightgreen", alpha = 0.5) +
        annotate("text", x = x_max * 0.9, y = y_max * 0.9, 
                 label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1, size = 5) +
        annotate("text", x = x_max * 0.9, y = y_min * 1.1, 
                 label = eqn, hjust = 1, vjust = 0, size = 5) +
        labs(x = "MTX/MGX ratio (cpn60)", y = "MTX/MGX ratio (humann)", title = "Comparison of MTX/MGX ratio") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_log10() + 
        scale_x_log10()
p
#check in oral
irep_oral = irep[irep$major_site == "oral",]
cor_test <- cor.test(irep_oral$MTX_cpn_RA, irep_oral$MTX_humann_RA, method = "spearman")
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)

lm_fit <- lm(MTX_humann_RA ~ MTX_cpn_RA, data = irep_oral)
slope <- round(coef(lm_fit)[2], 2)
intercept <- round(coef(lm_fit)[1], 2)
eqn <- paste("y =", slope, "*x +", intercept)

x_max <- max(irep_oral$MTX_cpn_RA, na.rm = TRUE)
y_max <- max(irep_oral$MTX_humann_RA, na.rm = TRUE)
y_min <- min(irep_oral$MTX_humann_RA, na.rm = TRUE)

p <- ggplot(irep_oral, aes(x = MTX_cpn_RA, y = MTX_humann_RA)) +
        geom_point(color = "dark green", alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, fill = "lightgreen", alpha = 0.5) +
        annotate("text", x = x_max * 0.9, y = y_max * 0.9, 
                 label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1, size = 5) +
        annotate("text", x = x_max * 0.9, y = y_min * 1.1, 
                 label = eqn, hjust = 1, vjust = 0, size = 5) +
        labs(x = "MTX cpn RA", y = "MTX humann RA", title = "Comparison of MTX_cpn_RA and MTX_humann_RA") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_sqrt() + 
        scale_x_sqrt()
ggsave("check_oral.pdf", plot = p, width = 5, height = 5)

###randoly select 300 HMP samples, nonIBD, UC and CD each 100
meta = read.csv("hmp2_metadata.csv", sep = ",", header = T)
mtx = meta[meta$data_type == "metatranscriptomics",]
mgx = meta[meta$data_type == "metagenomics",]
mtx_sample = unique(mtx$External.ID)
mgx_sample = unique(mgx$External.ID)
sample = mtx_sample[mtx_sample %in% mgx_sample]

tax_file = read.csv("/n/holystore01/LABS/huttenhower_lab/Users/leawang0705/cpn60/cpn_gene_23545_utaxID.tsv", sep = ",", header = T)
uniref90 = read.csv("uniref90_total_23230.csv", sep = ",", header = T)
who = c(uniref90$who_total_23230)

set.seed(123)  
selected_samples <- meta %>%
        filter(External.ID %in% sample) %>%  # Keep only samples in `sample`
        filter(diagnosis %in% c("nonIBD", "CD", "UC")) %>%  # Ensure correct diagnoses
        group_by(diagnosis) %>%
        slice_sample(n = 50) %>%
        ungroup() %>%
        select(External.ID) %>%
        distinct()
selected_samples$External.ID <- paste(selected_samples$External.ID, "Abundance-RPKs", sep = "_")
selected_samples <- unique(selected_samples$External.ID)

library(dplyr)
selected_samples_clean <- gsub("_Abundance-RPKs", "", selected_samples) ###check the final number of samples subset
meta$External.ID <- as.character(meta$External.ID)
meta$diagnosis <- as.character(meta$diagnosis)  
selected_meta <- meta %>%
        filter(External.ID %in% selected_samples_clean) %>%
        distinct(External.ID, diagnosis) %>%  
        group_by(diagnosis) %>%
        summarise(count = n(), .groups = "drop")  
print(selected_meta)

diagnosis count
<chr>     <int>
        1 CD           46
2 UC           45
3 nonIBD       46

MGX_300 = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", selected_samples))
MGX_300 = arrange(MGX_300, `# Gene Family`)
MGX_300 = MGX_300[grepl("\\|", MGX_300$`# Gene Family`),] #selecting only the rows where the # Gene Family column contains the "|" character

MGX_300$uniref90 = MGX_300$`# Gene Family`
MGX_300$uniref90 = gsub("\\|.*", "",MGX_300$uniref90)
MGX_300$taxa = MGX_300$`# Gene Family`
MGX_300$uniref90 = NULL
MGX_300 = MGX_300 %>% separate(taxa, sep = "\\|", c("uniref90", "taxa"))

MGX_300 <- MGX_300 %>%select(-c(`# Gene Family`, uniref90))
MGX_300 = MGX_300 %>% group_by(taxa) %>% summarise_all(mean)
MGX_300 = as.data.frame(MGX_300)
MGX_300 = melt(MGX_300) 

MGX_300 <- MGX_300 %>%
        group_by(variable) %>%
        mutate(MGX_RA_humann = value / sum(value, na.rm = TRUE)) %>% #relative abundance of each species
        ungroup()

check_sums <- MGX_300 %>% group_by(variable) %>% summarise(sum_RA = sum(MGX_RA_humann, na.rm = TRUE))
check_sums

MTX_300 = fread("MTX/genefamilies.tsv", select = c("# Gene Family", selected_samples))
MTX_300 = arrange(MTX_300, `# Gene Family`)
MTX_300 = MTX_300[grepl("\\|", MTX_300$`# Gene Family`),]
MTX_300$uniref90 = MTX_300$`# Gene Family`
MTX_300$uniref90 = gsub("\\|.*", "",MTX_300$uniref90)
MTX_300$taxa = MTX_300$`# Gene Family`
MTX_300$uniref90 = NULL
MTX_300 = MTX_300 %>% separate(taxa, sep = "\\|", c("uniref90", "taxa"))

MTX_300 <- MTX_300 %>%select(-c(`# Gene Family`, uniref90))
MTX_300 = MTX_300 %>% group_by(taxa) %>% summarise_all(mean)
MTX_300 = as.data.frame(MTX_300)
MTX_300 = melt(MTX_300) 

MTX_300 <- MTX_300 %>%
        group_by(variable) %>%
        mutate(MTX_RA_humann = value / sum(value, na.rm = TRUE)) %>% #relative abundance of each species
        ungroup()

check_sums <- MTX_300 %>% group_by(variable) %>% summarise(sum_RA = sum(MTX_RA, na.rm = TRUE))
check_sums

MGX_300$sample_S = paste(MGX_300$variable, MGX_300$taxa, sep = "-")
MTX_300$sample_S = paste(MTX_300$variable, MTX_300$taxa, sep = "-")
MGX_300$MTX_RA_humann = MTX_300[match(MGX_300$sample_S, MTX_300$sample_S), "MTX_RA_humann"]
MGX_300$sample_S <- gsub("_Abundance-RPKs-g__[^.]+\\.", "-", MGX_300$sample_S)
MGX_300$variable = gsub("_Abundance-RPKs", "", MGX_300$variable)
MGX_300$taxa <- gsub(".*s__", "s__", MGX_300$taxa)

#do the same for cpn60 UniRep90s
MGX_300_cpn = fread("MGX/HMP2_humann3.6_genefamilies.tsv", select = c("# Gene Family", selected_samples))
MGX_300_cpn = MGX_300_cpn[grepl("\\|", MGX_300_cpn$`# Gene Family`),] #selecting only the rows where the # Gene Family column contains the "|" character
MGX_300_cpn = arrange(MGX_300_cpn, `# Gene Family`)
MGX_300_cpn$uniref90 = MGX_300_cpn$`# Gene Family`
MGX_300_cpn$uniref90 = gsub("\\|.*", "",MGX_300_cpn$uniref90)
MGX_300_cpn = MGX_300_cpn[MGX_300_cpn$uniref90 %in% who,]

MGX_300_cpn$taxa = MGX_300_cpn$`# Gene Family`
MGX_300_cpn$uniref90 = NULL
MGX_300_cpn = MGX_300_cpn %>% separate(taxa, sep = "\\|", c("uniref90", "taxa"))
MGX_300_cpn <- MGX_300_cpn %>%select(-c(`# Gene Family`, uniref90))
MGX_300_cpn = MGX_300_cpn %>% group_by(taxa) %>% summarise_all(sum)
MGX_300_cpn = as.data.frame(MGX_300_cpn)
MGX_300_cpn = melt(MGX_300_cpn) 

MGX_300_cpn <- MGX_300_cpn %>%
        group_by(variable) %>%
        mutate(MGX_RA_cpn = value / sum(value, na.rm = TRUE)) %>% #relative abundance of each species
        ungroup()

check_sums <- MGX_300_cpn %>% group_by(variable) %>% summarise(sum_RA = sum(MGX_RA_cpn, na.rm = TRUE))
check_sums

MTX_300_cpn = fread("MTX/genefamilies.tsv", select = c("# Gene Family", selected_samples))
MTX_300_cpn = MTX_300_cpn[grepl("\\|", MTX_300_cpn$`# Gene Family`),]
MTX_300_cpn = arrange(MTX_300_cpn, `# Gene Family`)
MTX_300_cpn$uniref90 = MTX_300_cpn$`# Gene Family`
MTX_300_cpn$uniref90 = gsub("\\|.*", "",MTX_300_cpn$uniref90)
MTX_300_cpn = MTX_300_cpn[MTX_300_cpn$uniref90 %in% who,]
MTX_300_cpn$taxa = MTX_300_cpn$`# Gene Family`
MTX_300_cpn$uniref90 = NULL
MTX_300_cpn = MTX_300_cpn %>% separate(taxa, sep = "\\|", c("uniref90", "taxa"))

MTX_300_cpn <- MTX_300_cpn %>%select(-c(`# Gene Family`, uniref90))
MTX_300_cpn = MTX_300_cpn %>% group_by(taxa) %>% summarise_all(sum)
MTX_300_cpn = as.data.frame(MTX_300_cpn)
MTX_300_cpn = melt(MTX_300_cpn) 

MTX_300_cpn <- MTX_300_cpn %>%
        group_by(variable) %>%
        mutate(MTX_RA_cpn = value / sum(value, na.rm = TRUE)) %>% #relative abundance of each species
        ungroup()

check_sums <- MTX_300_cpn %>% group_by(variable) %>% summarise(sum_RA = sum(MTX_RA_cpn, na.rm = TRUE))
check_sums

MGX_300_cpn$sample_S = paste(MGX_300_cpn$variable, MGX_300_cpn$taxa, sep = "-")
MTX_300_cpn$sample_S = paste(MTX_300_cpn$variable, MTX_300_cpn$taxa, sep = "-")
MGX_300_cpn$sample_S <- gsub("_Abundance-RPKs-g__[^.]+\\.", "-", MGX_300_cpn$sample_S)
MGX_300_cpn$variable = gsub("_Abundance-RPKs", "", MGX_300_cpn$variable)
MGX_300_cpn$taxa <- gsub(".*s__", "s__", MGX_300_cpn$taxa)

MGX_300$MGX_RA_cpn = MGX_300_cpn[match(MGX_300$sample_S, MGX_300_cpn$sample_S), "MGX_RA_cpn"]
MGX_300$MTX_RA_cpn = MTX_300_cpn[match(MGX_300$sample_S, MTX_300_cpn$sample_S), "MTX_RA_cpn"]

MGX_300$MGX_RA_cpn = MGX_300$MGX_RA_cpn$MGX_RA_cpn
MGX_300$MTX_RA_cpn = MGX_300$MTX_RA_cpn$MTX_RA_cpn

MGX_300 <- MGX_300 %>%
        mutate(across(where(is.list), as.character)) %>%  # Convert list columns to characters
        mutate(across(where(is.factor), as.character))  # Convert factors to characters
MGX_300[is.na(MGX_300)] <- 0


# Extract correlation values
MGX_300_filtered <- MGX_300 %>% filter(MTX_RA_humann != 0 & MTX_RA_cpn != 0 & MTX_RA_humann != 1 & MTX_RA_cpn != 1)

cor_test <- cor.test(MGX_300_filtered$MTX_RA_cpn, MGX_300_filtered$MTX_RA_humann, method = "spearman")
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(MTX_RA_humann ~ MTX_RA_cpn, data = MGX_300_filtered)
slope <- round(coef(lm_fit)[2], 2)
intercept <- round(coef(lm_fit)[1], 2)
eqn <- paste("y =", slope, "*x +", intercept)

# Create scatter plot with updated axes
p <- ggplot(MGX_300_filtered, aes(x = MTX_RA_cpn, y = MTX_RA_humann)) +
        geom_point(color = "dark green", alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, fill = "lightgreen", alpha = 0.3) +
        annotate("text", x = max(MGX_300_filtered$MTX_RA_cpn, na.rm = TRUE) * 0.9, 
                 y = max(MGX_300_filtered$MTX_RA_humann, na.rm = TRUE) * 0.9, 
                 label = paste("Spearman r =", cor_val, "\np-value =", p_val), 
                 color = "black", hjust = 1, size = 5) +
        annotate("text", x = max(MGX_300_filtered$MTX_RA_cpn, na.rm = TRUE) * 0.9, 
                 y = min(MGX_300_filtered$MTX_RA_humann, na.rm = TRUE) * 1.1, 
                 label = eqn, color = "black", hjust = 1, size = 5) +
        labs(x = "MTX RA Cpn", y = "MTX RA Humann", title = "Updated Correlation: MTX_RA_cpn vs. MTX_RA_humann") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_sqrt() + scale_y_sqrt()

ggsave("MTX_RA_cpn_vs_MTX_RA_humann.pdf", plot = p, width = 5, height = 5)