all_1 = fread("/Users/leawang/Dropbox (Harvard University)/hutlab/Lea/Manuscript/Cpn60_manuscript/data/Supple_5_iRep.csv")

all_100 = all_1[!all_1$overall_depth < 100, ]
cor_test <- cor.test(all_100$ratio, all_100$humann_MTX_MGX_ratio)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
# Linear regression for equation
lm_fit <- lm(humann_MTX_MGX_ratio ~ ratio, data = all_100)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
#cpn to humann
ggplot(all_100, aes(x = ratio, y = humann_MTX_MGX_ratio)) +
        geom_point(color = "brown", alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_100$ratio), y = max(all_100$humann_MTX_MGX_ratio),
                 label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +
        annotate("text", x = max(all_100$ratio), y = min(all_100$humann_MTX_MGX_ratio),
                 label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MTX/MGX ratio (cpn60)", y = "MTX/MGX ratio (humann)", title = "") +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_sqrt() + scale_x_sqrt()

#cpn tp bPTR
all_100$`bPTR_GC` = as.numeric(as.character(all_100$`bPTR_GC`))
cor_test <- cor.test(all_100$ratio, all_100$`bPTR_GC`)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(`bPTR_GC` ~ ratio, data = all_100)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
p = ggplot(all_100, aes(x = ratio, y = `bPTR_GC`)) +geom_point(color = "brown", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_100$ratio), y = max(all_100$`bPTR-GC`), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +annotate("text", x = max(all_100$ratio), y = min(all_100$`bPTR-GC`),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MTX/MGX ratio (cpn60)", y = "bPTR (by GC skew)", title = "") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
#bPTR to humann
all_100$`bPTR_GC` = as.numeric(as.character(all_100$`bPTR_GC`))
cor_test <- cor.test(all_100$`bPTR_GC`, all_100$humann_MTX_MGX_ratio)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(humann_MTX_MGX_ratio ~ `bPTR_GC`, data = all_100)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(all_100, aes(x = `bPTR_GC`, y = humann_MTX_MGX_ratio)) +geom_point(color = "brown", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_100$`bPTR-GC`), y = max(all_100$humann_MTX_MGX_ratio), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +annotate("text", x = max(all_100$`bPTR-GC`), y = min(all_100$humann_MTX_MGX_ratio),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "bPTR (by GC skew)", y = "MTX/MGX ratio (humann)", title = "") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
#bPTR to cpn
all_100$`bPTR_GC` = as.numeric(as.character(all_100$`bPTR_GC`))
cor_test <- cor.test(all_100$`bPTR_GC`, all_100$ratio)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(ratio ~ `bPTR_GC`, data = all_100)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(all_100, aes(x = `bPTR_GC`, y = ratio)) +geom_point(color = "brown", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_100$`bPTR_GC`), y = max(all_100$ratio), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +annotate("text", x = max(all_100$`bPTR_GC`), y = min(all_100$ratio),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "bPTR (by GC skew)", y = "MTX/MGX ratio (cpn60)", title = "") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()

all_50 = all_1[!all_1$overall_depth< 50, ]
all_50$`bPTR_GC` = as.numeric(as.character(all_50$`bPTR_GC`))
cor_test <- cor.test(all_50$`bPTR_GC`, all_50$ratio)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(ratio ~ `bPTR_GC`, data = all_50)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(all_50, aes(x = `bPTR_GC`, y = ratio)) +geom_point(color = "brown", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_50$`bPTR_GC`), y = max(all_50$ratio), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +annotate("text", x = max(all_50$`bPTR_GC`), y = min(all_50$ratio),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "bPTR (by GC skew)", y = "MTX/MGX ratio (cpn60)", title = "") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()


#humann to cpn
all_100$`bPTR-GC` = as.numeric(as.character(all_100$`bPTR-GC`))
cor_test <- cor.test(all_100$humann_MTX_MGX_ratio, all_100$ratio)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(ratio ~ humann_MTX_MGX_ratio, data = all_100)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(all_100, aes(x = humann_MTX_MGX_ratio, y = ratio)) +geom_point(color = "brown", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_100$humann_MTX_MGX_ratio), y = max(all_100$ratio), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +annotate("text", x = max(all_100$humann_MTX_MGX_ratio), y = min(all_100$ratio),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MTX/MGX ratio (humann)", y = "MTX/MGX ratio (cpn60)", title = "") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
#humann to bPTR
cor_test <- cor.test(all_100$humann_MTX_MGX_ratio, all_100$`bPTR_GC`)
cor_val <- round(cor_test$estimate, 2)
p_val <- round(cor_test$p.value, 4)
lm_fit <- lm(`bPTR_GC` ~ humann_MTX_MGX_ratio, data = all_100)
eqn <- paste("y =", round(coef(lm_fit)[2], 2), "*x +", round(coef(lm_fit)[1], 2))
ggplot(all_100, aes(x = humann_MTX_MGX_ratio, y = `bPTR_GC`)) +geom_point(color = "brown", alpha = 0.5) + geom_smooth(method = "lm", se = TRUE, fill = "lightblue") +
        annotate("text", x = max(all_100$humann_MTX_MGX_ratio), y = max(all_100$`bPTR_GC`), label = paste("r =", cor_val, "\np-value =", p_val), hjust = 1, vjust = 1) +annotate("text", x = max(all_100$humann_MTX_MGX_ratio), y = min(all_100$`bPTR_GC`),label = eqn, hjust = 1, vjust = 0) +
        labs(x = "MTX/MGX ratio (humann)", y = "bPTR (by GC skew)", title = "") +  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +scale_y_sqrt() + scale_x_sqrt()
#barplot #export in 6x 3
ggplot(all_100, aes(x = Phylum, y = ratio, color = Phylum)) +
        geom_jitter(aes(color = Phylum), position=position_jitter(0.2), size = 1, alpha = 0.5) +
        geom_boxplot(aes(color = Phylum), position=position_dodge(), alpha = 0.5)  +
        labs(x = "Phyla", y = "MTX/MGX ratio (cpn60 Uniref90)") +
        theme_bw() +scale_y_sqrt() +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggplot(all_100, aes(x = Phylum, y = humann_MTX_MGX_ratio, color = Phylum)) +
        geom_jitter(aes(color = Phylum), position=position_jitter(0.2), size = 1, alpha = 0.5) +
        geom_boxplot(aes(color = Phylum), position=position_dodge(), alpha = 0.5)  +
        labs(x = "Phyla", y = "MTX/MGC ratio (humann)") +
        theme_bw() +scale_y_sqrt() +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
