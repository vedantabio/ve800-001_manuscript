# =============================================================================
# 04_serum_SCFA_analysis.R
# -----------------------------------------------------------------------------
# Serum short-chain fatty acids (SCFA) in melanoma patients: change over time and differences
# by Week-8 response (SD vs PD) and by VE800 engraftment (HE vs LE).
#
# Statistics (per analyte, nlme::lme with a random intercept per subject):
#   Section 1:   log10(Abun) ~ Time                 all subjects, all visits
#   Sections 2-3: log10(Abun) ~ Time * Group        Baseline, D1, D14, D30
#     emmeans contrasts: (a) each visit vs Baseline within each group
#                        (b) group difference at each visit
#   P-values are BH-adjusted across analytes.
#
# Heatmaps show "global baseline" z-scores: log(Abun) centred/scaled by the
# BASELINE mean and SD of that analyte (pooled over all subjects).
# Hexanoic acid is excluded from the group models (sections 2-3).
#
# Input:  Processed_data/SCFA_serum_melanoma.csv
# Output: Results/04_serum_SCFA_analysis/  (see README.md)
#   SCFA_overtime_boxplot.pdf, Tab_LME_SCFA_Time.csv
#   SD_PD/  timeline, LME tables, Heatmap_SCFA_SD_PD*.pdf
#   HE_LE/  timeline, LME tables, Heatmap_SCFA_HE_LE*.pdf
# =============================================================================

source(here::here("R/utils.R"))
suppressPackageStartupMessages({
  library(phyloseq)
  library(purrr)
  library(vegan)
  library(tidyverse)
  library(yingtools2)
  library(ggprism)
  library(ComplexHeatmap)
  library(nlme)
  library(lemon)
  library(emmeans)
  library(circlize)
  library(scales)
})

# Output folder
results_folder <- create_results_dir("Results/04_serum_SCFA_analysis")

# =============================================================================
# 0. Data
# =============================================================================
# Processed serum SCFA data (long format, one row per subject x visit x analyte)
scfa_dt  <-  read.csv(here("Processed_data/SCFA_serum_melanoma.csv"))

# =============================================================================
# 1. Effect of time (all subjects, each visit vs Baseline)
# =============================================================================
# One LME per analyte; each visit compared with Baseline

lme_dt <-  scfa_dt

result_lme <- list()
for (scfa_var in unique(scfa_dt$Analyte)){
  sel_dt <-  lme_dt[lme_dt$Analyte %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- log10(sel_dt$Abun)
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D21","D30","D58","D86"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  mod_bac <- lme(t_Abundance ~  Time  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Analyte <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only Visit
result_lme_dt <- result_lme_dt[grep("Time",result_lme_dt$Var),]

# Group by adjust
result_lme_dt <-  result_lme_dt %>%
  group_by(Var)%>%
  mutate(p_adj = p.adjust(pval, method = "BH"))


# #result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
# result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]

# Save the LME results
write.csv(result_lme_dt,paste(results_folder,paste0("Tab_LME_SCFA_Time.csv"),sep = "/"))

# ---- Box plot of every SCFA over time -----------------------------------------
plot_dt <-  scfa_dt
plot_dt$Time <- factor(plot_dt$Time , levels = c("Baseline","D1","D14","D21","D30","D58","D86"))
plot_dt$log10_Abun <- log10(plot_dt$Abun)

p_scfa <- ggplot(plot_dt, aes(x = Time, y = log10_Abun)) +
  geom_boxplot(aes(y = log10_Abun), alpha = 0.3, outlier.colour = NA,
               position=position_dodge2(width=0.75,preserve = "single"))+
  geom_point(aes(y = log10_Abun),
             position=position_dodge(width=0.75,preserve = "total"),
             alpha = 0.3,size = 2)+
  facet_rep_wrap(~Analyte , scales = "free",ncol = 3)+
  theme_prism()+
  theme(axis.text.x=element_text(size=14,angle = 45, hjust=1,vjust = 1,face="bold"),
        axis.text.y=element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14,face="bold"))+
  xlab("Day")+
  ylab("Log10(Abundance)")

pdf(paste(results_folder,paste0("SCFA_overtime_boxplot.pdf"),sep="/"),height =10, width = 10)
print(p_scfa)
dev.off()



# =============================================================================
# 2. SD vs PD (Week-8 response)
# =============================================================================
results_folder <- create_results_dir("Results/04_serum_SCFA_analysis/SD_PD")

# ---- Sampling timeline coloured by response ---------------------------------
col_resp <- c("#ff5e1e","#0099b0")
names(col_resp) <- c("PD", "SD")

timeline_dt <-  scfa_dt %>% select(Subject_ID,Time,Response_Wk8,Cohort,Engraftment)%>% unique()
timeline_dt$Subject_ID <- as.character(timeline_dt$Subject_ID)

timeline_dt$Timepoint <-  factor(timeline_dt$Time, levels = c("Baseline","D1","D14","D30","D58","D86"))
#levels(timeline_dt$VISIT_NAME) <-  c("Baseline","Wk_2","Wk_4","Wk_8")

p_time <- ggplot(timeline_dt, aes(x = Timepoint,y = Subject_ID,
                                  group = Subject_ID,color = Response_Wk8))+
  geom_point(size = 3,alpha = 0.7)+
  geom_line()+
  scale_color_manual(name = "Response",values = col_resp)+
  #ggforce::facet_col(Cohort ~ Response_Wk8, scales = 'free', space = 'free')+
  facet_rep_wrap(Cohort ~ Response_Wk8 , scales = "free_y",ncol = 2,repeat.tick.labels = "x" )+
  theme_bw()+
  ylab("Subjects")+
  theme_prism()+
  theme(axis.text.x=element_text(size=14,angle = 45, hjust=1,vjust = 1,face="bold"),
        axis.text.y=element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14,face="bold"))

#print(p_time)

pdf(paste(results_folder,paste0("Timeline_SCFA_SD_PD",".pdf"),sep="/"),height =4, width = 6)
print(p_time)
dev.off()



# ---- 2a. LME: Time x Response, emmeans contrasts -----------------------------

lme_dt <-  scfa_dt %>% filter(!Analyte %in% c("Hexanoic acid"))%>%
  filter(Time %in% c("Baseline","D1","D14","D30"))

result_lme <- list()
effect_time_lme <-  list()
group_time_lme <-  list()
for (scfa_var in unique(lme_dt$Analyte)){
  sel_dt <-  lme_dt[lme_dt$Analyte %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- log10(sel_dt$Abun)
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D30"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  sel_dt$Response_Wk8 <- factor(sel_dt$Response_Wk8, levels = c("SD","PD"))
  mod_bac <- lme(t_Abundance ~  Time*Response_Wk8  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Analyte <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
  
  # Time effect compared to baseline in each group
  time_emm <- emmeans(mod_bac, ~ Time | Response_Wk8)
  eff_dt <-  data.frame(pairs(time_emm))
  eff_dt <-  eff_dt[grep("Baseline",eff_dt$contrast),]
  eff_dt$Analyte <- scfa_var
  effect_time_lme[[scfa_var]] <- eff_dt
  
  # Cross-sectional comparison across each time point between two groups:
  group_emm <- emmeans(mod_bac, ~ Response_Wk8 | Time)
  group_eff_dt <-  data.frame(pairs(group_emm))
  group_eff_dt$Analyte <- scfa_var
  group_time_lme[[scfa_var]] <- group_eff_dt
}

# Time comparison compared to baseline across each group
time_comp_dt <-  do.call("rbind",effect_time_lme)
# Adjust pvalues per group per time comparison across all metabolites
time_comp_dt <-  time_comp_dt %>%
  group_by(contrast,Response_Wk8)%>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

# Group comparison across each time point:
group_comp_dt <-  do.call("rbind",group_time_lme)
# Adjust pvalues per time point across all metabolites
group_comp_dt <-  group_comp_dt %>%
  group_by(Time)%>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only interaction
result_lme_dt <- result_lme_dt[grep(":",result_lme_dt$Var),]
# Group by adjust
result_lme_dt <-  result_lme_dt %>%
  group_by(Var)%>%
  mutate(p_adj = p.adjust(pval, method = "BH"))

write.csv(result_lme_dt,paste(results_folder,"Tab_LME_SCFA_Int_SD_PD.csv",sep = "/"))
write.csv(group_comp_dt,paste(results_folder,"Tab_LME_SCFA_SD_Vs_PD_Time.csv",sep = "/"))
write.csv(time_comp_dt,paste(results_folder,"Tab_LME_SCFA_Time_Effect_SD_PD.csv",sep = "/"))

# No significant interactions across all time points compared to baseline


# ---- 2b. Heatmap: z-scores split by visit, SD vs PD statistics ---------------

mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))


ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  pivot_wider(id_cols = c(Analyte,Time),names_from = Subject_ID ,values_from = s_Abun)

# metadata on pvalue:
meta_dt <-  ht_mat %>% select(Analyte,Time) %>% left_join(group_comp_dt)


# Make a heatmap:
mat <-  ht_mat[,3:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte
# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Response_Wk8,Engraftment)%>% unique()
col_ann_dt <-  col_ann_dt[match(colnames(mat), col_ann_dt$Subject_ID),]   # align annotation to heatmap columns

# Row annotation:
row_ann_dt <-  meta_dt
pval <- meta_dt$p.value
padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pvalue = anno_numeric(pval, 
                                               rg = c(min(pval), 1),
                                               x_convert = function(x) ifelse(pval<= .05,1,0), 
                                               labels_format = function(x) sprintf("%.2e", pval),
                                               labels_gp = gpar(fontsize = 8)),
                         padj = anno_numeric(padj, 
                                             rg = c(min(padj), 1),
                                             x_convert = function(x) ifelse(padj<= .05,1,0), 
                                             labels_format = function(x) sprintf("%.2e", padj),
                                             labels_gp = gpar(fontsize = 8)),
                         annotation_name_rot = 0,
                         annotation_name_gp =gpar(fontsize = 10))

ha_left = rowAnnotation(numeric = anno_numeric(meta_dt$estimate, 
                                               labels_format = function(x) sprintf("%.2f",meta_dt$estimate),
                                               bg_gp = gpar(fill = c("#ff5e1e","#0099b0"))), 
                        annotation_name_rot = 0)
col_coef <- colorRamp2(c(min(meta_dt$estimate), 0, max(meta_dt$estimate)), c(c("#ff5e1e","white","#0099b0"))) 
ha_left =  rowAnnotation(Estimate = meta_dt$estimate,
                         col = list(Estimate =  col_coef),
                         show_legend = c(TRUE),
                         show_annotation_name = c(TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        legend_height = unit(4, "cm")))

max_min <-  max(abs(min(meta_dt$estimate)),abs(max(meta_dt$estimate)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(0.5*c(-4, 0, 4),c("#0a75ad","white","#ff5353"))
z_score_min <- round(min(mat,na.rm = T))
z_score_max <-  round(max(mat,na.rm = T))
col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))

split_cols <-  factor(col_ann_dt$Response_Wk8)
split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))
ht  =  Heatmap(mat,name = "Z-score",
               column_split = split_cols,
               row_split =  split_rows,
               #         top_annotation = ha_column,
               right_annotation =  ha_right,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 14,fontface = "bold"),
               column_names_gp = gpar(fontsize = 8,fontface = "bold"),
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = F,
               cluster_row_slices = F,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = T,
               cluster_column_slices = F,
               show_column_dend = F,
               border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                                      legend_direction = "horizontal"),
               row_names_max_width = max_text_width(rownames(mat),
                                                    gp = gpar(fontsize = 10)),
               column_names_max_height = max_text_width(colnames(mat),
                                                        gp = gpar(fontsize = 10)))
pdf(paste(results_folder,paste0("Heatmap_SCFA_SD_PD.pdf"),sep = "/"), height = 6, width = 8,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# ---- 2c. Heatmap: visit vs Baseline within SD (page 1) and PD (page 2) -------



mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))

ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  filter(Response_Wk8 %in% c("SD"))%>%
  pivot_wider(id_cols = c(Analyte), names_from = c(Subject_ID,Time) ,values_from = s_Abun)

# metadata on pvalue:
meta_dt <-  time_comp_dt %>% filter(Response_Wk8 %in% "SD")

# Make a heatmap:
mat <-  ht_mat[,2:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte
# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Time,Response_Wk8,Engraftment) %>%
  filter(Response_Wk8 %in% c("SD"))%>%
  unique() %>%
  mutate(colID = paste0(Subject_ID,"_",Time))
col_ann_dt <-  col_ann_dt[match(colnames(mat), col_ann_dt$colID),]   # align annotation to heatmap columns


split_cols <-  factor(col_ann_dt$Time, levels = c("Baseline","D1","D14","D30"))

colnames(mat) <- gsub("_.*","",colnames(mat))
#split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))

# Row annotation:
row_ann_dt <-  meta_dt %>% filter(Response_Wk8 %in% "SD") %>% data.frame()
row_ann_dt$Time <- gsub("Baseline - ","",row_ann_dt$contrast)
row_ann_dt <- row_ann_dt %>%
  select(Time,estimate,Analyte,p.value,p_adj)%>%
  pivot_wider(id_cols = Analyte,names_from = Time,values_from = c(estimate,p.value,p_adj))

pval_D1 <- row_ann_dt$p.value_D1
pval_D14 <- row_ann_dt$p.value_D14
pval_D30 <- row_ann_dt$p.value_D30
padj_D1 <- row_ann_dt$p_adj_D1
padj_D14 <- row_ann_dt$p_adj_D14
padj_D30 <- row_ann_dt$p_adj_D30



#padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pval_D1 = anno_numeric(pval_D1,
                                                rg = c(min(pval_D1), 1),
                                                x_convert = function(x) ifelse(pval_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", pval_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         pval_D14 = anno_numeric(pval_D14,
                                                 rg = c(min(pval_D14), 1),
                                                 x_convert = function(x) ifelse(pval_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         pval_D30 = anno_numeric(pval_D30,
                                                 rg = c(min(pval_D30), 1),
                                                 x_convert = function(x) ifelse(pval_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D1 = anno_numeric(padj_D1,
                                                rg = c(min(padj_D1), 1),
                                                x_convert = function(x) ifelse(padj_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", padj_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         padj_D14 = anno_numeric(padj_D14,
                                                 rg = c(min(padj_D14), 1),
                                                 x_convert = function(x) ifelse(padj_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D30 = anno_numeric(padj_D30,
                                                 rg = c(min(padj_D30), 1),
                                                 x_convert = function(x) ifelse(padj_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         
                         
                         annotation_name_rot = 45,
                         annotation_name_gp =gpar(fontsize = 5),
                         width = unit(4, "cm"))

col_coef <- colorRamp2(c(min(meta_dt$estimate), 0, max(meta_dt$estimate)), c(muted("red"),"white",muted("blue")))
ha_left =  rowAnnotation(Estimate_D1 = anno_simple(row_ann_dt$estimate_D1,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D1< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D14 =anno_simple(row_ann_dt$estimate_D14,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D14< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D30 = anno_simple(row_ann_dt$estimate_D30,
                                                    col = col_coef,
                                                    pch = ifelse(row_ann_dt$p_adj_D30< 0.05,"*",NA),
                                                    pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(F,F,F),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         width = unit(1, "cm"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 8,fontface = "bold" ),
                                                        legend_height = unit(0.5, "cm")))
lgd_pvalue = Legend(title = "Coeff", col_fun = col_coef)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")
max_min <-  max(abs(min(meta_dt$estimate)),abs(max(meta_dt$estimate)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(2*c(-4, 0, 4),c("#0a75ad","white","#ff5353"))
col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))
ht1 <- Heatmap(mat,name = "Z-score",
        column_split = split_cols,
        #row_split =  split_rows,
        #         top_annotation = ha_column,
        #right_annotation =  ha_right,
        right_annotation = ha_left,
        col = col_mat,
        row_names_side = "left",
        row_gap = unit(2, "mm"),
        row_title_gp = gpar(fontsize = 14,fontface = "bold"),
        row_title_rot = 0,
        column_title_gp = gpar(fontsize = 14,fontface = "bold"),
        column_names_gp = gpar(fontsize = 8,fontface = "bold"),
        row_names_gp = gpar(fontsize = 10),
        column_title_rot = 0,
        cluster_rows = T,
        cluster_row_slices = F,
        show_parent_dend_line = F,
        show_row_dend = F,
        cluster_columns = F,
        cluster_column_slices = F,
        show_column_dend = F,
        show_column_names = F,
        border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                               legend_direction = "horizontal"),
        row_names_max_width = max_text_width(rownames(mat),
                                             gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(mat),
                                                 gp = gpar(fontsize = 10)))
# pdf(paste(results_folder,paste0("Heatmap_SCFA_SD_Time.pdf"),sep = "/"), height = 3, width = 8,useDingbats = F)
if (interactive()) draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",
     annotation_legend_list = list(lgd_pvalue, lgd_sig),
     legend_grouping = "original")
# dev.off()

# PD panel (page 2):
mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))

ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  filter(Response_Wk8 %in% c("PD"))%>%
  pivot_wider(id_cols = c(Analyte), names_from = c(Subject_ID,Time) ,values_from = s_Abun)

# metadata on pvalue:
meta_dt <-  time_comp_dt %>% filter(Response_Wk8 %in% "PD")

# Make a heatmap:
mat <-  ht_mat[,2:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte
# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Time,Response_Wk8,Engraftment) %>%
  filter(Response_Wk8 %in% c("PD"))%>%
  unique() %>%
  mutate(colID = paste0(Subject_ID,"_",Time))
col_ann_dt <-  col_ann_dt[match(colnames(mat), col_ann_dt$colID),]   # align annotation to heatmap columns


split_cols <-  factor(col_ann_dt$Time, levels = c("Baseline","D1","D14","D30"))

colnames(mat) <- gsub("_.*","",colnames(mat))
#split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))

# Row annotation:
row_ann_dt <-  meta_dt %>% filter(Response_Wk8 %in% "PD") %>% data.frame()
row_ann_dt$Time <- gsub("Baseline - ","",row_ann_dt$contrast)
row_ann_dt <- row_ann_dt %>%
  select(Time,estimate,Analyte,p.value,p_adj)%>%
  pivot_wider(id_cols = Analyte,names_from = Time,values_from = c(estimate,p.value,p_adj))

pval_D1 <- row_ann_dt$p.value_D1
pval_D14 <- row_ann_dt$p.value_D14
pval_D30 <- row_ann_dt$p.value_D30
padj_D1 <- row_ann_dt$p_adj_D1
padj_D14 <- row_ann_dt$p_adj_D14
padj_D30 <- row_ann_dt$p_adj_D30



#padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pval_D1 = anno_numeric(pval_D1,
                                                rg = c(min(pval_D1), 1),
                                                x_convert = function(x) ifelse(pval_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", pval_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         pval_D14 = anno_numeric(pval_D14,
                                                 rg = c(min(pval_D14), 1),
                                                 x_convert = function(x) ifelse(pval_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         pval_D30 = anno_numeric(pval_D30,
                                                 rg = c(min(pval_D30), 1),
                                                 x_convert = function(x) ifelse(pval_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D1 = anno_numeric(padj_D1,
                                                rg = c(min(padj_D1), 1),
                                                x_convert = function(x) ifelse(padj_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", padj_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         padj_D14 = anno_numeric(padj_D14,
                                                 rg = c(min(padj_D14), 1),
                                                 x_convert = function(x) ifelse(padj_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D30 = anno_numeric(padj_D30,
                                                 rg = c(min(padj_D30), 1),
                                                 x_convert = function(x) ifelse(padj_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         
                         
                         annotation_name_rot = 45,
                         annotation_name_gp =gpar(fontsize = 5),
                         width = unit(4, "cm"))

col_coef <- colorRamp2(c(min(meta_dt$estimate), 0, max(meta_dt$estimate)), c(muted("red"),"white",muted("blue")))
ha_left =  rowAnnotation(Estimate_D1 = anno_simple(row_ann_dt$estimate_D1,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D1< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D14 =anno_simple(row_ann_dt$estimate_D14,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D14< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D30 = anno_simple(row_ann_dt$estimate_D30,
                                                    col = col_coef,
                                                    pch = ifelse(row_ann_dt$p_adj_D30< 0.05,"*",NA),
                                                    pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(F,F,F),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         width = unit(1, "cm"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 8,fontface = "bold" ),
                                                        legend_height = unit(1, "cm")))

lgd_pvalue = Legend(title = "Coeff", col_fun = col_coef)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")


max_min <-  max(abs(min(meta_dt$estimate)),abs(max(meta_dt$estimate)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(0.5*c(-4, 0, 4),c("#0a75ad","white","#ff5353"))

col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))

ht2  =  Heatmap(mat,name = "Z-score",
                column_split = split_cols,
                #row_split =  split_rows,
                #         top_annotation = ha_column,
                #right_annotation =  ha_right,
                right_annotation = ha_left,
                col = col_mat,
                row_names_side = "left",
                row_gap = unit(2, "mm"),
                row_title_gp = gpar(fontsize = 14,fontface = "bold"),
                row_title_rot = 0,
                column_title_gp = gpar(fontsize = 14,fontface = "bold"),
                column_names_gp = gpar(fontsize = 8,fontface = "bold"),
                row_names_gp = gpar(fontsize = 10),
                column_title_rot = 0,
                cluster_rows = T,
                cluster_row_slices = F,
                show_parent_dend_line = F,
                show_row_dend = F,
                cluster_columns = F,
                cluster_column_slices = F,
                show_column_dend = F,
                show_column_names = F,
                border = T,
                heatmap_legend_param = list(legend_width = unit(4, "cm"),
                                            legend_direction = "horizontal"),
                row_names_max_width = max_text_width(rownames(mat),
                                                     gp = gpar(fontsize = 15)),
                column_names_max_height = max_text_width(colnames(mat),
                                                         gp = gpar(fontsize = 10)))
# pdf(paste(results_folder,paste0("Heatmap_SCFA_PD_Time.pdf"),sep = "/"), height = 2.5, width = 12,useDingbats = F)
if (interactive()) draw(ht2,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original",
     annotation_legend_list = list(lgd_pvalue, lgd_sig))


pdf(paste(results_folder,paste0("Heatmap_SCFA_SD_PD_Time_norm_baseline_zscore.pdf"),sep = "/"), height = 2.5, width = 8,useDingbats = F)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original",
     annotation_legend_list = list(lgd_pvalue, lgd_sig))
#ht <-  ht1 %v% ht2
draw(ht2,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original",
     annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()


# ---- 2d. Heatmap: Time x Response interaction coefficients -------------------

# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Time,Response_Wk8,Engraftment) %>%
  unique() %>%
  mutate(colID = paste0(Subject_ID,"_",Time))

col_ann_dt <- col_ann_dt %>% arrange(Time,Response_Wk8)


mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))


ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  pivot_wider(id_cols = c(Analyte), names_from = c(Subject_ID,Time) ,values_from = s_Abun)


ht_mat <- ht_mat[,c("Analyte",col_ann_dt$colID),]

# metadata on pvalue:
meta_dt <-  result_lme_dt 

# Make a heatmap:
mat <-  ht_mat[,2:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte

split_cols <-  factor(col_ann_dt$Time, levels = c("Baseline","D1","D14","D30"))

colnames(mat) <- gsub("_.*","",colnames(mat))
#split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))

# Row annotation:
row_ann_dt <-  meta_dt %>% data.frame()
names(row_ann_dt)[c(1,5)] <- c("estimate","p.value")
row_ann_dt$Time <- gsub(":.*","",row_ann_dt$Var)
row_ann_dt$Time <- gsub("Time","",row_ann_dt$Time)
row_ann_dt <- row_ann_dt %>%
  select(Time,estimate,Analyte,p.value,p_adj)%>%
  pivot_wider(id_cols = Analyte,names_from = Time,values_from = c(estimate,p.value,p_adj))

row_ann_dt <- row_ann_dt[match(ht_mat$Analyte, row_ann_dt$Analyte),]   # align annotation to heatmap rows


pval_D1 <- row_ann_dt$p.value_D1
pval_D14 <- row_ann_dt$p.value_D14
pval_D30 <- row_ann_dt$p.value_D30
padj_D1 <- row_ann_dt$p_adj_D1
padj_D14 <- row_ann_dt$p_adj_D14
padj_D30 <- row_ann_dt$p_adj_D30





#padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pval_D1 = anno_numeric(pval_D1,
                                                rg = c(min(pval_D1), 1),
                                                x_convert = function(x) ifelse(pval_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", pval_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         pval_D14 = anno_numeric(pval_D14,
                                                 rg = c(min(pval_D14), 1),
                                                 x_convert = function(x) ifelse(pval_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         pval_D30 = anno_numeric(pval_D30,
                                                 rg = c(min(pval_D30), 1),
                                                 x_convert = function(x) ifelse(pval_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D1 = anno_numeric(padj_D1,
                                                rg = c(min(padj_D1), 1),
                                                x_convert = function(x) ifelse(padj_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", padj_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         padj_D14 = anno_numeric(padj_D14,
                                                 rg = c(min(padj_D14), 1),
                                                 x_convert = function(x) ifelse(padj_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D30 = anno_numeric(padj_D30,
                                                 rg = c(min(padj_D30), 1),
                                                 x_convert = function(x) ifelse(padj_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         
                         annotation_name_rot = 45,
                         annotation_name_gp =gpar(fontsize = 5),
                         width = unit(4, "cm"))

ha_col <- HeatmapAnnotation(Resp = col_ann_dt$Response_Wk8,col = list(Resp = col_resp))

col_coef <- colorRamp2(c(min(row_ann_dt$estimate_D1), 0, max(row_ann_dt$estimate_D1)), c(muted("red"),"white",muted("blue")))
ha_left =  rowAnnotation(Estimate_D1 = anno_simple(row_ann_dt$estimate_D1,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D1< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D14 =anno_simple(row_ann_dt$estimate_D14,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D14< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D30 = anno_simple(row_ann_dt$estimate_D30,
                                                    col = col_coef,
                                                    pch = ifelse(row_ann_dt$p_adj_D30< 0.05,"*",NA),
                                                    pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(F,F,F),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         width = unit(1, "cm"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 8,fontface = "bold" ),
                                                        legend_height = unit(0.5, "cm")))

lgd_pvalue = Legend(title = "Coeff", col_fun = col_coef)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")
max_min <-  max(abs(min(meta_dt$Value)),abs(max(meta_dt$Value)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(0.5*c(-4, 0, 4),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(c(min(mat,na.rm = T), 0, max(mat,na.rm = T)),c("#0a75ad","white","#ff5353"))
col_mat <-  colorRamp2(c(round(min(mat,na.rm = T)), 0, round(max(mat,na.rm = T))),c("#0a75ad","white","#ff5353"))

ht1  =  Heatmap(mat,name = "Z-score",
                column_split = split_cols,
                #row_split =  split_rows,
                top_annotation = ha_col,
                right_annotation =  ha_left,
                #left_annotation = ha_left,
                col = col_mat,
                row_names_side = "left",
                row_gap = unit(2, "mm"),
                row_title_gp = gpar(fontsize = 14,fontface = "bold"),
                row_title_rot = 0,
                column_title_gp = gpar(fontsize = 14,fontface = "bold"),
                column_names_gp = gpar(fontsize = 8,fontface = "bold"),
                row_names_gp = gpar(fontsize = 10),
                column_title_rot = 0,
                cluster_rows = T,
                cluster_row_slices = F,
                show_parent_dend_line = F,
                show_row_dend = F,
                cluster_columns = F,
                cluster_column_slices = F,
                show_column_dend = F,
                show_column_names = F,
                border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                                       legend_direction = "horizontal"),
                row_names_max_width = max_text_width(rownames(mat),
                                                     gp = gpar(fontsize = 15)),
                column_names_max_height = max_text_width(colnames(mat),
                                                         gp = gpar(fontsize = 10)))
pdf(paste(results_folder,paste0("Heatmap_SCFA_SD_PD_Interaction.pdf"),sep = "/"), height = 3, width = 12,useDingbats = F)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",
     annotation_legend_list = list(lgd_pvalue, lgd_sig),
     legend_grouping = "original")
dev.off()


# =============================================================================
# 3. HE vs LE (VE800 engraftment) - same analyses as section 2
# =============================================================================
results_folder <- create_results_dir("Results/04_serum_SCFA_analysis/HE_LE")


# ---- Sampling timeline coloured by engraftment ------------------------------
col_eng <- c("#009031","#7f1f98")
names(col_eng) <- c("HE", "LE")

p_time <- ggplot(timeline_dt, aes(x = Timepoint,y = Subject_ID,
                                  group = Subject_ID,color = Engraftment))+
  geom_point(size = 3,alpha = 0.7)+
  geom_line()+
  scale_color_manual(name = "Engraftment",values = col_eng)+
  #ggforce::facet_col(Cohort ~ Response_Wk8, scales = 'free', space = 'free')+
  facet_rep_wrap(Cohort ~ Engraftment , scales = "free_y",ncol = 2,repeat.tick.labels = "x" )+
  theme_bw()+
  ylab("Subjects")+
  theme_prism()+
  theme(axis.text.x=element_text(size=14,angle = 45, hjust=1,vjust = 1,face="bold"),
        axis.text.y=element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14,face="bold"))

#print(p_time)

pdf(paste(results_folder,paste0("Timeline_SCFA_HE_LE",".pdf"),sep="/"),height =4, width = 6)
print(p_time)
dev.off()


# ---- 3a. LME: Time x Engraftment, emmeans contrasts --------------------------

lme_dt <-  scfa_dt %>% filter(!Analyte %in% c("Hexanoic acid"))%>%
  filter(Time %in% c("Baseline","D1","D14","D30"))

result_lme <- list()
effect_time_lme <-  list()
group_time_lme <-  list()
for (scfa_var in unique(lme_dt$Analyte)){
  sel_dt <-  lme_dt[lme_dt$Analyte %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- log10(sel_dt$Abun)
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D30"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  sel_dt$Engraftment <- factor(sel_dt$Engraftment, levels = c("LE","HE"))
  mod_bac <- lme(t_Abundance ~  Time*Engraftment  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Analyte <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
  
  # Time effect compared to baseline in each group
  time_emm <- emmeans(mod_bac, ~ Time | Engraftment)
  eff_dt <-  data.frame(pairs(time_emm))
  eff_dt <-  eff_dt[grep("Baseline",eff_dt$contrast),]
  eff_dt$Analyte <- scfa_var
  effect_time_lme[[scfa_var]] <- eff_dt
  
  # Cross-sectional comparison across each time point between two groups:
  group_emm <- emmeans(mod_bac, ~ Engraftment | Time)
  group_eff_dt <-  data.frame(pairs(group_emm))
  group_eff_dt$Analyte <- scfa_var
  group_time_lme[[scfa_var]] <- group_eff_dt
}

# Time comparison compared to baseline across each group
time_comp_dt <-  do.call("rbind",effect_time_lme)
# Adjust pvalues per group per time comparison across all metabolites
time_comp_dt <-  time_comp_dt %>%
  group_by(contrast,Engraftment)%>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

# Group comparison across each time point:
group_comp_dt <-  do.call("rbind",group_time_lme)
# Adjust pvalues per time point across all metabolites
group_comp_dt <-  group_comp_dt %>%
  group_by(Time)%>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only interaction
result_lme_dt <- result_lme_dt[grep(":",result_lme_dt$Var),]
# Group by adjust
result_lme_dt <-  result_lme_dt %>%
  group_by(Var)%>%
  mutate(p_adj = p.adjust(pval, method = "BH"))

write.csv(result_lme_dt,paste(results_folder,"Tab_LME_SCFA_Int_LE_HE.csv",sep = "/"))
write.csv(group_comp_dt,paste(results_folder,"Tab_LME_SCFA_LE_Vs_HE_Time.csv",sep = "/"))
write.csv(time_comp_dt,paste(results_folder,"Tab_LME_SCFA_Time_Effect_LE_HE.csv",sep = "/"))



# ---- 3b. Heatmap: z-scores split by visit, HE vs LE statistics ---------------
mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))


ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  pivot_wider(id_cols = c(Analyte,Time),names_from = Subject_ID ,values_from = s_Abun)

# metadata on pvalue:
meta_dt <-  ht_mat %>% select(Analyte,Time) %>% left_join(group_comp_dt)


# Make a heatmap:
mat <-  ht_mat[,3:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte
# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Response_Wk8,Engraftment)%>% unique()
col_ann_dt <-  col_ann_dt[match(colnames(mat), col_ann_dt$Subject_ID),]   # align annotation to heatmap columns

# Row annotation:
row_ann_dt <-  meta_dt
pval <- meta_dt$p.value
padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pvalue = anno_numeric(pval, 
                                         rg = c(min(pval), 1),
                                         x_convert = function(x) ifelse(pval<= .05,1,0), 
                                         labels_format = function(x) sprintf("%.2e", pval),
                                         labels_gp = gpar(fontsize = 8)),
                   padj = anno_numeric(padj, 
                                         rg = c(min(padj), 1),
                                         x_convert = function(x) ifelse(padj<= .05,1,0), 
                                         labels_format = function(x) sprintf("%.2e", padj),
                                          labels_gp = gpar(fontsize = 8)),
                   annotation_name_rot = 0,
                   annotation_name_gp =gpar(fontsize = 10))

ha_left = rowAnnotation(numeric = anno_numeric(meta_dt$estimate, 
                                               labels_format = function(x) sprintf("%.2f",meta_dt$estimate),
                                          bg_gp = gpar(fill = c("#009031","#7f1f98"))), 
                   annotation_name_rot = 0)
col_coef <- colorRamp2(c(min(meta_dt$estimate), 0, max(meta_dt$estimate)), c(c("#009031","white","#7f1f98"))) 
ha_left =  rowAnnotation(Estimate = meta_dt$estimate,
                         col = list(Estimate =  col_coef),
                         show_legend = c(TRUE),
                         show_annotation_name = c(TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 12,fontface = "bold" ),
                                                        legend_height = unit(4, "cm")))

max_min <-  max(abs(min(meta_dt$estimate)),abs(max(meta_dt$estimate)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
z_score_min <- round(min(mat,na.rm = T))
z_score_max <-  round(max(mat,na.rm = T))
col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))

split_cols <-  factor(col_ann_dt$Engraftment)
split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))
ht  =  Heatmap(mat,name = "Z-score",
               column_split = split_cols,
               row_split =  split_rows,
               #         top_annotation = ha_column,
               right_annotation =  ha_right,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 14,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 14,fontface = "bold"),
               column_names_gp = gpar(fontsize = 8,fontface = "bold"),
               row_names_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = F,
               cluster_row_slices = F,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = T,
               cluster_column_slices = F,
               show_column_dend = F,
               border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                                      legend_direction = "horizontal"),
               row_names_max_width = max_text_width(rownames(mat),
                                     gp = gpar(fontsize = 10)),
              column_names_max_height = max_text_width(colnames(mat),
                                         gp = gpar(fontsize = 10)))
pdf(paste(results_folder,paste0("Heatmap_SCFA_HE_LE.pdf"),sep = "/"), height = 6, width = 8,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
dev.off()


# ---- 3c. Heatmap: visit vs Baseline within HE (page 1) and LE (page 2) -------


mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))



ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  filter(Engraftment %in% c("HE"))%>%
  pivot_wider(id_cols = c(Analyte), names_from = c(Subject_ID,Time) ,values_from = s_Abun)

# metadata on pvalue:
meta_dt <-  time_comp_dt %>% filter(Engraftment %in% "HE")

# Make a heatmap:
mat <-  ht_mat[,2:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte
# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Time,Response_Wk8,Engraftment) %>%
  filter(Engraftment %in% c("HE"))%>%
  unique() %>%
  mutate(colID = paste0(Subject_ID,"_",Time))
col_ann_dt <-  col_ann_dt[match(colnames(mat), col_ann_dt$colID),]   # align annotation to heatmap columns


split_cols <-  factor(col_ann_dt$Time, levels = c("Baseline","D1","D14","D30"))

colnames(mat) <- gsub("_.*","",colnames(mat))
#split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))

# Row annotation:
row_ann_dt <-  meta_dt %>% filter(Engraftment %in% "HE") %>% data.frame()
row_ann_dt$Time <- gsub("Baseline - ","",row_ann_dt$contrast)
row_ann_dt <- row_ann_dt %>%
  select(Time,estimate,Analyte,p.value,p_adj)%>%
  pivot_wider(id_cols = Analyte,names_from = Time,values_from = c(estimate,p.value,p_adj))

pval_D1 <- row_ann_dt$p.value_D1
pval_D14 <- row_ann_dt$p.value_D14
pval_D30 <- row_ann_dt$p.value_D30
padj_D1 <- row_ann_dt$p_adj_D1
padj_D14 <- row_ann_dt$p_adj_D14
padj_D30 <- row_ann_dt$p_adj_D30



#padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pval_D1 = anno_numeric(pval_D1,
                                                rg = c(min(pval_D1), 1),
                                                x_convert = function(x) ifelse(pval_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", pval_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         pval_D14 = anno_numeric(pval_D14,
                                                 rg = c(min(pval_D14), 1),
                                                 x_convert = function(x) ifelse(pval_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         pval_D30 = anno_numeric(pval_D30,
                                                 rg = c(min(pval_D30), 1),
                                                 x_convert = function(x) ifelse(pval_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D1 = anno_numeric(padj_D1,
                                                rg = c(min(padj_D1), 1),
                                                x_convert = function(x) ifelse(padj_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", padj_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         padj_D14 = anno_numeric(padj_D14,
                                                 rg = c(min(padj_D14), 1),
                                                 x_convert = function(x) ifelse(padj_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D30 = anno_numeric(padj_D30,
                                                 rg = c(min(padj_D30), 1),
                                                 x_convert = function(x) ifelse(padj_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         
                         
                         annotation_name_rot = 45,
                         annotation_name_gp =gpar(fontsize = 5),
                         width = unit(4, "cm"))

col_coef <- colorRamp2(c(min(meta_dt$estimate), 0, max(meta_dt$estimate)), c(muted("red"),"white",muted("blue")))
ha_left =  rowAnnotation(Estimate_D1 = anno_simple(row_ann_dt$estimate_D1,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D1< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D14 =anno_simple(row_ann_dt$estimate_D14,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D14< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D30 = anno_simple(row_ann_dt$estimate_D30,
                                                    col = col_coef,
                                                    pch = ifelse(row_ann_dt$p_adj_D30< 0.05,"*",NA),
                                                    pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(F,F,F),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         width = unit(1, "cm"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 8,fontface = "bold" ),
                                                        legend_height = unit(0.5, "cm")))
lgd_pvalue = Legend(title = "Coeff", col_fun = col_coef)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")
max_min <-  max(abs(min(meta_dt$estimate)),abs(max(meta_dt$estimate)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(0.5*c(-4, 0, 4),c("#0a75ad","white","#ff5353"))
col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))

ht1  =  Heatmap(mat,name = "Z-score",
                column_split = split_cols,
                #row_split =  split_rows,
                #         top_annotation = ha_column,
                #right_annotation =  ha_right,
                right_annotation = ha_left,
                col = col_mat,
                row_names_side = "left",
                row_gap = unit(2, "mm"),
                row_title_gp = gpar(fontsize = 14,fontface = "bold"),
                row_title_rot = 0,
                column_title_gp = gpar(fontsize = 14,fontface = "bold"),
                column_names_gp = gpar(fontsize = 8,fontface = "bold"),
                row_names_gp = gpar(fontsize = 10),
                column_title_rot = 0,
                cluster_rows = T,
                cluster_row_slices = F,
                show_parent_dend_line = F,
                show_row_dend = F,
                cluster_columns = F,
                cluster_column_slices = F,
                show_column_dend = F,
                show_column_names = F,
                border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                                       legend_direction = "horizontal"),
                row_names_max_width = max_text_width(rownames(mat),
                                                     gp = gpar(fontsize = 15)),
                column_names_max_height = max_text_width(colnames(mat),
                                                         gp = gpar(fontsize = 10)))
# pdf(paste(results_folder,paste0("Heatmap_SCFA_HE_Time.pdf"),sep = "/"), height = 2.5, width = 12,useDingbats = F)
if (interactive()) draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",
     annotation_legend_list = list(lgd_pvalue, lgd_sig),
     legend_grouping = "original")
# dev.off()


# LE panel (page 2):
mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))

ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  filter(Engraftment %in% c("LE"))%>%
  pivot_wider(id_cols = c(Analyte), names_from = c(Subject_ID,Time) ,values_from = s_Abun)

# metadata on pvalue:
meta_dt <-  time_comp_dt %>% filter(Engraftment %in% "LE")

# Make a heatmap:
mat <-  ht_mat[,2:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte
# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Time,Response_Wk8,Engraftment) %>%
  filter(Engraftment %in% c("LE"))%>%
  unique() %>%
  mutate(colID = paste0(Subject_ID,"_",Time))
col_ann_dt <-  col_ann_dt[match(colnames(mat), col_ann_dt$colID),]   # align annotation to heatmap columns


split_cols <-  factor(col_ann_dt$Time, levels = c("Baseline","D1","D14","D30"))

colnames(mat) <- gsub("_.*","",colnames(mat))
#split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))

# Row annotation:
row_ann_dt <-  meta_dt %>% filter(Engraftment %in% "LE") %>% data.frame()
row_ann_dt$Time <- gsub("Baseline - ","",row_ann_dt$contrast)
row_ann_dt <- row_ann_dt %>%
  select(Time,estimate,Analyte,p.value,p_adj)%>%
  pivot_wider(id_cols = Analyte,names_from = Time,values_from = c(estimate,p.value,p_adj))

pval_D1 <- row_ann_dt$p.value_D1
pval_D14 <- row_ann_dt$p.value_D14
pval_D30 <- row_ann_dt$p.value_D30
padj_D1 <- row_ann_dt$p_adj_D1
padj_D14 <- row_ann_dt$p_adj_D14
padj_D30 <- row_ann_dt$p_adj_D30



#padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pval_D1 = anno_numeric(pval_D1,
                                                rg = c(min(pval_D1), 1),
                                                x_convert = function(x) ifelse(pval_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", pval_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         pval_D14 = anno_numeric(pval_D14,
                                                 rg = c(min(pval_D14), 1),
                                                 x_convert = function(x) ifelse(pval_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         pval_D30 = anno_numeric(pval_D30,
                                                 rg = c(min(pval_D30), 1),
                                                 x_convert = function(x) ifelse(pval_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D1 = anno_numeric(padj_D1,
                                                rg = c(min(padj_D1), 1),
                                                x_convert = function(x) ifelse(padj_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", padj_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         padj_D14 = anno_numeric(padj_D14,
                                                 rg = c(min(padj_D14), 1),
                                                 x_convert = function(x) ifelse(padj_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D30 = anno_numeric(padj_D30,
                                                 rg = c(min(padj_D30), 1),
                                                 x_convert = function(x) ifelse(padj_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         
                         
                         annotation_name_rot = 45,
                         annotation_name_gp =gpar(fontsize = 5),
                         width = unit(4, "cm"))
col_coef <- colorRamp2(c(min(meta_dt$estimate), 0, max(meta_dt$estimate)), c(muted("red"),"white",muted("blue")))
ha_left =  rowAnnotation(Estimate_D1 = anno_simple(row_ann_dt$estimate_D1,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D1< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D14 =anno_simple(row_ann_dt$estimate_D14,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D14< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D30 = anno_simple(row_ann_dt$estimate_D30,
                                                    col = col_coef,
                                                    pch = ifelse(row_ann_dt$p_adj_D30< 0.05,"*",NA),
                                                    pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(F,F,F),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         width = unit(1, "cm"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 8,fontface = "bold" ),
                                                        legend_height = unit(1, "cm")))

lgd_pvalue = Legend(title = "Coeff", col_fun = col_coef)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")


max_min <-  max(abs(min(meta_dt$estimate)),abs(max(meta_dt$estimate)))
max_min <- round(max_min,digits = 1)
#col_mat <-  colorRamp2(10*c(-max_min, 0, max_min),c("#0a75ad","white","#ff5353"))
#col_mat <-  colorRamp2(0.5*c(-4, 0, 4),c("#0a75ad","white","#ff5353"))
col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))

ht2  =  Heatmap(mat,name = "Z-score",
                column_split = split_cols,
                #row_split =  split_rows,
                #         top_annotation = ha_column,
                #right_annotation =  ha_right,
                right_annotation = ha_left,
                col = col_mat,
                row_names_side = "left",
                row_gap = unit(2, "mm"),
                row_title_gp = gpar(fontsize = 14,fontface = "bold"),
                row_title_rot = 0,
                column_title_gp = gpar(fontsize = 14,fontface = "bold"),
                column_names_gp = gpar(fontsize = 8,fontface = "bold"),
                row_names_gp = gpar(fontsize = 10),
                column_title_rot = 0,
                cluster_rows = T,
                cluster_row_slices = F,
                show_parent_dend_line = F,
                show_row_dend = F,
                cluster_columns = F,
                cluster_column_slices = F,
                show_column_dend = F,
                show_column_names = F,
                border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                                       legend_direction = "horizontal"),
                row_names_max_width = max_text_width(rownames(mat),
                                                     gp = gpar(fontsize = 15)),
                column_names_max_height = max_text_width(colnames(mat),
                                                         gp = gpar(fontsize = 10)))
# pdf(paste(results_folder,paste0("Heatmap_SCFA_HE_Time.pdf"),sep = "/"), height = 2.5, width = 8,useDingbats = F)
#draw(ht2,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original",
#      annotation_legend_list = list(lgd_pvalue, lgd_sig))
# dev.off()

pdf(paste(results_folder,paste0("Heatmap_SCFA_HE_LE_Time.pdf"),sep = "/"), height = 2.5, width = 8,useDingbats = F)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original",
     annotation_legend_list = list(lgd_pvalue, lgd_sig))
#ht <-  ht1 %v% ht2
draw(ht2,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original",
     annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()



# ---- 3d. Heatmap: Time x Engraftment interaction coefficients ----------------

# Column annotation:
col_ann_dt <-  lme_dt %>% select(Subject_ID,Time,Response_Wk8,Engraftment) %>%
  unique() %>%
  mutate(colID = paste0(Subject_ID,"_",Time))

col_ann_dt <- col_ann_dt %>% arrange(Time,Engraftment)


mu_dt <- lme_dt %>%
  filter(Time %in% "Baseline")%>%
  group_by(Analyte)%>%
  summarise(mu_baseline =  mean(log(Abun)),
            sd_baseline = sd(log(Abun)))


ht_mat <-  lme_dt %>%
  left_join(mu_dt)%>%
  group_by(Analyte)%>%
  #mutate(s_Abun = as.numeric(scale(log(Abun))))%>%
  mutate(s_Abun = (log(Abun) - mu_baseline)/sd_baseline)%>%
  pivot_wider(id_cols = c(Analyte), names_from = c(Subject_ID,Time) ,values_from = s_Abun)

ht_mat <- ht_mat[,c("Analyte",col_ann_dt$colID),]

# metadata on pvalue:
meta_dt <-  result_lme_dt 

# Make a heatmap:
mat <-  ht_mat[,2:ncol(ht_mat)] %>% as.matrix
rownames(mat) <- ht_mat$Analyte

split_cols <-  factor(col_ann_dt$Time, levels = c("Baseline","D1","D14","D30"))

colnames(mat) <- gsub("_.*","",colnames(mat))
#split_rows <- factor(meta_dt$Time, levels = c("Baseline","D1","D14","D30"))

# Row annotation:
row_ann_dt <-  meta_dt %>% data.frame()
names(row_ann_dt)[c(1,5)] <- c("estimate","p.value")
row_ann_dt$Time <- gsub(":.*","",row_ann_dt$Var)
row_ann_dt$Time <- gsub("Time","",row_ann_dt$Time)
row_ann_dt <- row_ann_dt %>%
  select(Time,estimate,Analyte,p.value,p_adj)%>%
  pivot_wider(id_cols = Analyte,names_from = Time,values_from = c(estimate,p.value,p_adj))

row_ann_dt <- row_ann_dt[match(ht_mat$Analyte, row_ann_dt$Analyte),]   # align annotation to heatmap rows


pval_D1 <- row_ann_dt$p.value_D1
pval_D14 <- row_ann_dt$p.value_D14
pval_D30 <- row_ann_dt$p.value_D30
padj_D1 <- row_ann_dt$p_adj_D1
padj_D14 <- row_ann_dt$p_adj_D14
padj_D30 <- row_ann_dt$p_adj_D30





#padj <-  meta_dt$p_adj
ha_right = rowAnnotation(pval_D1 = anno_numeric(pval_D1,
                                                rg = c(min(pval_D1), 1),
                                                x_convert = function(x) ifelse(pval_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", pval_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         pval_D14 = anno_numeric(pval_D14,
                                                 rg = c(min(pval_D14), 1),
                                                 x_convert = function(x) ifelse(pval_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         pval_D30 = anno_numeric(pval_D30,
                                                 rg = c(min(pval_D30), 1),
                                                 x_convert = function(x) ifelse(pval_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", pval_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D1 = anno_numeric(padj_D1,
                                                rg = c(min(padj_D1), 1),
                                                x_convert = function(x) ifelse(padj_D1<= .05,1,0),
                                                labels_format = function(x) sprintf("%.2e", padj_D1),
                                                labels_gp = gpar(fontsize = 3)),
                         padj_D14 = anno_numeric(padj_D14,
                                                 rg = c(min(padj_D14), 1),
                                                 x_convert = function(x) ifelse(padj_D14<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D14),
                                                 labels_gp = gpar(fontsize = 3)),
                         padj_D30 = anno_numeric(padj_D30,
                                                 rg = c(min(padj_D30), 1),
                                                 x_convert = function(x) ifelse(padj_D30<= .05,1,0),
                                                 labels_format = function(x) sprintf("%.2e", padj_D30),
                                                 labels_gp = gpar(fontsize = 3)),
                         
                         annotation_name_rot = 45,
                         annotation_name_gp =gpar(fontsize = 5),
                         width = unit(4, "cm"))

ha_col <- HeatmapAnnotation(Eng = col_ann_dt$Engraftment,col = list(Eng = col_eng))

col_coef <- colorRamp2(c(min(row_ann_dt$estimate_D1), 0, max(row_ann_dt$estimate_D1)), c(muted("red"),"white",muted("blue")))
ha_left =  rowAnnotation(Estimate_D1 = anno_simple(row_ann_dt$estimate_D1,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D1< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D14 =anno_simple(row_ann_dt$estimate_D14,
                                                   col = col_coef,
                                                   pch = ifelse(row_ann_dt$p_adj_D14< 0.05,"*",NA),
                                                   pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         Estimate_D30 = anno_simple(row_ann_dt$estimate_D30,
                                                    col = col_coef,
                                                    pch = ifelse(row_ann_dt$p_adj_D30< 0.05,"*",NA),
                                                    pt_size = unit(rep(5,nrow(meta_dt)), "mm")),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(F,F,F),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         width = unit(1, "cm"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 8,fontface = "bold" ),
                                                        legend_height = unit(0.5, "cm")))

lgd_pvalue = Legend(title = "Coeff", col_fun = col_coef)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")
max_min <-  max(abs(min(meta_dt$Value)),abs(max(meta_dt$Value)))
max_min <- round(max_min,digits = 1)
z_score_min <- round(min(mat,na.rm = T))
z_score_max <-  round(max(mat,na.rm = T))
col_mat <-  colorRamp2(c(z_score_min, 0, z_score_max),c("#0a75ad","white","#ff5353"))



ht1  =  Heatmap(mat,name = "Z-score",
                column_split = split_cols,
                #row_split =  split_rows,
                top_annotation = ha_col,
                #right_annotation =  ha_right,
                right_annotation = ha_left,
                col = col_mat,
                row_names_side = "left",
                row_gap = unit(2, "mm"),
                row_title_gp = gpar(fontsize = 14,fontface = "bold"),
                row_title_rot = 0,
                column_title_gp = gpar(fontsize = 14,fontface = "bold"),
                column_names_gp = gpar(fontsize = 8,fontface = "bold"),
                row_names_gp = gpar(fontsize = 10),
                column_title_rot = 0,
                cluster_rows = T,
                cluster_row_slices = F,
                show_parent_dend_line = F,
                show_row_dend = F,
                cluster_columns = F,
                cluster_column_slices = F,
                show_column_dend = F,
                show_column_names = F,
                border = T,heatmap_legend_param = list(legend_width = unit(5, "cm"),
                                                       legend_direction = "horizontal"),
                row_names_max_width = max_text_width(rownames(mat),
                                                     gp = gpar(fontsize = 15)),
                column_names_max_height = max_text_width(colnames(mat),
                                                         gp = gpar(fontsize = 10)))
pdf(paste(results_folder,paste0("Heatmap_SCFA_HE_LE_Interaction.pdf"),sep = "/"), height = 3, width = 8,useDingbats = F)
draw(ht1,heatmap_legend_side = c("top"),annotation_legend_side="right",
     annotation_legend_list = list(lgd_pvalue, lgd_sig),
     legend_grouping = "original")
dev.off()

