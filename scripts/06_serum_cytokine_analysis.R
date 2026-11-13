# =============================================================================
# 06_serum_cytokine_analysis.R
# -----------------------------------------------------------------------------
# Serum cytokines in melanoma patients. All models use the per-subject log2
# fold change from baseline (Abundance_frac), which also removes plate effects.
#
# Statistics (per cytokine, nlme::lme with a random intercept per subject):
#   Section 1: FC ~ Time                     all visits, each vs Baseline
#   Section 2: FC ~ Time * Response_Wk8      Baseline, D1, D14, D30 (SD vs PD)
#   Section 3: FC ~ Time * Engraftment       Baseline, D1, D14, D30 (HC vs LC)
#   emmeans contrasts; BH adjustment across cytokines.
#
# NOTE: engraftment is labelled HC / LC (high / low colonization) in this
# script, versus HE / LE in the other scripts. Same groups.
#
# Input:  Processed_data/serum_cytokines_melanoma.csv
# Output: Results/06_serum_cytokine_analysis/   (see README.md)
#   Line_plot_Cytokines_VE800_LogAbun.pdf, Line_plot_Cytokines_VE800_FC.pdf,
#   Tab_LME_Serum_Cytokines_Time_FC.csv, Sig_Serum_Cytokines_overtime_FC_boxplot.pdf
#   SD_PD/  timeline, LME tables, Volcano_Cytokines_SD_PD.pdf, LME_SD_PD_INT_Cytokines.pdf
#   HE_LE/  timeline, LME tables, Volcano_Cytokines_HC_LC.pdf, LME_HE_LE_INT_Cytokines.pdf
# =============================================================================

source(here::here("R/utils.R"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(zoo)
  library(lemon)
  library(ggprism)
  library(nlme)
  library(rstatix)
  library(emmeans)
  library(ggplot2)
  library(ggrepel)
})

# Output folder
results_folder <- create_results_dir("Results/06_serum_cytokine_analysis")

# =============================================================================
# 0. Data
# =============================================================================
cyt_dt_final <-  read.csv(here("Processed_data/serum_cytokines_melanoma.csv"))
cyt_dt_final$Engraftment <- factor(cyt_dt_final$Engraftment)
levels(cyt_dt_final$Engraftment) <- c("HC","LC")

# ---- QC: every cytokine per plate (log10 abundance, red line = LLOQ) --------
pdf(paste(results_folder,paste0("Line_plot_Cytokines_VE800_LogAbun",".pdf"),sep="/"),height =4, width = 6)

for (cyt in unique(cyt_dt_final$Cytokines)){
  
  p_cyto <- ggplot(cyt_dt_final[cyt_dt_final$Cytokines == cyt,], aes(x = Time,y = log10(Abundance), group = Subject_ID))+
    geom_point(size = 2,alpha = 0.4)+
    geom_line(alpha = 0.3)+
    facet_rep_grid(~Plate,scales = "free",repeat.tick.labels = T)+
    #facet_grid(Treatment ~ ., scales = "free_y", space = "free_y")+
    #scale_x_continuous(breaks = c( seq(0, 260, by=10)))+
    theme_prism()+
    geom_hline(aes(yintercept = log10(LLOQ)),linetype="dashed",color = "red")+
    #ylab("Subjects")+
    theme(axis.text.x=element_text(size=10,angle = 45, hjust=1,vjust = 1,face="bold"),
          axis.text.y=element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"))+ggtitle(cyt)
  
  
  
  print(p_cyto)
}
dev.off()



# ---- Every cytokine per plate as log2 fold change vs baseline ---------------
pdf(paste(results_folder,paste0("Line_plot_Cytokines_VE800_FC",".pdf"),sep="/"),height =3, width = 8)
for (cyt in unique(cyt_dt_final$Cytokines)){
  
  p_cyto <- ggplot(cyt_dt_final[cyt_dt_final$Cytokines == cyt,], aes(x = Time,y = Abundance_frac, group = Subject_ID))+
    geom_point(size = 2,alpha = 0.4)+
    geom_line(alpha = 0.3)+
    facet_rep_grid(~Plate,scales = "free",repeat.tick.labels = T)+
    #facet_grid(Treatment ~ ., scales = "free_y", space = "free_y")+
    #scale_x_continuous(breaks = c( seq(0, 260, by=10)))+
    theme_prism()+
    #geom_hline(aes(yintercept = log10(LLOQ)),linetype="dashed",color = "red")+
    #ylab("Subjects")+
    ylab("Log2(Abundance/Baseline)")+
    theme(axis.text.x=element_text(size=10,angle = 45, hjust=1,vjust = 1,face="bold"),
          axis.text.y=element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"))+ggtitle(cyt)
  print(p_cyto)
  
}
dev.off()



# =============================================================================
# 1. Effect of time (all subjects): log2 FC ~ Time, each visit vs Baseline
# =============================================================================

lme_dt <-  cyt_dt_final
result_lme <- list()
for (scfa_var in unique(lme_dt$Cytokines)){
  sel_dt <-  lme_dt[lme_dt$Cytokines %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- sel_dt$Abundance_frac
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D21","D30","D58","D86"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  mod_bac <- lme(t_Abundance ~  Time  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Cytokines <- scfa_var
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
write.csv(result_lme_dt,paste(results_folder,paste0("Tab_LME_Serum_Cytokines_Time_FC.csv"),sep = "/"))

result_lme_dt$Time <-  gsub("Time","",result_lme_dt$Var)
result_lme_dt$group2 <- result_lme_dt$Time

# Filter results with pval
p_cutoff <-  0.1
result_sig <- result_lme_dt[result_lme_dt$p_adj < p_cutoff ,]

if(nrow(result_sig)>0){
  sel_dt <-  lme_dt %>% filter(Cytokines %in% unique(result_sig$Cytokines))
  sel_dt$L2FC  <- sel_dt$Abundance_frac
  # Make pval dataframe for display and use the LME results pval
  # Add pvalue
  df_p_val <- sel_dt %>%
    group_by(Cytokines)%>%
    wilcox_test(L2FC ~ Time,ref.group = "Baseline") %>%
    add_xy_position(fun = "mean_sd") %>% 
    select(!c(p))%>% 
    inner_join(result_sig)%>%
    mutate(p_label = sprintf("%.4f",p_adj) )%>%
    add_significance(p.col = "p_adj")%>%
    mutate(Engraftment = "")
  
  # Make a boxplot:
  p_sig <- ggplot(sel_dt, aes(x = Time, y = L2FC)) +
    geom_boxplot( alpha = 0.3, outlier.colour = NA)+
    geom_point(aes(y = L2FC),
               alpha = 0.5,size = 2)+
    facet_rep_wrap(~Cytokines , scales = "free",ncol = 3)+
    theme_prism()+
    theme(axis.text.x=element_text(size=14,angle = 45, hjust=1,vjust = 1,face="bold"),
          axis.text.y=element_text(size=14,face="bold"),
          axis.title.y = element_text(size=14,face="bold"),
          axis.title.x = element_text(size=14,face="bold"))+
    xlab("Day")+
    ylab("Log2(Abundance/Baseline)")+
    add_pvalue(df_p_val,label = "padj={p_label}",tip.length = 0.01,)
  
  pdf(paste(results_folder,paste0("Sig_Serum_Cytokines_overtime_FC_boxplot.pdf"),sep="/"),height =20, width = 12)
  print(p_sig)
  dev.off()
  
} 


# =============================================================================
# 2. SD vs PD (Week-8 response)
# =============================================================================
results_folder <- create_results_dir("Results/06_serum_cytokine_analysis/SD_PD")

# ---- Sampling timeline (shape = assay plate) ---------------------------------

col_resp <- c("#ff5e1e","#0099b0")
names(col_resp) <- c("PD", "SD")

timeline_dt <-  cyt_dt_final %>% select(Subject_ID,Time,Response_Wk8,Cohort,Engraftment,Plate)%>% unique()
timeline_dt$Subject_ID <- as.character(timeline_dt$Subject_ID)

timeline_dt$Timepoint <-  factor(timeline_dt$Time, levels = c("Baseline","D1","D14","D30","D58","D86"))
#levels(timeline_dt$VISIT_NAME) <-  c("Baseline","Wk_2","Wk_4","Wk_8")

p_time <- ggplot(timeline_dt, aes(x = Timepoint,y = Subject_ID,
                                  group = Subject_ID,color = Response_Wk8))+
  geom_point(size = 3,alpha = 0.7,aes(shape = Plate))+
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

pdf(paste(results_folder,paste0("Timeline_Serum_Cytokines_SD_PD",".pdf"),sep="/"),height =4, width = 7)
print(p_time)
dev.off()


# ---- 2a. LME: Time x Response, emmeans contrasts -----------------------------
lme_dt <-  cyt_dt_final %>% 
  filter(Time %in% c("Baseline","D1","D14","D30"))
result_lme <- list()
effect_time_lme <-  list()
group_time_lme <-  list()
for (scfa_var in unique(lme_dt$Cytokines)){
  sel_dt <-  lme_dt[lme_dt$Cytokines %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- sel_dt$Abundance_frac
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D30"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  sel_dt$Response_Wk8 <- factor(sel_dt$Response_Wk8, levels = c("SD","PD"))
  mod_bac <- lme(t_Abundance ~  Time*Response_Wk8  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Cytokines <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
  
  # Time effect compared to baseline in each group
  time_emm <- emmeans(mod_bac, ~ Time | Response_Wk8)
  eff_dt <-  data.frame(pairs(time_emm))
  eff_dt <-  eff_dt[grep("Baseline",eff_dt$contrast),]
  eff_dt$Cytokines <- scfa_var
  effect_time_lme[[scfa_var]] <- eff_dt
  
  # Cross-sectional comparison across each time point between two groups:
  group_emm <- emmeans(mod_bac, ~ Response_Wk8 | Time)
  group_eff_dt <-  data.frame(pairs(group_emm))
  group_eff_dt$Cytokines <- scfa_var
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

write.csv(result_lme_dt,paste(results_folder,"Tab_LME_Int_SD_PD.csv",sep = "/"))
write.csv(group_comp_dt,paste(results_folder,"Tab_LME_SD_Vs_PD_Time.csv",sep = "/"))
write.csv(time_comp_dt,paste(results_folder,"Tab_LME_Time_Effect_SD_PD.csv",sep = "/"))


# ---- 2b. Volcano: SD vs PD at each post-baseline visit (unadjusted p) --------

vol_dt <-  group_comp_dt %>% filter(!Time %in% c("Baseline"))


# Volcano plots:
res_dt <- vol_dt
names(res_dt)[7] <- "pval"
# NOTE: the response is already a log2 fold change, so `estimate` is on the log2
res_dt$L2FC <-  res_dt$estimate
res_dt$dir <-  ifelse(res_dt$L2FC > 0, "SD","PD")

# Not recommended but we are only looking for trends
q_cut_off <- 0.05  
data_pvalue_less <- res_dt[which(res_dt$pval< q_cut_off),]

res_dt$dir[res_dt$pval > q_cut_off] <- "NS"

res_dt$alpha <-  ifelse(res_dt$pval <= q_cut_off,0.8,0.1)

col_vol <-  c(col_resp,"NS" = "black")

res_dt$y <- -log10(res_dt$pval)
res_dt$label <- res_dt$Cytokines
res_dt$label[res_dt$pval > q_cut_off] <- "" 

p_vol  <- ggplot(data=res_dt,
                 aes(y=y,x=L2FC, fill = dir, alpha = alpha)) + 
  geom_point(shape = 21,size = 3,color="black")+
  geom_text_repel(aes(label = label))+
  scale_fill_manual(name="Response",values = col_vol) +
  geom_hline(yintercept = -log10(q_cut_off),linewidth = .1,linetype = "dashed")+
  geom_vline(xintercept = 0,linewidth = .1,linetype = "dashed")+
  theme_prism(border = TRUE)+
  facet_wrap(~Time,nrow = 1)+
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 7)))+
  theme(legend.position =  "none")+
  xlab("Log2 Fold Change (SD/PD)")+
  ylab("-Log10 pval")
# +#
#   xlim(-1,10)+
#   ylim(0,20)
#  geom_magnify(from = from,to = to,recompute = T)



pdf(paste0(results_folder,"/Volcano_Cytokines_SD_PD.pdf"), height = 4,width = 10,useDingbats = F)
print(p_vol)
dev.off()



# ---- 2c. Box plots: cytokines with a Time x Response interaction (p < 0.05) --
write.csv(result_lme_dt,paste(results_folder,"Tab_LME_Cytokines_Int_SD_PD.csv",sep = "/"))
result_lme_dt$Time <-  gsub("Time","",result_lme_dt$Var)
result_lme_dt$Time <-  gsub(":.*","",result_lme_dt$Time)
result_lme_dt$group2 <- result_lme_dt$Time
# Filter results with pval
p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$pval < p_cutoff ,]
#  unique(result_sig$Bac)

if(nrow(result_sig)>0){
  sel_dt <-  lme_dt %>% filter(Cytokines %in% unique(result_sig$Cytokines))
  #sel_dt$LogAbun  <- log10(sel_dt$Abun)
  # Make pval dataframe for display and use the LME results pval
  # Add pvalue
  df_p_val <- sel_dt %>%
    group_by(Cytokines)%>%
    wilcox_test(Abundance_frac ~ Time,ref.group = "Baseline") %>%
    add_xy_position(fun = "mean_sd") %>% 
    select(!c(p))%>% 
    inner_join(result_sig)%>%
    mutate(p_adj_label = sprintf("%.4f",p_adj) )%>%
    mutate(p_unadj_label = sprintf("%.4f",pval) )%>%
    add_significance(p.col = "p_adj")%>%
    mutate(Response_Wk8 = "")
  
  # Make a boxplot:
  p_sig <- ggplot(sel_dt, aes(x = Time,y =Abundance_frac ))+
    geom_boxplot(outlier.colour = NA,aes(fill = Response_Wk8),alpha = 0.2)+ 
    #scale_color_manual(values = col_resp)+
    scale_fill_manual(values = col_resp)+
    # new_scale_color()+  
    geom_point(position=position_dodge(width=0.75),aes(group = Response_Wk8,fill = Response_Wk8),shape = 21,color = "black",size = 2,alpha = 0.5)+
    facet_wrap(~Cytokines, scales = "free", ncol = 3)+
    theme_bw()+
    #geom_hline(aes(yintercept = 0),color = "red",linetype="dashed")+
    #ylab("Subjects")+
    #ylab("FC = Log2(Abundance/Baseline)")+
    theme_prism()+
    theme(axis.text.x=element_text(size=10,angle = 45, hjust=1,vjust = 1,face="bold"),
          axis.text.y=element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"))+
    theme(legend.position="right")+
    add_pvalue(df_p_val,label = "p_unadj = {p_unadj_label}",tip.length = 0.01,)
  pdf(paste(results_folder,paste0("LME_SD_PD_INT_Cytokines.pdf"),sep="/"),height =8, width = 10)
  print(p_sig)
  dev.off()
} 


# =============================================================================
# 3. HC vs LC (VE800 engraftment) - same analyses as section 2
# =============================================================================

results_folder <- create_results_dir("Results/06_serum_cytokine_analysis/HE_LE")


# ---- Sampling timeline coloured by engraftment ------------------------------
col_eng <- c("#009031","#7f1f98")
names(col_eng) <- c("HC", "LC")

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

pdf(paste(results_folder,paste0("Timeline_Serum_Cytokines_HE_LE",".pdf"),sep="/"),height =4, width = 6)
print(p_time)
dev.off()


# ---- 3a. LME: Time x Engraftment, emmeans contrasts --------------------------
lme_dt <-  cyt_dt_final %>% 
  filter(Time %in% c("Baseline","D1","D14","D30"))
result_lme <- list()
effect_time_lme <-  list()
group_time_lme <-  list()
for (scfa_var in unique(lme_dt$Cytokines)){
  sel_dt <-  lme_dt[lme_dt$Cytokines %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- sel_dt$Abundance_frac
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D30"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  sel_dt$Response_Wk8 <- factor(sel_dt$Response_Wk8, levels = c("SD","PD"))
  #sel_dt$Engraftment <- factor(sel_dt$Engraftment,levels = c("LE","HE"))
  mod_bac <- lme(t_Abundance ~  Time*Engraftment  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Cytokines <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
  
  # Time effect compared to baseline in each group
  time_emm <- emmeans(mod_bac, ~ Time | Engraftment)
  eff_dt <-  data.frame(pairs(time_emm))
  eff_dt <-  eff_dt[grep("Baseline",eff_dt$contrast),]
  eff_dt$Cytokines <- scfa_var
  effect_time_lme[[scfa_var]] <- eff_dt
  
  # Cross-sectional comparison across each time point between two groups:
  group_emm <- emmeans(mod_bac,~ Engraftment | Time)
  group_eff_dt <-  data.frame(pairs(group_emm))
  group_eff_dt$Cytokines <- scfa_var
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

write.csv(result_lme_dt,paste(results_folder,"Tab_LME_Int_HE_LE.csv",sep = "/"))
write.csv(group_comp_dt,paste(results_folder,"Tab_LME_HE_Vs_LE_Time.csv",sep = "/"))
write.csv(time_comp_dt,paste(results_folder,"Tab_LME_Time_Effect_HE_LE.csv",sep = "/"))


# ---- 3b. Volcano: HC vs LC at each post-baseline visit (unadjusted p) --------

vol_dt <-  group_comp_dt %>% filter(!Time %in% c("Baseline"))


# Volcano plots:
res_dt <- vol_dt
names(res_dt)[7] <- "pval"
# NOTE: the response is already a log2 fold change, so `estimate` is on the log2
res_dt$L2FC <-  res_dt$estimate
res_dt$dir <-  ifelse(res_dt$L2FC > 0, "HC","LC")

# Not recommended but we are only looking for trends
q_cut_off <- 0.05  
data_pvalue_less <- res_dt[which(res_dt$pval< q_cut_off),]

res_dt$dir[res_dt$pval > q_cut_off] <- "NS"

res_dt$alpha <-  ifelse(res_dt$pval <= q_cut_off,0.8,0.1)

col_vol <-  c(col_eng,"NS" = "black")

res_dt$y <- -log10(res_dt$pval)
res_dt$label <- res_dt$Cytokines
res_dt$label[res_dt$pval > q_cut_off] <- "" 

p_vol  <- ggplot(data=res_dt,
                 aes(y=y,x=L2FC, fill = dir, alpha = alpha)) + 
  geom_point(shape = 21,size = 3,color="black")+
  geom_text_repel(aes(label = label))+
  scale_fill_manual(name="Engraftment",values = col_vol) +
  geom_hline(yintercept = -log10(q_cut_off),linewidth = .1,linetype = "dashed")+
  geom_vline(xintercept = 0,linewidth = .1,linetype = "dashed")+
  theme_prism(border = TRUE)+
  facet_wrap(~Time,nrow = 1)+
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 7)))+
  theme(legend.position =  "none")+
  xlab("Log2 Fold Change (HC/LC)")+
  ylab("-Log10 pval")
# +#
#   xlim(-1,10)+
#   ylim(0,20)
#  geom_magnify(from = from,to = to,recompute = T)

pdf(paste0(results_folder,"/Volcano_Cytokines_HC_LC.pdf"), height = 4,width = 10,useDingbats = F)
print(p_vol)
dev.off()



# ---- 3c. Interaction model with LC as reference + box plots (p < 0.05) -------
lme_dt <-  cyt_dt_final %>% 
  filter(Time %in% c("Baseline","D1","D14","D30"))

result_lme <- list()
for (scfa_var in unique(lme_dt$Cytokines)){
  sel_dt <-  lme_dt[lme_dt$Cytokines %in% scfa_var,]
  sel_dt$PID <- factor(sel_dt$Subject_ID)
  sel_dt$t_Abundance  <- sel_dt$Abundance_frac
  sel_dt$Time <- factor(sel_dt$Time , levels = c("Baseline","D1","D14","D30"))
  sel_dt$Time <-  relevel(sel_dt$Time,ref = "Baseline")
  sel_dt$Engraftment <- factor(sel_dt$Engraftment, levels = c("LC","HC"))
  mod_bac <- lme(t_Abundance ~  Time*Engraftment  ,random = ~1 |PID, data = sel_dt)
  # Optimize if it doesnt work
  sum_mod <-  summary(mod_bac)
  sum_mod_dt <- data.frame(sum_mod$tTable)
  sum_mod_dt$Cytokines <- scfa_var
  sum_mod_dt$Var <-  rownames(sum_mod_dt)
  result_lme[[scfa_var]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
# Keep only interaction
result_lme_dt <- result_lme_dt[grep(":",result_lme_dt$Var),]

# Group by adjust
result_lme_dt <-  result_lme_dt %>%
  group_by(Var)%>%
  mutate(p_adj = p.adjust(pval, method = "BH"))

write.csv(result_lme_dt,paste(results_folder,"Tab_LME_Cytokines_Int_HE_LE.csv",sep = "/"))
result_lme_dt$Time <-  gsub("Time","",result_lme_dt$Var)
result_lme_dt$Time <-  gsub(":.*","",result_lme_dt$Time)
result_lme_dt$group2 <- result_lme_dt$Time
# Filter results with pval
p_cutoff <-  0.05
result_sig <- result_lme_dt[result_lme_dt$pval < p_cutoff ,]
#  unique(result_sig$Bac)

if(nrow(result_sig)>0){
  sel_dt <-  lme_dt %>% filter(Cytokines %in% unique(result_sig$Cytokines))
  #sel_dt$LogAbun  <- log10(sel_dt$Abun)
  # Make pval dataframe for display and use the LME results pval
  # Add pvalue
  df_p_val <- sel_dt %>%
    group_by(Cytokines)%>%
    wilcox_test(Abundance_frac ~ Time,ref.group = "Baseline") %>%
    add_xy_position(fun = "mean_sd") %>% 
    select(!c(p))%>% 
    inner_join(result_sig)%>%
    mutate(p_adj_label = sprintf("%.4f",p_adj) )%>%
    mutate(p_unadj_label = sprintf("%.4f",pval) )%>%
    add_significance(p.col = "p_adj")%>%
    mutate(Engraftment = "")
  
  # Make a boxplot:
  p_sig <- ggplot(sel_dt, aes(x = Time,y =Abundance_frac ))+
    geom_boxplot(outlier.colour = NA,aes(fill = Engraftment),alpha = 0.2)+ 
    #scale_color_manual(values = col_resp)+
    scale_fill_manual(values = col_eng)+
    # new_scale_color()+  
    geom_point(position=position_dodge(width=0.75),aes(group = Engraftment,fill = Engraftment),shape = 21,color = "black",size = 2,alpha = 0.5)+
    facet_wrap(~Cytokines, scales = "free", ncol = 3)+
    theme_bw()+
    #geom_hline(aes(yintercept = 0),color = "red",linetype="dashed")+
    #ylab("Subjects")+
    #ylab("FC = Log2(Abundance/Baseline)")+
    theme_prism()+
    theme(axis.text.x=element_text(size=10,angle = 45, hjust=1,vjust = 1,face="bold"),
          axis.text.y=element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"))+
    theme(legend.position="right")+
    add_pvalue(df_p_val,label = "p_unadj = {p_unadj_label}",tip.length = 0.01,)
  pdf(paste(results_folder,paste0("LME_HE_LE_INT_Cytokines.pdf"),sep="/"),height =8, width = 10)
  print(p_sig)
  dev.off()
} 
