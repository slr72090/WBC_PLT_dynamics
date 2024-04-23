

rm(list=ls())
require(dplyr)
require(tidyverse)
require(ggplot2)
require(magrittr)
require(lubridate)
require(MASS) 
require(reshape2)
library(ggpubr)
require(scales)
select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

library(ggplot2)
library(cowplot)
library(ggcorrplot)
library(gslnls)

## Specify files and directory --------------------------------
working.directory <- "/Users/sranjeva/Desktop/Research/Projects/Higgins_lab/R_code/experiments_3_29_2024"
plotting.directory <- "/Users/sranjeva/Desktop/Research/Projects/Higgins_lab/R_code/experiments_3_29_2024/plots"
setwd(working.directory)

# Specify whether or not to save plots --------------------------------
save_plots = F
fit_recovery_models = F

textSize = 11
source("./Plotting_scripts/plot_themes.R")
time_bounds = c(0,5) # time bounds for recovery models 
escalation_start = -4 # Days prior to WBC peeak for start of escalation phasee 
GOF_threshold = 0.8 # Specify threshold for R2
delta_plt_threshold = 20 #Exclude patients with absolute platelet change under this threshold 
n_min_observations = 3 # minimum number of pre-WBC-peak observations
## Load in and format data ## ----------------------------------------------------------
source("./Data_formatting_scripts/format_data_combined.R")

## Choose patients for analysis and format data ## --------------------------------------------------
fit_recovery_models = F
sub_population = "survived" # patients that survived and had good model fit 
source("data_from_model_results.R")


# Plot mean phase plane trajectories ---------------------------------------------------------------------
this_title = "WBC-PLT phase plane for escalation phase"
p_traj <- df_phase_plane %>% 
  filter(time > -5) %>%  
  arrange(.,time) %>%
  ggplot(aes(x = mean_delta_WBC, y = mean_delta_PLT)) + 
  #ggplot(aes(x = mean_delta_w, y = mean_delta_p)) + 
  geom_point(aes(color = as.factor(Diagnosis))) +
  geom_path(aes(color = as.factor(Diagnosis)),
            arrow = arrow(type = 'open', 
                          angle = 30, 
                          length = unit(0.5, "cm"))) + 
  geom_errorbar(aes(ymin = (mean_delta_PLT - se_delta_PLT/2),
                    ymax = (mean_delta_PLT  + se_delta_PLT/2), 
                    color = as.factor(Diagnosis)), linetype = 2, width = 0.2) + 
  #facet_wrap(~Diagnosis, scales = "free") +
  plot_themes + 
  ggtitle(this_title)+
  labs(color = "Diagnosis") + 
  xlab("WBC change from first observation") + ylab("PLT change from first observation")


df_stats_phase_plane <- df_pre_combined %>% select(MRN, time, p_w_ratio, Diagnosis) %>% 
  left_join(.,df_phase_plane %>% select(time, LQ_p_w, UQ_p_w, n,Diagnosis)) %>% 
  filter(is.finite(p_w_ratio)) %>% 
  group_by(time, Diagnosis) %>% 
  summarise(frac = sum(p_w_ratio > LQ_p_w& p_w_ratio < UQ_p_w)/n(),
            n = n()
            )

p_stats <- ggplot(df_stats_phase_plane, aes(x = time, y = frac)) + 
  geom_bar(aes(fill = as.factor(Diagnosis)), stat = "identity", position = "dodge") +
  plot_themes + 
  labs(fill = "N observed") + 
  xlab("Time prior to WBC peak") +
  ylab("Fraction of patients within IQR of median scaled PLT:WBC") + 
  ylim(0,1)


# PLT dynamics  ## ------------------------------------------

# Trends in PLT before inferred lag point 
PLT_class <- df_pre_best_PLT %>% filter(t_rel_plt_lag == -1) %>% 
  mutate(PLT_classification = as.numeric(p_scaled <0)) %>% 
  select(MRN, Diagnosis, PLT_classification) %>% 
  left_join(.,dat_demog %>% select(MRN, Diagnosis, Mortality)) %>% 
  as.data.frame()

p_pre_PLT <-ggplot(df_sum_PLT %>%  filter(t_rel_plt_lag < 2 & t_rel_plt_lag > -4 ), aes(x = t_rel_plt_lag, y = median)) +
  geom_line(data = df_pre_best_PLT %>% left_join(.,PLT_class) %>% filter(t_rel_plt_lag < 2 & t_rel_plt_lag > -4 ) , aes( x = t_rel_plt_lag, y = p_scaled, group = as.factor(MRN)), linetype = 2, alpha = 0.05) + 
  geom_line(data = df_sum_PLT %>% filter(t_rel_plt_lag < 2 & t_rel_plt_lag > -4 ),  aes(x = t_rel_plt_lag, y = median, color = Diagnosis)) + 
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = Diagnosis), alpha = 0.2) + 
  xlim(-3,1) + 
  plot_themes + 
  ylim(-1,2) + 
  xlab("Time relative to inferred PLT lag (following WBC peak)") + 
  ylab("scaled PLT") + 
  facet_wrap(~Diagnosis)


PLT_sum <- df_pre_best_PLT %>% 
  left_join(.,PLT_class) %>% 
  group_by(Diagnosis, t_rel_plt_lag, PLT_classification) %>% 
  summarise(mean = mean(p_scaled, na.rm = T),
            median = median(p_scaled, na.rm =T),
            LCI = quantile(p_scaled, c(0.25,0.75))[1],
            UCI = quantile(p_scaled, c(0.25,0.75))[2])

p_pre_PLT_class <- ggplot(data = df_pre_best_PLT %>% 
                  left_join(.,PLT_class) %>% 
                  filter(p_scaled > -5 & p_scaled < 5 & t_rel_plt_lag < 2 & t_rel_plt_lag > -4 ) , 
                aes( x = as.factor(t_rel_plt_lag), y = p_scaled, group = as.factor(MRN))) +
  geom_line(linetype = 2, alpha = 0.3, aes(color = as.factor(PLT_classification))) +
  plot_themes + facet_wrap(~Diagnosis, scales = "free") + ylim(-2,2)

p_pre_PLT_sum <- ggplot(data = PLT_sum %>% filter(t_rel_plt_lag < 2 & t_rel_plt_lag > -4 & !is.na(PLT_classification))) + 
  geom_line(aes(x = t_rel_plt_lag, y = median, color = as.factor(PLT_classification))) + 
  # geom_line(data = df_pre_best_PLT %>% 
  #             left_join(.,PLT_class) %>% 
  #             filter(p_scaled > -5 & p_scaled < 5 & t_rel_plt_lag < 2 & t_rel_plt_lag > -4 ), 
  #           linetype = 2, alpha = 0.4,
  #           aes( x = t_rel_plt_lag, y = p_scaled, color = as.factor(PLT_classification))) + 
  geom_ribbon(aes(x = t_rel_plt_lag, ymin = LCI, ymax = UCI, fill = as.factor(PLT_classification)), alpha = 0.2) +
  plot_themes + ylim(-1,1) +
  facet_wrap(~Diagnosis) + xlab("time relative to inferred platelet lag"
  ) + ylab("scaled PLT") +
  labs(fill = "pre-lag PLT classification", color = "pre-lag PLT classification")
  

p_combined = plot_grid(p_pre_WBC, p_pre_PLT, ncol = 2)


# PLT stability prior to WBC peak ##----------
source("data_for_risk_analysis.R")
my_comparisons = list(c("cdiff", "MI"), 
                      c("cdiff", "sepsis"), 
                      c("cdiff", "stroke"),
                      c("MI", "sepsis"),
                      c("MI", "stroke"),
                      c("sepsis","stroke") 
)


df_sub_COV = df_sub %>% filter(MRN %in% these_MRNs)

p_stability <- ggplot(df_sub_COV, aes(x = Diagnosis, y = CoV_PLT)) + 
  geom_boxplot(aes(fill = as.factor(Diagnosis)), alpha = 0.2) + 
  geom_jitter(aes(color = as.factor(Diagnosis)), alpha = 0.05) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("") + ylab("PLT coefficient of variation") + 
  plot_themes + labs(color = "", fill = "")

my_comparisons = as.list(all_comparisons)

p_mortality_COV <- ggplot(df_mortality_COV %>% filter(time == 0), aes(x = injury_type, y = rel_mort, fill = as.factor(q_binned))) + 
  geom_bar(aes(fill = as.factor(q_binned)), alpha = 0.2, stat= "identity" , position = "dodge") + 
  plot_themes +  
  stat_compare_means(comparisons = my_comparisons, aes(label = paste0("p = ",scales::label_pvalue()(..p..)))) +
  xlab("Day pre-WBC peak") + 
  ylab("Relative mortality compared to mean") + 
  labs(fill = "Quantile of PLT CoV prior to WBC peak") 














# Make individual plots ## ------------------------------------------------
# Look at individuals with >= 3 observations before peak
ind <- df_pre_best_WBC %>% filter(time > escalation_start) %>% 
  group_by(MRN, Diagnosis, Mortality) %>% 
  summarise(n = n()) %>% 
  filter(n >= n_min_observations & Mortality == 1)

for(i in c(1:nrow(ind))){
  this_MRN = ind[i,]$MRN
p_pre_WBC_ind <- ggplot() +
  geom_point(data = df_pre_best_WBC %>% filter(MRN == this_MRN), aes( x = time, y = w_scaled, group = as.factor(MRN))) + 
  geom_line(data = df_pre_best_WBC %>% filter(MRN == this_MRN), aes( x = time, y = w_scaled, group = as.factor(MRN)), linetype = 2, alpha = 0.5) + 
 # geom_line(data = df_sum_WBC,  aes(x = time, y = median, color = Diagnosis)) + 
  #geom_ribbon(data = df_sum_WBC,aes(x = time, ymin = LCI, ymax = UCI, fill = Diagnosis), alpha = 0.2) + 
  xlim(-4,0) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  ylab("WBC") +
  ylim(0,1)

p_pre_PLT_ind <- ggplot() +
  #geom_line(data = df_sum_PLT,  aes(x = time, y = median, color = Diagnosis)) + 
  geom_point(data = df_pre_best_PLT %>% filter(MRN == this_MRN), aes(x = time, y = p_scaled)) + 
  geom_line(data = df_pre_best_PLT %>% filter(MRN == this_MRN), aes(x = time, y = p_scaled), linetype =2, alpha = 0.5) + 
  #geom_ribbon(data = df_sum_PLT, aes(x = time, ymin = LCI, ymax = UCI, fill = Diagnosis), alpha = 0.2) + 
  xlim(-4,2) + 
  plot_themes + 
  xlab("Time relative to WBC peak") + 
  ylab("PLT")

p_combined = plot_grid(p_pre_WBC_ind, p_pre_PLT_ind, ncol = 2)
save_plot(paste0("./plots_mortality/plots_escalation_", this_MRN, ".pdf"), p_combined)
}