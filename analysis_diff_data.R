

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
working.directory <- "/Users/sranjeva/Desktop/Research/Projects/Higgins_lab/R_code/experiments_4_15_2024"
plotting.directory <- "/Users/sranjeva/Desktop/Research/Projects/Higgins_lab/R_code/experiments_4_15_2024/plots"
setwd(working.directory)
textSize = 11
save_plots = T
time_bounds = c(-4:1)
source("./Plotting_scripts/plot_themes.R")

dat_diff <- read.csv("./Diff_data/tDiffsFull.2024.04.15.csv") 

dat_MRN <- dat_diff %>% 
  select(EPIC_PMRN, MRN, MRN_Type) %>% 
  filter(MRN_Type == "MGH") %>% 
  unique()

data_lab<- dat_diff %>% 
  select(EPIC_PMRN, Seq_Date_Time, Group_Id, Test_Id, Result, Reference_Range) %>% 
  mutate(Date = as.Date(Seq_Date_Time, format =  c("%m/%d/%y"))) %>% 
  left_join(.,dat_MRN %>% select(EPIC_PMRN, MRN)) %>% 
  filter(!is.na(MRN))

##

data_WBC_max <- read.csv("max_WBC_all.csv") %>% select(-X) %>% 
  mutate(Date = as.Date(Date),
         Diag_date = as.Date(Diag_date))

data_diag = data_WBC_max %>% select(MRN, Diag_date, Diagnosis)
data <- data_lab %>% 
              select(
                MRN, 
                variable = Group_Id,
                Result,
                Date
              ) %>% 
  left_join(., full_cohort %>% select(MRN, Diag_date, Diagnosis)) %>% 
  filter(as.numeric(Date - Diag_date) %in% c(-10:20)) %>% 
  #filter(grepl(c("IMMGRAN|NEUT|LYM|MON|BAN" ),variable)) %>% 
  filter(variable %in% c("LYMP", "MON", "NEUT","IMMGRAN%", "BAND", "WBC")) %>% 
  group_by(MRN, Diagnosis, Diag_date, variable) %>% 
  mutate(t_rel = as.numeric(Date - min(Date))) %>% 
  as.data.frame() 




data2 <- data %>% group_by(MRN, Diag_date, Diagnosis,Date, variable, t_rel) %>% 
  summarize(result = max(Result,na.rm = T)) %>% 
  filter(is.finite(result)) %>% 
  as.data.frame() %>% 
  mutate(method = 2)

data2_all <- data2 %>% 
  group_by(MRN, Diag_date, Diagnosis, Date, variable, t_rel) %>% 
  filter(method == min(method)) %>% 
  as.data.frame() %>% 
  filter(is.finite(result))

df_max_wbc <- data2_all %>% select(MRN,  Diagnosis, Diag_date) %>% as.data.frame() %>% unique() %>% 
  left_join(.,data_WBC_max %>% select(MRN, Date, Diagnosis, Diag_date, max_wbc = Result)) %>% 
  mutate(t_max_wbc = as.numeric(Date - Diag_date)) %>% 
  select(MRN,max_wbc, t_max_wbc) %>% 
  as.data.frame()

data3 <- data2_all %>% 
  left_join(., df_max_wbc) %>% 
  mutate(t_rel_WBC = t_rel - t_max_wbc)

data_all = data3 %>% 
  filter(!(MRN %in% data_all_old$MRN)) %>% 
  bind_rows(.,data_all_old) %>% 
  left_join(.,dat_demog %>% select(MRN, Diagnosis, Mortality))

df_stats <- data_all %>% 
  group_by(variable, Diagnosis, t_rel_WBC) %>% 
  summarize(n = n())

dat_sum = data_all %>% 
  select(MRN,Diagnosis,Mortality, variable, result, t_rel_WBC) %>% 
  as.data.frame() %>% 
  group_by(Diagnosis,Mortality, variable, t_rel_WBC) %>% 
  summarize(
    n = n(),
    mean = mean(result, na.rm = T),
    median = median(result, na.rm =T),
    LQ = quantile(result,na.rm = T)[2],
    UQ = quantile(result, na.rm =T)[4]
  ) 


p_N <- ggplot(df_stats %>% filter(variable!="WBC"), aes(x = t_rel_WBC, y = n)) +
  geom_line( aes(color = Diagnosis)) + 
  geom_point( aes(color = Diagnosis)) +
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  facet_wrap(~variable, scales = "free_y") + 
  xlim(-4,0)

p <- ggplot(dat_sum %>% filter(t_rel_WBC %in% c(-4:0)), aes(x = t_rel_WBC, y = median)) +
  geom_line( aes(x = t_rel_WBC, y = median, color = Diagnosis)) + 
  geom_point( aes(x = t_rel_WBC, y = median, color = Diagnosis)) +
  #geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = Diagnosis), alpha = 0.2) + 
  geom_text(aes(label = n)) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  facet_wrap(~variable, scales = "free_y") + 
  xlim(-4,0)

p <- ggplot(dat_sum %>% filter(t_rel_WBC %in% c(-4:0) & Diagnosis == "stroke"), aes(x = t_rel_WBC, y = median)) +
  geom_line( aes(x = t_rel_WBC, y = median, color = as.factor(Mortality))) + 
  geom_point( aes(x = t_rel_WBC, y = median, color = as.factor(Mortality))) +
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = as.factor(Mortality)), alpha = 0.2) + 
  geom_text(aes(label = n)) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  facet_wrap(~variable, scales = "free_y") + 
  xlim(-4,0)

my_comparisons = list(c("cdiff", "MI"), 
                      c("cdiff", "sepsis"), 
                      c("cdiff", "stroke"),
                      c("MI", "sepsis"),
                      c("MI", "stroke"),
                      c("sepsis","stroke") 
)

for(i in c(1:length(unique(data_all$variable)))){
this_variable = unique(data_all$variable)[i]

p_var1 <- ggplot(dat_sum   %>% filter(variable == this_variable & Mortality == 1 & t_rel_WBC %in% time_bounds), aes(x = t_rel_WBC, y = median)) +
  geom_line( aes(x = t_rel_WBC, y = median, color = Diagnosis)) + 
  geom_point( aes(x = t_rel_WBC, y = median, color = Diagnosis)) +
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = Diagnosis), alpha = 0.2) + 
  plot_themes + 
  ggtitle(paste0(this_variable)) + 
  xlab("Time relative to  WBC peak") +
  xlim(min(time_bounds), max(time_bounds))

p_var2 <- ggplot(data_all %>% filter(variable == this_variable &  Mortality == 0 &t_rel_WBC %in% time_bounds), aes(x = as.factor(t_rel_WBC), y = log(result), fill = Diagnosis)) +
  geom_boxplot(aes(fill = Diagnosis), alpha = 0.4) +
  geom_point(aes(color = as.factor(Diagnosis)), alpha = 0.2, position=position_jitterdodge()) +
  ylab("Log value") + 
  plot_themes + 
  xlab("Time relative to  WBC peak")  + 
  ggtitle(paste0(this_variable))

p_var3 <- ggplot(data_all %>% filter(variable == this_variable &   Mortality == 1 &t_rel_WBC %in% time_bounds), aes(x = as.factor(Diagnosis), y = log(result), fill = Diagnosis)) +
  geom_boxplot(aes(fill = Diagnosis)) +
  geom_jitter(aes(color = Diagnosis), alpha = 0.5) +
  stat_compare_means(comparisons = my_comparisons, size = 2) +
  plot_themes + 
  facet_wrap(~t_rel_WBC, scales ="free_y") +
  ylab("Log value") +
  xlab("Time relative to  WBC peak")  + 
  ggtitle(paste0(this_variable)) + ylim(3,7)

p_var4 <- ggplot(data_all %>% filter(variable == this_variable & is.finite(result)), 
             aes(x = t_rel_WBC, y = result, group = as.factor(MRN))) +
  geom_point(aes(x = t_rel_WBC, y = result, color = Diagnosis)) + 
  geom_line(aes(x = t_rel_WBC, y = result, color = Diagnosis)) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  ggtitle(paste0(this_variable))

if(save_plots){
  save_plot(paste0("./plots_diff_data/",strsplit(this_variable, c("%"," "))[[1]][1],"_1.pdf"), p_var1, base_width = 8, base_height = 6)
  save_plot(paste0("./plots_diff_data/",strsplit(this_variable, c("%"," "))[[1]][1],"_2.pdf"), p_var2, base_width = 8, base_height = 6)
  save_plot(paste0("./plots_diff_data/",strsplit(this_variable, c("%"," "))[[1]][1],"_3.pdf"), p_var3, base_width = 8, base_height = 6)
  save_plot(paste0("./plots_diff_data/",strsplit(this_variable, c("%"," "))[[1]][1],"_4.pdf"), p_var4, base_width = 8, base_height = 6)
  
}

}


MRN_vec = unique(data_all$MRN)

for(i in c(1:length(MRN_vec))){
  df = data_all %>% filter(MRN == MRN_vec[i] & !(variable %in% c("WBC", "WBC_rel")) & t_rel_WBC > -4 & t_rel_WBC <2)
  p <- ggplot(df, aes(x = t_rel_WBC, y = result)) + 
    geom_line(aes(color = variable)) + 
    geom_point(aes(color = variable)) + 
    xlab("Time relative to  WBC peak") + 
    ggtitle(unique(df$diagnosis)) + 
    plot_themes
  save_plot(paste0(plotting.directory ,"/plot_",unique(df$MRN), "_", unique(df$diagnosis), ".pdf"), p, base_width = 8, base_height = 6)
}

