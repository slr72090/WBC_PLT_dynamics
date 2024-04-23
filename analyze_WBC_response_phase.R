

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
require(RSQLite)
select <- dplyr::select
spread <- tidyr::spread
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

library(ggplot2)
library(cowplot)
library(ggcorrplot)
library(gslnls)
library(nls.multstart)

## Specify files and directory --------------------------------
working.directory <- "/Users/sranjeva/Desktop/Research/Projects/Higgins_lab/R_code/experiments_4_22_2024"
plotting.directory <- "/Users/sranjeva/Desktop/Research/Projects/Higgins_lab/R_code/experiments_4_22_2024/plots"
setwd(working.directory)

# Specify whether or not to save plots --------------------------------
save_plots = F
fit_recovery_models = F

textSize = 11
source("./Plotting_scripts/plot_themes.R")
time_bounds = c(0,5) # time bounds for recovery models 
escalation_start = -4 # Days prior to WBC peeak for start of escalation phasee 
GOF_threshold = 0.85 # Specify threshold for R2
delta_plt_threshold = 20 #Exclude patients with absolute platelet change under this threshold 
n_min_observations = 3 # minimum number of pre-WBC-peak observations
generate_data = F #pull formatted data from SQLite Database
fit_models = F # load fits 

# Diagnoses for comparison 
my_comparisons = list(c("cdiff", "MI"), 
                      c("cdiff", "sepsis"), 
                      c("cdiff", "stroke"),
                      c("MI", "sepsis"),
                      c("MI", "stroke"),
                      c("sepsis","stroke") 
)



## Load in data ## ----------------------------------------------------------
if(generate_data){
  source("./Data_formatting_scripts/format_data_combined.R")
}
if(!generate_data){
  conn <- dbConnect(RSQLite::SQLite(), "./Data/Data.db")
  df_all <- dbGetQuery(conn, "SELECT * FROM WBC_data") %>% select(-c(t_rel_WBC, Test))
  dat_PLT <- dbGetQuery(conn, "SELECT * FROM PLT_data")
  dat_demog <- dbGetQuery(conn, "SELECT * FROM Demographic_data")
  dat_date_compare <- dbGetQuery(conn, "SELECT * FROM Date_comparisons")
  dat_diff <- dbGetQuery(conn, "SELECT * FROM Differential_data_all") %>% 
    mutate(Date = as.Date(Date),
           Diag_date = as.Date(Diag_date))
  dbDisconnect(conn)
}
## Choose patients for analysis and format data ## -------------------------------------------------

df_subset <- df_all %>% 
  group_by(MRN,Diagnosis) %>% 
  mutate(n_obs_pre = sum(time < 0),
         max_time = max(time,na.rm=T)) %>% 
  ungroup() %>% 
  filter(n_obs_pre > n_min_observations & max_time >= max(time_bounds)) %>% 
  select(MRN, Diagnosis) %>% unique() %>% 
  left_join(.,dat_demog %>% select(MRN, Diagnosis, Mortality))




if(!fit_models){
  df_fits <- read.csv("all_fits_scaled.csv") %>%
    left_join(.,(dat_demog %>% select(MRN,Mortality))) %>% 
    left_join(.,read.csv("PLT_fits.csv") %>% select(MRN, Diagnosis, R2_plt = R2))
}

MRNs_survived= df_subset %>% filter(Mortality == 0) %>% 
  select(MRN) %>% unlist()
MRNs_good_recovery <- df_fits %>% filter(pars == "d" & R2 > GOF_threshold & R2_plt > GOF_threshold) %>% select(MRN) %>% unlist()
MRNs_good_response <- df_fits %>% filter(pars == "c" & R2 > GOF_threshold) %>% select(MRN) %>% unlist()

## Pick a subset of patients 
these_MRNs = unique(df_subset$MRN)
these_MRNs = MRNs_good_response
these_MRNs = intersect(MRNs_good_recovery, MRNs_survived)
these_MRNs =intersect(MRNs_good_recovery, MRNs_good_response)
these_MRNs = intersect(intersect(MRNs_good_recovery, MRNs_survived), MRNs_good_response)

## Time between diagnosis and peak WBC ---------------------------------------------------------------------------------------------

df_date_compare <- dat_date_compare %>% filter(MRN %in% these_MRNs) %>% 
  filter(t_rel == t_max_WBC) %>% 
  mutate(d_date = as.numeric(Date - Diag_date)) %>% unique()

df_date_compare_sum <- df_date_compare %>% 
  group_by(Diagnosis) %>% 
  summarize(mean = mean(d_date,na.rm =T),
            LCI = quantile(d_date,na.rm=T)[2],
            UCI = quantile(d_date,na.rm=T)[4],
            n = n(),
            n_neg = sum(d_date < 0, na.rm =T))

p <- ggplot(df_date_compare, aes(x = as.factor(Diagnosis), y = d_date)) +
  geom_boxplot(aes(fill = as.factor(Diagnosis)), alpha = 0.2) + 
  geom_jitter(aes(color = as.factor(Diagnosis)), alpha = 0.05) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("") + ylab("Time (days) between diagnosis and peak WBC") + 
  plot_themes + labs(color = "", fill = "")

these_MRNs = df_date_compare %>% filter(d_date > -4 & d_date < 10) %>% select(MRN) %>% unlist()

######  Empiric response phase rates   ##-----------------

df_WBC_pre <- df_all %>% 
  filter(time >= escalation_start & time <=0 & MRN %in% these_MRNs) %>% 
  group_by(MRN, Diagnosis) %>% 
  mutate(min_value = min(Result), 
         max_value = max(Result),
         first_value = Result[1]) %>% 
  ungroup() %>% 
  left_join(.,dat_demog %>%  select(MRN, Diagnosis, Mortality), by = c("MRN", "Diagnosis"), relationship = "many-to-many") %>% 
  mutate(w_scaled = (Result - min_value)/(max_value - min_value)) %>% 
  mutate(WBC = Result,
         delta_WBC = Result - first_value)  %>%  
  filter(is.finite(w_scaled))

df_sum_WBC <- df_WBC_pre %>% 
  group_by(time,Diagnosis)%>% 
  summarise(mean = mean(w_scaled, na.rm = T),
            median = median(w_scaled, na.rm =T),
            LCI = quantile(w_scaled, c(0.25,0.75))[1],
            UCI = quantile(w_scaled, c(0.25,0.75))[2], 
            n = n())


p_pre_WBC <- ggplot(df_sum_WBC , aes(x = time, y = median)) +
  geom_line(data = df_WBC_pre , aes( x = time, y = w_scaled, group = as.factor(MRN)), linetype = 2, alpha = 0.02) + 
  geom_point(data = df_sum_WBC,  aes(x = time, y = median, color = Diagnosis)) + 
  geom_text(data = df_sum_WBC, aes(label = n), nudge_y = 0.1) +
  geom_line(data = df_sum_WBC,  aes(x = time, y = median, color = Diagnosis)) + 
  geom_smooth(data = df_sum_WBC,
              method = 'nls',
              formula = y ~ a*exp(r*x),
              method.args = list(start = c(a = 10, r = 0.01)),#list(family = gaussian(link = 'log')),
              aes(color = Diagnosis),
              se = F,
              linetype = 2
  ) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = Diagnosis), alpha = 0.2) + 
  xlim(escalation_start,0) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  ylab("scaled WBC") +
  ylim(0,1.1) + 
  facet_wrap(~Diagnosis)

p_pre_WBC_spline <- ggplot(df_sum_WBC , aes(x = time, y = median)) + 
  geom_smooth(data = df_sum_WBC,
              method = 'nls',
              formula = y ~ a*exp(r*x),
              method.args = list(start = c(a = 10, r = 0.01)),#list(family = gaussian(link = 'log')),
              aes(color = Diagnosis),
              se = F,
              linetype = 2
  ) +
  xlim(escalation_start,0) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  ylab("scaled WBC") +
  ylim(0,1.1)

######  Quantify confidence in escalation start ###  ##-----------------
#Slope of the WBC curve at each point ## ---

df_escalation_slope <-df_all %>% 
  filter(MRN %in% these_MRNs & time >=-6 & time <=0) %>% 
  group_by(MRN, Diagnosis) %>% 
  mutate(min_value = min(Result), 
         max_value = max(Result),
         first_value = Result[1]) %>% 
  mutate(lag_result = lag(Result),
         percent_change = (Result - lag(Result))/Result) %>% 
  filter(!is.na(percent_change)) %>% 
  group_by(time, Diagnosis) %>% 
  summarize(median = median(percent_change, na.rm = T),
            LCI = quantile(percent_change, na.rm =T)[2],
            UCI = quantile(percent_change, na.rm =T)[4],
            n = n(),
            f_5percent = sum(percent_change <= 0.05)/n(),
            f_10percent = sum(percent_change <= 0.1)/n(),
            f_25percent = sum(percent_change <= 0.24)/n())

p_slopes <- ggplot(df_escalation_slope, aes(x= time, y = median)) + 
  geom_line() +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2) + 
  geom_hline(aes(yintercept = 0.05, color = "5%"), linetype = 2) +
  geom_hline(aes(yintercept = 0.1, color = "10%"), linetype = 2) +
  geom_hline(aes(yintercept = 0.25, color = "25%"), linetype = 2) +
  labs(color = "WBC slope") + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  ylab("WBC slope") +
  facet_wrap(~Diagnosis) + xlim(-5,0)

###### Analyze WBC fits ## ----------------------------------------------------------------------------------------
df_fits$Mortality_fctr = recode(as.factor(df_fits$Mortality), "1" = "Died" ,  "0" = "Survived")
df_fits$Diagnosis_fctr = recode(df_fits$Diagnosis, "cdiff" = 1 ,  "sepsis" = 0, "stroke" = 2, "MI" = 3)

my_comparisons = list(c("cdiff", "MI"), 
                      c("cdiff", "sepsis"), 
                      c("cdiff", "stroke"),
                      c("MI", "sepsis"),
                      c("MI", "stroke"),
                      c("sepsis","stroke"))

p_r2_growth <- ggplot(df_fits %>% filter( pars == "c"), aes(x = R2))  +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = GOF_threshold, color = "red", linetype = 2) +
  xlab("Adjusted R^2 for exponential response model")  + labs(fill = "")+ 
  plot_themes

p_r_diagnosis_model_cohort <- ggplot(df_fits %>% filter( pars == "c"& MRN %in% MRNs_good_recovery & MRN %in% MRNs_good_response & Mortality == 0) %>% unique(), 
                        aes(x = reorder(as.factor(Diagnosis), Diagnosis_fctr), y = log(MLE))) +
  geom_boxplot(aes(fill = Diagnosis), alpha = 0.4) + 
  geom_jitter(aes(color =Diagnosis), alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("")  + labs(fill = "", color = "") + ylab("Log of scaled WBC growth rate") + 
  plot_themes

p_f_model_cohort <- ggplot(df_fits %>% filter(model == "decay" & pars == "f" & MRN %in% MRNs_good_recovery & MRN %in% MRNs_good_response & Mortality == 0), 
                           aes(x = as.factor(Diagnosis), y = log(MLE))) +
  geom_boxplot(aes(fill = as.factor(Diagnosis)), alpha = 0.4) + 
  geom_jitter(aes(color = as.factor(Diagnosis)), alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons) + ylab("Log of scaled deecay rate") + xlab("") + 
  labs(fill = "", color = "") + 
  plot_themes

p_r_mortality <- ggplot(df_fits %>% filter(model == "growth" & pars == "c" & MRN %in% MRNs_good_response) %>% unique(), 
                        aes(x = reorder(as.factor(Diagnosis), Diagnosis_fctr), y = log(MLE))) +
  geom_boxplot(aes(fill = as.factor(Mortality_fctr)), alpha = 0.4) + 
  geom_point(aes(color = as.factor(Mortality_fctr)), alpha = 0.2, position=position_jitterdodge()) +
  facet_wrap(~Mortality_fctr) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("")  + labs(fill = "", color = "") + theme(legend.position = "none") + ylab("Log of scaled WBC growth rate") + 
  plot_themes

p_r_mortality_2 <- ggplot(df_fits %>% filter(model == "growth" & pars == "c") %>% unique(), 
                          aes(x =as.factor(Mortality_fctr), y = log(MLE))) +
  geom_boxplot(aes(fill = as.factor(Mortality_fctr)), alpha = 0.4) + 
  geom_point(aes(color = as.factor(Mortality_fctr)), alpha = 0.2, position=position_jitterdodge()) +
  facet_wrap(~Diagnosis, scales = 'free') +
  stat_compare_means() +
  xlab("")  + labs(color = "", fill = "") + theme(legend.position = "none") + ylab("Log of scaled WBC growth rate") + 
  plot_themes


## Association of peak WBC and starting WBC count with inferrred response phase characteristics ##----------------------------

df_WBC_at_peak <- df_all %>% filter(time == 0 & MRN %in% df_subset$MRN & Result < 80) %>% unique() %>% 
  left_join(.,df_fits %>% filter(pars == "c") %>% select(MRN, MLE, Diagnosis, Mortality)) %>% unique() %>% 
  mutate(Diagnosis_fctr = recode(Diagnosis, "cdiff" = 1 ,  "sepsis" = 0, "stroke" = 2, "MI" = 3),
         Mortality_fctr = recode(Mortality, "1" = "Died" ,  "0" = "Survived")) 

df_sum_peak <- df_WBC_at_peak %>% 
  group_by(Diagnosis, Mortality) %>% 
  summarize(mean = mean(Result),
            sd = sd(Result),
            metric = "Peak WBC")

df_WBC_at_start <- df_all %>% filter(time == -4 & MRN %in% df_subset$MRN & Result < 80) %>% unique() %>% 
  left_join(.,df_fits %>% filter(pars == "c") %>% select(MRN, MLE, Diagnosis, Mortality)) %>% unique() %>% 
  mutate(Diagnosis_fctr = recode(Diagnosis, "cdiff" = 1 ,  "sepsis" = 0, "stroke" = 2, "MI" = 3),
         Mortality_fctr = recode(Mortality, "1" = "Died" ,  "0" = "Survived")) 

df_sum_start <- df_WBC_at_peak %>% 
  group_by(Diagnosis, Mortality) %>% 
  summarize(mean = mean(Result),
            sd = sd(Result),
            metric = "starting WBC")

df_sum_WBC <- rbind(df_sum_peak, df_sum_start) %>% 
  mutate(Diagnosis_fctr = recode(Diagnosis, "cdiff" = 1 ,  "sepsis" = 0, "stroke" = 2, "MI" = 3),
         Mortality_fctr = recode(Mortality, "1" = "Died" ,  "0" = "Survived")) 

p_peak_model_cohort <- ggplot(df_WBC_at_peak  %>% filter(MRN %in% MRNs_good_recovery & MRN %in% MRNs_good_response & Mortality == 0),
                              aes(x = reorder(Diagnosis, Diagnosis_fctr), y = Result)) +
  geom_boxplot(aes(fill = as.factor(Diagnosis)), alpha = 0.4) + 
  geom_jitter(aes(color = as.factor(Diagnosis)), alpha = 0.1) + 
  stat_compare_means(comparisons = "my_comparisons")+
  xlab("")  + labs(color = "", fill = "") + theme(legend.position = "none") + ylab("Peak WBC count") + 
  plot_themes

p_start_model_cohort <- ggplot(df_WBC_at_start  %>% filter(MRN %in% MRNs_good_recovery & MRN %in% MRNs_good_response & Mortality == 0),
                              aes(x = reorder(Diagnosis, Diagnosis_fctr), y = Result)) +
  geom_boxplot(aes(fill = as.factor(Diagnosis)), alpha = 0.4) + 
  geom_jitter(aes(color = as.factor(Diagnosis)), alpha = 0.1) + 
  stat_compare_means(comparisons = "my_comparisons")+
  xlab("")  + labs(color = "", fill = "") + theme(legend.position = "none") + ylab("Starting WBC count") + 
  plot_themes

p_peak_mortality <- ggplot(df_WBC_at_peak  %>% filter(MRN %in%  MRNs_good_response),
                              aes(x = Mortality_fctr, y = Result)) +
  geom_boxplot(aes(fill = as.factor(Mortality)), alpha = 0.4) + 
  geom_point(aes(color = as.factor(Mortality)), alpha = 0.2, position=position_jitterdodge())+ 
  stat_compare_means()+ facet_wrap(~Diagnosis) +
  xlab("")  + labs(color = "", fill = "") + theme(legend.position = "none") + ylab("Peak WBC count") + 
  plot_themes

p_correlation_model_cohort <- ggplot(df_WBC_at_peak  %>% filter(MRN %in% MRNs_good_recovery & MRN %in% MRNs_good_response & Mortality == 0),
                                     aes(x = Result, y = log(MLE))) + 
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap(~Diagnosis, scales = "free") + xlab("Max WBC count") + ylab("Scaled exponential WBC growth constant, log scale") + 
  plot_themes
                          

test <- lm(MLE ~ Result + Diagnosis, data = df_WBC_at_peak )

## Plot WBC fit characteristics --------------------------------------------------------------
p <- ggplot(df_sum_WBC, aes(x = reorder(Diagnosis, Diagnosis_fctr), y = mean)) + 
  geom_point(aes(color = as.factor(Mortality))) + 
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd, color = as.factor(Mortality)),, width = 0.1) +
  facet_grid(metric ~.) + ylim(0,50) + labs(color = "") + xlab("") + ylab("Mean (sd)") +
  plot_themes
 
df_compare_rates <- df_risk_stratification %>%
  filter(pars %in% c("c","f") & R2 >GOF_threshold & R2_plt >GOF_threshold & Mortality == 0) %>% 
  select(MRN, pars, MLE, Mortality, Diagnosis) %>% unique() %>% 
  spread(key = pars, value= MLE) %>% 
  na.omit()


p <- ggplot(df_compare_rates , aes(x = log(c), y = log(f))) +
  geom_point() + 
  stat_smooth(method = "lm", linetype = 2) + 
  facet_wrap(~Diagnosis, scales = "free") + 
  plot_themes + xlab("Log of scaled growth rate") + ylab("Log of scaled decay rate")


make_individual_plot <- function(this_MRN, main_df, fits_df, save_plots = T){
  dat_mod <- main_df %>% filter(MRN == this_MRN) %>% select(MRN,Diagnosis, time, w_scaled) %>% unique()
  pars = fits_df %>% filter(MRN == this_MRN) %>% select(MLE) %>% unlist()
  names(pars) <- fits_df %>% filter(MRN == this_MRN) %>% select(pars) %>% unlist()
  times = seq(min(dat_mod$time), max(dat_mod$time), by = 0.01)
  p_patient = ggplot(dat_mod, aes(x = time,y =  w_scaled)) +
    geom_point() +
    geom_function(fun = function(x) (pars["b"]*(exp(pars["c"]*x))), aes(color = "response fit"), linetype = 2) + 
    geom_function(fun = function(x) (pars["d"]*(exp(-pars["f"]*x))), aes(color = "decay fit"), linetype = 2) + 
    labs(color = "") + ylim(0,1.25) + 
    xlab("Time") + ylab("scaled WBC") + 
    ggtitle(paste0("MRN: ", unique(dat_mod$MRN), ", Diagnosis: ", unique(dat_mod$Diagnosis))) + 
    plot_themes 
  
  if(save_plots == T){
    save_plot(paste0("./plots/",unique(dat_mod$Diagnosis), "_", this_MRN, ".pdf"), p_patient, base_width = 8, base_height = 6)
  }
  
}
plot_individuals = F
if(plot_individuals){
  n_plots
  sapply(sample(unique(model_cohort$MRN), n_plots), 
         make_individual_plot, 
         main_df = df_data_scaled, 
         fits_df = df_fits)
}

### Analyze differential data ## -------------------------------------------------
generate_diff_data = F
time_seq = c(escalation_start:1)
if(generate_diff_data){
  source("format_diff_data.R")
}
data_diff = dat_diff
dat_diff <- data_diff %>% filter(MRN %in% these_MRNs)

df_stats <- dat_diff %>% 
  group_by(variable, Diagnosis, t_rel_WBC) %>% 
  summarize(n = n())


dat_diff <- dat_diff %>% mutate(Mortality_fctr = recode(as.factor(Mortality), "0" = "Survived", "1" = "Died"))
dat_sum <- dat_sum %>% mutate(Mortality_fctr = recode(as.factor(Mortality), "0" = "Survived", "1" = "Died"))


dat_sum = dat_diff %>% 
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

p <- ggplot(dat_sum %>% filter(t_rel_WBC %in% time_seq), aes(x = t_rel_WBC, y = median)) +
  geom_line( aes(x = t_rel_WBC, y = median, color = Diagnosis)) + 
  geom_point( aes(x = t_rel_WBC, y = median, color = Diagnosis)) +
 # geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = Diagnosis), alpha = 0.2) + 
  geom_text(aes(label = n)) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  facet_wrap(~variable + Mortality, scales = "free_y") + 
  xlim(-4,0)

p <- ggplot(dat_sum %>% filter(t_rel_WBC %in% time_seq & Diagnosis == "MI"), aes(x = t_rel_WBC, y = median)) +
  geom_line( aes(x = t_rel_WBC, y = median, color = as.factor(Mortality_fctr))) + 
  geom_point( aes(x = t_rel_WBC, y = median, color = as.factor(Mortality_fctr))) +
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = as.factor(Mortality_fctr)), alpha = 0.2) + 
  geom_text(aes(label = n)) + 
  plot_themes + 
  xlab("Time relative to  WBC peak") + 
  facet_wrap(~variable, scales = "free_y") + 
  xlim(-4,1)

dat_sum_all = dat_sum 
dat_diff_all = dat_diff
#dat_sum = dat_sum_all %>% filter(Diagnosis == "MI")
#dat_diff = dat_diff_all %>% filter(Diagnosis == "MI")
for(i in c(1:length(unique(dat_diff$variable)))){
  this_variable = unique(dat_diff$variable)[i]
  
  p_var1 <- ggplot(dat_sum   %>% filter(variable == this_variable & t_rel_WBC %in% time_seq), aes(x = t_rel_WBC, y = median)) +
    geom_line( aes(x = t_rel_WBC, y = median, color = Mortality_fctr)) + 
    geom_point( aes(x = t_rel_WBC, y = median, color = Mortality_fctr)) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = Mortality_fctr), alpha = 0.2) + 
    plot_themes + 
    ggtitle(paste0(this_variable)) + 
    facet_wrap(~Diagnosis) +
    xlab("Time relative to  WBC peak") +
    xlim(min(time_seq), max(time_seq))
  
  p_var2 <- ggplot(dat_diff %>% filter(variable == this_variable &t_rel_WBC %in% time_seq), aes(x = as.factor(t_rel_WBC), y = log(result), fill = Diagnosis)) +
    geom_boxplot(aes(fill = as.factor(Diagnosis)), alpha = 0.4) +
    geom_point(aes(color = as.factor(Diagnosis)), alpha = 0.2, position=position_jitterdodge()) +
    ylab("Log value") + 
    plot_themes + 
    stat_compare_means(comparisons = my_comparisons, size = 2, aes(label = paste0("p = ", after_stat(p.format)))) + 
    facet_grid(Mortality_fctr~.) + 
    xlab("Time relative to  WBC peak")  + 
    ggtitle(paste0(this_variable)) + labs(color = "") + theme(legend.position = "none") 
  
  p_var5 <- ggplot(dat_diff %>% filter(variable == this_variable &t_rel_WBC %in% time_seq), aes(x = as.factor(t_rel_WBC), y = log(result), fill = Mortality_fctr)) +
    geom_boxplot(aes(fill = as.factor(Mortality_fctr)), alpha = 0.4) +
    geom_point(aes(color = as.factor(Mortality_fctr)), alpha = 0.2, position=position_jitterdodge()) +
    ylab("Log value") + 
    stat_compare_means(size = 3, aes(label = paste0("p = ", after_stat(p.format)))) + 
    plot_themes + 
    facet_grid(Diagnosis~.) + 
    xlab("Time relative to  WBC peak")  + 
    ggtitle(paste0(this_variable)) + 
    labs(color = "", fill = "") 
  
  
  p_var3 <- ggplot(dat_diff %>% filter(variable == this_variable &t_rel_WBC %in% time_seq), aes(x = as.factor(t_rel_WBC), y = log(result), fill = Mortality_fctr)) +
    geom_boxplot(aes(fill = Mortality_fctr), alpha = 0.4) +
    geom_point(aes(color = Mortality_fctr), alpha = 0.2, position=position_jitterdodge()) +
    ylab("Log value") + 
    plot_themes + 
    facet_grid(Diagnosis~., scales = "free_y") + 
    xlab("Time relative to  WBC peak")  + 
    ggtitle(paste0(this_variable)) + labs(fill = "", color = "")
  
  p_var3 <- ggplot(dat_diff %>% filter(variable == this_variable &   Mortality == 0 &t_rel_WBC %in% time_seq), aes(x = as.factor(Diagnosis), y = log(result), fill = Diagnosis)) +
    geom_jitter(aes(color = Diagnosis), alpha = 0.2) +
    geom_boxplot(aes(fill = Diagnosis)) +
    stat_compare_means(comparisons = my_comparisons, size = 2) +
    plot_themes + 
    facet_wrap(~t_rel_WBC, scales ="free_y") +
    ylab("Log value") +
    xlab("Time relative to  WBC peak")  + 
    ggtitle(paste0(this_variable, ", survivors")) 
  
  p_var4<- ggplot(dat_diff %>% filter(variable == this_variable &   Mortality == 1 &t_rel_WBC %in% time_seq), aes(x = as.factor(Diagnosis), y = log(result), fill = Diagnosis)) +
    geom_jitter(aes(color = Diagnosis), alpha = 0.2) +
    geom_boxplot(aes(fill = Diagnosis)) +
    stat_compare_means(comparisons = my_comparisons, size = 2) +
    plot_themes +
    facet_wrap(~t_rel_WBC, scales ="free_y") +
    ylab("Log value") +
    xlab("Time relative to  WBC peak")  +
    ggtitle(paste0(this_variable, ", non-survivors"))

  if(save_plots){
    save_plot(paste0("./plots_diff_data/all_",strsplit(this_variable, c("%"," "))[[1]][1],"_1.pdf"), p_var1, base_width = 8, base_height = 6)
    save_plot(paste0("./plots_diff_data/all_",strsplit(this_variable, c("%"," "))[[1]][1],"_2.pdf"), p_var2, base_width = 8, base_height = 6)
    save_plot(paste0("./plots_diff_data/all_",strsplit(this_variable, c("%"," "))[[1]][1],"_3.pdf"), p_var3, base_width = 8, base_height = 6)
    save_plot(paste0("./plots_diff_data/all_",strsplit(this_variable, c("%"," "))[[1]][1],"_4.pdf"), p_var4, base_width = 8, base_height = 6)
    save_plot(paste0("./plots_diff_data/all_",strsplit(this_variable, c("%"," "))[[1]][1],"_5.pdf"), p_var5, base_width = 8, base_height = 6)
 
  }
  
}

