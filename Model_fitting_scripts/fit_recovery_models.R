library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)


# Starting guess parameters for model fit## ------------------------------------------------------
c = 0.5
b = 10
a = 8

start <- list(a = a, b = b, c = c)
fit_models = F

# MI model fits ## ------------------------------------------------
df_mod_fits <- data.frame()
df_resid  <- data.frame()

these_patients = MRNs_MI$MRN
main_df = dat_WBC_MI
source("./Model_fitting_scripts/model_fits_WBC_recovery.R")

df_mod_fits_MI = df_mod_fits
df_resid_MI = df_resid
df_mod_fits_MI$Diagnosis = "MI"
write.csv(df_mod_fits_MI, "WBC_fits_MI.csv")
rm(df_mod_fits, df_resid)

# Stroke model fits ## ------------------------------------------------
df_mod_fits <- data.frame()
df_resid  <- data.frame()

these_patients = MRNs_stroke$MRN
main_df = dat_WBC_stroke
source("./Model_fitting_scripts/model_fits_WBC_recovery.R")

df_mod_fits_stroke = df_mod_fits
df_resid_stroke = df_resid

df_mod_fits_stroke$Diagnosis = "stroke"
write.csv(df_mod_fits_stroke, "WBC_fits_stroke.csv")
rm(df_mod_fits, df_resid)

# Sepsis model fits ## ---------------------------------------------
df_mod_fits <- data.frame()
df_resid  <- data.frame()

these_patients = MRNs_sepsis$MRN
main_df = dat_WBC_sepsis
source("./Model_fitting_scripts/model_fits_WBC_recovery.R")
df_mod_fits_sepsis = df_mod_fits
df_resid_sepsis = df_resid
df_mod_fits_sepsis$Diagnosis = "sepsis"
df_resid_sepsis$Diagnosis = "sepsis"
write.csv(df_mod_fits_sepsis, "WBC_fits_sepsis.csv")
rm(df_mod_fits, df_resid)

# Cdiff model fits ## ---------------------------------------------
# Data frame to store model fits and R^2 
df_mod_fits <- data.frame()
df_resid <- data.frame()
these_patients = MRNs_cdiff$MRN
main_df <- dat_WBC_cdiff

source("./Model_fitting_scripts/model_fits_WBC_recovery.R")

df_mod_fits_cdiff = df_mod_fits
df_resid_cdiff = df_resid
df_mod_fits_cdiff$Diagnosis = "cdiff"
df_resid_cdiff$Diagnosis = "cdiff"

write.csv(df_mod_fits_cdiff, "WBC_fits_cdiff.csv")
rm(df_mod_fits, df_resid)


df_mod_fits_WBC <- rbind(df_mod_fits_sepsis %>% mutate(Diagnosis = "sepsis"),
                     df_mod_fits_cdiff %>% mutate(Diagnosis = "cdiff"),
                     df_mod_fits_stroke %>% mutate(Diagnosis = "stroke"),
                     df_mod_fits_MI %>% mutate(Diagnosis = "MI"))

write.csv(df_mod_fits_WBC, "WBC_recovery_fits.csv")

## Combine fits 

MRN_best_fits <- df_mod_fits_WBC %>% filter(R2 > GOF_threshold) %>% select(MRN,Diagnosis)
dat_PLT_fits <- dat_PLT %>% filter(MRN %in% MRN_best_fits$MRN) %>% unique()

#Check how many pts included for each diagnosis 
dat_PLT_fits %>% group_by(Diagnosis) %>% 
  summarise(n = length(unique(MRN)))

# Fit for MI patients ## ---------------------------
these_patients = MRN_best_fits %>% filter(Diagnosis == "MI" & MRN %in% dat_PLT$MRN & MRN %in% MRNs_MI$MRN) %>% select(MRN) %>% unlist()
this_diagnosis = "MI"
window_length = 5
df_plt_fits = data.frame()

source("./Model_fitting_scripts/model_fits_PLT_recovery.R")

df_plt_fits_MI <- df_plt_fits
rm(df_plt_fits)

# Fit for stroke patients ## ---------------------------
these_patients = MRN_best_fits %>% filter(Diagnosis == "stroke" & MRN %in% dat_PLT$MRN) %>% select(MRN) %>% unlist()
this_diagnosis = "stroke"
window_length = 5
df_plt_fits = data.frame()

source("./Model_fitting_scripts/model_fits_PLT_recovery.R")

df_plt_fits_stroke <- df_plt_fits
rm(df_plt_fits)

# Fit for sepsis patients ## ---------------------------
these_patients = MRN_best_fits %>% filter(Diagnosis == "sepsis" & MRN %in% dat_PLT$MRN) %>% select(MRN) %>% unlist()
this_diagnosis = "sepsis"
window_length = 5
df_plt_fits = data.frame()

source("./Model_fitting_scripts/model_fits_PLT_recovery.R")

df_plt_fits_sepsis <- df_plt_fits
rm(df_plt_fits)

# Fit for Cdiff patients ## -----------------------------------------
these_patients = MRN_best_fits %>% filter(Diagnosis == "cdiff" & MRN %in% dat_PLT$MRN) %>% select(MRN) %>% unlist()
this_diagnosis = "cdiff"
window_length = 5
df_plt_fits = data.frame()
source("./Model_fitting_scripts/model_fits_PLT_recovery.R")
df_plt_fits_cdiff <- df_plt_fits
rm(df_plt_fits)

df_plt_fits <- rbind(df_plt_fits_sepsis %>% mutate(Diagnosis = "sepsis"),
                     df_plt_fits_cdiff %>% mutate(Diagnosis = "cdiff"),
                     df_plt_fits_stroke %>% mutate(Diagnosis = "stroke"),
                     df_plt_fits_MI %>% mutate(Diagnosis = "MI"))

write.csv(df_plt_fits, "PLT_fits.csv")
