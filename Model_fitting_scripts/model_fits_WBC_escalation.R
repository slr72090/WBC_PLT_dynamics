
## Define functions ---------------------------------------------
#Function to calculate the adjusted R2 from the model based on length of data and num params
r2_adj <-function(mdl,y,n,param = 3) 
{
  adj <- (sum(!is.na(y)) - 1)/(sum(!is.na(y)) - param)
  sum.sq <- (sum(!is.na(y)) - 1) * var(y, na.rm = TRUE)
  rsq <- 1 - (adj * (deviance(mdl)/sum.sq))
  return(rsq)
}


exp_growth <- function(a,b,c,time) {
  y = a + b*exp(c*time)
  return(y)
}

log_growth <- function(a,b,c,time) {
  y = a + b*(1-exp(-c*time))
  return(y)
}

power_law <- function(a,b,c,time){
  y =  a + b*(time^c)
  return(y)
}



## fit model for patient subset  ## ----------------------------------------------------------
these_patients = unique(df_pre_combined$MRN)
main_df = df_pre_combined

df_mod_fits_escalation = data.frame()
df_mod_compare_escalation = data.frame()
df_resid = data.frame()

for(i in c(1:length(these_patients))){
  print(i)
  dat_sub <- main_df %>% filter(MRN == these_patients[i])
  dat_mod <- dat_sub %>% filter(time >= escalation_start & time <= 0) 

  if(nrow(dat_mod) > n_min_observations){
    # mod.1 <- gsl_nls(
    #   fn = WBC ~ a + b*exp(c*time),
    #   data = dat_mod,
    #   start = start,
    #   algorithm = "lm",
    #   control = list(scale = "levenberg")
    # )
    
    mod.1 <- nls_multstart(WBC ~ exp_growth(a, b, c,time),
                           data = dat_mod,
                           iter = 1000,
                           start_lower = c(a = -100, b = 0, c = 0),
                           start_upper = c(a=100, b=100, c=10),
                           supp_errors = 'Y',
                           convergence_count = 1000,
                           algorithm = "lm",
                           control = list(scale = "levenberg"),
                           na.action = na.omit)
    
    mod.2 <- nls_multstart(WBC ~ log_growth(a, b, c,time),
                           data = dat_mod,
                           iter = 1000,
                           start_lower = c(a = -100, b = 0, c = 0),
                           start_upper = c(a=100, b=100, c=10),
                           supp_errors = 'Y',
                           convergence_count = 100,
                           algorithm = "lm",
                           control = list(scale = "levenberg"),
                           na.action = na.omit)
    
    mod.3 <- nls_multstart(WBC ~ power_law(a, b, c,time),
                           data = dat_mod %>% mutate(time = rev(abs(time))),
                           iter = 1000,
                           start_lower = c(a = -100, b = 0, c = 0),
                           start_upper = c(a=100, b=100, c=10),
                           supp_errors = 'Y',
                           convergence_count = 100,
                           algorithm = "lm",
                           control = list(scale = "levenberg"),
                           na.action = na.omit)
    
   
    this_converged = summary(mod.1)$convInfo$isConv
    df_sub <- data.frame(MRN = unique(dat_sub$MRN),
                         a= coef(mod.1)["a"],
                         b = coef(mod.1)["b"],
                         c = coef(mod.1)["c"],
                         R2 = this_R2,
                         Diagnosis = unique(dat_sub$Diagnosis))
    
    df_mod.compare <- data.frame(MRN = unique(dat_sub$MRN),
                                 Diagnosis = unique(dat_sub$Diagnosis),
                                 Model = c("exponential", "logistic", "power law"),
                                 AIC = c(AIC(mod.1), AIC(mod.2), AIC(mod.3)),
                                 R2 = c(r2_adj(mdl = mod.1, y = dat_mod$WBC, n = nrow(dat.mod)), 
                                        r2_adj(mdl = mod.2, y = dat_mod$WBC, n = nrow(dat.mod)),
                                        r2_adj(mdl = mod.3, y = dat_mod$WBC, n = nrow(dat.mod)))) %>% 
      mutate(dAIC = AIC - min(AIC,na.rm=T))
    
    df_mod_fits_escalation <- rbind(df_mod_fits_escalation, df_sub)
    df_mod_compare_escalation = rbind(df_mod_compare_escalation, df_mod.compare)
  }
}



df_mod_fits_escalation$Fit_type = "growth"
rm(df_sub, dat_sub, dat_mod,these_residuals)

my_comparisons = my_comparisons = list(c("exponential", "logistic"), 
                                       c("exponential", "power law"), 
                                       c("logistic", "power law")
)
p_compare_models_AIC = ggplot(df_mod_compare_escalation %>% 
                                filter(dAIC <5) %>% 
                                melt(., id.vars = c("MRN", "Diagnosis", "Model")) %>% 
                                filter(variable =="dAIC"), aes(x = as.factor(Model), y = value)) + 
  geom_boxplot(aes(fill = Model), alpha = 0.5) + 
  geom_jitter(aes(color = Model), alpha = 0.2) + 
  stat_compare_means(comparisons = my_comparisons, size = 3) +
  facet_wrap(~Diagnosis, scale = "free") + 
  plot_themes + 
  xlab("Model") + ylab("delta AIC") + ylim(0,5
                                           )

p_compare_models_R2 = ggplot(df_mod_compare_escalation %>% 
                               melt(., id.vars = c("MRN", "Diagnosis", "Model")) %>% 
                               filter(variable =="R2"), aes(x = as.factor(Model), y = value)) + 
  geom_boxplot(aes(fill = Model), alpha = 0.5) + 
  geom_jitter(aes(color = Model), alpha = 0.2) + 
  stat_compare_means(comparisons = my_comparisons, size = 3) +
  facet_wrap(~Diagnosis) + 
  plot_themes + 
  xlab("Model") + ylab("R2")


## Compare model fits with setpoints - exponential model fixing setpoint ## ------------------------
df_setpoints <- read.csv("../experiments_3_20_2024/setpoints_partial.csv") %>% filter(Marker == "WBC")
df_fit_recovery <- read.csv("../experiments_3_20_2024/WBC_revovery_fits_fixed_setpoints.csv")
these_patients = df_setpoints$MRN
df_mod_fits_setpoint = data.frame()
df_mod_compare_setpoint = data.frame()
for(i in c(18:length(these_patients))){
  print(i)
  dat_sub <- main_df %>% filter(MRN == these_patients[i])
  dat_mod <- dat_sub %>% filter(time >= escalation_start & time <= 0) 
  this_setpoint = df_setpoints[i,]$Setpoint
  decay_const = df_fit_recovery %>% filter(MRN == these_patients[i]) %>% select(c) %>% unlist()

  
  if(nrow(dat_mod) > n_min_observations){
    mod.1 <- gsl_nls(
      fn = WBC ~ a + b*exp(c*time),
      data = dat_mod,
      start = start,
      algorithm = "lm",
      control = list(scale = "levenberg")
    )
    
    mod.2 <- gsl_nls(
      fn = WBC ~ this_setpoint + b*exp(c*time),
      data = dat_mod,
      start = start3,
      algorithm = "lm",
      control = list(scale = "levenberg")
    )
    
    mod.3 <- gsl_nls(
      fn = WBC ~ this_setpoint + b*exp((c-decay_const)*time),
      data = dat_mod,
      start = start3,
      algorithm = "lm",
      control = list(scale = "levenberg")
    )
    
    this_R2 <- r2_adj(mdl = mod.1, y = dat_mod$WBC, n = nrow(dat.mod))
    
    these_residuals <- data.frame(MRN = unique(dat_sub$MRN),
                                  resid = resid(mod.1))
    this_converged = summary(mod.1)$convInfo$isConv
    df_sub <- data.frame(MRN = unique(dat_sub$MRN),
                         a= this_setpoint,
                         b = coef(mod.2)["b"],
                         c = coef(mod.2)["c"],
                         R2 = this_R2,
                         Diagnosis = unique(dat_sub$Diagnosis))
    
    df_mod.compare <- data.frame(MRN = unique(dat_sub$MRN),
                                 Diagnosis = unique(dat_sub$Diagnosis),
                                 Model = c("without setpoint", "with setpoint", "with setpoint and decay constant"),
                                 AIC = c(AIC(mod.1), AIC(mod.2), AIC(mod.3)),
                                 R2 = c(r2_adj(mdl = mod.1, y = dat_mod$WBC, n = nrow(dat.mod)), 
                                        r2_adj(mdl = mod.2, y = dat_mod$WBC, n = nrow(dat.mod)),
                                        r2_adj(mdl = mod.3, y = dat_mod$WBC, n = nrow(dat.mod)))) %>% 
      mutate(dAIC = AIC - min(AIC,na.rm=T))
    use_NA = F
    if(use_NA == T){
      df_mod.compare <- data.frame(MRN = unique(dat_sub$MRN),
                                   Diagnosis = unique(dat_sub$Diagnosis),
                                   Model = c("without setpoint", "with setpoint", "with setpoint and decay constant"),
                                   AIC = c(NA, AIC(mod.2), AIC(mod.3)),
                                   R2 = c(NA, 
                                          r2_adj(mdl = mod.2, y = dat_mod$WBC, n = nrow(dat.mod)),
                                          r2_adj(mdl = mod.3, y = dat_mod$WBC, n = nrow(dat.mod)))) %>% 
        mutate(dAIC = AIC - min(AIC,na.rm=T))
    }
    
    df_mod_fits_setpoint <- rbind(df_mod_fits_setpoint, df_sub)
    df_mod_compare_setpoint = rbind(df_mod_compare_setpoint, df_mod.compare)

  }
}

p_compare_1 = ggplot(df_mod_compare_setpoint %>% 
                     select(-c(AIC,dAIC)) %>% 
                     spread(., key = c("Model"), value = "R2"),
                     aes(x = `with setpoint`, y = `without setpoint`)) + 
  geom_point(aes(color = Diagnosis)) +  
  geom_abline(linetype = 2) + 
  xlim(0,1) + ylim(0,1) +
  plot_themes +
  xlab("R2 for Model with fixed setpoint") +
  ylab("R2 for Model with inferred setpoint")

p_compare_2 = ggplot(df_mod_compare_setpoint %>% 
                       select(-c(AIC,dAIC)) %>% 
                       spread(., key = c("Model"), value = "R2"),
                     aes(x = `with setpoint and decay constant`, y = `without setpoint`)) + 
  geom_point(aes(color = Diagnosis)) +  
  geom_abline(linetype = 2) + 
  xlim(0.5,1) + ylim(0.5,1) +
  plot_themes + 
  xlab("R2 for Model with fixed setpoint and decay constant") +
  ylab("R2 for Model with inferred setpoint")


## Plot residuals
p_resid <- ggplot(df_resid, aes(x = resid)) + 
  geom_histogram() + 
  plot_themes

# Mke a histogram of R2 values ## ----------------------- 
p_hist <- ggplot(df_mod_fits_escalation) + 
  geom_histogram(aes(x = R2), bins = 100, alpha = 0.5) + 
  geom_vline(aes(xintercept = 0.85, color = "red"), linetype = 2, size =1) + 
  theme_classic()

# Compare model fits by diagnosis 
my_comparisons = my_comparisons = list(c("exponential", "linear"), 
                                       c("exponential", "power law"), 
                                       c("linear", "power law")
)
p_compare_models_AIC = ggplot(df_mod_compare_escalation %>% 
                          melt(., id.vars = c("MRN", "Diagnosis", "Model")) %>% 
                          filter(variable =="dAIC"), aes(x = as.factor(Model), y = log(value))) + 
  geom_boxplot(aes(fill = Model), alpha = 0.5) + 
  geom_jitter(aes(color = Model), alpha = 0.2) + 
  stat_compare_means(comparisons = my_comparisons, size = 3) +
  facet_wrap(~Diagnosis) + 
  plot_themes + 
  xlab("Model") + ylab("Log delta AIC")

p_compare_models_R2 = ggplot(df_mod_compare_escalation %>% 
                                melt(., id.vars = c("MRN", "Diagnosis", "Model")) %>% 
                                filter(variable =="R2"), aes(x = as.factor(Model), y = value)) + 
  geom_boxplot(aes(fill = Model), alpha = 0.5) + 
  geom_jitter(aes(color = Model), alpha = 0.2) + 
  stat_compare_means(comparisons = my_comparisons, size = 3) +
  facet_wrap(~Diagnosis) + 
  plot_themes + 
  xlab("Model") + ylab("R2")

# Plot individual fits ## -----------
plot_fits = T
df_mod_fits = df_mod_fits_escalation
these_patients = unique(df_mod_fits$MRN)
if(plot_fits){
  for(i in c(1:length(these_patients))){
    df_sub <- df_pre_combined %>% filter(MRN == 5573492)#these_patients[i]) 
    p_patient = ggplot(df_sub, aes(time, WBC)) +
      geom_point() +
      geom_function(fun = function(x) coef(mod.1)["a"] +coef(mod.1)["b"]*exp(coef(mod.1)["c"]*x), aes(color = "NLS fit")) + 
      #geom_function(fun = function(x) df_mod_fits[i,]$a +df_mod_fits[i,]$b*exp(df_mod_fits[i,]$c*x), aes(color = "NLS fit")) + 
      labs(color = "") + 
      ggtitle(paste0("MRN: ", unique(df_sub$MRN), ", Diagnosis: ", unique(df_sub$Diagnosis))) + 
      plot_themes
    save_plot(p_patient, base_width = 6, base_height = 4, file = paste0("./plots_R_fits/",unique(df_sub$MRN),"_","fits.pdf"))
  }
}


### Now compare fit that includes subtraction term ## -------------------------------------------------------
df_mod_fits_escalation_2 = data.frame()
df_resid_2 = data.frame()

for(i in c(422:length(these_patients))){
  print(i)
  dat_sub <- main_df %>% filter(MRN == these_patients[i])
  dat_mod <- dat_sub %>% filter(time >= escalation_start & time <= 0)# & w_prime > 0) 
  decay_const = df_mod_fits_combined %>% filter(MRN == these_patients[i]) %>% select(c) %>% unlist()
  if(nrow(dat_mod) > n_min_observations & !is.na(decay_const)){
    mod.1 <- gsl_nls(
      fn = WBC ~ a + b*exp((c-decay_const)*time),
      data = dat_mod,
      start = start,
      algorithm = "lm",
      control = list(scale = "levenberg")
    )
    
    this_R2 <- r2_adj(mdl = mod.1, y = dat_mod$WBC, n = nrow(dat.mod))
    
    these_residuals <- data.frame(MRN = unique(dat_sub$MRN),
                                  resid = resid(mod.1))
    this_converged = summary(mod.1)$convInfo$isConv
    df_sub <- data.frame(MRN = unique(dat_sub$MRN),
                         a= coef(mod.1)["a"],
                         b = coef(mod.1)["b"],
                         c = coef(mod.1)["c"],
                         R2 = this_R2,
                         Diagnosis = unique(dat_sub$Diagnosis))
    
    df_mod_fits_escalation_2 <- rbind(df_mod_fits_escalation_2, df_sub)
    df_resid_2 = rbind(df_resid,
                     these_residuals)
  }
}
write.csv("WBC_fits_exponential_with_subtraction_all.csv")
df_mod_fits_compare <- df_mod_fits_escalation %>% 
  left_join(.,df_mod_fits_escalation_2 %>% select(MRN, Diagnosis, c2 = c, R2_2 = R2))
p <- ggplot(df_mod_fits_compare, aes(x = R2 , y = R2_2)) + 
  geom_point(aes(color = Diagnosis)) + 
  xlab("Adjusted R^2 growth process") + ylab("Adjusted R^2 growth-death process") + 
  geom_abline(linetype=2 ) + 
  plot_themes +
  ylim(-1,1)
