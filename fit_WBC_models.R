require(nls_multstart)
## Define functions ## -----------------------------------
r2_adj <-function(mdl,y,n,param = nparam) 
{
  adj <- (sum(!is.na(y)) - 1)/(sum(!is.na(y)) - param)
  sum.sq <- (sum(!is.na(y)) - 1) * var(y, na.rm = TRUE)
  rsq <- 1 - (adj * (deviance(mdl)/sum.sq))
  return(rsq)
}

exp_growth <- function(b,c,time) {
  y = b*exp(c*time)
  return(y)
}

exp_decay <- function(a,d,f,time){
  y =  a + d*(exp(-f*time))
  return(y)
}

exp_decay_scaled <- function(d,f,time){
  y =  d*(exp(-f*time))
  return(y)
}


## Choose patients for analysis and format data ## --------------------------------------------------
df_data_scaled <- df_all %>% 
  filter(time >= escalation_start & time <=0 & MRN %in% these_MRNs) %>% 
  group_by(MRN, Diagnosis) %>% 
  mutate(min_value = min(Result), 
         max_value = max(Result),
         first_value = Result[1]) %>% 
  ungroup() %>% 
   mutate(w_scaled = (Result - min_value)/(max_value - min_value)) %>% 
  filter(is.finite(w_scaled)) %>% 
  bind_rows(.,
            df_all %>% 
              filter(time >= 0 & time <= max(time_bounds) & MRN %in%these_MRNs) %>% 
              group_by(MRN, Diagnosis) %>% 
              mutate(min_value = min(Result), 
                     max_value = max(Result),
                     first_value = Result[1]) %>% 
              ungroup() %>% 
              mutate(w_scaled = (Result - min_value)/(max_value - min_value)) %>% 
              filter(is.finite(w_scaled)))



## Make the fits for each patient ## -------------------------------------------------------------
if(fit_models){
  df_fits = data.frame()
  MRNs_df = df_data_scaled %>% select(MRN, Diagnosis) %>% unique() %>% 
    filter(MRN %in% test_MRNs)
  main_df <- df_data_scaled
  for(i in c(1:length(MRNs_df$MRN))){
    cat(i)
    dat_mod <- main_df %>% filter(MRN == MRNs_df[i,]$MRN)
    this_diagnosis = unique(dat_mod$Diagnosis)
    mod.1 <- nls_multstart(w_scaled ~ exp_growth(b,c,time),
                           data = dat_mod %>% filter(time <=0),
                           iter = 100,
                           start_lower = c( b = 0, c = 0),
                           start_upper = c(b = 10, c = 10),
                           lower = c(b = 0, c = 0), 
                           supp_errors = 'N',
                           convergence_count = 1000,
                           algorithm = "lm",
                           control = list(scale = "levenberg"),
                           na.action = na.omit)
    
    # mod.2 <- nls_multstart(Result ~ exp_decay(a,d,f,time),
    #                        data = dat_mod %>% filter(time >=0),
    #                        iter = 100,
    #                        start_lower = c(a = 0, b = 0, c = 0),
    #                        start_upper = c(a=100, b = 10, c = 10),
    #                        lower = c(a = 0, b = 0, c = 0),
    #                        supp_errors = 'N',
    #                        convergence_count = 1000,
    #                        algorithm = "lm",
    #                        control = list(scale = "levenberg"),
    #                        na.action = na.omit)
    
    mod.2 <- nls_multstart(w_scaled ~ exp_decay_scaled(d,f,time),
                           data = dat_mod %>% filter(time >=0),
                           iter = 100,
                           start_lower = c(d = 0,f= 0),
                           start_upper = c(d = 10, f = 10),
                           lower = c(d = 0, f = 0),
                           supp_errors = 'N',
                           convergence_count = 1000,
                           algorithm = "lm",
                           control = list(scale = "levenberg"),
                           na.action = na.omit)
    
    # if(is.null(mod.2)){
    #   mod.2 <- nls_multstart(Result ~ exp_decay(a,d,f,time),
    #                          data = dat_mod %>% filter(time >=0),
    #                          iter = 1000,
    #                          start_lower = c(a = 0, b = 0, c = 0),
    #                          start_upper = c(a=100, b = 10, c = 10),
    #                          lower = c(a = 0, b = 0, d = 0),
    #                          supp_errors = 'N',
    #                          convergence_count = 1000,
    #                          algorithm = "lm",
    #                          control = list(scale = "levenberg"),
    #                          na.action = na.omit)
    # }
    
    r2_growth  = r2_adj(mdl = mod.1, y = dat_mod$w_scaled, n = nrow(dat.mod), param = length(coef(mod.1)))
    r2_decay  = r2_adj(mdl = mod.2, y = dat_mod$w_scaled, n = nrow(dat.mod), param = length(coef(mod.2)))
  
    df_fit_ind = data.frame(MRN = MRNs_df[i,]$MRN,
                            model = "growth",
                            pars = c("b", "c"),
                            #LCI = confint(mod.1)[,1],
                            #UCI = confint(mod.1)[,2],
                            MLE = coef(mod.1),
                            R2 = r2_growth) %>% 
      bind_rows(., data.frame(MRN = MRNs_df[i,]$MRN,
                              model = "decay",
                              pars = c("d", "f"),
                              #LCI = confint(mod.2)[,1],
                              #UCI = confint(mod.2)[,2],
                              MLE = coef(mod.2),
                              R2 = r2_decay) )%>% 
      mutate(Diagnosis = this_diagnosis)
    
    df_fits <- rbind(df_fits,
                     df_fit_ind) 
    
  }
  
  
  write.csv(df_fits, "all_fits_scaled_2.csv")
}
