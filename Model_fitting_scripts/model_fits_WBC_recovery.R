
## Define functions ---------------------------------------------

exp_decay <- function(a,b,c,time) {
  y = a + b*exp(-c*time)
  return(y)
}

#Function to calculate the adjusted R2 from the model based on length of data and num params
r2_adj <-function(mdl,y,n,param = 3) 
{
  adj <- (sum(!is.na(y)) - 1)/(sum(!is.na(y)) - param)
  sum.sq <- (sum(!is.na(y)) - 1) * var(y, na.rm = TRUE)
  rsq <- 1 - (adj * (deviance(mdl)/sum.sq))
  return(rsq)
}

## fit model for patient subset  ## ----------------------------------------------------------
for(i in c(1:length(these_patients))){
  print(i)
  dat_sub <- main_df %>% filter(MRN == these_patients[i])
  dat_mod <- dat_sub %>% filter(time >= time_bounds[1] & time <= time_bounds[2])# & w_prime > 0) 
  if(nrow(dat_mod) > n_min_observations){
    # mod.1 <- gsl_nls(
    #   fn = Result ~ a + b*exp(-c*time),
    #   data = dat_mod,
    #   start = start,
    #   algorithm = "lm",
    #   control = list(scale = "levenberg")
    # )
    # 

     mod.1 <- nls_multstart(Result ~ exp_decay(a, b, c,time),
                         data = dat_mod,
                         iter = 1000,
                         start_lower = c(a = -100, b = 0, c = 0),
                         start_upper = c(a=100, b=100, c=10),
                         supp_errors = 'Y',
                         convergence_count = 100,
                         algorithm = "lm",
                         control = list(scale = "levenberg"),
                         na.action = na.omit)
    
    
    this_R2 <- r2_adj(mdl = mod.1, y = dat_mod$Result, n = nrow(dat.mod))
    
    these_residuals <- data.frame(MRN = unique(dat_sub$MRN),
                                  resid = resid(mod.1))
    this_converged = summary(mod.1)$convInfo$isConv
    df_sub <- data.frame(MRN = unique(dat_sub$MRN),
                         a= coef(mod.1)["a"],
                         b = coef(mod.1)["b"],
                         c = coef(mod.1)["c"],
                         R2 = this_R2)
    
    df_mod_fits <- rbind(df_mod_fits, df_sub)
    df_resid = rbind(df_resid,
                     these_residuals)
  }
}

rm(df_sub, dat_sub, dat_mod,these_residuals)

