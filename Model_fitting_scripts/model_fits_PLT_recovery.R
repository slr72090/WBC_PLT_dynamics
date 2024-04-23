r2_adj <-function(mdl,y,n,param = 3) 
{
  adj <- (sum(!is.na(y)) - 1)/(sum(!is.na(y)) - param)
  sum.sq <- (sum(!is.na(y)) - 1) * var(y, na.rm = TRUE)
  rsq <- 1 - (adj * (deviance(mdl)/sum.sq))
  return(rsq)
}


for(i in c(1:length(these_patients))){
  print(i)
  dat_sub <- dat_PLT %>% filter(Diagnosis == this_diagnosis & MRN == these_patients[i])
  df_this_fit = data.frame()
  for(j in c(0:(nrow(dat_sub) - window_length))){
    shift = j
    dat_mod <- dat_sub %>% filter(time >= time_bounds[1] + j & time <= time_bounds[2] + j)# & w_prime > 0) 
    if(nrow(dat_mod) > 4){
      mod.lm <- lm(Result ~time, data = dat_mod)
      R2 = summary(mod.lm)$adj.r.squared
      df_this_fit <- rbind(df_this_fit, 
                           data.frame(MRN = unique(dat_sub$MRN),
                                      shift = j,
                                      m = coef(mod.lm)[2],
                                      b = coef(mod.lm)[1],
                                      R2 = R2))
      df_this_fit_sub <- df_this_fit
    }
    if(nrow(dat_mod) < 4){
      df_this_fit <- rbind(df_this_fit,
                           data.frame(MRN = unique(dat_sub$MRN),
                                      shift = j,
                                      m = NA,
                                      b = NA,
                                      R2 = 0))
      df_this_fit_sub <- df_this_fit
      
    }
  }
  df_plt_fits <- rbind(df_plt_fits, df_this_fit_sub %>% filter(R2 == max(R2)))
  
}
df_plt_fits$Diagnosis = "MI"
write.csv(df_plt_fits, "PLT_fits_MI.csv")