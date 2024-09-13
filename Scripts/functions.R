### FUNCTIONS --------------------------------------------------------------------------

# Make function to streamline fitting
fit.models <- function(dat){
  
  for(ii in 1:length(unique(dat$TS))){
    # filter data by TS
    dat %>%
      mutate(era = case_when((Year <=1989) ~ "early",
                             (Year >1989) ~ "late")) %>%
      filter(TS == unique(dat$TS)[ii]) %>%
      group_by(TS, era) %>%
      mutate(N = n()) %>%
      ungroup()-> TS.dat
    
    knts <- c(3, 4)
    
    for(kk in 1:length(knts)){
      
      # fit gams, record pval and AIC
      gam.1 <- gamm(log.recruitment ~ s(unsmoothed, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      gam.2 <- gamm(log.recruitment ~ s(twoyear, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      gam.3 <- gamm(log.recruitment ~ s(threeyear, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      
      p.gam.1 <- signif(summary(gam.1$gam)$s.table[,4],2)
      p.gam.2 <- signif(summary(gam.2$gam)$s.table[,4],2)
      p.gam.3 <- signif(summary(gam.3$gam)$s.table[,4],2)
      
      p.lme.1 <- signif(summary(gam.1$lme)$tTable[2,5],2)
      p.lme.2 <- signif(summary(gam.2$lme)$tTable[2,5],2)
      p.lme.3 <- signif(summary(gam.3$lme)$tTable[2,5],2)
      
      AIC.gam.1 <- AIC(gam.1)
      AIC.gam.2 <- AIC(gam.2)
      AIC.gam.3 <- AIC(gam.3)
      
      
      # Build summary table
      model.out <- rbind(model.out, data.frame(TS = unique(dat$TS)[ii],
                                               knots = knts[kk],
                                               #Model = c("gam.1", "gam.2", "gam.3"),
                                               sst = c("unsmoothed", "twoyear", "threeyear"),
                                               p_lme = c(p.lme.1, p.lme.2, p.lme.3),
                                               p_gam = c(p.gam.1, p.gam.2, p.gam.3),
                                               AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3)))
      
    } # close knot loop
  } # close timeseries loop
  
  # Label best models by timeseries
  model.out %>%
    group_by(TS) %>%
    mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) -> model.out
  
  return(model.out)
}

# Create function to calculate AR1 and CV/SD for 15 year windows for biology/SST timeseries
sum.fun <- function(dat, data.type){
  
  # if(data.type == "SSB"){
  #   dat %>%
  #     rename(Value = SSB) -> dat
  # }else if(data.type == "Recruitment"){
  #   dat %>%
  #     rename(Value = Recruitment) -> dat
  # }else if(data.type == "Catch"){
  #   dat %>% 
  #     rename(Value = Catch) -> dat
  # }else if(data.type == "SST"){
  #   dat %>%
  #     rename(Value = mean.sst) %>%
  #     mutate(TS = "sst") -> dat
  # }else{
  #   dat = dat
  # }
  # 
  if(data.type == "SSB"){
    dat %>%
      rename(Value = log.SSB) -> dat
  }else if(data.type == "Recruitment"){
    dat %>%
      rename(Value = log.recruitment) -> dat
  }else if(data.type == "Catch"){
    dat %>% 
      rename(Value = log.catch) -> dat
  }else if(data.type == "SST"){
    dat %>%
      rename(Value = mean.sst) %>%
      mutate(TS = "sst") -> dat
  }else{
    dat = dat
  }
  
  
  if(data.type != "Mature biomass"){
    for(ii in 1:length(unique(dat$TS))){
      
      dat %>%
        filter(TS == unique(dat$TS)[ii]) %>%
        na.omit() -> TS.dat
      
      # Specify sliding window width
      width = 15
      
      # Calculate rolling window AR1
      ar1 <- sapply(rollapply(TS.dat$Value, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
      
      # Calculate rolling window SD
      sd <-  rollapply(TS.dat$Value, width = width, FUN = sd)
      
      # Calculate rolling window mean
      avg <- rollapply(TS.dat$Value, width = width, FUN = mean)
      
      # Calculate CV
      cv <- sd/(avg*100)
      
      # Calculate windows
      window <- seq(min(TS.dat$Year)+(width-1), max(TS.dat$Year), by = 1)
      
      
      # Compile output
      sum.out <- rbind(sum.out, data.frame(TS = unique(dat$TS)[ii],
                                           ar1 = ar1,
                                           sd = sd,
                                           avg = avg,
                                           cv = cv,
                                           window = window))
      
    }
  } else{
    for(ii in 1:length(unique(dat$TS))){ # for crab
      
      dat %>%
        filter(TS == unique(dat$TS)[ii], Type == "mmb") %>%
        na.omit() -> TS.dat.mmb
      
      dat %>%
        filter(TS == unique(dat$TS)[ii], Type == "fmb") %>%
        na.omit() -> TS.dat.fmb
      
      
      # Specify sliding window width
      width = 15
      
      # Calculate AR1
      
      if(unique(dat$TS)[ii] != "Bristol Bay red king crab"){
        
        # Calculate rolling window AR1
        ar1.mmb <- sapply(rollapply(TS.dat.mmb$Value, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
        ar1.fmb <- sapply(rollapply(TS.dat.fmb$Value, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
        
        # Calculate rolling window SD
        sd.mmb <-  rollapply(TS.dat.mmb$Value, width = width, FUN = sd)
        sd.fmb <-  rollapply(TS.dat.fmb$Value, width = width, FUN = sd)
        
        # Calculate rolling window mean
        avg.mmb <- rollapply(TS.dat.mmb$Value, width = width, FUN = mean)
        avg.fmb <- rollapply(TS.dat.fmb$Value, width = width, FUN = mean)
        
        # Calculate CV
        cv.mmb <- sd.mmb/(avg.mmb*100)
        cv.fmb <- sd.fmb/(avg.fmb*100)
        
        # Calculate windows
        window <- seq(min(TS.dat.mmb$Year)+(width-1), max(TS.dat.mmb$Year), by = 1)
        
        
        
      } else{
        # Calculate rolling window AR1
        ar1.mmb <- sapply(rollapply(TS.dat.mmb$Value, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
        
        # Calculate rolling window SD
        sd.mmb <-  rollapply(TS.dat.mmb$Value, width = width, FUN = sd)
        
        # Calculate rolling window mean
        avg.mmb <- rollapply(TS.dat.mmb$Value, width = width, FUN = mean)
        
        # Calculate CV
        cv.mmb <- sd.mmb/(avg.mmb*100)
        
        # Calculate windows
        window <- seq(min(TS.dat.mmb$Year)+(width-1), max(TS.dat.mmb$Year), by = 1)
        
        ar1.fmb <- sd.fmb <- avg.fmb <- cv.fmb <- NA
        
      }
      
      
      # Compile output
      sum.out <- rbind(sum.out, data.frame(TS = unique(dat$TS)[ii],
                                           ar1.mmb = ar1.mmb, 
                                           ar1.fmb = ar1.fmb, 
                                           sd.mmb = sd.mmb, 
                                           sd.fmb = sd.fmb,
                                           avg.mmb = avg.mmb,
                                           avg.fmb = avg.fmb,
                                           cv.mmb = cv.mmb,
                                           cv.fmb = cv.fmb,
                                           window = window))
      
    }
  }
  
  
  return(sum.out)
}

# Function to detrend data, calculate ar1, SD, avg, and cv for 15 year windows
trend.fun <- function(TS.dat, data.type){
  
  # # Specify response variable based on data type
  # if(data.type == "Recruitment"){
  #   TS.dat %>%
  #     mutate(Value = log.recruitment) -> TS.dat
  # }else{
  #   TS.dat %>%
  #     mutate(Value = mean.sst) -> TS.dat
  # }
  
  if(data.type == "Recruitment"){
    names = unique(TS.dat$TS)
    
    out <- data.frame()
    summary <- data.frame()
    
    for (ii in 1:length(names)){
      
      dat <- TS.dat %>% filter(TS == names[ii]) %>%
        na.omit()
      
      # Detrend data
      detrend.dat <- loess(log.recruitment ~ Year, dat, span = 0.25, degree = 1)
      
      # Extract residuals
      detrend.dat.resid <- data.frame(TS = rep(names[ii], length(unique(dat$Year))),
                                      Year = dat$Year, log.recruitment = detrend.dat$residuals)
      
      # Calculate AR1 and CV/SD on detrended data
      sum.fun(detrend.dat.resid, data.type) -> out
      
      rbind(summary, out) -> summary
      
    }
  } else{
    
    dat <- TS.dat
    
    # Detrend data
    detrend.dat <- loess(mean.sst ~ Year, dat, span = 0.25, degree = 1)
    
    # Extract residuals
    detrend.dat.resid <- data.frame(TS = rep("SST", length(unique(dat$Year))),
                                    Year = dat$Year, mean.sst = detrend.dat$residuals)
    
    # Calculate AR1 and CV/SD on detrended data
    sum.fun(detrend.dat.resid, data.type) -> out
    
    summary <- out
  }
  return(summary)
}


