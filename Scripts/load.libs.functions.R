### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(zoo)
library(stats)
library(mgcv)
library(MuMIn)
library(roll)
library(TTR)
library(ncdf4)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(ggplot2)
library(oce)
library(pracma)
library(sf)
library(FactoMineR)
library(cowplot)
library(ggpubr)
library(ggtext)
library(patchwork)

### DIRECTORIES ------------------------------
dir <- "Y:/KOD_Research/Ryznar/Reddening/"


### FUNCTIONS --------------------------------------------------------------------------

# Make function to streamline model fitting
fit.r0.SST.models <- function(dat){
  
  for(ii in 1:length(unique(dat$TS))){
    # filter data by TS
    dat %>%
      mutate(era = case_when((Year <=1989) ~ "early",
                             (Year >1989) ~ "late")) %>%
      filter(TS == unique(dat$TS)[ii]) %>%
      group_by(TS, era) %>%
      mutate(N = n()) %>%
      ungroup()-> TS.dat
    
    knts <- c(3, 4, 5)
    
    for(kk in 1:length(knts)){
      
      # fit gams, record pval and AIC
      gam.1 <- gam(log.recruitment ~ s(unsmoothed, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      gam.2 <- gam(log.recruitment ~ s(twoyear, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      gam.3 <- gam(log.recruitment ~ s(threeyear, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      
      p.gam.1 <- signif(summary(gam.1)$s.table[,4],2)
      p.gam.2 <- signif(summary(gam.2)$s.table[,4],2)
      p.gam.3 <- signif(summary(gam.3)$s.table[,4],2)
      
      # p.lme.1 <- signif(summary(gam.1$lme)$tTable[2,5],2)
      # p.lme.2 <- signif(summary(gam.2$lme)$tTable[2,5],2)
      # p.lme.3 <- signif(summary(gam.3$lme)$tTable[2,5],2)
      
      rsq.gam.1 <- signif(summary(gam.1)$r.sq,2)
      rsq.gam.2 <- signif(summary(gam.2)$r.sq,2)
      rsq.gam.3 <- signif(summary(gam.3)$r.sq,2)
      
      AIC.gam.1 <- AIC(gam.1)
      AIC.gam.2 <- AIC(gam.2)
      AIC.gam.3 <- AIC(gam.3)
      
      
      # Build summary table
      model.out <- rbind(model.out, data.frame(TS = unique(dat$TS)[ii],
                                               knots = knts[kk],
                                               #Model = c("gam.1", "gam.2", "gam.3"),
                                               sst = c("unsmoothed", "twoyear", "threeyear"),
                                               #p_lme = c(p.lme.1, p.lme.2, p.lme.3),
                                               p_gam = c(p.gam.1, p.gam.2, p.gam.3),
                                               AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3),
                                               rsq = c(rsq.gam.1, rsq.gam.2, rsq.gam.3)))
      
    } # close knot loop
  } # close timeseries loop
  
  # Label best models by timeseries
  model.out %>%
    #group_by(TS) %>%
    mutate(sig = BH2(p_gam, alph = 0.05)$BHSig,
           p_gam = p_gam,
           padj = p.adjust(p_gam, method = "fdr")) %>%
    group_by(TS) %>%
    mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) %>%
    filter(BEST == "Y") -> model.out
  
  return(model.out)
}

# Create function to calculate AR1 and CV/SD for 15 year windows for biology/SST timeseries
calc.AR1.SD <- function(dat, data.type, wind){
  
  sum.out <- data.frame()
  
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
  }else if(data.type == "SLP"){
    dat %>%
      rename(Value = SLP.win.anom) %>%
      mutate(TS = "slp") -> dat
  }else{
    dat = dat
  }
  
  
  for(ii in 1:length(unique(dat$TS))){
      
      dat %>%
        filter(TS == unique(dat$TS)[ii]) %>%
        na.omit() -> TS.dat
      
      # Specify sliding window width
      width = wind
      
      # Calculate rolling window AR1
      ar1 <- sapply(rollapply(TS.dat$Value, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
      
      # Calculate rolling window SD
      sd <-  rollapply(TS.dat$Value, width = width, FUN = sd, fill = NA)
      
      # Calculate rolling window mean
      avg <- rollapply(TS.dat$Value, width = width, FUN = mean, fill = NA)
      
      # Calculate CV
      cv <- sd/(avg*100)
      
      # Make data frame of sd-cv
      data.frame(year = unique(TS.dat$Year), val = TS.dat$Value, sd = sd, avg = avg, cv = cv) -> win.dat
      
      # Calculate windows
      win.yr <- na.omit(win.dat) %>% pull(year)
      
      # Make data frame of ar1
      data.frame(year = win.yr, ar1 = ar1) -> ar1.dat
      
      # Join
      left_join(win.dat, ar1.dat) -> win.dat
      
      
      # Compile output
      sum.out <- rbind(sum.out, cbind(win.dat,TS = unique(dat$TS)[ii]))
      
    }
  return(sum.out)
}

# Function to detrend data, calculate ar1, SD, avg, and cv for 15 year windows
trend.fun <- function(TS.dat, data.type, wind){
  
  # # Specify response variable based on data type
  # if(data.type == "Recruitment"){
  #   TS.dat %>%
  #     mutate(Value = log.recruitment) -> TS.dat
  # }else{
  #   TS.dat %>%
  #     mutate(Value = mean.sst) -> TS.dat
  # }
  
  #sum.out <- data.frame()
  
 
  
  if(data.type == "Recruitment"){
    names = unique(TS.dat$TS)
    
    out <- data.frame()
    summary <- data.frame()
    
    for (ii in 1:length(names)){
      
      dat <- TS.dat %>% filter(TS == names[ii]) %>%
        na.omit()
      
      # Detrend data
      detrend.dat <- lm(log.recruitment ~ Lagged.Year, dat)
      
      # Extract residuals
      detrend.dat.resid <- data.frame(TS = rep(names[ii], length(unique(dat$Lagged.Year))),
                                      Year = dat$Lagged.Year, log.recruitment = detrend.dat$residuals)
      
     
      
      # Calculate AR1 and CV/SD on detrended data
      wind %>%
        purrr::map(~calc.AR1.SD(detrend.dat.resid, data.type, .x)) -> out
      
      rbind(summary, as.data.frame(out)) -> summary
      
    }
  } else if(data.type == "SST"){
    
    dat <- TS.dat
    
    # Detrend data
    detrend.dat <- lm(mean.sst ~ Year, dat)
    
    # Extract residuals
    
    
    # # Extract residuals
    # detrend.dat.resid <- data.frame(TS = rep("SST", length(unique(dat$Year))),
    #                                 Year = dat$Year, mean.sst = dat$mean.sst)
    # 
    detrend.dat.resid <- data.frame(TS = rep("SST", length(unique(dat$Year))),
                                    Year = dat$Year, mean.sst = detrend.dat$residuals)

    # Calculate AR1 and CV/SD on detrended data on different window lengths
    wind %>%
    purrr::map(~calc.AR1.SD(detrend.dat.resid, data.type, .x)) -> out
    
    summary <- out
  } else{
    dat <- TS.dat
    
    # Detrend data
    detrend.dat <- lm(SLP.win.anom ~ Year, dat)
    
    # Extract residuals
    detrend.dat.resid <- data.frame(TS = rep("SLP", length(unique(dat$Year))),
                                    Year = dat$Year, SLP.win.anom = dat$SLP.win.anom)
    
    detrend.dat.resid <- data.frame(TS = rep("SLP", length(unique(dat$Year))),
                                    Year = dat$Year, SLP.win.anom = detrend.dat$residuals)
    
    
    # Calculate AR1 and CV/SD on detrended data on different window lengths
    wind %>%
      purrr::map(~calc.AR1.SD(detrend.dat.resid, data.type, .x)) -> out
    
    summary <- out
  }
  #return(list(detrended.data = detrend.dat.resid, ar1.var.summary = summary))
  
  return(ar1.var.summary = summary)
  
}

# Make function to streamline model fitting
assess.trend <- function(dat, data.type){
  
  for(ii in 1:length(unique(dat$TS))){
    # filter data by TS
    dat %>%
      filter(TS == unique(dat$TS)[ii]) %>%
      group_by(TS) %>%
      mutate(N = n()) %>%
      ungroup()-> TS.dat
    
    knts <- c(3, 4, 5)
    
    # specify variance value by data type
    if(data.type == "Recruitment"){
      var.val = TS.dat$cv
    } else{
      var.val = TS.dat$sd
    }
    
    # fit gams, iterating through different knots
    for(kk in 1:length(knts)){
      
      gam.ar1 <- gam(ar1 ~ s(year, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      
      gam.var <- gam(var.val ~ s(year, k = knts[kk]), 
                      data= TS.dat, correlation = corAR1())
      
      
      p.gam.ar1 <- summary(gam.ar1)$s.table[,4]
      p.gam.var <- summary(gam.var)$s.table[,4]
      
      r.gam.ar1 <- summary(gam.ar1)$r.sq
      r.gam.var <- summary(gam.var)$r.sq
      
      # p.lme.ar1 <- signif(summary(gam.ar1$lme)$tTable[2,5],2)
      # p.lme.var <- signif(summary(gam.var$lme)$tTable[2,5],2)
      
      AIC.gam.ar1 <- AIC(gam.ar1)
      AIC.gam.var <- AIC(gam.var)
      
      # Build summary table
      model.out <- rbind(model.out, data.frame(TS = unique(dat$TS)[ii],
                                               knots = knts[kk],
                                               response = c("ar1", "var.val"),
                                               #p_lme = c(p.lme.ar1, p.lme.var),
                                               p_gam = c(p.gam.ar1, p.gam.var),
                                               rsq_gam = c(r.gam.ar1, r.gam.var),
                                               AIC = c(AIC.gam.ar1, AIC.gam.var)))
      
    } # close knot loop
  } # close timeseries loop
 
  # Label best models by timeseries
  model.out %>%
    mutate(sig = BH2(p_gam, alph = 0.05)$BHSig, # adjust alpha across ALL hypotheses test
           p_gam = p_gam,
           padj = p.adjust(p_gam, method = "fdr")) %>%
    group_by(TS, response) %>%
    mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) %>% # designate best models
    filter(BEST == "Y") -> model.out
  
  return(model.out)
}

# Make function for model prediction
model.predict <- function(model.dat, TS.dat, data.type){
  
  for(ii in 1:length(unique(model.dat$TS))){
    for(jj in 1:length(unique(model.dat$response))){
      
      model.dat %>%
        filter(TS == unique(model.dat$TS)[ii],
               response == unique(model.dat$response[jj])) -> model.dat2
      
      TS.dat %>%
        filter(TS == unique(model.dat$TS)[ii]) %>%
        dplyr::select(TS, ar1, cv, sd, year) %>%
        na.omit() -> TS.dat2
      
      
      if(unique(model.dat$response)[jj] == "var.val" & data.type == "Recruitment"){
        value = TS.dat2$cv
      } else if(unique(model.dat$response)[jj] == "var.val" & data.type == "sst"){
        value = TS.dat2$sd
      } else{
        value = TS.dat2$ar1
      }
      
      gam.best <- gam(value~ s(year, k = model.dat2$knots), 
                      data= TS.dat2, correlation = corAR1())
      
      # Predict
      pred.vals <- rbind(pred.vals, data.frame(TS = model.dat2$TS,
                                               response = model.dat2$response,
                                               observed = value,
                                               pred = predict(gam.best, se.fit =TRUE)$fit,
                                               pred.CI = 1.96*(predict(gam.best, se.fit =TRUE)$se.fit),
                                               k = model.dat2$knots,
                                               rsq = model.dat2$rsq_gam,
                                               p_gam = model.dat2$p_gam,
                                               sig = model.dat2$sig,
                                               year= TS.dat2$year))
      
    } # close response loop
  } # close timeseries loop
  
  return(pred.vals)
} # close function


# BH takes p and an FDR as inputs and returns hypothesis tests and the threshold p value
BH2 <- function(p, alpha = 0.05){
  u <- sort(p)
  uThresh <- alpha * (1:length(p))/length(p)
  k <- max(c(which((u<=uThresh)), 0), na.rm = TRUE)
  pCrit <- u[k]
  return(list(BHSig = c(rep(TRUE, k), rep(FALSE, length(p)-k))[order(order(p))]
              , pCrit = pCrit)
  )
}

