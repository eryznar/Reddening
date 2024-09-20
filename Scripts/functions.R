### FUNCTIONS --------------------------------------------------------------------------

# Make function to streamline model fitting
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
    
    knts <- c(3, 4, 5)
    
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
      
      rsq.gam.1 <- signif(summary(gam.1$gam)$r.sq,2)
      rsq.gam.2 <- signif(summary(gam.2$gam)$r.sq,2)
      rsq.gam.3 <- signif(summary(gam.3$gam)$r.sq,2)
      
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
                                               AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3),
                                               rsq = c(rsq.gam.1, rsq.gam.2, rsq.gam.3)))
      
    } # close knot loop
  } # close timeseries loop
  
  # Label best models by timeseries
  model.out %>%
    group_by(TS) %>%
    mutate(sig = BH2(p_gam, alph = 0.05)$BHSig,
           p_gam = p_gam,
           BEST = ifelse(AIC == min(AIC), "Y", "N"),
           padj = p.adjust(p_gam, method = "fdr")) %>%
    filter(BEST == "Y") -> model.out
  
  return(model.out)
}

# Create function to calculate AR1 and CV/SD for 15 year windows for biology/SST timeseries
sum.fun <- function(dat, data.type, wind){
  
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
  }else{
    dat = dat
  }
  
  
  if(data.type != "Mature biomass"){
    for(ii in 1:length(unique(dat$TS))){
      
      dat %>%
        filter(TS == unique(dat$TS)[ii]) %>%
        na.omit() -> TS.dat
      
      # Specify sliding window width
      width = wind
      
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
                                           window = window,
                                           width =wind))
      
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
      width = wind
      
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
                                           window = window,
                                           width = wind))
      
    }
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
      detrend.dat <- loess(log.recruitment ~ Year, dat, span = 0.25, degree = 1)
      
      # Extract residuals
      detrend.dat.resid <- data.frame(TS = rep(names[ii], length(unique(dat$Year))),
                                      Year = dat$Year, log.recruitment = detrend.dat$residuals)
      
     
      
      # Calculate AR1 and CV/SD on detrended data
      wind %>%
        purrr::map(~sum.fun(detrend.dat.resid, data.type, .x)) -> out
      
      rbind(summary, as.data.frame(out)) -> summary
      
    }
  } else{
    
    dat <- TS.dat
    
    # Detrend data
    detrend.dat <- loess(mean.sst ~ Year, dat, span = 0.25, degree = 1)
    
    # Extract residuals
    detrend.dat.resid <- data.frame(TS = rep("SST", length(unique(dat$Year))),
                                    Year = dat$Year, mean.sst = detrend.dat$residuals)
    
    # Calculate AR1 and CV/SD on detrended data on different window lengths
    wind %>%
    purrr::map(~sum.fun(detrend.dat.resid, data.type, .x)) -> out
    
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
      
      # fit gams, record pval and AIC
      # gam.ar1 <- gam(ar1 ~ s(window, k = knts[kk]), 
      #                data= TS.dat) # removed ", correlation = corAR1()"
      
      gam.ar1 <- gamm(ar1 ~ s(window, k = knts[kk]), 
                    data= TS.dat, correlation = corAR1())
      
      # gam.var <- gam(var.val ~ s(window, k = knts[kk]), 
      #                data= TS.dat) # removed ", correlation = corAR1()"
      
      gam.var <- gamm(var.val ~ s(window, k = knts[kk]), 
                      data= TS.dat, correlation = corAR1())
      
      
      p.gam.ar1 <- summary(gam.ar1$gam)$s.table[,4]
      p.gam.var <- summary(gam.var$gam)$s.table[,4]
      
      r.gam.ar1 <- summary(gam.ar1$gam)$r.sq
      r.gam.var <- summary(gam.var$gam)$r.sq
      
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
    group_by(TS, response) %>%
    mutate(sig = BH2(p_gam, alph = 0.05)$BHSig,
           p_gam = p_gam,
           BEST = ifelse(AIC == min(AIC), "Y", "N"),
           padj = p.adjust(p_gam, method = "fdr")) %>%
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
        dplyr::select(TS, ar1, cv, sd, window) -> TS.dat2
      
      
      if(unique(model.dat$response)[jj] == "var.val" & data.type == "Recruitment"){
        value = TS.dat2$cv
      } else if(unique(model.dat$response)[jj] == "var.val" & data.type == "sst"){
        value = TS.dat2$sd
      } else{
        value = TS.dat2$ar1
      }
      
      gam.best <- gam(value~ s(window, k = model.dat2$knots), 
                      data= TS.dat2)
      
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
                                               window = TS.dat2$window))
      
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

# 
# BHA <- function(p){
#   u <- sort(p, decreasing = TRUE)
#   mult <- length(p)/(length(p):1)
#   pNew <- cummin(u*mult)[order(order(-p))]
#   return(pNew)
# }
# 
# control.FDR <- function(p, FDR = 0.05, arr = c("p.orig", "none")){
#   BHP <- BH2(p, FDR)
#   p.adj <- BHA(p)
#   p.orig <- p
#   message(paste0("For a ", FDR
#                  , "-level false discovery rate, Threshold p-value = ", BHP$pCrit))
#   out <- data.frame( p.orig
#                      , significant_after_FDR = BHP$BHSig
#                      , p.adj)
#   if(arr == "none"){
#     return(out)
#   } 
#   else {return(out[order(get(arr)),])}
# }
# 
# BHAdjust(Hedenfalk$x, FDR = 0.05, arr = "p.orig") -> tt

# 
# pvalues<-c(0.01,0.001, 0.05, 0.20, 0.15, 0.15)
# ranks<-rank(pvalues, ties.method = "last")
# p_m_over_k<-pvalues*length(pvalues)/ranks
# 
# for (r in length(pvalues):1) {
#   print(p_m_over_k[ranks>=r])
# }
# 
# pvalues_adj<-c()
# 
# for (i in 1:length(pvalues)) {
#   
#   # find the rank
#   tmp_rank<-ranks[2]
#   
#   # get all the p_m_over_k that are greater or equal to this rank
#   # and get the min value
#   pvalues_adj<-c(pvalues_adj, min(1,min(p_m_over_k[ranks>=tmp_rank])))
# }
