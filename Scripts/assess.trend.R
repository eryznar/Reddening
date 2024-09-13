#GAMs on detrendend sst and biological ts

#Purpose is to assess the following questions: 
# 1) Has the climate gotten redder? 
# 2) Has the biology gotten redder?
# To do so, this script will test for trends in sst AR1 and SD and recruitment AR1 and CV
# by fitting GAMs on detrended data

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")

source("./Scripts/functions.R")

### Run function to detrend data, calculate ar1, SD, avg, and cv for 15 year windows -----------------------------------------------------------------------------------------
  sum.out <- data.frame()
  
  # Recruitment
  trend.fun(bsai.r0, "Recruitment") -> bsai.r0.out
  
  trend.fun(goa.r0, "Recruitment") -> goa.r0.out
  
  # SST
  trend.fun(ebs.sst, "SST") -> ebs.sst.out
  
  trend.fun(goa.sst, "SST") -> goa.sst.out
  
### Test question 1: GAMs to predict sst AR1 and SST with time -----------------------------------------------------------
  # AR1 and SD for SST were calculated using 15-year windows
  
  # Make function to streamline model fitting
  assess.trend <- function(dat, data.type){
    
    for(ii in 1:length(unique(dat$TS))){
      # filter data by TS
      dat %>%
        filter(TS == unique(dat$TS)[ii]) %>%
        group_by(TS) %>%
        mutate(N = n()) %>%
        ungroup()-> TS.dat
      
      knts <- c(3, 4)
      
      # specify variance value by data type
      if(data.type == "Recruitment"){
        var.val = dat$cv
      } else{
        var.val = dat$sd
      }
      
      # fit gams, iterating through different knots
      for(kk in 1:length(knts)){
        
        # fit gams, record pval and AIC
        gam.ar1 <- gam(ar1 ~ s(window, k = knts[kk]), 
                      data= dat) # removed ", correlation = corAR1()"
        
        gam.var <- gam(var.val ~ s(window, k = knts[kk]), 
                        data= dat) # removed ", correlation = corAR1()"
        
        
        p.gam.ar1 <- summary(gam.ar1)$s.table[,4]
        p.gam.var <- summary(gam.var)$s.table[,4]
        
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
                                                 AIC = c(AIC.gam.ar1, AIC.gam.var)))
        
      } # close knot loop
    } # close timeseries loop
    
    # Label best models by timeseries
    model.out %>%
      group_by(TS, response) %>%
      mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) %>%
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
                                                 pred = predict(gam.best, se.fit =TRUE)$fit,
                                                 pred.CI = 1.96*(predict(gam.best, se.fit =TRUE)$se.fit),
                                                 k = model.dat2$knots,
                                                 window = TS.dat2$window))
        
      } # close response loop
    } # close timeseries loop
    
    return(pred.vals)
  } # close function
  
  # Run functions for SST
  model.out <- data.frame()
  
  assess.trend(ebs.sst.out, "sst") -> ebs.sst.best.mods
  
  model.out <- data.frame()
  
  assess.trend(goa.sst.out, "sst") -> goa.sst.best.mods
  
  pred.vals <- data.frame()
  
  model.predict(ebs.sst.best.mods, ebs.sst.out, "sst") -> ebs.pred.out
  
  pred.vals <- data.frame()
  
  model.predict(goa.sst.best.mods, goa.sst.out, "sst") -> goa.pred.out
  
  
  # bind
  rbind(ebs.pred.out %>% mutate(region = "Eastern Bering Sea"),
        goa.pred.out %>% mutate(region = "Gulf of Alaska")) -> plot.dat.sst
  
 
  # Plot
  plot.dat.sst %>%
    dplyr::select(TS, response, region, k)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") -> lab.dat
  
  labs <- paste0(lab.dat$region, " \n(k=", lab.dat$k, ")")
  names(labs) <- c("Eastern Bering Sea", "Gulf of Alaska")
  
  plot.dat.sst %>%
    filter(response == "ar1") -> plot.dat.sst2
  
  # Plot sst AR1 with time
  ggplot()+
    geom_ribbon(plot.dat.sst2, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI,
                              fill = region), alpha = 0.25)+
    geom_point(plot.dat.sst2, mapping=aes(x = window, y = pred, color = region))+
    geom_line(plot.dat.sst2, mapping = aes(x = window, y = pred, color = region), size = 1.25)+
    scale_color_manual(values = c("#6A6DB7", "#A34242"))+
    scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
    facet_wrap(~region, scales = "free", labeller = labeller(region = labs),
               ncol = 1)+
    theme_bw()+
    ggtitle("SST AR1 with time")+
    ylab("AR1")+
    xlab("Window")+
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16)) -> sst.ar1.window.plot
  
  ggsave(plot = sst.ar1.window.plot, "./Figures/sstAR1.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # SD
  plot.dat.sst %>%
    dplyr::select(TS, response, region, k)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "var.val") -> lab.dat
  
  labs <- paste0(lab.dat$region, " \n(k=", lab.dat$k, ")")
  names(labs) <- c("Eastern Bering Sea", "Gulf of Alaska")
  
  plot.dat.sst %>%
    filter(response == "var.val") -> plot.dat.sst2
  
  # Plot sst SD with time
  ggplot()+
    geom_ribbon(plot.dat.sst2, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI,
                              fill = region), alpha = 0.25)+
    geom_point(plot.dat.sst2, mapping=aes(x = window, y = pred, color = region))+
    geom_line(plot.dat.sst2, mapping = aes(x = window, y = pred, color = region), size = 1.25)+
    scale_color_manual(values = c("#6A6DB7", "#A34242"))+
    scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
    facet_wrap(~region, scales = "free", labeller = labeller(region = labs),
               ncol = 1)+
    theme_bw()+
    ggtitle("SST SD with time")+
    ylab("SD")+
    xlab("Window")+
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16)) -> sst.SD.window.plot
  
  ggsave(plot = sst.SD.window.plot, "./Figures/sstSD.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  
  # Predict and plot with best models GOA
  model.predict <- data.frame()
  
  for(ii in 1:length(unique(goa.r0.best$TS))){
    goa.r0.best %>%
      filter(TS == unique(goa.r0.best$TS)[ii]) -> model.dat
    
    goa.r0.sst %>%
      filter(TS == unique(goa.r0.best$TS)[ii]) %>%
      dplyr::select(TS, log.recruitment, grep(model.dat$sst, names(.))) -> TS.dat
    
    gamm.best <- gamm(log.recruitment ~ s(TS.dat[,3], k = model.dat$knots), 
                      data= TS.dat, correlation = corAR1())
    
    model.predict <- rbind(model.predict, data.frame(TS = unique(goa.r0.best$TS)[ii],
                                                     log.recruitment = TS.dat$log.recruitment,
                                                     pred.r0 = predict(gamm.best, se.fit =TRUE)$fit,
                                                     pred.CI = 1.96*(predict(gamm.best, se.fit =TRUE)$se.fit),
                                                     sst = TS.dat[,3],
                                                     sst.smooth = model.dat$sst,
                                                     k = model.dat$knots))
    
    saveRDS(gamm.best, paste0("./Output/", unique(goa.r0.best$TS)[ii], ".gamm.rda"))
    
  }
  
  model.predict %>%
    dplyr::select(TS, sst.smooth, k)%>%
    distinct() %>%
    arrange(., TS) %>%
    mutate(sst.smooth = case_when((sst.smooth == "unsmoothed") ~ "1-year",
                                  (sst.smooth == "twoyear") ~ "2-year",
                                  (sst.smooth == "threeyear") ~ "3-year")) -> lab.dat
  
  labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, ", sst=", lab.dat$sst.smooth, ")")
  names(labs) <- names(r0.labs.goa)
  
  # Plot BSAI r0/SST
  ggplot()+
    geom_ribbon(model.predict, mapping = aes(x = sst, ymin = pred.r0 - pred.CI, ymax= pred.r0+pred.CI),
                fill = "grey", alpha = 0.75)+
    geom_point(model.predict, mapping=aes(x = sst, y = log.recruitment))+
    geom_line(model.predict, mapping = aes(x = sst, y = pred.r0), size = 1.25, color = "#A34242")+
    facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("GOA groundfish/crab recruitment and SST")+
    ylab("log(millions of recruits)")+
    xlab("°C")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> goa.r0.sst.plot
  
  ggsave(plot = goa.r0.sst.plot, "./Figures/goa.r0.sst.plot2.png", width = 8.5, height = 11, units = "in")
  
