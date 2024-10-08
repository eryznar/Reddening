#GAMs on detrended sst and biological ts

#Purpose is to assess the following questions: 
# 1) Has the climate gotten redder?
# 2) Is climate related to the biology (recruitment)?
# 3) Has the biology gotten redder?
# 4) Does window length influence AR1 and variance?
# 5) Is AL variability related to SST reddening?
# 6) Are recruitment-SST relationships actually valid or spurious?
# 7) What are the impacts of increasing ar1?

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/load.libs.functions.R")
source("./Scripts/ts.processing.R")

### SET UNIVERSAL PARAMETERS --------------------------------------------------------
width = 15 # rolling window length

### Test question 1: Has the climate gotten redder? -----------------------------------------------------------
 # APPROACH: fit GAMs between SST AR1/SD and time

 # Load (trying anomaly data)  
  ebs.sst <- read.csv("./Output/SST.winter.anom.ebs.csv") %>%
    filter(Year %in% 1948:2024)
  goa.sst <- read.csv("./Output/SST.winter.anom.goa.csv")%>%
    filter(Year %in% 1948:2024)
  
  rbind(ebs.sst %>% mutate(region = "Eastern Bering Sea"),
        goa.sst %>% mutate(region = "Gulf of Alaska")) -> sst

  # Detrend data
  trend.fun(ebs.sst, "SST", width) -> ebs.sst.out
  
  trend.fun(goa.sst, "SST", width) -> goa.sst.out
  
  ar1var.EBS.sst <- as.data.frame(ebs.sst.out)
  ar1var.goa.sst <- as.data.frame(goa.sst.out)
  
  # Fit and select best models
  model.out <- data.frame()
  
  assess.trend(ar1var.EBS.sst, "sst") -> ebs.sst.best.mods
  
  model.out <- data.frame()
  
  assess.trend(ar1var.goa.sst, "sst") -> goa.sst.best.mods
  
  # Predict with best models
  pred.vals <- data.frame()
  
  model.predict(ebs.sst.best.mods, ar1var.EBS.sst %>% na.omit(), "sst") -> ebs.pred.out
  
  pred.vals <- data.frame()
  
  model.predict(goa.sst.best.mods, ar1var.goa.sst %>% na.omit(), "sst") -> goa.pred.out
  
  # Bind
  rbind(ebs.pred.out %>% mutate(region = "Eastern Bering Sea"),
        goa.pred.out %>% mutate(region = "Gulf of Alaska")) -> plot.dat.sst
  
  # Plot sst AR1 with time
   plot.dat.sst %>%
    dplyr::select(TS, response, region, k, rsq, p_gam, sig)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") %>%
     mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                             (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                             (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                             TRUE ~ paste0("p=", round(p_gam, 2))),
            psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                             TRUE ~ plab)) -> lab.dat
   
  labs <- paste0(lab.dat$region, " \n(k=", lab.dat$k, ", R2=", round(lab.dat$rsq, 2), ", ", lab.dat$psig, ")")
  names(labs) <- c("Eastern Bering Sea", "Gulf of Alaska")
  
  plot.dat.sst %>%
    filter(response == "ar1") -> plot.dat.sst2
  
  ggplot()+
    geom_ribbon(plot.dat.sst2, 
                mapping = aes(x = year, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    geom_point(plot.dat.sst2, mapping=aes(x =  year, y = observed), color = "black")+
    geom_line(plot.dat.sst2, mapping = aes(x =  year, y = pred, color = region), size = 1.25)+
    scale_color_manual(values = c("#6A6DB7", "#A34242"))+
    scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
    facet_wrap(~region, scales = "free_y", labeller = labeller(region = labs),
               ncol = 1)+
    theme_bw()+
    ggtitle("SST AR1 with time")+
    geom_vline(xintercept = 1988.5, linetype = "dashed")+
    ylab("AR1")+
    xlab("Year")+
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16)) -> sst.ar1.window.plot
  
  ggsave(plot = sst.ar1.window.plot, "./Figures/sst.AR1.x.time.GAM2.png", width = 8.5, height = 11, units = "in")
  
  # Plot sst SD with time
  plot.dat.sst %>%
    dplyr::select(TS, response, region, k, rsq, p_gam, sig)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "var.val") %>%
    mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                            (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                            (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                            TRUE ~ paste0("p=", round(p_gam, 2))),
           psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                            TRUE ~ plab)) -> lab.dat
   
  labs <- paste0(lab.dat$region, " \n(k=", lab.dat$k, ", R2=", round(lab.dat$rsq, 2), ", ", lab.dat$psig, ")")
  names(labs) <- c("Eastern Bering Sea", "Gulf of Alaska")
  
  plot.dat.sst %>%
    filter(response == "var.val") -> plot.dat.sst2
  
  ggplot()+
    geom_ribbon(plot.dat.sst2, 
                mapping = aes(x = year, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    geom_point(plot.dat.sst2, mapping=aes(x =  year, y = observed), color = "black")+
    geom_line(plot.dat.sst2, mapping = aes(x =  year, y = pred, color = region), size = 1.25)+
    scale_color_manual(values = c("#6A6DB7", "#A34242"))+
    scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
    facet_wrap(~region, scales = "free_y", labeller = labeller(region = labs),
               ncol = 1)+
    theme_bw()+
    ggtitle("SST SD with time")+
    geom_vline(xintercept = 1988.5, linetype = "dashed")+
    ylab("SD")+
    xlab("Year")+
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16)) -> sst.SD.window.plot
  
  ggsave(plot = sst.SD.window.plot, "./Figures/sst.SD.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  
### Test question 2: Is climate related to the biology (recruitment)? ----------------------------------------------
  # APPROACH: fit GAMS between BSAI and GOA recruitment (not detrended) and SST (not detrended)
  
  # Load (trying anomaly data)  
  ebs.sst <- read.csv("./Output/SST.winter.anom.ebs.csv") %>%
    filter(Year %in% 1948:2024)
  goa.sst <- read.csv("./Output/SST.winter.anom.goa.csv")%>%
    filter(Year %in% 1948:2024)
  
  rbind(ebs.sst %>% mutate(region = "Eastern Bering Sea"),
        goa.sst %>% mutate(region = "Gulf of Alaska")) -> sst
  
  # Evaluate candidate models based on unsmoothed, 2-year, and 3-year rolling mean SST
  sst.2.ebs <- rollapply((sst %>% filter(region == "Eastern Bering Sea"))$mean.sst, 2, mean, na.rm = T, fill = NA) #2-year
  sst.3.ebs <- rollapply((sst %>% filter(region == "Eastern Bering Sea"))$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
  
  sst.2.goa <- rollapply((sst %>% filter(region == "Gulf of Alaska"))$mean.sst, 2, mean, na.rm = T, fill = NA) #2-year
  sst.3.goa <- rollapply((sst %>% filter(region == "Gulf of Alaska"))$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
  
  sst.rollmeans <- data.frame(sst, twoyear = c(sst.2.ebs, sst.2.goa), 
                              threeyear = c(sst.3.ebs, sst.3.goa)) %>%
    rename(unsmoothed = mean.sst)

  # EBS 
      # Add rollmean sst to  r0
      bsai.r0 %>%
        right_join(., expand.grid(Year = min(bsai.r0$Year):max(bsai.r0$Year))) %>%
        right_join(., sst.rollmeans %>% filter(region == "Eastern Bering Sea") %>%
                     rename(Lagged.Year = Year), by = c("Lagged.Year")) %>%
        na.omit() %>%
        dplyr::select(!region) %>%
        filter(Year > 1987)-> bsai.r0.sst
      
      # Fit models for BSAI
      model.out <- data.frame()
      
      fit.r0.SST.models(bsai.r0.sst) -> bsai.r0.best
      
      # Calculate the proportion of significant results after controlling for FDR
      sum(bsai.r0.best$sig == TRUE)/nrow(bsai.r0.best) -> ebs.true.sig
      
      # Predict and plot with best models BSAI
      model.prediction <- data.frame()
      
      for(ii in 1:length(unique(bsai.r0.best$TS))){
        bsai.r0.best %>%
          filter(TS == unique(bsai.r0.best$TS)[ii]) -> model.dat
        
        bsai.r0.sst %>%
          filter(TS == unique(bsai.r0.best$TS)[ii]) %>%
          dplyr::select(TS, log.recruitment, grep(model.dat$sst, names(.))) -> TS.dat
        
        gam.best <- gam(log.recruitment ~ s(TS.dat[,3], k = model.dat$knots), 
                          data= TS.dat, correlation = corAR1())
        
        model.prediction <- rbind(model.prediction, data.frame(TS = unique(bsai.r0.best$TS)[ii],
                                                         log.recruitment = TS.dat$log.recruitment,
                                                         pred.r0 = predict(gam.best, se.fit =TRUE)$fit,
                                                         pred.CI = 1.96*(predict(gam.best, se.fit =TRUE)$se.fit),
                                                         sst = TS.dat[,3],
                                                         sst.smooth = model.dat$sst,
                                                         k = model.dat$knots,
                                                         p_gam = model.dat$p_gam,
                                                         #p_lme = model.dat$p_lme,
                                                         rsq = model.dat$rsq,
                                                         sig = model.dat$sig))
        
      }
      
      # Specify plotting labels
      model.prediction %>%
        dplyr::select(TS, sst.smooth, k, p_gam, rsq, sig)%>%
        distinct() %>%
        arrange(., TS) %>%
        mutate(sst.smooth = case_when((sst.smooth == "unsmoothed") ~ "1-year",
                                      (sst.smooth == "twoyear") ~ "2-year",
                                      (sst.smooth == "threeyear") ~ "3-year")) %>%
        mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                                (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                                (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                                TRUE ~ paste0("p=", round(p_gam, 2))),
               psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                                TRUE ~ plab)) -> lab.dat
      
      labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, ", sst=", lab.dat$sst.smooth, ", R2=", lab.dat$rsq, ", ",
                     lab.dat$psig, ")")
      names(labs) <- names(r0.labs.bsai)
      
      # Plot BSAI r0/SST
      ggplot()+
        geom_ribbon(model.prediction, mapping = aes(x = sst, ymin = pred.r0 - pred.CI, ymax= pred.r0+pred.CI),
                    fill = "grey", alpha = 0.75)+
        geom_point(model.prediction, mapping=aes(x = sst, y = log.recruitment))+
        geom_line(model.prediction, mapping = aes(x = sst, y = pred.r0), size = 1.25, color = "#6A6DB7")+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
                   ncol = 4)+
        theme_bw()+
        ggtitle("BSAI groundfish/crab recruitment and SST")+
        ylab("log(millions of recruits)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> bsai.r0.sst.plot
      
      ggsave(plot = bsai.r0.sst.plot, "./Figures/bsai.r0.sst.plot2.png", width = 11, height = 8.5, units = "in")
      
      
  # GOA 
      # Add rollmean sst to  r0
      goa.r0 %>%
        right_join(., expand.grid(Year = min(goa.r0$Year):max(goa.r0$Year))) %>%
        right_join(., sst.rollmeans %>% filter(region == "Gulf of Alaska") %>%
                     rename(Lagged.Year = Year), by = c("Lagged.Year")) %>%
        na.omit() %>%
        dplyr::select(!region) %>%
        filter(Year > 1987)-> goa.r0.sst
      
      # Fit models for GOA
      model.out <- data.frame()
      
      fit.r0.SST.models(goa.r0.sst) -> goa.r0.best
      
      # Calculate the proportion of significant results after controlling for FDR
      sum(goa.r0.best$sig == TRUE)/nrow(goa.r0.best) -> goa.true.sig
      
      # Predict and plot with best models BSAI
      model.prediction <- data.frame()
      
      for(ii in 1:length(unique(goa.r0.best$TS))){
        goa.r0.best %>%
          filter(TS == unique(goa.r0.best$TS)[ii]) -> model.dat
        
        goa.r0.sst %>%
          filter(TS == unique(goa.r0.best$TS)[ii]) %>%
          dplyr::select(TS, log.recruitment, grep(model.dat$sst, names(.))) -> TS.dat
        
        gam.best <- gam(log.recruitment ~ s(TS.dat[,3], k = model.dat$knots), 
                          data= TS.dat, correlation = corAR1())
        
        model.prediction <- rbind(model.prediction, data.frame(TS = unique(goa.r0.best$TS)[ii],
                                                         log.recruitment = TS.dat$log.recruitment,
                                                         pred.r0 = predict(gam.best, se.fit =TRUE)$fit,
                                                         pred.CI = 1.96*(predict(gam.best, se.fit =TRUE)$se.fit),
                                                         sst = TS.dat[,3],
                                                         sst.smooth = model.dat$sst,
                                                         k = model.dat$knots,
                                                         p_gam = model.dat$p_gam,
                                                         rsq = model.dat$rsq,
                                                         sig = model.dat$sig))
        
      }
      
      # Specify plotting labels
      model.prediction %>%
        dplyr::select(TS, sst.smooth, k, p_gam, rsq, sig)%>%
        distinct() %>%
        arrange(., TS) %>%
        mutate(sst.smooth = case_when((sst.smooth == "unsmoothed") ~ "1-year",
                                      (sst.smooth == "twoyear") ~ "2-year",
                                      (sst.smooth == "threeyear") ~ "3-year")) %>%
        mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                                (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                                (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                                TRUE ~ paste0("p=", round(p_gam, 2))),
               psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                                TRUE ~ plab)) -> lab.dat
      
      labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, ", sst=", lab.dat$sst.smooth, ", R2=", lab.dat$rsq, ", ",
                     lab.dat$psig, ")")
      names(labs) <- names(r0.labs.goa)
      
      # Plot GOA r0/SST
      ggplot()+
        geom_ribbon(model.prediction, mapping = aes(x = sst, ymin = pred.r0 - pred.CI, ymax= pred.r0+pred.CI),
                    fill = "grey", alpha = 0.75)+
        geom_point(model.prediction, mapping=aes(x = sst, y = log.recruitment))+
        geom_line(model.prediction, mapping = aes(x = sst, y = pred.r0), size = 1.25, color = "#A34242")+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
                   ncol = 4)+
        theme_bw()+
        ggtitle("GOA groundfish recruitment and SST")+
        ylab("log(millions of recruits)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> goa.r0.sst.plot
      
      ggsave(plot = goa.r0.sst.plot, "./Figures/goa.r0.sst.plot2.png", width = 11.5, height = 8.5, units = "in")
      
      
  
### Test question 3: Has the biology gotten redder?----------------------------------------------------
  # APPROACH: fitting GAMs to predict BSAI and GOA recruitment AR1 and CV with time
      
   # EBS
      # Detrend r0 data
      sum.out <- data.frame()
      
      trend.fun(bsai.r0, "Recruitment", width) -> ebs.r0.out
      
      ar1var.EBS.r0 <- as.data.frame(ebs.r0.out)
      
      # Fit and select best models
      model.out <- data.frame()
      
      assess.trend(ar1var.EBS.r0, "Recruitment") -> ebs.r0.best.mods
      
      # Predict with best models
      pred.vals <- data.frame()
      
      model.predict(ebs.r0.best.mods, ar1var.EBS.r0 %>% na.omit(), "Recruitment") -> ebs.r0.pred.out
      
     # Plot r0 AR1 with time
      ebs.r0.pred.out %>%
        dplyr::select(TS, response, k, rsq, p_gam, sig)%>%
        distinct() %>%
        arrange(., TS) %>%
        filter(response == "ar1") %>%
        mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                                (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                                (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                                TRUE ~ paste0("p=", round(p_gam, 2))),
               psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                                TRUE ~ plab)) -> lab.dat
                 
           
      labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ", ", lab.dat$psig, ")")
      names(labs) <- names(r0.labs.bsai)
      
      ebs.r0.pred.out %>%
        filter(response == "ar1") -> plot.dat
      
      ggplot()+
        geom_ribbon(plot.dat, 
                    mapping = aes(x = year, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
        geom_point(plot.dat, mapping=aes(x = year, y = observed), color = "black")+
        geom_line(plot.dat, mapping = aes(x = year, y = pred), color =  "#6A6DB7", size = 1.25)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
                   ncol = 4)+
        theme_bw()+
        ggtitle("EBS recruitment AR1 with time")+
        ylab("AR1")+
        xlab("Year")+
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> EBS.r0.ar1.window.plot
      
      ggsave(plot = EBS.r0.ar1.window.plot, "./Figures/EBS.r0.AR1.x.time.GAM.png", width = 11, height = 8.5, units = "in")
      
      # Plot r0 CV with time
      ebs.r0.pred.out %>%
        dplyr::select(TS, response, k, rsq, p_gam, sig)%>%
        distinct() %>%
        arrange(., TS) %>%
        filter(response == "var.val") %>%
        mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                                (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                                (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                                TRUE ~ paste0("p=", round(p_gam, 2))),
               psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                                TRUE ~ plab)) -> lab.dat
      
      
      labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ", ", lab.dat$psig, ")")
      names(labs) <- names(r0.labs.bsai)
      
      ebs.r0.pred.out %>%
        filter(response == "var.val") -> plot.dat
      
      
      ggplot()+
        geom_ribbon(plot.dat, 
                    mapping = aes(x = year, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
        geom_point(plot.dat, mapping=aes(x = year, y = observed), color = "black")+
        geom_line(plot.dat, mapping = aes(x = year, y = pred), color =  "#6A6DB7", size = 1.25)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
                   ncol = 4)+
        theme_bw()+
        ggtitle("EBS recruitment CV with time")+
        ylab("CV")+
        xlab("Year")+
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> EBS.r0.CV.window.plot
      
      ggsave(plot = EBS.r0.CV.window.plot, "./Figures/EBS.r0.CV.x.time.GAM.png", width = 11, height = 8.5, units = "in")
      
  # GOA 
      # Detrend data
      sum.out <- data.frame()
      
      trend.fun(goa.r0, "Recruitment", width) -> goa.r0.out
      
      ar1var.goa.r0 <- as.data.frame(goa.r0.out)
      
      # Fit and select best models
      model.out <- data.frame()
      
      assess.trend(ar1var.goa.r0, "Recruitment") -> goa.r0.best.mods
      
      # Predict with best models
      pred.vals <- data.frame()
      
      model.predict(goa.r0.best.mods, ar1var.goa.r0, "Recruitment") -> goa.r0.pred.out
      
      
      # Plot r0 AR1 with time
      goa.r0.pred.out %>%
        dplyr::select(TS, response, k, rsq, p_gam, sig)%>%
        distinct() %>%
        arrange(., TS) %>%
        filter(response == "ar1") %>%
        mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                                (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                                (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                                TRUE ~ paste0("p=", round(p_gam, 2))),
               psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                                TRUE ~ plab)) -> lab.dat
      
      
      labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ", ", lab.dat$psig, ")")
      names(labs) <- names(r0.labs.goa)
      
      goa.r0.pred.out %>%
        filter(response == "ar1") -> plot.dat
      
      ggplot()+
        geom_ribbon(plot.dat, 
                    mapping = aes(x = year, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
        geom_point(plot.dat, mapping=aes(x = year, y = observed), color = "black")+
        geom_line(plot.dat, mapping = aes(x = year, y = pred), color =  "#A34242", size = 1.25)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
                   ncol = 4)+
        theme_bw()+
        ggtitle("GOA recruitment AR1 with time")+
        ylab("AR1")+
        xlab("Year")+
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> GOA.r0.ar1.window.plot
      
      ggsave(plot = GOA.r0.ar1.window.plot, "./Figures/GOA.r0.AR1.x.time.GAM.png", width = 11, height = 8.5, units = "in")
    
      # Plot r0 CV with time
      goa.r0.pred.out %>%
        dplyr::select(TS, response, k, rsq, p_gam, sig)%>%
        distinct() %>%
        arrange(., TS) %>%
        filter(response == "var.val") %>%
        mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                                (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                                (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                                TRUE ~ paste0("p=", round(p_gam, 2))),
               psig = case_when((sig == TRUE) ~ paste0(plab, "*"),
                                TRUE ~ plab)) -> lab.dat
      
      
      labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ", ", lab.dat$psig, ")")
      names(labs) <- names(r0.labs.goa)
      
      goa.r0.pred.out %>%
        filter(response == "var.val") -> plot.dat
      
      ggplot()+
        geom_ribbon(plot.dat, 
                    mapping = aes(x = year, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
        geom_point(plot.dat, mapping=aes(x = year, y = observed), color = "black")+
        geom_line(plot.dat, mapping = aes(x = year, y = pred), color =  "#A34242", size = 1.25)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
                   ncol = 4)+
        theme_bw()+
        ggtitle("GOA recruitment CV with time")+
        ylab("CV")+
        xlab("Year")+
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> GOA.r0.CV.window.plot
      
      ggsave(plot = GOA.r0.CV.window.plot, "./Figures/GOA.r0.CV.x.time.GAM.png", width = 11, height = 8.5, units = "in")
      
      
### Test question 4: Does window length influence AR1 and variance? ---------------------------------------------------
  # SST 
 # Detrend data
 c(10, 15, 20, 25) %>%
   purrr::map_df(~trend.fun(ebs.sst, "SST", .x)) -> ebs.sst.out
 
 c(10, 15, 20, 25) %>%
   purrr::map_df(~trend.fun(goa.sst, "SST", .x)) -> goa.sst.out
 
 
 rbind(ebs.sst.out %>% mutate(region = "Eastern Bering Sea"), 
       goa.sst.out %>% mutate(region = "Gulf of Alaska")) -> plot.dat
 
 ggplot()+
   geom_line(plot.dat, mapping = aes(year, ar1, color = as.factor(width)), size = 1.5, group = 1) +
   facet_wrap(~region, scales = "free_y", ncol = 1)+
   theme_bw()+
   ggtitle("Detrended SST AR1 at different window lengths")+
   scale_color_manual(name = "Window length", values = c("#8B0069", "#A96C00", "#75C165","#B0F4FA"))+
   ylab("AR1")+
   xlab("Year")+
   theme(axis.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         strip.text = element_text(size = 16),
         legend.text = element_text(size = 14),
         title = element_text(size = 16)) -> sst.AR1.windowlength
 
 
 ggsave(plot = sst.AR1.windowlength, "./Figures/sst.AR1.windowlength.png", width = 10, height = 11, units = "in")
 
 ggplot()+
   geom_line(plot.dat, mapping = aes(window, sd, color = as.factor(width)), size = 1.5) +
   facet_wrap(~region, scales = "free_y",
              ncol = 1)+
   theme_bw()+
   ggtitle("Detrended SST SD at different window lengths")+
   scale_color_manual(name = "Window length", values = c("#8B0069", "#A96C00", "#75C165","#B0F4FA"))+
   ylab("SD")+
   xlab("Window")+
   theme(axis.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         strip.text = element_text(size = 16),
         legend.text = element_text(size = 14),
         title = element_text(size = 16)) -> sst.SD.windowlength
 
 ggsave(plot = sst.SD.windowlength, "./Figures/sst.SD.windowlength.png", width = 10, height = 11, units = "in")
 
 
### Test question 5: Is AL variability related to SST reddening? ----------------------------------------------------
 # APPROACH: fit models between sea level pressure SD/AR1 and SST AR1
 
 # Calculate SLP AR1 and SD 
     slp <- read.csv("./Output/SLP.winter.anom.csv") %>%
       #filter(year %in% ebs.sst$Year) %>%
       dplyr::select(!X) %>%
       rename(Year = year)
     
     # Run function
     trend.fun(slp, "SLP", width) -> out
     
     ar1var.slp = as.data.frame(out)
    
     # Plot AR1
     ggplot(ar1var.slp,  mapping=aes(x = year, y = ar1))+
       geom_point()+
       geom_line()+
       theme_bw()+
       ggtitle("SLP AR1")+
       #scale_x_continuous(limits = c(min(ar1var.slp$year), max(ar1var.slp$year)), breaks = seq(min(ar1var.slp$year), max(ar1var.slp$year), by = 10))+
       ylab("AR1")+
       xlab("Year")+
       geom_vline(xintercept = 1988.5, linetype = "dashed")+
       theme(legend.position = "none",
             axis.text = element_text(size = 14),
             axis.title = element_text(size = 16),
             strip.text = element_text(size = 16),
             legend.text = element_text(size = 14),
             title = element_text(size = 16)) -> slp.ar1.window.plot
     
     ggsave(plot = slp.ar1.window.plot, "./Figures/slp.ar1.window.plot.png", width = 11, height = 8.5, units = "in")
     
     
     # Plot SD
     ggplot(ar1var.slp,  mapping=aes(x = year, y = sd))+
       geom_point()+
       geom_line()+
       theme_bw()+
       ggtitle("SLP SD")+
       #scale_x_continuous(limits = c(min(ar1var.slp$year), max(ar1var.slp$year)), breaks = seq(min(ar1var.slp$year), max(ar1var.slp$year), by = 10))+
       ylab("SD")+
       xlab("Year")+
       geom_vline(xintercept = 1988.5, linetype = "dashed")+
       theme(legend.position = "none",
             axis.text = element_text(size = 14),
             axis.title = element_text(size = 16),
             strip.text = element_text(size = 16),
             legend.text = element_text(size = 14),
             title = element_text(size = 16)) -> slp.sd.window.plot
     
     ggsave(plot = slp.sd.window.plot, "./Figures/slp.sd.window.plot.png", width = 11, height = 8.5, units = "in")
     
 # Fit gams between sst and SLP 
   ar1var.slp %>%
     na.omit() %>%
     rename(slp.ar1 = ar1, slp.sd = sd) %>%
     dplyr::select(year, slp.sd, slp.ar1) -> slp.model.dat
   
   ar1var.EBS.sst %>%
     na.omit() %>%
     rename(ebs.sst.ar1 = ar1, ebs.sst.sd = sd) %>%
     dplyr::select(year, ebs.sst.sd, ebs.sst.ar1) -> ebs.sst.model.dat
   
   ar1var.goa.sst %>%
     na.omit() %>%
     rename(goa.sst.ar1 = ar1, goa.sst.sd = sd) %>%
     dplyr::select(year, goa.sst.sd, goa.sst.ar1) -> goa.sst.model.dat
   
   right_join(goa.sst.model.dat, right_join(ebs.sst.model.dat, slp.model.dat)) -> model.dat
   
   k = 6 # set knots
   
     # Model: EBS sst AR1 x slp AR1 
     rr <- model.dat$ebs.sst.ar1
     pp <- model.dat$slp.ar1
     
     mod <- gam(rr ~ s(pp, k = k), 
                   data= model.dat, correlation = corAR1())
     
     p.val <- summary(mod)$s.table[,4]
     p.val <- case_when((p.val < 0.001) ~ "p<0.001",
                        (p.val < 0.01 & p.val >= 0.001) ~ "p<0.01",
                        (p.val <0.05 & p.val >= 0.01) ~ "p<0.05",
                        TRUE ~ paste0("p=", round(p.val, 2)))
     
     r.sq <- round(summary(mod)$r.sq, 2)
     
     pred <- predict(mod, se=T)
     pred.CI = 1.96*(predict(mod, se.fit =TRUE)$se.fit)
     
     plot.dat <- data.frame(rr = rr, pp = pp, pred = pred$fit, pred.se = pred$se, pred.CI = pred.CI)
     
     ggplot()+
       geom_ribbon(plot.dat, 
                   mapping = aes(x = pp, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
       geom_point(plot.dat, mapping=aes(x = pp, y = rr), color = "black")+
       geom_line(plot.dat, mapping = aes(x = pp, y = pred), color =  "#6A6DB7", size = 1.25)+
       theme_bw()+
       ggtitle(paste0("SLP AR1 vs. EBS SST AR1 (", p.val, ", R2 = ", r.sq, ")"))+
       ylab("SST anomaly AR1")+
       xlab("SLP anomaly AR1") +
       theme(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 10),
             legend.text = element_text(size = 12)) -> plot.1
     
     # Model: EBS sst AR1 x slp SD 
     rr <- model.dat$ebs.sst.ar1
     pp <- model.dat$slp.sd
     
     mod <- gam(rr ~ s(pp, k = k), 
                data= model.dat, correlation = corAR1())
     
     p.val <- summary(mod)$s.table[,4]
     p.val <- case_when((p.val < 0.001) ~ "p<0.001",
                        (p.val < 0.01 & p.val >= 0.001) ~ "p<0.01",
                        (p.val <0.05 & p.val >= 0.01) ~ "p<0.05",
                        TRUE ~ paste0("p=", round(p.val, 2)))
     r.sq <- round(summary(mod)$r.sq, 2)
     
     pred <- predict(mod, se=T)
     pred.CI = 1.96*(predict(mod, se.fit =TRUE)$se.fit)
     
     plot.dat <- data.frame(rr = rr, pp = pp, pred = pred$fit, pred.se = pred$se, pred.CI = pred.CI)
     
     ggplot()+
       geom_ribbon(plot.dat, 
                   mapping = aes(x = pp, ymin = pred - pred.se, ymax= pred +pred.se), fill = "grey", alpha = 0.5)+
       geom_point(plot.dat, mapping=aes(x = pp, y = rr), color = "black")+
       geom_line(plot.dat, mapping = aes(x = pp, y = pred), color =  "#6A6DB7", size = 1.25)+
       theme_bw()+
       ggtitle(paste0("SLP SD vs. EBS SST AR1 (", p.val, ", R2 = ", r.sq, ")"))+
       ylab("SST anomaly AR1")+
       xlab("SLP anomaly SD") +
       theme(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 10),
             legend.text = element_text(size = 12)) -> plot.2
 
     # Model: GOA sst AR1 x slp AR1 
     rr <- model.dat$goa.sst.ar1
     pp <- model.dat$slp.ar1
     
     mod <- gam(rr ~ s(pp, k = k), 
                data= model.dat, correlation = corAR1())
     
     
     p.val <- summary(mod)$s.table[,4]
     p.val <- case_when((p.val < 0.001) ~ "p<0.001",
                        (p.val < 0.01 & p.val >= 0.001) ~ "p<0.01",
                        (p.val <0.05 & p.val >= 0.01) ~ "p<0.05",
                        TRUE ~ paste0("p=", round(p.val, 2)))
     r.sq <- round(summary(mod)$r.sq, 2)
     
     pred <- predict(mod, se=T)
     pred.CI = 1.96*(predict(mod, se.fit =TRUE)$se.fit)
     
     plot.dat <- data.frame(rr = rr, pp = pp, pred = pred$fit, pred.se = pred$se, pred.CI = pred.CI)
     
     ggplot()+
       geom_ribbon(plot.dat, 
                   mapping = aes(x = pp, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
       geom_point(plot.dat, mapping=aes(x = pp, y = rr), color = "black")+
       geom_line(plot.dat, mapping = aes(x = pp, y = pred), color =  "#A34242", size = 1.25)+
       theme_bw()+
       ggtitle(paste0("SLP AR1 vs. GOA SST AR1 (", p.val, ", R2 = ", r.sq, ")"))+
       ylab("SST anomaly AR1")+
       xlab("SLP anomaly AR1") +
       theme(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 10),
             legend.text = element_text(size = 12)) -> plot.3
     
     # Model: GOA sst AR1 x slp SD 
     rr <- model.dat$goa.sst.ar1
     pp <- model.dat$slp.sd
     
     mod <- gam(rr ~ s(pp, k = k), 
                data= model.dat, correlation = corAR1())
     
     p.val <- summary(mod)$s.table[,4]
     p.val <- case_when((p.val < 0.001) ~ "p<0.001",
                        (p.val < 0.01 & p.val >= 0.001) ~ "p<0.01",
                        (p.val <0.05 & p.val >= 0.01) ~ "p<0.05",
                        TRUE ~ paste0("p=", round(p.val, 2)))
     r.sq <- round(summary(mod)$r.sq, 2)
     
     pred <- predict(mod, se=T)
     pred.CI = 1.96*(predict(mod, se.fit =TRUE)$se.fit)
     
     plot.dat <- data.frame(rr = rr, pp = pp, pred = pred$fit, pred.se = pred$se, pred.CI = pred.CI)
     
     ggplot()+
       geom_ribbon(plot.dat, 
                   mapping = aes(x = pp, ymin = pred - pred.se, ymax= pred +pred.se), fill = "grey", alpha = 0.5)+
       geom_point(plot.dat, mapping=aes(x = pp, y = rr), color = "black")+
       geom_line(plot.dat, mapping = aes(x = pp, y = pred), color =  "#A34242", size = 1.25)+
       theme_bw()+
       ggtitle(paste0("SLP SD vs. GOA SST AR1 (", p.val, ", R2 = ", r.sq, ")"))+
       ylab("SST anomaly AR1")+
       xlab("SLP anomaly SD") +
       theme(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 10),
             legend.text = element_text(size = 12)) -> plot.4
     
    plot_grid(plot.1, plot.2, plot.3, plot.4) -> pp
    
    ggsave(plot = pp, "./Figures/SLP.vs.SST.png", width = 11, height = 8.5, units = "in")

 # Calculate correlations through time 
      year <- data.frame(year = min(slp$Year):max(slp$Year))
    
      right_join(model.dat, year) %>%
        mutate(cor.SD.ebsAR1 = NA,
               cor.SD.goaAR1 = NA) %>%
        arrange(year) -> cor.dat
    
      for(ii in (width+1):(nrow(cor.dat)-width)){ # 132 month rolling windows!
        
        cor.dat$cor.SD.ebsAR1[ii] <- cor(cor.dat$slp.sd[(ii-width):(ii+width)], cor.dat$ebs.sst.ar1[(ii-width):(ii+width)])
        cor.dat$cor.SD.goaAR1[ii] <- cor(cor.dat$slp.sd[(ii-width):(ii+width)], cor.dat$goa.sst.ar1[(ii-width):(ii+width)])
        
      }
      
      ggplot()+
        geom_line(na.omit(cor.dat), mapping = aes(year, cor.SD.ebsAR1, color = "salmon"), linewidth = 1.5)+
        geom_line(na.omit(cor.dat), mapping = aes(year, cor.SD.goaAR1,color = "steelblue"), linewidth = 1.5)+
        geom_vline(xintercept = 1988.5, linetype = "dashed", linewidth = 1.5)+
        scale_color_manual(name = "", values = c("salmon", "steelblue"), labels = c("SLP SD x EBS SST AR1", "SLP SD x GOA SST AR1"))+
        theme_bw()+
        ylab("Pearson's correlation")+
        xlab("Year")+
        theme(legend.position = "bottom",
              legend.direction= "horizontal",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16),
              legend.text = element_text(size = 14),
              title = element_text(size = 16)) -> cor.plot
      
      ggsave(plot = cor.plot, "./Figures/SLP.vs.SST.corplot.png", width = 11, height = 8.5, units = "in")
        
 
### Test question 6: Are recruitment-SST relationships actually valid or spurious? --------------------------------------------
  # APPROACH: simulate SST timeseries with same SD and AR1 of timeseries in each region, fit the same
  # GAM models as in question 2 on simulated timeseries
  
  # Load data
  read.csv("./Output/SST.winter.anom.ebs.csv") -> sst.ebs
  read.csv("./Output/SST.winter.anom.goa.csv") -> sst.goa
  
  # Calculate 2- and 3-year running means
  # Load (trying anomaly data)  
  ebs.sst <- read.csv("./Output/SST.winter.anom.ebs.csv") %>%
    filter(Year %in% 1948:2024)
  goa.sst <- read.csv("./Output/SST.winter.anom.goa.csv")%>%
    filter(Year %in% 1948:2024)
  
  rbind(ebs.sst %>% mutate(region = "Eastern Bering Sea"),
        goa.sst %>% mutate(region = "Gulf of Alaska")) -> sst
  
  # Evaluate candidate models based on unsmoothed, 2-year, and 3-year rolling mean SST
  sst.2.ebs <- rollapply((sst %>% filter(region == "Eastern Bering Sea"))$mean.sst, 2, mean, na.rm = T, fill = NA) #2-year
  sst.3.ebs <- rollapply((sst %>% filter(region == "Eastern Bering Sea"))$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
  
  sst.2.goa <- rollapply((sst %>% filter(region == "Gulf of Alaska"))$mean.sst, 2, mean, na.rm = T, fill = NA) #2-year
  sst.3.goa <- rollapply((sst %>% filter(region == "Gulf of Alaska"))$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
  
  
  # Calculate sd and ar1 for each sst ts
  # EBS
  sd((sst %>% filter(region == "Eastern Bering Sea"))$mean.sst) -> ebs.sd.unsmoothed
  sd(na.omit(sst.2.ebs)) -> ebs.sd.twoyear
  sd(na.omit(sst.3.ebs)) -> ebs.sd.threeyear
  
  acf((sst %>% filter(region == "Eastern Bering Sea"))$mean.sst, lag.max = 1, plot = FALSE)$acf[2] -> ebs.ar1.unsmoothed
  acf(na.omit(sst.2.ebs), lag.max = 1, plot = FALSE)$acf[2] -> ebs.ar1.twoyear
  acf(na.omit(sst.3.ebs), lag.max = 1, plot = FALSE)$acf[2] -> ebs.ar1.threeyear
  
  ebs.pars <- data.frame(smooth = c("unsmoothed", "twoyear", "threeyear"), 
                         sd = c(ebs.sd.unsmoothed, ebs.sd.twoyear, ebs.sd.threeyear),
                         ar1 = c(ebs.ar1.unsmoothed, ebs.ar1.twoyear, ebs.ar1.threeyear))
  
  
  # GOA
  sd((sst %>% filter(region == "Gulf of Alaska"))$mean.sst) -> goa.sd.unsmoothed
  sd(na.omit(sst.2.goa)) -> goa.sd.twoyear
  sd(na.omit(sst.3.goa)) -> goa.sd.threeyear
  
  acf((sst %>% filter(region == "Gulf of Alaska"))$mean.sst, lag.max = 1, plot = FALSE)$acf[2] -> goa.ar1.unsmoothed
  acf(na.omit(sst.2.goa), lag.max = 1, plot = FALSE)$acf[2] -> goa.ar1.twoyear
  acf(na.omit(sst.3.goa), lag.max = 1, plot = FALSE)$acf[2] -> goa.ar1.threeyear
  
  goa.pars <- data.frame(smooth = c("unsmoothed", "twoyear", "threeyear"), 
                         sd = c(goa.sd.unsmoothed, goa.sd.twoyear, goa.sd.threeyear),
                         ar1 = c(goa.ar1.unsmoothed, goa.ar1.twoyear, goa.ar1.threeyear))
  
  
  # Create function to generate random timeseries based on sst pars in each region and fit recruitment models
  sim.fun <- function(dat, pars, iter){
    # Run for loop
    for(ii in 1:length(unique(dat$TS))){
      dat %>%
        filter(TS == unique(dat$TS)[ii]) %>%
        group_by(TS) %>%
        mutate(N = n()) %>%
        ungroup() %>%
        filter(Lagged.Year > 1987) -> TS.dat
      
      sim.unsmoothed <- arima.sim(model = list(order = c(1,0,0), 
                                               ar = (pars %>% filter(smooth == "unsmoothed"))$ar1, 
                                               sd= (pars %>% filter(smooth == "unsmoothed"))$sd),
                                  n = nrow(TS.dat)) # joining by lagged year doesn't matter bc it is a random TS
      sim.twoyear <- arima.sim(model = list(order = c(1,0,0), 
                                            ar = (pars %>% filter(smooth == "twoyear"))$ar1, 
                                            sd= (pars %>% filter(smooth == "threeyear"))$sd),
                               n = nrow(TS.dat))
      sim.threeyear <- arima.sim(model = list(order = c(1,0,0), 
                                              ar = (pars %>% filter(smooth == "unsmoothed"))$ar1, 
                                              sd= (pars %>% filter(smooth == "unsmoothed"))$sd),
                                 n = nrow(TS.dat))
      
      TS.dat <- cbind(TS.dat, data.frame(unsmoothed = sim.unsmoothed,
                                         twoyear = sim.twoyear,
                                         threeyear = sim.threeyear))
      
      
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
                                                 rsq = c(rsq.gam.1, rsq.gam.2, rsq.gam.3),
                                                 iteration = iter))
        
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
      filter(BEST == "Y")  -> sim.out
    
    return(sim.out)
    
  }
  
  # Run function for EBS 1000 times
  model.out <- data.frame()
  sim.out <- data.frame
  
  # 1:1000 %>%
  # purrr::map_df(~sim.fun(bsai.r0, ebs.pars, .x)) -> ebs.out
  
  ebs.out %>%
    #group_by(iteration) %>%
    reframe(N_sig = sum(sig == TRUE),
            N = n(),
            prop_sig = N_sig/N) -> ebs.sum.out # calculate proportion significant for each interation
  
  
  mean(ebs.sum.out$prop_sig) -> ebs.null.sig # calculate mean proportion across all iterations
  
  # Run function for GOA 1000 times
  model.out <- data.frame()
  sim.out <- data.frame
  
  # 1:1000 %>%
  #   purrr::map_df(~sim.fun(goa.r0, goa.pars, .x)) -> goa.out
  
  goa.out %>%
    group_by(iteration) %>%
    reframe(N_sig = sum(sig == TRUE),
            N = n(),
            prop_sig = N_sig/N) -> goa.sum.out # calculate proportion significant for each interation
  
  
  mean(goa.sum.out$prop_sig) -> goa.null.sig # calculate mean proportion across all iterations
  
  
  # Make data frame of all outputs
  data.frame("Region" = c("Eastern Bering Sea", "Gulf of Alaska"), 
             "True significant proportion" = c(ebs.true.sig, goa.true.sig), # calculated in question 2 section above
             "Null significant proportion" = c(ebs.null.sig, goa.null.sig)) -> sum.df
  
  # Save output summary table
  write.csv(sum.df, "./Output/true.vs.null.significance.csv")
  
  
### Test question 7: What are the impacts of increasing ar1? -----------------------------------------------------------------
  
  # EBS 
  right_join(ar1var.EBS.r0, 
             ar1var.EBS.sst %>% 
               dplyr::select(ar1, sd, year) %>% 
               rename(sst_ar1 = ar1, sst_sd = sd)) %>%
    na.omit() -> r0.sst_AR1_EBS
  
  
  ar1.fits <- function(dat){
    for(ii in 1:length(unique(dat$TS))){
      
      knts <- c(3, 4, 5)
      for(kk in 1:length(knts)){
        dat %>%
          filter(TS == unique(.$TS)[ii]) -> model.dat
        
        
        # fit models
        gam.ar1 <- gam(ar1 ~ s(sst_ar1, k = knts[kk]), 
                       data= model.dat, correlation = corAR1())
        
        p.gam.ar1 <- summary(gam.ar1$gam)$s.table[,4]
        
        r.gam.ar1 <- summary(gam.ar1$gam)$r.sq
        
        AIC.gam.ar1 <- AIC(gam.ar1)
        
        # Build summary table
        model.out <- rbind(model.out, data.frame(TS = unique(dat$TS)[ii],
                                                 knots = knts[kk],
                                                 p_gam = c(p.gam.ar1),
                                                 rsq_gam = c(r.gam.ar1),
                                                 AIC = c(AIC.gam.ar1)))
        
      } # close knot loop
    } # close timeseries loop
    
    # Label best models by timeseries
    model.out %>%
      group_by(TS) %>%
      mutate(BEST = ifelse(AIC == min(AIC), "Y", "N"),
             padj = p.adjust(p_gam, method = "fdr")) %>%
      filter(BEST == "Y") -> model.out
  }
  
  
  ar1.fits(r0.sst_AR1_EBS) -> EBS.out
  
  model.predict <- data.frame()
  
  best.dat <- EBS.out
  ar1.dat <- r0.sst_AR1_EBS
  
  for(ii in 1:length(unique(EBS.out$TS))){
    best.dat %>%
      filter(TS == unique(best.dat$TS)[ii]) -> model.dat
    
    ar1.dat %>%
      filter(TS == unique(best.dat$TS)[ii]) %>%
      dplyr::select(TS, ar1, sst_ar1) -> TS.dat
    
    gamm.best <- gamm(ar1 ~ s(sst_ar1, k = model.dat$knots), 
                      data= TS.dat, correlation = corAR1())
    
    model.predict <- rbind(model.predict, data.frame(TS = unique(best.dat$TS)[ii],
                                                     pred.ar1 = predict(gamm.best, se.fit =TRUE)$fit,
                                                     pred.CI = 1.96*(predict(gamm.best, se.fit =TRUE)$se.fit),
                                                     sst_ar1 = TS.dat$sst_ar1,
                                                     TS_ar1 = TS.dat$ar1,
                                                     k = model.dat$knots,
                                                     p_gam = model.dat$padj,
                                                     rsq = model.dat$rsq_gam))
    
  }
  
  # Specify plotting labels
  model.predict %>%
    dplyr::select(TS, k, p_gam, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                            (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                            (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                            TRUE ~ paste0("p=", round(p_gam, 2)))) -> lab.dat
  
  labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, ", R2=", round(lab.dat$rsq, 2), ", ",
                 lab.dat$plab, ")")
  names(labs) <- names(r0.labs.bsai)
  
  # Plot BSAI r0/SST
  ggplot()+
    geom_ribbon(model.predict, mapping = aes(x = sst_ar1, ymin = pred.ar1 - pred.CI, ymax= pred.ar1+pred.CI),
                fill = "grey", alpha = 0.75)+
    geom_point(model.predict, mapping=aes(x = sst_ar1, y = TS_ar1))+
    geom_line(model.predict, mapping = aes(x = sst_ar1, y = pred.ar1), size = 1.25, color = "#6A6DB7")+
    facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
               ncol = 4)+
    theme_bw()+
    ggtitle("EBS sst AR1 vs. r0 AR1")+
    ylab("Recruitment AR1")+
    xlab("SST AR1")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> ebs.r0vs.sstAR1.plot
  
  ggsave(plot = ebs.r0vs.sstAR1.plot, "./Figures/ebs.r0vs.sstAR1.plot.png", width = 11, height = 8.5, units = "in")
  
  # GOA 
  right_join(ar1var.goa.r0, 
             ar1var.goa.sst %>% 
               dplyr::select(ar1, sd, window) %>% 
               rename(sst_ar1 = ar1, st_sd = sd)) %>%
    na.omit() -> r0.sst_AR1_goa
  
  
  ar1.fits(r0.sst_AR1_goa) -> goa.out
  
  model.predict <- data.frame()
  
  best.dat <- goa.out
  ar1.dat <- r0.sst_AR1_goa
  
  for(ii in 1:length(unique(best.dat$TS))){
    best.dat %>%
      filter(TS == unique(best.dat$TS)[ii]) -> model.dat
    
    ar1.dat %>%
      filter(TS == unique(best.dat$TS)[ii]) %>%
      dplyr::select(TS, ar1, sst_ar1) -> TS.dat
    
    gamm.best <- gamm(ar1 ~ s(sst_ar1, k = model.dat$knots), 
                      data= TS.dat, correlation = corAR1())
    
    model.predict <- rbind(model.predict, data.frame(TS = unique(best.dat$TS)[ii],
                                                     pred.ar1 = predict(gamm.best, se.fit =TRUE)$fit,
                                                     pred.CI = 1.96*(predict(gamm.best, se.fit =TRUE)$se.fit),
                                                     sst_ar1 = TS.dat$sst_ar1,
                                                     TS_ar1 = TS.dat$ar1,
                                                     k = model.dat$knots,
                                                     p_gam = model.dat$padj,
                                                     rsq = model.dat$rsq_gam))
    
  }
  
  # Specify plotting labels
  model.predict %>%
    dplyr::select(TS, k, p_gam, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    mutate(plab = case_when((p_gam < 0.001) ~ "p<0.001",
                            (p_gam < 0.01 & p_gam >= 0.001) ~ "p<0.01",
                            (p_gam <0.05 & p_gam >= 0.01) ~ "p<0.05",
                            TRUE ~ paste0("p=", round(p_gam, 2)))) -> lab.dat
  
  labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, ", R2=", round(lab.dat$rsq, 2), ", ",
                 lab.dat$plab, ")")
  names(labs) <- names(r0.labs.goa)
  
  # Plot BSAI r0/SST
  ggplot()+
    geom_ribbon(model.predict, mapping = aes(x = sst_ar1, ymin = pred.ar1 - pred.CI, ymax= pred.ar1+pred.CI),
                fill = "grey", alpha = 0.75)+
    geom_point(model.predict, mapping=aes(x = sst_ar1, y = TS_ar1))+
    geom_line(model.predict, mapping = aes(x = sst_ar1, y = pred.ar1), size = 1.25, color = "#A34242")+
    facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
               ncol = 4)+
    theme_bw()+
    ggtitle("GOA sst AR1 vs. r0 AR1")+
    ylab("Recruitment AR1")+
    xlab("SST AR1")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> goa.r0vs.sstAR1.plot
  
  ggsave(plot = goa.r0vs.sstAR1.plot, "./Figures/goa.r0vs.sstAR1.plot.png", width = 11, height = 8.5, units = "in")
  
  
  
  
  
  # Get ar1 for unsmoothed sst
  ar1.vals <- seq(0.1, 0.9, by = 0.1)
  
  # Simulate model 1000 times with increasing ar1
  iter <- 10

 
  
 sim.fun2 <- function(ar1.vals, iter){
   sim.out <- data.frame()  
   
   for(aa in 1:length(ar1.vals)){
     ts <- arima.sim(model = list(order = c(1,0,0),  ar = ar1.vals[aa]), n = 77) # length of ebs SST ts
     
     sim.out <- rbind(sim.out, data.frame(ar1 = rep(ar1.vals[aa], 77),
                                          iteration = rep(iter, 77),
                                          ts.sim = c(ts)))
   }
   return(sim.out)
 }
 
 
 1:1000 %>%
 purrr::map_df(~sim.fun2(ar1.vals, .x)) -> out
 
 out %>%
   group_by(ar1) %>%
   reframe(mean.ts = mean(ts.sim),
           SD = sd(ts.sim)) -> plot.dat
 
 ggplot(plot.dat, aes(as.factor(ar1), SD))+
   geom_line(group=1)+
   geom_point()+
   ylab("Standard deviation")+
   xlab("AR1")+
   ggtitle("Impact of increasing AR1")+
   theme_bw()+
   theme(legend.position = "none",
         axis.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         strip.text = element_text(size = 16),
         legend.text = element_text(size = 14),
         title = element_text(size = 16)) -> sd.ar1.plot
 