#Model exploration

#PURPOSE: To fit models GOA and EBS TS and climate timeseries

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")
source("./Scripts/functions.R")


### MODEL EXPLORATION -------------------------------------------------------------------------------

  # Evaluate candidate models based on unsmoothed, 2-year, and 3-year rolling mean SST ----
    sst.2 <- rollapply(sst$mean.sst, 2, mean, na.rm = T, fill = NA) #2-year
    sst.3 <- rollapply(sst$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
    
    sst.rollmeans <- data.frame(sst, twoyear = sst.2, threeyear = sst.3) %>%
                          rename(unsmoothed = mean.sst)
  
    # RECRUITMENT ----
      # Add rollmean sst to BSAI r0
      bsai.r0 %>%
        right_join(., expand.grid(Year = min(bsai.r0$Year):max(bsai.r0$Year))) %>%
        right_join(., sst.rollmeans %>% filter(region == "Eastern Bering Sea") %>%
                      rename(Lagged.Year = Year), by = c("Lagged.Year")) %>%
        na.omit() %>%
        dplyr::select(!region) -> bsai.r0.sst
      
      # Add rollmean sst to GOA r0
      goa.r0 %>%
        right_join(., expand.grid(Year = min(goa.r0$Year):max(goa.r0$Year))) %>%
        right_join(., sst.rollmeans %>% filter(region == "Gulf of Alaska") %>%
                     rename(Lagged.Year = Year), by = c("Lagged.Year")) %>%
        na.omit() %>%
        dplyr::select(!region) -> goa.r0.sst
      
      
      # Fit models for BSAI
        model.out <- data.frame()
        
        fit.models(bsai.r0.sst) -> bsai.r0.out
        
      
      # Fit models for GOA
        model.out <- data.frame()
        
        fit.models(dat = goa.r0.sst) -> goa.r0.out
        
      # Isolate best models
       bsai.r0.best <- bsai.r0.out %>% filter(BEST == "Y")
       
       goa.r0.best <- goa.r0.out %>% filter(BEST == "Y")
       
       
      # Predict and plot with best models BSAI
      model.predict <- data.frame()
      
      for(ii in 1:length(unique(bsai.r0.best$TS))){
        bsai.r0.best %>%
          filter(TS == unique(bsai.r0.best$TS)[ii]) -> model.dat
        
        bsai.r0.sst %>%
          filter(TS == unique(bsai.r0.best$TS)[ii]) %>%
                   dplyr::select(TS, log.recruitment, grep(model.dat$sst, names(.))) -> TS.dat
        
        gamm.best <- gamm(log.recruitment ~ s(TS.dat[,3], k = model.dat$knots), 
                      data= TS.dat, correlation = corAR1())
        
        model.predict <- rbind(model.predict, data.frame(TS = unique(bsai.r0.best$TS)[ii],
                                                         log.recruitment = TS.dat$log.recruitment,
                                                         pred.r0 = predict(gamm.best, se.fit =TRUE)$fit,
                                                         pred.CI = 1.96*(predict(gamm.best, se.fit =TRUE)$se.fit),
                                                         sst = TS.dat[,3],
                                                         sst.smooth = model.dat$sst,
                                                         k = model.dat$knots))
        
        saveRDS(gamm.best, paste0("./Output/", unique(bsai.r0.best$TS)[ii], ".gamm.rda"))
        
      }
      
      model.predict %>%
        arrange(., TS) -> lab.dat
      
      model.predict %>%
        dplyr::select(TS, sst.smooth, k)%>%
        distinct() %>%
        arrange(., TS) %>%
        mutate(sst.smooth = case_when((sst.smooth == "unsmoothed") ~ "1-year",
                                      (sst.smooth == "twoyear") ~ "2-year",
                                      (sst.smooth == "threeyear") ~ "3-year")) -> lab.dat
      
      labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, ", sst=", lab.dat$sst.smooth, ")")
      names(labs) <- names(r0.labs.bsai)
      
      # Plot BSAI r0/SST
       ggplot()+
         geom_ribbon(model.predict, mapping = aes(x = sst, ymin = pred.r0 - pred.CI, ymax= pred.r0+pred.CI),
                     fill = "grey", alpha = 0.75)+
         geom_point(model.predict, mapping=aes(x = sst, y = log.recruitment))+
         geom_line(model.predict, mapping = aes(x = sst, y = pred.r0), size = 1.25, color = "#6A6DB7")+
         facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
                    ncol = 3)+
        theme_bw()+
         ggtitle("BSAI groundfish/crab recruitment and SST")+
         ylab("log(millions of recruits)")+
         xlab("°C")+
         theme(axis.text = element_text(size = 10),
               axis.title = element_text(size = 12),
               strip.text = element_text(size = 10),
               legend.text = element_text(size = 12)) -> bsai.r0.sst.plot
       
       ggsave(plot = bsai.r0.sst.plot, "./Figures/bsai.r0.sst.plot2.png", width = 8.5, height = 11, units = "in")
       
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
       
    # SALMON CATCH ----
       # Add rollmean sst to salmon catch
       rbind(salmon.catch %>%
               right_join(., expand.grid(Year = min(salmon.catch$Year):max(salmon.catch$Year))) %>%
               filter(Region == "GOA") %>%
               right_join(., sst.rollmeans %>% filter(region == "Gulf of Alaska") %>%
                            rename(Lagged.Year = Year), by = c("Lagged.Year")) %>%
               na.omit() %>%
               dplyr::select(!region),
             salmon.catch %>%
               right_join(., expand.grid(Year = min(salmon.catch$Year):max(salmon.catch$Year))) %>%
               filter(Region == "BSAI") %>%
               right_join(., sst.rollmeans %>% filter(region == "Eastern Bering Sea") %>%
                            rename(Lagged.Year = Year), by = c("Lagged.Year")) %>%
               na.omit() %>%
               dplyr::select(!region))-> salmon.catch.sst
       
       # For loops forever!!
        model.out <- data.frame()
       
        for(ii in 1:length(unique(salmon.catch.sst$TS))){
             # filter data by TS
             salmon.catch.sst %>%
               mutate(era = case_when((Year <=1989) ~ "early",
                                      (Year >1989) ~ "late")) %>%
               filter(TS == unique(salmon.catch.sst$TS)[ii]) %>%
               group_by(TS, era) %>%
               mutate(N = n(), era = ifelse(N <8, "late", era)) %>%
               ungroup() %>%
               na.omit()-> TS.dat
             
             knts <- c(3, 4)
             
             for(kk in 1:length(knts)){
               for(jj in 1:length(unique(TS.dat$era))){
               
               # fit gams, record pval and AIC
               gam.1 <- gamm(log.catch ~ s(unsmoothed, k = knts[kk]), 
                             data= TS.dat %>% filter(era == unique(TS.dat$era)[jj]), correlation = corAR1())
               gam.2 <- gamm(log.catch ~ s(twoyear, k = knts[kk]), 
                             data= TS.dat %>% filter(era == unique(TS.dat$era)[jj]), correlation = corAR1())
               gam.3 <- gamm(log.catch ~ s(threeyear, k = knts[kk]), 
                             data= TS.dat %>% filter(era == unique(TS.dat$era)[jj]), correlation = corAR1())
               
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
               model.out <- rbind(model.out, data.frame(TS = unique(TS.dat$TS),
                                                        knots = knts[kk],
                                                        #Model = c("gam.1", "gam.2", "gam.3"),
                                                        sst = c("unsmoothed", "twoyear", "threeyear"),
                                                        era = unique(TS.dat$era)[jj],
                                                        p_lme = c(p.lme.1, p.lme.2, p.lme.3),
                                                        p_gam = c(p.gam.1, p.gam.2, p.gam.3),
                                                        AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3)))
            
              } # close era loop   
            } # close knot loop
           } # close timeseries loop
           
        # Label best models by timeseries
           model.out %>%
             group_by(TS, era) %>%
             mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) -> model.out.best
           
           
        # Isolate best models
        salm.catch.best <- model.out.best %>% filter(BEST == "Y")
           
        # Fit and save best models BSAI
        model.predict <- data.frame()
           
        for(ii in 1:length(unique(salm.catch.best$TS))){
          salm.catch.best %>%
            filter(TS == unique(salm.catch.best$TS)[10]) -> model.dat
          
          for(jj in 1:nrow(model.dat)){
            
            
            salm.catch.sst %>%
              mutate(era = case_when((Year <=1989) ~ "early",
                                     (Year >1989) ~ "late")) %>%
              filter(TS == unique(model.dat$TS), era == model.dat$era[jj]) %>%
              dplyr::select(TS, log.catch, grep(model.dat$sst[jj], names(.)), era) -> TS.dat
            
            
            gamm.best <- gamm(log.catch ~ s(TS.dat[,3], k = model.dat$knots[jj]), 
                              data= TS.dat, correlation = corAR1())
            
            model.predict <- rbind(model.predict, data.frame(TS = unique(TS.dat$TS),
                                                             log.catch = TS.dat$log.catch,
                                                             pred.catch = predict(gamm.best, se.fit =TRUE)$fit,
                                                             pred.CI = 1.96*(predict(gamm.best, se.fit =TRUE)$se.fit),
                                                             sst = TS.dat[,3],
                                                             sst.smooth = model.dat$sst[jj],
                                                             k = model.dat$knots[jj],
                                                             era = model.dat$era[jj]))
            
            saveRDS(gamm.best, paste0("./Output/", unique(salm.catch.best$TS)[ii], ".", model.dat$era[jj], ".gamm.rda"))
            
          }
        }
     
        
        model.predict %>%
          arrange(., TS) -> lab.dat
        
        model.predict %>%
          dplyr::select(TS, sst.smooth, k)%>%
          distinct() %>%
          arrange(., TS) %>%
          mutate(sst.smooth = case_when((sst.smooth == "unsmoothed") ~ "1-year",
                                        (sst.smooth == "twoyear") ~ "2-year",
                                        (sst.smooth == "threeyear") ~ "3-year")) -> lab.dat
        
        labs <- paste0(salm.labs, " \n(k=", lab.dat$k, ", sst=", lab.dat$sst.smooth, ")")
        names(labs) <- names(salm.labs)
        
      # Plot
       ggplot()+
          # geom_ribbon(model.predict, mapping = aes(x = sst, ymin = pred.catch - pred.CI, ymax= pred.catch+pred.CI),
          #             fill = "grey", alpha = 0.75)+
          geom_point(model.predict, mapping=aes(x = sst, y = log.catch))+
          #geom_line(model.predict, mapping = aes(x = sst, y = pred.catch), size = 1.25, color = "#6A6DB7")+
          facet_wrap(~TS, scales = "free", labeller = labeller(TS = labs),
                     ncol = 3)+
          theme_bw()+
          ggtitle("BSAI groundfish/crab recruitment and SST")+
          ylab("log(millions of recruits)")+
          xlab("°C")+
          theme(axis.text = element_text(size = 10),
                axis.title = element_text(size = 12),
                strip.text = element_text(size = 10),
                legend.text = element_text(size = 12)) -> bsai.r0.sst.plot
        
  # Calculate AR1 and CV/SD for 15 year windows for biology/SST timeseries ----
  
     
    # BSAI SSB ----
    sum.out <- data.frame()
    
    sum.fun(bsai.ssb, "SSB") -> sum.ssb.bsai 
    
    sum.ssb.bsai %>%
      pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.bsai), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(bsai.ssb$Year), max(bsai.ssb$Year), by = 15),
                         limits = c(min(bsai.ssb$Year), max(bsai.ssb$Year)))+     
      ylab("AR1")+
      ggtitle("BSAI SSB AR1")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/bsai.ssb.AR1.png", width =11, height = 8.5,
           units = "in")
    
    
    # Plot CV
    ggplot(long.dat %>% filter(Type %in% c("cv")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.bsai), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(bsai.ssb$Year), max(bsai.ssb$Year), by = 15),
                         limits = c(min(bsai.ssb$Year), max(bsai.ssb$Year)))+      ylab("CV")+
      ggtitle("BSAI SSB CV")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/bsai.ssb.CV.png", width =11, height = 8.5,
           units = "in")
    
    # GOA SSB ----
    sum.out <- data.frame()
    
    sum.fun(goa.ssb, "SSB") -> sum.ssb.goa
    
    sum.ssb.goa %>%
      pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#A34242")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.goa), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(goa.ssb$Year), max(goa.ssb$Year), by = 15),
                         limits = c(min(goa.ssb$Year), max(goa.ssb$Year)))+      ylab("AR1")+
      ggtitle("GOA SSB AR1")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/goa.ssb.AR1.png", width =11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("cv")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#A34242")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.goa), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(goa.ssb$Year), max(goa.ssb$Year), by = 15),
                         limits = c(min(goa.ssb$Year), max(goa.ssb$Year)))+
      ylab("CV")+
      ggtitle("GOA SSB CV")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/goa.ssb.CV.png", width =11, height = 8.5, units = "in")
    
    # BSAI R0 ----
    sum.out <- data.frame()
    
    sum.fun(bsai.r0, "Recruitment") -> sum.r0.bsai
    
    sum.r0.bsai %>%
      pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.bsai), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(bsai.r0$Year), max(bsai.r0$Year), by = 15),
                         limits = c(min(bsai.r0$Year), max(bsai.r0$Year)))+     
      ylab("AR1")+
      ggtitle("BSAI Recruitment AR1")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/bsai.r0.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("cv")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.bsai), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(bsai.r0$Year), max(bsai.r0$Year), by = 15),
                         limits = c(min(bsai.r0$Year), max(bsai.r0$Year)))+     
      ylab("CV")+
      ggtitle("BSAI Recruitment CV")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/bsai.r0.CV.png", width = 11, height = 8.5, units = "in")
    
    
    # GOA R0 ----
    sum.out <- data.frame()
    
    sum.fun(goa.r0, "Recruitment") -> sum.r0.goa
    
    sum.r0.goa %>%
      pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#A34242")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.goa), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(goa.r0$Year), max(goa.r0$Year), by = 15),
                         limits = c(min(goa.r0$Year), max(goa.r0$Year)))+     
      ylab("AR1")+
      ggtitle("GOA Recruitment AR1")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/goa.r0.AR1.png", width = 11, height = 8.5, units = "in")
    
   
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("cv")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#A34242")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.goa), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(goa.r0$Year), max(goa.r0$Year), by = 15),
                         limits = c(min(goa.r0$Year), max(goa.r0$Year)))+     
      ylab("CV")+
      ggtitle("GOA Recruitment CV")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/goa.r0.CV.png", width = 11, height = 8.5, units = "in")
    
    # CRAB MB ----
    sum.out <- data.frame()
    
    sum.fun(crab.mb, "Mature biomass") -> sum.mb.crab
    
    sum.mb.crab %>%
      pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") %>%
      mutate(type = ifelse(grepl("mmb", Type) == TRUE, 
                           "Mature male biomass", "Mature female biomass")) -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(grepl("ar1", Type) == TRUE),
           mapping = aes(window, Value, linetype = type))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", nrow = 3)+
      theme_bw()+
      scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("Female", "Male"))+
      ylab("AR1")+
      ggtitle("Crab mature biomass AR1")+
      scale_x_continuous(breaks = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5),
                         labels = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5),
                         limits = c(min(crab.mb$Year), max(crab.mb$Year)))+
      theme(legend.position = "bottom",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12))
    
    ggsave("./Figures/crab.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(grepl("cv", Type) == TRUE),
           mapping = aes(window, Value, linetype = type))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", nrow = 3)+
      theme_bw()+
      scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("Female", "Male"))+
      ylab("CV")+
      ggtitle("Crab mature biomass CV")+
      scale_x_continuous(breaks = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5),
                         labels = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5),
                         limits = c(min(crab.mb$Year), max(crab.mb$Year)))+
      theme(legend.position = "bottom",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12))
    
    ggsave("./Figures/crab.CV.png", width = 11, height = 8.5, units = "in")
    
    
    # BSAI and GOA salmon ----
    sum.out <- data.frame()
    
    sum.fun(bsai.salmon, "Catch") -> sum.catch.bsai
    
    sum.out <- data.frame()
    
    sum.fun(goa.salmon, "Catch") -> sum.catch.goa
    
    rbind(sum.catch.goa %>%
          pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") %>%
          mutate(Region = "GOA"),
        sum.catch.bsai %>%
          pivot_longer(!c(TS,window), values_to = "Value", names_to = "Type") %>%
          mutate(Region = "BSAI")) -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1")), mapping = aes(window, Value, color = Region))+
      #ggplot(salmon.catch %>% filter(Year > 1959), aes(x = Year, y = log.catch))+
      geom_line(linewidth = 1)+
      scale_color_manual(values = c("#6A6DB7", "#A34242"))+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = salm.labs),
                 ncol = 3)+
      scale_x_continuous(breaks = seq(min(salmon.catch$Year), max(salmon.catch$Year), by = 30),
                         limits = c(min(salmon.catch$Year), max(salmon.catch$Year)))+
      theme_bw()+
      ggtitle("BSAI and GOA salmon catch AR1")+
      ylab("AR1")+
      theme(axis.text = element_text(size = 12),
            legend.position  = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12))
    
    ggsave("./Figures/salm.catch.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot CV
    ggplot(long.dat %>% filter(Type %in% c("cv")), mapping = aes(window, Value, color = Region))+
      #ggplot(salmon.catch %>% filter(Year > 1959), aes(x = Year, y = log.catch))+
      geom_line(linewidth = 1)+
      scale_color_manual(values = c("#6A6DB7", "#A34242"))+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = salm.labs),
                 ncol = 3)+
      scale_x_continuous(breaks = seq(min(salmon.catch$Year), max(salmon.catch$Year), by = 30),
                         limits = c(min(salmon.catch$Year), max(salmon.catch$Year)))+
      theme_bw()+
      ggtitle("BSAI and GOA salmon catch CV")+
      ylab("CV")+
      theme(axis.text = element_text(size = 12),
            legend.position  = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12))
    
    ggsave("./Figures/salm.catch.CV.png", width = 11, height = 8.5, units = "in")
    
    # SST ----
      # EBS SST
      sum.out <- data.frame()
      
      sum.fun(ebs.sst, "SST") %>%
        mutate(region = "BSAI") -> sum.sst.bsai
    
      # GOA SST
      sum.out <- data.frame()
      
      sum.fun(goa.sst, "SST") %>%
        mutate(region = "GOA") -> sum.sst.goa
      
      rbind(sum.sst.bsai, sum.sst.goa) %>%
          pivot_longer(!c(TS, region, window), values_to = "Value", names_to = "Type") -> long.dat
      
      # Plot AR1
      ggplot(long.dat %>% filter(Type %in% c("ar1")), aes(x = window, y = Value, color = region))+
        facet_wrap(~region, scales = "free_y", nrow = 2)+
        scale_color_manual(values = c("#6A6DB7", "#A34242"))+
        geom_line(size = 1.5)+
        #geom_point(size = 1.5)+
        scale_x_continuous(breaks = seq(min(sst$Year), max(sst$Year), by = 10),
                           limits = c(min(sst$Year), max(sst$Year)))+
        theme_bw()+
        ggtitle("ERA5 SST AR1")+
        ylab("AR1")+
        theme(legend.position = "none",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16),
              legend.text = element_text(size = 14),
              title = element_text(size = 16))
      
      ggsave("./Figures/SST.AR1.png", width = 11, height = 8.5, units = "in")
      
      # Plot SD
      ggplot(long.dat %>% filter(Type %in% c("sd")), aes(x = window, y = Value, color = region))+
        facet_wrap(~region, scales = "free_y", nrow = 2)+
        scale_color_manual(values = c("#6A6DB7", "#A34242"))+
        geom_line(size = 1.5)+
        #geom_point(size = 1.5)+
        scale_x_continuous(breaks = seq(min(sst$Year), max(sst$Year), by = 10),
                           limits = c(min(sst$Year), max(sst$Year)))+
        theme_bw()+
        ggtitle("ERA5 SST SD")+
        ylab("SD")+
        theme(legend.position = "none",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16),
              legend.text = element_text(size = 14),
              title = element_text(size = 16))
      
      ggsave("./Figures/SST.SD.png", width = 11, height = 8.5, units = "in")
      
    
 
 
  