#Model exploration

#PURPOSE: To fit models GOA and EBS TS and climate timeseries

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")

### MODEL EXPLORATION -------------------------------------------------------------------------------

  # Evaluate candidate models based on unsmoothed, 2-year, and 3-year rolling mean SST ----
    sst.2 <- rollapply(sst$mean.sst, 2, mean, na.rm = T, fill = NA) #2-year
    sst.3 <- rollapply(sst$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
    
    sst.rollmeans <- data.frame(sst, sst.2 = sst.2, sst.3 = sst.3) %>%
                          rename(sst = mean.sst)
  
     # RECRUITMENT
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
            for(jj in 1:length(unique(TS.dat$era))){
              
              # fit gams, record pval and AIC
              gam.1 <- gamm(log.recruitment ~ s(sst, k = knts[kk]), 
                           data= TS.dat %>% filter(era == unique(TS.dat$era)[jj]), correlation = corAR1())
              gam.2 <- gamm(log.recruitment ~ s(sst.2, k = knts[kk]), 
                            data= TS.dat %>% filter(era == unique(TS.dat$era)[jj]), correlation = corAR1())
              gam.3 <- gamm(log.recruitment ~ s(sst.3, k = knts[kk]), 
                            data= TS.dat %>% filter(era == unique(TS.dat$era)[jj]), correlation = corAR1())
              
              p.gam.1 <- signif(summary(gam.1$gam)$s.table[,4],2)
              p.gam.2 <- signif(summary(gam.2$gam)$s.table[,4],2)
              p.gam.3 <- signif(summary(gam.3$gam)$s.table[,4],2)
              
              AIC.gam.1 <- AICc(gam.1)
              AIC.gam.2 <- AICc(gam.2)
              AIC.gam.3 <- AICc(gam.3)
  
              
              # Build summary table
              model.out <- rbind(model.out, data.frame(TS = unique(dat$TS)[ii],
                                                       knots = knts[kk],
                                                       era = unique(TS.dat$era)[jj],
                                                       #Model = c("gam.1", "gam.2", "gam.3"),
                                                       smooth = c("unsmoothed", "2-year", "3-year"),
                                                       p_val = c(p.gam.1, p.gam.2, p.gam.3),
                                                       AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3),
                                                       N = TS.dat %>% 
                                                            filter(era ==unique(TS.dat$era)[jj]) %>%
                                                         dplyr::select(N) %>%
                                                         pull() %>%
                                                         unique()))
            } # close era loop
          } # close knot loop
        } # close timeseries loop
        
        # Label best models by timeseries
        model.out %>%
          group_by(TS, era) %>%
          mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) -> model.out
        
        return(model.out)
      }
      
    
      # Fit models for BSAI
        model.out <- data.frame()
        
        fit.models(bsai.r0.sst) -> bsai.r0.out
        
      
      # Fit models for GOA
        model.out <- data.frame()
        
        fit.models(dat = goa.r0.sst) -> goa.r0.out
        
        
          
    
  # Calculate AR1 and SD for 15 year windows for biology/SST timeseries ----
  
    # Create function to calculate values
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
        
        # Calculate windows
        window <- seq(min(TS.dat$Year)+(width-1), max(TS.dat$Year), by = 1)
        
        
        # Compile output
        sum.out <- rbind(sum.out, data.frame(TS = unique(dat$TS)[ii],
                                             ar1 = ar1,
                                             sd = sd,
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
          
          
          # Calculate windows
          window <- seq(min(TS.dat.mmb$Year)+(width-1), max(TS.dat.mmb$Year), by = 1)
          
   
          
        } else{
          # Calculate rolling window AR1
          ar1.mmb <- sapply(rollapply(TS.dat.mmb$Value, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 

          # Calculate rolling window SD
          sd.mmb <-  rollapply(TS.dat.mmb$Value, width = width, FUN = sd)

          
          # Calculate windows
          window <- seq(min(TS.dat.mmb$Year)+(width-1), max(TS.dat.mmb$Year), by = 1)
          
          ar1.fmb <- sd.fmb <- NA
          
        }
        
        
        # Compile output
        sum.out <- rbind(sum.out, data.frame(TS = unique(dat$TS)[ii],
                                             ar1.mmb = ar1.mmb, 
                                             ar1.fmb = ar1.fmb, 
                                             sd.mmb = sd.mmb, 
                                             sd.fmb = sd.fmb,
                                             window = window))
        
      }
    }
    
      
      return(sum.out)
    }
    
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
    
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.bsai), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(bsai.ssb$Year), max(bsai.ssb$Year), by = 15),
                         limits = c(min(bsai.ssb$Year), max(bsai.ssb$Year)))+      ylab("SD")+
      ggtitle("BSAI SSB SD")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/bsai.ssb.SD.png", width =11, height = 8.5,
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
    ggplot(long.dat %>% filter(Type %in% c("sd")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#A34242")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.goa), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(goa.ssb$Year), max(goa.ssb$Year), by = 15),
                         limits = c(min(goa.ssb$Year), max(goa.ssb$Year)))+
      ylab("SD")+
      ggtitle("GOA SSB SD")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/goa.ssb.SD.png", width =11, height = 8.5, units = "in")
    
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
    ggplot(long.dat %>% filter(Type %in% c("sd")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.bsai), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(bsai.r0$Year), max(bsai.r0$Year), by = 15),
                         limits = c(min(bsai.r0$Year), max(bsai.r0$Year)))+     
      ylab("SD")+
      ggtitle("BSAI Recruitment SD")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/bsai.r0.SD.png", width = 11, height = 8.5, units = "in")
    
    
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
    ggplot(long.dat %>% filter(Type %in% c("sd")), mapping = aes(window, Value))+
      geom_line(linewidth = 1, color = "#A34242")+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.goa), ncol = 3)+
      theme_bw()+
      scale_x_continuous(breaks = seq(min(goa.r0$Year), max(goa.r0$Year), by = 15),
                         limits = c(min(goa.r0$Year), max(goa.r0$Year)))+     
      ylab("SD")+
      ggtitle("GOA Recruitment SD")+
      theme(axis.text = element_text(size = 12),
            legend.position = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) 
    
    
    ggsave("./Figures/goa.r0.SD.png", width = 11, height = 8.5, units = "in")
    
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
    ggplot(long.dat %>% filter(grepl("sd", Type) == TRUE),
           mapping = aes(window, Value, linetype = type))+
      geom_line(linewidth = 1, color = "#6A6DB7")+
      facet_wrap(~TS, scales = "free_y", nrow = 3)+
      theme_bw()+
      scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("Female", "Male"))+
      ylab("SD")+
      ggtitle("Crab mature biomass SD")+
      scale_x_continuous(breaks = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5),
                         labels = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5),
                         limits = c(min(crab.mb$Year), max(crab.mb$Year)))+
      theme(legend.position = "bottom",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12))
    
    ggsave("./Figures/crab.SD.png", width = 11, height = 8.5, units = "in")
    
    
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
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd")), mapping = aes(window, Value, color = Region))+
      #ggplot(salmon.catch %>% filter(Year > 1959), aes(x = Year, y = log.catch))+
      geom_line(linewidth = 1)+
      scale_color_manual(values = c("#6A6DB7", "#A34242"))+
      facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = salm.labs),
                 ncol = 3)+
      scale_x_continuous(breaks = seq(min(salmon.catch$Year), max(salmon.catch$Year), by = 30),
                         limits = c(min(salmon.catch$Year), max(salmon.catch$Year)))+
      theme_bw()+
      ggtitle("BSAI and GOA salmon catch SD")+
      ylab("AR1")+
      theme(axis.text = element_text(size = 12),
            legend.position  = "none",
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12))
    
    ggsave("./Figures/salm.catch.SD.png", width = 11, height = 8.5, units = "in")
    
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
      
    
 
 