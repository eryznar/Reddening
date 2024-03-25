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
            filter(TS == unique(dat$TS)[ii]) -> TS.dat
          
          # fit gams, record pval and AIC
          gam.1 <- gam(log.recruitment ~ s(sst, k = 3), data= TS.dat, method = "REML")
          gam.2 <- gam(log.recruitment ~ s(sst.2, k = 3), data= TS.dat, method = "REML")
          gam.3 <- gam(log.recruitment ~ s(sst.3, k = 3), data= TS.dat, method = "REML")
          
          p.gam.1 <- signif(summary(gam.1)$s.table[,4], 2)
          p.gam.2 <- signif(summary(gam.2)$s.table[,4],2)
          p.gam.3 <- signif(summary(gam.3)$s.table[,4],2)
          
          AIC.gam.1 <- AICc(gam.1)
          AIC.gam.2 <- AICc(gam.2)
          AIC.gam.3 <- AICc(gam.3)
          
          # fit gls, record pval and AIC
          gls.1 <- gls(log.recruitment ~ sst, data = TS.dat, correlation = corAR1())
          gls.2 <- gls(log.recruitment ~ sst.2, data = TS.dat, correlation = corAR1())
          gls.3 <- gls(log.recruitment ~ sst.3, data = TS.dat, correlation = corAR1())
          
          p.gls.1 <- signif(summary(gls.1)$tTable[2,4], 2)
          p.gls.2 <- signif(summary(gls.2)$tTable[2,4], 2)
          p.gls.3 <- signif(summary(gls.3)$tTable[2,4], 2)
          
          AIC.gls.1 <- AICc(gls.1)
          AIC.gls.2 <- AICc(gls.2)
          AIC.gls.3 <- AICc(gls.3)
          
          # Build summary table
          model.out <- rbind(model.out, data.frame(TS = unique(dat$TS)[ii],
                                                   Model = c("gam.1", "gam.2", "gam.3", "gls.1", "gls.2", "gls.3"),
                                                   p_val = c(p.gam.1, p.gam.2, p.gam.3, p.gls.1, p.gls.2, p.gls.3),
                                                   AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3, AIC.gls.1, AIC.gls.2, AIC.gls.3)))
          
        }
        
        return(model.out)
      }
      
      # Fit models for BSAI
        model.out <- data.frame()
        
        fit.models(bsai.r0.sst) -> bsai.r0.out
        
        bsai.r0.out %>%
          group_by(Model) %>%
          reframe(mean.AIC = mean(AIC)) -> bsai.sum
        
        as.data.frame(bsai.r0.out) %>%
          filter(grepl("gam", Model) == TRUE) -> out
        
        ggplot(out, mapping = aes(TS, abs(AIC), fill = Model))+
          geom_bar(stat = "identity", position = "dodge")+
          theme_bw() +
          scale_fill_manual(values = c("#0E3F5C","#579C97","#D1FBD4"), 
                            labels = c("Unsmoothed", "2-year mean", "3-year mean"),
                            name = "")+
          geom_text(out, mapping = aes(x = TS, y = abs(AIC)+2), label = ifelse(out$p_val < 0.05, "*", ""),
                    position = position_dodge(width = 1), size = 8)+
          scale_x_discrete(breaks = names(r0.labs.bsai), labels = function(x)
            str_wrap(r0.labs.bsai, width = 17))+
          ylab("abs(AIC)")+
          ggtitle("BSAI SST vs. recruitment AIC")+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                axis.title.x = element_blank(),
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 16),
                legend.position = "bottom",
                legend.text =  element_text(size = 14),
                title = element_text(size = 16))
        
        ggsave("./Figures/bsai.r0.SST.png", width = 11, height = 8.5, units = "in")
        
      
      # Fit models for GOA
        model.out <- data.frame()
        
        fit.models(dat = goa.r0.sst) -> goa.r0.out
        
          goa.r0.out %>%
            group_by(Model) %>%
            reframe(mean.AIC = mean(AIC)) -> goa.sum
          
          as.data.frame(goa.r0.out) %>%
            filter(grepl("gam", Model) == TRUE) -> out
          
          ggplot(out, mapping = aes(TS, abs(AIC), fill = Model))+
            geom_bar(stat = "identity", position = "dodge")+
            theme_bw() +
            scale_fill_manual(values = c("#0E3F5C","#579C97","#D1FBD4"), 
                              labels = c("Unsmoothed", "2-year mean", "3-year mean"),
                              name = "")+
            geom_text(out, mapping = aes(x = TS, y = abs(AIC)+2), label = ifelse(out$p_val < 0.05, "*", ""),
                      position = position_dodge(width = 1), size = 8)+
            scale_x_discrete(breaks = names(r0.labs.goa), labels = function(x)
              str_wrap(r0.labs.goa, width = 17))+
            ylab("abs(AIC)")+
            ggtitle("GOA SST vs. recruitment AIC")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  axis.title.x = element_blank(),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  legend.position = "bottom",
                  legend.text =  element_text(size = 14),
                  title = element_text(size = 16))
          
          ggsave("./Figures/goa.r0.SST.png", width = 11, height = 8.5, units = "in")
          
    
  # Calculate AR1 and SD for 15 year windows for biology/SST timeseries ----
    # Specify 15-year windows
    win.1 <- seq(min(years), min(years)+14, by= 1)
    win.2 <- seq(max(win.1)+1, max(win.1)+15)
    
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
          filter(TS == unique(dat$TS)[ii], Year %in% years) %>%
          na.omit() -> TS.dat
        
        win1.dat <- TS.dat %>%
                      filter(Year %in% win.1)
        
        win2.dat <- TS.dat %>%
                      filter(Year %in% win.2)
        
        # Calculate AR1
        ar1.win1 <- acf(win1.dat$Value)$acf[2]
        ar1.win2 <- acf(win2.dat$Value)$acf[2]
        
        # Calculate SD
        sd.win1 <- sd(win1.dat$Value)
        sd.win2 <- sd(win2.dat$Value)
        
        # Compile output
        sum.out <- rbind(sum.out, data.frame(TS = unique(dat$TS)[ii],
                                             ar1.win1 = ar1.win1, 
                                             ar1.win2 = ar1.win2,
                                             sd.win1 = sd.win1, 
                                             sd.win2 = sd.win2))
        
      }
    } else{
      for(ii in 1:length(unique(dat$TS))){ # for crab
        
        dat %>%
          filter(TS == unique(dat$TS)[ii], Year %in% years, Type == "mmb") %>%
          na.omit() -> TS.dat.mmb
      
        dat %>%
          filter(TS == unique(dat$TS)[ii], Year %in% years, Type == "fmb") %>%
          na.omit() -> TS.dat.fmb
        
        # Filter by windows
        win1.dat.mmb <- TS.dat.mmb %>%
          filter(Year %in% win.1)
        
        win2.dat.mmb <- TS.dat.mmb %>%
          filter(Year %in% win.2)
        
        win1.dat.fmb <- TS.dat.fmb %>%
          filter(Year %in% win.1)
        
        win2.dat.fmb <- TS.dat.fmb %>%
          filter(Year %in% win.2)
        
        # Calculate AR1
        
        if(unique(dat$TS)[ii] != "Bristol Bay red king crab"){
          ar1.win1.mmb <- acf(win1.dat.mmb$Value)$acf[2]
          ar1.win2.mmb <- acf(win2.dat.mmb$Value)$acf[2]
          
          ar1.win1.fmb <- acf(win1.dat.fmb$Value)$acf[2]
          ar1.win2.fmb <- acf(win2.dat.fmb$Value)$acf[2]
          
          # Calculate SD
          sd.win1.mmb <- sd(win1.dat.mmb$Value)
          sd.win2.mmb <- sd(win2.dat.mmb$Value)
          
          sd.win1.fmb <- sd(win1.dat.fmb$Value)
          sd.win2.fmb <- sd(win2.dat.fmb$Value)
          
        } else{
          ar1.win1.mmb <- acf(win1.dat.mmb$Value)$acf[2]
          ar1.win2.mmb <- acf(win2.dat.mmb$Value)$acf[2]
          
          sd.win1.mmb <- sd(win1.dat.mmb$Value)
          sd.win2.mmb <- sd(win2.dat.mmb$Value)
          
          ar1.win1.fmb <- ar1.win2.fmb <- sd.win1.fmb <- sd.win2.fmb <- NA
          
        }
        
        
        # Compile output
        sum.out <- rbind(sum.out, data.frame(TS = unique(dat$TS)[ii],
                                             ar1.win1.mmb = ar1.win1.mmb, 
                                             ar1.win2.mmb = ar1.win2.mmb,
                                             ar1.win1.fmb = ar1.win1.fmb, 
                                             ar1.win2.fmb = ar1.win2.fmb,
                                             sd.win1.mmb = sd.win1.mmb, 
                                             sd.win2.mmb = sd.win2.mmb,
                                             sd.win1.fmb = sd.win1.fmb, 
                                             sd.win2.fmb = sd.win2.fmb))
        
      }
    }
    
      
      return(sum.out)
    }
    
    # BSAI SSB
    sum.out <- data.frame()
    
    sum.fun(bsai.ssb, "SSB") -> sum.ssb.bsai 
    
    sum.ssb.bsai %>%
      pivot_longer(!TS, values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1.win1", "ar1.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = names(ssb.labs.bsai), labels = function(x)
        str_wrap(ssb.labs.bsai, width = 17))+
      ylab("AR1")+
      ggtitle("BSAI SSB AR1")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/bsai.ssb.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd.win1", "sd.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ggtitle("BSAI SSB SD")+
      scale_x_discrete(breaks = names(ssb.labs.bsai), labels = function(x)
        str_wrap(ssb.labs.bsai, width = 17))+
      ylab("SD")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/bsai.ssb.SD.png", width = 11, height = 8.5, units = "in")
    
    # GOA SSB
    sum.out <- data.frame()
    
    sum.fun(goa.ssb, "SSB") -> sum.ssb.goa
    
    sum.ssb.goa %>%
      pivot_longer(!TS, values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1.win1", "ar1.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = names(ssb.labs.goa), labels = function(x)
        str_wrap(ssb.labs.goa, width = 17))+
      ylab("AR1")+
      ggtitle("GOA SSB AR1")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/goa.ssb.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd.win1", "sd.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ggtitle("GOA SSB SD")+
      scale_x_discrete(breaks = names(ssb.labs.goa), labels = function(x)
        str_wrap(ssb.labs.goa, width = 17))+
      ylab("SD")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/goa.ssb.SD.png", width = 11, height = 8.5, units = "in")
    
    # BSAI R0
    sum.out <- data.frame()
    
    sum.fun(bsai.r0, "Recruitment") -> sum.r0.bsai
    
    sum.r0.bsai %>%
      pivot_longer(!TS, values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1.win1", "ar1.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = names(r0.labs.bsai), labels = function(x)
        str_wrap(r0.labs.bsai, width = 17))+
      ylab("AR1")+
      ggtitle("BSAI recruitment AR1")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/bsai.r0.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd.win1", "sd.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ggtitle("BSAI recruitment SD")+
      scale_x_discrete(breaks = names(r0.labs.bsai), labels = function(x)
                       str_wrap(r0.labs.bsai, width = 17))+
      ylab("SD")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/bsai.r0.SD.png", width = 11, height = 8.5, units = "in")
    
    
    # GOA R0
    sum.out <- data.frame()
    
    sum.fun(goa.r0, "Recruitment") -> sum.r0.goa
    
    sum.r0.goa %>%
      pivot_longer(!TS, values_to = "Value", names_to = "Type") -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1.win1", "ar1.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = names(r0.labs.goa), labels = function(x)
        str_wrap(r0.labs.goa, width = 17))+
      ylab("AR1")+
      ggtitle("GOA recruitment AR1")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/goa.r0.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd.win1", "sd.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ggtitle("GOA recruitment SD")+
      scale_x_discrete(breaks = names(r0.labs.goa), labels = function(x)
        str_wrap(r0.labs.goa, width = 17))+
      ylab("SD")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/goa.r0.SD.png", width = 11, height = 8.5, units = "in")
    
    # CRAB MB
    sum.out <- data.frame()
    
    sum.fun(crab.mb, "Mature biomass") -> sum.mb.crab
    
    sum.mb.crab %>%
      pivot_longer(!TS, values_to = "Value", names_to = "Type") %>%
      mutate(window = ifelse(grepl("win1", Type) == TRUE, "win1", "win2"),
             type = ifelse(grepl("mmb", Type) == TRUE, 
                           "Mature male biomass", "Mature female biomass")) -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(grepl("ar1", Type) == TRUE),
           mapping = aes(TS, Value, fill = window))+
      geom_bar(stat = "identity", position = "dodge")+
      facet_wrap(~type, strip.position = "bottom", scales = "free_x")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = c("Bristol Bay red king crab", "Snow crab",
                                  "Tanner crab"), labels = function(x)
        str_wrap(c("Bristol Bay red king crab", "Snow crab",
                   "Tanner crab"), width = 17))+
      ylab("AR1")+
      ggtitle("Crab mature biomass AR1")+
      theme(
            axis.title.x = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 16),
            strip.background = element_blank(),
            strip.text = element_text(size = 14),
            strip.placement = "outside",
            legend.position = "bottom",
            title = element_text(size = 16),
            legend.text =  element_text(size = 14))
    
    ggsave("./Figures/crab.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(grepl("sd", Type) == TRUE),
           mapping = aes(TS, Value, fill = window))+
      geom_bar(stat = "identity", position = "dodge", na.rm = TRUE)+
      facet_wrap(~type, strip.position = "bottom", scales = "free_x")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = c("Bristol Bay red king crab", "Snow crab",
                                  "Tanner crab"), labels = function(x)
                                    str_wrap(c("Bristol Bay red king crab", "Snow crab",
                                               "Tanner crab"), width = 17))+
      ylab("SD")+
      ggtitle("Crab mature biomass SD")+
      theme(
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        strip.placement = "outside",
        legend.position = "bottom",
        title = element_text(size = 16),
        legend.text =  element_text(size = 14))
    
    ggsave("./Figures/crab.SD.png", width = 11, height = 8.5, units = "in")
    
    
    # BSAI and GOA salmon
    sum.out <- data.frame()
    
    sum.fun(bsai.salmon, "Catch") -> sum.catch.bsai
    
    sum.out <- data.frame()
    
    sum.fun(goa.salmon, "Catch") -> sum.catch.goa
    
    rbind(sum.catch.goa %>%
          pivot_longer(!TS, values_to = "Value", names_to = "Type"),
        sum.catch.bsai %>%
          pivot_longer(!TS, values_to = "Value", names_to = "Type")) -> long.dat
    
    # Plot AR1
    ggplot(long.dat %>% filter(Type %in% c("ar1.win1", "ar1.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      scale_x_discrete(breaks = names(salm.labs), labels = function(x)
        str_wrap(salm.labs, width = 17))+
      ylab("AR1")+
      ggtitle("BSAI and GOA salmon catch AR1")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/salm.catch.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(Type %in% c("sd.win1", "sd.win2")), mapping = aes(TS, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ggtitle("BSAI and GOA salmon catch SD")+
      scale_x_discrete(breaks = names(salm.labs), labels = function(x)
        str_wrap(salm.labs, width = 17))+
      ylab("SD")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.position = "bottom",
            legend.text =  element_text(size = 14),
            title = element_text(size = 16))
    
    ggsave("./Figures/salm.catch.SD.png", width = 11, height = 8.5, units = "in")
    
    

    
    # EBS SST
    sum.out <- data.frame()
    
    sum.fun(ebs.sst, "SST") %>%
      mutate(region = "BSAI") -> sum.sst.bsai
  
    # GOA SST
    sum.out <- data.frame()
    
    sum.fun(goa.sst, "SST") %>%
      mutate(region = "GOA") -> sum.sst.goa
    
    rbind(sum.sst.bsai, sum.sst.goa) %>%
        pivot_longer(!c(TS, region), values_to = "Value", names_to = "Type") -> long.dat
    
    
    # Plot AR1
    ggplot(long.dat %>% filter(grepl("ar1", Type) == TRUE),
           mapping = aes(region, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ylab("AR1")+
      ggtitle("SST AR1")+
      theme(
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        strip.placement = "outside",
        legend.position = "bottom",
        title = element_text(size = 16),
        legend.text =  element_text(size = 16))
    
    ggsave("./Figures/SST.AR1.png", width = 11, height = 8.5, units = "in")
    
    # Plot SD
    ggplot(long.dat %>% filter(grepl("sd", Type) == TRUE),
           mapping = aes(region, Value, fill = Type))+
      geom_bar(stat = "identity", position = "dodge")+
      theme_bw() +
      scale_fill_manual(values = c("cadetblue", "turquoise"), labels = c("1995:2009", "2010:2021"),
                        name = "")+
      ylab("SD")+
      ggtitle("SST SD")+
      theme(
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        strip.placement = "outside",
        legend.position = "bottom",
        title = element_text(size = 16),
        legend.text =  element_text(size = 16))
    
    ggsave("./Figures/SST.SD.png", width = 11, height = 8.5, units = "in")
    
  
 
  # Rolling SD ------------------
    salmon.catch %>%
      filter(Year %in% years) %>%
      na.omit() %>%
      mutate(roll.SD.2 = roll_sd(log.catch, width = 2),
             roll.SD.3 = roll_sd(log.catch, width = 3)) %>%
      pivot_longer(c(log.catch, roll.SD.2, roll.SD.3), names_to = "window", values_to = "value") -> salm.roll 
    
    ggplot(salm.roll, aes(x = Year, roll.SD, color = window))+
      facet_wrap(~TS, scales = "free_y")+
      geom_line(size = 1)+
      geom_point(size = 1)
    
    bsai.r0 %>%
      filter(Year %in% years) %>%
      na.omit() %>%
      mutate(roll.SD.2 = roll_sd(log.recruitment, width = 2),
             roll.SD.3 = roll_sd(log.recruitment, width = 3)) %>%
      pivot_longer(c(roll.SD.2, roll.SD.3), names_to = "window", values_to = "roll.SD") -> bsai.r0.roll 
    
    ggplot(bsai.r0.roll %>% filter(TS == "bsai.opi.r0"), aes(x = Year, roll.SD, color = window))+
      facet_wrap(~TS, scales = "free_y")+
      geom_line(size = 1)+
      scale_x_continuous(breaks = seq(min(bsai.r0.roll$Year), max(bsai.r0.roll$Year), by = 2),
                         labels = seq(min(bsai.r0.roll$Year), max(bsai.r0.roll$Year), by = 2))
      geom_point(size = 1)
    