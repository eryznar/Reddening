#Model exploration

#PURPOSE: To fit models GOA and EBS stock and climate timeseries

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
       
        for(ii in 1:length(unique(dat$Stock))){
          # filter data by stock
          dat %>%
            filter(Stock == unique(dat$Stock)[ii]) -> stock.dat
          
          # fit gams, record pval and AIC
          gam.1 <- gam(Recruitment ~ s(sst, k = 3), data= stock.dat, method = "REML")
          gam.2 <- gam(Recruitment ~ s(sst.2, k = 3), data= stock.dat, method = "REML")
          gam.3 <- gam(Recruitment ~ s(sst.3, k = 3), data= stock.dat, method = "REML")
          
          p.gam.1 <- signif(summary(gam.1)$s.table[,4], 2)
          p.gam.2 <- signif(summary(gam.2)$s.table[,4],2)
          p.gam.3 <- signif(summary(gam.3)$s.table[,4],2)
          
          AIC.gam.1 <- AIC(gam.1)
          AIC.gam.2 <- AIC(gam.2)
          AIC.gam.3 <- AIC(gam.3)
          
          # fit gls, record pval and AIC
          gls.1 <- gls(Recruitment ~ sst, data = stock.dat, correlation = corAR1())
          gls.2 <- gls(Recruitment ~ sst.2, data = stock.dat, correlation = corAR1())
          gls.3 <- gls(Recruitment ~ sst.3, data = stock.dat, correlation = corAR1())
          
          p.gls.1 <- signif(summary(gls.1)$tTable[2,4], 2)
          p.gls.2 <- signif(summary(gls.2)$tTable[2,4], 2)
          p.gls.3 <- signif(summary(gls.3)$tTable[2,4], 2)
          
          AIC.gls.1 <- AIC(gls.1)
          AIC.gls.2 <- AIC(gls.2)
          AIC.gls.3 <- AIC(gls.3)
          
          # Build summary table
          model.out <- rbind(model.out, data.frame(Stock = unique(dat$Stock)[ii],
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
      
      # Fit models for GOA
        model.out <- data.frame()
        
        fit.models(dat = goa.r0.sst) -> goa.r0.out
        
          goa.r0.out %>%
            group_by(Model) %>%
            reframe(mean.AIC = mean(AIC)) -> goa.sum
    
  # Calculate AR1 and SD for 15 year windows for biology/SST timeseries ----
    # Specify 15-year windows
    win.1 <- seq(min(years), min(years)+14, by= 1)
    win.2 <- seq(max(win.1)+1, max(win.1)+15)
    
    
    # Create function to calculate values
    sum.fun <- function(dat, data.type){
      
      if(data.type == "SSB"){
        dat %>%
          rename(Value = SSB) -> dat
      }else if(data.type == "Recruitment"){
        dat %>%
          rename(Value = Recruitment) -> dat
      }else if(data.type == "Catch"){
        dat %>% 
          rename(Value = Catch) -> dat
      }else if(data.type == "SST"){
        dat %>%
          rename(Value = mean.sst) -> dat
      }else{
        dat = dat
      }
      
      for(ii in 1:length(unique(dat$Stock))){
        dat %>%
          filter(Stock == unique(dat$Stock)[ii], Year %in% years) %>%
          na.omit() -> stock.dat
        
        win1.dat <- stock.dat %>%
                      filter(Year %in% win.1)
        
        win2.dat <- stock.dat %>%
                      filter(Year %in% win.2)
        
        # Calculate AR1
        ar1.win1 <- acf(win1.dat$Value)$acf[2]
        ar1.win2 <- acf(win2.dat$Value)$acf[2]
        
        # Calculate SD
        sd.win1 <- sd(win1.dat$Value)
        sd.win2 <- sd(win2.dat$Value)
        
        # Compile output
        sum.out <- rbind(sum.out, data.frame(Stock = unique(dat$Stock)[ii],
                                             ar1.win1 = ar1.win1, 
                                             ar1.win2 = ar1.win2,
                                             sd.win1 = sd.win1, 
                                             sd.win2 = sd.win2))
        
      }
      return(sum.out)
    }
    
    # BSAI SSB
    sum.out <- data.frame()
    
    sum.fun(bsai.ssb, "SSB") -> sum.ssb.bsai 
    
    # GOA SSB
    sum.out <- data.frame()
    
    sum.fun(goa.ssb, "SSB") -> sum.ssb.goa
    
    # BSAI R0
    sum.out <- data.frame()
    
    sum.fun(bsai.r0, "Recruitment") -> sum.r0.bsai
    
    # GOA R0
    sum.out <- data.frame()
    
    sum.fun(goa.r0, "Recruitment") -> sum.r0.goa
    
    # CRAB MB
    sum.out <- data.frame()
    
    sum.fun(crab.mb, "Mature biomass") -> sum.mb.crab
    
    # BSAI salmon
    sum.out <- data.frame()
    
    sum.fun(bsai.salmon, "Catch") -> sum.catch.bsai
    
    # GOA salmon
    sum.out <- data.frame()
    
    sum.fun(goa.salmon, "Catch") -> sum.catch.goa
    
    # EBS SST
    sum.out <- data.frame()
    
    sum.fun(ebs.sst, "SST") -> sum.sst.bsai
    
    