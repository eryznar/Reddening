#GAMs on detrended sst and biological ts

#Purpose is to assess the following questions: 
# 1) Has the climate gotten redder? 
# 2) Has the biology gotten redder?
# To do so, this script will test for trends in sst AR1 and SD and recruitment AR1 and CV
# by fitting GAMs on detrended data

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")

source("./Scripts/functions.R")

### Test question 1: GAMs to predict sst AR1 and SD with time -----------------------------------------------------------
  # Detrend data
  trend.fun(ebs.sst, "SST") -> ebs.sst.out
  
  trend.fun(goa.sst, "SST") -> goa.sst.out
  
  # Fit and select best models
  model.out <- data.frame()
  
  assess.trend(ebs.sst.out, "sst") -> ebs.sst.best.mods
  
  model.out <- data.frame()
  
  assess.trend(goa.sst.out, "sst") -> goa.sst.best.mods
  
  # Predict with best models
  pred.vals <- data.frame()
  
  model.predict(ebs.sst.best.mods, ebs.sst.out, "sst") -> ebs.pred.out
  
  pred.vals <- data.frame()
  
  model.predict(goa.sst.best.mods, goa.sst.out, "sst") -> goa.pred.out
  
  # Bind
  rbind(ebs.pred.out %>% mutate(region = "Eastern Bering Sea"),
        goa.pred.out %>% mutate(region = "Gulf of Alaska")) -> plot.dat.sst
  
  # Plot sst AR1 with time
   plot.dat.sst %>%
    dplyr::select(TS, response, region, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") -> lab.dat
  
  labs <- paste0(lab.dat$region, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- c("Eastern Bering Sea", "Gulf of Alaska")
  
  plot.dat.sst %>%
    filter(response == "ar1") -> plot.dat.sst2
  
  ggplot()+
    geom_ribbon(plot.dat.sst2, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat.sst2, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat.sst2, mapping = aes(x = window, y = pred, color = region), size = 1.25)+
    scale_color_manual(values = c("#6A6DB7", "#A34242"))+
    scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
    facet_wrap(~region, scales = "free_y", labeller = labeller(region = labs),
               ncol = 1)+
    theme_bw()+
    ggtitle("Detrended SST AR1 with time")+
    ylab("AR1")+
    xlab("Window")+
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16)) -> sst.ar1.window.plot
  
  ggsave(plot = sst.ar1.window.plot, "./Figures/sst.AR1.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # Plot sst SD with time
  plot.dat.sst %>%
    dplyr::select(TS, response, region, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "var.val") -> lab.dat
  
  labs <- paste0(lab.dat$region, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- c("Eastern Bering Sea", "Gulf of Alaska")
  
  plot.dat.sst %>%
    filter(response == "var.val") -> plot.dat.sst2
  
  ggplot()+
    geom_ribbon(plot.dat.sst2, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat.sst2, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat.sst2, mapping = aes(x = window, y = pred, color = region), size = 1.25)+
    scale_color_manual(values = c("#6A6DB7", "#A34242"))+
    scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
    facet_wrap(~region, scales = "free_y", labeller = labeller(region = labs),
               ncol = 1)+
    theme_bw()+
    ggtitle("Detrended SST SD with time")+
    ylab("SD")+
    xlab("Window")+
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16)) -> sst.SD.window.plot
  
  ggsave(plot = sst.SD.window.plot, "./Figures/sst.SD.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  
### Test question 2: GAMs to predict BSAI and GOA AR1 and CV with time ----------------------------------------------------
  # EBS ----
  # Detrend data
  sum.out <- data.frame()
  
  trend.fun(bsai.r0, "Recruitment") -> ebs.r0.out
  
  # Fit and select best models
  model.out <- data.frame()
  
  assess.trend(ebs.r0.out, "Recruitment") -> ebs.r0.best.mods
  
  # Predict with best models
  pred.vals <- data.frame()
  
  model.predict(ebs.r0.best.mods, ebs.r0.out, "Recruitment") -> ebs.r0.pred.out
  
 # Plot r0 AR1 with time
  ebs.r0.pred.out %>%
    dplyr::select(TS, response, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") -> lab.dat
  
  labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- names(r0.labs.bsai)
  
  ebs.r0.pred.out %>%
    filter(response == "ar1") -> plot.dat
  
  ggplot()+
    geom_ribbon(plot.dat, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat, mapping = aes(x = window, y = pred), color =  "#6A6DB7", size = 1.25)+
    facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("EBS detrended recruitment AR1 with time")+
    ylab("AR1")+
    xlab("Window")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> EBS.r0.ar1.window.plot
  
  ggsave(plot = EBS.r0.ar1.window.plot, "./Figures/EBS.r0.AR1.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # Plot r0 CV with time
  ebs.r0.pred.out %>%
    dplyr::select(TS, response, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "var.val") -> lab.dat
  
  labs <- paste0(r0.labs.bsai, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- names(r0.labs.bsai)
  
  ebs.r0.pred.out %>%
    filter(response == "var.val") -> plot.dat
  
  ggplot()+
    geom_ribbon(plot.dat, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat, mapping = aes(x = window, y = pred), color =  "#6A6DB7", size = 1.25)+
    facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("EBS detrended recruitment CV with time")+
    ylab("CV")+
    xlab("Window")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> EBS.r0.CV.window.plot
  
  ggsave(plot = EBS.r0.CV.window.plot, "./Figures/EBS.r0.CV.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # GOA ----
  # Detrend data
  sum.out <- data.frame()
  
  trend.fun(goa.r0, "Recruitment") -> goa.r0.out
  
  # Fit and select best models
  model.out <- data.frame()
  
  assess.trend(goa.r0.out, "Recruitment") -> goa.r0.best.mods
  
  # Predict with best models
  pred.vals <- data.frame()
  
  model.predict(goa.r0.best.mods, goa.r0.out, "Recruitment") -> goa.r0.pred.out
  
  # Plot r0 AR1 with time
  goa.r0.pred.out %>%
    dplyr::select(TS, response, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") -> lab.dat
  
  labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- names(r0.labs.goa)
  
  goa.r0.pred.out %>%
    filter(response == "ar1") -> plot.dat
  
  ggplot()+
    geom_ribbon(plot.dat, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat, mapping = aes(x = window, y = pred), color =  "#A34242", size = 1.25)+
    facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("GOA detrended recruitment AR1 with time")+
    ylab("AR1")+
    xlab("Window")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> GOA.r0.ar1.window.plot
  
  ggsave(plot = GOA.r0.ar1.window.plot, "./Figures/GOA.r0.AR1.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # Plot r0 AR1 with time
  goa.r0.pred.out %>%
    dplyr::select(TS, response, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") -> lab.dat
  
  labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- names(r0.labs.goa)
  
  goa.r0.pred.out %>%
    filter(response == "ar1") -> plot.dat
  
  ggplot()+
    geom_ribbon(plot.dat, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat, mapping = aes(x = window, y = pred), color =  "#A34242", size = 1.25)+
    facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("GOA detrended recruitment AR1 with time")+
    ylab("AR1")+
    xlab("Window")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> GOA.r0.ar1.window.plot
  
  ggsave(plot = GOA.r0.ar1.window.plot, "./Figures/GOA.r0.AR1.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # Plot r0 AR1 with time
  goa.r0.pred.out %>%
    dplyr::select(TS, response, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "ar1") -> lab.dat
  
  labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- names(r0.labs.goa)
  
  goa.r0.pred.out %>%
    filter(response == "ar1") -> plot.dat
  
  ggplot()+
    geom_ribbon(plot.dat, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat, mapping = aes(x = window, y = pred), color =  "#A34242", size = 1.25)+
    facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("GOA detrended recruitment AR1 with time")+
    ylab("AR1")+
    xlab("Window")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> GOA.r0.ar1.window.plot
  
  ggsave(plot = GOA.r0.ar1.window.plot, "./Figures/GOA.r0.AR1.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  # Plot r0 SD with time
  goa.r0.pred.out %>%
    dplyr::select(TS, response, k, rsq)%>%
    distinct() %>%
    arrange(., TS) %>%
    filter(response == "var.val") -> lab.dat
  
  labs <- paste0(r0.labs.goa, " \n(k=", lab.dat$k, " , R2=", round(lab.dat$rsq, 2), ")")
  names(labs) <- names(r0.labs.goa)
  
  goa.r0.pred.out %>%
    filter(response == "var.val") -> plot.dat
  
  ggplot()+
    geom_ribbon(plot.dat, 
                mapping = aes(x = window, ymin = pred - pred.CI, ymax= pred +pred.CI), fill = "grey", alpha = 0.5)+
    #geom_point(plot.dat, mapping=aes(x = window, y = observed), color = "black")+
    geom_line(plot.dat, mapping = aes(x = window, y = pred), color =  "#A34242", size = 1.25)+
    facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = labs),
               ncol = 3)+
    theme_bw()+
    ggtitle("GOA detrended recruitment CV with time")+
    ylab("CV")+
    xlab("Window")+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 12)) -> GOA.r0.CV.window.plot
  
  ggsave(plot = GOA.r0.CV.window.plot, "./Figures/GOA.r0.CV.x.time.GAM.png", width = 8.5, height = 11, units = "in")
  
  