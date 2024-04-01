#Timeseries processing and visualization

#PURPOSE: To load, process, and visualize GOA and EBS TS and climate timeseries

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")

### PROCESS/PLOT COMMUNITY AND SST DATA ----------------------------------------------------------------------------------------------
  # BSAI GROUNDFISH SSB -----------
   
    # Plot and save
    ggplot(bsai.ssb, aes(x = Year, y = log.SSB))+
        geom_line(linewidth = 1, color = "#6A6DB7")+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.bsai), ncol = 3)+
        theme_bw()+
        scale_x_continuous(breaks = seq(min(bsai.ssb$Year), max(bsai.ssb$Year), by = 15))+
        ggtitle("BSAI groundfish spawning stock biomass")+
        ylab("log(kilotons)")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.position = "none",
              strip.text = element_text(size = 10))-> bsai.ssb.plot
      
      ggsave(plot = bsai.ssb.plot, "./Figures/bsai.ssb.plot.png", width = 11, height = 8.5, units = "in")
      
    
  # GOA GROUNDFISH SSB -----------
    # Plot and save
      ggplot(goa.ssb, aes(x = Year, y = log.SSB))+
      #ggplot(goa.ssb %>% filter(Year > 1959), aes(x = Year, y = log.SSB))+
        geom_line(linewidth = 1, color = "#A34242")+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.goa), ncol =3)+
        theme_bw()+
        scale_x_continuous(breaks = seq(min(goa.ssb$Year), max(goa.ssb$Year), by = 15))+
        ggtitle("GOA groundfish spawning stock biomass")+
        ylab("log(kilotons)")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.position = "none",
              strip.text = element_text(size = 10)) -> goa.ssb.plot
      
      ggsave(plot = goa.ssb.plot, "./Figures/goa.ssb.plot.png", width = 11, height = 8.5, units = "in")
      
  # SALMON CATCH ----------------
  
    # Plot and save
      ggplot(salmon.catch, aes(x = Year, y = log.catch, color = Region))+
      #ggplot(salmon.catch %>% filter(Year > 1959), aes(x = Year, y = log.catch))+
        geom_line(linewidth = 1)+
        scale_color_manual(values = c("#6A6DB7", "#A34242"))+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = salm.labs),
                   ncol = 3)+
        scale_x_continuous(breaks = seq(min(salmon.catch$Year), max(salmon.catch$Year), by = 30))+
        theme_bw()+
        ggtitle("BSAI and GOA salmon catch")+
        ylab("log(number of fish)")+
        theme(axis.text = element_text(size = 12),
              legend.position  = "none",
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> salmon.catch.plot
      
      ggsave(plot = salmon.catch.plot, "./Figures/salmon.catch.plot.png", width = 11, height = 8.5, units = "in")
      
  # CRAB MB ----------------
    # Plot and save
      ggplot(crab.mb, aes(x = Year, y = log.value, linetype = Type))+
      #ggplot(crab.mb %>% filter(Year > 1953), aes(x = Year, y = log.value, linetype = Type))+
        geom_line(linewidth = 1, color = "#6A6DB7")+
        facet_wrap(~TS, scales = "free_y", nrow = 3)+
        theme_bw()+
        scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("Female", "Male"))+
        ggtitle("BSAI crab mature biomass")+
        ylab("log(kilotons)")+
        scale_x_continuous(breaks = seq(min(crab.mb$Year), max(crab.mb$Year), by = 5))+
        theme(legend.position = "bottom",
              axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14),
                    strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> crab.mb.plot
      
      ggsave(plot = crab.mb.plot, "./Figures/crab.mb.plot.png", width = 11, height = 8.5, units = "in")
      
  # GROUNDFISH AND CRAB RECRUITMENT ---------------
   
    # Plot and save BSAI
      ggplot(bsai.r0, aes(x = Year, y = log.recruitment))+
      #ggplot(bsai.r0 %>% filter(Year ==), aes(x = Lagged.Year, y = log.recruitment))+
        geom_line(linewidth = 1, color = "#6A6DB7")+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.bsai), ncol = 3)+
        theme_bw()+
        ggtitle("BSAI groundfish and crab recruitment")+
        scale_x_continuous(breaks = seq(min(bsai.r0$Year), max(bsai.r0$Year), by = 15))+
        ylab("log(millions of recruits)") + 
        xlab("Year")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> bsai.r0.plot
      
      ggsave(plot = bsai.r0.plot, "./Figures/bsai.r0.plot.png", width = 11, height = 8.5, units = "in")
      
  
      
    # Plot and save GOA
      ggplot(goa.r0, aes(x = Year, y = log.recruitment))+
      #ggplot(goa.r0 %>% filter(Year > 1959), aes(x = Lagged.Year, y = log.recruitment))+
        geom_line(linewidth = 1, color = "#A34242")+
        facet_wrap(~TS, scales = "free_y", ncol= 3, labeller = labeller(TS=r0.labs.goa))+
        theme_bw()+
        ggtitle("GOA groundfish recruitment")+
        scale_x_continuous(breaks = seq(min(goa.r0$Year), max(goa.r0$Year), by = 15))+
        ylab("log(millions of recruits)") + 
        xlab("Year")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> goa.r0.plot
      
      ggsave(plot = goa.r0.plot, "./Figures/goa.r0.plot.png", width = 11, height = 8.5, units = "in")
      
  # SST ---------------
    # Plot and save
      ggplot(sst %>% filter(Year < 2024), aes(x = Year, y = mean.sst, color = region))+
        facet_wrap(~region, scales = "free_y", nrow = 2)+
        scale_color_manual(values = c("#6A6DB7", "#A34242"))+
        geom_point(size = 2)+
        geom_line(size = 1.5)+
        scale_x_continuous(breaks = seq(min(sst$Year), max(sst$Year), by = 10))+
        theme_bw()+
        ggtitle("ERA5 SST")+
        ylab("°C")+
        theme(legend.position = "none",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16),
              legend.text = element_text(size = 14),
              title = element_text(size = 16)) -> sst.plot
      
      ggsave(plot= sst.plot, "./Figures/ERA5sst.png", width = 11, height = 8.5, units = "in")
    
  # SSB/MB/R0/CATCH DATA WITH SST ------------
  
    # Plot and save BSAI SSB/SST
      ggplot(bsai.ssb.sst %>% filter(Year > 1952), aes(x = mean.sst, y = log.SSB))+
        geom_point()+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = ssb.labs.bsai),
                   ncol =3)+
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
        theme_bw()+
        ggtitle("BSAI groundfish spawning stock biomass and SST")+
        ylab("log(kilotons)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> bsai.ssb.sst.plot
      
      ggsave(plot= bsai.ssb.sst.plot, "./Figures/bsai.ssb.sst.plot.png", height = 11, width = 8.5, units = "in")
    
    # Plot and save GOA SSB/SST
      ggplot(goa.ssb.sst %>% filter(Year > 1959), aes(x = mean.sst, y = log.SSB))+
        geom_point()+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = ssb.labs.goa),
                   ncol =3)+
        theme_bw()+
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
        ggtitle("GOA groundfish spawning stock biomass and SST")+
        ylab("log(kilotons)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> goa.ssb.sst.plot
      
      ggsave(plot= goa.ssb.sst.plot, "./Figures/goa.ssb.sst.plot.png", height = 11, width = 8.5, units = "in")
   
    # Plot and save salmon catch/SST
      ggplot(salmon.catch.sst, aes(x = Lagged.sst, y = log.catch))+
        geom_point()+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = salm.labs),
                   ncol = 3)+
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
        theme_bw()+
        ggtitle("BSAI and GOA salmon catch")+
        ylab("log(number of fish)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> salmon.catch.sst.plot
      
      ggsave(plot = salmon.catch.sst.plot, "./Figures/salmon.catch.sst.plot.png", width = 8.5, height = 11, units = "in")
      
    # Plot and save crab mb/SST
      ggplot(crab.mb.sst, aes(x = mean.sst, y = log.value, color = Type))+
        geom_point(size = 2)+
        facet_wrap(~TS, scales = "free", nrow = 3)+
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
        theme_bw()+
        ggtitle("BSAI crab mature biomass and SST")+
        scale_color_manual(values = c("salmon", "turquoise"), labels = c("Female", "Male"),
                           name = "")+
        ylab("log(kilotons)")+
        xlab("°C")+
        theme_bw()+
        theme(legend.position = "bottom",
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> crab.mb.sst.plot
      
      ggsave(plot = crab.mb.sst.plot, "./Figures/crab.mb.sst.plot.png", width = 8.5, height = 11, units = "in")

    # Plot BSAI r0/SST
      ggplot(bsai.r0.sst, aes(x = Lagged.sst, y = log.recruitment))+
        geom_point()+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = r0.labs.bsai),
                   ncol = 3)+
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
        theme_bw()+
        ggtitle("BSAI groundfish/crab recruitment and SST")+
        ylab("log(millions of recruits)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> bsai.r0.sst.plot
      
      ggsave(plot = bsai.r0.sst.plot, "./Figures/bsai.r0.sst.plot.png", width = 8.5, height = 11, units = "in")
      
    
    # Plot GOA r0/SST
      ggplot(goa.r0.sst, aes(x = Lagged.sst, y = log.recruitment))+
        geom_point()+
        facet_wrap(~TS, scales = "free", labeller = labeller(TS = r0.labs.goa),
                   ncol = 3)+
        theme_bw()+
        ggtitle("GOA groundfish recruitment and SST")+
        geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
        ylab("log(millions of recruits)")+
        xlab("°C")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> goa.r0.sst.plot
      
      ggsave(plot = goa.r0.sst.plot, "./Figures/goa.r0.sst.plot.png", width = 8.5, height = 11, units = "in")
      
    