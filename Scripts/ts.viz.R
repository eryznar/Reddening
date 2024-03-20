#Timeseries processing and visualization

#PURPOSE: To load, process, and visualize GOA and EBS TS and climate timeseries

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")

### PROCESS/PLOT COMMUNITY AND SST DATA ----------------------------------------------------------------------------------------------
  # BSAI GROUNDFISH SSB -----------
    # Specify labels
      ssb.labs.bsai <- c("Aleutian Islands cod", "Aleutian Islands pollock", "Alaska plaice", "Arrowtooth flounder",
                    "Atka mackerel", "Blackspotted/rougheye rockfish", "Eastern Bering Sea cod",
                    "Eastern Bering Sea pollock", "Kamchatka flounder", "Northern rockfish", "Northern rock sole",
                    "Pacific ocean perch", "Sablefish", "Skate", "Greenland turbot", "Yellowfin sole")
      
      names(ssb.labs.bsai) <- c("bsai.ai.cod.ssb", "bsai.ai.pol.ssb", "bsai.apl.ssb", "bsai.atf.ssb", "bsai.atk.ssb",
                           "bsai.brr.ssb", "bsai.ebs.cod.ssb", "bsai.ebs.pol.ssb", "bsai.kam.ssb", "bsai.nrf.ssb",
                           "bsai.nrs.ssb", "bsai.pop.ssb", "bsai.sab.ssb", "bsai.ska.ssb", "bsai.turb.ssb", 
                           "bsai.yfs.ssb")
    # Plot and save
      ggplot(bsai.ssb %>% filter(Year > 1952), aes(x = Year, y = SSB))+
        geom_line(linewidth = 1)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.bsai), ncol = 3)+
        theme_bw()+
        ggtitle("BSAI groundfish spawning stock biomass")+
        ylab("log(kilotons)")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> bsai.ssb.plot
      
      ggsave(plot = bsai.ssb.plot, "./Figures/bsai.ssb.plot.png", width = 8.5, height = 11, units = "in")
      
    
  # GOA GROUNDFISH SSB -----------
    # Specify labels
      ssb.labs.goa <- c("Arrowtooth flounder", "Cod", "Rougheye/blackspotted rockfish", "Dover sole", "Dusky rockfish",
                    "Flathead sole", "Northern rockfish", "Northern rock sole", "Pollock", "Pacific ocean perch",
                    "Rex sole", "Sablefish", "Southern rock sole")
      
      names(ssb.labs.goa) <- c("goa.atf.ssb", "goa.cod.ssb", "goa.dbr.ssb", "goa.dov.ssb", "goa.drf.ssb", "goa.fhs.ssb",
                           "goa.nrf.ssb", "goa.nrs.ssb", "goa.pol.ssb", "goa.pop.ssb", "goa.rex.ssb", "goa.sab.ssb",
                           "goa.srs.ssb")
      
    # Plot and save
      ggplot(goa.ssb %>% filter(Year > 1959), aes(x = Year, y = SSB))+
        geom_line(linewidth = 1)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = ssb.labs.goa), ncol =3)+
        theme_bw()+
        ggtitle("GOA groundfish spawning stock biomass")+
        ylab("log(kilotons)")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10)) -> goa.ssb.plot
      
      ggsave(plot = goa.ssb.plot, "./Figures/goa.ssb.plot.png", width = 8.5, height = 11, units = "in")
      
  # SALMON CATCH ----------------
    # Specify labels
      salm.labs <- c("Bristol Bay sockeye", "Chignik chum", "Chignik pink", "Chignik sockeye", "Cook Inlet chum",
                        "Cook Inlet pink", "Cook Inlet sockeye", "Kodiak chum", "Kodiak pink", "Kodiak sockeye",
                        "Peninsula chum", "Peninsula pink", "Peninsula sockeye", "PWS chum",
                        "PWS pink", "PWS sockeye", "Southeast chum", "Southeast pink",
                        "Southeast sockeye")
      
      names(salm.labs) <- c("bsai.sck.catch", "chig.chm.catch", "chig.pnk.catch", "chig.sck.catch",
                               "cook.chm.catch", "cook.pnk.catch", "cook.sck.catch", "kod.chm.catch",
                               "kod.pnk.catch", "kod.sck.catch", "pen.chm.catch", "pen.pnk.catch",
                               "pen.sck.catch", "pws.chm.catch", "pws.pnk.catch", "pws.sck.catch", "se.chm.catch",
                               "se.pnk.catch", "se.sck.catch")
    
    # Plot and save
      ggplot(salmon.catch %>% filter(Year > 1959), aes(x = Year, y = Catch))+
        geom_line(linewidth = 0.75)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = salm.labs),
                   ncol = 3)+
        theme_bw()+
        ggtitle("BSAI and GOA salmon catch")+
        ylab("log(number of fish)")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> salmon.catch.plot
      
      ggsave(plot = salmon.catch.plot, "./Figures/salmon.catch.plot.png", width = 8.5, height = 11, units = "in")
      
  # CRAB MB ----------------
    # Plot and save
      ggplot(crab.mb %>% filter(Year > 1953), aes(x = Year, y = Value, linetype = Type))+
        geom_line(linewidth = 0.75)+
        facet_wrap(~TS, scales = "free_y", nrow = 3)+
        theme_bw()+
        scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("Female", "Male"))+
        ggtitle("BSAI crab mature biomass")+
        ylab("log(kilotons)")+
        theme(legend.position = "bottom",
              axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14),
                    strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> crab.mb.plot
      
      ggsave(plot = crab.mb.plot, "./Figures/crab.mb.plot.png", width = 8.5, height = 11, units = "in")
      
  # GROUNDFISH AND CRAB RECRUITMENT ---------------
    # Specify BSAI labels
      r0.labs.bsai <- c("Aleutian Islands cod", "Aleutian Islands pollock", "Alaska plaice", "Arrowtooth flounder",
                        "Atka mackerel", "Bristol Bay red king crab", "Blackspotted/rougheye rockfish", "Eastern Bering Sea cod",
                        "Eastern Bering Sea pollock", "Kamchatka flounder", "Northern rockfish", "Northern rock sole",
                        "Snow crab", "Pacific ocean perch", "Sablefish", "Skate", 
                        "Tanner crab", "Greenland turbot", "Yellowfin sole")
      
      names(r0.labs.bsai) <- c("bsai.ai.cod.r0", "bsai.ai.pol.r0", "bsai.apl.r0", "bsai.atf.r0", "bsai.atk.r0",
                               "bsai.bbrkc.r0", "bsai.brr.r0", "bsai.ebs.cod.r0", "bsai.ebs.pol.r0", "bsai.kam.r0", "bsai.nrf.r0",
                               "bsai.nrs.r0", "bsai.opi.r0", "bsai.pop.r0", "bsai.sab.r0", "bsai.ska.r0", 
                               "bsai.tanner.r0", "bsai.turb.r0",  "bsai.yfs.r0")
      
    # Plot and save BSAI
      ggplot(bsai.r0 %>% filter(Year > 1950), aes(x = Lagged.Year, y = Recruitment))+
        geom_line(linewidth = 0.75)+
        facet_wrap(~TS, scales = "free_y", labeller = labeller(TS = r0.labs.bsai), ncol = 3)+
        theme_bw()+
        ggtitle("BSAI groundfish and crab recruitment")+
        ylab("log(millions of recruits)") + 
        xlab("Year")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> bsai.r0.plot
      
      ggsave(plot = bsai.r0.plot, "./Figures/bsai.r0.plot.png", width = 8.5, height = 11, units = "in")
      
    # Specify GOA labels
      r0.labs.goa <- c("Arrowtooth flounder", "Cod", "Rougheye/blackspotted rockfish", "Dover sole", "Dusky rockfish",
                        "Flathead sole", "Northern rockfish", "Northern rock sole", "Pollock", "Pacific ocean perch",
                        "Rex sole", "Sablefish", "Southern rock sole")
      
      names(r0.labs.goa) <- c("goa.atf.r0", "goa.cod.r0", "goa.dbr.r0", "goa.dov.r0", "goa.drf.r0", "goa.fhs.r0",
                               "goa.nrf.r0", "goa.nrs.r0", "goa.pol.r0", "goa.pop.r0", "goa.rex.r0", "goa.sab.r0",
                               "goa.srs.r0")
      
    # Plot and save GOA
      ggplot(goa.r0 %>% filter(Year > 1959), aes(x = Lagged.Year, y = Recruitment))+
        geom_line(linewidth = 0.75)+
        facet_wrap(~TS, scales = "free_y", ncol= 3, labeller = labeller(TS=r0.labs.goa))+
        theme_bw()+
        ggtitle("BSAI groundfish and crab recruitment")+
        ylab("log(millions of recruits)") + 
        xlab("Year")+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10))-> goa.r0.plot
      
      ggsave(plot = goa.r0.plot, "./Figures/goa.r0.plot.png", width = 8.5, height = 11, units = "in")
      
  # SST ---------------
    # Plot and save
      ggplot(sst %>% filter(Year < 2024), aes(x = Year, y = mean.sst))+
        facet_wrap(~region, scales = "free_y", nrow = 2)+
        geom_line(linewidth = 0.75)+
        geom_point(size = 1.5)+
        scale_x_continuous(breaks = seq(min(sst$Year), max(sst$Year), by = 10))+
        theme_bw()+
        ggtitle("ERA5 SST")+
        ylab("°C")+
        theme(legend.position = "bottom",
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 12)) -> sst.plot
      
      ggsave(plot= sst.plot, "./Figures/ERA5sst.png", width = 6, height = 8, units = "in")
    
  # SSB/MB/R0/CATCH DATA WITH SST ------------
  
    # Plot and save BSAI SSB/SST
      ggplot(bsai.ssb.sst %>% filter(Year > 1952), aes(x = mean.sst, y = SSB))+
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
      ggplot(goa.ssb.sst %>% filter(Year > 1959), aes(x = mean.sst, y = SSB))+
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
      ggplot(salmon.catch.sst, aes(x = Lagged.sst, y = Catch))+
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
      ggplot(crab.mb.sst, aes(x = mean.sst, y = Value, color = Type))+
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
      ggplot(bsai.r0.sst, aes(x = Lagged.sst, y = Recruitment))+
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
      ggplot(goa.r0.sst, aes(x = Lagged.sst, y = Recruitment))+
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
      
    