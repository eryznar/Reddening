# LOAD SOURCE FILE --------------------------------------------------------------
source("./Scripts/ts.processing.R")
source("./Scripts/load.libs.functions.R")


# LOAD DATA ---------------------------------------------------------------------
sst.lme <- read.csv("./Output/lme.sst_df.csv") 

# PLOT DATA ---------------------------------------------------------------------
ggplot(sst.lme, aes(year, mean.sst))+
  geom_line(size = 1)+
  facet_wrap(~LME, scales = "free_y",  labeller = label_wrap_gen())+
  theme_bw()+
  ylab("°C")+
  xlab("Year")

ggsave("./Figures/lme.sst.png", height = 8.5, width = 11, units = "in")


# CALCULATE AR1 AND SD ----------------------------------------------------------
lme <- unique(sst.lme$LME)
sum.out <- data.frame()

for(ii in 1:length(lme)){
  sst.lme %>%
    filter(LME == lme[ii]) %>%
    na.omit() -> TS.dat
  
  # Specify sliding window width
  width = 15
  
  # Calculate rolling window AR1
  ar1 <- sapply(rollapply(TS.dat$mean.sst, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
  
  # Calculate rolling window SD
  sd <-  rollapply(TS.dat$mean.sst, width = width, FUN = sd)
  
  # Calculate windows
  window <- seq(min(TS.dat$year)+(width-1), max(TS.dat$year), by = 1)
  
  
  # Compile output
  sum.out <- rbind(sum.out, data.frame(LME = lme[ii],
                                       ar1 = ar1,
                                       sd = sd,
                                       window = window))
}


# PLOT AR1 AND SD ---------------------------------------------------------------
ggplot(sum.out, aes(window, ar1))+
  geom_line(size = 1)+
  facet_wrap(~LME, scales = "free_y",  labeller = label_wrap_gen())+
  theme_bw()+
  ylab("AR1")+
  xlab("Window")

ggsave("./Figures/lme.sstAR1.png", height = 8.5, width = 11, units = "in")

ggplot(sum.out, aes(window, sd))+
  geom_line(size = 1)+
  facet_wrap(~LME, scales = "free_y",  labeller = label_wrap_gen())+
  theme_bw()+
  ylab("SD")+
  xlab("Window")

ggsave("./Figures/lme.sstSD.png", height = 8.5, width = 11, units = "in")
