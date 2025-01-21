# Script to calculate AR1 and SD for monthly SST anomalies, PDO index values, and SLP anomalies
source("./Scripts/load.libs.functions.R")

# AR1 and SD for monthly SST anomalies ----
# Load 
ebs.sst <- read.csv("./Data/ebs.monthlySSTanomalies.csv") %>%
  filter(Year %in% 1948:2024) %>%
  mutate(dec.yr = Year + (Month - 0.5)/12)
goa.sst <- read.csv("./Data/goa.monthlySSTanomalies.csv")%>%
  filter(Year %in% 1948:2024) %>%
  mutate(dec.yr = Year + (Month - 0.5)/12)

rbind(ebs.sst %>% mutate(region = "Eastern Bering Sea"),
      goa.sst %>% mutate(region = "Gulf of Alaska")) -> sst


width = nrow(ebs.sst)/2 # per Boulton and Lenton

  # EBS SST ----
  # Detrend data
  detrend.dat <- lm(month.anom ~ dec.yr, ebs.sst)
  
  # Extract residuals
  resid <- data.frame(Year = ebs.sst$Year, month = ebs.sst$Month, month.anom = detrend.dat$residuals)
  
  # Calculate rolling window AR1
  ar1 <- sapply(rollapply(resid$month.anom, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
  
  # Calculate rolling window SD
  sd <-  rollapply(resid$month.anom, width = width, FUN = sd, fill = NA)
  
  
  # Make data frame of sd-cv
  data.frame(year = ebs.sst$Year, month = ebs.sst$Month, dec.yr = ebs.sst$dec.yr, sd = sd) -> win.dat
  
  # Calculate windows
  win.yr <- na.omit(win.dat) %>% pull(year)
  win.month <- na.omit(win.dat) %>% pull(month)
  
  # Make data frame of ar1
  data.frame(year = win.yr, month = win.month, ar1 = ar1) -> ar1.dat
  
  # Join
  left_join(win.dat, ar1.dat) -> win.dat
  
  plot(win.dat$ar1, type = "l")
  plot(win.dat$sd, type = "l")
  
  win.dat %>%
    mutate(date = paste0(year, "-", month)) -> win.dat.ebs
  
  tau.ar1 <- cor.test(win.dat$ar1, win.dat$dec.yr, method = "kendall")
  
  ggplot(win.dat.ebs, aes(dec.yr, ar1))+
    geom_line(group = 1)+
    ggtitle("EBS AR1")+
    xlab("Month")+
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16))
  
  ggplot(win.dat.ebs, aes(dec.yr, sd))+
    geom_line(group = 1)+
    ggtitle("EBS SD")+
    xlab("Month")+
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16))
  
  # GOA SST ----
  # Detrend data
  detrend.dat <- lm(month.anom ~ dec.yr, goa.sst)
  
  # Extract residuals
  resid <- data.frame(Year = goa.sst$Year, month = goa.sst$Month, month.anom = detrend.dat$residuals)
  
  # Calculate rolling window AR1
  ar1 <- sapply(rollapply(resid$month.anom, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
  
  # Calculate rolling window SD
  sd <-  rollapply(resid$month.anom, width = width, FUN = sd, fill = NA)
  
  
  # Make data frame of sd-cv
  data.frame(year = goa.sst$Year, month = goa.sst$Month, sd = sd, dec.yr = goa.sst$dec.yr) -> win.dat
  
  # Calculate windows
  win.yr <- na.omit(win.dat) %>% pull(year)
  win.month <- na.omit(win.dat) %>% pull(month)
  
  # Make data frame of ar1
  data.frame(year = win.yr, month = win.month, ar1 = ar1) -> ar1.dat
  
  # Join
  left_join(win.dat, ar1.dat) -> win.dat
  
  plot(win.dat$ar1, type = "l")
  plot(win.dat$sd, type = "l")
  
  win.dat %>%
    mutate(date = paste0(year, "-", month)) -> win.dat.goa
  
  tau.ar1 <- cor.test(win.dat$ar1, win.dat$dec.yr, method = "kendall")
  tau.sd <- cor.test(win.dat$sd, win.dat$dec.yr, method = "kendall")
  
  
  ggplot(win.dat.goa, aes(dec.yr, ar1))+
    geom_line(group = 1)+
    ggtitle("GOA AR1")+
    xlab("Month")+
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          title = element_text(size = 16))
    
    ggplot(win.dat.goa, aes(1:nrow(win.dat.goa), sd))+
      geom_line(group = 1)+
      ggtitle("GOA SD")+
    xlab("Month")+
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 14),
            title = element_text(size = 16))
    
    

# AR1 and SD for monthly PDO index values ----
    pdo <- read.csv("./Data/pdo.timeseries.ersstv5.csv")
    names(pdo) <- c("Date", "index")
    
    pdo %>%
      mutate(year = as.numeric(substr(Date, 1, 4)),
             month = as.numeric(substr(Date, 6, 7)),
             dec.yr = year + (month - 0.5)/12) %>%
      filter(year %in% 1948:2024, index > -2000)-> pdo2 
    
    
    plot(pdo2$dec.yr, pdo2$index, type = "l")
    
    width = 460 # per Boulton and Lenton
    
    
    # Calculate rolling window AR1
    ar1 <- sapply(rollapply(pdo2$index, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
    
    # Calculate rolling window SD
    sd <-  rollapply(pdo2$index, width = width, FUN = sd, fill = NA)
    
    
    # Make data frame of sd-cv
    data.frame(year = pdo2$year, month = pdo2$month, dec.yr = pdo2$dec.yr, sd = sd) -> win.dat
    
    # Calculate windows
    win.yr <- na.omit(win.dat) %>% pull(year)
    win.month <- na.omit(win.dat) %>% pull(month)
    
    # Make data frame of ar1
    data.frame(year = win.yr, month = win.month, ar1 = ar1) -> ar1.dat
    
    # Join
    left_join(win.dat, ar1.dat, relationship = "many-to-many") -> win.dat.pdo
    
    ggplot(win.dat.pdo, aes(dec.yr, ar1))+
      ggtitle("PDO")+
      geom_line()+
      xlab("Year")+
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 14),
            title = element_text(size = 16))
    
    ggplot(win.dat.pdo, aes(dec.yr, sd))+
      ggtitle("PDO SD")+
      geom_line()+
      xlab("Year")+
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 14),
            title = element_text(size = 16))
    
# AR1 and SD for monthly SLP anomalies ----
# Load 
  slp <- read.csv("./Data/monthlySLPanomalies.csv") %>%
      filter(Year %in% 1948:2024) %>%
      mutate(dec.yr = Year + (Month - 0.5)/12)
    
    width = 460 # per Boulton and Lenton
    
    
    # Calculate rolling window AR1
    ar1 <- sapply(rollapply(slp$month.anom, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
    
    # Calculate rolling window SD
    sd <-  rollapply(slp$month.anom, width = width, FUN = sd, fill = NA)
    
    
    # Make data frame of sd-cv
    data.frame(year = slp$Year, month = slp$month, dec.yr = slp$dec.yr, sd = sd) -> win.dat
    
    # Calculate windows
    win.yr <- na.omit(win.dat) %>% pull(year)
    win.month <- na.omit(win.dat) %>% pull(month)
    
    # Make data frame of ar1
    data.frame(year = win.yr, month = win.month, ar1 = ar1) -> ar1.dat
    
    # Join
    left_join(win.dat, ar1.dat, relationship = "many-to-many") -> win.dat.slp
    
    ggplot(win.dat.slp, aes(dec.yr, ar1))+
      ggtitle("SLP AR1")+
      geom_line()+
      xlab("Year")+
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 14),
            title = element_text(size = 16))
    
    ggplot(win.dat.slp, aes(dec.yr, sd))+
      ggtitle("SLP SD")+
      geom_line()+
      xlab("Year")+
      geom_vline(xintercept = 1988.5, linetype = "dashed")+
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 14),
            title = element_text(size = 16))
  
  
  
# Compare all values ----
rbind(win.dat.ebs %>% mutate(type = "ebs sst") %>% dplyr::select(!date), 
      win.dat.goa %>% mutate(type = "goa sst") %>% dplyr::select(!date), 
      win.dat.pdo %>% mutate(type = "pdo"),
      win.dat.slp %>% mutate(type = "slp")) -> plot.dat

plot.dat %>% filter(type == "ebs sst") -> v1
plot.dat %>% filter(type == "pdo") -> v2

cor(na.omit(scale(v1$ar1))[1:461], na.omit(scale(v2$sd))[1:461])

ggplot()+
  geom_line(plot.dat %>% filter(type == "ebs sst"),  mapping = aes(dec.yr, scale(ar1), color = "salmon"))+
  geom_line(plot.dat %>% filter(type == "pdo"),  mapping = aes(dec.yr, scale(sd), color = "steelblue"))+
  scale_color_manual(values = c("salmon", "steelblue"), labels = c("EBS SST AR1", "PDO SD"), name = "")+
  ylab("Value")+
  geom_vline(xintercept = 1988.5, linetype = "dashed") -> plot.1
  

plot.dat %>% filter(type == "ebs sst") -> v1
plot.dat %>% filter(type == "pdo") -> v2

cor(na.omit(scale(v1$ar1))[1:461], na.omit(scale(v2$ar1))[1:461])

ggplot()+
  geom_line(plot.dat %>% filter(type == "ebs sst"), mapping = aes(dec.yr, scale(ar1), color = "salmon"))+
  geom_line(plot.dat %>% filter(type == "pdo"),  mapping = aes(dec.yr, scale(ar1), color = "steelblue"))+
  scale_color_manual(values = c("salmon", "steelblue"), labels = c("EBS SST AR1", "PDO AR1"), name = "")+
  ylab("Value")+
  geom_vline(xintercept = 1988.5, linetype = "dashed") -> plot.2



plot.dat %>% filter(type == "goa sst") -> v1
plot.dat %>% filter(type == "pdo") -> v2

cor(na.omit(scale(v1$ar1))[1:461], na.omit(scale(v2$sd))[1:461])

ggplot()+
  geom_line(plot.dat %>% filter(type == "goa sst"),  mapping = aes(dec.yr, scale(ar1), color = "salmon"))+
  geom_line(plot.dat %>% filter(type == "pdo"),  mapping = aes(dec.yr, scale(sd), color = "steelblue"))+
  scale_color_manual(values = c("goldenrod", "steelblue"), labels = c("GOA SST AR1", "PDO SD"), name = "")+
  ylab("Value")+
  geom_vline(xintercept = 1988.5, linetype = "dashed") -> plot.3



plot.dat %>% filter(type == "goa sst") -> v1
plot.dat %>% filter(type == "pdo") -> v2

cor(na.omit(scale(v1$ar1))[1:461], na.omit(scale(v2$ar1))[1:461])

ggplot()+
  geom_line(plot.dat %>% filter(type == "goa sst"), mapping = aes(dec.yr, scale(ar1), color = "salmon"))+
  geom_line(plot.dat %>% filter(type == "pdo"), mapping = aes(dec.yr, scale(ar1), color = "steelblue"))+
  scale_color_manual(values = c("goldenrod", "steelblue"), labels = c("GOA SST AR1", "PDO AR1"), name = "")+
  ylab("Value")+
  geom_vline(xintercept = 1988.5, linetype = "dashed") -> plot.4

cowplot::plot_grid(plot.1, plot.2, plot.3, plot.4, nrow = 2)
    