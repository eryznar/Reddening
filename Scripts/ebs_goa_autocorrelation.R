# ACF for GOA/EBS SST

source("./Scripts/load.libs.functions.R")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set theme
theme_set(theme_bw())

goa_sst <- read.csv(paste0(dir, "Data/goa.monthlySSTanomalies.csv"), row.names = 1) %>%
  filter(Year %in% 1948:2024) 

ebs_sst <- read.csv(paste0(dir, "Data/ebs.monthlySSTanomalies.csv"), row.names = 1) %>%
  filter(Year %in% 1948:2024) 

goa.acf <- acf(goa_sst$month.anom, lag.max = 60)
ebs.acf <- acf(ebs_sst$month.anom, lag.max = 60)
ccf(goa_sst$month.anom, ebs_sst$month.anom)

plot_month_acf <- data.frame(lag = 1:60,
                             GOA = goa.acf$acf[2:61],
                             EBS = ebs.acf$acf[2:61]) %>%
  pivot_longer(cols = -lag)

ggplot(plot_month_acf, aes(lag, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_x_continuous(breaks = c(1,12,24,36,48,60))

ggsave("./Figures/goa_ebs_monthly_acf.png", units = 'in', width = 6, height = 4)

# seasonal!

goa.ann <- tapply(goa_sst$month.anom, goa_sst$Year, mean)
goa.123 <- tapply(goa_sst$month.anom[goa_sst$Month %in% 1:3], goa_sst$Year[goa_sst$Month %in% 1:3], mean)
goa.456 <- tapply(goa_sst$month.anom[goa_sst$Month %in% 4:6], goa_sst$Year[goa_sst$Month %in% 4:6], mean)
goa.789 <- tapply(goa_sst$month.anom[goa_sst$Month %in% 7:9], goa_sst$Year[goa_sst$Month %in% 7:9], mean)
goa.101112 <- tapply(goa_sst$month.anom[goa_sst$Month %in% 10:12], goa_sst$Year[goa_sst$Month %in% 10:12], mean)

ebs.ann <- tapply(ebs_sst$month.anom, ebs_sst$Year, mean)
ebs.123 <- tapply(ebs_sst$month.anom[ebs_sst$Month %in% 1:3], ebs_sst$Year[ebs_sst$Month %in% 1:3], mean)
ebs.456 <- tapply(ebs_sst$month.anom[ebs_sst$Month %in% 4:6], ebs_sst$Year[ebs_sst$Month %in% 4:6], mean)
ebs.789 <- tapply(ebs_sst$month.anom[ebs_sst$Month %in% 7:9], ebs_sst$Year[ebs_sst$Month %in% 7:9], mean)
ebs.101112 <- tapply(ebs_sst$month.anom[ebs_sst$Month %in% 10:12], ebs_sst$Year[ebs_sst$Month %in% 10:12], mean)

goa_seas <- data.frame(lag = 0:6,
                       annual = acf(goa.ann, lag.max = 6)$acf, 
                       JFM = acf(goa.123, lag.max = 6)$acf,
                       AMJ = acf(goa.456, lag.max = 6)$acf,
                       JAS = acf(goa.789, lag.max = 6)$acf,
                       OND = acf(na.omit(goa.101112), lag.max = 6)$acf,
                       region = "goa")

ebs_seas <- data.frame(lag = 0:6,
                       annual = acf(ebs.ann, lag.max = 6)$acf, 
                       JFM = acf(ebs.123, lag.max = 6)$acf,
                       AMJ = acf(ebs.456, lag.max = 6)$acf,
                       JAS = acf(ebs.789, lag.max = 6)$acf,
                       OND = acf(na.omit(ebs.101112), lag.max = 6)$acf,
                       region = "ebs")

seasonal <- rbind(goa_seas, ebs_seas) %>%
  pivot_longer(cols = c(-lag, -region)) %>%
  filter(lag > 0)

ggplot(seasonal, aes(lag, value, fill = region)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cb[c(2,6)]) +
  facet_wrap(~name)

ggsave("./Figures/seasonal_autocorrelation_ebs_goa.png", width = 6, height = 6, units = 'in')
