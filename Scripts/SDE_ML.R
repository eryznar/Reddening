# evaluate time-dependent slp-sst relationships affecting Bering Sea ecosystem

library(tidyverse)
library(sde)
library(rstan)
library(ggpubr)
library(pracma)
library(nlme)

# plot settings
theme_set(theme_bw())

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load slp & sst data
slp <- read.csv("./Data/monthlySLPanomalies.csv", row.names = 1) %>%
  group_by(Month) %>%
  mutate(scaled_slp = scale(month.anom)[,1]) %>%
  ungroup()

goa.sst <- read.csv("./Data/goa.monthlySSTanomalies.csv", row.names = 1) %>%
  filter(Year >= 1948) %>%
  mutate(dec.yr = Year + (Month - 0.5)/12) %>%
  group_by(Month) %>%
  mutate(scaled_sst = scale(month.anom)[,1]) %>%
  ungroup() %>%
  mutate(detr_scaled_sst = resid(lm(scaled_sst ~ dec.yr)))

ebs.sst <- read.csv("./Data/ebs.monthlySSTanomalies.csv", row.names = 1) %>%
  filter(Year >= 1948) %>%
  mutate(dec.yr = Year + (Month - 0.5)/12) %>%
  group_by(Month) %>%
  mutate(scaled_sst = scale(month.anom)[,1]) %>%
  ungroup() %>%
  mutate(detr_scaled_sst = resid(lm(scaled_sst ~ dec.yr)))

# plot
plot_slp <- slp %>%
  mutate(dec.yr = Year + (Month - 0.5)/12,
         variable = "SLP") %>%
  rename(anomaly = scaled_slp) %>%
  dplyr::select(dec.yr, anomaly, variable)

plot_goa <- goa.sst %>%
  mutate(variable = "GOA SST") %>%
  rename(anomaly = scaled_sst) %>%
  dplyr::select(dec.yr, anomaly, variable)

plot_ebs <- ebs.sst %>%
  mutate(variable = "EBS SST") %>%
  rename(anomaly = scaled_sst) %>%
  dplyr::select(dec.yr, anomaly, variable)

plot_dat <- rbind(plot_slp, plot_goa, plot_ebs)

ggplot(plot_dat, aes(dec.yr, anomaly)) +
  geom_line() +
  facet_wrap(~variable, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2, col = "red") # looks like SST needs to be detrended - will do that in t




# define functions
# ar_ls calculates the process deviations after
# accounting for forcing variables and autocorrelation,
# (1-gamma)
ar_ls = function(time,forcing,gamma) {
  #S(t+1) = (1-GAMMA*DT)*S(t) + F(t)*DT
  forcing = c(forcing - mean(forcing))
  T=length(forcing)
  sig = 0
  
  for(t in 1:(T-1)) {
    #sig[t+1] = -theta*sig[t] + forcing[t]
    sig[t+1] = (1-gamma)*sig[t] + forcing[t]
  }
  
  # next estimates are linearly de-trended
  #s.sig = sig
  sig = sig - lm(sig ~ time)$fitted.values
  # interpolate output on the original time grid
  s.sig=(sig[-1]+sig[-T])/2 # midpoint
  # final step is normalize
  s.sig=s.sig/sd(s.sig)
  return(s.sig)
}

# vector of decorrelation scales
decor <- 1:12

# object to catch results
cor.out <- p.out <-  data.frame()

# create sst df
sst <- plot_dat %>%
  filter(variable != "SLP") %>%
  rename(date = dec.yr,
         sst = anomaly,
         region = variable)

regions <- unique(sst$region)

for(r in 1:length(regions)){ # loop through regions
# r <- 1  
# set up data   
  
temp.sst <- sst %>%
  filter(region == regions[r])

  dat <- data.frame(date = temp.sst$date,
                    sst = temp.sst[,2],
                    slp.0 = slp$month.anom,
                    slp.1 = c(NA, slp$month.anom[1:919]),
                    slp.2 = c(NA, NA, slp$month.anom[1:918]),
                    slp.3 = c(NA, NA, NA, slp$month.anom[1:917]),
                    slp.4 = c(NA, NA, NA, NA, slp$month.anom[1:916]),
                    slp.5 = c(NA, NA, NA, NA, NA, slp$month.anom[1:915]),
                    slp.6 = c(NA, NA, NA, NA, NA, NA, slp$month.anom[1:914]))
  
  
  # and drop NAs
  dat <- na.omit(dat)

for(l in 3:ncol(dat)){ # loop through lags
# l <- 1
for(i in 1:length(decor)){ # loop through decorrelation scale
  # i <- 1
pred_ts = ar_ls(1:nrow(dat), forcing=dat[,l],
                gamma = 1/decor[i])


pred.sst = data.frame(t = dat$date,
                      sst = dat$sst,
                      integrated.slp = c(0,-as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP


cor.out <- rbind(cor.out, 
                 data.frame(region = regions[r],
                            lag = l - 3,
                            decor = decor[i],
                            cor = cor(pred.sst$sst, pred.sst$integrated.slp)))

# and p-values
mod <- nlme::gls(sst ~ integrated.slp, correlation = corAR1(), data = pred.sst)

p.out <- rbind(p.out, 
                 data.frame(region = regions[r],
                            lag = l - 3,
                            decor = decor[i],
                            p_value = summary(mod)$tTable[2,4]))

}

}

}

cor.out



ggplot(cor.out, aes(decor, cor, color = as.factor(lag))) + 
  geom_line() +
  geom_point() +
  facet_wrap(~region) # very different decorrelation scales!

ggsave("./Figures/sst-slp_lag_decorrelation_by_region.png", width = 9, height = 6, units = 'in')

# plot EBS and GOA for report

# 
# ggplot(filter(cor.out, region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska")), aes(decor, cor, color = as.factor(lag))) + 
#   geom_line() +
#   geom_point() +
#   facet_wrap(~region) + # very different decorrelation scales!
#   labs(x = "Decorrelation scale (months)",
#        y = "Correlation coefficient",
#        color = "Lag (months)") +
#   scale_x_continuous(breaks = 1:12)
#   
# ggsave("./Figures/sst-slp_lag_decorrelation_by_region_EBS_GOA.png", width = 9, height = 6, units = 'in')

decor.use <- cor.out %>%
  group_by(region) %>%
  summarise(decor = decor[which.max(cor)])

# check SST decor scale for GOA and EBS
# report.sst <- sst


# decor.EBS <- acf(report.sst$monthly.anom[report.sst$region == "EBS SST"])

# decor.GOA <- acf(report.sst$monthly.anom[report.sst$region == "GOA SST"])


# now loop through and fit at the best decorrelation scale for each region
predicted.sst <- data.frame()

for(i in 1:nrow(decor.use)){

  # i <- 1
  
  temp.sst <- sst %>%
    filter(region == decor.use$region[i])
  
  dat <- data.frame(date = temp.sst$date,
                    sst = temp.sst[,2],
                    slp.0 = slp$month.anom)
  
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp.0,
                  gamma = 1/decor.use$decor[i])
  
  
  predicted.sst = rbind(predicted.sst,
                        data.frame(region = decor.use$region[i],
                        t = temp.sst$date,
                        sst = temp.sst$sst,
                        integrated.slp = c(0,-as.numeric(pred_ts))))

}
  
predicted.sst <- predicted.sst %>%
  pivot_longer(cols = c(-region, -t))

ggplot(predicted.sst, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP", "SST")) +
  facet_wrap(~region)
  # could add correlations for each!
  
# # plot EBS and GOA for report
# ggplot(filter(predicted.sst, region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska")), aes(t, value, color = name)) +
#   geom_hline(yintercept = 0) +
#   geom_line() +
#   scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP", "SST")) +
#   facet_wrap(~region, scales = "free_y", ncol = 1) +
#   theme(legend.title = element_blank(),
#         axis.title.x = element_blank()) +
#   ylab("Anomaly")
# 
# ggsave("./figs/EBS_GOA_SST_integrated_SLP_time_series.png", width = 7, height = 5)



# get statistics to report
## first ebs
temp.ebs <- predicted.sst %>%
  filter(region == "EBS SST") %>%
  dplyr::select(t, name, value) %>%
  pivot_wider(names_from = name)

cor(temp.ebs$sst, temp.ebs$integrated.slp) # r = 0.18

mod <- nlme::gls(sst ~ integrated.slp, corAR1(), data = temp.ebs)
summary(mod)$tTable[2,4] # p = 0.08

## now goa
temp.goa <- predicted.sst %>%
  filter(region == "GOA SST") %>%
  dplyr::select(t, name, value) %>%
  pivot_wider(names_from = name)

cor(temp.goa$sst, temp.goa$integrated.slp) # r = 0.365

mod <- nlme::gls(sst ~ integrated.slp, corAR1(), data = temp.goa)
summary(mod)$tTable[2,4] # p = 2.31748e-07

# set window length
# calculate AR(1) and SD
out_ar <- data.frame(region = c(rep("GOA", 461*2), rep("EBS", 461*2)),
                     time_series = c(rep("SST", 461), rep("int_SLP", 461), rep("SST", 461), rep("int_SLP", 461)),
                  AR1 = c(sapply(rollapply(temp.goa$sst, width = 460, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
                          sapply(rollapply(temp.goa$integrated.slp, width = 460, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
                          sapply(rollapply(temp.ebs$sst, width = 460, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
                          sapply(rollapply(temp.ebs$integrated.slp, width = 460, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2)))
                  
                  
 out_sd <- data.frame(t = rep(temp.goa$t, 4), 
                      region = c(rep("GOA", 920*2), rep("EBS", 920*2)),
                      time_series = c(rep("SST", 920), rep("int_SLP", 920), rep("SST", 920), rep("int_SLP", 920)),
                      SD = c(rollapply(temp.goa$sst, width = 460, FUN = sd, fill = NA),
                             rollapply(temp.goa$integrated.slp, width = 460, FUN = sd, fill = NA),
                             rollapply(temp.ebs$sst, width = 460, FUN = sd, fill = NA),
                             rollapply(temp.ebs$integrated.slp, width = 460, FUN = sd, fill = NA)))


 
 # Calculate windows
 t <- na.omit(out_sd) %>% pull(t)
 
 # Make data frame of ar1
 cbind(t, out_ar) -> ar.dat
 
 # Join
 left_join(out_sd, ar.dat, relationship = "many-to-many") -> ar.sd.dat
   
   
ggplot(ar.sd.dat, aes(t, AR1, color = time_series)) +
  facet_wrap(~region, scales = "free_y", nrow= 2)+
  geom_line()

ggplot(ar.sd.dat, aes(t, SD, color = time_series)) +
  facet_wrap(~region, scales = "free_y", nrow = 2)+
  geom_line()

ggplot() +
  geom_line(ar.sd.dat %>% filter(time_series == "int_SLP"), mapping = aes(t, SD, color = time_series))+
  geom_line(ar.sd.dat %>% filter(time_series == "SST"), mapping = aes(t, AR1, color = time_series))+
  facet_wrap(~region, scales = "free_y", nrow = 2)+
  ylab("Value")+
  scale_color_manual(values = c("salmon", "steelblue"), labels = c("int_SLP SD", "SST AR1"))


## integrate slp to recreate PDO

# load PDO