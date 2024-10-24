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
  group_by(Month) %>%
  mutate(scaled_sst = scale(month.anom)[,1]) %>%
  ungroup() %>%
  mutate(scaled_sst = detrend(scaled_sst))

# Recreating AR1 that I sent
    # Detrend data
    detrend.dat <- loess(month.anom ~ Year, goa.sst, span = 0.25, degree = 1)
    
    # Extract residuals
    TS.dat <- data.frame(Year = goa.sst$Year, month.anom = detrend.dat$residuals)
    pp <- sapply(rollapply(TS.dat$month.anom, width = 460, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2)
    plot(pp, type = "l")


ebs.sst <- read.csv("./Data/ebs.monthlySSTanomalies.csv", row.names = 1) %>%
  filter(Year >= 1948) %>%
  group_by(Month) %>%
  mutate(scaled_sst = scale(month.anom)[,1]) %>%
  ungroup() %>%
  mutate(scaled_sst = detrend(scaled_sst))

# plot
plot_slp <- slp %>%
  mutate(dec.yr = Year + (Month - 0.5)/12,
         variable = "SLP") %>%
  rename(anomaly = scaled_slp) %>%
  dplyr::select(dec.yr, anomaly, variable)

plot_goa <- goa.sst %>%
  mutate(dec.yr = Year + (Month - 0.5)/12,
         variable = "GOA SST") %>%
  rename(anomaly = scaled_sst) %>%
  dplyr::select(dec.yr, anomaly, variable)

plot_ebs <- ebs.sst %>%
  mutate(dec.yr = Year + (Month - 0.5)/12,
         variable = "EBS SST") %>%
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

# reset plotting order for regions
reg.ord <- data.frame(region = regions,
                      order = 1:6)

cor.out <- left_join(cor.out, reg.ord)

cor.out$region <- reorder(cor.out$region, cor.out$order)

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

cor(temp.ebs$sst, temp.ebs$integrated.slp) # r = 0.219

mod <- nlme::gls(sst ~ integrated.slp, corAR1(), data = temp.ebs)
summary(mod)$tTable[2,4] # p = 0.0175

## now goa
temp.goa <- predicted.sst %>%
  filter(region == "GOA SST") %>%
  dplyr::select(t, name, value) %>%
  pivot_wider(names_from = name)

cor(temp.goa$sst, temp.goa$integrated.slp) # r = 0.365

mod <- nlme::gls(sst ~ integrated.slp, corAR1(), data = temp.goa)
summary(mod)$tTable[2,4] # p = 2.31748e-07


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



# Calculate rolling window SD
sd <-  rollapply(temp.goa$sst, width = 460, FUN = sd)


# something is wrong!!


## below is old
                  
  goa.ar.sst <- sapply(rollapply(temp.goa, width = 460, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
  





# now loop through on 240-month (20 year) rolling windows

# first, fit EOF to each rolling window of SLP to use for all regions
slp.pc1.time.series <- data.frame()

for(i in 240:nrow(slp)){ # loop through windows
# i <- 240

# parcel out the 240 months ending at time == i
temp.slp <- slp[(i-239):(i),] 

# and fit EOF 
pca.temp <- svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]

# scale
pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!

# save
slp.pc1.time.series <- rbind(slp.pc1.time.series, 
                             data.frame(i = i,
                                        pc1.temp = pc1.temp))

}

slp.pc1.time.series <- slp.pc1.time.series %>%
  rename(window = i)

# make object to capture correlation output
cor.out <- p.out <- data.frame()

for(r in 2:(length(regions)-1)){ # excluding N. Pac. and SCC!
# r <- 2
  # limit to regions of interest
  temp.sst <- sst %>%
    filter(region == regions[r])
  
for(i in 240:nrow(temp.sst)){ # loop through windows
  # i <- 240
  
  # now select sst

  temp.sst.window <- temp.sst[(i-239):i,]
  
  temp.slp.window <- slp.pc1.time.series %>%
    filter(window == i)
  
  # fit AR(1) model
  pred_ts = ar_ls(1:nrow(temp.slp.window), forcing=temp.slp.window$pc1.temp,
                  gamma = 1/decor.use$decor[(r-1)])
  

  cor.out <- rbind(cor.out,
                   data.frame(regions = regions[r],
                              end.date = as.Date(temp.sst.window[nrow(temp.sst.window),3]),
                              r = cor(temp.sst.window[,2], c(0, as.numeric(pred_ts)))))
  
  # and p-values
  temp.dat <- data.frame(sst = temp.sst.window[,2],
                         integrated.slp = c(0, as.numeric(pred_ts))
  )
  
  mod <- nlme::gls(sst ~ integrated.slp, correlation = corAR1(), data = temp.dat)
  
  p.out <- rbind(p.out, 
                 data.frame(region = regions[r],
                            end.date = as.Date(temp.sst.window[nrow(temp.sst.window),3]),
                            p_value = summary(mod)$tTable[2,4]))
  

} # close window loop

} # close region loop

cor.out <- cor.out %>%
  pivot_longer(cols = c(-end.date, -regions)) %>%
  rename(region = regions) %>%
  left_join(., reg.ord) %>%
  mutate(region = reorder(region, order))

ggplot(cor.out, aes(end.date, value)) +
  geom_line() +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = as.Date("1989-01-01"), lty=3) +
  theme(axis.title.x = element_blank()) +
  ylab("correlation")
  

ggsave("./figs/20_year_SLP_SDE_vs_SST_by region.png", width = 6, height = 9, units = 'in')

# plot EBS and GOA for report
ebs.goa.cor.plot <- ggplot(filter(cor.out, region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska")), aes(end.date, value)) +
  geom_line() +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = as.Date("1989-01-01"), lty=3) +
  theme(axis.title.x = element_blank()) +
  ylab("Correlation coefficient")

ebs.goa.cor.plot

ggsave("./figs/20_year_SLP_SDE_vs_SST_by_region_EBS_GOA.png", width = 6, height = 5, units = 'in')

# plot p-values
p.out <- p.out %>%
  pivot_longer(cols = c(-end.date, -region)) %>%
  left_join(., reg.ord) %>%
  mutate(region = reorder(region, order))

ebs.goa.p.plot <- ggplot(filter(p.out, region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska")), 
                         aes(end.date, log(value, 10))) +
  geom_line() +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = as.Date("1989-01-01"), lty=3) +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0),
                     labels = c("0.001", "0.01", "0.1", "1")) +
  geom_hline(yintercept = -2, lty = 2) +
  ylab("p-value")

ebs.goa.p.plot

# combined plot
png("./figs/combined_cor_p_plots.png", width = 8, height = 4, units = 'in', res = 300)

ggpubr::ggarrange(ebs.goa.cor.plot,
                  ebs.goa.p.plot,
                  ncol = 2, 
                  labels = "auto")

dev.off()

# and we really just want to plot r values
r.plot <- cor.out %>%
  filter(name == "r")

rolling.window <- ggplot(r.plot, aes(end.date, value)) +
  geom_line() +
  labs(x = "End year", 
       y = "r")

rolling.window

ggsave("./figs/20_year_SLP_SDE_vs_SST_correlation.png", width = 6, height = 4, units = 'in')


# combine SDE panels

png("./figs/Bering_SDE_plots.png", 
    width=5.5, height=9, units='in', res=300)

ggpubr::ggarrange(sst.reconstruct, sst.reconstruct.sm, rolling.window, nrow = 3, ncol = 1,
                  labels = "auto")

dev.off()




## fit Bayesian regression to correlation time series -------------------------
# (can we distinguish the trend from 0?)
library(plyr)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./scripts/stan_utils.R")

# subset data
dat <- cor.out %>%
  dplyr::filter(name == "r") %>%
  dplyr::select(-name) %>%
  dplyr::rename(r = value) %>%
  mutate(year = lubridate::year(end.date),
         month = lubridate::month(end.date)) 

# now identify a winter year corresponding to Oct - Sept
dat$winter.year <- if_else(dat$month %in% 10:12, dat$year + 1, dat$year)

dat <- dat %>%
  dplyr::group_by(winter.year) %>%
  dplyr::summarise(r = mean(r))

# check to plot
# dat$dec.yr <- dat$year + (dat$month-0.5)/12
# 
# ggplot(dat, aes(dec.yr, r)) +
#   geom_line() +
#   geom_point() 

ggplot(dat, aes(winter.year, r)) +
  geom_line() +
  geom_point()
  
# # define model formula
# dat$month_fac <- as.factor(dat$month)

# slp_sst_cor_formula <-  bf(r ~ s(year, k = 5) + (1 | month_fac/year) + ar(p = 1)) # limiting k to ~ number of decades in time series
# 
# # show default priors
# get_prior(slp_sst_cor_formula, dat)
# 
# # set priors
# priors <- c(set_prior("normal(0, 3)", class = "ar"),
#             set_prior("normal(0, 3)", class = "b"),
#             set_prior("normal(0, 3)", class = "Intercept"),
#             set_prior("student_t(3, 0, 3)", class = "sd"),
#             set_prior("student_t(3, 0, 3)", class = "sds"),
#             set_prior("student_t(3, 0, 3)", class = "sigma"))

slp_sst_cor_formula <-  bf(r ~ s(winter.year, k = 4) + ar(p = 1)) # limiting k to ~ number of decades in time series

# show default priors
get_prior(slp_sst_cor_formula, dat)

# set priors
priors <- c(set_prior("normal(0, 3)", class = "ar"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"))


# fit model 
sst_slp_cor_brm <- brm(slp_sst_cor_formula,
                    prior = priors,
                    data = dat,
                    cores = 4, chains = 4, iter = 5000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99999, max_treedepth = 10))

codR_dfa_brm  <- add_criterion(codR_dfa_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_dfa_brm, file = "output/codR_dfa_brm.rds")

codR_dfa_brm <- readRDS("./output/codR_dfa_brm.rds")
check_hmc_diagnostics(sst_slp_cor_brm$fit)
neff_lowest(sst_slp_cor_brm$fit)
rhat_highest(sst_slp_cor_brm$fit)
summary(sst_slp_cor_brm)
bayes_R2(sst_slp_cor_brm)
plot(sst_slp_cor_brm$criteria$loo, "k")
plot(conditional_effects(sst_slp_cor_brm), ask = FALSE)

 ## fit model without ar()
slp_sst_cor_formula <-  bf(r ~ s(winter.year, k = 4)) # limiting k to ~ number of decades in time series

# show default priors
get_prior(slp_sst_cor_formula, dat)

# set priors
priors <- c(set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"))


# fit model 
sst_slp_cor_brm <- brm(slp_sst_cor_formula,
                       prior = priors,
                       data = dat,
                       cores = 4, chains = 4, iter = 3000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.9999, max_treedepth = 12))


check_hmc_diagnostics(sst_slp_cor_brm$fit)
neff_lowest(sst_slp_cor_brm$fit)
rhat_highest(sst_slp_cor_brm$fit)
summary(sst_slp_cor_brm)
bayes_R2(sst_slp_cor_brm)
plot(sst_slp_cor_brm$criteria$loo, "k")
plot(conditional_effects(sst_slp_cor_brm), ask = FALSE)

# so these models only show a trend if ar(1) term isn't used!
# either that means a trend can't be demonstrated given the autocorrelation in the data
# or I'm asking too much of the model to effectively use the same 
# term as a covariate and an AR term
# might be a bad idea to fit any statistical models to these highly overlapping windows!



## plot sst-slp regression maps for contrasting eras -----------------------------
# 1950-1968, 1969-1988, 1989-2008

library(maps)
library(maptools)
library(mapdata)
library(fields)

# load slp and cell weights
slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)

# get a dates vector to deal with!
d <- chron::dates(rownames(slp))
yr <- as.numeric(as.character(chron::years(d)))

# fix bad years
change <- yr > 2030
yr[change] <- yr[change]-100

# get y and z
x <- seq(152.5, 230, 2.5)
y <- seq(40, 62.5, 2.5)

# load mean EBS sst anomaly
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)
sst$year <- as.numeric(as.character(chron::years(sst$date)))

# fix bad year values
fix <- sst$year > 2030
sst$year[fix] <- sst$year[fix] - 100

# select the three sst windows (response)
sst1 <- sst %>%
  filter(year %in% 1950:1968)

sst2 <- sst %>%
  filter(year %in% 1969:1988)

sst3 <- sst %>%
  filter(year %in% 1989:2008)

# select the corresponding slp windows (1 year behind to allow smoothing!)
slp1 <- slp[yr %in% 1950:1968,]

slp2 <- slp[yr %in% 1969:1988,]

slp3 <- slp[yr %in% 1988:2008,]

# smooth and scale slp
ff <- function(x) as.vector(scale(rollmean(x, 6, fill = NA, align = "right")))
slp.sm.1 <- apply(slp1, 2, ff)
d1 <- chron::dates(rownames(slp))[yr %in% 1950:1968]
sst1$date

# remove first sst value to accommodate lag
sst1 <- sst1[2:nrow(sst1),]

# and last row of slp1
slp.sm.1 <- slp.sm.1[1:(nrow(slp.sm.1)-1),]

# now loop through slp columns, fit lm, and keep coef
coef1 <- NA

for(i in 1:ncol(slp.sm.1)){
  # i <- 1
  mod <- lm(sst1$anom ~ slp.sm.1[,i])
  coef1[i] <- mod$coefficients[2]
}

## second window
# smooth and scale slp
slp.sm.2 <- apply(slp2, 2, ff)
d2 <- chron::dates(rownames(slp))[yr %in% 1969:1988]
sst2$date

# # trim slp.sm.2
# slp.sm.2 <- slp.sm.2[12:(nrow(slp.sm.2)-1),]

# remove first sst value to accommodate lag
sst2 <- sst2[2:nrow(sst2),]

# and last row of slp1
slp.sm.2 <- slp.sm.2[1:(nrow(slp.sm.2)-1),]


# now loop through slp columns, fit lm, and keep coef
coef2 <- NA

for(i in 1:ncol(slp.sm.2)){
  # i <- 1
  mod <- lm(sst2$anom ~ slp.sm.2[,i])
  coef2[i] <- mod$coefficients[2]
}

## third window
# smooth and scale slp
slp.sm.3 <- apply(slp3, 2, ff)
d3 <- chron::dates(rownames(slp))[yr %in% 1988:2008]
sst3$date
d3
# trim slp.sm.2
slp.sm.3 <- slp.sm.3[12:(nrow(slp.sm.3)-1),]

# now loop through slp columns, fit lm, and keep coef
coef3 <- NA

for(i in 1:ncol(slp.sm.3)){
  # i <- 1
  mod <- lm(sst3$anom ~ slp.sm.3[,i])
  coef3[i] <- mod$coefficients[2]
}

#plot

# get range
range <- range(coef1, coef2, coef3)
zlim <- c(range[1], -range[1])

# draw sebs box for sst area
box.x <- c(187, 199, 199, 203, 203, 199, 199, 193, 193, 187)
box.y <- c(61, 61, 59, 59, 57, 57, 55, 55, 53, 53)

# get x and y for plotting
x <- read.csv("./output/slp_x.csv")
x <- x[,1]

y <- read.csv("./output/slp_y.csv")
y <- y[,1]


png("./figs/EBS era sst-slp regressions.png", 3, 6, units="in", res=300)

par(mfcol=c(3,1), mar=c(0,0.5,2,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

z <- coef1 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("a) 1950-1968",  cex=1, side=3, adj=0)

z <- coef2 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("b) 1969-1988",  cex=1, side=3, adj=0)

z <- coef3 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("c) 1989-2008",  cex=1, side=3, adj=0)

dev.off()

## now same maps for GOA----------
# load slp and cell weights
slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)

# get a dates vector to deal with!
d <- chron::dates(rownames(slp))
yr <- as.numeric(as.character(chron::years(d)))

# fix bad years
change <- yr > 2030
yr[change] <- yr[change]-100

# load mean GOA sst anomaly
sst <- read.csv("./data/regional_monthly_sst.csv") %>%
  filter(region == "Gulf_of_Alaska") %>%
  dplyr::select(-region) %>%
  rename(anom = monthly.anom)

sst$year <- as.numeric(as.character(chron::years(sst$date)))

# fix bad year values
fix <- sst$year > 2030
sst$year[fix] <- sst$year[fix] - 100

# select the three sst windows (response)
sst1 <- sst %>%
  filter(year %in% 1950:1968)

sst2 <- sst %>%
  filter(year %in% 1969:1988)

sst3 <- sst %>%
  filter(year %in% 1989:2008)

# select the corresponding slp windows (1 year behind to allow smoothing!)
slp1 <- slp[yr %in% 1950:1968,]

slp2 <- slp[yr %in% 1969:1988,]

slp3 <- slp[yr %in% 1988:2008,]

# smooth and scale slp
ff <- function(x) as.vector(scale(rollmean(x, 6, fill = NA, align = "right")))
slp.sm.1 <- apply(slp1, 2, ff)
d1 <- chron::dates(rownames(slp))[yr %in% 1950:1968]
sst1$date

# remove first sst value to accommodate lag
sst1 <- sst1[2:nrow(sst1),]

# and last row of slp1
slp.sm.1 <- slp.sm.1[1:(nrow(slp.sm.1)-1),]

# now loop through slp columns, fit lm, and keep coef
coef1_GOA <- NA

for(i in 1:ncol(slp.sm.1)){
  # i <- 1
  mod <- lm(sst1$anom ~ slp.sm.1[,i])
  coef1_GOA[i] <- mod$coefficients[2]
}

## second window
# smooth and scale slp
slp.sm.2 <- apply(slp2, 2, ff)
d2 <- chron::dates(rownames(slp))[yr %in% 1969:1988]
sst2$date

# # trim slp.sm.2
# slp.sm.2 <- slp.sm.2[12:(nrow(slp.sm.2)-1),]

# remove first sst value to accommodate lag
sst2 <- sst2[2:nrow(sst2),]

# and last row of slp1
slp.sm.2 <- slp.sm.2[1:(nrow(slp.sm.2)-1),]


# now loop through slp columns, fit lm, and keep coef
coef2_GOA <- NA

for(i in 1:ncol(slp.sm.2)){
  # i <- 1
  mod <- lm(sst2$anom ~ slp.sm.2[,i])
  coef2_GOA[i] <- mod$coefficients[2]
}

## third window
# smooth and scale slp
slp.sm.3 <- apply(slp3, 2, ff)
d3 <- chron::dates(rownames(slp))[yr %in% 1988:2008]
sst3$date
d3
# trim slp.sm.2
slp.sm.3 <- slp.sm.3[12:(nrow(slp.sm.3)-1),]

# now loop through slp columns, fit lm, and keep coef
coef3_GOA <- NA

for(i in 1:ncol(slp.sm.3)){
  # i <- 1
  mod <- lm(sst3$anom ~ slp.sm.3[,i])
  coef3_GOA[i] <- mod$coefficients[2]
}

#plot

# get range
range <- range(coef1_GOA, coef2_GOA, coef3_GOA)
zlim <- c(range[1], -range[1])

# draw goa box for sst area
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

# get x and y for plotting
x <- read.csv("./output/slp_x.csv")
x <- x[,1]

y <- read.csv("./output/slp_y.csv")
y <- y[,1]


png("./figs/GOA era sst-slp regressions.png", 3, 6, units="in", res=300)

par(mfcol=c(3,1), mar=c(0,0.5,2,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

z <- coef1_GOA # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(goa.x, goa.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("a) 1950-1968",  cex=1, side=3, adj=0)

z <- coef2_GOA # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(goa.x, goa.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("b) 1969-1988",  cex=1, side=3, adj=0)

z <- coef3_GOA # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(goa.x, goa.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("c) 1989-2008",  cex=1, side=3, adj=0)

dev.off()

# and combine both regions in one plot
range <- range(coef1_GOA, coef2_GOA, coef3_GOA, coef1, coef2, coef3)
zlim <- c(range[1], -range[1])

png("./figs/GOA EBS era sst-slp regressions.png", 6.2, 6, units="in", res=300)

par(mfcol=c(3,2), mar=c(0,0.5,2,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

z <- coef1 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("a) EBS 1950-1968",  cex=1, side=3, adj=0)

z <- coef2 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("b) EBS 1969-1988",  cex=1, side=3, adj=0)

z <- coef3 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("c) EBS 1989-2008",  cex=1, side=3, adj=0)

z <- coef1_GOA # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(goa.x, goa.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("d) GOA 1950-1968",  cex=1, side=3, adj=0)

z <- coef2_GOA # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(goa.x, goa.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("e) GOA 1969-1988",  cex=1, side=3, adj=0)

z <- coef3_GOA # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(goa.x, goa.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("f) GOA 1989-2008",  cex=1, side=3, adj=0)

dev.off()
