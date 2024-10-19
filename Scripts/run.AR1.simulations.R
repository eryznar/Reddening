# APPROACH: fit GAMs between SST AR1/SD and time

# Load 
ebs.sst <- read.csv("./Output/SST.anom.ebs.csv") %>%
  filter(Year %in% 1948:2024)
goa.sst <- read.csv("./Output/SST.anom.goa.csv")%>%
  filter(Year %in% 1948:2024)

# Fit linear trend to EBS
ebs.mod <- lm(mean.sst ~ Year, ebs.sst)

ar1.ebs <- acf(ebs.sst$mean.sst, lag.max = 1, plot = FALSE)$acf[2]
ar1.goa <- acf(goa.sst$mean.sst, lag.max = 1, plot = FALSE)$acf[2]

# Extract residuals
resid.ebs <- data.frame(TS = rep("SST", length(unique(ebs.sst$Year))),
                        Year = ebs.sst$Year, resid = ebs.mod$residuals)

resid.ebs %>% filter(Year %in% 2018:2019) %>%
  pull(resid) -> ebs.outliers

# Fit linear trend to GOA
goa.mod <- lm(mean.sst ~ Year, goa.sst)

# Extract residuals
resid.goa <- data.frame(TS = rep("SST", length(unique(ebs.sst$Year))),
                        Year = goa.sst$Year, resid = goa.mod$residuals)

resid.goa %>% filter(Year %in% 2014:2016) %>%
  pull(resid) -> goa.outliers

# Calculate AR1 and SD
trend.fun(ebs.sst, "SST", 15) -> ebs.sst.out

trend.fun(goa.sst, "SST", 15) -> goa.sst.out

ar1var.EBS.sst <- as.data.frame(ebs.sst.out)
ar1var.goa.sst <- as.data.frame(goa.sst.out)

max(na.omit(ar1var.EBS.sst)$ar1)
max(na.omit(ar1var.goa.sst)$ar1)
min(na.omit(ar1var.EBS.sst)$ar1)
min(na.omit(ar1var.goa.sst)$ar1)

max(na.omit(ar1var.EBS.sst)$sd)
max(na.omit(ar1var.goa.sst)$sd)
min(na.omit(ar1var.EBS.sst)$sd)
min(na.omit(ar1var.goa.sst)$sd)


# plot time series with variable levels of first-order autocorrelation
# set ar(1) and SD for BSAI and GOA to observed mins and maxes levels

conditions <- data.frame(system = rep(c("EBS", "GOA"), each = 2),
                         state = rep(c("White noise", "Red noise"), 2),
                         ar = c(-0.7, 0.64, -0.51, 0.47),
                         sd = c(0.24, 0.76, 0.25, 0.58), 
                         trend = rep(c(0.021, 0.011), each = 2))

# set ar(1) and SD for BSAI and GOA to observed mins and maxes levels

conditions <- data.frame(system = rep(c("EBS", "GOA"), each = 2),
                         state = rep(c("White noise", "Red noise"), 2),
                         ar = c(0.0001, 0.64, 0.0001, 0.47), # from 15 year running mean
                         sd = c(0.24, 0.76, 0.25, 0.58), 
                         trend = rep(c(0.021, 0.011), each = 2))
output <- data.frame()


set.seed(999)
for(i in 1:nrow(conditions)){
  
  # Set parameters
  N <- 10000  # Number of time points
  phi <- conditions$ar[i]  # AR(1) coefficient
  sd <- conditions$sd[i] # Standard deviation of error term
  beta <- conditions$trend[i]  # Trend coefficient
  
  # Generate time vector and trend component
  time <- 1:N
  trend <- beta * time
  
  # Simulate AR(1) process
  ar_component <- arima.sim(model = list(ar = phi), n = N, sd = sd)
  
  # Combine trend and AR(1) components
  simulated_series <- trend + ar_component
  
  output <- rbind(output,
                  data.frame(system = conditions$system[i],
                             state = conditions$state[i],
                             time = time,
                             temperature = as.vector(simulated_series),
                             diff_temp = c(NA, diff(simulated_series))))  
  
}


ggplot(filter(output, time <=100), aes(time, temperature)) +
  geom_line() +
  facet_grid(system~state)+
  geom_smooth(method = "lm", se = F)

# plot first differences
ggplot(filter(output, time <=100), aes(time, diff_temp)) +
  geom_line() +
  facet_grid(system~state) 

# calculate change of > 1 degree change in 1 year
sum <- na.omit(output) %>%
  group_by(system, state) %>%
  reframe(prop_1_degree = sum(abs(diff_temp) > 1)/length(diff_temp),
          prop_out_degree = case_when((system == "EBS") ~ sum(abs(diff_temp) > max(ebs.outliers))/length(diff_temp),
                                      TRUE ~ sum(abs(diff_temp) > max(goa.outliers))/length(diff_temp))) %>%
  distinct()

sum

# detrend

output_detrend <- data.frame()

for(i in 1:nrow(sum)){
  
  temp <- output %>%
    filter(system == sum$system[i],
           state == sum$state[i])
  
  mod <- lm(temperature ~ time, temp)
  
  temp_out <- data.frame(system = sum$system[i],
                         state = sum$state[i],
                         time = 1:N,
                         detrended_temp = mod$residuals)
  
  
  #temp_out$detrended_temp_5 <- zoo::rollmean(temp_out$detrended_temp, 15, "right")
  
  output_detrend <- rbind(output_detrend, temp_out)
  
}


# calculate change of > 1 degree change in 1 year
sum_detrend <- na.omit(output_detrend) %>%
  group_by(system, state) %>%
  reframe(prop_1_degree = sum(detrended_temp > 1)/length(detrended_temp),
          prop_out_degree = case_when((system == "EBS") ~ sum(detrended_temp > max(ebs.outliers))/length(detrended_temp),
                                      TRUE ~ sum(detrended_temp > max(goa.outliers))/length(detrended_temp))) %>%
  distinct()

sum_detrend



## RUN simulation
ar = seq(0.01, 0.99, by = 0.1)
sd = seq(0.01, 0.99, length.out = length(ar))
system = c("EBS", "GOA")

conditions <- expand.grid(ar, sd, system) %>%
  rename(ar = Var1, sd = Var2, system = Var3) %>%
  # mutate(state = case_when((ar<0) ~ "white noise",
  #                          TRUE ~ "red noise"),
  mutate(trend = case_when((system == "EBS") ~ 0.021,
                           TRUE ~ 0.011))

output <- data.frame()


for(i in 1:nrow(conditions)){
  
  # Set parameters
  N <- 10000  # Number of time points
  phi <- conditions$ar[i]  # AR(1) coefficient
  sd <- conditions$sd[i] # Standard deviation of error term
  beta <- conditions$trend[i]  # Trend coefficient
  
  # Generate time vector and trend component
  time <- 1:N
  trend <- beta * time
  
  # Simulate AR(1) process
  ar_component <- arima.sim(model = list(ar = phi), n = N, sd = sd)
  
  # Combine trend and AR(1) components
  simulated_series <- trend + ar_component
  
  output <- rbind(output,
                  data.frame(system = conditions$system[i],
                             ar = conditions$ar[i],
                             sd = conditions$sd[i],
                             #state = conditions$state[i],
                             time = time,
                             temperature = as.vector(simulated_series),
                             diff_temp = c(NA, diff(simulated_series))))  
  
}

# calculate change of > 1 degree change in 1 year
summary <- na.omit(output) %>%
  group_by(system, ar, sd) %>%
  reframe(prop_1_degree = sum(abs(diff_temp) > 1)/length(diff_temp),
          prop_out_degree = case_when((system == "EBS") ~ sum(abs(diff_temp) > max(ebs.outliers))/length(diff_temp),
                                      TRUE ~ sum(abs(diff_temp) > max(goa.outliers))/length(diff_temp))) %>%
  distinct()


output_detrend <- data.frame()

for(i in 1:nrow(conditions)){
  
  temp <- output %>%
    filter(system == conditions$system[i],
           ar == conditions$ar[i],
           sd == conditions$sd[i])
  
  
  mod <- lm(temperature ~ time, temp)
  
  temp_out <- data.frame(system = conditions$system[i],
                         ar = conditions$ar[i],
                         sd = conditions$sd[i],
                         time = 1:N,
                         detrended_temp = mod$residuals)
  
  
  #temp_out$detrended_temp_5 <- zoo::rollmean(temp_out$detrended_temp, 15, "right")
  
  output_detrend <- rbind(output_detrend, temp_out)
  
}

summary <- na.omit(output_detrend) %>%
  group_by(system, ar, sd) %>%
  reframe(prop_1_degree = sum(detrended_temp > 1)/length(detrended_temp),
          prop_out_degree = case_when((system == "EBS") ~ sum(detrended_temp > max(ebs.outliers))/length(detrended_temp),
                                      TRUE ~ sum(abs(detrended_temp) > max(goa.outliers))/length(detrended_temp))) %>%
  distinct()

ggplot()+
  geom_tile(summary, mapping = aes(ar, sd, fill = prop_out_degree))+
  scale_fill_viridis_c(option = "turbo", name = "Proportion \nanomalous conditions")+
  facet_wrap(~system)+
  # geom_hline(data.frame(system = c("EBS", "GOA"), y = c(0.76, 0.58)), mapping = aes(yintercept = y),
  #            color = "white", linewidth = 1.5, linetype = "dashed")+
  # geom_vline(data.frame(system = c("EBS", "GOA"), x = c(0.64, 0.47)), mapping = aes(xintercept = x),
  #            color = "white", linewidth = 1.5, linetype = "dashed")+
  geom_segment(data_frame(x = c(-Inf, 0.64, -Inf, 0.47),
                          y = c(0.75, 0.75, 0.58, 0.58),
                          xend = c(0.64, 0.64, 0.47, 0.47),
                          yend = c(0.75, -Inf, 0.58, -Inf),
                          system = c("EBS", "EBS", "GOA", "GOA")),
               mapping = aes(x=x, xend=xend, y=y, yend=yend), linewidth = 1.5,
               color = "white", alpha = 0.5, linetype = "dashed")+
  ylab("Standard deviation")+
  xlab("Autocorrelation")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
