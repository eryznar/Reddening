# plot time series with variable levels of first-order autocorrelation

library(tidyverse)
theme_set(theme_bw())
# set ar(1) and SD for BSAI and GOA to observed mins and maxes levels

conditions <- data.frame(system = rep(c("EBS", "GOA"), each = 2),
                         state = rep(c("White noise", "Red noise"), 2),
                         ar = c(-0.5, 0.73, -0.16, 0.67),
                         sd = c(0.23, 0.77, 0.34, 0.7), 
                         trend = rep(c(0.034, 0.018), each = 2))

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
  summarise(prop_1_degree = sum(abs(diff_temp) > 1)/length(diff_temp),
            prop_2_degree = sum(abs(diff_temp) > 2)/length(diff_temp),
            prop_3_degree = sum(abs(diff_temp) > 3)/length(diff_temp))

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
  
  
  temp_out$detrended_temp_5 <- zoo::rollmean(temp_out$detrended_temp, 5, "right")
  
  output_detrend <- rbind(output_detrend, temp_out)
  
}


# calculate change of > 1 degree change in 1 year
sum_detrend <- na.omit(output_detrend) %>%
  group_by(system, state) %>%
  summarise(prop_1_degree = sum(detrended_temp_5 > 1)/length(detrended_temp_5))

sum_detrend