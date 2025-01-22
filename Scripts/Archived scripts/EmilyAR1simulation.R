# goa cod, 4 knots, 3 year sst smooth
goa.r0.sst %>%
  filter(TS == "goa.cod.r0") -> dat

temp.ar1 <- acf(dat$threeyear, lag.max = 1, plot = FALSE)$acf[2]
temp.sd <- sd(dat$threeyear)

cod.ar1 <- acf(dat$log.recruitment, lag.max = 1, plot = FALSE)$acf[2]
cod.sd <- sd(dat$log.recruitment)

N <- nrow(dat)


# Best model

cod.sim <- arima.sim(n = N, model = list(order = c(1,0,0),  ar = cod.ar1, sd = cod.sd)) # length of ebs SST ts

data.frame(dat, cod.sim = cod.sim) -> dat

mod <- gam(cod.sim ~ s(threeyear, k=4), correlation = corAR1(), data = dat)


# Get ar1 for threeyear sst
ar1.vals <- seq(0.1, 0.9, by = 0.03)
sd.vals <- seq(0.1, 2, length.out = length(ar1.vals))

# Simulate model 1000 times with increasing ar1

sim.fun2 <- function(ar1.vals, iter){
  sim.out <- data.frame()  
 for(ss in 1:length(sd.vals)){
   for(aa in 1:length(ar1.vals)){
     ts <- arima.sim(model = list(order = c(1,0,0),  ar = ar1.vals[aa]), sd = sd.vals[ss], n = N) # length of ebs SST ts
     data.frame(threeyear = ts) -> dat2
     
     predict(mod, newdata = dat2) -> pred.vals
     
     # sim.out <- rbind(sim.out, data.frame(ar1 = ar1.vals[aa],
     #                                      temp.sd = sd.vals[ss],
     #                                      iteration = iter,
     #                                      temp.mean = mean(ts),
     #                                      #temp.sd = sd(ts),
     #                                      rec.mean = mean(pred.vals),
     #                                      rec.sd = sd(pred.vals),
     #                                      rec.cv = sd(pred.vals)/(mean(pred.vals))))
     
     sim.out <- rbind(sim.out, data.frame(ar1 = rep(ar1.vals[aa], length(ts)),
                                          temp.sd = rep(sd.vals[ss], length(ts)),
                                          iteration = rep(iter, length(ts)),
                                          ts = c(ts),
                                          pred.vals = pred.vals))
   }
 }
  
  return(sim.out)
}

1:100 %>%
purrr::map_df(~sim.fun2(ar1.vals, .x)) -> out

# out %>%
#   group_by(ar1, temp.sd, iteration) %>%
#  mutate(rec.sd = sd(pred.vals),
#         rec.mean = mean(pred.vals),
#         N = n(),
#         var = case_when((pred.vals >= pred.vals+rec.sd | pred.vals <= pred.vals-rec.sd) ~ "TRUE",
#                       TRUE ~ "FALSE")) %>%
#   ungroup() %>%
#   group_by(ar1, temp.sd) %>%
#   reframe(prop.var = sum(var == TRUE)/N) -> out2

out %>%
  group_by(ar1, temp.sd) %>%
  reframe(rec.sd = sd(pred.vals)) -> out2


ggplot()+
  geom_tile(out2, mapping = aes(ar1, temp.sd, fill = rec.sd))+
  scale_fill_viridis_c(option = "turbo", name = "Recruitment variation")+
  ylab("SST variation")+
  xlab("AR1")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

ggsave("./Figures/reddening.impacts.r0.simmodel.png", width = 11, height = 8.5, units = "in")
