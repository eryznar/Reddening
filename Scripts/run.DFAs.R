# EBS RAW RECRUITMENT ----
value <- "val"
years = 1987:2016
data <- ar1var.EBS.r0
labs <- r0.labs.bsai
region <- "Eastern Bering Sea"

# process data
data %>%
  rename(resp = grep(value, colnames(.))) %>%
  filter(year %in% years) %>%
  mutate(scaled = scale(resp)) %>%
  dplyr::select(TS, year, scaled) %>%
  pivot_wider(., names_from = TS, values_from = scaled) %>%
  dplyr::select(!year) %>%
  as.data.frame() %>%
  t() -> t.dat

# Specify colnames and change to matrix
colnames(t.dat) <- years

dfa.dat <- as.matrix(t.dat)

# find best error structure for 1-trend model

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices

levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS::MARSS(dfa.dat, model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare

model.data$dAICc <- model.data$AICc-min(model.data$AICc)

model.data <- model.data %>%
  arrange(dAICc)

model.data

# diagonal and equal and diagonal and uequal best (unconstrained did not fit)
cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="diagonal and equal")
recruit.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

recruit.CI <- MARSSparamCIs(recruit.mod)

recruit.plot.CI <- data.frame(names=rownames(dfa.dat),
                              mean=recruit.CI$par$Z,
                              upCI=recruit.CI$par.upCI$Z,
                              lowCI=recruit.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

recruit.plot.CI$names <- reorder(recruit.plot.CI$names, recruit.CI$par$Z)

# Plot
cc <- ifelse(region == "Eastern Bering Sea", "#6A6DB7", "#A34242")
value <- ifelse(value == "val", "raw", value)

recruit.loadings.plot <- ggplot(recruit.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cc) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  ggtitle(region)+
  scale_x_discrete(labels = labs)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 0)

# plot trend
recruit.trend <- data.frame(t=years,
                            estimate=as.vector(recruit.mod$states),
                            conf.low=as.vector(recruit.mod$states)-1.96*as.vector(recruit.mod$states.se),
                            conf.high=as.vector(recruit.mod$states)+1.96*as.vector(recruit.mod$states.se))

recruit.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, recruit.trend = estimate) -> r0.trend

write.csv(r0.trend, paste0("./Output/r0.trend.", value, ".", region, ".csv"))

recruit.trend.plot <- ggplot(recruit.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cc) +
  geom_hline(yintercept = 0) +
  #ggtitle(paste0(region, " ", value, " DFA trend"))+
  scale_x_continuous(breaks = seq(min(recruit.trend$t), max(recruit.trend$t), by = 2))+
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cc) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))+
  ylab("Trend")+
  xlab("")


ggpubr::ggarrange(recruit.loadings.plot, recruit.trend.plot, nrow = 2, #labels = "auto"
                  widths = c(0.4, 0.6))

ggsave(paste0("./Figures/", region, ".", value, ".DFAplot.png"), width = 8.5, height = 11, units = "in")

# GOA RAW RECRUITMENT ----
value <- "val"
years = 1987:2016
data <- ar1var.goa.r0
labs <- r0.labs.goa
region <- "Gulf of Alaska"

# process data
data %>%
  rename(resp = grep(value, colnames(.))) %>%
  filter(year %in% years) %>%
  mutate(scaled = scale(resp)) %>%
  dplyr::select(TS, year, scaled) %>%
  pivot_wider(., names_from = TS, values_from = scaled) %>%
  dplyr::select(!year) %>%
  as.data.frame() %>%
  t() -> t.dat

# Specify colnames and change to matrix
colnames(t.dat) <- years

dfa.dat <- as.matrix(t.dat)

# find best error structure for 1-trend model

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices

levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS::MARSS(dfa.dat, model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare

model.data$dAICc <- model.data$AICc-min(model.data$AICc)

model.data <- model.data %>%
  arrange(dAICc)

model.data

# diagonal and equal and diagonal and uequal best (unconstrained did not fit)
cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="equalvarcov")
recruit.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

recruit.CI <- MARSSparamCIs(recruit.mod)

recruit.plot.CI <- data.frame(names=rownames(dfa.dat),
                              mean=recruit.CI$par$Z,
                              upCI=recruit.CI$par.upCI$Z,
                              lowCI=recruit.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

recruit.plot.CI$names <- reorder(recruit.plot.CI$names, recruit.CI$par$Z)

# Plot
cc <- ifelse(region == "Eastern Bering Sea", "#6A6DB7", "#A34242")
value <- ifelse(value == "val", "raw", value)

recruit.loadings.plot <- ggplot(recruit.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cc) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  ggtitle(region)+
  scale_x_discrete(labels = labs)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 0)

# plot trend
recruit.trend <- data.frame(t=years,
                            estimate=as.vector(recruit.mod$states),
                            conf.low=as.vector(recruit.mod$states)-1.96*as.vector(recruit.mod$states.se),
                            conf.high=as.vector(recruit.mod$states)+1.96*as.vector(recruit.mod$states.se))

recruit.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, recruit.trend = estimate) -> r0.trend

write.csv(r0.trend, paste0("./Output/r0.trend.", value, ".", region, ".csv"))

recruit.trend.plot <- ggplot(recruit.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cc) +
  geom_hline(yintercept = 0) +
  #ggtitle(paste0(region, " ", value, " DFA trend"))+
  scale_x_continuous(breaks = seq(min(recruit.trend$t), max(recruit.trend$t), by = 2))+
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cc) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))+
  ylab("Trend")+
  xlab("")


ggpubr::ggarrange(recruit.loadings.plot, recruit.trend.plot, nrow = 2, #labels = "auto"
                  widths = c(0.4, 0.6))

ggsave(paste0("./Figures/", region, ".", value, ".DFAplot.png"), width = 8.5, height = 11, units = "in")

# EBS R0 AR1 ----
value <- "ar1"
years = 1987:2016
data <- ar1var.EBS.r0
labs <- r0.labs.bsai
region <- "Eastern Bering Sea"

# process data
data %>%
  rename(resp = grep(value, colnames(.))) %>%
  filter(year %in% years) %>%
  mutate(scaled = scale(resp)) %>%
  dplyr::select(TS, year, scaled) %>%
  pivot_wider(., names_from = TS, values_from = scaled) %>%
  dplyr::select(!year) %>%
  as.data.frame() %>%
  t() -> t.dat

# Specify colnames and change to matrix
colnames(t.dat) <- years

dfa.dat <- as.matrix(t.dat)

# find best error structure for 1-trend model

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices

levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS::MARSS(dfa.dat, model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare

model.data$dAICc <- model.data$AICc-min(model.data$AICc)

model.data <- model.data %>%
  arrange(dAICc)

model.data

# diagonal and unequal and equalvarcov (unconstrained did not fit)
cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="diagonal and unequal")
recruit.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

recruit.CI <- MARSSparamCIs(recruit.mod)

recruit.plot.CI <- data.frame(names=rownames(dfa.dat),
                              mean=recruit.CI$par$Z,
                              upCI=recruit.CI$par.upCI$Z,
                              lowCI=recruit.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

recruit.plot.CI$names <- reorder(recruit.plot.CI$names, recruit.CI$par$Z)

# Plot
cc <- ifelse(region == "Eastern Bering Sea", "#6A6DB7", "#A34242")

recruit.loadings.plot <- ggplot(recruit.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cc) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  ggtitle(region)+
  scale_x_discrete(labels = labs)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 0)

# plot trend
recruit.trend <- data.frame(t=years,
                            estimate=as.vector(recruit.mod$states),
                            conf.low=as.vector(recruit.mod$states)-1.96*as.vector(recruit.mod$states.se),
                            conf.high=as.vector(recruit.mod$states)+1.96*as.vector(recruit.mod$states.se))

recruit.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, recruit.trend = estimate) -> r0.trend

write.csv(r0.trend, paste0("./Output/r0.trend.", value, ".", region, ".csv"))

recruit.trend.plot <- ggplot(recruit.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cc) +
  geom_hline(yintercept = 0) +
  #ggtitle(paste0(region, " ", value, " DFA trend"))+
  scale_x_continuous(breaks = seq(min(recruit.trend$t), max(recruit.trend$t), by = 2))+
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cc) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))+
  ylab("Trend")+
  xlab("")


ggpubr::ggarrange(recruit.loadings.plot, recruit.trend.plot, nrow = 2, #labels = "auto"
                  widths = c(0.4, 0.6))

ggsave(paste0("./Figures/", region, ".", value, ".DFAplot.png"), width = 8.5, height = 11, units = "in")
       

# GOA R0 AR1 ----
value <- "ar1"
years = 1987:2016
data <- ar1var.goa.r0
labs <- r0.labs.goa
region <- "Gulf of Alaska"

# process data
data %>%
  rename(resp = grep(value, colnames(.))) %>%
  filter(year %in% years) %>%
  mutate(scaled = scale(resp)) %>%
  dplyr::select(TS, year, scaled) %>%
  pivot_wider(., names_from = TS, values_from = scaled) %>%
  dplyr::select(!year) %>%
  as.data.frame() %>%
  t() -> t.dat

# Specify colnames and change to matrix
colnames(t.dat) <- years

dfa.dat <- as.matrix(t.dat)

# find best error structure for 1-trend model

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices

levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS::MARSS(dfa.dat, model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare

model.data$dAICc <- model.data$AICc-min(model.data$AICc)

model.data <- model.data %>%
  arrange(dAICc)

model.data

# diagonal and unequal and equalvarcov (unconstrained did not fit)
cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="diagonal and unequal")
recruit.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

recruit.CI <- MARSSparamCIs(recruit.mod)

recruit.plot.CI <- data.frame(names=rownames(dfa.dat),
                              mean=recruit.CI$par$Z,
                              upCI=recruit.CI$par.upCI$Z,
                              lowCI=recruit.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

recruit.plot.CI$names <- reorder(recruit.plot.CI$names, recruit.CI$par$Z)

# Plot
cc <- ifelse(region == "Eastern Bering Sea", "#6A6DB7", "#A34242")

recruit.loadings.plot <- ggplot(recruit.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cc) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  ggtitle(region)+
  scale_x_discrete(labels = labs)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 0)

# plot trend
recruit.trend <- data.frame(t=years,
                            estimate=as.vector(recruit.mod$states),
                            conf.low=as.vector(recruit.mod$states)-1.96*as.vector(recruit.mod$states.se),
                            conf.high=as.vector(recruit.mod$states)+1.96*as.vector(recruit.mod$states.se))

recruit.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, recruit.trend = estimate) -> r0.trend

write.csv(r0.trend, paste0("./Output/r0.trend.", value, ".", region, ".csv"))

recruit.trend.plot <- ggplot(recruit.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cc) +
  geom_hline(yintercept = 0) +
  #ggtitle(paste0(region, " ", value, " DFA trend"))+
  scale_x_continuous(breaks = seq(min(recruit.trend$t), max(recruit.trend$t), by = 2))+
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cc) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))+
  ylab("Trend")+
  xlab("")


ggpubr::ggarrange(recruit.loadings.plot, recruit.trend.plot, nrow = 2, #labels = "auto"
                  widths = c(0.4, 0.6))

ggsave(paste0("./Figures/", region, ".", value, ".DFAplot.png"), width = 8.5, height = 11, units = "in")




# DFA AR1 trends vs. sst AR1 ----
# EBS
dfa.ar1.ebs <- read.csv("./Output/r0.trend.ar1.Eastern Bering Sea.csv")

ar1var.EBS.sst %>%
  filter(year %in% dfa.ar1.ebs$year) -> sst.ar1.ebs

gam(dfa.ar1.ebs$recruit.trend ~ s(sst.ar1.ebs$ar1), correlation = corAR1()) -> mod.1

pred.r0 = predict(mod.1, se.fit =TRUE)$fit
pred.CI = 1.96*(predict(mod.1, se.fit =TRUE)$se.fit)

data.frame(dfa.ar1 = dfa.ar1.ebs$recruit.trend, sst.ar1 = sst.ar1.ebs$ar1, pred.ar1 = pred.r0, pred.C1 = pred.CI) -> plot.dat.ebs

plot.1 <- ggplot()+
            geom_ribbon(plot.dat.ebs, mapping = aes(x = sst.ar1, ymin = pred.ar1 - pred.CI, ymax= pred.ar1+pred.CI),
                        fill = "grey", alpha = 0.75)+
            geom_point(plot.dat.ebs, mapping=aes(x = sst.ar1, y = dfa.ar1))+
            geom_line(plot.dat.ebs, mapping = aes(x = sst.ar1, y = pred.ar1), size = 1.25, color = "#6A6DB7")+
            ylab("Recruitment AR1 trend (DFA)")+
            xlab("SST AR1")+
            ggtitle("Eastern Bering Sea (p<0.001, R2 = 0.69)")+
            theme_bw()+
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12),
                  strip.text = element_text(size = 10),
                  legend.text = element_text(size = 12)) 

ggsave(plot= plot.1, "./Figures/EBS.r0trend.vs.sst.png", width = 8.5, height = 5.5, units = "in")

# GOA
dfa.ar1.goa <- read.csv("./Output/r0.trend.ar1.Gulf of Alaska.csv")

ar1var.goa.sst %>%
  filter(year %in% dfa.ar1.goa$year) -> sst.ar1.goa

gam(dfa.ar1.goa$recruit.trend ~ s(sst.ar1.goa$ar1), correlation = corAR1()) -> mod.1

pred.r0 = predict(mod.1, se.fit =TRUE)$fit
pred.CI = 1.96*(predict(mod.1, se.fit =TRUE)$se.fit)

data.frame(dfa.ar1 = dfa.ar1.goa$recruit.trend, sst.ar1 = sst.ar1.goa$ar1, pred.ar1 = pred.r0, pred.C1 = pred.CI) -> plot.dat.goa

plot.2 <- ggplot()+
            geom_ribbon(plot.dat.goa, mapping = aes(x = sst.ar1, ymin = pred.ar1 - pred.CI, ymax= pred.ar1+pred.CI),
                        fill = "grey", alpha = 0.75)+
            geom_point(plot.dat.goa, mapping=aes(x = sst.ar1, y = dfa.ar1))+
            geom_line(plot.dat.goa, mapping = aes(x = sst.ar1, y = pred.ar1), size = 1.25, color = "#A34242")+
            ylab("Recruitment AR1 trend (DFA)")+
            xlab("SST AR1")+
            ggtitle("Gulf of Alaska (p=0.84, R2 = 0.034)")+
            theme_bw()+
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 12),
                  strip.text = element_text(size = 10),
                  legend.text = element_text(size = 12))

ggsave(plot= plot.2, "./Figures/GOA.r0trend.vs.sst.png", width = 8.5, height = 5.5, units = "in")

# DFA raw r0 trends vs. sst ----
# EBS 
read.csv("./Output/r0.trend.raw.Eastern Bering Sea.csv") %>%
  dplyr::select(!X) -> dfa.trend

dfa.trend %>%
  right_join(., sst.rollmeans %>% filter(region == "Eastern Bering Sea") %>%
               rename(year = Year), by = c("year")) %>%
  na.omit() %>%
  dplyr::select(!c(region, X)) %>%
  filter(year > 1987) -> dfa.sst

TS.dat <- dfa.sst

knts <- c(3, 4, 5)

model.out <- data.frame()

for(kk in 1:length(knts)){
  
  # fit gams, record pval and AIC
  gam.1 <- gam(recruit.trend ~ s(unsmoothed, k = knts[kk]), 
               data= TS.dat, correlation = corAR1())
  gam.2 <- gam(recruit.trend ~ s(twoyear, k = knts[kk]), 
               data= TS.dat, correlation = corAR1())
  gam.3 <- gam(recruit.trend ~ s(threeyear, k = knts[kk]), 
               data= TS.dat, correlation = corAR1())
  
  p.gam.1 <- signif(summary(gam.1)$s.table[,4],2)
  p.gam.2 <- signif(summary(gam.2)$s.table[,4],2)
  p.gam.3 <- signif(summary(gam.3)$s.table[,4],2)
  
  # p.lme.1 <- signif(summary(gam.1$lme)$tTable[2,5],2)
  # p.lme.2 <- signif(summary(gam.2$lme)$tTable[2,5],2)
  # p.lme.3 <- signif(summary(gam.3$lme)$tTable[2,5],2)
  
  rsq.gam.1 <- signif(summary(gam.1)$r.sq,2)
  rsq.gam.2 <- signif(summary(gam.2)$r.sq,2)
  rsq.gam.3 <- signif(summary(gam.3)$r.sq,2)
  
  AIC.gam.1 <- AIC(gam.1)
  AIC.gam.2 <- AIC(gam.2)
  AIC.gam.3 <- AIC(gam.3)
  
  
  # Build summary table
  model.out <- rbind(model.out, data.frame(knots = knts[kk],
                                           #Model = c("gam.1", "gam.2", "gam.3"),
                                           sst = c("unsmoothed", "twoyear", "threeyear"),
                                           #p_lme = c(p.lme.1, p.lme.2, p.lme.3),
                                           p_gam = c(p.gam.1, p.gam.2, p.gam.3),
                                           AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3),
                                           rsq = c(rsq.gam.1, rsq.gam.2, rsq.gam.3)))
  
  } # close knot loop



# Label best models by timeseries
model.out %>%
  #group_by(TS) %>%
  mutate(sig = BH2(p_gam, alph = 0.05)$BHSig,
         p_gam = p_gam,
         padj = p.adjust(p_gam, method = "fdr")) %>%
  mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) %>%
  filter(BEST == "Y") -> model.out


# Fit best model
mod.1 <- gam(recruit.trend ~ s(threeyear, k = 5), 
             data= TS.dat, correlation = corAR1())


pred.r0 = predict(mod.1, se.fit =TRUE)$fit
pred.CI = 1.96*(predict(mod.1, se.fit =TRUE)$se.fit)

data.frame(dfa.r0 = dfa.sst$recruit.trend, sst = dfa.sst$threeyear, pred = pred.r0, pred.C1 = pred.CI) -> plot.dat.ebs

ggplot()+
  geom_ribbon(plot.dat.ebs, mapping = aes(x = sst, ymin = pred - pred.CI, ymax= pred+pred.CI),
              fill = "grey", alpha = 0.75)+
  geom_point(plot.dat.ebs, mapping=aes(x = sst, y = dfa.r0))+
  geom_line(plot.dat.ebs, mapping = aes(x = sst, y = pred), size = 1.25, color = "#6A6DB7")+
  ylab("Recruitment trend (DFA)")+
  xlab("SST (three year smooth)")+
  ggtitle("Eastern Bering Sea (p<0.05*, R2 = 0.17)")+
  theme_bw()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) 

ggsave("./Figures/EBS.r0trend.vs.sst.png", width = 8.5, height = 5.5, units = "in")


# GOA
read.csv("./Output/r0.trend.raw.Gulf of Alaska.csv") %>%
  dplyr::select(!X) -> dfa.trend

dfa.trend %>%
  right_join(., sst.rollmeans %>% filter(region == "Gulf of Alaska") %>%
               rename(year = Year), by = c("year")) %>%
  na.omit() %>%
  dplyr::select(!c(region, X)) %>%
  filter(year > 1987) -> dfa.sst

TS.dat <- dfa.sst

knts <- c(3, 4, 5)

model.out <- data.frame()

for(kk in 1:length(knts)){
  
  # fit gams, record pval and AIC
  gam.1 <- gam(recruit.trend ~ s(unsmoothed, k = knts[kk]), 
               data= TS.dat, correlation = corAR1())
  gam.2 <- gam(recruit.trend ~ s(twoyear, k = knts[kk]), 
               data= TS.dat, correlation = corAR1())
  gam.3 <- gam(recruit.trend ~ s(threeyear, k = knts[kk]), 
               data= TS.dat, correlation = corAR1())
  
  p.gam.1 <- signif(summary(gam.1)$s.table[,4],2)
  p.gam.2 <- signif(summary(gam.2)$s.table[,4],2)
  p.gam.3 <- signif(summary(gam.3)$s.table[,4],2)
  
  # p.lme.1 <- signif(summary(gam.1$lme)$tTable[2,5],2)
  # p.lme.2 <- signif(summary(gam.2$lme)$tTable[2,5],2)
  # p.lme.3 <- signif(summary(gam.3$lme)$tTable[2,5],2)
  
  rsq.gam.1 <- signif(summary(gam.1)$r.sq,2)
  rsq.gam.2 <- signif(summary(gam.2)$r.sq,2)
  rsq.gam.3 <- signif(summary(gam.3)$r.sq,2)
  
  AIC.gam.1 <- AIC(gam.1)
  AIC.gam.2 <- AIC(gam.2)
  AIC.gam.3 <- AIC(gam.3)
  
  
  # Build summary table
  model.out <- rbind(model.out, data.frame(knots = knts[kk],
                                           #Model = c("gam.1", "gam.2", "gam.3"),
                                           sst = c("unsmoothed", "twoyear", "threeyear"),
                                           #p_lme = c(p.lme.1, p.lme.2, p.lme.3),
                                           p_gam = c(p.gam.1, p.gam.2, p.gam.3),
                                           AIC = c(AIC.gam.1, AIC.gam.2, AIC.gam.3),
                                           rsq = c(rsq.gam.1, rsq.gam.2, rsq.gam.3)))
  
} # close knot loop



# Label best models by timeseries
model.out %>%
  #group_by(TS) %>%
  mutate(sig = BH2(p_gam, alph = 0.05)$BHSig,
         p_gam = p_gam,
         padj = p.adjust(p_gam, method = "fdr")) %>%
  mutate(BEST = ifelse(AIC == min(AIC), "Y", "N")) %>%
  filter(BEST == "Y") -> model.out


# Fit best model
mod.1 <- gam(recruit.trend ~ s(twoyear, k = 5), 
             data= TS.dat, correlation = corAR1())


pred.r0 = predict(mod.1, se.fit =TRUE)$fit
pred.CI = 1.96*(predict(mod.1, se.fit =TRUE)$se.fit)

data.frame(dfa.r0 = dfa.sst$recruit.trend, sst = dfa.sst$twoyear, pred = pred.r0, pred.C1 = pred.CI) -> plot.dat.goa

ggplot()+
  geom_ribbon(plot.dat.goa, mapping = aes(x = sst, ymin = pred - pred.CI, ymax= pred+pred.CI),
              fill = "grey", alpha = 0.75)+
  geom_point(plot.dat.goa, mapping=aes(x = sst, y = dfa.r0))+
  geom_line(plot.dat.goa, mapping = aes(x = sst, y = pred), size = 1.25, color = "#A34242")+
  ylab("Recruitment trend (DFA)")+
  xlab("SST (two year smooth)")+
  ggtitle("Gulf of Alaska (p<0.001*, R2 = 0.32)")+
  theme_bw()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) 

ggsave("./Figures/GOA.r0trend.vs.sst.png", width = 8.5, height = 5.5, units = "in")


