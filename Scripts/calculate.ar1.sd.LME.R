# LOAD SOURCE FILE --------------------------------------------------------------
source("./Scripts/load.libs.functions.R")
source("./Scripts/ts.processing.R")

# LOAD DATA ---------------------------------------------------------------------
sst.lme <- read.csv("./Output/lme.sst_df.csv") %>%
  filter(!LME %in% c("East Bering Sea", "Gulf of Alaska"),
         year > 1960 & year < 2024) %>%
  dplyr::select(year, LME, mean.anom)


ar1var.EBS.sst %>%
  mutate(LME = "Eastern Bering Sea") %>%
  dplyr::select(year, val, LME) %>%
  rename(mean.anom = val) -> ebs

ar1var.goa.sst %>%
  mutate(LME = "Gulf of Alaska") %>%
  dplyr::select(year, val, LME) %>%
  rename(mean.anom = val) -> goa

 rbind(sst.lme, ebs, goa)  -> sst.lme

# PLOT DATA ---------------------------------------------------------------------
ggplot(sst.lme, aes(year, mean.anom))+
  geom_line(size = 1)+
  facet_wrap(~LME, scales = "free_y",  labeller = label_wrap_gen())+
  theme_bw()+
  ylab("SSTa (°C)")+
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
  
  # Detrend data
  detrend.dat <- loess(mean.anom ~ year, TS.dat, span = 0.25, degree = 1)
  
  # Extract residuals
  TS.dat <- data.frame(year = TS.dat$year, mean.anom = detrend.dat$residuals)
  
  
  # Calculate rolling window AR1
  ar1 <- sapply(rollapply(TS.dat$mean.anom, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
  
  # Calculate rolling window SD
  sd <-  rollapply(TS.dat$mean.anom, width = width, FUN = sd, fill = NA)
  
  # Make data frame of sd-cv
  data.frame(year = unique(TS.dat$year), val = TS.dat$mean.anom, sd = sd) -> win.dat
  
  # Calculate windows
  win.yr <- na.omit(win.dat) %>% pull(year)
  
  # Make data frame of ar1
  data.frame(year = win.yr, ar1 = ar1) -> ar1.dat
  
  # Join
  left_join(win.dat, ar1.dat) %>%
    mutate(LME = lme[ii]) -> win.dat
  
  # Calculate tau on ar1 and sd
  tau.dat <- na.omit(win.dat)
  tau.ar1 <- cor.test(tau.dat$ar1, tau.dat$year, method = "kendall")
  tau.sd <- cor.test(tau.dat$sd, tau.dat$year, method = "kendall")
  
  tt <- data.frame(tau.ar1 = round(tau.ar1$estimate, 2), tau.ar1.p = round(tau.ar1$p.value, 2),
                   tau.sd = round(tau.sd$estimate, 2),  tau.sd.p = round(tau.sd$p.value, 2))
  
  cbind(win.dat, tt) -> win.dat
  
  # Compile output
  sum.out <- rbind(sum.out, win.dat)
}

# Adjust tau for FDR
sum.out %>%
  dplyr::select(!c(year, val, sd, ar1)) %>%
  distinct() %>%
  mutate(tau.ar1.sig = BH2(tau.ar1.p, alph = 0.05)$BHSig,
         tau.ar1.padj = p.adjust(tau.ar1.p, method = "fdr"),
         tau.sd.sig = BH2(tau.sd.p, alph = 0.05)$BHSig,
         tau.sd.padj = p.adjust(tau.sd.p, method = "fdr")) %>%
  right_join(sum.out, .) -> sum.out2

# PLOT AR1 AND SD ---------------------------------------------------------------
sum.out2 %>%
  dplyr::select(LME, tau.ar1, tau.ar1.p, tau.ar1.sig, tau.sd, tau.sd.p, tau.sd.sig)%>%
  arrange(., LME) %>%
  distinct() %>%
  mutate(p.ar1 = case_when((tau.ar1.p < 0.001) ~ "p<0.001",
                          (tau.ar1.p < 0.01 & tau.ar1.p >= 0.001) ~ "p<0.01",
                          (tau.ar1.p <0.05 & tau.ar1.p >= 0.01) ~ "p<0.05",
                          TRUE ~ paste0("p=", tau.ar1.p)),
         p.sd = case_when((tau.sd.p < 0.001) ~ "p<0.001",
                          (tau.sd.p < 0.01 & tau.sd.p >= 0.001) ~ "p<0.01",
                          (tau.sd.p <0.05 & tau.sd.p >= 0.01) ~ "p<0.05",
                          TRUE ~ paste0("p=", tau.sd.p))) %>%
  mutate(p.ar1 = case_when((tau.ar1.sig == TRUE) ~ paste0(p.ar1, "*"),
                           TRUE ~ paste(p.ar1)),
         p.sd =  case_when((tau.sd.sig == TRUE) ~ paste0(p.sd, "*"),
                           TRUE ~ paste(p.sd))) -> lab.dat

labs <- paste0(lab.dat$LME, " \n(tau=", lab.dat$tau.ar1, ", ", lab.dat$p.ar1, ")")
names(labs) <- lab.dat$LME

ggplot(sum.out2, aes(year, ar1))+
  geom_line(size = 1)+
  facet_wrap(~LME, scales = "free_y",  labeller = labeller(LME = labs))+
  theme_bw()+
  ylab("AR1")+
  xlab("Year")

ggsave("./Figures/lme.sstAR1.png", height = 8.5, width = 11, units = "in")

library(forcats)
sum.out2 %>% mutate(LME = fct_reorder(LME, tau.ar1)) %>%
  dplyr::select(LME, tau.ar1, tau.ar1.sig) %>% distinct()-> data
ggplot(data, aes(LME, tau.ar1))+
  geom_bar(stat = "identity", color = "black", fill = "goldenrod2")+
  geom_text(data = data[data$tau.ar1.sig ==TRUE, ], aes(x = LME, y = tau.ar1, label = "*",
                vjust = ifelse(tau.ar1 >= 0, tau.ar1-0.5 ,tau.ar1 +2)), size = 10) +
  ylim(c(-0.45, 0.55))+
  theme_bw()+
  ylab("Kendall's tau")+
  ggtitle("AR1")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=16), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16),
        title = element_text(size = 18))

ggsave("./Figures/lme.ar1.tau.png", height = 8.5, width = 10.5, units = "in")



labs <- paste0(lab.dat$LME, " \n(tau=", lab.dat$tau.sd, ", ", lab.dat$p.sd, ")")
names(labs) <- lab.dat$LME

ggplot(sum.out, aes(year, sd))+
  geom_line(size = 1)+
  facet_wrap(~LME, scales = "free_y",  labeller = labeller(LME = labs))+
  theme_bw()+
  ylab("SD")+
  xlab("Year")

ggsave("./Figures/lme.sstSD.png", height = 8.5, width = 11, units = "in")

sum.out2 %>% mutate(LME = fct_reorder(LME, tau.ar1)) %>%
  dplyr::select(LME, tau.sd, tau.sd.sig) %>% distinct()-> data
ggplot(data, aes(LME, tau.sd))+
  geom_bar(stat = "identity", color = "black", fill = "goldenrod2")+
  geom_text(data = data[data$tau.sd.sig ==TRUE, ], aes(x = LME, y = tau.sd, label = "*",
                                                        vjust = ifelse(tau.sd >= 0, tau.sd-0.5 ,tau.sd +2)), size = 10) +
  ylim(c(-0.53, 0.8))+
  theme_bw()+
  ylab("Kendall's tau")+
  ggtitle("SD")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=16), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 16),
        title = element_text(size = 18))

ggsave("./Figures/lme.sd.tau.png", height = 8.5, width = 11, units = "in")


sum.out2%>% mutate(LME = fct_reorder(LME, tau.ar1)) %>%
  dplyr::select(!c(year, val, sd, ar1)) %>% distinct()-> data

data %>% 
  dplyr::select(!c(tau.ar1.padj, tau.sd.padj, tau.ar1.p, tau.sd.p)) %>%
pivot_longer(., c(tau.ar1, tau.sd), values_to = "value", names_to = "name") %>%
  mutate(sig = case_when((name == "tau.ar1" & tau.ar1.sig == TRUE) ~ TRUE,
                         (name == "tau.sd" & tau.sd.sig == TRUE) ~ TRUE,
                         TRUE ~ FALSE)) %>%
  dplyr::select(!c(tau.ar1.sig, tau.sd.sig)) -> data2

labs <- c("AR1", "Standard deviation")
names(labs) <- c("tau.ar1", "tau.sd")

ggplot()+
  geom_bar(data2, mapping = aes(LME, value), stat = "identity", color = "black", fill = "goldenrod2", position= "dodge")+
  facet_wrap(~name, nrow = 1, labeller = labeller(name = labs))+
  geom_text(data = data2[data2$sig ==TRUE, ], aes(x = LME, y = value, label = "*", group = name,
                                                       hjust = ifelse(value >= 0, value-0.8 ,value+2)),
            vjust = 0.75, size = 10)+
  theme_bw()+
  ylab("Kendall's tau")+
  xlab("")+
  coord_flip()+
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=16), legend.title = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        title = element_text(size = 18))


ggsave("./Figures/lme.tau.combined.png", height = 9, width = 15, units = "in")

  

# DFA on AR1 --------------------------------------------------------------------
# process data
sum.out %>%
  filter(year %in% 1967:2016) %>%
  rename(resp = ar1) %>%
  mutate(scaled = scale(resp)) %>%
  dplyr::select(LME, year, scaled) %>%
  pivot_wider(., names_from = LME, values_from = scaled) %>%
  dplyr::select(!year) %>%
  as.data.frame() %>%
  t() -> t.dat

# Specify colnames and change to matrix
colnames(t.dat) <- 1967:2016

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

# equalvarcov and diagonal and unequal (unconstrained did not fit)
cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="equalvarcov")
lme.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

lme.CI <- MARSSparamCIs(lme.mod)

lme.plot.CI <- data.frame(names=rownames(dfa.dat),
                              mean=lme.CI$par$Z,
                              upCI=lme.CI$par.upCI$Z,
                              lowCI=lme.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

lme.plot.CI$names <- reorder(lme.plot.CI$names, lme.CI$par$Z)

# Plot
lme.loadings.plot <- ggplot(lme.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill="goldenrod2") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  ggtitle("AR1")+
  #scale_x_discrete(labels = labs)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 0)

# plot trend
lme.trend <- data.frame(t=1967:2016,
                            estimate=as.vector(lme.mod$states),
                            conf.low=as.vector(lme.mod$states)-1.96*as.vector(lme.mod$states.se),
                            conf.high=as.vector(lme.mod$states)+1.96*as.vector(lme.mod$states.se))

lme.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, lme.trend = estimate) -> lme.trend2

write.csv(lme.trend2, "./Output/lme.trend.AR1.csv")

lme.trend.plot <- ggplot(lme.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color="goldenrod2") +
  geom_hline(yintercept = 0) +
  #ggtitle("AR1")+
  scale_x_continuous(breaks = seq(min(lme.trend$t), max(lme.trend$t), by = 5))+
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill="goldenrod2") + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))+
  ylab("Trend")+
  xlab("")


ggpubr::ggarrange(lme.loadings.plot, lme.trend.plot, nrow = 2, #labels = "auto"
                  widths = c(0.4, 0.6))

ggsave(paste0("./Figures/lme.AR1.DFA.png"), width = 8.5, height = 11, units = "in")


# DFA on SD --------------------------------------------------------------------
# process data
sum.out %>%
  filter(year %in% 1967:2016) %>%
  rename(resp = sd) %>%
  mutate(scaled = scale(resp)) %>%
  dplyr::select(LME, year, scaled) %>%
  pivot_wider(., names_from = LME, values_from = scaled) %>%
  dplyr::select(!year) %>%
  as.data.frame() %>%
  t() -> t.dat

# Specify colnames and change to matrix
colnames(t.dat) <- 1967:2016

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
lme.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

lme.CI <- MARSSparamCIs(lme.mod)

lme.plot.CI <- data.frame(names=rownames(dfa.dat),
                          mean=lme.CI$par$Z,
                          upCI=lme.CI$par.upCI$Z,
                          lowCI=lme.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

lme.plot.CI$names <- reorder(lme.plot.CI$names, lme.CI$par$Z)

# Plot

lme.loadings.plot <- ggplot(lme.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill="goldenrod2") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  ggtitle("SD")+
  #scale_x_discrete(labels = labs)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top',
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  geom_hline(yintercept = 0)

# plot trend
lme.trend <- data.frame(t=1967:2016,
                        estimate=as.vector(lme.mod$states),
                        conf.low=as.vector(lme.mod$states)-1.96*as.vector(lme.mod$states.se),
                        conf.high=as.vector(lme.mod$states)+1.96*as.vector(lme.mod$states.se))

lme.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, lme.trend = estimate) -> lme.trend2

write.csv(lme.trend2, "./Output/lme.trend.SD.csv")

lme.trend.plot <- ggplot(lme.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color="goldenrod2") +
  geom_hline(yintercept = 0) +
  #ggtitle("SD")+
  scale_x_continuous(breaks = seq(min(lme.trend$t), max(lme.trend$t), by = 5))+
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill="goldenrod2") + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))+
  ylab("Trend")+
  xlab("")


ggpubr::ggarrange(lme.loadings.plot, lme.trend.plot, nrow = 2, #labels = "auto"
                  widths = c(0.4, 0.6))

ggsave(paste0("./Figures/lme.SD.DFA.png"), width = 8.5, height = 11, units = "in")

