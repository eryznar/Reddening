### LOAD PACKAGES --------------------------------------------------------------
source("./Scripts/load.libs.functions.R")
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

### LOAD/PROCESS DATA ----------------------------------------------------------

# NPH sd in winter and spring
nph <- read.csv(paste0(dir, "Output/NPH_sd.csv")) %>%
  dplyr::select(!X)

# SLP winter anomalies in high activity center
slp <- read.csv(paste0(dir, "Output/SLP.winter.anom.csv")) %>%
  #filter(year %in% ebs.sst$Year) %>%
  dplyr::select(!X) %>%
  rename(Year = year) %>%
  filter(Year %in% nph$year)

# Run function to detrend/calculate ar1 and sd on rolling windows
trend.fun(slp, "SLP", 15) -> out

slp2 = as.data.frame(out)


# Scale both
slp.norm <- na.omit(slp2) %>%
  mutate(sd.norm = scale(sd))

nph.norm <- na.omit(nph) %>%
  group_by(season) %>%
  mutate(sd = sd*100, # convert to pascals
         sd.norm = scale(sd))

nph.norm <- nph.norm[-44,]

### PLOT -----------------------------------------------------------------------
# Winter
nph.win <- nph.norm %>% filter(season == "Winter")

mod <- lm(slp.norm$sd.norm ~ nph.win$sd.norm)

p <- data.frame(x = 2010, y = -1.8, plab = "<0.001", rlab = round(sqrt(summary(mod)$r.squared),2))


ggplot()+
  geom_line(nph.win, mapping = aes(year, sd.norm, color = "1"), linewidth = 1.5)+
  geom_point(nph.win, mapping = aes(year, sd.norm, color = "1"), size = 2.25)+
  geom_line(slp.norm, mapping = aes(year, sd.norm, color = "2"), linewidth = 1.5)+
  geom_point(slp.norm, mapping = aes(year, sd.norm, color = "2"), size = 2.25)+
  theme_bw()+
  scale_color_manual(name = "", values = c("darkgoldenrod", "darkviolet"), 
                     labels = c("North Pacific High (pa)\nstandard deviation", 
                                "Aleutian Low (pa)\nstandard deviation"))+
  theme_bw()+
  ylab("Normalized value")+
  xlab("Year")+
  ggtitle("Winter NPH")+
  #geom_richtext(data = p, aes(x = x,  y = y, label = lab))+
  geom_richtext(data = p, aes(x = x,  y = y,
                              label = paste0(plab, "<br>\nr = ", rlab)), size = 4)+
  theme(legend.position = "bottom",
        legend.direction= "horizontal",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16)) -> plot.1

# Spring
nph.sp <- nph.norm %>% filter(season == "Spring")

mod <- lm(slp.norm$sd.norm ~ nph.sp$sd.norm)

p <- data.frame(x = 2010, y = -1.8, plab = "<0.001", rlab = round(sqrt(summary(mod)$r.squared),2))


ggplot()+
  geom_line(nph.sp, mapping = aes(year, sd.norm, color = "1"), linewidth = 1.5)+
  geom_point(nph.sp, mapping = aes(year, sd.norm, color = "1"), size = 2.25)+
  geom_line(slp.norm, mapping = aes(year, sd.norm, color = "2"), linewidth = 1.5)+
  geom_point(slp.norm, mapping = aes(year, sd.norm, color = "2"), size = 2.25)+
  theme_bw()+
  scale_color_manual(name = "", values = c("darkgoldenrod", "darkviolet"), 
                     labels = c("North Pacific High (pa)\nstandard deviation", 
                                "Aleutian Low (pa)\nstandard deviation"))+
  theme_bw()+
  ylab("Normalized value")+
  xlab("Year")+
  ggtitle("Spring NPH")+
  #geom_richtext(data = p, aes(x = x,  y = y, label = lab))+
  geom_richtext(data = p, aes(x = x,  y = y,
                              label = paste0(plab, "<br>\nr = ", rlab)), size = 4)+
  theme(legend.position = "bottom",
        legend.direction= "horizontal",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16)) -> plot.2

plot.1 + plot.2 + plot_layout(nrow =2, guides = "collect") & theme(legend.position = "bottom")

ggsave("./Figures/ALvsNPH_sd.png", width = 5, height = 8, units = "in")
