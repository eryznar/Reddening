source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")


nc.early <- stack("./Data/ERA5_sst_1960-1982.nc")
nc.mid <- brick("./Data/ERA5_sst_1983-2005.nc")
#nc.late.pre24 <- rast(raster::stack("./Data/ERA5_sst_2006-2024.nc"))
nc.late <- brick("./Data/ERA5_sst_2006-2023.nc")

lme <- st_read("./Data/LME shapefiles/LMEs66.shp")
north40 <- c(1, 2, 8, 9, 18, 19, 20, 21, 22, 23, 24, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64)
subset(lme, lme$ARCTIC=="Arctic") -> arcticlme
subset(lme, lme$LME_NUMBER %in% north40) ->  north40lme

unique(north40lme$LME_NAME)
drop <- c("Laptev Sea", "Kara Sea", "Canadian Eastern Arctic - West Greenland", 
          "Beaufort Sea", "Central Arctic", "East Siberian Sea")

subset(north40lme, !north40lme$LME_NAME %in% drop) -> north40lme


unique(north40lme$LME_NAME) -> lme_names
unique(north40lme$LME_NUMBER) -> lme_num

sst_df <- data.frame()
yrs <- 1960:2023

for(ii in 1:length(lme_names)){
    subset(north40lme, north40lme$LME_NAME == lme_names[ii]) -> shp
    
    terra::extract(nc.early, shp) %>%
      as.data.frame() -> early
    terra::extract(nc.mid, shp) %>%
      as.data.frame() -> mid
    terra::extract(nc.late, shp) %>%
      as.data.frame() -> late
    
    cbind(early, mid, late) %>%
      na.omit() -> tt
    
    tt[!duplicated(colnames(tt))] -> dat
    
    dat %>%
      mutate(cell.mean = rowMeans(.)) %>%
      pivot_longer(., cols = 1:768, names_to = "time", values_to = "sst") %>%
      mutate(anom = sst-cell.mean,
             year = str_sub(time, 2, 5),
             month = str_sub(time, 7, 8),
             LME = lme_names[ii])  %>%
      group_by(year, month, LME) %>%
      reframe(mean.sst = mean(sst-273.15),
              mean.anom = mean(anom)) -> df
    
    rbind(sst_df, df) -> sst_df
    
}


sst_df %>%
  right_join(., data.frame("LME" = lme_names, "LME_NUM"= lme_num), by = "LME") -> sst_df

sst_df %>%
  mutate(month = as.numeric(as.character(month)),
         year = as.numeric(as.character(year)),
         year = case_when((month %in% c(11, 12)) ~ (year +1),
                          TRUE ~ year)) %>%
  filter(month %in% c(11,12, 1, 2, 3)) %>%
  group_by(year, LME, LME_NUM) %>%
  reframe(mean.sst = mean(mean.sst),
          mean.anom = mean(mean.anom)) -> sst_df2
  


# sst_df <- read.csv("./Output/lme.sst_df.csv") %>%
#   mutate(mean.sst = mean.sst-273.15)
write.csv(sst_df2, "./Output/lme.sst_df.csv")
