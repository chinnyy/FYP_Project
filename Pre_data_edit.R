## Setup
###Load packages & plot setup


library(plyr) # data wrangling
library(dplyr) # data wrangling
library(tidyverse) # data wrangling and management
library(data.table) # To read large csv
library(tidyverse) # data wrangling and management
library(ggplot2) # visualization


## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

asv_sample_wide<-read.table(file = paste0(wd,'/DATA/ALL/metapr2_ASVs_selected_abundance_Eukaryota.tsv'), sep = '\t', header = TRUE)

asv_sample_wide <- as.data.frame(unclass(asv_sample_wide),                     # Convert all columns to factor
                       stringsAsFactors = TRUE)

woa23_04deg_temp <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_t00mn04.csv'))

# Prepare data for sea surface temperature

woa23_04deg_ss_temp <- cbind(woa23_04deg_temp[,1:2], coalesce(!!!woa23_04deg_temp[,3:23])) # Fill with non-missing vector until 100m 

colnames(woa23_04deg_ss_temp) <- c("latitude","longitude","ss_temperature") # Renaming the columns to match sample data

woa23_04deg_ss_temp <- cbind(latitude= round_any(woa23_04deg_ss_temp$latitude, 1,f = ceiling),
                             longitude=round_any(woa23_04deg_ss_temp$longitude, 1,f = ceiling),
                             woa23_04deg_ss_temp[,3]) # rounding off 

# finding average of repeated locations 
woa23_04deg_ss_temp<- woa23_04deg_ss_temp %>%
  #group_by(latitude,longitude) %>%
  replace(is.na(.), 0) %>%
  aggregate( . ~ latitude+longitude, mean)%>%
as.data.frame()

asv_sample_wide_2dp<- cbind(asv_sample_wide[,1:19],
                            latitude=round_any(asv_sample_wide$latitude,1,f = ceiling),
                            longitude=round_any(asv_sample_wide$longitude,1,f = ceiling),
                            asv_sample_wide[,22:47]) # rounding off 

# Adding sst column to asv_sample_wide dataset
asv_sample_wide_sst<-merge(asv_sample_wide_2dp, woa23_04deg_ss_temp, by = c("latitude", "longitude"), all.x = T)

#THIS CODE IS TO CHECK VALIDITY OF SST DATA
# Merging down to 3181 samples
asv_sample_wide_merge <- asv_sample_wide_sst %>%
  group_by(file_code,ss_temperature) %>%
  replace(is.na(.), 0) %>%
  aggregate(. ~ file_code, mean)%>%
  as.data.frame()

# replacing climate column using annual temp data

asv_sample_wide_sst <- mutate(asv_sample_wide_sst, climate = case_when((ss_temperature > 18 ) ~ "tropical",
                                                                       (ss_temperature < 10 ) ~ "polar",
                                                                       (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) 

asv_sample_wide_merge <- mutate(asv_sample_wide_merge, climate = case_when((ss_temperature > 18 ) ~ "tropical",
                                                                       (ss_temperature < 10 ) ~ "polar",
                                                                       (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) 


# Checking the distribution of climate type

# Create base map
worldmap <- map_data ("world")# From the tidyverse package

base_world  <- ggplot()+ 
  coord_fixed() +
  xlab("") + ylab("") + 
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), 
               colour="dark green", fill="light green")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'light blue', colour = 'light blue'), 
        axis.line = element_line(colour = "white"), #legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

# Save plot 
png(file=paste0(wd,"/PLOTS/Pre/Climate_type_map.jpeg"),width=1536,height=802) 


base_world +
  geom_point(data = asv_sample_wide_merge, aes(x = longitude, y = latitude,fill = climate), size = 4, 
             shape = 21)+
  scale_fill_brewer(palette = "RdYlBu", direction=-1)

dev.off()

# Prepare data for salinity 
woa23_04deg_sali <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_s00mn04.csv'))

# Converting from wide to long data 
woa23_04deg_sali_long <- gather(woa23_04deg_sali, depth, salinity, 3:104, factor_key=TRUE)

# Rounding off 
woa23_04deg_sali_long <- cbind(latitude= round_any(woa23_04deg_sali_long$LATITUDE, 1,f = ceiling),
                             longitude=round_any(woa23_04deg_sali_long$LONGITUDE, 1,f = ceiling),
                             woa23_04deg_sali_long[,3:4]) # rounding off


# finding average of repeated locations 
woa23_04deg_sali_long_round<- woa23_04deg_sali_long %>%
  #group_by(latitude,longitude) %>%
  replace(is.na(.), 0) %>%
  aggregate( . ~ latitude+longitude+depth, mean)%>%
  as.data.frame()

# Adding WOA_salinity column to asv_sample_wide dataset
asv_sample_wide_sali<-merge(asv_sample_wide_2dp, woa23_04deg_sali_long_round, by = c("latitude", "longitude","depth"), all.x = T)

# Replacing extreme values of salinity with NA
asv_sample_wide_sali$salinity.x[asv_sample_wide_sali$salinity.x < 30 | asv_sample_wide_sali$salinity.x > 40] <- NA

# Fill missing salinity data with WOA salinity data 

coalesce_salinity_df <- as.data.frame(coalesce(asv_sample_wide_sali$salinity.x,
                                                              asv_sample_wide_sali$salinity.y))

asv_sample_wide_sali<- cbind(asv_sample_wide_sali,
                             coalesce_salinity = coalesce_salinity_df)

colnames(asv_sample_wide_sali)[49]  <- "coalesce_salinity"    # change column name for x column

## Add coalesce_salinity to asv_sample_wide_sst
asv_sample_wide_sst_sali<- cbind (asv_sample_wide[,1:21], # For orignal lat and long values
                                  asv_sample_wide_sst[,22:48], # For corrected climate types and sst
                                  asv_sample_wide_sali$coalesce_salinity) # For salinity

colnames(asv_sample_wide_sst_sali)[49]  <- "coalesce_salinity"    # change column name for x column

# Replacing extreme values of salinity with NA
asv_sample_wide_sst_sali$salinity[asv_sample_wide_sst_sali$salinity < 30 | asv_sample_wide_sst_sali$salinity > 40] <- NA

asv_sample_wide_sst_sali$coalesce_salinity[asv_sample_wide_sst_sali$coalesce_salinity < 30 | asv_sample_wide_sst_sali$coalesce_salinity > 40] <- NA

## Export asv_sample_wide_sst_sali to csv 
write.csv(asv_sample_wide_sst_sali, paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'), row.names=FALSE)

# TO VALIDATE DATA
# finding average of repeated locations 
asv_sample_wide_sali_merge<- asv_sample_wide_sali %>%
  group_by(latitude,longitude,depth,coalesce_salinity) %>%
  replace(is.na(.), 0) %>%
  aggregate( . ~ latitude+longitude+depth, mean)%>%
  as.data.frame()
