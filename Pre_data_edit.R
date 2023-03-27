## Setup
###Load packages & plot setup


library(plyr) # data wrangling
library(dplyr) # data wrangling
library(tidyverse) # data wrangling and management
library(data.table) # To read large csv
library(tidyverse) # data wrangling and management
library(ggplot2) # visualization
library(sf) 
library(viridis)
library(raster)


## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

asv_sample_wide<-read.table(file = paste0(wd,'/DATA/ALL/metapr2_ASVs_selected_abundance_Eukaryota.tsv'), sep = '\t', header = TRUE)

asv_sample_wide <- as.data.frame(unclass(asv_sample_wide),  # Convert all columns to factor
                       stringsAsFactors = TRUE)

woa23_temp <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_t00mn01.csv'))

# Prepare data for sea surface temperature

woa23_ss_temp <- cbind(woa23_temp[,1:2], coalesce(!!!woa23_temp[,3:23])) # Fill with non-missing vector until 100m 

colnames(woa23_ss_temp) <- c("latitude","longitude","ss_temperature") # Renaming the columns to match sample data

round_05<- function(value){
  value_1<- value + 0.5 
  value_2<- round(value_1,0)
  value_3<- value_2 -0.5
  
  return(value_3)
}

# Merging down to 3128 samples
asv_sample_wide_merge_location <- asv_sample_wide %>%
  dplyr::select(file_code,longitude,latitude) %>%
  dplyr::group_by(file_code,longitude,latitude) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  as.data.frame() 

# To find out the 0.5 rounded off lat and long for each sample
asv_sample_wide_merge_round<- asv_sample_wide_merge_location
asv_sample_wide_merge_round[,2:3] <- lapply(asv_sample_wide_merge_location[2:3], round_05)

# Adding sst column to asv_sample_wide_merge_location dataset
asv_sample_wide_sst<-merge(asv_sample_wide_merge_round, woa23_ss_temp, by = c("latitude", "longitude"), all.x = T)

# Filtering for NA SST
asv_sample_wide_sst_na <- asv_sample_wide_sst[rowSums(is.na(asv_sample_wide_sst)) > 0,]

# Plotting location of NA values
ggplot() + geom_raster(data = woa23_ss_temp, aes(x=longitude, y = latitude, fill=ss_temperature)) +
  coord_fixed(ratio = 1) +
  geom_jitter(data=asv_sample_wide_sst_na,aes(x=longitude, y = latitude),size=1,shape=4,color="red")+
  labs(fill="Sea surface \ntemperature (°c)") +
  scale_fill_viridis() +
  theme_bw()

ggsave(filename=paste0(wd,"/PLOTS/Pre/sst_Na_location.png"),width=40,height=19,units=c("cm"))


# To solve the issue, find which grid closest to the point has a filled in value and extract that value
woa23_ss_temp <- woa23_ss_temp[,c("longitude", "latitude", "ss_temperature")]

min_dist<- cbind(asv_sample_wide_merge_location[,2:3],index=data.frame(apply(raster::pointDistance(as.matrix(asv_sample_wide_merge_location[,2:3]), 
                            as.matrix(woa23_ss_temp[,1:2]), 
                            lonlat = T), 1,which.min)))

names(min_dist)<- c("longitude",'latitude',"index")

woa23_ss_temp$index <- 1:nrow(woa23_ss_temp)

min_dist<- merge(min_dist, woa23_ss_temp, by = 'index', all.x = TRUE)

# Combining sst to file code and aggregating repeats
asv_sample_woa_sst<- asv_sample_wide_merge_location %>%
  merge(., min_dist[c(2,3,6)], by.x = c('longitude','latitude'),
        by.y = c("longitude.x","latitude.x"))%>%
  dplyr::select(file_code,longitude,latitude,ss_temperature) %>%
  dplyr::group_by(file_code,longitude,latitude,ss_temperature) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  as.data.frame() 

# Adding sst to original dataset
asv_sample_wide<-merge(asv_sample_wide, asv_sample_woa_sst[,c("file_code","ss_temperature")], by = c("file_code"), all.x = T)


# replacing climate column using annual temp data
asv_sample_wide <- mutate(asv_sample_wide, climate = case_when((ss_temperature > 18 ) ~ "tropical",
                                                               (ss_temperature < 10 ) ~ "polar",
                                                               (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) 

asv_sample_merge_sst<- asv_sample_wide %>%
  dplyr::select(file_code,latitude,longitude,climate,ss_temperature) %>%
  dplyr::group_by(file_code,latitude,longitude,climate,ss_temperature) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  as.data.frame() 

asv_sample_wide_merge_climate <- asv_sample_merge_sst%>%
  mutate( climate = case_when((ss_temperature > 18 ) ~ "tropical",
                              (ss_temperature < 10 ) ~ "polar",
                              (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) %>%
  na.omit()



# Checking the distribution of climate type

# Create base map
worldmap <- map_data ("world")# From the tidyverse package

base_world  <- ggplot()+ 
  coord_fixed() +
  xlab("") + ylab("") + 
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group), 
                fill="grey80")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), #legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

# Plot distribution of climate type 
base_world +
  geom_jitter(data = asv_sample_wide_merge_climate, aes(x = longitude, y = latitude,color = climate), size = 3, 
             shape = 16)+
  guides(color=guide_legend(title="Climate"))+
  scale_color_viridis_d()+
  theme_minimal()+
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        # legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()
        )

ggsave(filename=paste0(wd,"/PLOTS/Pre/Climate_type_map.png"),width=40,height=19,units=c("cm"),bg="white")

# Prepare data for salinity 
woa23_sali <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_s00mn01.csv'))

# Coalesce WOA data, maximum depth of 100m and extract the salinty at 0m
woa23_sali_coalesce<- cbind(woa23_sali[,1:2],
                                 salinity = as.data.frame(coalesce(!!!woa23_sali[,3:length(woa23_sali)])))

names(woa23_sali_coalesce)[3] <- "salinity"

# Find which grid closest to the point has a filled in value and extract that value
woa23_sali_coalesce <- woa23_sali_coalesce[,c("longitude", "latitude", "salinity")]

s_min_dist<- cbind(asv_sample_wide_merge_location[,2:3],index=data.frame(apply(raster::pointDistance(as.matrix(asv_sample_wide_merge_location[,2:3]), 
                                                                                                   as.matrix(woa23_sali_coalesce[,1:2]), 
                                                                                                   lonlat = T), 1,which.min)))

names(s_min_dist)<- c("longitude",'latitude',"index")

woa23_sali_coalesce$index <- 1:nrow(woa23_sali_coalesce)

s_min_dist<- merge(s_min_dist, woa23_sali_coalesce, by = 'index', all.x = TRUE)

# Combining sst to file code and aggregating repeats
asv_sample_woa_sali<- asv_sample_wide_merge_location %>%
  merge(., s_min_dist[c(2,3,6)], by.x = c('longitude','latitude'),
        by.y = c("longitude.x","latitude.x"))%>%
  dplyr::select(file_code,longitude,latitude,salinity) %>%
  dplyr::group_by(file_code,longitude,latitude,salinity) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  dplyr::rename("salinity_woa"="salinity")%>%
  as.data.frame() 

# Adding sali to original dataset
asv_sample_wide<-merge(asv_sample_wide, asv_sample_woa_sali[,c("file_code","salinity_woa")], by = c("file_code"), all.x = T)

# Filling NA values in salinity using WOA data
asv_sample_wide<- cbind(asv_sample_wide,
                        as.data.frame(coalesce(!!!asv_sample_wide[,c("salinity","salinity_woa")])))

names(asv_sample_wide)[50] <- "coalesce_salinity"

# Removing extreme values of salinity with NA >40
asv_sample_wide$salinity[asv_sample_wide$salinity > 40] <- NA
asv_sample_wide$salinity_woa[asv_sample_wide$salinity_woa > 40] <- NA
asv_sample_wide$coalesce_salinity[asv_sample_wide$coalesce_salinity > 40] <- NA


## Export updated asv_sample_wide dataset to csv 
write.csv(asv_sample_wide, paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'), row.names=FALSE)

####  Adding edited values to samples dataset

# Load samples dataset
samples_df <- fread(paste0(wd,'/DATA/ALL/samples.csv'))

# Adding sst to original dataset
samples_df<-merge(samples_df, asv_sample_woa_sst[,c("file_code","ss_temperature")], by = c("file_code"), all.x = T)


# replacing climate column using annual temp data
samples_df <- mutate(samples_df, climate = case_when((ss_temperature > 18 ) ~ "tropical",
                                                               (ss_temperature < 10 ) ~ "polar",
                                                               (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) 

samples_df<-merge(samples_df, asv_sample_woa_sali[,c("file_code","salinity_woa")], by = c("file_code"), all.x = T)


# Filling NA values in salinity using WOA data
samples_df<- cbind(samples_df,
                        as.data.frame(coalesce(!!!samples_df[,c("salinity","salinity_woa")])))

names(samples_df)[50] <- "coalesce_salinity"

# Removing extreme values of salinity with NA >40
samples_df$salinity[samples_df$salinity > 40] <- NA
samples_df$salinity_woa[samples_df$salinity_woa > 40] <- NA
samples_df$coalesce_salinity[samples_df$coalesce_salinity > 40] <- NA

## Export updated sample dataset to csv 
write.csv(samples_df, paste0(wd,'/DATA/ALL/samples_edit.csv'), row.names=FALSE)
