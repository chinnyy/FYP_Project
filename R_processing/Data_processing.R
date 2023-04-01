## Setup
###Load packages & plot setup


library(plyr) # data wrangling
library(dplyr) # data wrangling
library(tidyverse) # data wrangling and management
library(data.table) # To read large csv
library(tidyverse) # data wrangling and management
library(ggplot2) # visualization
library(viridis) # Color palette visualization
library(raster) # Managing spatial data
library(readxl) # Read excel files

## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

# Navigate up to parent directory
setwd("..")
wd <- getwd()

# Load samples dataset and WOA temperature and salinity dataset
samples_df <- read_excel(paste0(wd,'/DATA/METAPR2/samples.xlsx'))

woa23_temp <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_t00mn01.csv'))

woa23_sali <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_s00mn01.csv'))

colnames(woa23_temp)[1:3] <- c("latitude", "longitude","0")

colnames(woa23_sali)[1:3] <- c("latitude", "longitude","0")

# Prepare data for sea surface temperature

woa23_ss_temp <- cbind(woa23_temp[,1:2], coalesce(!!!woa23_temp[,3:23])) # Fill with non-missing vector until 100m 

colnames(woa23_ss_temp) <- c("latitude","longitude","ss_temperature") # Renaming the columns to match sample data

round_05<- function(value){
  value_1<- value + 0.5 
  value_2<- round(value_1,0)
  value_3<- value_2 -0.5
  
  return(value_3)
}

# To find out the 0.5 rounded off lat and long for each sample
samples_df_round<- samples_df%>%
  dplyr::select(file_code,longitude,latitude) %>%
  as.data.frame()

samples_df_round[,c("latitude","longitude")] <- lapply(samples_df[,c("latitude","longitude")], round_05)

# Adding sst column to asv_sample_wide_merge_location dataset
samples_df_sst<-merge(samples_df_round, woa23_ss_temp, by = c("latitude", "longitude"), all.x = T)

# Filtering for NA SST
samples_df_sst_na <- samples_df_sst[rowSums(is.na(samples_df_sst)) > 0,]

# Create file directory for plots from R_processing
dir.create(paste0(wd,'/PLOTS/R_processing'))

# Plotting location of NA values
ggplot() + geom_raster(data = woa23_ss_temp, aes(x=longitude, y = latitude, fill=ss_temperature)) +
  coord_fixed(ratio = 1) +
  geom_jitter(data=samples_df_sst_na,aes(x=longitude, y = latitude),size=1,shape=4,color="red")+
  labs(fill="Sea surface \ntemperature (°c)") +
  scale_fill_viridis() +
  theme_bw()

ggsave(filename=paste0(wd,"/PLOTS/R_processing/sst_Na_location.pdf"),width=40,height=19,units=c("cm"))


# To solve the issue, find which grid closest to the point has a filled in value and extract that value
woa23_ss_temp <- woa23_ss_temp[,c("longitude", "latitude", "ss_temperature")]
  

min_dist<- cbind(samples_df[,c("longitude","latitude")],index=data.frame(apply(raster::pointDistance(as.matrix(samples_df[,c("longitude","latitude")]), 
                            as.matrix(woa23_ss_temp[,c("longitude","latitude")]), 
                            lonlat = T), 1,which.min)))

names(min_dist)<- c("longitude",'latitude',"index")

woa23_ss_temp$index <- 1:nrow(woa23_ss_temp)

# Merge WOA data to lat and long data using indexing
min_dist_sst<- min_dist %>%
  merge(woa23_ss_temp, by = 'index', all.x = TRUE)

# Create a dataframe with only file_code, lat and long data
samples_df_location<- samples_df %>%
  dplyr::select(file_code,latitude,longitude)

# Combining sst to sample data using lat and long data and aggregate
asv_sample_woa_sst_location<- samples_df_location %>%
  merge(., min_dist_sst[,c(2,3,6)], by.x = c('longitude','latitude'),
        by.y = c("longitude.x","latitude.x"),all.x = T,allow.cartesian=TRUE)%>%
  dplyr::group_by(file_code,longitude,latitude,ss_temperature) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  as.data.frame() 

# Adding sst to original dataset
samples_df<-merge(samples_df, asv_sample_woa_sst_location[,c("file_code","ss_temperature")], by = c("file_code"), all.x = T)


# replacing climate column using annual temp data
samples_df <- mutate(samples_df, climate = case_when((ss_temperature > 18 ) ~ "tropical",
                                                               (ss_temperature < 10 ) ~ "polar",
                                                               (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) 

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
  geom_jitter(data = subset(samples_df, !is.na(climate)), aes(x = longitude, y = latitude,color = climate), size = 3, 
             shape = 16)+
  guides(color=guide_legend(title="Climate"))+
  scale_color_viridis_d()+
  theme_minimal()+
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

ggsave(filename=paste0(wd,"/PLOTS/R_processing/Climate_type_map.pdf"),width=40,height=19,units=c("cm"),bg="white")

# To edit salinity data

# Coalesce WOA data, maximum depth of 100m and extract the salinty at 0m
woa23_sali_coalesce<- cbind(woa23_sali[,1:2],
                                 salinity = as.data.frame(coalesce(!!!woa23_sali[,3:length(woa23_sali)])))

names(woa23_sali_coalesce)[3] <- "salinity"

# Find which grid closest to the point has a filled in value and extract that value
woa23_sali_coalesce <- woa23_sali_coalesce[,c("longitude", "latitude", "salinity")]

s_min_dist<- cbind(samples_df[,c("longitude","latitude")],
                   index=data.frame(apply(raster::pointDistance(as.matrix(samples_df[,c("longitude","latitude")]), 
                                                                as.matrix(woa23_sali_coalesce[,c("longitude","latitude")]), 
                                                                lonlat = T), 1,which.min)))

names(s_min_dist)<- c("longitude",'latitude',"index")

woa23_sali_coalesce$index <- 1:nrow(woa23_sali_coalesce)

s_min_dist<- s_min_dist%>%
  merge(woa23_sali_coalesce, by = 'index', all.x = TRUE)

# Combining salinity to sample data using lat and long data and aggregate
asv_sample_woa_sali_location<- samples_df_location %>%
  merge(., s_min_dist[c(2,3,6)], by.x = c('longitude','latitude'),
        by.y = c("longitude.x","latitude.x"),all.x = T,allow.cartesian=TRUE)%>%
  dplyr::group_by(file_code,longitude,latitude,salinity) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  dplyr::rename("salinity_woa"="salinity")%>%
  as.data.frame() 

# Adding sali to original dataset
samples_df<-merge(samples_df, asv_sample_woa_sali_location[,c("file_code","salinity_woa")], by = c("file_code"), all.x = T)


# Filling NA values in salinity using WOA data
samples_df<- cbind(samples_df,
                        as.data.frame(coalesce(!!!samples_df[,c("salinity","salinity_woa")])))

names(samples_df)[50] <- "coalesce_salinity"

# Removing extreme values of salinity with NA >40
samples_df$salinity[samples_df$salinity > 40] <- NA
samples_df$salinity_woa[samples_df$salinity_woa > 40] <- NA
samples_df$coalesce_salinity[samples_df$coalesce_salinity > 40] <- NA


## Export updated sample dataset to csv 
write.csv(samples_df, paste0(wd,'/DATA/METAPR2/samples_edit.csv'), row.names=FALSE)
