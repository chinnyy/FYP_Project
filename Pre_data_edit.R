## Setup
###Load packages & plot setup


library(plyr) # data wrangling
library(dplyr) # data wrangling
library(tidyverse) # data wrangling and management
library(data.table) # To read large csv
library(tidyverse) # data wrangling and management
library(ggplot2) # visualization
library(sf) 

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
# Merging down to 3128 samples
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
                fill="grey80")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), #legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

# Save plot 
png(file=paste0(wd,"/PLOTS/Pre/Climate_type_map.jpeg"),width=1536,height=802) 


base_world +
  geom_jitter(data = asv_sample_wide_merge, aes(x = longitude, y = latitude,color = climate), size = 3, 
             shape = 16)+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(title="Climate"))+
  scale_color_viridis_d()

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

coalesce_salinity_df<- asv_sample_wide_sali[,c("file_code","salinity.x","salinity.y")]

coalesce_salinity_df<-cbind(file_code = asv_sample_wide_sali$file_code,
                            as.data.frame(coalesce(!!!coalesce_salinity_df[,2:3])))
  
colnames(coalesce_salinity_df)[2]  <- "coalesce_salinity" 

# Merging down to 3128 samples
coalesce_salinity_df_merge <- coalesce_salinity_df %>%
  select(file_code,coalesce_salinity) %>%
  group_by(file_code,coalesce_salinity) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  select( -total_count) %>%
  as.data.frame() 

## Add coalesce_salinity to asv_sample_wide_sst


asv_sample_wide_sst_sali<- cbind (asv_sample_wide[,1:21], # For orignal lat and long values
                                  asv_sample_wide_sst[,22:48]) # For corrected climate types and sst # For salinity

asv_sample_wide_sst_sali<- merge(asv_sample_wide_sst_sali, coalesce_salinity_df_merge, by = c("file_code"), all.x = T)


# Replacing extreme values of salinity with NA
asv_sample_wide_sst_sali$salinity[asv_sample_wide_sst_sali$salinity < 30 | asv_sample_wide_sst_sali$salinity > 40] <- NA

asv_sample_wide_sst_sali$coalesce_salinity[asv_sample_wide_sst_sali$coalesce_salinity < 30 | asv_sample_wide_sst_sali$coalesce_salinity > 40] <- NA

## Export asv_sample_wide_sst_sali to csv 
write.csv(asv_sample_wide_sst_sali, paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'), row.names=FALSE)

# For sample.xlsx data

samples_data<- read.csv(file = paste0(wd,'/DATA/ALL/samples.csv'))

#Add sst to the samples_data

sub_asv_sample_wide_merge<-asv_sample_wide_merge[,c('file_code','ss_temperature')]

samples_data_1<-merge(samples_data, sub_asv_sample_wide_merge, by = c("file_code"), all.x = T)

# Change climate type
samples_data_1 <- mutate(samples_data_1, climate = case_when((ss_temperature > 18 ) ~ "tropical",
                                                                           (ss_temperature < 10 ) ~ "polar",
                                                                           (ss_temperature >= 10 & ss_temperature <= 18) ~ "temperate")) 

# Add coalesce_salinity to samples_data

# Merging down to 3128 samples (removing extreme)
asv_sample_wide_sst_sali_merge <- asv_sample_wide_sst_sali %>%
  select(file_code,coalesce_salinity) %>%
  group_by(file_code,coalesce_salinity) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  select( -total_count) %>%
  as.data.frame() 

samples_data_1<- merge(samples_data_1, asv_sample_wide_sst_sali_merge, by = c("file_code"), all.x = T)

## Export asv_sample_wide_sst_sali to csv 
write.csv(samples_data_1, paste0(wd,'/DATA/ALL/samples1.csv'), row.names=FALSE)

#######################
## Sort samples by Longhurst codes 

# Load the shapefile
lh_shape <- read_sf(paste0(wd,'/DATA/LONGH/Longhurst_world_v4_2010.shp'))

# Extract out all the samples by file_code, lat and long 3128 samples
asv_sample<- asv_sample_wide %>%
  group_by(file_code,latitude,longitude) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>% # Find the total count for each location 
  select(-total_count)%>%
  as.data.frame()


asv_sample$region <- apply(asv_sample[,2:3], 1, function(row) {  
  # transformation to palnar is required, since sf library assumes planar projection 
  lh_shape_pl <- st_transform(lh_shape, 2163)   
  coords <- as.data.frame(matrix(row, nrow = 1, 
                                 dimnames = list("", c("latitude", "longitude"))))   
  pnt_sf <- st_transform(st_sfc(st_point(row),crs = 4326), 2163)
  # st_intersects with sparse = FALSE returns a logical matrix
  # with rows corresponds to argument 1 (points) and 
  # columns to argument 2 (polygons)
  
  lh_shape_pl[which(st_intersects(pnt_sf, lh_shape_pl, sparse = FALSE)), ]$ProvCode
})

# Rearrange values 
# Split characters
asv_sample_split <- asv_sample %>% separate(region, c('Region_1', 'Region_2','Region_3'))

# Remove redundant SANT (SANT is seen in many double reported regions)
asv_sample_split <- replace(asv_sample_split, asv_sample_split=="c", NA)
asv_sample_split <- replace(asv_sample_split, asv_sample_split=="character", NA)
asv_sample_split <- replace(asv_sample_split, asv_sample_split==0, NA)
asv_sample_split["Region_2"][asv_sample_split["Region_2"] == "SANT"] <- NA
asv_sample_split["Region_3"][asv_sample_split["Region_3"] == "SANT"] <- NA
asv_sample_split_2 <- cbind(asv_sample_split[,1:3], region = coalesce(!!!asv_sample_split[,4:6])) 

# Save plot 
png(file=paste0(wd,"/PLOTS/Pre/Region_type_map.jpeg"),width=1536,height=802) 

base_world +
  geom_jitter(data = asv_sample_split_2, aes(x = longitude, y = latitude,color = region), size = 3, 
              shape = 16)+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(title="Region"))+
  scale_color_viridis_d()

dev.off()

# Save plot after removing NA
asv_sample_split_3 <- na.omit(replace(asv_sample_split_2, asv_sample_split_2=="", NA))

png(file=paste0(wd,"/PLOTS/Pre/Region_type_NAomit_map.jpeg"),width=1536,height=802) 

base_world +
  geom_jitter(data = asv_sample_split_3, aes(x = longitude, y = latitude,color = region), size = 2, 
              shape = 16)+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(title="Region"))+
  scale_color_viridis_d()

dev.off()

# Plotting NA sample points on longhurst shp
asv_sample_split_4 <-subset(replace(asv_sample_split_2, asv_sample_split_2=="", NA),is.na(region))

png(file=paste0(wd,"/PLOTS/Pre/Region_type_NAomit_map.jpeg"),width=1536,height=802) 

ggplot() + 
  geom_polygon(data = lh_shape, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

lh_bp<- ggplot(data = lh_shape) +
  geom_sf(aes(fill = ProvCode))+
  geom_sf_text( aes(label = ProvCode))+
  scale_fill_viridis_d()+
  theme(legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), #legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

png(file=paste0(wd,"/PLOTS/Pre/Longhurst_NA_map.jpeg"),width=1536,height=802) 

lh_bp+
  geom_jitter(data = asv_sample_split_4, aes(x = longitude, y = latitude), size = 2, 
             shape = 16, color = "red")

dev.off()

################
# # Extract out all the samples by file_code, lat and long 3128 samples
# asv_sample_region<- asv_sample_wide %>%
#   group_by(file_code,latitude,longitude,region) %>%
#   dplyr::summarise(total_count=n(),.groups = 'drop')%>% # Find the total count for each location 
#   select(-total_count)%>%
#   as.data.frame()
# 
# base_world +
#   geom_jitter(data = asv_sample_region, aes(x = longitude, y = latitude,fill = region), size = 3, 
#               shape = 21)+
#   scale_fill_viridis_d()



