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
library(gstat)


## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

asv_sample_wide<-read.table(file = paste0(wd,'/DATA/ALL/metapr2_ASVs_selected_abundance_Eukaryota.tsv'), sep = '\t', header = TRUE)

asv_sample_wide <- as.data.frame(unclass(asv_sample_wide),  # Convert all columns to factor
                       stringsAsFactors = TRUE)

woa23_1deg_temp <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_t00mn01.csv'))

# Prepare data for sea surface temperature

woa23_1deg_ss_temp <- cbind(woa23_1deg_temp[,1:2], coalesce(!!!woa23_1deg_temp[,3:23])) # Fill with non-missing vector until 100m 

colnames(woa23_1deg_ss_temp) <- c("latitude","longitude","ss_temperature") # Renaming the columns to match sample data

round_05<- function(value){
  value_1<- value + 0.5 
  value_2<- round(value_1,0)
  value_3<- value_2 -0.5
  
  return(value_3)
}

# Merging down to 3128 samples
asv_sample_wide_merge_location <- asv_sample_wide %>%
  dplyr::select(file_code,latitude,longitude) %>%
  dplyr::group_by(file_code,latitude,longitude) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  as.data.frame() 

# To find out the 0.5 rounded off lat and long for each sample

asv_sample_wide_merge_location[2:3] <- lapply(asv_sample_wide_merge_location[2:3], round_05)

# Adding sst column to asv_sample_wide_merge_location dataset
asv_sample_wide_sst<-merge(asv_sample_wide_merge_location, woa23_1deg_ss_temp, by = c("latitude", "longitude"), all.x = T)


# Plotting location of NA values
asv_sample_wide_sst_na <- asv_sample_wide_sst[rowSums(is.na(asv_sample_wide_sst)) > 0,]

png(file=paste0(wd,"/PLOTS/Pre/sst_Na_location.jpeg"),width=1536,height=802) 

ggplot() + geom_raster(data = woa23_1deg_ss_temp, aes(x=longitude, y = latitude, fill=ss_temperature)) +
  coord_fixed(ratio = 1) +
  geom_jitter(data=asv_sample_wide_sst_na,aes(x=longitude, y = latitude),size=1,shape=4)+
  scale_fill_viridis() +
  theme_bw()

dev.off()

# To solve the issue, interpolate WOA ss values
# Create an empty raster where the extend is that of 
longitude <- seq(min(woa23_1deg_ss_temp$longitude), max(woa23_1deg_ss_temp$longitude), 1)
latitude <- seq(min(woa23_1deg_ss_temp$latitude), max(woa23_1deg_ss_temp$latitude), 1)

e_ras<- expand.grid(longitude,latitude)
names(e_ras)       <- c("X", "Y")
coordinates(e_ras) <- c("X", "Y")
gridded(e_ras)     <- TRUE  # Create SpatialPixel object
fullgrid(e_ras)    <- TRUE  # Create SpatialGrid object

woa23_1deg_ss_temp<-na.omit(woa23_1deg_ss_temp)

xy <- woa23_1deg_ss_temp[,c(2,1)]

spdf <- SpatialPointsDataFrame(coords = xy, data = woa23_1deg_ss_temp,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Add P's projection information to the empty grid
proj4string(spdf) <- proj4string(spdf) # Temp fix until new proj env is adopted
proj4string(e_ras) <- proj4string(spdf)

P.idw <- gstat::idw(ss_temperature ~ 1, spdf, newdata=e_ras, idp=2.0)

woa_idw_raster<- raster(P.idw)

woa_idw_df<- as.data.frame(rasterToPoints(woa_idw_raster))

names(woa_idw_df)<- c("longitude", "latitude","ss_temperature")

# Plotting interpolated values
png(file=paste0(wd,"/PLOTS/Pre/sst_interpolate.jpeg"),width=1536,height=802) 

ggplot() + geom_raster(data = woa_idw_df, aes(x=longitude, y = latitude, fill=ss_temperature)) +
  coord_fixed(ratio = 1) +
  geom_jitter(data=asv_sample_wide_sst_na,aes(x=longitude, y = latitude),size=1,shape=4)+
  scale_fill_viridis() +
  theme_bw()

dev.off()


# Now merging with the new interpolated dataset
asv_sample_sst_interp<-merge(asv_sample_wide_merge_location, woa_idw_df, by = c("latitude", "longitude"), all.x = T)

# Adding sst to original dataset
asv_sample_wide<-merge(asv_sample_wide, asv_sample_sst_interp[,c(3,4)], by = c("file_code"), all.x = T)


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

asv_sample_wide_merge_climate <- mutate(asv_sample_merge_sst, climate = case_when((ss_temperature > 18 ) ~ "tropical",
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
  geom_jitter(data = asv_sample_wide_merge_climate, aes(x = longitude, y = latitude,color = climate), size = 3, 
             shape = 16)+
  guides(color=guide_legend(title="Climate"))+
  scale_color_viridis_d()+
  theme_minimal()+
  theme(legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        legend.position = "bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

dev.off()

# Prepare data for salinity 
woa23_1deg_sali <- fread(paste0(wd,'/DATA/WOA/woa23_decav91C0_s00mn01.csv'))

# Find depth of NA samples

# Merging down to 3128 samples
asv_sample_merge_depth <- asv_sample_wide %>%
  dplyr::select(file_code,latitude,longitude,depth,salinity) %>%
  dplyr::group_by(file_code,latitude,longitude,depth,salinity) %>%
  dplyr::summarise(total_count=n(),.groups = 'drop')%>%
  dplyr::select( -total_count) %>%
  as.data.frame()

# # Plot distribution of depth in samples
# ggplot(asv_sample_merge_depth, aes(x=log(depth))) + 
#   geom_histogram(binwidth=1)

# Coalesce WOA data, maximum depth of 100m and extract the salinty at 0m
woa23_1deg_sali_coalesce<- cbind(woa23_1deg_sali[,1:2],
                                 salinity = as.data.frame(coalesce(!!!woa23_1deg_sali[,3:length(woa23_1deg_sali)])))

names(woa23_1deg_sali_coalesce)[3] <- "salinity"

# Interpolate WOA salinity values

woa23_1deg_sali_coalesce<-na.omit(woa23_1deg_sali_coalesce)

xy_sali <- woa23_1deg_sali_coalesce[,c(2,1)]

spdf_sali <- SpatialPointsDataFrame(coords = xy_sali, data = woa23_1deg_sali_coalesce,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Add P's projection information to the empty grid
proj4string(spdf_sali) <- proj4string(spdf_sali) # Temp fix until new proj env is adopted
proj4string(e_ras) <- proj4string(spdf_sali)

P.idw_sali <- gstat::idw(salinity ~ 1, spdf_sali, newdata=e_ras, idp=2.0)

woa_idw_sali_raster<- raster(P.idw_sali)

woa_idw_sali_df<- as.data.frame(rasterToPoints(woa_idw_sali_raster))

names(woa_idw_sali_df)<- c("longitude", "latitude","salinity_woa")

# Merging new salinity data with merged dataset 
asv_sample_sali_interp<-merge(asv_sample_wide_merge_location, woa_idw_sali_df, by = c("latitude", "longitude"), all.x = T)

# Adding WOA_salinity column to asv_sample_wide dataset
asv_sample_wide<-merge(asv_sample_wide, asv_sample_sali_interp[,3:4], by = c("file_code"), all.x = T)

# Filling NA values in salinity using WOA data
asv_sample_wide<- cbind(asv_sample_wide,
                        as.data.frame(coalesce(!!!asv_sample_wide[,c("salinity","salinity_woa")])))

names(asv_sample_wide)[50] <- "coalesce_salinity"

# Removing extreme values of salinity with NA >40
asv_sample_wide$salinity[asv_sample_wide$salinity.x > 40] <- NA
asv_sample_wide$salinity_woa[asv_sample_wide$salinity.x > 40] <- NA
asv_sample_wide$coalesce_salinity[asv_sample_wide$salinity.x > 40] <- NA


## Export asv_sample_wide_sst_sali to csv 
write.csv(asv_sample_wide, paste0(wd,'/DATA/ALL/samples1.csv'), row.names=FALSE)

####################### Codes above runs smoothly and are accurate ############

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




