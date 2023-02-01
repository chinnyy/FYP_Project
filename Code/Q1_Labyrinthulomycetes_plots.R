## Setup
###Load packages & plot setup

wd<- setwd("C:/Users/chiny/Desktop/FYP/R")

library(plyr) # data wrangling
library(dplyr) # data wrangling
library(tidyverse) # data wrangling and management
library(ggplot2) # visualization
library(phyloseq) # metabarcoding
library(ape) # phylogentic tree creation 
library(cluster) # Find clusters in data
library(openxlsx) # Read xlsx files
library(vegan) # Ecological statistic tools 
library(tabula) # Tests and measures of diversity
library(goeveg) # NMDS screeplot
library(scatterpie) # Scatterpie plot

## Load dataset
asv_sample_wide<-read.table(file = paste0(wd,'/DATA/ALL/metapr2_ASVs_selected_abundance_Eukaryota.tsv'), sep = '\t', header = TRUE)

# #Filter for only biotic and abiotic variables
# envi_var<- c("fraction_name","date","depth_level","depth","substrate","latitude","longitude","climate","temperature","salinity","ice_coverage","region","ecosystem","ecosystem_climate","substrate_type")
# avs_sample_wide_f <- select(asv_sample_wide,envi_var)

# Filter for class laby 
laby_asv_wide <- filter(asv_sample_wide,class == "Labyrinthulomycetes")

# Preparing laby data: abundance within class
laby_order <- laby_asv_wide %>%
  group_by(latitude,longitude,order) %>%
  dplyr::summarise(total_count=n(), # To count the number of occurrence for each order at each location
            .groups = 'drop')%>%
  spread(key = order, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2))))%>% # Add another row that counts the abundance of all 
  as.data.frame()

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

# Plot abundance map of laby only
png(file=paste0(wd,"/PLOTS/Laby_abundance_map.jpeg"),width=800,height=600) 

base_world+
  geom_point(data=laby_order, 
             aes(x=longitude, y=latitude, 
                 colour =log(sum),  # Edit out the log() if you dont want a logged version
                 fill =log(sum),
                 size = log(sum)), 
             pch=21, alpha=I(0.7))+
  scale_color_gradient(low="yellow", high="red")+
  scale_fill_gradient(low="yellow", high="red")

dev.off()

# Plot scatterpie map of laby only
png(file=paste0(wd,"/PLOTS/Laby_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude, r=log(sum)), 
                  data=laby_order,cols=colnames(laby_order[,c(3:8)]), color=NA)+
  geom_scatterpie_legend(log(laby_order$sum), x=-160, y=-55)

dev.off()

# Preparing laby and protist data: abundance within Sagenista division
laby_stram <- asv_sample_wide %>%
  filter(division == "Sagenista")%>%
  group_by(latitude,longitude,class) %>%
  dplyr::summarise(total_count=n(), # To count the number of occurance for each order at each location
                   .groups = 'drop')%>%
  spread(key = class, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2))))%>% # Add another row that counts the abundance of all 
  as.data.frame()

# Plot scatterpie map of laby and protist data: abundance within Sagenista division
png(file=paste0(wd,"/PLOTS/Laby_w_Sagenista_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude, r=log(sum)), 
                  data=laby_stram,cols=colnames(laby_stram[,c(3:(ncol(laby_stram)-1))]), color=NA)+
  geom_scatterpie_legend(log(laby_stram$sum), x=-160, y=-55)

dev.off()

# Arrange the order of Stramenopiles supergroup in decreasing order (Probably could be combined with the chunk below)
laby_stram_ord <- asv_sample_wide %>%
  filter(supergroup == "Stramenopiles")%>%
  group_by(class) %>%
  dplyr::summarise(total_count=n(), # To count the number of occurrence for each order at each location
                   .groups = 'drop') %>%
  arrange(desc(total_count)) 

# Preparing laby and protist data: abundance within Stramenopiles supergroup (Class level)
laby_stram <- asv_sample_wide %>%
  filter(supergroup == "Stramenopiles")%>%
  group_by(latitude,longitude,class) %>%
  dplyr::summarise(total_count=n(), # To count the number of occurrence for each order at each location
                   .groups = 'drop')%>%
  filter(class %in% laby_stram_ord$class[1:10]) %>% # Keep the top 10 counts
  spread(key = class, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2))))%>% # Add another row that counts the abundance of all 
  as.data.frame()

# Plot scatterpie map of laby and protist data: abundance within Stramenopiles supergroup
png(file=paste0(wd,"/PLOTS/Laby_w_10_stramenopiles_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude, r=log(sum)), 
                  data=laby_stram,cols=colnames(laby_stram[,c(3:(ncol(laby_stram)-1))]), color=NA)+
  geom_scatterpie_legend(log(laby_stram$sum), x=-160, y=-55)

dev.off()

######### Abundance and diversity index   #########

# Preparing laby data: abundance within class 
laby_order_wide <- laby_asv_wide %>%
  group_by(latitude,longitude,date,fraction_name,depth_level,depth,substrate,climate,temperature,salinity,ice_coverage,ecosystem,ecosystem_climate,order) %>%
  dplyr::summarise(total_count=n(), # To count the number of occurrence for each order at each location
                   .groups = 'drop')%>%
  spread(key = order, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., 14:19)))%>% # Add another row that counts the abundance of all 
  as.data.frame()

# Calculating shannon's diversity index(Species Heterogeneity & evenness)
laby_order_wide$shannond_h <- heterogeneity(laby_order_wide[,14:19], method = "shannon") # Assumes random sampling
laby_order_wide$shannond_e<- evenness(laby_order_wide[,14:19], method = "shannon") # Assumes random sampling

# Calculating Brillouin diversity index(Species Heterogeneity & evenness)
laby_order_wide$brillouin_h <- heterogeneity(laby_order_wide[,14:19], method = "brillouin") # Does not assume random sampling
laby_order_wide$brillouin_e <- evenness(laby_order_wide[,14:19], method = "brillouin") # Does not assume random sampling

# Calculating Chao1 index (Species Rarefaction)
laby_order_wide$Chao1 <- composition(laby_order_wide[,14:19], method = "chao1")

# Calculating Margalef index (Species Richness)
laby_order_wide$Margalef <- richness(laby_order_wide[,14:19], method = "margalef")

# Plot diversity, rarefaction and Richness index map of laby 

# Make list of variable names to loop over.
plot_name_list = combn(names(laby_order_wide)[21:26], 1, simplify=FALSE)

# Make plots.
plot_list = list()
for (i in 1:length(plot_name_list)) {
  p = base_world+
    geom_point(data=laby_order_wide, 
                 aes(x=longitude, y=latitude, 
                     colour =shannond_h,  # Edit out the log() if you dont want a logged version
                     fill =shannond_h,
                     size = shannond_h), 
                 pch=21, alpha=I(0.7))+
    scale_color_gradient(low="yellow", high="red")+
    scale_fill_gradient(low="yellow", high="red")+
    ggtitle(paste("Scatter map of Labyrinthulomycetes according to",plot_name_list[i],"index"))+
    theme(plot.title = element_text(hjust = 0.5))
  plot_list[[i]] = p
}

# Save plots to .jpeg. Makes a separate file for each plot.
for (i in 1:length(plot_name_list)) {
#  file_name = paste("Laby_", i, ".jpeg", sep="")
#  tiff(file_name)
  png(file=paste0(wd,"/PLOTS/Laby_",plot_name_list[i],"_scatter_map.jpeg"),width=1536,height=802) 
  print(plot_list[[i]])
  dev.off()
}


