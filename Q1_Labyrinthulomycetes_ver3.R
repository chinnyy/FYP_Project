## Setup
###Load packages & plot setup


# Library is super messy but i will clean it up :")
library(plyr) # data wrangling
library(dplyr) # data wrangling
library(tidyverse) # data wrangling and management
library(ggplot2) # visualization
library(phyloseq) # metabarcoding
library(ape) # phylogentic tree creation 
library(cluster) # Find clusters in data
library(vegan) # Ecological statistic tools 
library(tabula) # Tests and measures of diversity
library(goeveg) # NMDS screeplot
library(scatterpie) # Scatterpie plot
library(smplot2) # Statistical data visualization that complements ggplot2
library(gridExtra) # Extensions to the grid system in ggplot2
library(treemapify) # Extensions to create treemaps in ggplot2
library(viridis) # Colourblind friendly map
library(rstudioapi)
library(scales)
library("xlsx")
library(ggpubr)
library(ggallin)
library(cowplot)


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

## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

asv_sample<-read.csv(file = paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'))

sample_info<- read.csv(file = paste0(wd,'/DATA/ALL/samples_edit.csv'))

## Part 0: Removing unwanted data and rearranging categories

# Remove the non-protist
asv_sample<-asv_sample[!(asv_sample$division == "Metazoa" | asv_sample$division == "Fungi"),]

# rearrange depth_level
asv_sample<-asv_sample %>%
  mutate(depth_level = case_when((depth_level == "mesopelagic" ) ~ "pelagic",
                                 (depth_level == "bathypelagic" ) ~ "pelagic",
                                 TRUE ~ depth_level))
sample_info<-sample_info %>%
  mutate(depth_level = case_when((depth_level == "mesopelagic" ) ~ "pelagic",
                                 (depth_level == "bathypelagic" ) ~ "pelagic",
                                 TRUE ~ depth_level))
## Part 1: Distribution of laby (surface)

# Create file directory for plots from part 1
dir.create(paste0(wd,'/PLOTS/Q1'))
dir.create(paste0(wd,'/PLOTS/Q1/P1'))

# For global distribution plots, only surface communities will be plotted (100m) and exclude time series data
asv_sample_surface_nts<-asv_sample %>%
  subset(depth_level == 'surface' | depth <= 100 )%>%
  subset(dataset_id!="16" & dataset_id!="18"& dataset_id!="19"& dataset_id!="20"& dataset_id!="28"& dataset_id!="29"& dataset_id!="36"& dataset_id!="49"& dataset_id!="384")


# Preparing Presence/absence data
laby_PAD <- asv_sample_surface_nts %>%
  group_by(file_code, latitude, longitude, class, n_reads) %>%
  dplyr::summarise(total_count = n() * n_reads, .groups = 'drop') %>% # Find the total count for each location
  dplyr::select(-n_reads) %>% group_by(file_code, latitude, longitude, class) %>%
  dplyr::summarise(total_count = sum(total_count), .groups = 'drop') %>% # Merging values of the same location
  spread(key = class, value = total_count) %>% # Convert to wide data
  dplyr::select(file_code, latitude, longitude,Labyrinthulomycetes)%>%
  mutate(P_A = case_when((Labyrinthulomycetes > 0) ~ "present",
                         is.na(Labyrinthulomycetes) ~ "absent"))%>%
  as.data.frame()

# Prepare data for Laby dominance map 
# Filter for class laby 
laby_asv_sample<- filter(asv_sample,class == "Labyrinthulomycetes")
laby_asv_sample_surface_nts <- filter(asv_sample_surface_nts,class == "Labyrinthulomycetes")

dominant_taxon <- laby_asv_sample_surface_nts %>%
  group_by(file_code, order,n_reads_pct) %>%
  arrange(file_code, desc(n_reads_pct)) %>%
  group_by(file_code) %>%
  dplyr::slice(1) %>%
  mutate(dominant_taxon = .data[["order"]]) %>%
  dplyr::select(file_code,latitude,longitude ,dominant_taxon,n_reads_pct)

# Merge information on P/A and dominance
dominant_merge <- laby_PAD %>%
  merge(dominant_taxon[,c("file_code","dominant_taxon","n_reads_pct")], 
        by=c('file_code'),all.x = T)

# Separate plots will be created and merged together on a separate software affinity 
main_map<- base_world +
  geom_point(data=subset(dominant_merge, is.na(dominant_taxon)), 
             aes(x = longitude, 
                 y = latitude,
                 shape='NA'), stroke=0.8, shape=4,color = "red") +
  geom_jitter(data = dominant_merge, aes(x = longitude, 
                                         y = latitude,
                                         color=dominant_taxon,
                                         size = n_reads_pct*2),
              alpha=0.6)+
  guides(size = guide_legend(reverse=T))+
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 40), 
                        breaks = c(1, 5, 20),
                        trans=pseudolog10_trans)+
  scale_color_viridis_d()+
  labs(size="% of microbial \n eukaryotes", color="Dominant taxon")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.box.background = element_rect(colour = NA))

main_map +
  geom_polygon(aes(x=c(-7,-7,36,36),y=c(30,60,60,30)),
               fill = "NA", colour = "black",linewidth = 0.4)+
  geom_polygon(aes(x=c(-98,-98,-55,-55),y=c(13,48,48,13)),
               fill = "NA", colour = "black",linewidth = 0.4)+
  geom_polygon(aes(x=c(164,164,185,185),y=c(-60,-30,-30,-60)),
               fill = "NA", colour = "black",linewidth = 0.4)+
  geom_polygon(aes(x=c(-80,-80,-50,-50),y=c(-70,-60,-60,-70)),
               fill = "NA", colour = "black",linewidth = 0.4)

ggsave(filename=paste0(wd,"/PLOTS/Q1/P1/main_map.pdf"),width=40,height=19,units=c("cm"),bg="white")

# Using these coordinates, create small individual boxes of maps and save them individually 

# Create a dataframe to store the plots coordinates
cord_df<- data.frame(plot_num = c("A","B","C","D"),
                     x1 = c(-7,-98,164,-80),
                     x2 = c(36,-55,185,-50),
                     y1 = c(30,13,-60,-70),
                     y2 = c(60,48,-30,-60))

for (i in 1:nrow(cord_df)){
  
  small_p<- main_map+
    coord_cartesian(xlim =c(cord_df[i,"x1"], cord_df[i,"x2"]), ylim = c(cord_df[i,"y1"], cord_df[i,"y2"]))+
    theme(legend.position="none")
  
  print(small_p)
  
  ggsave(filename=paste0(wd,"/PLOTS/Q1/P1/small_plot",cord_df[i,"plot_num"],".pdf"),dpi=300,bg="white")
}

# Part 2: Treemap+ scatterpie maps of each environmental variable (surface)
# Create file directory for plots from part 2
dir.create(paste0(wd,'/PLOTS/Q1/P2'))


# Function to plot treemap +scatterpie map
p2_plot<- function(variable){
  
  if (variable == "ecosystem"|variable == "climate"|variable =="coalesce_salinity"){
    df<- asv_sample %>%
      subset(depth_level == 'surface' | depth <= 100 )%>% # Keep surface points only
      filter(class == "Labyrinthulomycetes") # Filter for laby
  }  else if (variable == "depth_level"){
    df<- asv_sample %>%  # for depth_level 
      mutate(across(depth_level, factor, levels=c("surface","euphotic","pelagic")))%>%
      filter(ecosystem == "oceanic" & class == "Labyrinthulomycetes") # Filter for laby and oceanic samples 
  } 
  
  
  ## Part 2.1: Scatterpie map of laby (orders)
  
  # Preparing laby scatterpie data: Relative abundance within class
  laby_order <- df %>%
    dplyr::select(file_code,any_of(c(variable)),latitude,longitude,order,n_reads_pct) %>%
    dplyr::group_by(file_code,across(any_of(variable)),latitude,longitude,order)%>%
    dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
    dplyr::arrange(order, desc(n_reads_pct))%>%
    pivot_wider(names_from = order, values_from = n_reads_pct,values_fill = 0) %>%
    na.omit()%>%
    as.data.frame()
  
  if (variable == "coalesce_salinity"){
    laby_order <- laby_order %>% 
      mutate( coalesce_salinity = cut( coalesce_salinity, breaks = c(0,5,30,40) ))%>%
      drop_na()
  }
  
  ### Plotting scatterpie map of orders in class Laby 
  
  # Combine colour palette to dataframe
  
  for (i in unique(laby_order[,variable])){
    laby_order_f <- laby_order%>%
      subset(get(variable) == i)
  
    if (i == "polar"|i == "euphotic"|i == "pelagic"){
      laby_scatter<-base_world+
        geom_scatterpie(aes(x=longitude, y=latitude), 
                        data=laby_order_f,cols=colnames(laby_order_f[,c(5:10)]), 
                        color=NA,alpha=0.8,pie_scale = 0.3)+
        theme(legend.position = "bottom")+
        guides(fill=guide_legend(title="Order"))+
        scale_fill_manual(values=c("#440154","#2a788e","#22a884","#7ad151","#fde725"))
    } else if (i == "surface"){
      laby_scatter<-base_world+
        geom_scatterpie(aes(x=longitude, y=latitude), 
                        data=laby_order_f,cols=colnames(laby_order_f[,c(5:10)]), 
                        color=NA,alpha=0.8,pie_scale = 0.3)+
        theme(legend.position = "bottom")+
        guides(fill=guide_legend(title="Order"))+
        scale_fill_manual(values=c("#414487","#2a788e","#22a884","#7ad151","#fde725"))
    }
    else {
        laby_scatter<-base_world+
    geom_scatterpie(aes(x=longitude, y=latitude), 
                    data=laby_order_f,cols=colnames(laby_order_f[,c(5:10)]), 
                    color=NA,alpha=0.8,pie_scale = 0.3)+
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(title="Order"))+
    scale_fill_viridis_d()
    }
    


  
  print(laby_scatter)
  
  ggsave(filename=paste0(wd,"/PLOTS/Q1/P2/Laby_scatterpie_",variable,"_",i,".pdf"),width=40,height=22,units=c("cm"),bg="white")
  
  }
  ## Part 2.2: Tree map of laby order subgroup (species)
  
  # Preparing laby treemap data: order and subgroup species 
  laby_order_sub <- df %>%
    dplyr::select(any_of(c(variable)),latitude,longitude,order,species,n_reads_pct) %>%
    group_by(across(any_of(variable)),order,species)%>%
    dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
    drop_na()%>%
    as.data.frame()
  
  if (variable == "coalesce_salinity"){
    laby_order_sub <- laby_order_sub %>% 
      mutate( coalesce_salinity = cut(coalesce_salinity, breaks = c(0,5,30,40) ))%>%
      drop_na()
  }
  
  # Finding average laby % and number of samples for each variable
  f_labels <- df %>%
    group_by(file_code,across(any_of(variable)))%>%
    dplyr::summarise(file_code_sum = sum(n_reads_pct),.groups = 'drop')%>% 
    group_by(across(any_of(variable)))%>%
    dplyr::summarise(avg_reads = mean(file_code_sum),total_count=n(),.groups = 'drop')%>%
    drop_na()%>%
    as.data.frame()
  
  ### Plotting treemap of orders and subgroup in class laby 

  for (i in unique(laby_order_sub[,variable])){
    laby_order_sub_f <- laby_order_sub%>%
      subset(get(variable) == i)
    
    if (i == "polar"|i == "euphotic"|i == "pelagic"){
      treemap<- ggplot(laby_order_sub_f, aes(area = n_reads_pct, 
                               fill = order,
                               subgroup = order, 
                               label = species)) +
      treemapify::geom_treemap() +
      treemapify::geom_treemap_text(colour = "white", place = "centre", grow = T) +
      treemapify::geom_treemap_subgroup_border() +
      treemapify::geom_treemap_subgroup_text(place = "topleft", grow = T, 
                                             alpha = 0.5, colour = "black", 
                                             min.size = 0) +
      theme_bw() +
      scale_fill_manual(values=c("#440154","#2a788e","#22a884","#7ad151","#fde725"))+
      coord_cartesian(xlim = c(50, 350), ylim = c(10, 35), clip = "off")+
      guides(fill = "none") 
    } else if (i == "surface") {
      treemap<- ggplot(laby_order_sub_f, aes(area = n_reads_pct, 
                                             fill = order,
                                             subgroup = order, 
                                             label = species)) +
        treemapify::geom_treemap() +
        treemapify::geom_treemap_text(colour = "white", place = "centre", grow = T) +
        treemapify::geom_treemap_subgroup_border() +
        treemapify::geom_treemap_subgroup_text(place = "topleft", grow = T, 
                                               alpha = 0.5, colour = "black", 
                                               min.size = 0) +
        theme_bw() +
        scale_fill_manual(values=c("#414487","#2a788e","#22a884","#7ad151","#fde725"))+
        coord_cartesian(xlim = c(50, 350), ylim = c(10, 35), clip = "off")+
        guides(fill = "none") 
    }
    else{      
      treemap<- ggplot(laby_order_sub_f, aes(area = n_reads_pct, 
                                                        fill = order,
                                                        subgroup = order, 
                                                        label = species)) +
      treemapify::geom_treemap() +
      treemapify::geom_treemap_text(colour = "white", place = "centre", grow = T) +
      treemapify::geom_treemap_subgroup_border() +
      treemapify::geom_treemap_subgroup_text(place = "topleft", grow = T, 
                                             alpha = 0.5, colour = "black", 
                                             min.size = 0) +
      theme_bw() +
      scale_fill_viridis_d()+
      coord_cartesian(xlim = c(50, 350), ylim = c(10, 35), clip = "off")+
      guides(fill = "none") 
    }
    
    
    
    print(treemap)
    
    ggsave(filename=paste0(wd,"/PLOTS/Q1/P2/Laby_treemap_",variable,"_",i,".pdf"),width=40,height=40,units=c("cm"),bg="white")
    
  }
  return(as.data.frame(f_labels))
}

# Using the function, plot the different treescatterpie map according to the variables 
variable_df <- data.frame(variable=c("ecosystem","depth_level","climate"))

f_labels_combi <- list() #  New dataframe to store information of plots

for (i in variable_df[,1]){
  f_lab<- p2_plot(i)
  f_labels_combi[[i]] <- f_lab # add it to your list
}

## Part 2.3: NMDS 
# Create file directory for plots from part 2
dir.create(paste0(wd,'/PLOTS/Q1/P3'))

# Preparing NDMS data: Relative abundance within class (surface)
laby_nmds_surface <- laby_asv_sample %>%
  subset(depth_level == 'surface' | depth <= 100 )%>%
  dplyr::select(file_code,asv_code,n_reads_pct) %>%
  group_by(file_code,asv_code)%>%
  pivot_wider(names_from = asv_code, values_from = n_reads_pct,values_fill = 0) %>%
  mutate_at(vars("69fe46f018":"b95e00fe3c"),as.numeric)%>%
  as.data.frame()

# Run NMDs
NMDS_laby <- metaMDS(laby_nmds_surface[,-c(1)], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# Extract NMDS scores for plotting with ggplot2
data.scores = as.data.frame(scores(NMDS_laby)$sites)

data.scores$file_code = laby_nmds_surface$file_code
data.scores<- merge(data.scores,sample_info[,c("file_code","climate","ecosystem")],by="file_code")


ggplot(data.scores, aes(x = NMDS1, y = NMDS2,label = file_code)) + 
  geom_jitter(size = 4, aes( colour = climate,shape=ecosystem))+ 
  # coord_cartesian(xlim =c(-1,1), ylim = c(-1, 1))+ # SETTING AXIS LIMIT 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.key=element_blank()) + 
  scale_color_viridis_d()

ggsave(filename=paste0(wd,"/PLOTS/Q1/P3/NMDS_surface.pdf"),width=40,height=40,units=c("cm"),bg="white")

# Running anosim
surface_com = as.matrix(laby_nmds_surface[,2:ncol(laby_nmds_surface)])#make community matrix

surface_ano = anosim(surface_com, data.scores$climate, distance = "bray", permutations = 20)

# Attempt 2 after normalizing data

# Run NMDs after normalizing data 

# Function to normalize data
fun_norm<- function(data) {
  
  # Normalize all columns
  laby_norm<- cbind(data["file_code"],round(data[,2:(ncol(data)-2)]/data[,"factor"],digit=0))
  return(laby_norm)
}

# Creating dataframe with n_reads instead of n_reads_pct 
laby_nmds_surface_norm <- laby_asv_sample %>%
  subset(depth_level == 'surface' | depth <= 100 )%>%
  dplyr::select(file_code,asv_code,n_reads) %>%
  group_by(file_code,asv_code)%>%
  pivot_wider(names_from = asv_code, values_from = n_reads,values_fill = 0) %>%
  as.data.frame()

# Normalizing data
laby_nmds_surface_norm<-   laby_nmds_surface_norm %>%
  mutate(sum = rowSums(select(., -c(1))))%>% # Add another row that counts the abundance of all
  mutate(factor = sum / mean(sum))%>%
  fun_norm(.)

# Run NMDs on normalised data
NMDS_laby_norm <- metaMDS(laby_nmds_surface_norm[,-c(1)], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# Extract NMDS scores for plotting with ggplot2
data.scores_norm = as.data.frame(scores(NMDS_laby_norm)$sites)

data.scores_norm$file_code = laby_nmds_surface_norm$file_code
data.scores_norm<- merge(data.scores_norm,sample_info[,c("file_code","climate","ecosystem")],by="file_code")


ggplot(data.scores_norm, aes(x = NMDS1, y = NMDS2,label = file_code)) + 
  geom_jitter(size = 4, aes( colour = climate,shape=ecosystem))+ 
  # coord_cartesian(xlim =c(-1,1), ylim = c(-1, 1))+ # SETTING AXIS LIMIT 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.key=element_blank()) + 
  scale_color_viridis_d()


ggsave(filename=paste0(wd,"/PLOTS/Q1/P3/NMDS_surface_norm.pdf"),width=40,height=40,units=c("cm"),bg="white")

#######################
# Data below not yet tested do not run

# Preparing NDMS data: Relative abundance within class (depth)
laby_order_depth <- laby_asv_sample %>%
  dplyr::select(file_code,depth_level,order,family,genus,species,asv_code,n_reads_pct) %>%
  group_by(file_code,depth_level,order,family,genus,species,asv_code)%>%
  dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
  arrange(order, desc(n_reads_pct))%>%
  pivot_wider(names_from = asv_code, values_from = n_reads_pct,values_fill = 0) %>%
  na.omit()%>%
  as.data.frame()

# Run NMDs
NMDS_laby_depth <- metaMDS(laby_order_depth[,-c(1:6)], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# Extract NMDS scores for plotting with ggplot2
data.scores_depth = as.data.frame(scores(NMDS_laby_depth)$sites)

data.scores_depth$order = laby_order_depth$order
data.scores_depth$family = laby_order_depth$family
data.scores_depth$genus = laby_order_depth$genus
data.scores_depth$species = laby_order_depth$species
data.scores_depth$depth_level = laby_order_depth$depth_level


data.scores_depth.f<- filter(data.scores_depth,order == "Thraustochytrida")

ggplot(data.scores_depth.f, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( colour = depth_level))+ 
  # stat_ellipse(aes(x = NMDS1,
  #                  y = NMDS2,
  #                  colour = order),
  #              type = "norm")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  #  labs(x = "NMDS1", colour = "climate", y = "NMDS2", shape = "ecosystem")  + 
  facet_wrap(~species)+
  scale_color_viridis_d()


