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
library(rstudioapi) #Access the RStudio API
library(scales) # Customize axis and legend labels
library(ggpubr)
library(ggallin)
library(cowplot) # Supplement ggplot2

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

# Navigate up to parent directory
setwd("..")
wd <- getwd()

asv_sample_raw<-read_tsv(file = paste0(wd,'/DATA/METAPR2/woa_metapr2_ASVs_selected_abundance_Eukaryota.tsv')) # Change the name of your MetaPR2 file accordingly

sample_info<- read.csv(file = paste0(wd,'/DATA/METAPR2/samples_edit.csv'))

## Part 0: Data processing

# Remove the non-mircobial eukaryotes
asv_sample_raw<-asv_sample_raw[!(asv_sample_raw$division == "Metazoa" | asv_sample_raw$division == "Fungi"),]

# Merging information from updated samples data set to the asv data set by matching file_code

asv_sample<- asv_sample_raw %>%
  dplyr::select(-c("climate"))%>%
  merge(sample_info[,c("file_code","climate","ss_temperature","salinity_woa","coalesce_salinity")],by=c('file_code'),all.x = T)


# Reassign depth level to avoid overlap
asv_sample<-asv_sample %>%
  mutate(depth_level = case_when((depth >=  0 & depth <= 10) ~ "surface",
                                 (depth >  10 & depth <= 200 ) ~ "euphotic",
                                 (depth >  200) ~ "pelagic",
                                 TRUE ~ depth_level))
sample_info<-sample_info %>%
  mutate(depth_level = case_when((depth >=  0 & depth <= 10) ~ "surface",
                                 (depth >  10 & depth <= 200 ) ~ "euphotic",
                                 (depth >  200) ~ "pelagic",
                                 TRUE ~ depth_level))

# Correct dataset ID 385 as it was mistakenly classified as oceanic
asv_sample<- asv_sample %>%
  mutate(ecosystem = case_when((dataset_id ==  385) ~ "coastal",
                                 TRUE ~ ecosystem))

sample_info<- sample_info %>%
  mutate(ecosystem = case_when((dataset_id ==  385) ~ "coastal",
                               TRUE ~ ecosystem))

## Part 1: Distribution of laby (surface)

# Create file directory for plots from part 1
dir.create(paste0(wd,'/PLOTS/R'))
dir.create(paste0(wd,'/PLOTS/R/P1'))

# For global distribution plots, only surface communities will be plotted and exclude time series data
asv_sample_surface_nts<-asv_sample %>%
  subset(depth_level == 'surface')%>%
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

# Creating colour palette for the 6 orders
vir_order<- c("Amphifilida" = "#440154", "Amphitremida" = "#414487","Labyrinthulida" = "#2a788e","Labyrinthulomycetes_X" = "#22a884","Oblongichytrida" = "#7ad151","Thraustochytrida" = "#fde725")

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
  scale_color_manual(values=vir_order)+
  labs(size="% of microbial \n eukaryotes", color="Dominant taxon")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.box.background = element_rect(colour = NA))

main_map_plot<- main_map +
  geom_polygon(aes(x=c(-7,-7,36,36),y=c(30,60,60,30)),
               fill = "NA", colour = "black",linewidth = 0.4)+
  geom_polygon(aes(x=c(-98,-98,-55,-55),y=c(13,48,48,13)),
               fill = "NA", colour = "black",linewidth = 0.4)+
  geom_polygon(aes(x=c(164,164,185,185),y=c(-60,-30,-30,-60)),
               fill = "NA", colour = "black",linewidth = 0.4)+
  geom_polygon(aes(x=c(-80,-80,-50,-50),y=c(-70,-60,-60,-70)),
               fill = "NA", colour = "black",linewidth = 0.4)


ggsave(filename=paste0(wd,"/PLOTS/R/P1/main_map.pdf"),plot=main_map_plot,width=40,height=19,units=c("cm"),bg="white")

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
    theme(legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))

  ggsave(filename=paste0(wd,"/PLOTS/R/P1/small_plot",cord_df[i,"plot_num"],".pdf"),plot = small_p,dpi=300,bg="white")
}

# Part 2: Treemap+ scatterpie maps of each environmental variable
# Create file directory for plots from part 2
dir.create(paste0(wd,'/PLOTS/R/P2'))


# Function to plot treemap +scatterpie map
p2_plot<- function(variable){
  
  if (variable == "ecosystem"|variable == "climate"|variable =="coalesce_salinity"){
    df<- laby_asv_sample %>%
      subset(depth_level == 'surface' )# Keep surface points only
  }  else if (variable == "depth_level"){
    df<- laby_asv_sample %>%  # for depth_level 
      mutate(across(depth_level, factor, levels=c("surface","euphotic","pelagic")))%>%
      filter(ecosystem == "oceanic") # Filter for laby and oceanic samples 
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
  for (i in unique(laby_order[,variable])){
    laby_order_f <- laby_order%>%
      subset(get(variable) == i)
  
    laby_scatter<-base_world+
      geom_scatterpie(aes(x=longitude, y=latitude), 
                      data=laby_order_f,cols=colnames(laby_order_f[,c(5:10)]), 
                      color=NA,alpha=0.8,pie_scale = 0.3)+
      theme(legend.position = "bottom")+
      guides(fill=guide_legend(title="Order"))+
      scale_fill_manual(values=vir_order)
  
  ggsave(filename=paste0(wd,"/PLOTS/R/P2/Laby_scatterpie_",variable,"_",i,".pdf"),plot = laby_scatter,width=40,height=22,units=c("cm"),bg="white")
  
  }
  ## Part 2.2: Tree map of laby order subgroup (species)
  
  # Preparing laby treemap data: order and subgroup species 
  laby_order_sub <- df %>%
    dplyr::select(any_of(c(variable)),latitude,longitude,order,species,n_reads_pct) %>%
    group_by(across(any_of(variable)),order,species)%>%
    dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
    drop_na()%>%
    as.data.frame()
  
  # Finding how many distinct ASV are there in each species 
  laby_species_distinct <- df %>%
    dplyr::select(any_of(c(variable)),latitude,longitude,order,species,asv_code) %>%
    group_by(across(any_of(variable)),order,species)%>%
    dplyr::summarize(distinct_asv = n_distinct(asv_code),.groups = 'drop')%>%
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
      scale_fill_manual(values=vir_order)+
      coord_cartesian(xlim = c(50, 350), ylim = c(10, 35), clip = "off")+
      guides(fill = "none") 
    
    ggsave(filename=paste0(wd,"/PLOTS/R/P2/Laby_treemap_",variable,"_",i,".pdf"),plot=treemap,width=40,height=40,units=c("cm"),bg="white")
    
  }
  output_list<- list(f_labels,laby_order_sub,laby_species_distinct)
  
  return(output_list)
}

# # Using the function, plot the different treescatterpie map according to the variables 
# variable_df <- data.frame(variable=c("ecosystem","depth_level","climate"))
# 
# f_labels_combi <- list() #  New dataframe to store information of plots
# 
# for (i in variable_df[,1]){
#   f_lab<- p2_plot(i)
#   f_labels_combi[[i]] <- f_lab # add it to your list
# }

# # Prepare data for depth plot
# laby_depth<- laby_asv_sample %>%  # for depth_level
#   mutate(across(depth_level, factor, levels=c("surface","euphotic","pelagic")))%>%
#   filter(ecosystem == "oceanic") %>% # Filter for laby and oceanic samples
#   mutate(depth = fct_rev(cut(depth, breaks = c(0, 2.5, 5, 7.5, 10, 50, 100, 150,200, 1000,2000,3000,9000))))%>%
#   drop_na(depth)%>%
#   group_by(depth,order) %>%
#   summarize(n_reads_pct = sum(n_reads_pct)) %>%
#   group_by(depth) %>% 
#   mutate(n_reads_pct = n_reads_pct/sum(n_reads_pct)*100)
# 
# depth_plot<- ggplot(laby_depth)+
#   geom_col(aes (x= depth,
#                 y=n_reads_pct, 
#                 fill=order)) +
#   scale_fill_viridis_d()+
#   theme_minimal()+
#   theme(panel.border = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")+
#   scale_y_continuous(breaks = seq(0,100,20)) +
#   ylab("% of reads") + xlab("Depth (m)") +
#   coord_flip()
# 
# ggsave(filename=paste0(wd,"/PLOTS/R/P2/Laby_depth_plot.pdf"),plot=depth_plot,width=40,height=40,units=c("cm"),bg="white")
# 
# ## Part 3: NMDS 
# # Create file directory for plots from part 2
# dir.create(paste0(wd,'/PLOTS/R/P3'))

# Preparing NDMS community data (surface)

# Function to normalize data
fun_norm<- function(data) {
  
  # Find the factor to normalise each sample by
  data$factor<- data$sum/mean(data$sum)
  
  # Normalize all columns
  laby_norm<- cbind(data["file_code"],round(data[,2:(ncol(data)-2)]/data[,"factor"],digit=0))
  return(laby_norm)
}


# Function to create a dataframe of community matrix and plot NMDS 

# For variable, you can only choose surface or depth 

nmds_plot <- function (variable,data){
  
  # Create dataframe with top 100 asv for all samples 
  if (variable == "surface"){
    top_asv <- data %>%
      subset(depth_level == 'surface')%>%
      drop_na(climate)%>%
      dplyr::select(asv_code,n_reads) %>%
      group_by(asv_code)%>%
      dplyr::summarise(total = sum(n_reads),.groups = 'drop')%>% # merging values of the same asv
      arrange(desc(total))%>%
      dplyr::filter(total >= 200) 
  } else if (variable == "depth") {
    top_asv <- data %>%
      subset(ecosystem == 'oceanic')%>%
      dplyr::select(asv_code,n_reads) %>%
      group_by(asv_code)%>%
      dplyr::summarise(total = sum(n_reads),.groups = 'drop')%>% # merging values of the same asv
      arrange(desc(total))%>%
      dplyr::filter(total >= 200) 
  }
  
  # Create community matrix with top 100 asv
  if (variable == "surface"){
    laby_nmds_norm <- data %>%
      subset(depth_level == 'surface' )%>%
      drop_na(climate)%>%
      dplyr::select(file_code,asv_code,n_reads) %>%
      merge(top_asv,., by="asv_code")%>%
      dplyr::select(-total) %>%
      pivot_wider(names_from = asv_code, values_from = n_reads,values_fill = 0) %>%
      as.data.frame()
  } else if (variable == "depth"){
    laby_nmds_norm <- data %>%
      subset(ecosystem == 'oceanic')%>%
      dplyr::select(file_code,asv_code,n_reads) %>%
      merge(top_asv,., by="asv_code")%>%
      dplyr::select(-total) %>%
      pivot_wider(names_from = asv_code, values_from = n_reads,values_fill = 0) %>%
      as.data.frame()
  }
  
  # Normalizing data and filter for samples with less than a total of 50 reads
  laby_nmds_norm<-   laby_nmds_norm %>%
    mutate(sum = rowSums(dplyr::select(., -c(1))))%>% # Add another row that counts the abundance of all
    filter(sum>=50)%>%
    fun_norm(.)
  
  # Run NMDs on normalised data
  NMDS_laby_norm <- metaMDS(laby_nmds_norm[,-c(1)], k = 2, trymax = 40, 
                            trace = F, autotransform = FALSE,
                            distance = "bray")
  
  # Extract NMDS scores for plotting with ggplot2
  data.scores_norm = as.data.frame(scores(NMDS_laby_norm)$sites)
  data.scores_norm$file_code = laby_nmds_norm$file_code
  if (variable == "surface"){
    data.scores_norm<- merge(data.scores_norm,sample_info[,c("file_code","climate","ecosystem")],by="file_code")
  }else if (variable == "depth"){
    data.scores_norm<- merge(data.scores_norm,sample_info[,c("file_code","depth_level")],by="file_code")
  }
  
  # Remove extreme values from NMDs plots
  # if (variable == "surface"){
  # data.scores_norm <- data.scores_norm %>%
  #   filter(NMDS2 >= -0.4 & NMDS2 <= -0.1) %>%
  #   filter(NMDS1 >= -0.3 & NMDS1 <= 0) 
  # } else if (variable =="depth"){
  #   data.scores_norm <- data.scores_norm %>%
  #     filter(NMDS2 >= -0.4 & NMDS2 <= 0) %>%
  #     filter(NMDS1 >= -0.5 & NMDS1 <= 0) 
  # }
  NMDS1_quantiles <- quantile(data.scores_norm$NMDS1, c(0.05, 0.95)) 
  NMDS2_quantiles <- quantile(data.scores_norm$NMDS2, c(0.05, 0.95))
  data.scores_norm <- data.scores_norm[data.scores_norm$NMDS1 > NMDS1_quantiles[1] &   # Drop rows below/above percentiles
                                         data.scores_norm$NMDS1 < NMDS1_quantiles[2], ]
  data.scores_norm <- data.scores_norm[data.scores_norm$NMDS2 > NMDS2_quantiles[1] &   # Drop rows below/above percentiles
                                         data.scores_norm$NMDS2 < NMDS2_quantiles[2], ]
  # Plot ggplot
  nmds_plot<- ggplot(data.scores_norm, aes(x = NMDS1, y = NMDS2)) + #,label = file_code
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
          legend.key=element_blank()) +
    scale_color_viridis_d()
  
  
  if (variable == "surface"){
    facet_ecosystem<- nmds_plot+
      geom_jitter(size = 4, aes( colour = climate))+
      facet_wrap(~ecosystem)
    
    
    ggsave(filename=paste0(wd,"/PLOTS/R/P3/NMDS_surface_facet_ecosystem.pdf"),plot=facet_ecosystem,width=40,height=40,units=c("cm"),bg="white")
    
    facet_climate<-nmds_plot+
      geom_jitter(size = 4, aes( colour = ecosystem))+
      facet_wrap(~climate)
    
    ggsave(filename=paste0(wd,"/PLOTS/R/P3/NMDS_surface_facet_climate.pdf"),plot=facet_climate,width=40,height=40,units=c("cm"),bg="white")
    
  }else if (variable == "depth"){
    nmds_depth<-nmds_plot+
      geom_jitter(size = 4, aes( colour = depth_level))+
      guides(color=guide_legend(title="Depth level"))+ 
      scale_color_viridis_d(limits = c("surface", "euphotic", "pelagic"),direction = -1)
    
    ggsave(filename=paste0(wd,"/PLOTS/R/P3/NMDS_depth.pdf"),plot=nmds_depth,width=40,height=40,units=c("cm"),bg="white")
  }
}

# # Using the function, plot the NMDs plots according to the variables 
# nmds_variable_df <- data.frame(variable=c("surface","depth"))
# 
# for (i in nmds_variable_df[,1]){
#   nmds_com_mat(i,laby_asv_sample)
# }

############# CODE GRAVEYARD###########


# # Create dataframe with top 100 asv for all samples 
# top_asv <- laby_asv_sample %>%
#   subset(depth_level == 'surface')%>%
#   drop_na(climate)%>%
#   dplyr::select(asv_code,n_reads) %>%
#   group_by(asv_code)%>%
#   dplyr::summarise(total = sum(n_reads),.groups = 'drop')%>% # merging values of the same asv
#   arrange(desc(total))%>%
#   dplyr::slice(1:25) 
# 
# # Create community matrix with top 100 asv
# laby_nmds_surface_norm <- laby_asv_sample %>%
#   subset(depth_level == 'surface' )%>%
#   drop_na(climate)%>%
#   dplyr::select(file_code,asv_code,n_reads) %>%
#   merge(top_asv,., by="asv_code")%>%
#   dplyr::select(-total) %>%
#   pivot_wider(names_from = asv_code, values_from = n_reads,values_fill = 0) %>%
#   as.data.frame()
# 
# # Normalizing data and filter for samples with less than a total of 50 reads
# laby_nmds_surface_norm<-   laby_nmds_surface_norm %>%
#   mutate(sum = rowSums(select(., -c(1))))%>% # Add another row that counts the abundance of all
#   filter(sum>=50)%>%
#   fun_norm(.)
# 
# 
# # Obtaining environmental data
# envi_data <- laby_nmds_surface_norm %>%
#   merge(sample_info[,c("file_code","climate","ecosystem")],by = "file_code") %>%
#   select("file_code","climate","ecosystem")
# 
# 
# # Run NMDs on normalised data
# NMDS_laby_norm <- metaMDS(laby_nmds_surface_norm[,-c(1)], k = 2, trymax = 20, 
#                      trace = F, autotransform = FALSE,
#                      distance = "bray")
# 
# en = envfit(NMDS_laby_norm, envi_data[,2:3], permutations = 20, na.rm = TRUE)
# 
# 
# # Extract NMDS scores for plotting with ggplot2
# data.scores_norm = as.data.frame(scores(NMDS_laby_norm)$sites)
# 
# data.scores_norm$file_code = laby_nmds_surface_norm$file_code
# data.scores_norm<- merge(data.scores_norm,sample_info[,c("file_code","climate","ecosystem")],by="file_code")
# 
# 
# ggplot(data.scores_norm, aes(x = NMDS1, y = NMDS2)) + #,label = file_code
#   geom_jitter(size = 4, aes( colour = climate))+ 
#   # coord_cartesian(xlim =c(0.5225,0.5275), ylim = c(0.541, 0.546))+ # for combine plot filter 100
#   # coord_cartesian(xlim =c(-0.355,-0.345), ylim = c(-0.045, -0.035))+ # for combine plot filter 50
#   coord_cartesian(xlim =c(-0.51,-0.45), ylim = c(-0.135, -0.115))+ # for facet plot filter 50
#   # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
#   #              data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
#   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
#         axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
#         legend.text = element_text(size = 12, face ="bold", colour ="black"), 
#         legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
#         axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
#         legend.title = element_text(size = 14, colour = "black", face = "bold"), 
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#         legend.key=element_blank()) + 
#   facet_wrap(~ecosystem)+
#   scale_color_viridis_d()
# 
# 
# ggsave(filename=paste0(wd,"/PLOTS/R/P3/NMDS_surface_norm.pdf"),width=40,height=40,units=c("cm"),bg="white")
# 
# # Running ANOSIM
# surface_com = as.matrix(laby_nmds_surface_norm[,2:ncol(laby_nmds_surface_norm)])#make community matrix
# 
# anosim(surface_com, data.scores_norm$climate, distance = "bray", permutations = 20)
# 
# #######################
# # For depth 
# 
# # Create dataframe with top 100 asv for all samples 
# top_asv_depth <- laby_asv_sample %>%
#   subset(ecosystem == 'oceanic')%>%
#   dplyr::select(asv_code,n_reads) %>%
#   group_by(asv_code)%>%
#   dplyr::summarise(total = sum(n_reads),.groups = 'drop')%>% # merging values of the same asv
#   arrange(desc(total))%>%
#   dplyr::slice(1:25) 
# 
# # Create community matrix with top 100 asv
# laby_nmds_depth_norm <- laby_asv_sample %>%
#   subset(ecosystem == 'oceanic')%>%
#   dplyr::select(file_code,asv_code,n_reads) %>%
#   merge(top_asv_depth,., by="asv_code")%>%
#   dplyr::select(-total) %>%
#   pivot_wider(names_from = asv_code, values_from = n_reads,values_fill = 0) %>%
#   as.data.frame()
# 
# # Normalizing data
# laby_nmds_depth_norm<-   laby_nmds_depth_norm %>%
#   mutate(sum = rowSums(select(., -c(1))))%>% # Add another row that counts the abundance of all
#   filter(sum>=50)%>%
#   fun_norm(.)
# 
# # Run NMDs on normalised data
# NMDS_laby_depth_norm <- metaMDS(laby_nmds_depth_norm[,-c(1)], k = 2, trymax = 20, 
#                           trace = F, autotransform = FALSE,
#                           distance = "bray")
# 
# 
# # Extract NMDS scores for plotting with ggplot2
# data.scores_depth_norm = as.data.frame(scores(NMDS_laby_depth_norm)$sites)
# 
# data.scores_depth_norm$file_code = laby_nmds_depth_norm$file_code
# data.scores_depth_norm<- merge(data.scores_depth_norm,sample_info[,c("file_code","depth_level")],by="file_code")
# 
# 
# ggplot(data.scores_depth_norm, aes(x = NMDS1, y = NMDS2)) + #,label = file_code
#   geom_jitter(size = 4, aes( colour = depth_level))+ 
#   coord_cartesian(xlim =c(-0.05,0.04), ylim = c(-0.025, 0.025))+ # for depth filter 50
#   # coord_cartesian(xlim =c(-0.08,0.06), ylim = c(-0.06, 0.06))+ # for depth filter 100
#   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
#         axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
#         legend.text = element_text(size = 12, face ="bold", colour ="black"), 
#         legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
#         axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
#         legend.title = element_text(size = 14, colour = "black", face = "bold"), 
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#         legend.key=element_blank()) + 
#   guides(color=guide_legend(title="Depth level"))+ 
#   scale_color_viridis_d(limits = c("surface", "euphotic", "pelagic"),direction = -1)
# 
# # Obtaining environmental data
# envi_data_depth <- laby_nmds_depth_norm %>%
#   merge(sample_info[,c("file_code","depth_level")],by = "file_code") %>%
#   select("file_code","depth_level")
# 
# 
# en_depth = envfit(NMDS_laby_depth_norm, envi_data_depth[,2], permutations = 20, na.rm = TRUE)
# 
