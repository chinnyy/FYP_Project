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

asv_sample_wide<-read.csv(file = paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'))

## Part 0: Removing unwanted data and rearranging categories

# Remove the non-protist
asv_sample_wide<-asv_sample_wide[!(asv_sample_wide$division == "Metazoa" | asv_sample_wide$division == "Fungi"),]

# rearrange depth_level
asv_sample_wide<-mutate(asv_sample_wide, 
                        depth_level = case_when((depth_level == "mesopelagic" ) ~ "pelagic",
                                                (depth_level == "bathypelagic" ) ~ "pelagic",
                                                TRUE ~ depth_level))

## Part 1: Distribution of laby (surface)

# Create file directory for plots from part 1
dir.create(paste0(wd,'/PLOTS/Q1'))
dir.create(paste0(wd,'/PLOTS/Q1/P1'))

# For global distribution plots, only surface communities will be plotted (100m)
asv_sample_wide_surface<-subset(asv_sample_wide, depth_level == 'surface' | depth <= 100 )

# Preparing Presence/absence data
laby_PAD <- asv_sample_wide_surface %>%
  group_by(file_code, latitude, longitude, class, n_reads) %>%
  dplyr::summarise(total_count = n() * n_reads, .groups = 'drop') %>% # Find the total count for each location
  select(-n_reads) %>% group_by(file_code, latitude, longitude, class) %>%
  dplyr::summarise(total_count = sum(total_count), .groups = 'drop') %>% # Merging values of the same location
  spread(key = class, value = total_count) %>% # Convert to wide data
  select(file_code, latitude, longitude,Labyrinthulomycetes)%>%
  mutate(P_A = case_when((Labyrinthulomycetes > 0) ~ "present",
                         is.na(Labyrinthulomycetes) ~ "absent"))%>%
  as.data.frame()

# Prepare data for Laby dominance map 
# Filter for class laby 
laby_asv_wide<- filter(asv_sample_wide,class == "Labyrinthulomycetes")
laby_asv_wide_surface <- filter(asv_sample_wide_surface,class == "Labyrinthulomycetes")

dominant_taxon <- laby_asv_wide_surface %>%
  arrange(file_code, desc(n_reads_pct)) %>%
  group_by(file_code) %>%
  dplyr::slice(1) %>%
  mutate(dominant_taxon = .data[["order"]]) %>%
  select(file_code,latitude,longitude ,dominant_taxon,n_reads_pct)

# Merge information on P/A and dominance
dominant_merge <- laby_PAD %>%
  merge(dominant_taxon[,c(1,4,5)], 
        by=c('file_code'),all.x = T)

png(file=paste0(wd,"/PLOTS/Q1/P1/Laby_global_map.jpeg"),width=1536,height=802) 

base_world +
  geom_jitter(data = dominant_merge, aes(x = longitude, 
                                         y = latitude,
                                         color=dominant_taxon,
                                         size = log(n_reads_pct)),
              alpha=0.5)+
  geom_point(data=subset(dominant_merge, is.na(dominant_taxon)), 
             aes(x = longitude, 
                 y = latitude,
                 shape='NA'), stroke=1.5, shape=4,color = "red") +
  guides(size = guide_legend(reverse=T))+
  scale_size(range = c(.1, 10))+
  scale_color_viridis_d()+
  labs(size="logged % of eukaryotes", color="Dominant_taxon")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.box.background = element_rect(colour = "black"))

dev.off()

# Part 2: Treemap+ scatterpie maps of each environmental variable (surface)
# Create file directory for plots from part 2
dir.create(paste0(wd,'/PLOTS/Q1/P2'))


# Function to plot treemap +scatterpie map
p2_plot<- function(variable){
  
  if (variable == "ecosystem"|variable == "climate"|variable =="coalesce_salinity"){
    df= laby_asv_wide_surface
  }  else if (variable == "depth_level"){
    df<- laby_asv_wide %>%  # for depth_level 
          mutate(across(depth_level, factor, levels=c("surface","euphotic","pelagic")))
  } 
  
  
  ## Part 2.1: Scatterpie map of laby (orders)
  
  # Preparing laby scatterpie data: Relative abundance within class
  laby_order <- df %>%
    select(file_code,any_of(c(variable)),latitude,longitude,order,n_reads_pct) %>%
    group_by(file_code,across(any_of(variable)),latitude,longitude,order)%>%
    dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
    arrange(order, desc(n_reads_pct))%>%
    pivot_wider(names_from = order, values_from = n_reads_pct) %>%
    replace(is.na(.), 0) %>%
    as.data.frame()
  
  if (variable == "coalesce_salinity"){
    laby_order <- laby_order %>% 
      mutate( coalesce_salinity = cut( coalesce_salinity, breaks = c(0,5,30,40) ))%>%
      drop_na()
  }
  
  ### Plotting scatterpie map of orders in class Laby (surface)
  laby_scatter<-base_world+
    geom_scatterpie(aes(x=longitude, y=latitude), 
                    data=laby_order,cols=colnames(laby_order[,c(5:10)]), 
                    color=NA,alpha=0.8,pie_scale = 0.3)+
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(title="Order"))+
    facet_wrap(~get(variable))+
    scale_fill_viridis_d()
  
  
  ## Part 2.2: Tree map of laby order subgroup (species)
  
  # Preparing laby treemap data: order and subgroup species (surface)
  laby_order_sub <- df %>%
    select(file_code,any_of(c(variable)),latitude,longitude,order,species,n_reads_pct) %>%
    group_by(across(any_of(variable)),order,species)%>%
    dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
    drop_na()%>%
    as.data.frame()
  
  if (variable == "coalesce_salinity"){
    laby_order_sub <- laby_order_sub %>% 
      mutate( coalesce_salinity = cut( coalesce_salinity, breaks = c(0,5,30,40) ))%>%
      drop_na()
  }
  
  # Finding average laby % for each variable
  f_labels <- laby_order_sub %>%
    group_by(across(any_of(variable)))%>%
    dplyr::summarise(avg_reads = mean(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
    drop_na()%>%
    as.data.frame()
  
  ### Plotting treemap of orders and subgroup in class laby (surface)
  
  laby_tree<- ggplot(laby_order_sub, aes(area = n_reads_pct, 
                                         fill = order,
                                         subgroup = order, 
                                         label = species)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white", place = "centre", grow = FALSE) +
    treemapify::geom_treemap_subgroup_border() +
    treemapify::geom_treemap_subgroup_text(place = "topleft", grow = F, 
                                           alpha = 0.5, colour = "black", 
                                           min.size = 0) +
    theme_bw() +
    title()
    scale_fill_viridis_d()+
    facet_wrap(~get(variable))+
    coord_cartesian(xlim = c(50, 350), ylim = c(10, 35), clip = "off")+
    guides(fill = "none") 
  
  
  png(file=paste0(wd,"/PLOTS/Q1/P2/Laby_",variable,"_treescatter_map.jpeg"),width=938,height=794) 
  
  treescatter<- ggarrange(laby_tree,laby_scatter, nrow = 2) 
  
  print(treescatter)
  
  dev.off()
  
}

# Using the function, plot the different treescatterpie map according to the variables 
variable_df <- data.frame(variable=c("ecosystem","depth_level","climate","coalesce_salinity"))
for (i in variable_df[,1]){
  p2_plot(i)
}

## Part 2.3: NMDS 
# Create file directory for plots from part 2
dir.create(paste0(wd,'/PLOTS/Q1/P3'))

# Attempt 1

# Preparing NDMS data: Relative abundance within class (surface)
laby_order_surface <- laby_asv_wide_surface %>%
  select(file_code,ecosystem,climate,coalesce_salinity,latitude,longitude,order,n_reads_pct) %>%
  group_by(file_code,ecosystem,climate,coalesce_salinity,latitude,longitude,order)%>%
  dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
  arrange(order, desc(n_reads_pct))%>%
  pivot_wider(names_from = order, values_from = n_reads_pct) %>%
  replace(is.na(.), 0) %>%
  mutate( coalesce_salinity = as.character(cut( coalesce_salinity, breaks = c(0,5,30,40) )))%>%
  na.omit()%>%
  as.data.frame()

# Run NMDs
NMDS_laby <- metaMDS(laby_order_surface[,-c(1:6)], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# Extract NMDS scores for plotting with ggplot2
## species scores
species.scores <- as.data.frame(scores(NMDS_laby, "species"))
species.scores$species <- rownames(species.scores) 
species.scores$type <- "Species"

## site scores
site.scores <- as.data.frame(scores(NMDS_laby, "sites"))
site.scores <- cbind (site.scores, 
                      climate = laby_order_surface$climate, 
                      coalesce_salinity = laby_order_surface$coalesce_salinity, 
                      ecosystem = laby_order_surface$ecosystem) # Add environmental factors

plot_nmds <- function(species.scores,site.scores,variable,round){
  plot<-ggplot() + 
    geom_point(data = species.scores,
               aes(x = NMDS1,
                   y = NMDS2),
               size = 2) +
    geom_point(data = site.scores,
               aes(x = NMDS1,
                   y = NMDS2,
                   colour = get(variable)),
               shape = 17,
               size = 2) +
    # classify sites by variables and display using 95% confidence ellipses
    stat_ellipse(data = site.scores,
                 aes(x = NMDS1,
                     y = NMDS2,
                     colour = get(variable)),
                 type = "norm") +
    labs(color=variable) +
    scale_color_viridis_d()+
    ggtitle(paste0("NMDS plot of laby and ",variable)) +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size=22))
  
  # Saving plots
  png(file=paste0(wd,"/PLOTS/Q1/P3/NMDS_plot_",round,"_laby_",variable,".jpeg"),width=1536,height=802) 
  print(plot)
  dev.off()
}

# Plot NMDS
for (i in c("climate","ecosystem","coalesce_salinity")){
  plot_nmds(species.scores,site.scores,i,1)
}

# Preparing NDMS data: Relative abundance within class (depth)
laby_order_depth <- laby_asv_wide %>%
  select(file_code,depth_level,latitude,longitude,order,n_reads_pct) %>%
  group_by(file_code,depth_level,latitude,longitude,order)%>%
  dplyr::summarise(n_reads_pct = sum(n_reads_pct),.groups = 'drop')%>% # merging values of the same file_code
  arrange(order, desc(n_reads_pct))%>%
  pivot_wider(names_from = order, values_from = n_reads_pct) %>%
  replace(is.na(.), 0) %>%
  as.data.frame()

# Run NMDs
NMDS_laby_depth <- metaMDS(laby_order_depth[,-c(1:4)], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# Extract NMDS scores for plotting with ggplot2
## species scores
species.scores.depth <- as.data.frame(scores(NMDS_laby_depth, "species"))
species.scores.depth$species <- rownames(species.scores.depth) 
species.scores.depth$type <- "Species"

## site scores
species.scores.depth <- as.data.frame(scores(NMDS_laby_depth, "sites"))
species.scores.depth <- cbind (species.scores.depth, 
                      depth_level = laby_order_depth$depth_level) # Add environmental factors

# Plot NMDS
plot_nmds(species.scores.depth,species.scores.depth,"depth_level",1)

# Attempt 2

# Preparing NDMS data: Relative abundance within class (surface)
laby_order_2 <- laby_asv_wide_surface %>%
  select(file_code,ecosystem,climate,coalesce_salinity,latitude,longitude,order,n_reads) %>%
  group_by(file_code,ecosystem,climate,coalesce_salinity,latitude,longitude,order)%>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location
  group_by(file_code,ecosystem,climate,coalesce_salinity,latitude,longitude,order)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = order, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate( coalesce_salinity = as.character(cut( coalesce_salinity, breaks = c(0,5,30,40) )))%>%
  mutate(sum = rowSums(select(., -c(1:6))))%>% # Add another row that counts the abundance of all
  na.omit()%>%
  as.data.frame()

# Function to remove bottom 5% of the dataset and normalize all columns
fun_rem_5<- function(data,col_start) {
  p05 <- quantile(data$sum, 0.05)
  laby_95<- data[which(data$sum >= p05),]
  
  # Normalize all columns
  laby_95$factor<- laby_95$sum/mean(laby_95$sum)
  
  laby_norm<- cbind(laby_95[,1:col_start-1],round(laby_95[,(col_start+2):length(laby_95)-2]/laby_95$factor, digits = 0))
  
  return(laby_norm)
}

# Normalize data
laby_order_norm_2<- fun_rem_5(laby_order_2,7)

# Run NMDs
NMDS_laby_2 <- metaMDS(laby_order_norm_2[,-c(1:6)], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# Extract NMDS scores for plotting with ggplot2
## species scores
species.scores_2 <- as.data.frame(scores(NMDS_laby_2, "species"))
species.scores_2$species <- rownames(species.scores_2) 
species.scores_2$type <- "Species"

## site scores
site.scores_2 <- as.data.frame(scores(NMDS_laby_2, "sites"))
site.scores_2 <- cbind (site.scores_2, 
                      climate = laby_order_norm_2$climate, 
                      coalesce_salinity = laby_order_norm_2$coalesce_salinity, 
                      ecosystem = laby_order_norm_2$ecosystem) # Add environmental factors

# Plot NMDS
for (i in c("climate","ecosystem","coalesce_salinity")){
  plot_nmds(species.scores_2,site.scores_2,i,2)
}

# Preparing NDMS data: Relative abundance within class (depth)
laby_order_depth_2 <- laby_asv_wide_surface %>%
  select(file_code,depth_level,latitude,longitude,order,n_reads) %>%
  group_by(file_code,depth_level,latitude,longitude,order)%>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location
  group_by(file_code,depth_level,latitude,longitude,order)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = order, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1:4))))%>% # Add another row that counts the abundance of all
  as.data.frame()

# Normalize data
laby_order_norm_depth_2<- fun_rem_5(laby_order_depth_2,5)

# Run NMDs
NMDS_laby_depth_2 <- metaMDS(laby_order_norm_depth_2[,-c(1:4)], k = 2, trymax = 20, 
                       trace = F, autotransform = FALSE,
                       distance = "bray")

# Extract NMDS scores for plotting with ggplot2
## species scores
species.scores_depth_2 <- as.data.frame(scores(NMDS_laby_depth_2, "species"))
species.scores_depth_2$species <- rownames(species.scores_depth_2) 
species.scores_depth_2$type <- "Species"

## site scores
site.scores_depth_2 <- as.data.frame(scores(NMDS_laby_depth_2, "sites"))
site.scores_depth_2 <- cbind (site.scores_depth_2, 
                        depth_level = laby_order_norm_depth_2$depth_level) # Add environmental factors


# Plot NMDS

plot_nmds(species.scores_depth_2,site.scores_depth_2,"depth_level",2)


##################### CODES BELOW ARE FOR ARCHIVE #####################

# Find the correlation between classes within the kingdom of eukaryotes and laby and find the top 10 R scores

# Preparing data: count total count of each class under each sampling occurrence 

laby_all_euk_fc_wide <- asv_sample_wide %>%
  group_by(file_code,latitude,longitude,class,n_reads) %>% 
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each label 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same class
  spread(key = class, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2,3))))%>% # Add another row that counts the abundance of all 
  as.data.frame()

# Normalize values
laby_all_euk_fc_norm <- fun_rem_5(laby_all_euk_fc_wide,4) 

# Removing 0 values after logging
laby_all_euk_fc_logged <- laby_all_euk_fc_norm
laby_all_euk_fc_logged[, 4:length(laby_all_euk_fc_logged)] <- log(laby_all_euk_fc_norm[4:length(laby_all_euk_fc_logged)])
laby_all_euk_fc_logged[laby_all_euk_fc_logged == -Inf] <- NA

# Loop to obtain r values
cor_res_df <- data.frame() # Create empty list

all_class_name<- data.frame(name= unique(asv_sample_wide$class))%>% 
  filter(!name == "Labyrinthulomycetes")

for (i in all_class_name[,1]){
  cor_res <- cor.test(log(laby_all_euk_fc_norm$Labyrinthulomycetes), log(laby_all_euk_fc_norm[,i]), method = "spearman",exact=FALSE,formula = ~ x + y,alternative = "two.sided")
  cor_res_col<- cbind(class = i, r = round(cor_res$estimate[["rho"]],2),p_value = cor_res$p.value)
  cor_res_df <- rbind(cor_res_df,cor_res_col)
}

# Select top 10 R values
cor_res_df <-cor_res_df%>%
  arrange(desc(r))%>%
  slice(1:10)

# Plotting statistics from a paired correlation
m2_cor_plot_laby_vs_10_euk = list()
for (i in cor_res_df[,1]) {
  plot = ggplot( mapping = aes(x = laby_all_euk_fc_logged$Labyrinthulomycetes, y = laby_all_euk_fc_logged[,i])) +
    geom_point( color = '#440154', size = 1) +
    sm_statCorr(corr_method = 'spearman',
                fit.params = list(linetype = 'dashed'),
                borders = FALSE,
                lm_method = lm)+
    ggtitle(paste("Correlation plot of ",i))+
    labs(y= i, x = "Labyrinthulomycetes")+
    theme(plot.title = element_text(hjust = 0.5))
  m2_cor_plot_laby_vs_10_euk[[i]] <- ggplot_gtable(ggplot_build(plot))
}

png(file=paste0(wd,"/PLOTS/Q1/top_10_cooccurance_plot.jpeg"),width=1536,height=802) 

do.call("grid.arrange", c(m2_cor_plot_laby_vs_10_euk, ncol=3))## display plot

dev.off()

## Part 2: Oceanic vs coastal Laby

# Create file directory for plots from part 2
dir.create(paste0(wd,'/PLOTS/Q1/P2'))

## Part 2.1: Maps comparing coastal and oceanic laby

# Preparing Presence/absence data
laby_PAD_p2.1 <- asv_sample_wide %>%
  group_by(file_code, latitude, longitude,ecosystem, class, n_reads) %>%
  dplyr::summarise(total_count = n() * n_reads, .groups = 'drop') %>% # Find the total count for each location
  select(-n_reads) %>% group_by(file_code, latitude, longitude,ecosystem, class) %>%
  dplyr::summarise(total_count = sum(total_count), .groups = 'drop') %>% # Merging values of the same location
  spread(key = class, value = total_count) %>% # Convert to wide data
  select(file_code, latitude, longitude,ecosystem,Labyrinthulomycetes)%>%
  mutate(P_A = case_when((Labyrinthulomycetes > 0) ~ "present",
                         is.na(Labyrinthulomycetes) ~ "absent"))%>%
  as.data.frame()

# Plotting presence/absence map
png(file=paste0(wd,"/PLOTS/Q1/P2/Laby_ovsc_presence_absence_map.jpeg"),width=1536,height=802) 

base_world +
  geom_jitter(data = laby_PAD_p2.1, aes(x = longitude, y = latitude,color = P_A), size = 3, 
              shape = 16,width = 0.5, height = 0.5)+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#440154","#21918c"))+
  facet_wrap(~ ecosystem)+
  guides(color=guide_legend(title="Present or absent"))

dev.off()


### Preparing laby data: Relative abundance within class (Coastal vs Oceanic)
laby_order_p2.1 <- laby_asv_wide %>%
  group_by(file_code,latitude,longitude,ecosystem,order,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location
  select( -n_reads) %>% group_by(file_code,latitude,longitude,ecosystem,order)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = order, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2,3,4))))%>% # Add another row that counts the abundance of all
  as.data.frame()

laby_order_norm_p2.1<- fun_rem_5(laby_order_p2.1,5)

### Plotting scatterpie map of orders in class Laby
png(file=paste0(wd,"/PLOTS/Q1/P2/Laby_ovsc_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude), 
                  data=laby_order_norm_p2.1,cols=colnames(laby_order_norm_p2.1[,c(5:10)]), 
                  color=NA,alpha=0.8,pie_scale = 0.3)+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(title="Order"))+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~ ecosystem)+
  scale_fill_viridis_d()

dev.off()

## Part 2.2: Treemap to compare coastal and oceanic laby

# Preparing data for treemap
treemap_p2.2 <- laby_order_norm_p2.1 %>%
  select(ecosystem,Amphifilida,Amphitremida,Labyrinthulida,Labyrinthulomycetes_X,Oblongichytrida,Thraustochytrida) %>%
  gather(order, count, Amphifilida:Thraustochytrida, factor_key=TRUE)%>%
  aggregate(count ~ ecosystem+order, ., FUN=sum)%>%
  as.data.frame()

# Plot tree map 
png(file=paste0(wd,"/PLOTS/Q1/P2/Laby_ovsc_treemap.jpeg"),width=1536,height=802) 

ggplot(treemap_p2.2, aes(area = count, fill = order, label = order)) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15) +
  scale_fill_viridis_d()+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(title="Order"))+
  facet_wrap(~ ecosystem)

dev.off()
## Part 3: Zooming into family/ genus level for each environmental variable

# Create file directory for plots from part 1
dir.create(paste0(wd,'/PLOTS/Q1/P3'))

# Part 3.1: Loop to generate bar plots according to the taxa and factors facetwrap(~ecosystem)
dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.1'))
dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.1/Order'))

envi_var<- data.frame(var = c("depth_level","climate","latitude","temperature","coalesce_salinity","region"))

## Function for barplot

barplot <- function(df, variable, pre_taxo_level,pre_taxo_level_n,taxo_level) {
  
  
  # Preparing data for barplot
  df <- df %>% 
    select(file_code, ecosystem ,any_of(c(taxo_level, variable)), n_reads_pct)

  
  # Changing continuous data into discrete data
  if (variable == "temperature") {
    if (length(unique(df$temperature)) > 1) {
      df <- df %>%
        mutate(temperature =  fct_rev(cut_width(temperature, width=5, boundary=0)))
    } else {
      df <- df %>%
        mutate(temperature =  as.factor(temperature))
    }
  }
  
  if (variable == "coalesce_salinity") {
    if (length(unique(df$coalesce_salinity)) > 1) {
      df <- df %>%
        mutate(coalesce_salinity =  fct_rev(cut_width(coalesce_salinity, width=5, boundary=0)))
    } else {
      df <- df %>%
        mutate(coalesce_salinity =  as.factor(coalesce_salinity))
    }
  }
  
  if (variable == "latitude") {  
    if (length(unique(df$latitude)) > 1) {
      df <- df %>%
        mutate(latitude =  fct_rev(cut_width(latitude, width=20, boundary=0)))
    } else {
      df <- df %>%
        mutate(latitude =  as.factor(latitude))
    }
  }
  
  # Obtaining sample number for barplot
  samples <- df %>% 
    select(file_code,ecosystem, any_of(c(variable))) %>% 
    distinct() %>% 
    group_by(ecosystem,across(any_of(variable))) %>% 
    count() 
  
  df <- df %>% 
    group_by(ecosystem,across(any_of(c(taxo_level, variable)))) %>%
    summarize(n_reads_pct = sum(n_reads_pct),.groups = 'drop') %>%
    group_by(ecosystem,across(any_of(variable))) %>% 
    mutate(n_reads_pct = n_reads_pct/sum(n_reads_pct)*100) 
  
  # Plotting barplot
  plot<- ggplot(df) +
    geom_col(aes(y= fct_rev(.data[[variable]]),
                 x=n_reads_pct, 
                 fill=.data[[taxo_level]])) +
    geom_text(data = samples, 
              aes(x = 100, 
                  y = fct_rev(.data[[variable]]), 
                  label = glue::glue("n = {n}"),
              ),
              nudge_x = 10, 
              size = 3) +
    xlab("% of reads") + ylab("") +
    facet_wrap(~ecosystem)+
    scale_fill_viridis_d()+
    ggtitle(paste0("Bar plot of ",variable," in ",pre_taxo_level," ",pre_taxo_level_n)) +
    scale_x_continuous(breaks = seq(0,100,20))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}

## Loop for Class Laby (order level)


for (i in envi_var[,1]){
  
  laby_asv_wide_NAo<- laby_asv_wide %>% 
    drop_na(i)
  
  plot<- barplot(laby_asv_wide_NAo, i, "class","Labyrinthulomycetes","order")
  
  png(file=paste0(wd,"/PLOTS/Q1/P3/P3.1/Barplot_of_Labyrinthulomycetes_",i,".jpeg"),width=1536,height=802) 
  
  print(plot)
  
  dev.off()
}

# Dataframe containing laby orders
laby_order_name<- data.frame(unique(laby_asv_wide$order))

# Replace all "/" with "_
laby_asv_wide <- laby_asv_wide %>% 
  mutate_at(c("family","genus","species"),str_replace_all, "/", "_")

# Loop for all other levels 
for (order_n in laby_order_name[,1]){
  dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.1/Order/',order_n))
  
  df<- laby_asv_wide%>%
    filter(order == order_n)
  
  # Down to family level
  laby_family <- data.frame(unique(df$family))
  dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.1/Order/',order_n,'/family'))
  
  # Plot bar graph for order level 
  
  for (i in envi_var[,1]){
    
    df_NAo<- df %>% 
    drop_na(i)
    
    if(nrow(df_NAo) == 0){
      print("This data frame is empty")
    }else{
      plot<- barplot(df_NAo, i, "order",order_n,"family")
      
      png(file=paste0(wd,"/PLOTS/Q1/P3/P3.1/Order/",order_n,"/Barplot_of_",order_n,"_",i,".jpeg"),width=1536,height=802) 
      
      print(plot)
      
      dev.off()
    }
  }
  
  for (family_n in laby_family[,1]){
    dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.1/Order/',order_n,'/family/',family_n))

    df_2<- df%>%
      filter(family == family_n)

    # Down to genus level
    laby_genus <- data.frame(unique(df_2$genus))
    dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.1/Order/',order_n,'/family/',family_n,'/genus'))

    # Plot bar graph for family level

    for (i in envi_var[,1]){
      
      df_2_NAo<- df_2 %>% 
        drop_na(i)
      
      if(nrow(df_2_NAo) == 0){
        print("This data frame is empty")
      }else{
        plot<- barplot(df_2_NAo, i, "family",family_n,"genus")
        
        png(file=paste0(wd,"/PLOTS/Q1/P3/P3.1/Order/",order_n,"/family/Barplot_of_",family_n,"_",i,".jpeg"),width=1536,height=802)
        
        print(plot)
        
        dev.off()
      }
    }

    # Down to genus level
    for (genus_n in laby_genus[,1]){

      df_3<- df_2%>%
        filter(genus == genus_n)

      # Plot bar graph for order level
      for (i in envi_var[,1]){
        
        df_3_NAo<- df_3 %>% 
          drop_na(i)
        
        if(nrow(df_3_NAo) == 0){
          print("This data frame is empty")
        }else{
          plot<- barplot(df_3_NAo, i, "genus",genus_n,"species")
          
          png(file=paste0(wd,"/PLOTS/Q1/P3/P3.1/Order/",order_n,"/family/",family_n,"/genus/Barplot_of_",genus_n,"_",i,".jpeg"),width=1536,height=802)
          
          print(plot)
          
          dev.off()
        }
      }
    }
  }
}

############## Code below is tentative


# Part 3.2: Plot how genetic diversity (ASV) varies across temperature/latitude 
dir.create(paste0(wd,'/PLOTS/Q1/P3/P3.2'))

# Load necessary files
fasta_data<-read.csv(file = paste0(wd,'/DATA/ALL/asv.csv'))

samples_data<- read.csv(file = paste0(wd,'/DATA/ALL/samples1.csv'))

# Function to create phyloseq file 

make_phyloseq <- function(samples, df, fasta){
  
  # samples = 1000, 192 Mb
  # samples = 2000, 822 Mb
  
  cols_to_include <- c("ecosystem","depth_level","depth","climate","latitude","temperature","coalesce_salinity","region")
  
  # 1. samples table : row names are labeled by file_code
  samples_df <- samples %>%
    select(file_code, any_of(cols_to_include)) %>%
    tibble::column_to_rownames(var = "file_code")
  
  
  # 2. otu table :
  otu <- df %>%
    select(asv_code, file_code, n_reads) %>%
    tidyr::pivot_wider(names_from=file_code,
                       values_from = n_reads,
                       values_fill=list(n_reads=0),
                       values_fn = mean) %>%
    tibble::column_to_rownames(var = "asv_code")
  
  # 3. Taxonomy table
  
  tax <-  fasta %>%
    select(asv_code, kingdom:species) %>%
    distinct(asv_code, .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "asv_code")
  
  ## Create and save to phyloseq object
  
  # Transform into matrixes
  otu_mat <- as.matrix(otu)
  tax_mat <- as.matrix(tax)
  
  # Transform to phyloseq object and save to Rdata file
  OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = phyloseq::tax_table(tax_mat)
  samples = phyloseq::sample_data(samples_df)
  
  
  cat("Make phyloseq done \n")
  
  ps <- phyloseq::phyloseq(OTU, TAX, samples)
  
}

laby_phyloseq<- make_phyloseq(samples_data,laby_asv_wide,fasta_data)
  
# Discretize the data (must make sure that there is more than one value)

sample_data(laby_phyloseq)[["latitude"]] =  cut_width(phyloseq::get_variable(laby_phyloseq,c("latitude")), width=20, boundary=0) # For latitude

phyloseq::sample_data(laby_phyloseq)[["temperature"]] =  fct_rev(cut_width(phyloseq::get_variable(laby_phyloseq,c("temperature")), width=5, boundary=0)) # For temperature

phyloseq::sample_data(laby_phyloseq)[["coalesce_salinity"]] =  fct_rev(cut_width(phyloseq::get_variable(laby_phyloseq,c("coalesce_salinity")), width=5, boundary=0)) # For salinity

phyloseq::sample_data(laby_phyloseq)[["depth"]] =  fct_rev(cut_width(phyloseq::get_variable(laby_phyloseq,c("depth")), width=1000, boundary=0)) # For depth

measures = data.frame(c("Shannon","Chao1","Simpson"))

for (i in measures[,1]){
  diversity <- phyloseq::estimate_richness(laby_phyloseq, split = TRUE, measures = i) %>%
    tibble::rownames_to_column(var = "file_code")
  
  
  samples <- data.frame(phyloseq::sample_data(laby_phyloseq)) 
  
  samples <- samples %>% 
    tibble::rownames_to_column(var = "file_code") %>% 
    left_join(diversity) %>% 
    pivot_longer(cols = i, values_to = "diversity", names_to = "measures")%>%
    drop_na()
  
  png(file=paste0(wd,"/PLOTS/Q1/P3/P3.2/Plot_",i,"_temp_lat_trend.jpeg"),width=1536,height=802)
  
  plot<- ggplot(data = samples, aes(x= latitude, y = diversity)) +
    scale_y_continuous(limits = c(0.01, NA)) +
    geom_smooth(method = "gam") +
    geom_jitter(aes(color = temperature),size=2, alpha=0.85 )+
    facet_grid(rows = vars(ecosystem))+
    ggtitle(paste0("Plotting ",i," diversity trend across latitude and temperature")) +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_viridis_d()+
    coord_flip()
  
  print(plot)
  
  dev.off()
}

