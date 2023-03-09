## Setup
###Load packages & plot setup

# NOTE PLEASE CHANGE DIRECTORY ACCORDINGLY
# BE AWARE OF WHERE THE CODE AND THE DATA IS STORED


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
library(reshape2)


## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

asv_sample_wide<-read.csv(file = paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'))


## Part 0: Removing unwanted data
asv_sample_wide<-asv_sample_wide %>%
  select(-c("NCBI_SRA_accession","NCBI_SRA_biosample","sample_code","dataset_id","substrate","DNA_RNA","fraction_name","sample_fixation","project","elevation","biome","vegetation","ice_coverage","dataset_code","dataset_name","gene","gene_region","organelle","ecosystem_climate","substrate_type","label"))%>%
  filter(depth_level=='euphotic'|depth_level=='mesopelagic'|depth_level=='surface')%>%
  as.data.frame()


#Preparing data for RDA

# Filter for class laby 
laby_asv_wide <- filter(asv_sample_wide,class == "Labyrinthulomycetes")

# Preparing laby data with envi factors: Relative abundance within class
laby_order_wide <- laby_asv_wide %>%
  group_by(file_code,latitude,longitude,date,depth_level,climate,temperature,salinity,ecosystem,ss_temperature,coalesce_salinity,order,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,date,depth_level,climate,temperature,salinity,ecosystem,ss_temperature,coalesce_salinity,order)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = order, value = total_count)%>% # Convert to wide data
  mutate_at(c(12:17), ~replace_na(.,0))%>%
  mutate(sum = rowSums(select(., -c(1:11))))%>% # Add another row that counts the abundance of all
  as.data.frame()

laby_order_wide$date<- months(as.Date(laby_order_wide$date))
laby_order_wide$date = factor(laby_order_wide$date, levels=month.name)


# Make each collection site unique according to their labels 
laby_order_wide$file_code <- ave(as.character(laby_order_wide$file_code), laby_order_wide$file_code, FUN=function(x) if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])

# Ensure that they are all factors
laby_order_wide<- laby_order_wide%>%
  mutate_if(is.character,as.factor)

#Remove rows with NA
laby_order_wide<-na.omit(laby_order_wide)

# Check for collinearity in environmental factors

# Change categorical data into numerical data 
laby_edited_unclass <- cbind(laby_order_wide[,-c(4:6,9,12:18)],data.frame(sapply(laby_order_wide[,c(4:6,9)], unclass)) )

# creating correlation matrix
corr_mat <- round(abs(cor(laby_edited_unclass[,2:11])),2)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr_mat){
  corr_mat[lower.tri(corr_mat)]<- NA
  return(corr_mat)
}

upper_tri <- get_upper_tri(corr_mat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

# plotting the correlation heatmap
png(file=paste0(wd,"/PLOTS/Q2/Correlation_heatmap_of_environmental_variables.jpeg"),width=1536,height=802) 

ggplot(data = melted_cormat, aes(x=Var2, y=Var1,
                                   fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black", size = 10)+
  scale_fill_distiller(palette = "YlOrRd", direction=1) +
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))

dev.off()

# Remove ss_temperature,salinity and longitude

# Change categorical data into numerical data 
laby_edited_unclass_2 <- cbind(laby_order_wide[,-c(3:6,8:10,12:18)],data.frame(sapply(laby_order_wide[,c(4:6,9)], unclass)) )

# creating correlation matrix
corr_mat_2 <- round(abs(cor(laby_edited_unclass_2[,2:8])),2)

upper_tri_2 <- get_upper_tri(corr_mat_2)

melted_cormat_2 <- melt(upper_tri_2, na.rm = TRUE)

# plotting the correlation heatmap
png(file=paste0(wd,"/PLOTS/Q2/Correlation_heatmap_of_environmental_variables_reduce.jpeg"),width=1536,height=802) 

ggplot(data = melted_cormat_2, aes(x=Var2, y=Var1,
                                 fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black", size = 10)+
  scale_fill_distiller(palette = "YlOrRd", direction=1) +
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))

dev.off()

# Redundancy analysis in class Labyrinthulomycetes

# Remove ss_temperature, salinity and longitude
laby_reduce_wide<- laby_order_wide[,-c(3,8,10)]

# Normalize data (i realise my previous function is too hard coded, will make a new version for this)
p05 <- quantile(laby_reduce_wide$sum, 0.05)
laby_95<- laby_reduce_wide[which(laby_reduce_wide$sum >= p05),]

# Normalize all columns
laby_95$factor<- laby_95$sum/mean(laby_95$sum)

laby_norm<- cbind(laby_95[1:8],round(laby_95[,11:length(laby_95)-2]/laby_95$factor, digits = 0)) # Indexing is weird need to +2

# Subset environmental factors 
laby_norm.envi<- laby_norm[,1:8]

# Subsetting population data
laby_norm.pop<- laby_norm[,-(2:8)]

# Since there are alot of zeros seen there maybe a double zero problem, hence I will be using the Hellinger transformation which expresses abundances as the square-root of their relative abundance at each site 
laby_hel_wide.pop <- decostand(laby_norm.pop[,-1], method = "hellinger")

laby.rda <- rda(laby_hel_wide.pop ~ ., data = laby_norm.envi[,-1])

# Selecting variables 

fwd.sel <- ordiR2step(rda(laby_hel_wide.pop ~ 1, data = laby_norm.envi[,-1]), # lower model limit 
                      laby.rda, # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE)

# Check the new model with forward-selected variables
fwd.sel$call

# Test again with new model
laby.rda.signif <- rda(formula = laby_hel_wide.pop ~ ecosystem + temperature + climate + 
                         date + coalesce_salinity + latitude, data = laby_norm.envi[,-1])

# Check the adjusted r2
RsquareAdj(laby.rda.signif)

# Test significance of RDA
anova.cca(laby.rda.signif, step = 1000, by = "term")

# RDA plot

# Type 1 scaling: Distances among objects reflect their similarities
png(file=paste0(wd,"/PLOTS/Q2/RDA_plot_type_1_scaling.jpeg"),width=1536,height=802) 

ordiplot(laby.rda.signif, scaling = 1, type = "text")

dev.off()
# Type 2 scaling: Angles between variables reflect their correlation
png(file=paste0(wd,"/PLOTS/Q2/RDA_plot_type_2_scaling.jpeg"),width=1536,height=802) 

ordiplot(laby.rda.signif, scaling = 2, type = "text")

dev.off()

# To test for collinearity in data 
sqrt(vif.cca(laby.rda))

# Non-metric Multidimensional Scaling (NMDS) in class Labyrinthulomycetes

# Run scree plot
png(file=paste0(wd,"/PLOTS/Q2/Scree_plot_NMDs.jpeg"),width=781,height=762) 

dimcheckMDS(laby_norm.pop[,-1], distance = "bray",
            k = 6, trymax = 20,
            autotransform = FALSE) 

dev.off()

# Run NMDs
NMDS_laby <- metaMDS(laby_norm.pop[,-1], k = 2, trymax = 20, 
                     trace = F, autotransform = FALSE,
                     distance = "bray")

# View NMDS results
NMDS_laby

# Plot Stress plot
png(file=paste0(wd,"/PLOTS/Q2/Stress_plot_NMDs.jpeg"),width=781,height=762) 

stressplot(NMDS_laby)

dev.off()

# Plot NMDs
# Extract NMDS scores for plotting with ggplot2
## species scores
species.scores <- as.data.frame(scores(NMDS_laby, "species"))
species.scores$species <- rownames(species.scores) 
species.scores$type <- "Species"

## site scores
site.scores <- as.data.frame(scores(NMDS_laby, "sites"))
site.scores <- cbind (site.scores, 
                      latitude = laby_norm.envi[,-1]$latitude, 
                      date = laby_norm.envi[,-1]$date, 
                      depth_level = laby_norm.envi[,-1]$depth_level,
                      climate = laby_norm.envi[,-1]$climate, 
                      temperature = laby_norm.envi[,-1]$temperature,
                      salinity = laby_norm.envi[,-1]$coalesce_salinity, 
                      ecosystem = laby_norm.envi[,-1]$ecosystem) # Add environmental factors

## Cut continuous variable into discrete variable
site.scores$temperature<- as.factor(cut(site.scores$temperature,breaks = seq(min(site.scores$temperature, na.rm = TRUE), 
                                                                                     max(site.scores$temperature, na.rm = TRUE), 
                                                                                     by = 5)))


site.scores$salinity<- as.factor(cut(site.scores$salinity,breaks = seq(min(site.scores$salinity, na.rm = TRUE), 
                                                                               max(site.scores$salinity, na.rm = TRUE), 
                                                                               by = 2)))

site.scores$latitude<- as.factor(cut(site.scores$latitude,breaks = seq(min(site.scores$latitude, na.rm = TRUE), 
                                                                                                 max(site.scores$latitude, na.rm = TRUE), 
                                                                                                 by = 30)))


# Plot NMDS using ggplot2
envi_var<- data.frame(var = c("depth_level","climate","ecosystem","temperature","salinity","latitude","date"))

for (i in envi_var[,1]){
  plot<-ggplot() + 
  geom_point(data = species.scores,
             aes(x = NMDS1,
                 y = NMDS2),
             size = 2) +
  geom_point(data = site.scores,
             aes(x = NMDS1,
                 y = NMDS2,
                 colour = get(i)),
             shape = 17,
             size = 2) +
  # classify sites by variables and display using 95% confidence ellipses
  stat_ellipse(data = site.scores,
               aes(x = NMDS1,
                   y = NMDS2,
                   colour = get(i)),
               type = "norm") +
  labs(color=i) +
  scale_color_viridis_d()+
  ggtitle(paste0("NMDS plot of laby and ",i)) +
  theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,size=22))
  
  # Saving plots
  png(file=paste0(wd,"/PLOTS/Q2/NMDS_plot_laby_",i,".jpeg"),width=1536,height=802) 
  print(plot)
  dev.off()

}
