## Setup
###Load packages & plot setup

# NOTE PLEASE CHANGE DIRECTORY ACCORDINGLY
# BE AWARE OF WHERE THE CODE AND THE DATA IS STORED

wd<- setwd("C:/Users/chiny/Desktop/FYP/R") 

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
asv_sample_wide<-read.table(file = paste0(wd,'/DATA/ALL/metapr2_ASVs_selected_abundance_Eukaryota.tsv'), sep = '\t', header = TRUE)

#Preparing data for RDA

# Filter for class laby 
laby_asv_wide <- filter(asv_sample_wide,class == "Labyrinthulomycetes")


### Preparing laby data: abundance within class 
laby_order_wide <- laby_asv_wide %>%
  group_by(label,latitude,date,depth_level,substrate,climate,temperature,salinity,ecosystem,order,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each label
  select( -n_reads) %>% group_by(label,latitude,date,depth_level,substrate,climate,temperature,salinity,ecosystem,order)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same order
  spread(key = order, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  as.data.frame()

## Filling in the gap for climate variable
## Noted: I will be dropping the variable ecosystem_climate as it is extremely similar to climate
laby_edited_wide <- mutate(laby_order_wide, climate = case_when((latitude <= 90 & latitude > 66.5) ~ "polar",
                                                               (latitude >= -90 & latitude < -66.5) ~ "polar",
                                                               (latitude <= 66.5 & latitude > 23.5) ~ "temperate",
                                                               (latitude >= -66.5 & latitude < -23.5) ~ "temperate",
                                                               (latitude <= 23.5 & latitude >= -23.5) ~ "equatorial")) 


## Cut continuous variable into discrete variable
laby_edited_wide$temperature<- as.factor(cut(laby_edited_wide$temperature,breaks = seq(min(laby_edited_wide$temperature, na.rm = TRUE),
                                                                                     max(laby_edited_wide$temperature, na.rm = TRUE),
                                                                                     by = 5)))
laby_edited_wide$salinity<- as.factor(cut(laby_edited_wide$salinity,breaks = seq(min(laby_edited_wide$salinity, na.rm = TRUE),
                                                                               max(laby_edited_wide$salinity, na.rm = TRUE),
                                                                               by = 5)))

laby_edited_wide$latitude<- as.factor(cut(laby_edited_wide$latitude,breaks = seq(min(-90, na.rm = TRUE), 
                                                                               max(90, na.rm = TRUE), 
                                                                               by = 30)))

## Converting dates into months
laby_edited_wide$date<- months(as.Date(laby_edited_wide$date))

# Make each collection site unique according to their labels 
laby_edited_wide$label <- ave(as.character(laby_edited_wide$label), laby_edited_wide$label, FUN=function(x) if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])

# Ensure that they are all factors
laby_edited_wide<- laby_edited_wide%>%
  mutate_if(is.character,as.factor)

#Remove rows with NA
laby_edited_wide<-na.omit(laby_edited_wide)

# Check for collinearity in environmental factors

# Change categorical data into numerical data
laby_edited_unclass <- data.frame(sapply(laby_edited_wide[,-1], unclass)) 

# creating correlation matrix
corr_mat <- round(abs(cor(laby_edited_unclass[,1:8])),2)
corr_mat <- melt(corr_mat)

# plotting the correlation heatmap
png(file=paste0(wd,"/PLOTS/Q2/Correlation_heatmap_of_environmental_variables.jpeg"),width=1536,height=802) 

ggplot(data = corr_mat, aes(x=Var1, y=Var2,
                                   fill=value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black", size = 10)+
  scale_fill_distiller(palette = "YlOrRd", direction=1) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))

dev.off()

# Redundancy analysis in class Labyrinthulomycetes

# Subset environmental factors 
laby_edited_wide.envi<- laby_edited_wide[,1:9]

# Subsetting population data
laby_edited_wide.pop<- laby_edited_wide[,-(2:9)]

# Since there are alot of zeros seen there maybe a double zero problem, hence I will be using the Hellinger transformation which expresses abundances as the square-root of their relative abundance at each site 
laby_hel_wide.pop <- decostand(laby_edited_wide.pop[,-1], method = "hellinger")

laby.rda <- rda(laby_hel_wide.pop ~ ., data = laby_edited_wide.envi[,-1])

# Selecting variables 

fwd.sel <- ordiR2step(rda(laby_hel_wide.pop ~ 1, data = laby_edited_wide.envi[,-1]), # lower model limit 
                      laby.rda, # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE)

# Check the new model with forward-selected variables
fwd.sel$call

# Test again with new model
laby.rda.signif <- rda(laby_hel_wide.pop ~latitude + date + salinity + temperature + ecosystem + depth_level + climate, data = laby_edited_wide.envi[,-1])

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

dimcheckMDS(laby_edited_wide.pop[,-1], distance = "bray",
            k = 6, trymax = 20,
            autotransform = FALSE) 

dev.off()

# Run NMDs
NMDS_laby <- metaMDS(laby_edited_wide.pop[,-1], k = 2, trymax = 20, 
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
                      latitude = laby_edited_wide.envi[,-1]$latitude, 
                      date = laby_edited_wide.envi[,-1]$date, 
                      depth_level = laby_edited_wide.envi[,-1]$depth_level,
                      substrate = laby_edited_wide.envi[,-1]$substrate, 
                      climate = laby_edited_wide.envi[,-1]$climate, 
                      temperature = laby_edited_wide.envi[,-1]$temperature,
                      salinity = laby_edited_wide.envi[,-1]$salinity, 
                      ecosystem = laby_edited_wide.envi[,-1]$ecosystem) # Add environmental factors

# Plot NMDS using ggplot2
envi_var<- data.frame(var = c("depth_level","substrate","climate","ecosystem","temperature","salinity","latitude","date"))

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
  theme_minimal() 
  
  # Saving plots
  png(file=paste0(wd,"/PLOTS/Q2/NMDS_plot_laby_",i,".jpeg"),width=1536,height=802) 
  print(plot)
  dev.off()

}
