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

## Load dataset
wd<- setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory

asv_sample_wide<-read.csv(file = paste0(wd,'/DATA/ALL/woa_metapr2_ASVs_selected_abundance_Eukaryota.csv'))

## Part 0: Removing unwanted data
asv_sample_wide<-asv_sample_wide %>%
  select(-c("NCBI_SRA_accession","NCBI_SRA_biosample","sample_code","substrate","DNA_RNA","fraction_name","sample_fixation","project","elevation","biome","vegetation","ice_coverage","dataset_code","dataset_name","gene","gene_region","organelle","ecosystem_climate","substrate_type","label"))%>%
  filter(depth_level=='euphotic'|depth_level=='mesopelagic'|depth_level=='surface')%>%
  as.data.frame()


# Part 1: Global abundance map Laby vs other protists

# Filter for class laby 
laby_asv_wide <- filter(asv_sample_wide,class == "Labyrinthulomycetes")

# Preparing laby data: Relative abundance within class
laby_order <- laby_asv_wide %>%
  group_by(file_code,latitude,longitude,order,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,order)%>%
   dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
   spread(key = order, value = total_count)%>% # Convert to wide data
   replace(is.na(.), 0) %>%
   mutate(sum = rowSums(select(., -c(1,2,3))))%>% # Add another row that counts the abundance of all 
  as.data.frame()

# Function to remove bottom 5% of the dataset and normalize all columns
fun_rem_5<- function(data) {
  p05 <- quantile(data$sum, 0.05)
  laby_95<- data[which(data$sum >= p05),]
  
  # Normalize all columns
  laby_95$factor<- laby_95$sum/mean(laby_95$sum)
  
  laby_norm<- cbind(laby_95[1:3],round(laby_95[,6:length(laby_95)-2]/laby_95$factor, digits = 0))


  return(laby_norm)
}

laby_order_norm<- fun_rem_5(laby_order)

### Plotting scatterpie map of orders in class Laby
png(file=paste0(wd,"/PLOTS/Q1/Laby_global_abun_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude), 
                  data=laby_order_norm,cols=colnames(laby_order_norm[,c(4:9)]), 
                  color=NA,alpha=0.8,pie_scale = 0.75)+
  theme(legend.position = "bottom")+
  scale_fill_viridis_d()

dev.off()

# Preparing laby vs protist data: abundance with other class (Kingdom level)
laby_euk<- asv_sample_wide %>%
  group_by(file_code,latitude,longitude,class,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  mutate(class = replace(class, class != "Labyrinthulomycetes", "Not Labyrinthulomycetes"))%>%
  group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = class, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>% 
  mutate(sum = rowSums(select(., -c(1,2,3))))%>% # Add another row that counts the abundance of all 
  as.data.frame()

# Normalize data
laby_euk_norm<- fun_rem_5(laby_euk) 

# Preparing laby vs protist data: abundance with other class (Supergroup level)
laby_stram<- asv_sample_wide %>%
  filter(supergroup == "Stramenopiles")%>%
 group_by(file_code,latitude,longitude,class,n_reads) %>%
 dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  mutate(class = replace(class, class != "Labyrinthulomycetes", "Not Labyrinthulomycetes"))%>%
  group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = class, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2,3))))%>% # Add another row that counts the abundance of all
  as.data.frame()

# Normalize data
laby_stram_norm<- fun_rem_5(laby_stram) 

# Preparing laby vs protist data: abundance with other class (Division level)
laby_sagen<- asv_sample_wide %>%
  filter(division == "Sagenista")%>%
  group_by(file_code,latitude,longitude,class,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  mutate(class = replace(class, class != "Labyrinthulomycetes", "Not Labyrinthulomycetes"))%>%
  group_by(file_code,latitude,longitude,class)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = class, value = total_count)%>% # Convert to wide data
  replace(is.na(.), 0) %>%
  mutate(sum = rowSums(select(., -c(1,2,3))))%>% # Add another row that counts the abundance of all
  as.data.frame()

# Normalize data
laby_sagen_norm<- fun_rem_5(laby_sagen) 

# Plotting scatterpie map of laby vs protist data: abundance with other class (Kingdom level)
png(file=paste0(wd,"/PLOTS/Q1/Laby_vs_eukaryotes_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude), 
                  data=laby_euk_norm,cols=colnames(laby_euk_norm[,c(4:5)]), 
                  color=NA,alpha=0.8,pie_scale = 0.75)+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("#440154","#21918c"))

dev.off()

# Plotting scatterpie map of laby vs protist data: abundance with other class (Supergroup level)
png(file=paste0(wd,"/PLOTS/Q1/Laby_vs_stramenopiles_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude), 
                  data=laby_stram_norm,cols=colnames(laby_stram_norm[,c(4:5)]), 
                  color=NA,alpha=0.8,pie_scale = 0.75)+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("#440154","#21918c"))

dev.off()

# Plotting scatterpie map of laby vs protist data: abundance with other class (Division level)
png(file=paste0(wd,"/PLOTS/Q1/Laby_vs_sagenista_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude), 
                  data=laby_sagen_norm,cols=colnames(laby_sagen_norm[,c(4:5)]), 
                  color=NA,alpha=0.8,pie_scale = 0.75)+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("#440154","#21918c"))

dev.off()


# Part 2: Co-occurrence plots of Laby vs other protist classes

# METHOD 2: Find the correlation between classes within the kingdom of eukaryotes and laby and find the top 10 R scores

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
laby_all_euk_fc_norm <- fun_rem_5(laby_all_euk_fc_wide) 

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
  plot = ggplot( mapping = aes(x = log(laby_all_euk_fc_norm$Labyrinthulomycetes), y = log(laby_all_euk_fc_norm[,i]))) +
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

# Save plots to .jpeg. Makes a separate file for each plot.
for (i in 1:10) {
  png(file=paste0(wd,"/PLOTS/Q1/Method_2_Laby_vs_",cor_res_df[i,1],"_corr_plot.jpeg"),width=1536,height=802) 
  gridExtra::grid.arrange(m2_cor_plot_laby_vs_10_euk[[i]])
  dev.off()
}

################# All parts above are normalized and running #####################

# Part 3: Distribution of order within laby according to environmental variables

# Preparing laby data with envi factors: Relative abundance within class
laby_order_wide <- laby_asv_wide %>%
  group_by(file_code,latitude,longitude,dataset_id,date,depth_level,climate,temperature,salinity,ecosystem,ss_temperature,coalesce_salinity,order,n_reads) %>%
  dplyr::summarise(total_count=n()*n_reads,.groups = 'drop')%>% # Find the total count for each location 
  select( -n_reads) %>% group_by(file_code,latitude,longitude,dataset_id,date,depth_level,climate,temperature,salinity,ecosystem,ss_temperature,coalesce_salinity,order)%>%
  dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same location
  spread(key = order, value = total_count)%>% # Convert to wide data
  mutate_at(c(13:18), ~replace_na(.,0))%>%
  mutate(sum = rowSums(select(., -c(1:12))))%>% # Add another row that counts the abundance of all
  as.data.frame()

## Cut continuous variable into discrete variable
laby_order_wide$temperature<- as.factor(cut(laby_order_wide$temperature,breaks = seq(min(laby_order_wide$temperature, na.rm = TRUE), 
                                                                                     max(laby_order_wide$temperature, na.rm = TRUE), 
                                                                                     by = 5)))

laby_order_wide$ss_temperature<- as.factor(cut(laby_order_wide$ss_temperature,breaks = seq(min(laby_order_wide$ss_temperature, na.rm = TRUE), 
                                                                                     max(laby_order_wide$ss_temperature, na.rm = TRUE), 
                                                                                     by = 5)))

laby_order_wide$salinity<- as.factor(cut(laby_order_wide$salinity,breaks = seq(min(laby_order_wide$salinity, na.rm = TRUE), 
                                                                               max(laby_order_wide$salinity, na.rm = TRUE), 
                                                                               by = 2)))

laby_order_wide$coalesce_salinity<- as.factor(cut(laby_order_wide$coalesce_salinity,breaks = seq(min(laby_order_wide$coalesce_salinity, na.rm = TRUE), 
                                                                               max(laby_order_wide$coalesce_salinity, na.rm = TRUE), 
                                                                               by = 2)))

laby_order_wide$latitude<- as.factor(cut(laby_order_wide$latitude,breaks = seq(min(-90, na.rm = TRUE), 
                                                                               max(90, na.rm = TRUE), 
                                                                               by = 30)))
## Converting dates into months
laby_order_wide$date<- months(as.Date(laby_order_wide$date))
laby_order_wide$date = factor(laby_order_wide$date, levels=month.name)

# Normalize data (i realise my previous function is too hard coded, will make a new version for this)
p05 <- quantile(laby_order_wide$sum, 0.05)
laby_95<- laby_order_wide[which(laby_order_wide$sum >= p05),]

# Normalize all columns
laby_95$factor<- laby_95$sum/mean(laby_95$sum)

laby_norm<- cbind(laby_95[1:12],round(laby_95[,15:length(laby_95)-2]/laby_95$factor, digits = 0)) # Indexing is weird need to +2

# Plotting Treemap of environmental factors affecting abundance in class Labyrinthulomycetes 

envi_var<- data.frame(var = c("depth_level","climate","ecosystem","latitude","temperature","ss_temperature","salinity","coalesce_salinity"))

for (i in envi_var[,1]){
  # Preparing data
  laby_order_envi <- laby_norm %>%
    gather(.,order, total_count,13:18 , factor_key=TRUE)%>% # Converting to long data
    group_by(across(all_of(c(i,"order","total_count"))))%>%
    dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>% #Merging all data according to ecosystem
    group_by(across(all_of(c(i,"order")))) %>%
    dplyr::summarise(total_count = sum(total_count),.groups = 'drop')%>%# Merging values of the same order
    as.data.frame()

  
  # Treemap of laby
  plot<- ggplot(laby_order_envi, aes(fill = order, x= get(i),y=total_count)) +
    geom_bar(position="fill", stat="identity")+
    scale_fill_viridis_d()+
    ylab("% of reads")+
    scale_y_continuous(labels = percent)+
    ggtitle(paste0("Bar plot of Labyrinthulomycetes distribution according to ",i)) +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5,size=22))+
    coord_flip()
  
  # Saving plots
  png(file=paste0(wd,"/PLOTS/Q1/Laby_envi_factor_",i,"_treemap.jpeg"),width=1536,height=802) 
  print(plot)
  dev.off()
}


############## ALL CODES ABOVE RUN SMOOTHLY ############### MAY CONSIDER SHORTENING THE PREPING DATA AND PLOTTING PART

########### CODE STORAGE SPACE ######################
#USELESS CODES I DONT WANT TO DELETE


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
  png(file=paste0(wd,"/PLOTS/Laby_",plot_name_list[i],"_scatter_map.jpeg"),width=1536,height=802) 
  print(plot_list[[i]])
  dev.off()
}


