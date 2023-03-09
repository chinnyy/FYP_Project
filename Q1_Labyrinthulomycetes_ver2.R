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

## Part 0: Removing unwanted data

# Remove the non-protist
asv_sample_wide<-asv_sample_wide[!(asv_sample_wide$division == "Metazoa" | asv_sample_wide$division == "Fungi"),]


## Part 1: Distribution of laby

# Create file directory for plots from part 1.1
dir.create(paste0(wd,'/PLOTS/Q1'))

## Part 1.1: Presence/absence map

# Create file directory for plots from part 1.1
dir.create(paste0(wd,'/PLOTS/Q1/P1.1'))

# Preparing Presence/absence data
laby_PAD <- asv_sample_wide %>%
  group_by(file_code, latitude, longitude, class, n_reads) %>%
  dplyr::summarise(total_count = n() * n_reads, .groups = 'drop') %>% # Find the total count for each location
  select(-n_reads) %>% group_by(file_code, latitude, longitude, class) %>%
  dplyr::summarise(total_count = sum(total_count), .groups = 'drop') %>% # Merging values of the same location
  spread(key = class, value = total_count) %>% # Convert to wide data
  select(file_code, latitude, longitude,Labyrinthulomycetes)%>%
  mutate(P_A = case_when((Labyrinthulomycetes > 0) ~ "present",
                         is.na(Labyrinthulomycetes) ~ "absent"))%>%
  as.data.frame()

# Plotting presence/absence map
png(file=paste0(wd,"/PLOTS/Q1/P1.1/Laby_presence_absence_map.jpeg"),width=1536,height=802) 

base_world +
  geom_jitter(data = laby_PAD, aes(x = longitude, y = latitude,color = P_A), size = 2, 
              shape = 16)+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#440154","#21918c"))+
  guides(color=guide_legend(title="Present or absent"))

dev.off()

## Part 1.2: Scatterpie map of laby (orders)

# Create file directory for plots from part 1.2
dir.create(paste0(wd,'/PLOTS/Q1/P1.2'))

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
fun_rem_5<- function(data,col_start) {
  p05 <- quantile(data$sum, 0.05)
  laby_95<- data[which(data$sum >= p05),]
  
  # Normalize all columns
  laby_95$factor<- laby_95$sum/mean(laby_95$sum)
  
  laby_norm<- cbind(laby_95[,1:col_start-1],round(laby_95[,(col_start+2):length(laby_95)-2]/laby_95$factor, digits = 0))
  
  return(laby_norm)
}

laby_order_norm<- fun_rem_5(laby_order,4)

### Plotting scatterpie map of orders in class Laby
png(file=paste0(wd,"/PLOTS/Q1/P1.2/Laby_global_abun_scatterpie_map.jpeg"),width=1536,height=802) 

base_world+
  geom_scatterpie(aes(x=longitude, y=latitude), 
                  data=laby_order_norm,cols=colnames(laby_order_norm[,c(4:9)]), 
                  color=NA,alpha=0.8,pie_scale = 0.3)+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(title="Order"))+
  scale_fill_viridis_d()

dev.off()


## Part 1.3: Map of the laby contribution metabarcodes and dominant laby group [TENTATIVE]

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

