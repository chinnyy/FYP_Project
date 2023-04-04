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
  
  NMDS1_quantiles <- quantile(data.scores_norm$NMDS1, c(0.10, 0.90)) 
  NMDS2_quantiles <- quantile(data.scores_norm$NMDS2, c(0.10, 0.90))
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
