# Set Working Directory ---------------------------------------------------

setwd("")
getwd()

# Set color palette ----------------------------------------------

CINARpalette=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288")

# Load required packages --------------------------------------------------

library(tidyverse)
library(ggpubr)
library(data.table)
library(vegan)
library(ggcorrplot)
library(RColorBrewer)
library(lubridate)
library(ggforce)
library(Hmisc)
library(reshape2)
library(readxl)
library(ggthemes)
library(viridis)
library(ggdendro)
library(grid)
library(FactoMineR)
library(factoextra)
library(rcartocolor)
library(corrplot)
library(stats)
library(broom)
library(stringr)
library(ggpubr)
library(cowplot)
library(sf)
library(here)
library(ggspatial)
library(ggrepel)
library(rcartocolor)
library(UpSetR)


# Load data ---------------------------------------------------------------

data <- read_csv("mtab_files/CINARhabitat_Quant_table.csv")
data <- as_tibble(data,rownames="rowID_full")
data$rep_group<-as.factor(data$rep_group)
seq <- read_csv("mtab_files/CINAThabitat_Seq_Sheet.csv")
limits <- read_csv("mtab_files/CINARhabitat_LODLOQ_limits.csv")
benthic <- read.delim("mtab_files/CINARhabitat_Benthic_cover.txt", sep = "\t", header = TRUE)

# Merge in file name ------------------------------------------------------

seq_pos <- seq%>%
  filter(ionMode == "pos", goodData == 1, sType == "rep" )%>%
  select(File_Name,sMatlabName)%>%
  rename(File_Name_Pos = File_Name)

seq_neg <- seq%>%
  filter(ionMode == "neg", goodData == 1, sType == "rep" )%>%
  select(File_Name,sMatlabName)%>%
  rename(File_Name_Neg = File_Name)

seq_mode_sidebyside <- merge(seq_neg,seq_pos,by = "sMatlabName")%>%
  rename(adaptedDervLabel = sMatlabName)

data_wFileName <- merge(seq_mode_sidebyside,data,by = "adaptedDervLabel",all.y = TRUE)

rm(seq_pos,seq_neg, seq_mode_sidebyside,seq, data)

# Subset metabolite data -------------------------------------------------------------

head(data_wFileName)

habitat <- data_wFileName %>%
  filter(sample_type=="peak_photo" & Date <= "1/23/2021" & Site != "LB_seagrass")%>%
  rename(NPOC = NPOC_uM_, TN = TN_uM_)%>%
  mutate(Site = factor(Site, levels = c("Ditliff", "Coco_loba", "Joel's Shoal", "Yawzi", "Tektite")))%>%
  mutate(Site_abv =
           case_when(Site =="Ditliff" ~ "DL", 
                     Site =="Coco_loba" ~ "CO",
                     Site =="Joel's Shoal" ~ "JS",
                     Site =="Yawzi" ~ "YZ",
                     Site =="Tektite" ~ "TK"),.after = Site)%>%
  mutate(Site_abv = factor(Site_abv, levels = c("DL", "CO", "JS", "YZ", "TK")))%>%
  mutate(grouping =
           case_when(Site =="Ditliff" ~ "high_coral", 
                     Site =="Coco_loba" ~ "low_coral",
                     Site =="Joel's Shoal" ~ "high_coral",
                     Site =="Yawzi" ~ "low_coral",
                     Site =="Tektite" ~ "high_coral"),.after = Site_abv)%>%
  mutate(grouping = factor(grouping, levels = c("high_coral", "low_coral")))%>%
  mutate(bay_group =
           case_when(Site =="Ditliff" ~ "Fish Bay", 
                     Site =="Coco_loba" ~ "Fish Bay",
                     Site =="Joel's Shoal" ~ "Fish Bay",
                     Site =="Yawzi" ~ "Lameshur Bay",
                     Site =="Tektite" ~ "Lameshur Bay"),.after = grouping)%>%
  mutate(bay_group = factor(bay_group, levels = c("Fish Bay", "Lameshur Bay")))%>%
  arrange(desc(Site_abv))


str(habitat)

habitat.mat <- habitat
rownames(habitat.mat) <- habitat$adaptedDervLabel
habitat.mat <- habitat.mat[,19:ncol(habitat.mat)]
habitat.mat <- habitat.mat[, ! apply(habitat.mat , 2 , function(x) 
  sd(x, na.rm = TRUE)==0 ) ] #remove metabolites with a SD of zero 

x <- colnames(habitat.mat) #pull out metabolite names

habitat_long <- habitat %>%
  pivot_longer(cols=19:ncol(habitat),names_to="metabolites",values_to="concentration",values_drop_na=TRUE)%>%
  filter(metabolites %in% x) #create long format with only metabolites SD > 0

x <- unique(habitat_long$metabolites)

habitat_wide <- habitat_long %>%
  pivot_wider(names_from = metabolites,values_from = concentration)

limits_habitat <- limits%>%
  filter(metNames_filtered %in% x)%>%
  rename(metabolites = metNames_filtered)

habitat_long <- merge(x = limits_habitat, y = habitat_long, by = "metabolites", all = TRUE)%>%
mutate(QuantFlag = concentration >= LOQ_filtered_nM, .after = "LOQ_filtered_nM")%>%
  mutate(LODflag = concentration >= LOD_filtered_nM, .after = "LOD_filtered_nM")%>%
  rename(LOD = LOD_filtered_nM, LOQ = LOQ_filtered_nM)

rm(limits,data_wFileName, habitat_wide)

# Subset benthic data -----------------------------------------------------
fivereefs <- benthic %>%
  filter(Site %in% c("DL", "CO", "JS", "YZ", "TK")) %>%
  mutate(Site = factor(Site, levels = c("DL", "CO", "JS", "YZ", "TK")))%>%
  mutate(bay_group =
           case_when(Site =="DL" ~ "fish_bay", 
                     Site =="CO" ~ "fish_bay",
                     Site =="JS" ~ "fish_bay",
                     Site =="YZ" ~ "lameshur_bay",
                     Site =="TK" ~ "lameshur_bay"))%>%
  arrange(desc(Site))

fivereefs_long <- fivereefs%>%
  pivot_longer(cols=4:14,names_to="Category",values_to="prop")%>%
  mutate(counts = prop*100)%>%
  filter(prop > 0)%>%
  mutate(Category = case_when(Category == "DiseasedCoral" ~ "Diseased coral", 
                              Category == "HardCoral" ~ "Hard coral", 
                              Category == "SoftCoral" ~ "Soft coral", 
                              Category == "TurfAlgae" ~ "Turf algae", TRUE ~ Category))

# Kruskal Wallis & Pairwise Wilcoxon Test of FCM & Water Chemistry by Site --------

metadata.mat <- habitat[,13:18]

set.seed(1)
kw <- apply(metadata.mat,2,kruskal.test, g = habitat$Site_abv)

kw

kw_results<- rbindlist(kw,idcol="metadata")%>%
  mutate(sig=p.value <= 0.05)

set.seed(1)
pwrst <- apply(metadata.mat,2,pairwise.wilcox.test, g = habitat$Site_abv, p.adjust.method="BH",paired= FALSE)

pwrst

pwrst_results<- rbindlist(pwrst,use.names=TRUE,idcol="metadata")%>%
  mutate("sig <= 0.05" = p.value <= 0.05,
         "sig <= 0.1" = p.value <= 0.1)

rm(kw, kw_results,pwrst,pwrst_results)

# Kruskal Wallis & Pairwise Wilcoxon Test of FCM & Water Chemistry by Bay --------

metadata.mat <- habitat[,13:18]

set.seed(1)
kw <- apply(metadata.mat,2,kruskal.test, g = habitat$bay_group)

kw

kw_results<- rbindlist(kw,idcol="metadata")%>%
  mutate(sig=p.value <= 0.05)

set.seed(1)
pwrst <- apply(metadata.mat,2,pairwise.wilcox.test, g = habitat$bay_group, p.adjust.method="BH",paired= FALSE)

pwrst

pwrst_results<- rbindlist(pwrst,use.names=TRUE,idcol="metadata")%>%
  mutate("sig <= 0.05" = p.value <= 0.05,
         "sig <= 0.1" = p.value <= 0.1)

rm(kw, kw_results,pwrst,pwrst_results)

# Figure 1. St.John,USVI Map ----------------------------------------------

# Raw shapefiles
usa <- st_read("site_map/stanford-vt021tk4894-shapefile/","vt021tk4894")
sttstj <- st_read("site_map/stsj_fin")
nps <- st_read("site_map/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service")

# site metadata

reefs <- read.table("site_map/CINARhabitat_reef_coords.txt", sep = "\t", header = TRUE)
reefs$Site <- factor(reefs$Site, levels = c("Ditliff", "Cocoloba", "Joels Shoal", "Yawzi", "Tektite"))
usvi <- usa %>% filter(state == "United States Virgin Islands")

reefZONE <- sttstj %>%
  filter(ZONE %in% c("Forereef", "Reef Crest", "Backreef"))

vinp <- nps %>% filter(PARKNAME == "Virgin Islands")

#plot
Figure1 <- ggplot() +
  geom_sf(data = usvi, fill = "#e5dccd",color = "#6699CC") +
  geom_sf(data = reefZONE, fill = "#cc6677", color = "#882255",alpha=0.5) +
  geom_sf(data = vinp, fill = NA, color = "#117733", linewidth = 0.6) +
  coord_sf(xlim = c(-64.8040, -64.655), ylim = c(18.2895, 18.378364), expand = FALSE) +
  geom_point(data = reefs, mapping = aes(x = Lon, y = Lat, fill = Site),colour="#333333", pch = 21, size = 3,alpha=0.9) +
  scale_fill_manual(values = CINARpalette) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "#6699CC"),legend.position = c(0.902, 0.21))+
  labs(x = "Longitude",y="Latitude")+
  geom_text() +
  annotate("text", label = "St. John,\nU.S. Virgin Islands",x = -64.74, y = 18.34, size = 4, colour = "black",fontface=2)+
  annotate("text", label = "Coral reef habitat",x = -64.785, y = 18.307, size = 3, colour = "#882255",fontface=2)+
  annotate("text", label = "Virgin Islands \nNational Park",x = -64.715, y = 18.3585, size = 3, colour = "#117733",fontface=2)+
  annotate("text", label = "Fish Bay",x = -64.757, y = 18.309, size = 3, colour = "black",fontface=2)+
  annotate("text", label = "Lameshur Bay",x = -64.735, y = 18.305, size = 3, colour = "black",fontface=2)

Figure1 + geom_segment(aes(x = -64.735, y = 18.307, xend = -64.728, yend = 18.312),
                 arrow = arrow(length = unit(0.01, "npc")))

# Figure 3. Heatmap of Metabolite Concentration by Site w/Clustering by Concentration----------------------------------------------------

#Create a dataframe in a wide format with the median concentration of each metabolite by Site. 
df<-habitat_long%>%
  mutate(metabolites = case_when(metabolites %in%
                                   c("lysine 2","ornithine 2","cysteine 2",
                                     "glutathione 2","putrescine 2",
                                     "spermidine 3") 
                                 ~ str_sub(metabolites,end=-2), 
                                 TRUE ~ metabolites))%>%
  mutate(metabolites = case_when(metabolites == "ectoine plus H2O" ~ "ectoine",
                                 TRUE ~ metabolites))%>%
  select(Site_abv,metabolites,concentration)%>%
  group_by(Site_abv,metabolites)%>%
  summarise(median_conc = median(concentration))%>%
  pivot_wider(names_from = Site_abv, values_from = median_conc)

# Run clustering
data_matrix <- as.matrix(df[, -c(1)]) #remove metabolite column
rownames(data_matrix) <- df$metabolites #set metabolite names as rownames
dendro <- as.dendrogram(hclust(d = dist(x = data_matrix)))

# Create dendrogram plot
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE) + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

print(dendro_plot)

# Heatmap

# Convert dataframe to long format
df_long <- pivot_longer(data = df,
                        cols = -c(metabolites),
                        names_to = "Site_abv",
                        values_to = "median_conc")

# Extract the order of the tips in the dendrogram
order <- order.dendrogram(dendro)

#Merge in the LOD values by matching metabolite names

limits_heatmap <- limits_habitat%>%
  mutate(metabolites = case_when(metabolites %in%
                                   c("lysine 2","ornithine 2","cysteine 2",
                                     "glutathione 2","putrescine 2",
                                     "spermidine 3") 
                                 ~ str_sub(metabolites,end=-2), 
                                 TRUE ~ metabolites))%>%
  mutate(metabolites = case_when(metabolites == "ectoine plus H2O" ~ "ectoine",
                                 TRUE ~ metabolites))

df_long <- merge(limits_heatmap,df_long,by="metabolites")

#Replace values below the LOD with NA so they will be white in the heatmap
for (i in 1:nrow(df_long)) {
  # Check if concentration is less than LOD
  if (df_long$median_conc[i] < df_long$LOD_filtered_nM[i]) {
    # Replace concentration with NA
    df_long$median_conc[i] <- NA
  }
}

# Order the levels according to their position in the cluster so the heatmap and dendogram are in the same order

df_long$metabolites <- factor(x = df_long$metabolites,
                              levels = df$metabolites[order], 
                              ordered = TRUE)

df_long$Site_abv <-   factor(df_long$Site_abv, levels = c("DL", "CO", "JS", "YZ", "TK"))


str(df_long) #make sure this is a factored ordered dataframe 

# Create heatmap plot
Figure3 <- df_long%>%
  ggplot(aes(x = Site_abv, y = metabolites)) +
  geom_tile(color="black",aes(fill = median_conc)) +
  scale_fill_viridis(trans = 'log', na.value = "white", direction = 1,
                     breaks = c(0.001, 0.01, 0.1, 1, 8), "Median Concentration (nM)") + 
  theme_bw(base_size = 14)+
  theme(legend.position = "top",text = element_text(colour = "black"))+
  labs(x = "Reef", y="Metabolite")

print(Figure3)

# All together (manually need to adjust the dendo_plot y and height values to get it to align with the heatmap, kind of tedious but it works)
grid.newpage() #need to run this each iteration or it will plot ontop of the old plot
print(dendro_plot, 
      vp = viewport(x = 0.85, y = 0.465, width = 0.3, height = 0.99))
print(Figure3, 
      vp = viewport(x = 0.35, y = 0.5, width = 0.75, height = 1.0))

# Figure 4. Metabolite Comparison -----------------------------------------

metcomp_palette <-c("#ca0020", "#f4a582", "#92c5de","#0571b0","#033F63")

# Load in datasets
d1 <- read_csv("met_comp/CINAR_20230606.csv")
d2 <- read_csv("met_comp/Weber2022.csv")
d3 <- read_csv("met_comp/Weber2020.csv")
d4 <- read_csv("met_comp/Fiore2017.csv")
d5 <- read_csv("met_comp/Becker2023.csv")

# Extract all unique InchiKey (metabolites) from the four respective datasets (d1-d4)
x <- as.data.frame(c(d1$`Chemical Name`,d2$`Chemical Name`,d3$`Chemical Name`,d4$`Chemical Name`,d5$`Chemical Name`))
y <- as.data.frame(c(d1$InChIKey,d2$InChIKey,d3$InChIKey,d4$InChIKey,d5$InChIKey))

mets <- cbind(x,y)%>%
  rename("Modified_chemical_name" = "c(d1$`Chemical Name`, d2$`Chemical Name`, d3$`Chemical Name`, d4$`Chemical Name`, d5$`Chemical Name`)", "InchiKeys" = "c(d1$InChIKey, d2$InChIKey, d3$InChIKey, d4$InChIKey, d5$InChIKey)")%>%
  drop_na()%>%
  distinct()

met_comp<- mets%>%
  mutate(Garcia_2023 = InchiKeys %in% d1$InChIKey, Weber_2022 = InchiKeys %in% d2$InChIKey,
         Weber_2020 = InchiKeys %in% d3$InChIKey,Fiore_2017 = InchiKeys %in% d4$InChIKey, Becker_2023 = InchiKeys %in% d5$InChIKey)%>%
  mutate(Garcia_2023 = as.integer(as.logical(Garcia_2023)),
                                Weber_2022 = as.integer(as.logical(Weber_2022)),
                                Weber_2020 = as.integer(as.logical(Weber_2020)),
                                Fiore_2017 = as.integer(as.logical(Fiore_2017)),
                                Becker_2023 = as.integer(as.logical(Becker_2023)))

str(met_comp)
rm(x,y,d1,d2,d3,d4,d5)

Figure4 <- upset(met_comp, nsets=5, nintersects = 100,
                    sets = c("Garcia_2023", "Weber_2020", "Weber_2022","Fiore_2017","Becker_2023"),show.numbers = "yes")

print(Figure4)

# PERMANOVA on metabolite and benthic composition ---------------------------------------------------------------

#  PERMANOVA on metabolite data

set.seed(100)
site <- adonis2(habitat.mat ~ Site_abv, data = habitat, method = "bray")
site

set.seed(100)
bay <- adonis2(habitat.mat ~ bay_group, data = habitat, method = "bray")
bay

#  PERMANOVA on benthic data

fivereefs.mat <- fivereefs[,c(4:14)]

set.seed(100)
site <- adonis2(fivereefs.mat ~ Site, data = fivereefs, method = "bray")
site

set.seed(100)
bay <- adonis2(fivereefs.mat ~ bay_group, data = fivereefs, method = "bray")
bay

rm(fivereefs.mat)
# Figure 5. Spearman Correlation between Benthic Survey & Metabolites --------------------------------------------------------------------

fivereefs_cor <- benthic %>%
  filter(Site %in% c("DL", "CO", "JS", "YZ", "TK")) %>%
  unite("row_id",Site:Transect, remove = FALSE)%>%
  mutate(Site = factor(Site, levels = c("DL", "CO", "JS", "YZ", "TK")))%>%
  select(Site,CCA:TurfAlgae)%>%
  group_by(Site)%>%
  summarise_all(median)%>%
  rename("Diseased Coral" = "DiseasedCoral",
         "Hard Coral" = "HardCoral",
         "Soft Coral" = "SoftCoral",
         "Turf Algae" = "TurfAlgae")

fivereefs_mat <- as.matrix(fivereefs_cor[,-1])
rownames(fivereefs_mat) <- as.matrix(fivereefs_cor[,1])

x<- colnames(habitat.mat)

habitat_cor <- habitat%>%
  group_by(Site_abv)%>%
  select(Site_abv,all_of(x))%>%
  summarise_all(median)%>%
  rename("cysteine" = "cysteine 2",
         "ectoine" = "ectoine plus H2O",
         "glutathione" = "glutathione 2",
         "ornithine" = "ornithine 2",
         "spermidine" =  "spermidine 3",
         "lysine" = "lysine 2",
         "putrescine" = "putrescine 2")

habitat_cor_mat <- as.matrix(habitat_cor[,-1])
rownames(habitat_cor_mat) <- as.matrix(habitat_cor[,1])  

x <- colnames(habitat_cor_mat)

habitat_cor_mat <- habitat_cor_mat[, ! apply(habitat_cor_mat , 2 , function(x) 
  sd(x, na.rm = TRUE)==0 ) ] #remove benthic categories with a SD of zero 

#perform correlation - rcorr is nice because it returns your deg of freedom, and pvalue
# output is a list, so just keeping the correlation result for now since we aren't using the other info

set.seed(100)
df.cor.list <- rcorr(habitat_cor_mat,fivereefs_mat, type = "spearman")
df.cor <- df.cor.list[[1]]
df.cor.p <- df.cor.list[[3]]

melted_cor.p <- melt(df.cor.p)%>%
  mutate(Var1_type =
           case_when(Var1 %in% x ~ "metabolite"),
         Var2_type =
           case_when(Var2 %in% x ~ "metabolite"))%>%
  mutate(Var1_type = ifelse(is.na(Var1_type), "benthic", Var1_type))%>%
  mutate(Var2_type = ifelse(is.na(Var2_type), "benthic", Var2_type))%>%
  mutate(type = paste(Var1_type, Var2_type, sep = "-"))%>%
  filter(type == "benthic-metabolite")%>%
  rename(pvalue = value)%>%
  unite("Pair",Var1:Var2,remove=FALSE)

#this function will calculate the distance matrix, cluster, and reorder your correlation matrix
reorder_cormat <- function(cormat){
  
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]}

# Reorder the correlation matrix
cormat <- reorder_cormat(df.cor)
dd <- as.dist((1-df.cor)/2)
hc <- hclust(dd)

# Melt the correlation matrix
melted_cormat <- melt(cormat)

melted_cormat$Var1 <- as.character(melted_cormat$Var1)
melted_cormat$Var2 <- as.character(melted_cormat$Var2)

melted_cormat <- melted_cormat%>%
  mutate(Var1_type =
           case_when(Var1 %in% x ~ "metabolite"),
         Var2_type =
           case_when(Var2 %in% x ~ "metabolite"))%>%
  mutate(Var1_type = ifelse(is.na(Var1_type), "benthic", Var1_type))%>%
  mutate(Var2_type = ifelse(is.na(Var2_type), "benthic", Var2_type))%>%
  mutate(type = paste(Var1_type, Var2_type, sep = "-"))%>%
  filter(type == "benthic-metabolite")%>%
  rename(corr_value = value)%>%
  unite("Pair",Var1:Var2,remove=FALSE)

idx <- which(! hc$labels[hc$order] %in% x)
Var1.levels <- hc$labels[hc$order]
Var1.levels <- Var1.levels[idx]

idx <- which(hc$labels[hc$order] %in% x)
Var2.levels <- hc$labels[hc$order]
Var2.levels <- Var2.levels[idx]

corrPlot <- melted_cormat %>%
  merge(melted_cor.p,by="Pair")%>%
  select(Pair,Var1.x,Var2.x,corr_value,pvalue)%>%
  rename(Var1 = Var1.x, Var2 = Var2.x)%>%
  mutate(Var1 = factor(Var1, levels = Var1.levels),
         Var2 = factor(Var2, levels = Var2.levels)) %>%
  filter(Var1 != "Other")%>%
  mutate(sig = pvalue < 0.05)

corrPlot$Flag <- ifelse(corrPlot$pvalue < 0.05, "*", "")

Figure5 <- corrPlot%>%
  ggplot(aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = corr_value), colour = "grey45") + 
  coord_equal() + 
  scale_fill_viridis(na.value = "white", direction = 1,
                     breaks = c(-1.0, -0.5, -0.0, 0.5, 1.0), "Spearman's Correlation")+
  theme_bw(base_size = 14)+
  theme(legend.position = "top", axis.text.y = element_text(color="black"),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, color="black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank())+
  geom_text(aes(label = Flag), color = "darkorange2",fontface = "bold",size=8,vjust = 0.7) +  # Add text labels
  labs(x= "Metabolites", y = "Benthic Categories", fill = "Spearman's Correlation")

print(Figure5)

corrSummary <- merge(melted_cor.p,melted_cormat, by="Pair")
colnames(corrSummary)

corrSummary <- corrSummary%>%
  select(Pair,pvalue,corr_value,type.x)%>%
  rename(pair_type = type.x)%>%
  mutate(sig = pvalue < 0.05, 
         strong = abs(corr_value) >= 0.8)

# Figure 6. Kruskall-Wallis test and Pairwise Wilcoxon Rank Sum Test of Metabolites by Site ----------------------------------------------------

#Kruskal Wallis 
rm(kw, kw_results, pwrst, pwrst_results, pwrst_sig)

set.seed(100)
kw <- apply(habitat.mat,2,kruskal.test, g = habitat$Site_abv)

kw

kw_results<- rbindlist(kw,idcol="metabolite")%>%
  mutate(sig=p.value <= 0.05)

sum(kw_results$sig == TRUE,na.rm=TRUE)
index <- which(kw_results$sig == TRUE)
df <- data.frame(metabolite=kw_results$metabolite[index],pvalue_kw = kw_results$p.value[index])

x <- colnames(habitat.mat) #pull out metabolite names

kw_sig_Site <- habitat.mat %>%
  select(x[index])

#Pairwise Wilcoxon Rank Sum Test 

set.seed(100)
pwrst <- apply(kw_sig_Site,2,pairwise.wilcox.test, g = habitat$Site_abv, p.adjust.method="BH",paired= FALSE)

pwrst

pwrst_results<- rbindlist(pwrst,idcol="metabolite")%>%
  mutate(sig=p.value <= 0.1)

pwrst_sig <- pwrst_results%>%
  filter(sig==TRUE)

length(unique(pwrst_sig$metabolite))

id <- unique(pwrst_sig$metabolite)

Figure6 <- habitat_long%>%
  filter(metabolites %in%  id)%>%
  mutate(mtab_bc = metabolites, mtab_bc = case_when(mtab_bc == "homoserine betaine" ~ "homoserine betaine*", TRUE ~ mtab_bc),.after = metabolites)%>%
  ggplot(aes(x=Site_abv, y=concentration)) +
  geom_boxplot(lwd=0.25,fatten=1.5)+
  geom_jitter(aes(color= Site_abv, shape = bay_group), size=2.5, alpha=0.7) +
  scale_color_manual("Reef",values=CINARpalette)+
  facet_wrap(~mtab_bc,scales="free_y")+
  labs(x="Reef", y = "Concentration (nM)", shape = "Bay")+
  theme_bw(base_size = 14)+
  theme(legend.position = "top")

print(Figure6)

# Figure 7. Kruskall-Wallis test and Pairwise Wilcoxon Rank Sum Test of Metabolites by Bay ----------------------------------------------------

#Kruskal Wallis 
rm(kw, kw_results, pwrst, pwrst_results, pwrst_sig)

set.seed()
kw <- apply(habitat.mat,2,kruskal.test, g = habitat$bay_group)

kw

kw_results<- rbindlist(kw,idcol="metabolite")%>%
  mutate(sig=p.value <= 0.05)

sum(kw_results$sig == TRUE,na.rm=TRUE)
index <- which(kw_results$sig == TRUE)

kw_sig_Bay <- habitat %>%
  select(x[index])

#Pairwise Wilcoxon Rank Sum Test 

set.seed(100)
pwrst <- apply(kw_sig_Bay,2,pairwise.wilcox.test, g = habitat$bay_group, p.adjust.method="BH",paired= FALSE)

pwrst

pwrst_results<- rbindlist(pwrst,idcol="metabolite")%>%
  mutate(sig=p.value <= 0.05)

pwrst_sig <- pwrst_results%>%
  filter(sig==TRUE)

length(unique(pwrst_sig$metabolite))

id <- unique(pwrst_sig$metabolite)

Figure7 <- habitat_long%>%
  filter(metabolites %in%  id)%>%
  mutate(metabolites = case_when(metabolites == "glutathione 2" ~ "glutathione",TRUE ~ metabolites))%>%
  mutate(mtab_bc = metabolites, mtab_bc = case_when(mtab_bc == "cysteate" ~ "cysteate*",
                                          mtab_bc == "glutamine" ~"glutamine*",
                                          mtab_bc == "isethionate" ~ "isethionate*", TRUE ~ mtab_bc),.after = metabolites)%>%
  mutate(mtab_bc = factor(mtab_bc, levels = c("cysteate*","glutamine*","isethionate*",
                                    "glutamic acid","glutathione","guanosine","sn-glycerol 3-phosphate",
                                    "tryptophan")))%>%
  ggplot(aes(x=bay_group, y=concentration)) +
  geom_boxplot(lwd=0.25,fatten=1.5)+
  geom_jitter(aes(color= Site_abv,shape= bay_group), size=2.5, alpha=0.7) +
  scale_color_manual("Reef",values = CINARpalette)+
  facet_wrap(~mtab_bc,scales="free_y",nrow=2)+
  labs(x="Bay", y = "Concentration (nM)", shape = "Bay")+
  theme_bw(base_size = 14)+
  theme(legend.position = "top")

print(Figure7)

# Figure S1. Metabolite Comparison Bubbleplot -----------------------------

met_comp_longer <- as_tibble(met_comp)%>%
  mutate(sum = rowSums(.[3:7]))%>%
  pivot_longer(cols = c("Garcia_2023", "Weber_2022", "Weber_2020", "Fiore_2017","Becker_2023"),
               names_to="study",
               values_to="value",values_drop_na=FALSE)%>%
  mutate(study = factor(study,level = c("Fiore_2017","Weber_2022","Garcia_2023","Weber_2020","Becker_2023")))

tmp <- met_comp_longer%>%
  arrange(desc(sum))%>%
  select(Modified_chemical_name)

FigureS1 <- met_comp_longer%>%
  mutate(Modified_chemical_name = factor(Modified_chemical_name, levels = unique(tmp$Modified_chemical_name)))%>%
  mutate(sum = factor(sum, levels = c(1,2,3,4,5)))%>%
  ggplot(aes(y=Modified_chemical_name,x=study,size= value,fill=study)) +
  geom_point(shape=21)+
  scale_size(range = c(0,3), guide="none")+
  geom_hline(yintercept =c(6.5,17.5,31.5,47.5),linewidth=.5)+
  theme(axis.text.y=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="Metabolite",x="Study")+
  scale_fill_manual("Study",values=metcomp_palette, labels = c("Fiore et al. 2017","Weber et al. 2022","*Garcia et al. 2023","Weber et al. 2020","Becker et al. 2023"))+
  coord_fixed(ratio = .3)

print(FigureS1)

# Figure S2-A. PCA of Benthic Data ------------------------------------------------

fivereefs <- benthic %>%
  filter(Site %in% c("DL", "CO", "JS", "YZ", "TK")) %>%
  unite("row_id",Site:Transect, remove = FALSE)%>%
  mutate(Site = factor(Site, levels = c("DL", "CO", "JS", "YZ", "TK")))%>%
  rename(Site_abv = Site)%>%
  mutate(bay_group =
           case_when(Site_abv =="DL" ~ "Fish Bay", 
                     Site_abv =="CO" ~ "Fish Bay",
                     Site_abv =="JS" ~ "Fish Bay",
                     Site_abv =="YZ" ~ "Lameshur Bay",
                     Site_abv =="TK" ~ "Lameshur Bay"),.after = Site_abv)%>%
  mutate(bay_group = factor(bay_group, levels = c("Fish Bay", "Lameshur Bay")))

pca_data_short_final <- as.matrix(fivereefs)                #change from tibble to matrix
rownames(pca_data_short_final) <- pca_data_short_final[,2]  #make rownames the unique sites
pca_data_short_final <- pca_data_short_final[,6:16]         #select numeric variables
class(pca_data_short_final) <- "numeric"                    #change from character to numeric values for pca

#Do the PCA
pca <- PCA(pca_data_short_final,scale.unit=TRUE) #performs the Principal component analysis. #scale.unit=TRUE then data are scaled to unit variance

# Graph
FigureS2A <- fviz_pca(pca,
                     geom=c("point"),
                     pointsize=3,
                     col.ind=fivereefs$Site_abv,
                     invisible="quali",
                     geom.var =c("arrow", "text"),
                     palette=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288"),
                     col.var = "black",
                     repel=TRUE,
                     legend.title="Reef",
                     title = "Benthic",ggtheme = theme_bw(base_size=12))+           
                    scale_shape_manual(values=c(16,16,16,17,17))

FigureS2A

# Figure S2-B. NMDS of Metabolite Data ------------------------------------------------
dist <- vegdist(habitat.mat,dist='bray')

set.seed(1)
nmds <-
  metaMDS(dist,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(nmds)
stressplot(nmds)

#need to calculate the distance between the points from the ordination itself
#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

#calculate data.scores as a Euclidean distance
distanceXY = dist(data.scores, method = "euclidean")
cf = mantel(dist,distanceXY,method = "spearman",permutations = 9, na.rm = TRUE)
r2_overall <- cf$statistic * cf$statistic
r2_overall

distanceXY = dist(data.scores[,1], method = "euclidean")
cfMatrices = mantel(dist,distanceXY,method = "spearman", permutations = 9, na.rm = TRUE)
r2_axis1 <- cfMatrices$statistic * cfMatrices$statistic
r2_axis1

distanceXY = dist(data.scores[,2], method = "euclidean")
cfMatrices = mantel(dist,distanceXY,method = "spearman", permutations = 9, na.rm = TRUE)
r2_axis2 <- cfMatrices$statistic * cfMatrices$statistic
r2_axis2

nmds_df <- data.frame(habitat$Site_abv,habitat$bay_group,nmds$points)

mat<- nmds_df%>%
  rename(Site_abv = habitat.Site_abv, bay_group = habitat.bay_group)%>%
  mutate(color_site = case_when(Site_abv == "DL" ~ "#88CCEE",Site_abv == "CO" ~  "#CC6677",
                                Site_abv == "JS" ~  "#DDCC77", Site_abv == "YZ" ~  "#117733",
                                Site_abv == "TK" ~ "#332288"))

FigureS2B <- ggscatter(mat, x = "MDS1", y = "MDS2", color = "Site_abv",size=3,shape="bay_group")+
  geom_hline(yintercept=0,linetype=2)+
  geom_vline(xintercept=0,linetype=2)+
  scale_color_manual("Reef",values=CINARpalette, guide = "none")+
  labs(title="Metabolites", shape = "Bay", 
       x = paste0("MDS1 (",round(r2_axis1*100,1), "%)"),
       y = paste0("MDS2 (",round(r2_axis2*100,1), "%)"))+
  theme_bw(base_size = 12)


print(FigureS2B)

# Figure S2. (Full). Create single plot with multi-panels for all PCA/NMDS  --------

FigureS2 <- ggarrange(FigureS2A,FigureS2B,heights = c(3,3),
                     labels = c("A","B"),
                     ncol = 2, nrow = 1,common.legend = FALSE, legend = "bottom")

FigureS2

rm(fivereefs,nmds,nmds_df,pca,pca_data_short_final,mat,dist)

# Figure S3. Barplot of FCM and Water Chemistry Data --------------------------------------

bp_data <- habitat%>%
  rename(`Unpigmented cells` = Unpigmented_cells,
         TOC = NPOC)%>%
  pivot_longer(cols=13:18,names_to="metadata",values_to="value",values_drop_na=TRUE)%>%
  mutate(metadata = factor(metadata, levels = c("Picoeukaryotes", "Prochlorococcus", "Synechococcus",
                                                   "Unpigmented cells", "TN", "TOC")))%>%
  group_by(Site_abv,metadata)%>%
  summarise(sd = sd(value, na.rm=TRUE),avg_value = mean(value, na.rm=TRUE))

FigureS3 <- ggplot(bp_data, aes(x=Site_abv, y=avg_value, fill=Site_abv)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=avg_value-sd, ymax=avg_value+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual("Reef",values=CINARpalette)+
  facet_wrap(~metadata,scales="free")+
  labs(x = "Reef",y = "Average value")+
  theme_bw(base_size=12)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

print(FigureS3)

# Figure S4. Proportion of hard+soft coral --------------------------------

FigureS4 <- fivereefs_long%>%
  filter(Category %in% c("Hard coral", "Soft coral"))%>%
  group_by(Site,Transect,Category)%>%
  summarise_at("prop",mean)%>%
  summarise_at("prop", sum)%>%
  rename(sum_prop_Coral = prop)%>%
  ggplot(aes(x=Site, y=sum_prop_Coral, fill=Site))+
  geom_boxplot(alpha=0.8)+
  geom_point()+
  scale_fill_manual("Reef",values=CINARpalette)+
  labs(x="Reef", y="Proportion of coral (hard + soft)")+
  theme_bw(base_size=14)

print(FigureS4)

# Figure S5. Benthic Composition ------------------------------------------

FigureS5 <- fivereefs_long%>%
  select(Site,Category,prop)%>%
  group_by(Site,Category)%>%
  summarise_all(mean)%>%
  ggplot(aes(fill=Category, y=prop, x=Site)) + 
  geom_bar(position="fill", stat="identity", color="grey13")+
  scale_fill_brewer(type = "div", palette = "RdYlBu")+
  labs(y="Proportion",x="Reef")+
  theme_bw(base_size=14)

print(FigureS5)

# Figure S6. Ratio hard coral to macroalgae -------------------------------

coral_vs_algae <- fivereefs_long%>%
  filter(Category %in% c("Hard coral", "Macroalgae"))%>%
  group_by(Site,Transect,Category)%>%
  summarise_at("prop",mean)%>%
  rename(avg_prop = prop)%>%
  pivot_wider(names_from = Category, values_from = avg_prop)%>%
  mutate(ratio = `Hard coral` / Macroalgae)%>%
  group_by(Site)%>%
  mutate(mean_ratio = mean(ratio, na.rm = TRUE), med_ratio = median(ratio, na.rm = TRUE))

FigureS6 <- coral_vs_algae %>%
  ggplot(aes(x=Site, y=ratio, fill=Site))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  geom_boxplot(alpha=0.8)+
  geom_point()+
  geom_point(aes(x=Site, y= mean_ratio), fill="white",color="grey12", shape = 23, size=3)+
  scale_fill_manual("Reef",values=CINARpalette)+
  labs(x="Reef", y= "Ratio (Hard coral:Macroalgae)")+
  theme_bw(base_size=14)

print(FigureS6)
