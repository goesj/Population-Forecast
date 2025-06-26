##### Script to create all graphics within the text ############################
pacman::p_load("tidyverse","reshape2","ggplot2", "cowplot","ggdist","ggrepel",
               "paletteer")


options(mc.cores = parallel::detectCores())
source("02_Functions.R")

load(file = file.path(getwd(),"Data/TotalData.RData"))



##### Population Projections ###################################################
load(file = file.path(getwd(),"Results/PopTrajectories_RaTotal_24.RData"))

PopTrajectories_Total <- PopTrajectories_Male + PopTrajectories_Female


RegionNames_English <- gsub("(Lkr)", "Region", unique(PopCounts$RegionName)) #Change Lkr to Region
RegionNames_English <- gsub("(Krfr.St)", "City", RegionNames_English) #Change Krfr. St to City
RegionNames_English[13] <- c("Wunsiedel (Region)")



#### 1. Plot the population of all regions together ###########################
#Summarise forecasts and transform into long format for plotting
PopLong_FC_Tot <- 
  PopTrajectories_Total %>% 
  reshape2::melt(value.name = "Population", as.is = TRUE) %>% 
  reframe("Population" = sum(Population), 
          .by = c(Year, RegionNumber, Iteration)) %>% #sum over age
  reframe(ggdist::median_qi(Population, .width = c(0.8, 0.95)), 
          .by = c(Year, RegionNumber)) %>% #calculate PI's
  mutate("Year" = as.numeric(Year)) %>% 
  filter(Year < 2048) %>% 
  mutate(".width" = as.factor(.width)) %>% 
  rename(Population = y) %>% 
  left_join(x = ., 
            y = select(PopCounts, c(RegionNumber, RegionName)), 
            multiple = "first") %>% #add RegionNumber 
  mutate("Type" = "FC") %>% 
  mutate(RegionID = match(RegionNumber, unique(RegionNumber)), #add english region Names
         RegionName = RegionNames_English[RegionID]) %>% 
  select(-RegionID)
  

PopFC_Total_Plot <- 
  PopCounts %>% 
  reframe(Population = sum(Population),  #sum over sex
          .by = c(Year, RegionNumber, RegionName)) %>% 
  mutate("Type" = "Obs") %>% 
  mutate(RegionID = match(RegionNumber, unique(RegionNumber)),
         RegionName = RegionNames_English[RegionID]) %>% 
  ggplot(aes(x = Year, y = Population))+
  geom_line(aes(linetype = Type),
            linewidth = 1)+
  ggdist::geom_lineribbon(data = PopLong_FC_Tot, 
                          aes(y = Population, ymin = ymin, ymax = ymax, 
                              fill = .width,
                              linetype = Type), alpha = 0.8,
                          linewidth = 1)+
  facet_wrap(~ RegionName, scales ="free_y")+
  theme_bw()+
  scale_linetype_manual(values = c("dashed","solid"))+
  scale_fill_manual(values = c("#4292C6","#97BCDD"))+
  scale_y_continuous(labels = scales::label_number())+ #no scientific notation
  theme(legend.key.width = unit(1,"cm"),
        axis.text = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),#rotate text
        axis.title = element_text(size = 20 , face = "bold"),
        legend.title = element_text(size = 20, face ="bold"))+ 
  guides(fill = guide_legend(title = "PI",
                             override.aes = list(linetype = 0)), #remove line 
         linetype = guide_legend(override.aes =   #change appearance in legend
                                   list(fill = NA, 
                                        linetype = c("dashed","solid"),
                                        color = "black")))

ggsave(filename = file.path(getwd(),"Pictures/PopTotal_OF.pdf"),
       plot = PopFC_Total_Plot, 
       device = "pdf", dpi = 500, 
       width = 16, height = 12)

###### 2. Plot Population Pyramid #############################################
# @Region should be number between 1 and 13
# @Year should be either 2029, 2034, 2039 or 2044

PP_Bayreuth <- 
  PopulationPyramid(dataF = PopTrajectories_Female,
                    dataM = PopTrajectories_Male, 
                    Region = 2, YearUpper = 2044) #select region and year


ggsave(filename = file.path(getwd(),"Pictures/PopPyramid.pdf"),
       plot = PP_Bayreuth, device = "pdf", dpi = 500, 
       width = 16, height = 12)


##### 4. Map plot of Difference in Population Estimates ########################
PopLong_FC_Diff <- 
  PopTrajectories_Total %>% 
  reshape2::melt(value.name = "Population", as.is = TRUE) %>% 
  reframe("Population" = sum(Population), 
          .by = c(Year, RegionNumber, Iteration)) %>% #sum over age
  reframe(ggdist::median_qi(Population, .width = c(0.8)), 
          .by = c(Year, RegionNumber)) %>% #calculate PI's
  mutate("Year" = as.numeric(Year)) %>% 
  mutate(".width" = as.factor(.width)) %>% 
  rename(Population = y) %>% 
  filter(Year %in% c(2024, 2044)) %>% 
  reframe("PopDiff" = diff(Population), 
          .by = RegionNumber) %>% 
  rename(ARS = RegionNumber)


NamesData <-
  UpperFranconia %>% 
  mutate("Region" = c(rep(NA, 4),tail(RegionNames_English,9))) %>% 
  mutate("City" = c(RegionNames_English[1:4], rep(NA, 9))) %>% 
  sf::st_centroid(.) %>% #Get middle point of Area (for Plotting text)
  tidyr::extract(geometry, into = c('Lat', 'Lon'), '\\((.*),(.*)\\)', conv = TRUE)


## Merge Data with Res
pacman::p_load(paletteer)
Merged_MapData <- 
  left_join(x = UpperFranconia, 
            y = PopLong_FC_Diff) %>% 
  mutate("Region" = c(rep(NA, 4),
                      tail( #Remove (Region from name for nicer plotting)
                        unlist(strsplit(RegionNames_English, " (Region)", fixed = TRUE)),
                        9))) 

MyPalAdj <- c(paletteer_d("MexBrewer::Revolucion")[1:4],
              "#FEE08B", "#FFFFBF", #RdYlBu
              "#D1E5F0" ,"#92C5DE" ,"#4393C3", "#2166AC") #RdBu

PopDiff_Plot <- 
  ggplot(Merged_MapData)+
  geom_sf(aes(fill = PopDiff), linewidth = 0, alpha = 0.9,
          color = "gray50")+
  theme_void() +
  geom_sf_text(aes(label=Region, color = PopDiff < -8000),
               nudge_x = c(rep(0,4),-4000,0,-15000,-100,0,200,0,50,-50),
               nudge_y = c(rep(0,4),6000,-2000, 11000,0,-2000,0,0,-2000,0),
               size = 5.5)+
  ggrepel::geom_text_repel(data = NamesData, aes(x = Lat, y = Lon, label = City),
                           nudge_x = c(-20000, 35000, 6000, -4000, rep(0,9)),
                           nudge_y = c(-25000, -2000, 20000, 18000, rep(0,9)),
                           size = 5.5, 
                           force_pull = 0)+
  scale_fill_gradientn(str_wrap("Difference in population", width = 15), #line break 
                       colors=MyPalAdj,
                       limits = c(-11500, 10000), 
                       labels = scales::label_number())+
  scale_color_manual(guide = "none", values = c("black", "white"))+
  theme(legend.position = "bottom",
        legend.key.width = unit(1.75,"cm"),
        legend.key.height = unit(0.75,"cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))+
  guides(fill = guide_colorbar(frame.colour = "gray20", 
                               frame.linewidth = 0.8))

### Plotting share of old People ##
Share_OldPeople <- 
  PopCounts %>% 
  reframe(Population = sum(Population), 
          .by = c(Year, RegionNumber, AgeID_New, AgeGroup_New)) %>% 
  mutate("over65" = ifelse(AgeID_New > 12, 1, 0)) %>% 
  reframe("PopS" = sum(Population), .by = c(Year, RegionNumber, over65)) %>% 
  reframe("PopShare" = prop.table(PopS)*over65, 
                       .by = c(Year, RegionNumber)) %>% 
  filter(Year == 2023 & PopShare != 0) %>% 
  rename(ARS = RegionNumber)


## Merge Data with Res
Merged_MapData_Share <- 
  left_join(x = UpperFranconia, 
            y = Share_OldPeople) 

PopShare_Plot <- 
  ggplot(Merged_MapData_Share)+
  geom_sf(aes(fill = PopShare), linewidth = 0, alpha = 0.9,
          color = "gray50")+
  theme_void() +
  scale_fill_gradientn(str_wrap("Share of people above 60", width = 15), 
                       colors=MyPalAdj, 
                       limits = c(0.25, 0.36))+
  theme(legend.position = "bottom",
        legend.key.width = unit(1.75,"cm"),
        legend.key.height = unit(0.75,"cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))+
  guides(fill = guide_colorbar(frame.colour = "gray20", 
                               frame.linewidth = 0.8))


Joined_Map_Plot <- cowplot::plot_grid(PopDiff_Plot, 
                                      PopShare_Plot, 
                                      align = c("h"), 
                                      ncol = 2)


ggsave(filename = file.path(getwd(),"Pictures/PopDiff_Map.pdf"),
       plot = Joined_Map_Plot, 
       device = "pdf", dpi = 500, 
       width = 16, height = 7)


  
##### 5. Age Migration Schedule ################################################
### 5.1 Males ###
OutMigSchedule_Plot <- 
  AgeSpecific_Mig_Male %>% 
  mutate("AgeGroup_New" = factor(AgeGroup_New, levels = unique(AgeGroup_New))) %>% 
  mutate("Ra"=Out/sum(Out), 
         .by = c(Year, RegionNumber)) %>%
  mutate("RegionID" = match(RegionNumber, unique(RegionNumber)),
         "RegionName" = RegionNames_English[RegionID]) %>% 
  ggplot(aes(x = AgeGroup_New))+
  geom_line(aes(y = Ra, group = Year, col = Year))+
  facet_wrap(~RegionName, scales = "free_y")+
  scale_color_continuous(labels = scales::label_number(accuracy = 1, #no decimals points in legend
                                                       big.mark = ""))+
  xlab(label = "")+
  ylab(label ="Out-migration schedule")+
  theme_bw()+
  theme(legend.key.width = unit(3,"cm"),
        axis.text = element_text(size = 15, face ="bold"),
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),#rotate text
        axis.title = element_text(size = 20 , face = "bold"),
        legend.title = element_text(size = 20, face ="bold"),
        legend.position = "bottom")


InMigSchedulePlot <- 
  AgeSpecific_Mig_Male %>% 
  mutate("AgeGroup_New" = factor(AgeGroup_New, levels = unique(AgeGroup_New))) %>% 
  mutate("Ra_In"=In/sum(In), 
         .by = c(Year, RegionNumber)) %>%
  mutate("RegionID" = match(RegionNumber, unique(RegionNumber)),
         "RegionName" = RegionNames_English[RegionID]) %>% 
  ggplot(aes(x = AgeGroup_New))+
  geom_line(aes(y = Ra_In, group = Year, col = Year))+
  facet_wrap(~RegionName, scales = "free_y")+
  scale_color_continuous(labels = scales::label_number(accuracy = 1, #no decimals points in legend
                                                       big.mark = ""))+
  xlab(label = "")+
  ylab(label ="In-migration schedule")+
  theme_bw()+
  theme(legend.key.width = unit(3,"cm"),
        axis.text = element_text(size = 15, face ="bold"),
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),#rotate text
        axis.title = element_text(size = 20 , face = "bold"),
        legend.title = element_text(size = 20, face ="bold"),
        legend.position = "bottom")


### 5.2 Females ####
OutMigSchedule_Plot_F <- 
  AgeSpecific_Mig_Female %>% 
  mutate("AgeGroup_New" = factor(AgeGroup_New, levels = unique(AgeGroup_New))) %>% 
  mutate("Ra"=Out/sum(Out), 
         .by = c(Year, RegionNumber)) %>%
  mutate("RegionID" = match(RegionNumber, unique(RegionNumber)),
         "RegionName" = RegionNames_English[RegionID]) %>% 
  ggplot(aes(x = AgeGroup_New))+
  geom_line(aes(y = Ra, group = Year, col = Year))+
  facet_wrap(~RegionName, scales = "free_y")+
  scale_color_continuous(labels = scales::label_number(accuracy = 1, #no decimals points in legend
                                                       big.mark = ""))+
  xlab(label = "")+
  ylab(label ="Out-migration schedule")+
  theme_bw()+
  theme(legend.key.width = unit(3,"cm"),
        axis.text = element_text(size = 15, face ="bold"),
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),#rotate text
        axis.title = element_text(size = 20 , face = "bold"),
        legend.title = element_text(size = 20, face ="bold"),
        legend.position = "bottom")


InMigSchedulePlot_F <- 
  AgeSpecific_Mig_Female %>% 
  mutate("AgeGroup_New" = factor(AgeGroup_New, levels = unique(AgeGroup_New))) %>% 
  mutate("Ra_In"=In/sum(In), 
         .by = c(Year, RegionNumber)) %>%
  mutate("RegionID" = match(RegionNumber, unique(RegionNumber)),
         "RegionName" = RegionNames_English[RegionID]) %>% 
  ggplot(aes(x = AgeGroup_New))+
  geom_line(aes(y = Ra_In, group = Year, col = Year))+
  facet_wrap(~RegionName, scales = "free_y")+
  scale_color_continuous(labels = scales::label_number(accuracy = 1, #no decimals points in legend
                                                       big.mark = ""))+
  xlab(label = "")+
  ylab(label ="In-migration schedule")+
  theme_bw()+
  theme(legend.key.width = unit(3,"cm"),
        axis.text = element_text(size = 15, face ="bold"),
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),#rotate text
        axis.title = element_text(size = 20 , face = "bold"),
        legend.title = element_text(size = 20, face ="bold"),
        legend.position = "bottom")


ggsave(filename = file.path(getwd(),"Pictures/Ra_Out.pdf"),
       plot = OutMigSchedule_Plot, 
       device = "pdf", dpi = 500, 
       width = 16, height = 12)


ggsave(filename = file.path(getwd(),"Pictures/Ra_Out_F.pdf"),
       plot = OutMigSchedule_Plot_F, 
       device = "pdf", dpi = 500, 
       width = 16, height = 12)


#### 6. Net Migration Rates ####################################################
pacman::p_load("paletteer")
NetMigPlot <- 
  rbind(Male_Mig_Data, 
        Female_Mig_Data) %>% 
  ggplot(aes(x = Year, y = NetMig))+
  geom_line(aes(col = RegionNumber, group = RegionNumber))+
  facet_wrap(~Sex, scales = "free_y", 
             nrow = 2)+
  ylab(label ="Net migration counts")+
  scale_colour_paletteer_d("nbapalettes::suns_city")+
  theme_bw()+
  theme(axis.text = element_text(size = 20, face ="bold"),
        strip.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20 , face = "bold"),
        legend.title = element_text(size = 20, face ="bold"),
        legend.position = "none")

ggsave(filename = file.path(getwd(),"Pictures/NetMig_OF.pdf"), 
       plot = NetMigPlot,
       device = "pdf", dpi = 500, 
       width = 16, height = 12)


  