#### SCRIPT TO RUN MORTALITY MODEL ############################################
pacman::p_load("tidyverse","rstan","spdep","reshape2")


options(mc.cores = parallel::detectCores())
source("02_Functions.R") #load helper functions

load(file = file.path(getwd(),"Data/TotalData.RData")) #load data

###### Combine Death data with respective population data #########
DeathData <- left_join(x = DeathCounts, 
                         y = PopCounts, 
                         by = c("Year","RegionNumber","RegionName",
                                "AgeGroup_New","Sex","AgeID_New")) %>% 
  mutate("MortRate" = Deaths / Exposure) %>% 
  mutate("RegionID" = match(RegionName, unique(RegionName)),
         "YearID" = match(Year, unique(Year))) %>% 
  mutate("CohortID" = CohortIndex(agegroup = AgeID_New, time = YearID, 
                                  maxAge = max(AgeID_New),
                                  M = 5))

nm.adj <- spdep::poly2nb(UpperFranconia)
#Create Adjacency Matrix
AdjMatOF <- as(spdep::nb2mat(nm.adj, style = "B"), "Matrix") 
#Get Nodes for Stan
Nodes <- CARData4Stan(NeighborhoodMatrix = AdjMatOF) 
#scaling Factor BYM2 Model
scalingFac <- ScalingFacBYM2(Nodes = Nodes, AdjMat = AdjMatOF) 


### 01_ Females ################################################################
FemaleDeaths <- DeathData %>% #Get data
  filter(Sex == "female") 

#create stan data
StanData <- list("T" = max(FemaleDeaths$YearID),
                 "A" = max(FemaleDeaths$AgeID_New),
                 "R" = max(FemaleDeaths$RegionID), 
                 "y" = FemaleDeaths$Deaths, "E" = FemaleDeaths$Population, 
                 "M" = 5, "TFor" = 33,
                 "TInd" = FemaleDeaths$YearID, 
                 "RInd" = FemaleDeaths$RegionID,
                 "AInd" = FemaleDeaths$AgeID_New, 
                 "CInd" = FemaleDeaths$CohortID,
                 "N_edges" = Nodes$N_edges,
                 "node1" = Nodes$node1, "node2" = Nodes$node2,
                 "scaling_factor" = scalingFac)

#run model (attention, may take a while)
RH_BYM2_Female <- rstan::stan(file = file.path(getwd(),"StanCode/RH_BYM2.stan"),
                              data = StanData, 
                              chains = 4, iter = 4000, warmup = 2000, thin = 2, 
                              control = list("max_treedepth" = 12,
                                             "adapt_delta" = 0.82))

#save results 
save(RH_BYM2_Female, file = file.path(getwd(),"Results/MortRate_F_RH_BYM2.RData"))

# extract forecasts and transform results into array for population projection##
#Note, Results are not found in Git Repo since they are too big. 
load(file = file.path(getwd(),"Results/MortRate_F_RH_BYM2.RData"))

Mort_Rate_FC_F_RH <- rstan::extract(RH_BYM2_Female, pars ="mufor")$mufor

MortRateHelperData <- expand.grid("AgeGroup" = unique(PopCounts$AgeGroup_New),
                                  "Region" = unique(PopCounts$RegionNumber),
                                  "Year" = 2018:2050)

#Transform Mortality Results into Array for forecasting
Mort_Rate_Array_F <- 
  Mort_Rate_FC_F_RH %>% t() %>% as.data.frame() %>% 
  mutate("Year" = MortRateHelperData$Year, #add information
         "RegionNumber" = MortRateHelperData$Region, 
         "AgeGroup" = MortRateHelperData$AgeGroup, 
         .before = 1) %>% 
  pivot_longer(cols = 4:ncol(.), names_to ="Iteration", #into long format
               values_to = "mu_for") %>% 
  mutate("mu_for" = exp(mu_for)) %>% #exponentiate to get rates
  mutate(Iteration = as.numeric(gsub("V", "", Iteration))) %>% #extract the numbers
  reshape2::acast(formula =  AgeGroup ~ Year ~ RegionNumber ~ Iteration, 
                  value.var = "mu_for")


#### 02_Males ##################################################################
MaleDeaths <- DeathData %>% 
  filter(Sex == "male") 


StanDataM <- list("T" = max(MaleDeaths$YearID),
                 "A" = max(MaleDeaths$AgeID_New),
                 "R" = max(MaleDeaths$RegionID), 
                 "y" = MaleDeaths$Deaths, "E" = MaleDeaths$Population, 
                 "M" = 5, "TFor" = 33,
                 "TInd" = MaleDeaths$YearID, 
                 "RInd" = MaleDeaths$RegionID,
                 "AInd" = MaleDeaths$AgeID_New, 
                 "CInd" = MaleDeaths$CohortID,
                 "N_edges" = Nodes$N_edges,
                 "node1" = Nodes$node1, "node2" = Nodes$node2,
                 "scaling_factor" = 0.4366992)

RH_BYM2_Male <- rstan::stan(file = file.path(getwd(),"StanCode/RH_BYM2.stan"),
                              data = StanDataM, 
                              chains = 4, iter = 4000, warmup = 2000, thin = 2, 
                            control = list("max_treedepth" = 13,
                                           "adapt_delta" = 0.82))

save(RH_BYM2_Male, file = file.path(getwd(),"Results/MortRate_M_RH_BYM2.RData"))

## Extract forecasts and transform into array for population projection ########
load(file = file.path(getwd(),"Results/MortRate_M_RH_BYM2.RData"))

Mort_Rate_FC_M_RH <- rstan::extract(RH_BYM2_Male, pars ="mufor")$mufor

MortRateHelperData <- expand.grid("AgeGroup" = unique(PopCounts$AgeGroup_New),
                                  "Region" = unique(PopCounts$RegionNumber),
                                  "Year" = 2018:2050)

Mort_Rate_Array_M <- 
  Mort_Rate_FC_M_RH %>% t() %>% as.data.frame() %>% 
  mutate("Year" = MortRateHelperData$Year,#add information
         "RegionNumber" = MortRateHelperData$Region, 
         "AgeGroup" = MortRateHelperData$AgeGroup, 
         .before = 1) %>% 
  pivot_longer(cols = 4:ncol(.), names_to ="Iteration", 
               values_to = "mu_for") %>% 
  mutate("mu_for" = exp(mu_for)) %>% 
  mutate(Iteration = as.numeric(gsub("V", "", Iteration))) %>% #extract the numbers
  reshape2::acast(formula = AgeGroup ~ Year ~ RegionNumber ~ Iteration, 
                  value.var = "mu_for")

save(Mort_Rate_Array_M, 
     Mort_Rate_Array_F, 
     file = file.path(getwd(),"Results/MortalityForecasts.RData"))


