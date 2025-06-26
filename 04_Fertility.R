#### SCRIPT TO RUN Fertility MODEL ############################################
pacman::p_load("tidyverse","openxlsx","rstan","reshape2")


options(mc.cores = parallel::detectCores())
source("02_Functions.R")

load(file = file.path(getwd(),"Data/TotalData.RData"))

##### Getting Data #############################################################
BirthData <- 
  bind_rows(PopCounts,
            Population_Old_Fem) %>% #incl. Population Counts from 1995-1999
  filter(Sex =="female") %>% 
  inner_join(x = BirthCounts, 
             y = . , 
             by = c("Year","RegionNumber","AgeGroup_New")) %>% 
  select(-c(RegionName.y)) %>% #remove due to different naming structure
  rename(RegionName = RegionName.x) %>% 
  select(-c(AgeID_New, AgeGroup, Sex)) %>% 
  mutate("ASFR" = Births / Population) %>% #calculate asfr
  mutate("n_ASFR" = 5*ASFR) #for calculation of TFR

#Save data frame of TFR separately
TFR_Frame <- BirthData %>%
  reframe("TFR" = sum(n_ASFR), #calculate TFR
          .by = c("RegionNumber", "Year"))

#Calculate Age Proportion
BirthData <- 
  left_join(x = BirthData, 
            y = TFR_Frame, by =c("RegionNumber", "Year")) %>% #add TFR
  mutate("Prop_Age" = n_ASFR/TFR) %>% 
  arrange(Year, RegionNumber) #put into correct order

##### Running the LC type model ################################################
StanDat <- list("T" = length(unique(BirthData$Year)),
                "A" = length(unique(BirthData$AgeGroup_New)),
                "R" = length(unique(BirthData$RegionNumber)),
                "y" = BirthData$Births, 
                "E" = BirthData$Exposure,
                "M" = 5, TFor = 21) # Forecast until 2044

# run model
LC_Model <- rstan::stan(file = file.path(getwd(),"StanCode/LC_Birth_Poisson.stan"),
                        data =  StanDat, chains = 4, 
                        iter = 5000, warmup = 3000, 
                        thin = 2)
#save results
save(LC_Model, file = file.path(getwd(),"Results/ASFR_LC_Poi.RData"))

## Extract forecasts and transform into array for population projection ########
Fert_Rate_FC_F <- rstan::extract(LC_Model, pars ="mufor")$mufor

FertRateHelper <- expand.grid("AgeGroup" = unique(BirthCounts$AgeGroup_New),
                              "Region" = unique(BirthCounts$RegionNumber),
                              "Year" = 2024:2050)

Fert_Rate_Array_Direct <- 
  Fert_Rate_FC_F %>% t() %>% as.data.frame() %>% 
  mutate("Year" = FertRateHelper$Year,#add information
         "RegionNumber" = FertRateHelper$Region, 
         "AgeGroup" = FertRateHelper$AgeGroup, 
         .before = 1) %>% 
  pivot_longer(cols = 4:ncol(.), names_to ="Iteration", 
               values_to = "mu_for") %>% 
  mutate("mu_for" = exp(mu_for)) %>% 
  mutate(Iteration = as.numeric(gsub("V", "", Iteration))) %>% #extract the numbers
  reshape2::acast(formula = AgeGroup ~ Year ~ RegionNumber ~ Iteration,
                  value.var = "mu_for")

save(Fert_Rate_Array_Direct, 
     file = file.path(getwd(),"Results/FertilityForecasts.RData"))


####### 2 Out-of-sample evaluation ############################################
#### 2.1 Direct Approach #######################################################
BirthData_Train <- BirthData %>%
  filter(Year < 2015)

StanDat <- list("T" = length(unique(BirthData_Train$Year)),
                "A" = length(unique(BirthData_Train$AgeGroup)),
                "R" = length(unique(BirthData_Train$RegionNumber)),
                "y" = BirthData_Train$Births, 
                "E" = BirthData_Train$Population,
                "M" = 5, TFor = 8) #forecast till 2023

### 3.1.3 Poisson Model 
LC_Poi_95_15 <- rstan::stan(file = file.path(getwd(),"Stan/LC_BirthPoisson.stan"),
                            data =  StanDat, chains = 4, 
                            iter = 5000, warmup = 3000, 
                            thin = 2)

save(LC_Poi_95_15, file = file.path(getwd(),"Results/ASFR_LC_Poi_95_15.RData"))

#### 2.2 Indirect Approach ######################################################
### 2.2.1 Dirichlet Regression Prop Age 
DataArray_Train <- 
  BirthData_Train %>% 
  select(RegionNumber, Year, AgeGroup_New, Prop_Age) %>%
  reshape2::acast(formula = RegionNumber~Year~AgeGroup_New, 
                  value.var = "Prop_Age")  


#Data as Array, rows = Regions, Cols = Year, thrid dim = AgeGroup
StanDat_Dir_R <- list("T" = length(unique(OF_Data_Train$Year)), 
                      "A" = length(unique(OF_Data_Train$AgeGroup)),
                      "y" = DataArray_Train,
                      "X_t" = 1:dim(DataArray_Train)[2],
                      "R" = length(unique(OF_Data$RegionNumber)),
                      "H" = 8) #forecast till 2023


Dir_Reg_95_15 <- 
  rstan::stan(file = file.path(getwd(),"Stan/DirichletRegression_PropAge.stan"),
              data =  StanDat_Dir_R, chains = 4, 
              iter = 4000, warmup = 2000, 
              thin = 2, control = list("max_treedepth" = 12))

## Save results
save(Dir_Reg_Model_AgeRegionInt_95_15, 
     file = file.path(getwd(),"Results/Samples_Dirichlet_Reg_95_15.RData"))

#### 2.2.2 TFR Forecasts
TFR_Mat <- TFR_Frame %>% 
  filter(Year < 2016) %>% 
  select(RegionNumber, Year, TFR) %>% 
  reshape2::acast(formula = RegionNumber~Year, 
                  value.var = "TFR") 

StanList <- list("T" = ncol(TFR_Mat),
                 "R" = nrow(TFR_Mat), 
                 "y" = TFR_Mat,
                 "H" = 8)  #forecast till 2023

TFR_AR1_95_15 <- rstan::stan(file = file.path(getwd(),"Stan/AR1_Hier.stan"),
                             data =  StanList, chains = 4, 
                             iter = 5000, warmup = 3000, thin = 2)

save(TFR_AR1_95_15, file = file.path(getwd(),"Results/Samples_TFR_AR1_95_15.RData"))



