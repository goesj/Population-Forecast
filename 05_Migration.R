pacman::p_load("tidyverse","openxlsx","rstan",
               "lme4","MuMIn","ggdist")


options(mc.cores = parallel::detectCores())
source("02_Functions.R")

load(file = file.path(getwd(),"Data/TotalData.RData"))
age_groups_new <- c(paste0(0:16*5, "-", 1:17*5-1), "85+")


#### Calculate net migration rates ############################################

#First calculate total population counts 
PopCounts_Total <- 
  PopCounts %>% 
  reframe("Population" = sum(Population), 
          .by = c(Year, RegionNumber, RegionName, Sex))

### Combine Both Datasets ###
Female_Mig_Data <- 
  PopCounts_Total %>% 
  filter(Sex == "female") %>% 
  inner_join(x = ., 
             y =  select(MigCounts_Female_Tot, -RegionName),
             by = c("RegionNumber", "Year")) %>%
  mutate(PopCountsRisk = Population - NetMig) %>% 
  mutate(NetMigRate = NetMig / PopCountsRisk,
         InMigRate = InMig / PopCountsRisk,
         OutMigRate = OutMig / PopCountsRisk) %>% 
  mutate(RegionID = match(RegionNumber, unique(RegionNumber)),
         YearID = match(Year, unique(Year)))

Male_Mig_Data <- 
  PopCounts_Total %>% 
  filter(Sex == "male") %>% 
  inner_join(x = ., 
             y =  select(MigCounts_Male_Tot, -RegionName),
             by = c("RegionNumber", "Year")) %>% 
  mutate(PopCountsRisk = Population - NetMig) %>% 
  mutate(NetMigRate = NetMig / PopCountsRisk,
         InMigRate = InMig / PopCountsRisk,
         OutMigRate = OutMig / PopCountsRisk) %>% 
  mutate(RegionID = match(RegionNumber, unique(RegionNumber)),
         YearID = match(Year, unique(Year)))


#### Relationship Net migration rate to in-migration rate ###################
model_female_counts <- lmer(InMig ~ NetMig  + (1|RegionNumber) -1 ,
                          data = Female_Mig_Data)

model_male_counts <- lmer(InMig ~ NetMig + (1|RegionNumber) -1 ,
                          data = Male_Mig_Data)

###### 02_ Female  #############################################################
#### Forecasting Counts ####
NetMigMatF_Counts <- reshape2::acast(Female_Mig_Data, 
                                     formula = RegionID ~ Year, 
                                     value.var = "NetMig")

StanDat_F_Counts <- list("R" = nrow(NetMigMatF_Counts),
                  "T" = ncol(NetMigMatF_Counts), 
                  "y" = NetMigMatF_Counts, 
                  "H" = 27)


AzRa_F_Skew_N_Counts <- rstan::stan(file = file.path(getwd(),"StanCode/AzozeRaftery_SkewNormal_Total.stan"), 
                                    chains = 4, iter = 3000, warmup = 2000, 
                                    data = StanDat_F_Counts)



###### 03_ Male ################################################################
#### Forecasting Counts #########
NetMigMatM_Counts <- reshape2::acast(Male_Mig_Data, 
                                     formula = RegionID ~ Year, 
                                     value.var = "NetMig")

StanDat_M_Counts <- list("R" = nrow(NetMigMatM_Counts),
                          "T" = ncol(NetMigMatM_Counts), 
                          "y" = NetMigMatM_Counts, 
                          "H" = 27)


AzRa_M_Skew_N_Counts <- rstan::stan(file = file.path(getwd(),"StanCode/AzozeRaftery_SkewNormal_Total.stan"), 
                                    chains = 4, iter = 3000, warmup = 2000, 
                                    data = StanDat_M_Counts)

#Hierarchical sigma
save(AzRa_F_Skew_N_Counts,
     AzRa_M_Skew_N_Counts,
     file = file.path(getwd(),"Results/AzRa_Skew_Counts.RData"))


######## 04_Calculate Out Migration from relationship ##########################
load(file.path(getwd(),"Results/AzRa_Skew_Counts.RData"))


NetMig_Count_For_M <- rstan::extract(AzRa_M_Skew_N_Counts_H, pars ="y_for")$y_for %>% 
  aperm(perm = c(3,2,1))

NetMig_Count_For_F <- rstan::extract(AzRa_F_Skew_N_Counts_H, pars ="y_for")$y_for %>% 
  aperm(perm = c(3,2,1))

NetMig_Count_For_F_Long <- 
  NetMig_Count_For_F %>% 
  reshape2::melt(value.name = "NetMigCounts", varnames = c("Year","Region","Iteration")) %>% 
  reframe(ggdist::mean_qi(NetMigCounts), .by = c(Year,Region)) %>% 
  mutate("Year" = Year + 2023,
         "RegionNumber" = unique(Female_Mig_Data$RegionNumber)[Region]) %>% 
  rename(NetMig = y)

NetMig_Count_For_M_Long <- 
  NetMig_Count_For_M %>% 
  reshape2::melt(value.name = "NetMigCounts", varnames = c("Year","Region","Iteration")) %>% 
  reframe(ggdist::mean_qi(NetMigCounts), .by = c(Year,Region)) %>% 
  mutate("Year" = Year + 2023,
         "RegionNumber" = unique(Female_Mig_Data$RegionNumber)[Region]) %>% 
  rename(NetMig = y)

save(NetMig_Count_For_F_Long, 
     NetMig_Count_For_M_Long, 
     NetMig_Count_For_F,
     NetMig_Count_For_M,
     file = file.path(getwd(),"Results/MigrationForecasts_Counts.RData"))

####### 04_02 Forecast Migration Counts 20 - 24 ################################
Male_Mig_Data_00_19 <- 
  filter(Male_Mig_Data, Year < 2020)

NetMigMatM_Counts_00_19 <- reshape2::acast(Male_Mig_Data_00_19, 
                                           formula = RegionID ~ Year, 
                                           value.var = "NetMig")

StanDat_M_Counts_00_19 <- list("R" = nrow(NetMigMatM_Counts_00_19),
                               "T" = ncol(NetMigMatM_Counts_00_19), 
                               "y" = NetMigMatM_Counts_00_19, 
                               "H" = 5)


AzRa_M_Skew_N_Counts_00_19 <- rstan::stan(file = file.path(getwd(),"StanCode/AzozeRaftery_SkewNormal_Total.stan"), 
                                    chains = 4, iter = 3000, warmup = 2000, 
                                    data = StanDat_M_Counts_00_19)

#Females
Female_Mig_Data_00_19 <- 
  filter(Female_Mig_Data, Year < 2020)

NetMigMatF_Counts_00_19 <- reshape2::acast(Female_Mig_Data_00_19, 
                                     formula = RegionID ~ Year, 
                                     value.var = "NetMig")

StanDat_F_Counts_00_19 <- list("R" = nrow(NetMigMatF_Counts_00_19),
                               "T" = ncol(NetMigMatF_Counts_00_19), 
                               "y" = NetMigMatF_Counts_00_19, 
                               "H" = 5)

AzRa_F_Skew_N_Counts_00_19 <- 
  rstan::stan(file = file.path(getwd(),"StanCode/AzozeRaftery_SkewNormal_Total.stan"), 
              chains = 4, iter = 3000, warmup = 2000, 
              data = StanDat_F_Counts_00_19)

####### 05_ Age Migration Schedule #############################################
MigPop_Data_Female <- 
  inner_join(x = AgeSpecific_Mig_Female, 
             y = filter(PopCounts, Sex == "female"), 
             by = c("Year", "RegionNumber", "AgeID_New","AgeGroup_New", "RegionName")) %>% 
  mutate("PopCountsRisk" = Population - Net) %>% 
  mutate("OutMigRate" = Out/PopCountsRisk,
         "InMigRate" = In/PopCountsRisk) %>% 
  mutate("AgeGroup_New" = factor(AgeGroup_New, levels = unique(AgeGroup_New)))


MigPop_Data_Male <- 
  inner_join(x = AgeSpecific_Mig_Male, 
             y = filter(PopCounts, Sex == "male"),
             by = c("Year", "RegionNumber", "AgeID_New","AgeGroup_New", "RegionName")) %>% 
  mutate("PopCountsRisk" = Population - Net) %>% 
  mutate("OutMigRate" = Out/PopCountsRisk,
         "InMigRate" = In/PopCountsRisk) %>% 
  mutate("AgeGroup_New" = factor(AgeGroup_New, levels = unique(AgeGroup_New)))


############ Create data of Migration Schedules ########################
RaData_Male <- 
  MigPop_Data_Male %>% 
  mutate("Ra_Out" = Out/sum(Out),
         "Ra_In" = In/sum(In),
         .by = c(Year, RegionNumber)) 

#Matrix of out-migration Schedue
Ra_Out_Matrix_Male <- 
  RaData_Male %>% 
  reshape2::acast(formula = RegionNumber ~ Year ~ AgeID_New, value.var = "Ra_Out") #transform into array

#Problem, Ra cannot be zero since Dirichlet likelihood does not allow that. 
#Therefore add small value to Ra and renormalize
Ra_In_Matrix_Male <- 
  RaData_Male %>% 
  mutate(Ra_In = ifelse(Ra_In == 0, 1e-6, Ra_In)) %>% 
  reframe(Ra_In = Ra_In / sum(Ra_In), 
          AgeID_New = AgeID_New, 
          .by = c("Year","RegionNumber")) %>% 
  reshape2::acast(formula = RegionNumber ~ Year ~ AgeID_New, 
                  value.var = "Ra_In") #transform into array

### Running Dirichlet Regression to estimate migration Schedules ###############
RC_List_Male_Out <- 
  list("y" = Ra_Out_Matrix_Male, 
       "T" = dim(Ra_Out_Matrix_Male)[2],
       "R" = dim(Ra_Out_Matrix_Male)[1],
       "A" = dim(Ra_Out_Matrix_Male)[3])

RC_List_Male_In <- 
  list("y" = Ra_In_Matrix_Male, 
       "T" = dim(Ra_In_Matrix_Male)[2],
       "R" = dim(Ra_In_Matrix_Male)[1],
       "A" = dim(Ra_In_Matrix_Male)[3])

RogersCastroStan_Out_M <- 
  rstan::stan(file = file.path(getwd(),"StanCode/RogersCastro_Dirichlet.stan"),
              data = RC_List_Male_Out, chains = 4, iter = 4000, warmup = 2000, 
              thin = 2)

RogersCastroStan_In_M <- 
  rstan::stan(file = file.path(getwd(),"StanCode/RogersCastro_Dirichlet.stan"),
              data = RC_List_Male_In, chains = 4, iter = 4000, warmup = 2000, 
              thin = 2)

##### Same for females, first create data then run Dirichlet Regression ########
RaData_Female <- 
  MigPop_Data_Female %>% 
  mutate("Ra_Out" = Out/sum(Out),
         "Ra_In" = In/sum(In),
         .by = c(Year, RegionNumber)) 

Ra_Out_Matrix_Female <- 
  RaData_Female %>% 
  reshape2::acast(formula = RegionNumber ~ Year ~ AgeID_New, value.var = "Ra_Out") #transform into array

Ra_In_Matrix_Female <- 
  RaData_Female %>% 
  reshape2::acast(formula = RegionNumber ~ Year ~ AgeID_New, value.var = "Ra_In") #transform into array

RC_List_Female_Out <- 
  list("y" = Ra_Out_Matrix_Female, 
       "T" = dim(Ra_Out_Matrix_Female)[2],
       "R" = dim(Ra_Out_Matrix_Female)[1],
       "A" = dim(Ra_Out_Matrix_Female)[3])

RC_List_Female_In <- 
  list("y" = Ra_In_Matrix_Female, 
       "T" = dim(Ra_In_Matrix_Female)[2],
       "R" = dim(Ra_In_Matrix_Female)[1],
       "A" = dim(Ra_In_Matrix_Female)[3])

RogersCastroStan_Out_F <- 
  rstan::stan(file = file.path(getwd(),"StanCode/RogersCastro_Dirichlet.stan"),
              data = RC_List_Female_Out, chains = 4, iter = 4000, warmup = 2000, 
              thin = 2)

RogersCastroStan_In_F <- 
  rstan::stan(file = file.path(getwd(),"StanCode/RogersCastro_Dirichlet.stan"),
              data = RC_List_Female_In, chains = 4, iter = 4000, warmup = 2000, 
              thin = 2)

### Save Results
save(RogersCastroStan_Out_M,
     RogersCastroStan_In_M, 
     RogersCastroStan_Out_F,
     RogersCastroStan_In_F,
     file = file.path(getwd(),"Results/RogersCastro_Totals.RData"))


############# 5.2. Analysis of Results ########################################
load(file.path(getwd(),"Results/RogersCastro_Totals.RData"))

Ra_Out_M <- rstan::extract(RogersCastroStan_Out_M, pars ="Y_Rep")$Y_Rep
Ra_In_M <- rstan::extract(RogersCastroStan_In_M, pars ="Y_Rep")$Y_Rep

Ra_Out_F <- rstan::extract(RogersCastroStan_Out_F, pars ="Y_Rep")$Y_Rep
Ra_In_F <- rstan::extract(RogersCastroStan_In_F, pars = "Y_Rep")$Y_Rep

save(Ra_Out_M, Ra_In_F, Ra_Out_F, Ra_In_M, 
     file = file.path(getwd(),"Results/Ra_Samples_Totals.RData"))


