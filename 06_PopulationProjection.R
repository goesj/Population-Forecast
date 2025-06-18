### Population Total ###
pacman::p_load("tidyverse","openxlsx","rstan","reshape2",
               "ggdist","lme4")


options(mc.cores = parallel::detectCores())
source("02_Functions.R")

load(file = file.path(getwd(),"Data/TotalData.RData"))

age_groups_new <- c(paste0(0:16*5, "-", 1:17*5-1), "85+")
AgesUpper <- 0:17*5 
## Combine all of migration, fertility and mortality data to obtain population estimates##

################################################################################
####              --- NET MIGRATION TOTALS ---                            ######
################################################################################
load(file = file.path(getwd(),"Results/MortalityForecasts.RData")) #Mortality Rates
load(file = file.path(getwd(),"Results/FertilityForecasts.RData")) #Mortality Rates
load(file = file.path(getwd(),"Results/MigrationForecasts_Counts.RData")) # Migration Rates
load(file = file.path(getwd(),"Results/Ra_Samples_Totals.RData")) #Estimates of Age-Migration Schedule

#### POPULATION PROJECTION with YEARLY migration and yearly Ra #################
InMig_Count_Array_M <- InMig_Count_Array_F <- OutMig_Count_Array_F <- 
  OutMig_Count_Array_M <- array(0, 
                                dim = dim(NetMig_Count_For_F),
                                dimnames = dimnames(NetMig_Count_For_F))


Iter <- 4000
n_AgeGroups <- length(unique(PopCounts$AgeGroup_New))
H <- length(seq(2024,2044,5))
R <- 13

PopTrajectories_Male <- PopTrajectories_Female <- 
  array(data = 0, dim = c(n_AgeGroups,  H, 13, Iter),
        dimnames=list("age_group" = age_groups_new,
                      "Year" = seq(2024,2044,5),
                      "RegionNumber" = unique(PopCounts$RegionNumber),
                      "Iteration" = 1:Iter))

#Natural Population (w/o Migration)
PopTrajectories_Male_Nat <- 
  PopTrajectories_Female_Nat <- 
  array(data = 0, dim = c(n_AgeGroups,  H, 13, Iter),
        dimnames=list("age_group" = age_groups_new,
                      "Year" = seq(2024,2044,5),
                      "RegionNumber" = unique(PopCounts$RegionNumber),
                      "Iteration" = 1:Iter))


## Last Observed Population Data
Pop_Array_2024 <- PopCounts %>% 
  mutate(AgeGroup_New = factor(AgeGroup_New, levels = unique(AgeGroup_New))) %>% 
  filter(Year == 2024) %>% 
  reshape2::acast(formula = AgeGroup_New ~ RegionNumber~Sex, 
                  value.var = "Population")

# Add population of 2023 to Pop Trajectories array (as starting population for Leslie Matrix)
PopTrajectories_Female[,1,,] <- Pop_Array_2024[,,1] #Females
PopTrajectories_Male[,1,,] <- Pop_Array_2024[,,2] #Males

PopTrajectories_Female_Nat[,1,,] <- Pop_Array_2024[,,1] #Females
PopTrajectories_Male_Nat[,1,,] <- Pop_Array_2024[,,2] #Males


UniqueRegions <- unique(PopCounts$RegionNumber)

# Empty Array for age specific migration rates
AS_Mig_Counts_F <- AS_Mig_Counts_M <- 
  array(0, dim = c(n_AgeGroups, H-1, 13, Iter, 3),
        dimnames = list(age_groups_new, 
                        seq(2029,2044,5),
                        unique(PopCounts$RegionNumber), 
                        1:Iter, 
                        c("InMig", "OutMig", "NetMig"))) #one year less

### Calculate mixed effects model to transform net migration into in and out migraion
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

## Run mixed effects model
model_female_counts <- lmer(InMig ~ NetMig  + (1|RegionNumber) -1 ,
                            data = Female_Mig_Data)

model_male_counts <- lmer(InMig ~ NetMig + (1|RegionNumber) -1 ,
                          data = Male_Mig_Data)

Coef_F <- coef(model_female_counts)$RegionNumber
Coef_M <- coef(model_male_counts)$RegionNumber


pb <- txtProgressBar(min = 1, max = Iter, style = 3)

for(s in 1:Iter){ #all iterations
  setTxtProgressBar(pb, s)
  for(h in 1:(H-1)){ #entire FC Period
    j <- h*5 #helper for data frames with yearly values (e.g. Mort Rates, Fertility Rates)
    for(r in 1:R){
      
      ##### 1.) Calculate Migration Counts 
      
      ### 1.1 Females 
      #In Migration for the next five years
      InMig_Count_Array_F[(j-3):(j+1), r, s] <-  
        Coef_F[r,1]+ #scalar is added to each element of vector
        Coef_F[r,2]*NetMig_Count_For_F[(j-3):(j+1),r,s]
      
      OutMig_Count_Array_F[(j-3):(j+1), r, s] <- 
        InMig_Count_Array_F[(j-3):(j+1), r, s] - 
        NetMig_Count_For_F[(j-3):(j+1),r,s] #since Net MIg of 5 years
      
      
      #In Migration Counts
      #elementwise multiplication. Each element of the vector is multiplied by each row of the matrix
      AS_Mig_Counts_F[, h, r, s, 1] <- colSums(InMig_Count_Array_F[(j-3):(j+1), r, s] * 
                                                 Ra_In_F[s, r, sample(1:13, size = 5, replace = TRUE) , ])
      
      
      
      #Out Migration Counts
      AS_Mig_Counts_F[, h, r, s, 2] <- colSums(OutMig_Count_Array_F[(j-3):(j+1), r, s] * 
                                                 Ra_Out_F[s, r, sample(1:13, size = 5, replace = TRUE) , ])
      
      # Net Migtaion Count
      AS_Mig_Counts_F[, h, r, s, 3] <- AS_Mig_Counts_F[, h, r, s, 1] - 
        AS_Mig_Counts_F[, h, r, s, 2]
      
      
      ### 1.2 Males ###
      #Migration for the next five years
      InMig_Count_Array_M[(j-3):(j+1), r, s] <-    
        Coef_M[r,1]+
        Coef_M[r,2]*NetMig_Count_For_M[(j-3):(j+1),r,s]
      
      OutMig_Count_Array_M[(j-3):(j+1), r, s] <- InMig_Count_Array_M[(j-3):(j+1), r, s] - 
        NetMig_Count_For_M[(j-3):(j+1),r,s]
      
      #In Migration Counts
      #elementwise multiplikation. Each element of the vector is multiplied by each row of the matrix
      AS_Mig_Counts_M[, h, r, s, 1] <- colSums(InMig_Count_Array_M[(j-3):(j+1), r, s] * 
                                                 Ra_In_M[s, r, sample(1:13, size = 5, replace = TRUE), ]) #sample years randomly 
      
      
      #Out Migration Counts
      AS_Mig_Counts_M[, h, r, s, 2] <- colSums(OutMig_Count_Array_M[(j-3):(j+1), r, s] * 
                                                 Ra_Out_M[s, r, sample(1:13, size = 5, replace = TRUE), ])
      
      AS_Mig_Counts_M[, h, r, s, 3] <- AS_Mig_Counts_M[, h, r, s, 1] - 
        AS_Mig_Counts_M[, h, r, s, 2]
      
      
      
      
      ##### 2. ) Forecast the Population 
      
      #### 2.1 Females 
      LifeTable_F <- LifeTableFun(Age = AgesUpper,
                                  MortalityRate = Mort_Rate_Array_F[,j + 7, r, s], 
                                  radix = 1)
      #1.1.2) Get Fertility Rate 
      FertRates <- c(0,0,0, Fert_Rate_Array_Direct[, j+1, r, s], rep( 0, 9))
      
      
      #2.1 Population Forecast Females
      PopTrajectories_Female[,h + 1, r, s] <- 
        (Leslie(L = LifeTable_F[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_F[1,"lx"], 
                sex ="Female")%*%(PopTrajectories_Female[,h, r, s]+
                                    AS_Mig_Counts_F[,h, r, s, 3]/2) #half of net migration of age group before
        )+ 
        AS_Mig_Counts_F[,h, r, s, 3]/2 #half at the end of the period (current age group)
      
      #### 2.2 Males
      LifeTable_M <- LifeTableFun(Age = AgesUpper,
                                  MortalityRate = Mort_Rate_Array_M[, j + 7 , r, s], 
                                  radix = 1)
      
      ## 2.2.1 Population Forecast males
      PopTrajectories_Male[,h + 1, r, s] <- 
        (Leslie(L = LifeTable_M[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_M[1,"lx"], 
                sex ="male")%*%(PopTrajectories_Male[,h, r, s]+
                                  AS_Mig_Counts_M[, h, r, s, 3]/2)+
           AS_Mig_Counts_M[, h, r, s, 3]/2)
      
      BM <- (Leslie(L=LifeTable_F[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_F[1,"lx"], sex="Male")%*%
               (PopTrajectories_Female[,h, r, s])+ AS_Mig_Counts_F[,h, r, s, 3]/2) %>% first() #Male Births
      
      PopTrajectories_Male[1,h+1, r, s] <- BM + AS_Mig_Counts_M[1, h, r, s,3]/2
      
      
      
      #### Natural Population ###
      PopTrajectories_Female_Nat[,h + 1, r, s] <- 
        Leslie(L = LifeTable_F[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_F[1,"lx"], 
               sex ="Female")%*%PopTrajectories_Female_Nat[,h, r, s]
      
      PopTrajectories_Male_Nat[,h + 1, r, s] <- 
        Leslie(L = LifeTable_M[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_M[1,"lx"], 
               sex ="male")%*%PopTrajectories_Male_Nat[,h, r, s]
      
      PopTrajectories_Male_Nat[1,h+1, r, s] <- 
        (Leslie(L=LifeTable_F[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_F[1,"lx"], sex="Male")%*%
           PopTrajectories_Female_Nat[,h, r, s]) %>% first() #Male Births
      
    }
  }
}

save(PopTrajectories_Female,
     PopTrajectories_Male,
     PopTrajectories_Female_Nat,
     PopTrajectories_Male_Nat, #Natural Population
     file = file.path(getwd(),"Results/PopTrajectories_RaTotal_24.RData"))

PopTrajectories_Total <- PopTrajectories_Female + PopTrajectories_Male

################################################################################
#### OOS EVALUATION ############################################################
#### Obs of Mortality Rates till 2019, Forecast 2020 - 2024 ##################### 
load(file = file.path(getwd(),"Results/MortalityForecasts.RData")) #Mortality Rates
load(file = file.path(getwd(),"Results/FertilityForecasts_20_24.RData"))
load(file = file.path(getwd(),"Results/MigrationForecasts_Counts_20_24.RData")) # Migration Rates
load(file = file.path(getwd(),"Results/Ra_Samples_Totals.RData")) #Estimates of Age-Migration Schedule


InMig_Count_Array_M <- InMig_Count_Array_F <- OutMig_Count_Array_F <- 
  OutMig_Count_Array_M <- array(0, 
                                dim = dim(NetMig_Count_For_F_20_24),
                                dimnames = dimnames(NetMig_Count_For_F_20_24))

Iter <- 4000
n_AgeGroups <- length(unique(PopCounts$AgeGroup_New))
H <- length(seq(2019,2024,5))
R <- 13

PopTrajectories_Male_OOS <- 
  PopTrajectories_Female_OOS <- 
  array(data = 0, dim = c(n_AgeGroups,  H, 13, Iter),
        dimnames=list("age_group" = age_groups_new,
                      "Year" = seq(2019,2024,5),
                      "RegionNumber" = unique(PopCounts$RegionNumber),
                      "Iteration" = 1:Iter))

## Last Observed Population Data
Pop_Array_2019 <- PopCounts %>% 
  mutate(AgeGroup_New = factor(AgeGroup_New, levels = unique(AgeGroup_New))) %>% 
  filter(Year == 2019) %>% 
  reshape2::acast(formula = AgeGroup_New ~ RegionNumber~Sex, 
                  value.var = "Population")

# Add population of 2023 to Pop Trajectories array (as starting population for Leslie Matrix)
PopTrajectories_Female_OOS[,1,,] <- Pop_Array_2019[,,1] #Females
PopTrajectories_Male_OOS[,1,,] <- Pop_Array_2019[,,2] #Males

UniqueRegions <- unique(PopCounts$RegionNumber)

# Empty Array for age specific migration rates
AS_Mig_Counts_F <- AS_Mig_Counts_M <- 
  array(0, dim = c(n_AgeGroups, H-1, 13, Iter, 3)) #one year less


pb <- txtProgressBar(min = 1, max = Iter, style = 3)

model_female_counts <- lmer(InMig ~ NetMig  + (1|RegionNumber) -1 ,
                            data = Female_Mig_Data_00_19)

model_male_counts <- lmer(InMig ~ NetMig + (1|RegionNumber) -1 ,
                          data = Male_Mig_Data_00_19)


Coef_F <- coef(model_female_counts)$RegionNumber
Coef_M <- coef(model_male_counts)$RegionNumber

for(s in 1:Iter){ #all iterations
  setTxtProgressBar(pb, s)
  for(h in 1:(H-1)){ #entire FC Period
    j <- h*5 #helper for data frames with yearly values (e.g. Mort Rates, Fertility Rates)
    for(r in 1:R){
      
      ##### 1.) Calculate Migration Counts 
      
      ### 1.1 Females
      InMig_Count_Array_F[(j-4):j, r, s] <-  Coef_F[r,1]+
        Coef_F[r,2]*NetMig_Count_For_F_20_24[(j-4):j,r,s]
      
      OutMig_Count_Array_F[(j-4):j, r, s] <- 
        InMig_Count_Array_F[(j-4):j, r, s] - 
        NetMig_Count_For_F_20_24[(j-4):j,r,s] #since Net MIg of 5 years
      
      
      #In Migration Counts
      #elementwise multiplikation. Each element of the vector is multiplied by each row of the matrix
      AS_Mig_Counts_F[, h, r, s, 1] <- colSums(InMig_Count_Array_F[(j-4):j, r, s] * 
                                                 Ra_In_F[s, r, sample(1:13, size = 5, replace = TRUE) , ])
      
      
      
      #Out Migration Counts
      AS_Mig_Counts_F[, h, r, s, 2] <- colSums(OutMig_Count_Array_F[(j-4):j, r, s] * 
                                                 Ra_Out_F[s, r, sample(1:13, size = 5, replace = TRUE) , ])
      
      # Net Migtaion Count
      AS_Mig_Counts_F[, h, r, s, 3] <- AS_Mig_Counts_F[, h, r, s, 1] - 
        AS_Mig_Counts_F[, h, r, s, 2]
      
      
      ### 1.2 Males ###
      #Migration for the next five years
      InMig_Count_Array_M[(j-4):j, r, s] <-   Coef_M[r,1]+
        Coef_M[r,2]*NetMig_Count_For_M_20_24[(j-4):j,r,s]
      
      OutMig_Count_Array_M[(j-4):j, r, s] <- InMig_Count_Array_M[(j-4):j, r, s] - 
        NetMig_Count_For_M[(j-4):j,r,s]
      
      #In Migration Counts
      #elementwise multiplikation. Each element of the vector is multiplied by each row of the matrix
      AS_Mig_Counts_M[, h, r, s, 1] <- colSums(InMig_Count_Array_M[(j-4):j, r, s] * 
                                                 Ra_In_M[s, r, sample(1:13, size = 5, replace = TRUE), ]) #sample years randomly 
      
      
      #Out Migration Counts
      AS_Mig_Counts_M[, h, r, s, 2] <- colSums(OutMig_Count_Array_M[(j-4):j, r, s] * 
                                                 Ra_Out_M[s, r, sample(1:13, size = 5, replace = TRUE), ])
      
      AS_Mig_Counts_M[, h, r, s, 3] <- AS_Mig_Counts_M[, h, r, s, 1] - 
        AS_Mig_Counts_M[, h, r, s, 2]
      
      
      
      ##### 2. ) Forecast the Population 
      
      #### 2.1 Females 
      LifeTable_F <- LifeTableFun(Age = AgesUpper,
                                  MortalityRate = Mort_Rate_Array_F[,j + 2, r, s], 
                                  radix = 1)
      #1.1.2) Get Fertility Rate 
      FertRates <- c(0,0,0, Fert_Rate_Array_Direct[, j, r, s], rep( 0, 9))
      
      
      #2.1 Population Forecast Females
      PopTrajectories_Female_OOS[,h + 1, r, s] <- 
        (Leslie(L = LifeTable_F[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_F[1,"lx"], 
                sex ="Female")%*%(PopTrajectories_Female_OOS[,h, r, s]+
                                    AS_Mig_Counts_F[,h, r, s, 3]/2) #half of net migration of age group before
        )+ 
        AS_Mig_Counts_F[,h, r, s, 3]/2 #half at the end of the period (current age group)
      
      #### 2.2 Males
      LifeTable_M <- LifeTableFun(Age = AgesUpper,
                                  MortalityRate = Mort_Rate_Array_M[, j + 2 , r, s], 
                                  radix = 1)
      
      ## 2.2.1 Population Forecast males
      PopTrajectories_Male_OOS[,h + 1, r, s] <- 
        (Leslie(L = LifeTable_M[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_M[1,"lx"], 
                sex ="male")%*%(PopTrajectories_Male_OOS[,h, r, s]+
                                  AS_Mig_Counts_M[, h, r, s, 3]/2)+
           AS_Mig_Counts_M[, h, r, s, 3]/2)
      
      BM <- (Leslie(L=LifeTable_F[,"Lx"], f = FertRates, SRB = SRB, l0 = LifeTable_F[1,"lx"], sex="Male")%*%
               (PopTrajectories_Female_OOS[,h, r, s])+ AS_Mig_Counts_F[,h, r, s, 3]/2) %>% first() #Male Births
      
      PopTrajectories_Male_OOS[1,h+1, r, s] <- BM + AS_Mig_Counts_M[1, h, r, s,3]/2
    }
  }
}

save(PopTrajectories_Female_OOS, 
     PopTrajectories_Male_OOS, 
     file = file.path(getwd(),"Results/PopTrajectories_RaTotal_OOSEval.RData"))

load("Results/PopTrajectories_RaTotal_OOSEval.RData")

PopTrajectories_Total_OOS <- PopTrajectories_Female_OOS + 
  PopTrajectories_Male_OOS

### OOS Evaluation ####
Pop24 <- filter(PopCounts, Year == 2024) %>% 
  reframe(Population = sum(Population), #sum by sex
          .by = c(Year, RegionNumber, AgeGroup_New)) %>% 
  rename(age_group = AgeGroup_New)

## Evaluation
ResTable <- 
  PopTrajectories_Total_OOS %>% 
  reshape2::melt(value.name = "PopFC", 
                 as.is = TRUE) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  reframe(mean_qi(PopFC, .width = 0.95), 
          .by = c(age_group, Year, RegionNumber)) %>% 
  filter(Year == 2024) %>% 
  left_join(x = ., 
            y = Pop24) %>% 
  mutate("et" = y-Population, #absolut error
         "pt" = 100*et/Population) %>% #percentage error 
  reframe("MAE" = sqrt(mean(et^2)),
          "RMSE" = mean(abs(et)),
          "MAPE" = mean(abs(pt)), #see Hyndman, Koehler
          "RMSPE" = sqrt(mean(pt^2)),
          "Cov" = mean(dplyr::between(Population, ymin, ymax)))


