pacman::p_load("tidyverse","openxlsx")


options(mc.cores = parallel::detectCores())
### Loading all the dava and saving the results ###############################
#Helper function
age_groups_new <- c(paste0(0:16*5, "-", 1:17*5-1), "85+")

############## 01_Population Counts ############################################
PopulationM <- openxlsx::read.xlsx(xlsxFile = file.path(getwd(),
                                                        "Data/BY_12411(Population)_2000-17.xlsx"),
                                   sheet = "PopulationMale", colNames = TRUE)

PopulationW <- openxlsx::read.xlsx(xlsxFile = file.path(getwd(),
                                                        "Data/BY_12411(Population)_2000-17.xlsx"),
                                   sheet = "PopulationFemale", colNames = TRUE)

#Transform Data Into Long Format
PM <- pivot_longer(PopulationM, cols = -(1:5),
                   names_to = "AgeGroup", 
                   values_to="Population") %>% 
  mutate("Year" = as.numeric(Year)) %>% 
  mutate("Sex" = "male")

PW <- pivot_longer(PopulationW, cols = -c(1:5),
                   names_to = "AgeGroup",
                   values_to = "Population") %>% 
  mutate("Year" = as.numeric(Year)) %>% 
  mutate("Sex" = "female")



# calculate the sex ratio at birth
SRB <- full_join(PM, 
                 PW) %>% 
  filter(AgeGroup == "unter.1") %>% 
  filter(grepl("094",RegionNumber)) %>% 
  pivot_wider(., names_from = "Sex", values_from = "Population") %>% 
  reframe("SRB" = male/female , .by = c(Year, RegionNumber)) %>% 
  reframe("SRB" = mean(SRB)) %>% pull()



### Combine age groups together so they are a width of 5##################
# combine first two age groups and last three together
AgeID_Helper17 <- data.frame("AgeGroup" = unique(PW$AgeGroup),
                             "AgeID" = 1:length(unique(PW$AgeGroup)),
                             "AgeIDNew" = c(1,1,2:17,18,18,18))

PopCounts_17 <- 
  full_join(x = PW,
            y = PM) %>% 
  mutate("AgeID" = match(AgeGroup, unique(AgeGroup)),
         "AgeID_New" = AgeID_Helper17$AgeIDNew[AgeID]) %>% 
  reframe("Population" = sum(Population), 
          .by = c(Year, RegionNumber, RegionName, AgeID_New, Sex)) %>% 
  mutate("AgeGroup_New" = age_groups_new[AgeID_New])


## Loading in other Population Counts (different granularity)##################
PopCounts_24 <- openxlsx::read.xlsx(
  xlsxFile =  file.path(getwd(),"Data/BY_12411(Population)_2018-24.xlsx"),
  startRow = 5 
)

#Helper to create proper age groups as population is given for single ages
AgeIDHelper <- data.frame("AgeGroup" = unique(PopCounts_24$Age),
                          "AgeID" = 1:length(unique(PopCounts_24$Age)),
                          "AgeIDNew" = c(rep(1:17, each = 5), 18, 18, 18,19))
#Data Clensing
PopCounts_24 <- PopCounts_24 %>% 
  select(-c(Nummer, Kreis)) %>% 
  mutate(RegionName = stringr::str_trim(RegionName)) %>%  #Remove WhiteSpaces
  filter(nchar(RegionNumber) == 5) %>% #only regions
  filter(Age != "Insgesamt") %>% 
  select(-Total) %>% 
  mutate("AgeID" = match(Age, unique(Age))) %>% 
  mutate("AgeID_New" = AgeIDHelper$AgeIDNew[AgeID]) %>% 
  reframe("male" = sum(male), 
          "female" = sum(female), 
          .by = c(Year, RegionNumber, RegionName, AgeID_New)) %>% 
  pivot_longer(cols = 5:6, names_to = "Sex", values_to = "Population") %>% 
  mutate("AgeGroup_New" = age_groups_new[AgeID_New])


PopCounts_Total <- 
  full_join(x = PopCounts_17, 
            y = PopCounts_24) 

PopCounts <- 
  PopCounts_Total %>% 
   mutate("PopLag" = lag(Population),
          .by = c(RegionNumber, AgeGroup_New, Sex)) %>%
   #filter(Year > 2000) %>% 
   mutate("PopLag" = rnorm(n = n(), mean = PopLag, sd = 0.001)) %>% # add small noise, otherwise if pop stays the same division by 0 
   mutate("Exposure" = (Population - PopLag ) / log(Population/PopLag)) %>% 
   filter(grepl("094", RegionNumber)) %>% #only Upper Franconia
   select(-PopLag)
   
##### Population 1995 - 1999 ###################################################
# only needed for Birth Data
Population_Old_Fem <- 
  openxlsx::read.xlsx(xlsxFile = 
                        file.path(getwd(),"Data/BY_12411(Population)_1995-1999.xlsx"),
                                      sheet = "FemaleBavaria") %>% 
  mutate(Year = as.numeric(Year)) %>% 
  pivot_longer(cols = 4:ncol(.), names_to = "AgeGroup", values_to = "PopFemale") %>% 
  mutate("AgeID" = match(AgeGroup, unique(AgeGroup))) %>% 
  filter(AgeID > 4 & AgeID < 12) %>%  #only select women of reproductive age
  mutate("AgeID_Grouped" = ifelse(AgeID < 7, 1, AgeID -5)) %>%  #combine lowest two age groups
  reframe(PopFemale = sum(PopFemale), #sum over lowest two age groups
          .by = c(Year, RegionNumber, RegionName, AgeID_Grouped)) %>% 
  mutate("AgeGroup_New" = age_groups_new[AgeID_Grouped+3]) %>%  #create new age group
  mutate("Sex" = "female") %>% 
  rename(Population = PopFemale) 


################### 02_Deaths ##################################################
DeathM <- openxlsx::read.xlsx(xlsxFile = file.path(getwd(),
                                                   "Data/BY_12613(Deaths)_2000-17.xlsx"),
                              sheet = "SterbefälleMännlich", colNames = TRUE)

DeathW <- openxlsx::read.xlsx(xlsxFile = file.path(getwd(),
                                                   "Data/BY_12613(Deaths)_2000-17.xlsx"),
                              sheet = "SterbefälleWeiblich", colNames = TRUE)

DM <- DeathM %>% pivot_longer(., cols = -(1:5),
                              names_to = "AgeGroup", 
                              values_to="Deaths") %>% 
  mutate("Year" = as.numeric(Year)) %>% 
  mutate("Sex" = "male")

DF <- DeathW %>% pivot_longer(., cols = -(1:5), 
                              names_to = "AgeGroup", 
                              values_to="Deaths") %>% 
  mutate("Year" = as.numeric(Year)) %>% 
  mutate("Sex" = "female")


DeathCounts <- 
  full_join(x = DM, 
            y = DF) %>% 
  mutate("AgeID" = match(AgeGroup, unique(AgeGroup)),
         "AgeID_New" = AgeID_Helper17$AgeIDNew[AgeID]) %>% 
  reframe("Deaths" = sum(Deaths), 
          .by = c(Year, RegionNumber, RegionName, AgeID_New, Sex)) %>% 
  mutate("AgeGroup_New" = age_groups_new[AgeID_New]) %>% 
  filter(grepl("094", RegionNumber)) #only Upper Franconia


############ 03_Births #########################################################
Births <- openxlsx::read.xlsx(xlsxFile = 
                                file.path(getwd(),"Data/BY_12612(Births)_1983-2023.xlsx"),
                              startRow = 4, na.strings = ".") %>% 
  select(-c(Nummer, Name)) %>% 
  filter(nchar(RegionNumber) > 3) %>% #only select regions  
  arrange(RegionNumber) %>% #sort by region number
  mutate(Year = as.numeric(Year)) %>% 
  mutate(RegionName = trimws(RegionName)) %>% #removs whitespace characters
  pivot_longer(cols = 4:ncol(.), names_to = "AgeGroup",values_to = "Births") %>% 
  filter(AgeGroup != "Insgesamt")



#Change age group names of Births. 15-19 denotes under 20 and 45-45 denotes 40+
AgeGroupNames <- age_groups_new[4:9]

#change Age Groups
BirthCounts <- Births %>% 
  mutate("AgeID" = match(AgeGroup, unique(AgeGroup))) %>% 
  mutate("AgeGroup_New" = AgeGroupNames[AgeID]) %>% 
  filter(grepl("094", RegionNumber)) #only Upper Franconia


#Question, what to do with the NA's. NA's only below 3 birhts after 2021
Births %>% reframe("NA's" = sum(is.na(Births)))
# No missing value in Oberfranken
Births %>% filter(Year > 2000) %>% 
  filter(grepl("094",RegionNumber)) %>% 
  reframe("Prop_NA" = sum(is.na(Births))/n() * 100)


######## 04_Migration ##########################################################
AgeSpecific_MigData <- 
  openxlsx::read.xlsx(xlsxFile = 
                        file.path(getwd(),"Data/BY_NetMigration_2011-2023.xlsx"),
                      startRow = 10, na.strings = "NA")

#Helper function
In_Out_Selction <- function(x,sex ="männlich", type ="Zu"){
  out <- 
    x %>% select(c(1,2), matches({{sex}})) %>% 
    select(c(1,2),matches({{type}})) %>% 
    rename_with(~ str_replace(., paste0(type, sex), "")) 
  return(out)
}

## create respective male and female data of age-specific in/out migration
AgeSpecificMig_In_Male <- 
  In_Out_Selction(x = AgeSpecific_MigData, sex ="männlich", type ="Zu")  %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "AgeGroup", 
               values_to = "In") %>% 
  mutate(sex = "male", .before = 4)

AgeSpecificMig_Out_Male <- 
  In_Out_Selction(x = AgeSpecific_MigData, sex ="männlich", type ="Fo") %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "AgeGroup", 
               values_to = "Out") %>% 
  mutate(sex = "male", .before = 4)

AgeSpecific_Mig_Net_Male <- 
  In_Out_Selction(x = AgeSpecific_MigData, sex ="männlich", type ="Wa") %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "AgeGroup", 
               values_to = "Net") %>% 
  mutate(sex = "male", .before = 4)

AgeSpecificMig_In_Female <- 
  In_Out_Selction(x = AgeSpecific_MigData, sex ="weiblich", type ="Zu") %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "AgeGroup", 
               values_to = "In") %>% 
  mutate(sex = "female", .before = 4)

AgeSpecificMig_Out_Female <- 
  In_Out_Selction(x = AgeSpecific_MigData, sex ="weiblich", type ="Fo") %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "AgeGroup", 
               values_to = "Out") %>% 
  mutate(sex = "female", .before = 4)

AgeSpecific_Mig_Net_Female <- 
  In_Out_Selction(x = AgeSpecific_MigData, sex ="weiblich", type ="Wa") %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "AgeGroup", 
               values_to = "Net") %>% 
  mutate(sex = "female", .before = 4)

RegionNameHelper <- 
  data.frame("RegionNumber" = unique(PopCounts$RegionNumber),
             "RegionName" = unique(PopCounts$RegionName))


#combine both in and out migration and add net migration column
AgeSpecific_Mig_Male <- 
  left_join(x = AgeSpecificMig_In_Male,
            y = AgeSpecificMig_Out_Male) %>%
  left_join(x = ., 
            y = AgeSpecific_Mig_Net_Male) %>% 
  mutate("Year" = as.numeric(Year)) %>% 
  mutate("AgeID_New" = match(AgeGroup, unique(AgeGroup))) %>% 
  mutate("AgeGroup_New" = age_groups_new[AgeID_New]) %>% 
  mutate("RegionName" = RegionNameHelper$RegionName[
    match(RegionNumber, RegionNameHelper$RegionNumber)]) %>% 
  filter(grepl("094", RegionNumber))

#combine both in and out migration and add net migration column
AgeSpecific_Mig_Female <- 
  left_join(x = AgeSpecificMig_In_Female,
            y = AgeSpecificMig_Out_Female) %>% 
  left_join(x = ., 
            y = AgeSpecific_Mig_Net_Female) %>% 
  mutate("Year" = as.numeric(Year)) %>% 
  mutate("AgeID_New" = match(AgeGroup, unique(AgeGroup))) %>% 
  mutate("AgeGroup_New" = age_groups_new[AgeID_New]) %>% 
  mutate("RegionName" = RegionNameHelper$RegionName[
    match(RegionNumber, RegionNameHelper$RegionNumber)]) %>% 
  filter(grepl("094", RegionNumber))

### Imputation of Migration values using CART #################################
## Missing only in 2019 and 2020, potentially error in the data ... 
AgeSpecific_Mig_Female %>% 
  reframe("min" = min(In, na.rm = TRUE), .by = Year) 

AgeSpecific_Mig_Female %>% 
  reframe("Missings" = sum(is.na(Out)), .by = Year)

AgeSpecific_Mig_Male %>% 
  reframe("Missings" = sum(is.na(Out)), .by = Year)

#thus imputation with cart fine. 
MigFemaleData_Imp <- 
  AgeSpecific_Mig_Female %>% select(Year, RegionNumber, In, Out, Net, AgeID_New) %>% 
  mutate("RegID" = factor(match(RegionNumber, unique(RegionNumber)))) %>% 
  select(-RegionNumber) %>% 
  select(-In)


meth <- mice::make.method(MigFemaleData_Imp)
meth["Out"] <- "cart"

pred <- mice::make.predictorMatrix(MigFemaleData_Imp)

Im_Fem_Mig <- mice::mice(MigFemaleData_Imp, 
                         meth = meth, pred = pred, m = 1)

Im_Fem_Mig$imp

AgeSpecific_Mig_Female[is.na(AgeSpecific_Mig_Female$Out),"Out"] <- Im_Fem_Mig$imp$Out
AgeSpecific_Mig_Female[,"In"] <- AgeSpecific_Mig_Female$Net + AgeSpecific_Mig_Female$Out


#### Impuataion Model for Males ################################################
MigMaleData_Imp <- 
  AgeSpecific_Mig_Male %>% select(Year, RegionNumber, In, Out, Net, AgeID_New) %>% 
  mutate("RegID" = factor(match(RegionNumber, unique(RegionNumber)))) %>% 
  select(-RegionNumber) %>% 
  select(-In)


meth <- mice::make.method(MigMaleData_Imp)
meth["Out"] <- "cart"

pred <- mice::make.predictorMatrix(MigMaleData_Imp)

Im_Male_Mig <- mice::mice(MigMaleData_Imp, 
                         meth = meth, pred = pred, m = 1)

Im_Male_Mig$imp

AgeSpecific_Mig_Male[is.na(AgeSpecific_Mig_Male$Out),"Out"] <- Im_Male_Mig$imp$Out
AgeSpecific_Mig_Male[,"In"] <- AgeSpecific_Mig_Male$Net + AgeSpecific_Mig_Male$Out



### add other migration data (totals only from 2000 - 2010) ##
NetMigCounts_BY_Tot_2000_10 <- 
  openxlsx::read.xlsx(
    xlsxFile =file.path(getwd(), "Data/BY_NetMigration_Total_2000-2011.xlsx"),
    sheet = "Total", colNames = TRUE
  ) %>% 
  filter(grepl("094", RegionNumber)) %>% 
  filter(Year != 2011) #already in data


OutMigCounts_BY_Tot_2000_10 <- 
  openxlsx::read.xlsx(
    xlsxFile =file.path(getwd(), "Data/BY_OutMigration_Total_2000-2011.xlsx"),
    sheet = "Total", colNames = TRUE
  ) %>% 
  filter(grepl("094", RegionNumber))

InMigCounts_BY_Tot_2000_10 <-
  openxlsx::read.xlsx(
    xlsxFile =file.path(getwd(), "Data/BY_InMigration_Total_2000-2011.xlsx"),
    sheet = "Total", colNames = TRUE
  ) %>% 
  filter(grepl("094", RegionNumber))
  
MigCounts_Female_Tot_2000_10 <- 
  left_join(x = select(NetMigCounts_BY_Tot_2000_10, -NetMale),
            y = OutMigCounts_BY_Tot_2000_10) %>% 
  left_join(x = ., 
            y = InMigCounts_BY_Tot_2000_10) %>% 
  select(-c(OutMale,InMale)) %>% 
  rename("NetMig" = "NetFemale",
         "OutMig" = "OutFemale",
         "InMig" = "InFemale") 

MigCounts_Male_Tot_2000_10 <- 
  left_join(x = select(NetMigCounts_BY_Tot_2000_10, -NetFemale),
            y = OutMigCounts_BY_Tot_2000_10) %>% 
  left_join(x = ., 
            y = InMigCounts_BY_Tot_2000_10) %>% 
  select(-c(OutFemale,InFemale)) %>% 
  rename("NetMig" = "NetMale",
         "OutMig" = "OutMale",
         "InMig" = "InMale") 


###### CALCULATE TOTAL MIGRATION Net #####################################
#Calculation of total net migration for each region from 2000 - 2023
MigCounts_Female_Tot <- 
  AgeSpecific_Mig_Female %>% 
  reframe(InMig = sum(In),
          OutMig = sum(Out),
          NetMig = sum(Net),
          .by = c(RegionNumber, RegionName ,Year))%>% 
  select(Year, RegionNumber, NetMig, OutMig, InMig, RegionName) %>% 
  full_join(x = MigCounts_Female_Tot_2000_10, 
            y = ., 
            by =c("Year","RegionNumber","RegionName","OutMig","InMig","NetMig")) %>% 
  mutate(RegionName = stringr::str_trim(RegionName))  #Remove WhiteSpaces


MigCounts_Male_Tot <- 
  AgeSpecific_Mig_Male %>% 
  reframe(InMig = sum(In),
          OutMig = sum(Out),
          NetMig = sum(Net),
          .by = c(RegionNumber, RegionName ,Year))%>% 
  select(Year, RegionNumber, NetMig, OutMig, InMig, RegionName) %>% 
  full_join(x = MigCounts_Male_Tot_2000_10, 
            y = ., 
            by =c("Year","RegionNumber","RegionName","OutMig","InMig","NetMig")) %>% 
  mutate(RegionName = stringr::str_trim(RegionName))  #Remove WhiteSpaces


##### Loading geo Data #########################################################
GermanyData <- sf::read_sf(
  file.path(getwd(),"Data/Maps/vg250_01-01.tm32.shape.ebenen/vg250_ebenen_0101/VG250_KRS.shp"),
  as_tibble = FALSE)

#Filter Bavaria from German Data
OF <- which(grepl("094",GermanyData$AGS)) #Only Upper Franconia Regions
UpperFranconia <- GermanyData[OF,] #Dataset with only Bavaria


save(UpperFranconia,
     SRB,
     PopCounts, 
     Population_Old_Fem,
     DeathCounts, 
     BirthCounts,
     MigCounts_Female_Tot,
     AgeSpecific_Mig_Female,
     MigCounts_Male_Tot,
     AgeSpecific_Mig_Male,
     file = file.path(getwd(),"Data/TotalData.RData"))
