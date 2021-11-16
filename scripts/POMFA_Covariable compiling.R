
# clean 2015-2018 QU39 data: POMS SI, Chl+phaeopigments, Nutrients, CTD temp and salinity, rain + wind data

# Compile all data for the time period in one data set
# Compile all covariables for FA dates in another


library(tidyverse)
library(lubridate)
library(oce)


# Load FA data to extract dates later 
POMFAdata <- read.csv("POMFA data package 20210728.csv", stringsAsFactors = FALSE, na=c("", "NA"))
POMFAdata$Date <- as.Date(POMFAdata$Date, "%Y-%m-%d")
# Make data frame of Dates with FA data
FAdates <- (POMFAdata$Date)



# Isotopes and elements ---------------------------------------------------

# QU39 POMs isotope data compiling and filtering 

# Load POM Isotope Data 
POMisotopeDataQU39 <- read.csv("~amclaskey/Google Drive (akmclaskey@gmail.com)/UBC Work/Work/Fatty Acid Data/Data analysis POM FA/POMFA Data Files to Share/Covariates raw data/2021-07-27_190508_HakaiData_poms.csv", stringsAsFactors = FALSE, na=c("", "NA"))


POMisotopeDataQU39 <- read.csv("Covariates raw data/2021-07-27_190508_HakaiData_poms.csv", stringsAsFactors = FALSE, na="")
POMisotopeDataQU39$Date <- as.Date(POMisotopeDataQU39$Date, "%Y-%m-%d")

# use carbon isotopes from acidified samples, nitrogen and C:N from not acidified samples
# Need acidified and non-acidified ug.C, ug.N, corr_delta15n, corr_delta13c
# Also need Volume filtered, C.Flag, N.Flag, date, Line.Out.Depth
smallPOMisoData <- POMisotopeDataQU39 %>% select(Date, Line.Out.Depth, Acidified, Volume..ml.,  corr_delta15n, corr_delta13c, ug.C, ug.N, C.Flag, N.Flag) 

# Do all samples have acidified or not?
table(is.na(smallPOMisoData$Acidified))
# Four are missing Acidified. Remove them 
smallPOMisoData <- smallPOMisoData %>%filter(!is.na(Acidified))

# Check for sample dates+depths where both samples were not Acidified, or both were, so they don't have unique combo of keys.
smallPOMisoData[duplicated(smallPOMisoData[,1:3]),]
smallPOMisoData[duplicated(smallPOMisoData[,1:3], fromLast = T),]
# There are four dates+depths with duplicate samples

# There is only one sample from 0 or 5 m, looks like both samples really weren't Acidified. Just going to remove duplicates for now
smallPOMisoData.new <- smallPOMisoData %>% distinct(Date, Line.Out.Depth, Acidified, .keep_all = T)


# Check quality flags
table(smallPOMisoData.new$C.Flag)
table(is.na(smallPOMisoData.new$C.Flag))
table(smallPOMisoData.new$N.Flag)
table(is.na(smallPOMisoData.new$N.Flag))

# Acidified samples are either good for carbon or aren't good (i.e. flags on nitrogen don't matter)
# Non-Acidified samples  could be good for just nitrogen, or also C_N (i.e. flags on C and N both matter)

# Split the dataset in two to QC acidfied and not-Acidified separately
smallPOMisoData.acid <- smallPOMisoData.new %>% filter(Acidified=="TRUE")
smallPOMisoData.notacid <- smallPOMisoData.new %>% filter(Acidified=="FALSE")

# Only the carbon matters out of the Acidified sample
smallPOMisoData.acid.noCflag <- smallPOMisoData.acid %>% filter(is.na(C.Flag))
# both N and C matter out of the unAcidified sample, although I may want to separate these later to maximize the amount of data I have 
smallPOMisoData.notacid.noflags <- smallPOMisoData.notacid %>% filter(is.na(N.Flag) & is.na(C.Flag))
# I need them both for unAcidified samples because I'm calculated C:N, but only need the C.Flag for Acidified samples
smallPOMisoData.QC <- rbind(smallPOMisoData.acid.noCflag, smallPOMisoData.notacid.noflags)

# Before spreading the samples I need to remove N.Flag because it is different between the rows now
smallPOMisoData.QC <- smallPOMisoData.QC %>% select(-N.Flag)

# Calculate C:N 
smallPOMisoData.QC <- smallPOMisoData.QC %>% mutate(C_N =  ((ug.C*12.01) / (ug.N*14.001)) )

# Plot to check vol filtered
plot(smallPOMisoData.QC$Volume..ml. ~ smallPOMisoData.QC$Date)

# Calculate elemental concentrations based on volume filtered
smallPOMisoData.QC <- smallPOMisoData.QC %>% 
  mutate(C_ug.L =  (ug.C / (Volume..ml./1000))) %>% 
  mutate(N_ug.L =  (ug.N / (Volume..ml./1000))) 

# spread the acidified/non-acidified samples 
# Remove Vol filtered, ug.C, and ug.N which aren't meaningful to spread
smallPOMisoData.QC <- smallPOMisoData.QC %>% select(-c(Volume..ml., ug.C, ug.N, C.Flag))
smallPOMisoData_long <- pivot_wider(smallPOMisoData.QC, names_from = Acidified, values_from = c(C_N, corr_delta15n, corr_delta13c, C_ug.L, N_ug.L))

# I don't find true=acidified and false=not acidified to be very helpful. Change col names 
# CHECK TO MAKE SURE COLS ARE STILL IN SAME ORDER
colnames(smallPOMisoData_long)[3:12]<- c("C_N_acidified", "C_N_notacid",  "delta15n.acidified","delta15n.notacid",
                                         "delta13c.acidified", "delta13c.notacid", "C_ug.L.acidified", "C_ug.L.notacid", 
                                         "N_ug.L.acidified", "N_ug.L.notacid")

# Look at distribution of sampling depths
table(smallPOMisoData_long$Line.Out.Depth)


# compare acidified vs not acidified 
plot(smallPOMisoData_long$C_ug.L.acidified ~ smallPOMisoData_long$C_ug.L.notacid)
abline(0,1)

plot(smallPOMisoData_long$N_ug.L.acidified ~ smallPOMisoData_long$N_ug.L.notacid)
abline(0,1)

plot(smallPOMisoData_long$delta15n.acidified ~ smallPOMisoData_long$delta15n.notacid)
abline(0,1)
plot(smallPOMisoData_long$delta13c.acidified ~ smallPOMisoData_long$delta13c.notacid)
abline(0,1)


plot((smallPOMisoData_long$C_ug.L.notacid-smallPOMisoData_long$C_ug.L.acidified) ~ smallPOMisoData_long$C_ug.L.acidified)
abline(0, 0)
abline(v=100)

plot((smallPOMisoData_long$delta13c.notacid-smallPOMisoData_long$delta13c.acidified) ~ smallPOMisoData_long$C_ug.L.acidified)
abline(0, 0)
abline(v=100)




## filter out 0 and 5 m only
# Average the zero and 5 m samples 
POMisoData_0and5m <- smallPOMisoData_long %>% filter(Line.Out.Depth==5 | Line.Out.Depth==0)
POMisoData_0and5m.agg = aggregate(POMisoData_0and5m,
                                  by = list(POMisoData_0and5m$Date),
                                  FUN = mean, na.rm=TRUE)


# Isotope data all dates
POMisoData_all <- POMisoData_0and5m.agg %>% select(-Group.1, -Line.Out.Depth) %>% 
  filter(Date >= "2015-01-01")

# Isotope data for FA dates only
POMisoData_FAdates <- POMisoData_all %>% filter(Date %in% FAdates)
length(unique(POMisoData_FAdates$Date))
# Should be 133





# Chlorophyll -------------------------------------------------------------

chlData <- read.csv("~amclaskey/Google Drive (akmclaskey@gmail.com)/UBC Work/Work/Fatty Acid Data/Data analysis POM FA/POMFA Data Files to Share/Covariates raw data/2021-01-05_104110_HakaiData_chlorophyll.csv", stringsAsFactors = FALSE, na=c("", "NA"))

chlData <- read.csv("Covariates raw data/2021-01-05_104110_HakaiData_chlorophyll.csv")
chlData$date <- as.Date(chlData$date, "%m/%d/%y")
str(chlData)
# Remove 2019 values
chlData <- chlData %>% filter(date<as.Date("2019-01-01")) 

chlSumm.FA2 <- chlData %>% select(date, line_out_depth, volume, 
                                  hakai_id, filter_type, before_acid, after_acid, acid_flag,
                                  chla, chla_flag, chla_final, phaeo, phaeo_flag, phaeo_final,
                                  row_flag, quality_level, quality_log)


# Use Chla_final to avoid negative values. Those below detection limit have chla_flag= "DL"

# Looking at the acid ratio is an important QC step
plot(chlSumm.FA2$before_acid ~ chlSumm.FA2$after_acid)
abline(0,1)
# There are a couple outliers

unique(chlSumm.FA2$chla_flag)
# Let's start by just removing them. 

# There are three different type of quality flags, acid, chl, phaeopigment
# They don't always correspond
# Do I need to pay attention to quality_level ? Most are PI, some technicianmr or technichianm

# remove SVC and SVD, leave DL (detection limit)
unique(chlSumm.FA2$phaeo_flag)
chl.noflags <- chlSumm.FA2[chlSumm.FA2$chla_flag %in% c("AV", "DL"),]
phaeo.noflags <- chlSumm.FA2[chlSumm.FA2$phaeo_flag %in% c("AV", "DL"),]



# Although the units are not actually indicated, they are already ug/L
 chl.noflags <- chl.noflags %>% mutate(chla_ug.L = chla_final)
 phaeo.noflags <- phaeo.noflags %>% mutate(phaeo_ug.L = phaeo_final)

# there are some sampleS with negative phaeo measurements but marked AV. Going to remove for now. 
negative.phaeo <- phaeo.noflags[phaeo.noflags$phaeo_final<=0,]
phaeo.noflags <- phaeo.noflags[phaeo.noflags$phaeo_final>=0,]


# Will need to spread filter type. Make reduced dataset
colnames(chl.noflags)
chl.noflags <- select(chl.noflags, date, line_out_depth, filter_type, chla_ug.L)
# Check for duplicates
chl.noflags[duplicated(chl.noflags[1:3]),]
chl.noflags[duplicated(chl.noflags[1:3], fromLast = T),]
# Triplicate samples were taken on 2016-04-01
# There are two 10 m Bulk GF/F samples on 2016-04-14, but one is much lower than the other. 
chl.noflags[chl.noflags$date=="2016-04-14" & chl.noflags$line_out_depth==10,]
# I think this is a mistake and the low one is actually the GF/F size fraction

# UPDATE THIS -------------------------------------------------------------
# UPDATE THis MANUAL FIX IF NEW FILE DOWNLOADED
chl.noflags[1047,3] <- "GF/F"

# Remove triplicates
chl.noflags.2 <- chl.noflags %>% distinct(date, line_out_depth, filter_type, .keep_all = T)

chl.noflags.2[duplicated(chl.noflags.2[1:3]),]

# Spread
chl.noflags.wide <- chl.noflags.2 %>% pivot_wider(names_from = filter_type, values_from = chla_ug.L) 
colnames(chl.noflags.wide)[3:7] <- paste("chl", colnames(chl.noflags.wide)[3:7], sep = "_")


phaeo.noflags <- select(phaeo.noflags, date, line_out_depth, filter_type, phaeo_ug.L)

phaeo.noflags[duplicated(phaeo.noflags[1:3]),]
# There are two 10 m Bulk GF/F samples on 2016-04-14, but one is much lower than the other. 
phaeo.noflags[phaeo.noflags$date=="2016-04-14" & phaeo.noflags$line_out_depth==10,]
phaeo.noflags[1026,] 
phaeo.noflags[1026,3] <- "GF/F"

# Remove the triplicates
phaeo.noflags.2 <- phaeo.noflags %>% distinct(date, line_out_depth, filter_type, .keep_all = T)

phaeo.noflags.wide <- phaeo.noflags.2 %>% pivot_wider(names_from = filter_type, values_from = phaeo_ug.L) 
colnames(phaeo.noflags.wide)[3:7] <- paste("phaeo", colnames(phaeo.noflags.wide)[3:7], sep = "_")

# I want the four size classes of Chl and Phaeo (8 values total) for each date and depth. 
chl.phaeo.noflags <- full_join(chl.noflags.wide, phaeo.noflags.wide)

colnames(chl.phaeo.noflags)
colnames(chl.phaeo.noflags)[c(5,7,10,12)] <- c("chl_GF.F", "chl_Bulk.GF.F", "phaeo_GF.F", "phaeo_Bulk.GF.F")



# On the first two sampling dates 2 um filters were used rather than 3 um. Use values (these aren't FA dates i.e. 
# these are not in any models, just the plots. Then remove those cols
chl.phaeo.noflags[1:8,6] <- chl.phaeo.noflags[1:8,3]
chl.phaeo.noflags[1:8,11] <- chl.phaeo.noflags[1:8,8]
chl.phaeo.noflags <- chl.phaeo.noflags %>% select(-chl_2um, -phaeo_2um)



# Find mean of 0 and 5 m 
chl.phaeo.noflags.0_5m <- chl.phaeo.noflags %>% 
  filter(line_out_depth<=5)

# Find mean by filter size
chl.phaeo.noflags.agg.mean = aggregate(chl.phaeo.noflags.0_5m,
                                  by = list(chl.phaeo.noflags.0_5m$date),
                                  FUN = mean, na.rm=TRUE)
# remove extra cols
chl.phaeo_all <- chl.phaeo.noflags.agg.mean[-c(1)] 
# Now I have my file of all dates 
colnames(chl.phaeo_all)[1] <- "Date"

# Calculate ratio of phaeopigments/chl a 
chl.phaeo_all <- chl.phaeo_all %>%   mutate(phaeo.chla = phaeo_Bulk.GF.F / chl_Bulk.GF.F) 


# filter dataset to only include FA Dates
chl.phaeo_FAdates <- chl.phaeo_all %>% filter(Date %in% FAdates)
length(unique(chl.phaeo_FAdates$Date))
# Five FA dates are missing chl data
# Should be 128



# Calculate SumChl for days that have all three size fracs

chl.phaeo_FAdates <- chl.phaeo_FAdates %>%  mutate(SumChl = (chl_20um + chl_3um + chl_GF.F))
chl.phaeo_FAdates$SumChl




# Nutrients ---------------------------------------------------------------

# 4/24/20 Chris Mackenzie brought the nutrient replicates to my attention again.
# But it looks like when I filter out 2015-2018 0-30 m there are only a few replicates taken in Mar 2015

# All nutrients are in umol/L despite what the column header says
# nutQU39 <- read.csv("~amclaskey/Google Drive (akmclaskey@gmail.com)/UBC Work/Work/Fatty Acid Data/Data analysis POM FA/POMFA Data Files to Share/Covariates raw data/
#                     2021-01-08_091656_HakaiData_nutrients.csv", stringsAsFactors = FALSE, na=c("", "NA"))
# nutQU39 <- read.csv("Covariates raw data/2021-01-08_091656_HakaiData_nutrients.csv", na="")



# Updated data file 2021-06-10
nutQU39 <- read.csv("~amclaskey/Google Drive (akmclaskey@gmail.com)/UBC Work/Work/Fatty Acid Data/Data analysis POM FA/POMFA Data Files to Share/Covariates raw data/2021-06-15_133145_HakaiData_nutrients.csv", stringsAsFactors = FALSE, na=c("", "NA"))
nutQU39 <- read.csv("Covariates raw data/2021-06-15_133145_HakaiData_nutrients.csv", na="")

nutQU39$date <- as.Date(nutQU39$date, "%m/%d/%y")

# Subset data to remove 2019 samples and deep samples
nutQU39 <- nutQU39 %>% 
  filter(date<"2019-01-01" & line_out_depth<=30)






# Plot the data to visually check
plot(nutQU39$no2_no3_um[nutQU39$line_out_depth==0] ~ nutQU39$date[nutQU39$line_out_depth==0])
plot(nutQU39$no2_no3_um[nutQU39$line_out_depth==5] ~ nutQU39$date[nutQU39$line_out_depth==5])
plot(nutQU39$no2_no3_um[nutQU39$line_out_depth==10] ~ nutQU39$date[nutQU39$line_out_depth==10])
plot(nutQU39$no2_no3_um[nutQU39$line_out_depth==20] ~ nutQU39$date[nutQU39$line_out_depth==20])
plot(nutQU39$no2_no3_um[nutQU39$line_out_depth==30] ~ nutQU39$date[nutQU39$line_out_depth==30])
# Those all look fine
plot(nutQU39$po4[nutQU39$line_out_depth==0] ~ nutQU39$date[nutQU39$line_out_depth==0])
plot(nutQU39$po4[nutQU39$line_out_depth==5] ~ nutQU39$date[nutQU39$line_out_depth==5])
plot(nutQU39$po4[nutQU39$line_out_depth==10] ~ nutQU39$date[nutQU39$line_out_depth==10])
plot(nutQU39$po4[nutQU39$line_out_depth==20] ~ nutQU39$date[nutQU39$line_out_depth==20])
plot(nutQU39$po4[nutQU39$line_out_depth==30] ~ nutQU39$date[nutQU39$line_out_depth==30])
# Those all look fine
plot(nutQU39$sio2[nutQU39$line_out_depth==0] ~ nutQU39$date[nutQU39$line_out_depth==0])
plot(nutQU39$sio2[nutQU39$line_out_depth==5] ~ nutQU39$date[nutQU39$line_out_depth==5])
plot(nutQU39$sio2[nutQU39$line_out_depth==10] ~ nutQU39$date[nutQU39$line_out_depth==10])
plot(nutQU39$sio2[nutQU39$line_out_depth==20] ~ nutQU39$date[nutQU39$line_out_depth==20])
plot(nutQU39$sio2[nutQU39$line_out_depth==30] ~ nutQU39$date[nutQU39$line_out_depth==30])


# Detection limit for nitrate analysis is 0.036 umol L-1
nutQU39.lowNO3 <- nutQU39 %>% 
  filter(no2_no3_um<0.036)
# Detection limit for phosphate analysis is 0.008 umol L-1
nutQU39.lowPO4 <- nutQU39 %>% 
  filter(po4<0.008)
# Detection limit for silicate analysis is 0.1 umol L-1
nutQU39.lowSiO2 <- nutQU39 %>% 
  filter(sio2<0.1)


# Check for data flags
table(nutQU39$no2_no3_flag)
table(nutQU39$po4_flag)
table(nutQU39$sio2_flag)

# Remove row with flags
nutQU39.2 <- nutQU39 %>% filter(no2_no3_flag!="SVC" & no2_no3_flag!="SVD" | is.na(no2_no3_flag))


table(nutQU39$rn)
# There are only three replicates in this subset of the data. 

# rn stands for replicate number. This removes repeats
nutQU39 <- nutQU39[nutQU39$rn==1,]

# Check for other repeats
nutQU39[duplicated(nutQU39[,c(5,19)]),]
nutQU39[duplicated(nutQU39[,c(5,19)], fromLast = T),]
# 236, 237, 238, 816, 817

# The first one is different from the other two, drop first two
# Other ones are actual repeated data (both QNUT4904)
nutQU39[c(233:235,813:814),] # this is a manual deletion so make sure this is actually what you want to remove
nutQU39 <- nutQU39[-c(233:234,813),]

# Make smaller dataset to aggregate by date
nutrientsOnly <- nutQU39 %>% select(date, line_out_depth, no2_no3_um, po4, sio2) 



# Find mean of 0 and 5 m
nutrientsOnly.0_5 <- nutrientsOnly %>%  filter( line_out_depth<=5)
nutrientsOnly.agg = aggregate(nutrientsOnly.0_5,
                                   by = list(nutrientsOnly.0_5$date),
                                   FUN = mean, na.rm=TRUE)
# # When there are no measurements, returns -Inf
 nutrientsOnly.agg <- nutrientsOnly.agg %>%  mutate_if(is.numeric, ~ replace(., is.infinite(.), NA))

# # remove extra cols
 nutrients_all <- nutrientsOnly.agg[-c(1,3)]
 # Now I have my file of all dates
 colnames(nutrients_all) <- c("Date", "NO3", "PO4", "Si")

# # filter dataset to only include FA Dates
 nutrients_FAdates <- nutrients_all %>% filter(Date %in% FAdates)
 length(unique(nutrients_FAdates$Date))
# # Should be 133



# # 5 m only
# nutrientsOnly.5m <- nutrientsOnly %>% filter(line_out_depth==5)
# 
# # Remove depth column and rename Cols
# colnames(nutrientsOnly.5m)
# nutrientsOnly.5m <- nutrientsOnly.5m[-c(2)] 
# colnames(nutrientsOnly.5m) <- c("Date", "N", "P", "Si")
# 
# 
# # filter dataset to only include FA Dates
# nutrientsOnly.5m_FAdates <- nutrientsOnly.5m %>% filter(Date %in% FAdates)
# length(unique(nutrientsOnly.5m_FAdates$Date))
# # 133
# nutrientsOnly.5m_FAdates.comp <- nutrientsOnly.5m_FAdates[complete.cases(nutrientsOnly.5m_FAdates),]
# # 132 days w complete nutrients
# length(unique(nutrientsOnly.5m_FAdates.comp$date))
# 
# 
# nutrientsOnly.5m_FAdates
# 
# # Rename dataframe for compiling 
# nutrients_FAdates <- nutrientsOnly.5m_FAdates
# nutrients_all <- nutrientsOnly.5m



# CTD data ----------------------------------------------------------------

# This is CTD data downloaded for QU24 and QU39 from 2015-01-01 to 2018-12-31
# CTDdata <- read.csv("Covariates raw data/ctd-bulk-1610561316379.csv", stringsAsFactors = F)
# Up to date version 2021-06-14
 CTDdata <- read.csv("Covariates raw data/ctd-bulk-1623709310669.csv", stringsAsFactors = F)
 
str(CTDdata)
CTDdata$Start.timeDate <- as.Date(CTDdata$Start.time, "%Y-%m-%d  %H:%M:%S")
CTDdata$Start.time <- as.POSIXct(CTDdata$Start.time)

length(unique(CTDdata$Cast.PK))
length(unique(CTDdata$Start.timeDate))


casts <- unique(CTDdata$Cast.PK)
times <- unique(CTDdata$Start.time)

# Plot a profile
plot(CTDdata$Depth..m.[CTDdata$Cast.PK==casts[2]]~CTDdata$Fluorometry.Chlorophyll..ug.L.[CTDdata$Cast.PK==casts[2]], ylim=c(50,0))
# Plot several profiles
for(i in 201:278){
  plot(CTDdata$Depth..m.[CTDdata$Cast.PK==casts[i]]~CTDdata$Fluorometry.Chlorophyll..ug.L.[CTDdata$Cast.PK==casts[i]], ylim=c(50,0))
  abline(h=10)
  abline(h=5)
}


# *Calculate buoyancy flux --------------------------------------

CTDdata$N2 <- rep(NA, nrow(CTDdata))
for(i in 1:length(casts)){
  ctd <- as.ctd(salinity = CTDdata$Salinity..PSU.[CTDdata$Cast.PK==casts[i]], 
                temperature = CTDdata$Temperature..deg.C.[CTDdata$Cast.PK==casts[i]], 
                pressure = CTDdata$Pressure..dbar.[CTDdata$Cast.PK==casts[i]]) 
  CTDdata$N2[CTDdata$Cast.PK==casts[i]] <- swN2(ctd)
}

plot(CTDdata$N2)


# Calculate summary stats 

castSummary<- matrix(nrow=length(casts), ncol=12, byrow=TRUE)
castSummary<- as.data.frame(castSummary)
for(i in 1:length(casts)){
  castSummary[i,1] <- as.character(unique(CTDdata$Start.time[CTDdata$Cast.PK==casts[i]]))
  castSummary[i,2] <- (unique(CTDdata$Cast.PK[CTDdata$Cast.PK==casts[i]]))
# SST is mean temperature between 4-6 m inclusive
  castSummary[i,3] <- mean(CTDdata$Temperature..deg.C.[which(CTDdata$Cast.PK==casts[i] & CTDdata$Depth..m.<=6 & CTDdata$Depth..m.>=4)])
# SSS is mean salinity between 4-6 m inclusive
  castSummary[i,4] <- mean(CTDdata$Salinity..PSU.[which(CTDdata$Cast.PK==casts[i] & CTDdata$Depth..m.<=6 & CTDdata$Depth..m.>=4)])
  castSummary[i,5] <- max(CTDdata$Fluorometry.Chlorophyll..ug.L.[which(CTDdata$Cast.PK==casts[i])])
  castSummary[i,6] <- mean(CTDdata$Fluorometry.Chlorophyll..ug.L.[which(CTDdata$Cast.PK==casts[i] & CTDdata$Depth..m.<=30)])
  castSummary[i,7] <- min(CTDdata$Dissolved.O2..mL.L.[which(CTDdata$Cast.PK==casts[i] & CTDdata$Dissolved.O2..mL.L.>0)], na.rm = T)
  castSummary[i,8] <- (unique(CTDdata$Start.timeDate[CTDdata$Cast.PK==casts[i]]))
# Adding in min fluorescence as a QC check
  castSummary[i,9] <- min(CTDdata$Fluorometry.Chlorophyll..ug.L.[which(CTDdata$Cast.PK==casts[i])])
  castSummary[i,10] <- mean(CTDdata$Temperature..deg.C.[which(CTDdata$Cast.PK==casts[i] & CTDdata$Depth..m.<=30 & CTDdata$Depth..m.>=4)])
  castSummary[i,11] <- mean(CTDdata$N2[which(CTDdata$Cast.PK==casts[i] & CTDdata$Depth..m.<=30 & CTDdata$Depth..m.>=0.2)])
  castSummary[i,12] <- as.character(unique(CTDdata$Station[CTDdata$Cast.PK==casts[i]]))
  # There are some extra lines w/o data on a few profiles that give me NAs when I calculate buoyancy frequency with them included
}

table(castSummary$V12)

names(castSummary) <- c("date", "Cast.PK", "SST", "SSS", "Max.Chl","Mean.Chl", "Min.O2", "Start.timeDate", "Min.Chl", "30mTemp", "N2.30m", "Station")
str(castSummary)
castSummary$Start.timeDate <- as.Date.numeric(castSummary$Start.timeDate, origin = "1970-01-01")
castSummary$date <- as.POSIXct(castSummary$date)

# Calculate log buoyancy flux
castSummary$log10.N2 <- log10(castSummary$N2.30m)
plot(castSummary$log10.N2 ~ castSummary$date)


# Take out QU24 casts from 2015-03-18 onwards, clean columns
CTD_all <- castSummary %>% 
  filter(Station=="QU39" | Start.timeDate<"2015-03-18") %>% 
  select(-date, -Cast.PK, -Station)  %>% 
  dplyr::rename(Date = Start.timeDate)
length(unique(CTD_all$Date))

# Sometimes there were multiple casts per day 
CTD_all.dups <- CTD_all[duplicated(CTD_all$Date),]
CTD_all.dups2 <- CTD_all[duplicated(CTD_all$Date, fromLast = T),]
CTD_all.dups.comp <- full_join(CTD_all.dups, CTD_all.dups2)

dates.CTD <- unique(CTD_all.dups.comp$Date)
# Plot SSS and SST replicates
for(i in 1:67){
  plot(CTD_all.dups.comp$SSS[CTD_all.dups.comp$Date==dates.CTD[i]]~CTD_all.dups.comp$SST[CTD_all.dups.comp$Date==dates.CTD[i]], ylim=c(25,33), xlim=c(6,20))
}

# Remove duplicates
CTD_all <- CTD_all %>% distinct(Date, .keep_all = T)


# filter dataset to only include FA Dates
CTD_FAdates <- CTD_all %>% filter(Date %in% FAdates)
length(unique(CTD_FAdates$Date))
# Two FA sampling dates do not have CTD casts (and one more has NaN data)
# Should be 131




# Rain  -------------------------------------------------------------------

# Rain data are only calculated for the FA sampling dates only, not the whole continuous time series

QuadraRain<- read.csv("Covariates raw data/2021-01-11.1hourSamples_Rain.csv", stringsAsFactors = F)
# I think the last row of header makes the most sense for the col names
colnames(QuadraRain)  <- (QuadraRain[3,])
# Remove top three rows
QuadraRain<- QuadraRain[-c(1:3),]

# Need to reformat dates
QuadraRain$Date.time <- as.POSIXct(QuadraRain$measurementTime, tz = "", format="%Y-%m-%d %H:%M")
QuadraRain$Rain <- as.numeric(QuadraRain$Rain)

# For rain use total over a time period

# Calculate total Rain over periods prior to each FA date
Rain_FAdates <- as.data.frame(matrix(ncol=6, nrow=length(FAdates)))
for(i in 1:length(FAdates)){
  Rain_FAdates[i,1] <- FAdates[i]
  Rain_FAdates[i,2] <- sum(QuadraRain$Rain[QuadraRain$Date.time<=FAdates[i] & QuadraRain$Date.time>(FAdates[i] - days(1))], na.rm = T)
  Rain_FAdates[i,3] <- sum(QuadraRain$Rain[QuadraRain$Date.time<=FAdates[i] & QuadraRain$Date.time>(FAdates[i] - days(2))], na.rm = T)
  Rain_FAdates[i,4] <- sum(QuadraRain$Rain[QuadraRain$Date.time<=FAdates[i] & QuadraRain$Date.time>(FAdates[i] - days(3))], na.rm = T)
  Rain_FAdates[i,5] <- sum(QuadraRain$Rain[QuadraRain$Date.time<=FAdates[i] & QuadraRain$Date.time>(FAdates[i] - days(5))], na.rm = T)
  Rain_FAdates[i,6] <- sum(QuadraRain$Rain[QuadraRain$Date.time<=FAdates[i] & QuadraRain$Date.time>(FAdates[i] - days(10))], na.rm = T)
}
# I am taking the mean for the third col to deal with the NA that I'm not sure why it is happening
colnames(Rain_FAdates) <- c("Date", "Rain24hr", "Rain2Day", "Rain3Day","Rain5Day", "Rain10Day")
Rain_FAdates$Date <- as.Date(Rain_FAdates$Date, origin = "1970-01-01")



# Wind Speed --------------------------------------------------------------

WindSpAll <- read.csv("Covariates raw data/2021-01-11.1hourSamples_WindSpdAvg.csv", stringsAsFactors = F)
# I think the last row of header makes the most sense for the col names
colnames(WindSpAll)  <- (WindSpAll[3,])
# Remove top three rows
WindSpAll<- WindSpAll[-c(1:3),]

# Need to reformat dates
WindSpAll$Date.time <- as.POSIXct(WindSpAll$measurementTime, tz = "", format="%Y-%m-%d %H:%M")
WindSpAll$WindSpd_Avg <- as.numeric(WindSpAll$WindSpd_Avg)


# Calculate Mean wind speed over time intervals prior for each FA date
WindSp_FAdates <- as.data.frame(matrix(ncol=6, nrow=length(FAdates)))
for(i in 1:length(FAdates)){
  WindSp_FAdates[i,1] <- FAdates[i]
  WindSp_FAdates[i,2] <- mean(WindSpAll$WindSpd_Avg[WindSpAll$Date.time<=FAdates[i] & WindSpAll$Date.time>(FAdates[i] - days(1))], na.rm = T)
  WindSp_FAdates[i,3] <- mean(WindSpAll$WindSpd_Avg[WindSpAll$Date.time<=FAdates[i] & WindSpAll$Date.time>(FAdates[i] - days(2))], na.rm = T)
  WindSp_FAdates[i,4] <- mean(WindSpAll$WindSpd_Avg[WindSpAll$Date.time<=FAdates[i] & WindSpAll$Date.time>(FAdates[i] - days(3))], na.rm = T)
  WindSp_FAdates[i,5] <- mean(WindSpAll$WindSpd_Avg[WindSpAll$Date.time<=FAdates[i] & WindSpAll$Date.time>(FAdates[i] - days(5))], na.rm = T)
  WindSp_FAdates[i,6] <- mean(WindSpAll$WindSpd_Avg[WindSpAll$Date.time<=FAdates[i] & WindSpAll$Date.time>(FAdates[i] - days(10))], na.rm = T)
}
colnames(WindSp_FAdates) <- c("Date", "WindSp24hr", "WindSp2Day", "WindSp3Day", "WindSp5Day", "WindSp10Day")
WindSp_FAdates$Date <- as.Date(WindSp_FAdates$Date, origin = "1970-01-01")



# Wind Direction ---------------------------------------------------------------------

# These are mean wind direction in each hour period  
windDirecthourly<- read.csv("Covariates raw data/2021-01-11.1hourSamples_WindDirect.csv", stringsAsFactors = F)
# I think the last row of header makes the most sense for the col names
colnames(windDirecthourly)  <- (windDirecthourly[3,])
# Remove top three rows
windDirecthourly<- windDirecthourly[-c(1:3),]

# Need to reformat dates
windDirecthourly$Date.time <- as.POSIXct(windDirecthourly$measurementTime, tz = "", format="%Y-%m-%d %H:%M")
windDirecthourly$WindDir_Avg <- as.numeric(windDirecthourly$WindDir_Avg)


# Calculate Mean wind direction over time periods prior to each FA date
WindDir_FAdates <- as.data.frame(matrix(ncol=6, nrow=length(FAdates)))
for(i in 1:length(FAdates)){
  WindDir_FAdates[i,1] <- FAdates[i]
  WindDir_FAdates[i,2] <- mean(windDirecthourly$WindDir_Avg[windDirecthourly$Date.time<=FAdates[i] & windDirecthourly$Date.time>(FAdates[i] - days(1))], na.rm = T)
  WindDir_FAdates[i,3] <- mean(windDirecthourly$WindDir_Avg[windDirecthourly$Date.time<=FAdates[i] & windDirecthourly$Date.time>(FAdates[i] - days(2))], na.rm = T)
  WindDir_FAdates[i,4] <- mean(windDirecthourly$WindDir_Avg[windDirecthourly$Date.time<=FAdates[i] & windDirecthourly$Date.time>(FAdates[i] - days(3))], na.rm = T)
  WindDir_FAdates[i,5] <- mean(windDirecthourly$WindDir_Avg[windDirecthourly$Date.time<=FAdates[i] & windDirecthourly$Date.time>(FAdates[i] - days(5))], na.rm = T)
  WindDir_FAdates[i,6] <- mean(windDirecthourly$WindDir_Avg[windDirecthourly$Date.time<=FAdates[i] & windDirecthourly$Date.time>(FAdates[i] - days(10))], na.rm = T)
}
colnames(WindDir_FAdates) <- c("Date", "WindDir24hr", "WindDir2Day", "WindDir3Day", "WindDir5Day", "WindDir10Day")
WindDir_FAdates$Date <- as.Date(WindDir_FAdates$Date, origin = "1970-01-01")




# Microplankton microscopy ------------------------------------------------


# When I open this file in excel it messes up the dates so they are in two different formats
PP.QU39 <- read.csv("Covariates raw data/QU39_v1.2.csv", stringsAsFactors = FALSE, na="")

# Fix header and colnames
colnames(PP.QU39) <- PP.QU39[1,]
PP.QU39 <- PP.QU39[-c(1:3),]
colnames(PP.QU39)[1:2] <- c("group", "Taxa")

PP.QU39[3:123] <- as.numeric(unlist(PP.QU39[3:123]))
str(PP.QU39)


# Still need to change "2016-12-14" to "2016-12-15" and "2017-04-06" to "2017-04-07"
colnames(PP.QU39)[c(24,32)] <- c("2016-12-15", "2017-04-07")

# Calculate Shannon-Wiener diversity
library(vegan)
SWindex<- as.data.frame(matrix(nrow=121, ncol=2, byrow=TRUE))

for(i in 3:123){
  SWindex[i-2,1] <- colnames(PP.QU39)[i]
  SWindex[i-2,2] <- diversity(PP.QU39[,i][complete.cases(PP.QU39[,i])], index = "shannon", MARGIN = 1, base = exp(1))
}
colnames(SWindex) <- c("Date", "SW.Diversity")

#write.csv(SWindex, "Phyoplankton diversity SW index.csv", row.names = F)



# Transform to long form
PP.QU39.long <- PP.QU39 %>% pivot_longer("2016-02-23":"2019-07-24",names_to = "sampling.date", values_to = "Ind.per.L")
str(PP.QU39.long)


# Leave dates as characters to fix
# There are a few dates that are wrong in the phytoplankton community database (PP.QU39) that need to be changed

# 2017-04-06 in phyto community data is actually 2017-04-07 BUT NOT A FA DATE
# 2016-12-14 in phyto community data is actually 2016-12-15
# 2018-11-14 in phyto community data is actually 2018-11-19
# 2018-07-23 is not showing up but it was analyzed for community data I think it is miss-labeled as 2018-05-23
PP.QU39.long <- mutate(PP.QU39.long, sampling.date = ifelse(sampling.date == "14-12-2016", "15-12-2016", sampling.date))
PP.QU39.long <- mutate(PP.QU39.long, sampling.date = ifelse(sampling.date == "06-04-2017", "07-04-2017", sampling.date))
PP.QU39.long <- mutate(PP.QU39.long, sampling.date = ifelse(sampling.date == "14-11-2018", "19-11-2018", sampling.date))
PP.QU39.long <- mutate(PP.QU39.long, sampling.date = ifelse(sampling.date == "23-05-2018", "23-07-2018", sampling.date))

unique(PP.QU39.long$sampling.date)

PP.QU39.long$sampling.date <- as.Date(PP.QU39.long$sampling.date, "%Y-%m-%d")
PP.QU39.long$Ind.per.L <- as.numeric(PP.QU39.long$Ind.per.L)
str(PP.QU39.long)

# Filter out 2019 samples
PP.QU39.long <- PP.QU39.long %>% filter(sampling.date < "2019-01-01")

PP.QU39.long$group <- as.factor(PP.QU39.long$group)

# Make summary of Taxa for each date
taxaSummary <- PP.QU39.long %>%
      dplyr::group_by(sampling.date, group) %>%
  dplyr::summarise(sum = sum(Ind.per.L, na.rm = T)) %>% 
  pivot_wider(names_from = group, values_from = sum)
  

colnames(taxaSummary)[1] <- "Date"

# taxaSummary is the table of raw abundnaces of all taxa groups

# filter dataset to only include FA Dates
taxa_FAdates <- taxaSummary %>% filter(Date %in% FAdates)





# Bacterial cell counts ----------------------------------------------------

cellcounts <- read.csv("Covariates raw data/hakai_prokcell_counts.csv", stringsAsFactors = FALSE, na="")
cellcounts$Date <- as.Date(cellcounts$Date)
cellcounts$Counts <- as.numeric(cellcounts$Counts)

ggplot(cellcounts) + geom_point(aes(y=Counts, x=Date, color=as.factor(Depth)))

plot(cellcounts$Counts ~ cellcounts$Date)
plot(cellcounts$Counts[cellcounts$Depth==0] ~ cellcounts$Date[cellcounts$Depth==0])
plot(cellcounts$Counts[cellcounts$Depth==5] ~ cellcounts$Date[cellcounts$Depth==5])
plot(cellcounts$Counts[cellcounts$Depth==30] ~ cellcounts$Date[cellcounts$Depth==30])




# Assemble datasets -------------------------------------------------------


# covariates_FAdates
covariates_FAdates1 <- full_join(POMisoData_FAdates, chl.phaeo_FAdates)
covariates_FAdates2 <- full_join(covariates_FAdates1, nutrients_FAdates) # 133
covariates_FAdates3 <- full_join(covariates_FAdates2, CTD_FAdates)  # There were multiple casts per day which is adding rows
covariates_FAdates4 <- full_join(covariates_FAdates3, Rain_FAdates)
covariates_FAdates5 <- full_join(covariates_FAdates4, WindSp_FAdates)
covariates_FAdates6 <- full_join(covariates_FAdates5, WindDir_FAdates)
covariates_FAdates <- full_join(covariates_FAdates6, taxa_FAdates)

# write.csv(covariates_FAdates, "Covariates FA dates redone 7.28.21.csv", row.names = F)


# covariates all dates (I am including the weather sation vars I only calculated for FA dates)
covariates_all1 <- full_join(POMisoData_all, chl.phaeo_all)
covariates_all2 <- full_join(covariates_all1, nutrients_all)
covariates_all3 <- full_join(covariates_all2, CTD_all)
covariates_all4 <- full_join(covariates_all3, Rain_FAdates)
covariates_all5 <- full_join(covariates_all4, WindSp_FAdates)
covariates_all6 <- full_join(covariates_all5, WindDir_FAdates)
covariates_all <- full_join(covariates_all6, taxaSummary)

# write.csv(covariates_all, "Covariates all dates redone 7.28.21.csv", row.names = F)




# Write Data package covariates -------------------------------------------------


# All FA as percentages, total FA in all units
POMpackageCovariates <- covariates_FAdates %>% select(Date, C_N_notacid, delta15n.notacid, delta13c.acidified, 
                                                      C_ug.L.acidified, N_ug.L.notacid,
                                                      chl_20um, chl_GF.F, chl_3um, phaeo.chla, 
                                                      NO3, PO4, Si, SST, SSS, 
                                                      log10.N2, Rain2Day, WindSp2Day, WindDir2Day, 
                                                      Bacillariophyta, Chlorophyta, Ciliophora, Cryptophyta, Dinoflagellata)

colnames(POMpackageCovariates)[c(2:6)] <- c("C.N", "delta15n",   "delta13c", "POC_ug.L", "PN_ug.L")


#write.csv(POMpackageCovariates, "POMFA data package Covariates 20210728.csv", row.names = F)





