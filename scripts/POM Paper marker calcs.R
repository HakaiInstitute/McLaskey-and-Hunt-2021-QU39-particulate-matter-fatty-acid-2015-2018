# This script combines the POMFA data table with the assembled covariates,
# calculates fatty acid concentrations, FA trophic markers, and nutritional
# indices used in McLaskey et al. in prep

library(tidyverse)
library(Hmisc)


POMFAdata <- read.csv("POMFA Data Files to Share/POMFA data package 20210728.csv", stringsAsFactors = FALSE, na=c("", "NA"))
colnames(POMFAdata)



# *Calculate DW concentrations --------------------------------------------

# Multiply percentage by Total (ug) then divide by sample weight (converted to mg)
# Yields FA concentrations in ug/mg or mg/g DW
columns <- c(5:44)
colnames(POMFAdata)[columns]

FAweights <- (data.frame(matrix(ncol = 40, nrow = nrow(POMFAdata))))
for(i in 1:length(columns)){
  FAweights[i] <- ((POMFAdata[columns[i]]*POMFAdata$TotalFA_ug)/(POMFAdata$Sample.dryweight_g*1000))
}

colnames(FAweights) <- paste(str_remove(colnames(POMFAdata)[columns], "_PERCENT"), "mg.gDW", sep = "_")
FAweights$TotalFA_mg.gDW <- POMFAdata$TotalFA_ug/(POMFAdata$Sample.dryweight_g*1000)

POMFAdata <- cbind(POMFAdata, FAweights)

POMFAdata$TotalFA_ug

# *Calculate WW concentrations -------------------------------------------

# Multiply percentage by Total (ug) then divide by sample weight (converted to mg)
# Yields FA concentrations in ug/mg or mg/g WW

FAWETweights <- (data.frame(matrix(ncol = 40, nrow = nrow(POMFAdata))))
for(i in 1:length(columns)){
  FAWETweights[i] <- ((POMFAdata[columns[i]]*POMFAdata$TotalFA_ug)/(POMFAdata$Sample.wetweight_g*1000))
}

colnames(FAWETweights) <- paste(str_remove(colnames(POMFAdata)[columns], "_PERCENT"), "mg.g_WetWt", sep = "_")
FAWETweights$TotalFA_mg.g_WetWt <- POMFAdata$TotalFA_ug/(POMFAdata$Sample.wetweight_g*1000)

POMFAdata <- cbind(POMFAdata, FAWETweights)



# *Calculate /L concentrations --------------------------------------------

# Multiply percentage by Total (ug) then divide by volume filtered (converted to L)
# Yields FA concentrations in ug/L
FAliter <- (data.frame(matrix(ncol = 40, nrow = nrow(POMFAdata))))
for(i in 1:length(columns)){
  FAliter[i] <- (POMFAdata[columns[i]]*POMFAdata$TotalFA_ug)/(POMFAdata$Volume.filtered_mL/1000)
}

colnames(FAliter) <- paste(str_remove(colnames(POMFAdata)[columns], "_PERCENT"), "ug.L", sep = "_")
FAliter$TotalFA_ug.L <- POMFAdata$TotalFA_ug/(POMFAdata$Volume.filtered_mL/1000)

POMFAdata <- cbind(POMFAdata, FAliter)



# Calculate Markers -------------------------------------------------------

# 16:1n7/16:0
POMFAdata$Ratio16.1 <- (POMFAdata$C16.1n.7_PERCENT/ POMFAdata$C16.0_PERCENT)

# Diatom FAs / Flagellate FAs
POMFAdata$Diatom.Flag <- ((POMFAdata$C16.1n.7_PERCENT + POMFAdata$C20.5n.3_PERCENT +  POMFAdata$C16.2n.4_PERCENT + POMFAdata$C16.3n.4_PERCENT) 
                          / (POMFAdata$C22.6n.3_PERCENT + POMFAdata$C18.3n.3_PERCENT + POMFAdata$C18.3n.6_PERCENT + POMFAdata$C18.4n.3_PERCENT))

# DHA/EPA : Dinos/Diatoms
POMFAdata$DHA.EPA <- (POMFAdata$C22.6n.3_PERCENT/ POMFAdata$C20.5n.3_PERCENT)


# 18:2(n-6)+ 18:3(n-3) > 2.5 Terrestrial signal
POMFAdata$Terr_18.2n6_18.3n3 <- (POMFAdata$C18.2n.6c_PERCENT + POMFAdata$C18.3n.3_PERCENT)


# 22:0 and 24:0 indicate terrestrial input or riverine detritus
POMFAdata$River.C22C24 <- (POMFAdata$C22.0_PERCENT + POMFAdata$C24.0_PERCENT)


# 15:0 and 17:0 : Bacteria
POMFAdata$Bacteria_15_17 <- (POMFAdata$C15.0_PERCENT + POMFAdata$C17.0_PERCENT)


# 18:1n-9/18:1n-7 Carnivory or omnivory
POMFAdata$Carn18.1N9_n7 <- (POMFAdata$C18.1n.9c_PERCENT/ POMFAdata$C18.1n.7_PERCENT)


# Make 18:3n-3_18:4n-4 (pico-Chl)
POMFAdata$C18.3andC18.4 <- POMFAdata$C18.3n.3_PERCENT + POMFAdata$C18.4n.3_PERCENT



# Calculate proportions for each main class

percent.SFA <- vector(length=nrow(POMFAdata))
for(i in 1:nrow(POMFAdata)){
  percent.SFA[i]=sum(
    POMFAdata$C10.0_PERCENT[i], 
    POMFAdata$C11.0_PERCENT[i], 
    POMFAdata$C12.0_PERCENT[i], 
    POMFAdata$C13.0_PERCENT[i], 
    POMFAdata$C14.0_PERCENT[i], 
    POMFAdata$C15.0_PERCENT[i], 
    POMFAdata$C16.0_PERCENT[i], 
    POMFAdata$C17.0_PERCENT[i], 
    POMFAdata$C18.0_PERCENT[i], 
    POMFAdata$C20.0_PERCENT[i],
    POMFAdata$C22.0_PERCENT[i],
    POMFAdata$C23.0_PERCENT[i],
    POMFAdata$C24.0_PERCENT[i],
    na.rm=TRUE)
}
POMFAdata$percent.SFA <- percent.SFA


percent.MUFA <- vector(length=nrow(POMFAdata))
for(i in 1:nrow(POMFAdata)){
  percent.MUFA[i]=sum(POMFAdata$C14.1_PERCENT[i], 
                      POMFAdata$C16.1n.7_PERCENT[i], 
                      POMFAdata$C18.1n.7_PERCENT[i], 
                      POMFAdata$C18.1n.9c_PERCENT[i], 
                      POMFAdata$C20.1n.9_PERCENT[i],
                      POMFAdata$C22.1n.9_PERCENT[i],
                      POMFAdata$C24.1n.9_PERCENT[i],na.rm=TRUE)
}
POMFAdata$percent.MUFA <- percent.MUFA


percent.PUFA <- vector(length=nrow(POMFAdata))
for(i in 1:nrow(POMFAdata)){
  percent.PUFA[i]=sum(POMFAdata$C18.2n.6c_PERCENT[i], 
                      POMFAdata$C18.3n.3_PERCENT[i], 
                      POMFAdata$C18.3n.6_PERCENT[i], 
                      POMFAdata$C18.4n.3_PERCENT[i],
                      POMFAdata$C16.2n.4_PERCENT[i],
                      POMFAdata$C16.3n.4_PERCENT[i], 
                      POMFAdata$C20.3n.3_PERCENT[i], 
                      POMFAdata$C20.4n.6_PERCENT[i], 
                      POMFAdata$C20.5n.3_PERCENT[i],
                      POMFAdata$C22.2n.6_PERCENT[i], 
                      POMFAdata$C22.4n.6_PERCENT[i], 
                      POMFAdata$C22.5n.3_PERCENT[i], 
                      POMFAdata$C22.5n.6._PERCENT[i],
                      POMFAdata$C22.6n.3_PERCENT[i], na.rm=TRUE)
} 
POMFAdata$percent.PUFA <- percent.PUFA



# PUFA/SFA carnivory index
POMFAdata$PUFA.SFA <- POMFAdata$percent.PUFA/POMFAdata$percent.SFA



POMFAdata <- mutate(POMFAdata, PUFA_mg.g = (percent.PUFA*TotalFA_ug)/Sample.dryweight_g)
POMFAdata <- POMFAdata %>% mutate(PUFA.L = (percent.PUFA*TotalFA_ug)/(Volume.filtered_mL/1000))


percent.n3_PUFA <- vector(length=nrow(POMFAdata))
for(i in 1:nrow(POMFAdata)){
  percent.n3_PUFA[i]=sum(POMFAdata$C18.3n.3_PERCENT[i], 
                         POMFAdata$C18.4n.3_PERCENT[i],
                         POMFAdata$C20.3n.3_PERCENT[i], 
                         POMFAdata$C20.5n.3_PERCENT[i],
                         POMFAdata$C22.5n.3_PERCENT[i], 
                         POMFAdata$C22.6n.3_PERCENT[i], na.rm=TRUE)
} 
POMFAdata$percent.n3_PUFA <- percent.n3_PUFA

percent.n6_PUFA <- vector(length=nrow(POMFAdata))
for(i in 1:nrow(POMFAdata)){
  percent.n6_PUFA[i]=sum(POMFAdata$C18.2n.6_PERCENT[i], 
                         POMFAdata$C18.3n.6_PERCENT[i], 
                         POMFAdata$C20.4n.6_PERCENT[i], 
                         POMFAdata$C22.2n.6_PERCENT[i], 
                         POMFAdata$C22.4n.6_PERCENT[i], 
                         POMFAdata$C22.5n.6._PERCENT[i], na.rm=TRUE)
} 
POMFAdata$percent.n6_PUFA <- percent.n6_PUFA


# n-3 : n-6 ratio (too many samples w zero n-6)
POMFAdata <- POMFAdata %>% mutate(n3.n6 = percent.n3_PUFA/percent.n6_PUFA)

# EFA = 18C PUFA and longer; HUFA = 20C PUFA and longer
POMFAdata <- POMFAdata %>% mutate(EFA.ug.L = ((C18.3n.6_PERCENT + C18.3n.3_PERCENT + C18.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT)* TotalFA_ug)/(Volume.filtered_mL/1000))
POMFAdata <- POMFAdata %>% mutate(HUFA.ug.L = ((C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT)* TotalFA_ug)/(Volume.filtered_mL/1000))

POMFAdata <- POMFAdata %>% mutate(Prop.EFA = (C18.3n.6_PERCENT + C18.3n.3_PERCENT + C18.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT))
POMFAdata <- POMFAdata %>% mutate(Prop.HUFA = (C20.5n.3_PERCENT + C22.6n.3_PERCENT + C20.4n.6_PERCENT))



# Combine w covariates ----------------------------------------------------

covariatedata <- read.csv("POMFA Data Files to Share/POMFA data package Covariates 20210728.csv", stringsAsFactors = FALSE, na=c("", "NA"))
colnames(POMFAdata)
colnames(covariatedata)

allData <- full_join(POMFAdata, covariatedata)


#write.csv(allData, "POMFA Data Files to Share/POM Paper data table Final 7.28.21.csv", row.names = F)




