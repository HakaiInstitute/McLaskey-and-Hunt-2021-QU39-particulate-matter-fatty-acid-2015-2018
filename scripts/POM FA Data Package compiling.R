# POM FA data compiling 

# 2015-2018 POM samples from QU39 for data package and publication


# POM FA metadata compiling for data package

# 1) First DFO metadata and Hakai metadata are combined, data are checked, 
#    and combined sample info (from filters run together) is generated. 
# 2) Combined metadata are matched to raw GC data.  

require(tidyverse)


# 1) Metadata: Combine and checks -----------------------------------------

# Load Data files ---------------------------------------------------------

# Load DFO metadata (sample weights)
DFOmetadataPOM <- read.csv("POMFA Metadata for Data Package.csv", stringsAsFactors = FALSE, na="")
DFOmetadataPOM$Date.new <- as.Date(DFOmetadataPOM$Date, "%m/%d/%y")

# Load Hakai Metadata (date, station, filter weights, etc.)
HakaiFAMetadataPOM <- read.csv("2020-08-05_HakaiData_pomfas.csv", stringsAsFactors = FALSE)
HakaiFAMetadataPOM$date <- as.Date(HakaiFAMetadataPOM$date, "%Y-%m-%d")


# Bring in filter weights for IDs listed on DFO metadata ------------------

# Filter codes currently in columns 8:12
# I am specifying this manually so double check
colnames(DFOmetadataPOM)

# This matches the filter IDs from the DFOmetadata with the Hakai metadata to extract:
# 1) the number of filters run, 2) the combined weight of those filters, and 3) the combined volume filtered onto all

# Empty data frame to populate with new data
POMfilterweights <- (data.frame(matrix(ncol = 3, nrow = nrow(DFOmetadataPOM))))

for(i in 1:nrow(DFOmetadataPOM)){
  # Number of filters per sample
  POMfilterweights[i,1] <- length(which(!is.na(DFOmetadataPOM[i,8:12])))
  
# pre weights for all filters
  weight1<- HakaiFAMetadataPOM$pre_weight_mg[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,8]]
  if(!is.na(DFOmetadataPOM[i,9])){
    weight2 <- HakaiFAMetadataPOM$pre_weight_mg[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,9]]}
  else{weight2 <- 0}
  if(!is.na(DFOmetadataPOM[i,10])){
    weight3 <- HakaiFAMetadataPOM$pre_weight_mg[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,10]]}
  else{weight3 <- 0}
  if(!is.na(DFOmetadataPOM[i,11])){
    weight4 <- HakaiFAMetadataPOM$pre_weight_mg[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,11]]}
  else{weight4 <- 0}
  if(!is.na(DFOmetadataPOM[i,12])){
    weight5 <- HakaiFAMetadataPOM$pre_weight_mg[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,12]]}
  else{weight5 <- 0}
  # Sum of all filters
  POMfilterweights[i,2] <- sum(weight1, weight2, weight3, weight4, weight5)
  
# volume filtered for all filters
  vol1<- HakaiFAMetadataPOM$volume[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,8]]
  if(!is.na(DFOmetadataPOM[i,9])){
    vol2 <- HakaiFAMetadataPOM$volume[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,9]]}
  else{vol2 <- 0}
  if(!is.na(DFOmetadataPOM[i,10])){
    vol3 <- HakaiFAMetadataPOM$volume[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,10]]}
  else{vol3 <- 0}
  if(!is.na(DFOmetadataPOM[i,11])){
    vol4 <- HakaiFAMetadataPOM$volume[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,11]]}
  else{vol4 <- 0}
  if(!is.na(DFOmetadataPOM[i,12])){
    vol5 <- HakaiFAMetadataPOM$volume[HakaiFAMetadataPOM$hakai_id==DFOmetadataPOM[i,12]]}
  else{vol5 <- 0}
  # Sum of all filters
  POMfilterweights[i,3] <- sum(vol1, vol2, vol3, vol4, vol5)
}
colnames(POMfilterweights) <- c("Num.filters.analyzed","Filter.wts.analyzed", "volume.analyzed")
# Combine with full DFO metadata file
DFOmetadataPOM <- cbind(DFOmetadataPOM, POMfilterweights)



# Check filters that were combined ----------------------------------------

# Generate a list of filter IDs from the Hakai metadata that are from the same date, 
# and either should have been run together or will be combined post hoc, to compare to what was done. 

# Make reduced size dataset to work with of just date and Hakai ID
QU39POMFAsamples <- HakaiFAMetadataPOM[,c(4,12)]
# Remove NAs from strange formating in this version of the Hakai data
QU39POMFAsamples <- QU39POMFAsamples[complete.cases(QU39POMFAsamples),]

# List of Unique sampling dates
dates <- unique(QU39POMFAsamples$date)

# Empty data frame to fill 
POMfilters2combine <- setNames(data.frame(matrix(ncol = 6, nrow = 100)), c("Date", "hakai_id", "hakai_id", "hakai_id", "hakai_id", "hakai_id"))
# For each unique sampling date, pull together all filter numbers
for(i in 1:length(dates)){
  POMfilters2combine[i,1] <- dates[i]
  values <- QU39POMFAsamples$hakai_id[QU39POMFAsamples$date==dates[i]]
  POMfilters2combine[i,2:(1+length(values))] <- values
}
POMfilters2combine$Date <- as.Date.numeric(POMfilters2combine$Date, origin = "1970-01-01")

# This is the file to compare what was actually run to


# Were all filters run? --------------------------------------------------------

# For each ID in POMfilters2combine, is it found in the DFOmetadata POM file?
# DFOmetadataPOM Filter codes currently in columns 8:12
POMfiltersRun <- (data.frame(matrix(ncol = 7, nrow = nrow(POMfilters2combine))))

for(i in 1:nrow(POMfilters2combine)){
  POMfiltersRun[i,1] <- POMfilters2combine[i,1]
  POMfiltersRun[i,2] <- any(grepl(POMfilters2combine[i,2], DFOmetadataPOM[DFOmetadataPOM$Date.new==POMfilters2combine[i,1],c(8:12)]))
  POMfiltersRun[i,3] <- any(grepl(POMfilters2combine[i,3], DFOmetadataPOM[DFOmetadataPOM$Date.new==POMfilters2combine[i,1],c(8:12)]))
  POMfiltersRun[i,4] <- any(grepl(POMfilters2combine[i,4], DFOmetadataPOM[DFOmetadataPOM$Date.new==POMfilters2combine[i,1],c(8:12)]))
  POMfiltersRun[i,5] <- any(grepl(POMfilters2combine[i,5], DFOmetadataPOM[DFOmetadataPOM$Date.new==POMfilters2combine[i,1],c(8:12)]))
  POMfiltersRun[i,6] <- any(grepl(POMfilters2combine[i,6], DFOmetadataPOM[DFOmetadataPOM$Date.new==POMfilters2combine[i,1],c(8:12)]))
  POMfiltersRun[i,7] <- all(POMfiltersRun[i,c(2:6)]==TRUE, na.rm=T)
}
colnames(POMfiltersRun) <- c(colnames(POMfilters2combine[1:6]),"All.Run")

# Combine with main data file
POMfilters2combine <- cbind(POMfilters2combine, POMfiltersRun)

# This is searching in the whole database whether they've been run, 
# so includes some samples that were run separate and need to be combined

# Bring over info on whether all filters have been run to the DFO metadata file- matching based on date
allfilters <- (data.frame(matrix(ncol = 12, nrow = nrow(DFOmetadataPOM))))
for(i in 1:nrow(DFOmetadataPOM)){
  if(any(grepl(DFOmetadataPOM$Date.new[i], POMfilters2combine$Date)==TRUE)){
    allfilters[i,1:12] <- POMfilters2combine[POMfilters2combine$Date==DFOmetadataPOM$Date.new[i],c(1:6,8:13)]
  } else{allfilters[i,1:12] <- NA}
  
}
colnames(allfilters) <- colnames(POMfilters2combine[c(1:6,8:13)])

DFOmetadataPOM <- cbind(DFOmetadataPOM, allfilters)


# Filters not combined but should have? ----------------------

# search for filters that should have been combined but weren't
# DFOmetadataPOM Filter codes currently in columns 8:12
# Hakai_ids are in 30:34

# I want to make sure filter.IDs are the same as hakai_ids for each row
missing.filters <- (data.frame(matrix(ncol = 1, nrow = nrow(DFOmetadataPOM))))
for(i in 1:nrow(DFOmetadataPOM)){
  hakaiID1 <- any(grepl(DFOmetadataPOM[i,30], DFOmetadataPOM[i,8:12]))  
  hakaiID2 <- any(grepl(DFOmetadataPOM[i,31], DFOmetadataPOM[i,8:12]))  
  hakaiID3 <- any(grepl(DFOmetadataPOM[i,32], DFOmetadataPOM[i,8:12]))  
  hakaiID4 <- any(grepl(DFOmetadataPOM[i,33], DFOmetadataPOM[i,8:12]))  
  hakaiID5 <- any(grepl(DFOmetadataPOM[i,34], DFOmetadataPOM[i,8:12]))  
  missing.filters[i,1] <- any(c(hakaiID1, hakaiID2, hakaiID3, hakaiID4, hakaiID5)==FALSE, na.rm = T)
}

DFOmetadataPOM <- cbind(DFOmetadataPOM, missing.filters)
colnames(DFOmetadataPOM)[41] <- "Filters.not.included"



# Extra filters included? -------------------------------------------------

# Then I want to make sure hakai_ids are the same as filter.IDs for each row
extra.filters <- (data.frame(matrix(ncol = 1, nrow = nrow(DFOmetadataPOM))))
for(i in 1:nrow(DFOmetadataPOM)){
  filterID1 <- any(grepl(DFOmetadataPOM[i,8], DFOmetadataPOM[i,30:34]))  
  filterID2 <- any(grepl(DFOmetadataPOM[i,9], DFOmetadataPOM[i,30:34]))  
  filterID3 <- any(grepl(DFOmetadataPOM[i,10], DFOmetadataPOM[i,30:34]))  
  filterID4 <- any(grepl(DFOmetadataPOM[i,11], DFOmetadataPOM[i,30:34]))  
  filterID5 <- any(grepl(DFOmetadataPOM[i,12], DFOmetadataPOM[i,30:34]))  
  extra.filters[i,1] <- any(c(filterID1, filterID2, filterID3, filterID4, filterID5)==FALSE, na.rm = T)
}

DFOmetadataPOM <- cbind(DFOmetadataPOM, extra.filters)
colnames(DFOmetadataPOM)[42] <- "Extra.filters"



# Calculate sample weight -------------------------------------------------

### Make column of sample weights
Sample.wt.g <- (data.frame(matrix(ncol = 1, nrow = nrow(DFOmetadataPOM))))

for(i in 1:nrow(DFOmetadataPOM)){
  Sample.wt.g[i,1] <- DFOmetadataPOM$Dry.Weight..g.[i] - (DFOmetadataPOM$Filter.wts.analyzed[i]/1000)
}

DFOmetadataPOM <- cbind(DFOmetadataPOM, Sample.wt.g)
colnames(DFOmetadataPOM)[43] <- "Sample.wt.g"



# Calculate sample Wet Weight
Sample.WETwt.g <- (data.frame(matrix(ncol = 1, nrow = nrow(DFOmetadataPOM))))

for(i in 1:nrow(DFOmetadataPOM)){
  Sample.WETwt.g[i,1] <- DFOmetadataPOM$Wet.Disk.s..weight..g.[i] - (DFOmetadataPOM$Filter.wts.analyzed[i]/1000)
}

DFOmetadataPOM <- cbind(DFOmetadataPOM, Sample.WETwt.g)
colnames(DFOmetadataPOM)[44] <- "Sample.WETwt.g"


# The 2019 samples don't have filter weights so are not quantitative based on weight 
# But they are not included in this data package
plot(DFOmetadataPOM$Filter.wts.analyzed)



# Write combined metadata file --------------------------------------------
# write.csv(DFOmetadataPOM, "DFO metadata w filter info 1.4.21.csv", row.names = F, na="")


# Remaining quality decisions I made: 

# 5/5/2015 (P82) was run with three filters, but only two were collected - no record of QF771. 
# I am assuming that this is a record keeping issue and it was collected on the date it was labeled with. (I.e. I am leaving it in)

# 11/16/2018 (J319-P53) is missing a filter QF6154 (not run) I haven't found this filter so I am including it as is. 



# 2) GC data matched to combined metadata  -------------------------------------------------


# Combined Metadata file
str(DFOmetadataPOM)

# Load FA data from DFO: output from summarizer
# Add a column w date of run to look for patterns in analysis (and check I'm using the correct run for each sample)

FAData2020_02_19 <- read.csv("Raw GC Data/20200309_153433_SUMMARY 2020-02-19.csv", stringsAsFactors = FALSE, na="")
FAData2020_02_19[,2] <- rep("2020_02_19", nrow(FAData2020_02_19))
FAData2020_02_18 <- read.csv("Raw GC Data/20200221_135408_SUMMARY 2020-02-18.csv", stringsAsFactors = FALSE, na="")
FAData2020_02_18[,2] <- rep("2020_02_18", nrow(FAData2020_02_18))

FAData2019_12_05<- read.csv("Raw GC Data/20200221_161302_SUMMARY 2019-12-05 min height 75.csv", stringsAsFactors = FALSE, na="")
FAData2019_12_05[,2] <- rep("2019_12_05", nrow(FAData2019_12_05))
# FAData2019_12_03<- read.csv("Raw GC Data/20200221_161313_SUMMARY 2019-12-03 min height 75 original.csv", stringsAsFactors = FALSE, na="")
FAData2019_12_03<- read.csv("Raw GC Data/20200221_161313_SUMMARY 2019-12-03 min height 75.csv", stringsAsFactors = FALSE, na="")
FAData2019_12_03[,2] <- rep("2019_12_03", nrow(FAData2019_12_03))

FAData2018_06_26<- read.csv("Raw GC Data/20191217_172155_SUMMARY_2018-06-26.csv", stringsAsFactors = FALSE, na="")
FAData2018_06_26[,2] <- rep("2018_06_26", nrow(FAData2018_06_26))
FAData2018_06_11<- read.csv("Raw GC Data/20191217_172628_SUMMARY_2018-06-11.csv", stringsAsFactors = FALSE, na="")
FAData2018_06_11[,2] <- rep("2018_06_11", nrow(FAData2018_06_11))
FAData2018_06_04<- read.csv("Raw GC Data/20191217_171543_SUMMARY_2018-06-04.csv", stringsAsFactors = FALSE, na="")
FAData2018_06_04[,2] <- rep("2018_06_04", nrow(FAData2018_06_04))
FAData2018_05_30<- read.csv("Raw GC Data/20191219_144434_SUMMARY_2018-05-30.csv", stringsAsFactors = FALSE, na="")
FAData2018_05_30[,2] <- rep("2018_05_30", nrow(FAData2018_05_30))
FAData2018_03_23<- read.csv("Raw GC Data/20191219_112055_SUMMARY_2018-03-23.csv", stringsAsFactors = FALSE, na="")
FAData2018_03_23[,2] <- rep("2018_03_23", nrow(FAData2018_03_23))


# When I bind these in reverse order they search and select for newer runs over older runs, if there are doubles for samples
FAData1 <- full_join(FAData2020_02_19, FAData2020_02_18)
FAData2 <- full_join(FAData1, FAData2019_12_05)
FAData3 <- full_join(FAData2, FAData2019_12_03)
FAData4<- full_join(FAData3, FAData2018_06_26)
FAData5 <- full_join(FAData4, FAData2018_06_11)
FAData6 <- full_join(FAData5, FAData2018_06_04)
FAData7 <- full_join(FAData6, FAData2018_05_30)
FADataRaw <- full_join(FAData7, FAData2018_03_23)


# Clean up peak names
colnames(FADataRaw)[5:64] <- str_remove(colnames(FADataRaw)[5:64], "X1_")



# Add Collection Data -----------------------------------------------------

### For each sample in the DFO metadata, match up output from GC

# Make vector of empty values for samples that are on the metadata list but not in these data files
FAdataEmpty <- (data.frame(matrix(ncol = ncol(FADataRaw), nrow = 1, data= NA)))

# Empty data frame to populate with new data
FAdataTemp <- (data.frame(matrix(ncol = ncol(DFOmetadataPOM)+ncol(FADataRaw), nrow = nrow(DFOmetadataPOM))))
# Search for sample names from metadata in the file path names from data and line up. If no match use empty data vector
for(i in 1:nrow(DFOmetadataPOM)){
  if(any(grepl(DFOmetadataPOM[i,6], FADataRaw[,1]))==TRUE){
    row <- grep(DFOmetadataPOM[i,6], FADataRaw[,1])
    FAdataTemp[i,] <- data.frame(DFOmetadataPOM[i,], FADataRaw[row,])}
  else(FAdataTemp[i,] <- data.frame(DFOmetadataPOM[i,], FAdataEmpty))
}
colnames(FAdataTemp) <- c(colnames(DFOmetadataPOM) , colnames(FADataRaw))
# Warnings are for samples that have an older run the data as well. It takes the first one only. Run date from GC should match date of best run from metadata. 

# Do manual check of run matching with best day
# write.csv(FAdataTemp, "FA compiled w metadata 1.4.21.csv", row.names = F)

# This file still includes samples that I don't have good runs for. Remove them
FAdataTemp <- FAdataTemp[!is.na(FAdataTemp$file),]



# Combine filters from same sample --------------------------------------

# Calculate Total FA Area without the standard C19:0 (which was added to samples)
FAdataTemp$AREA.Total.STD <- (FAdataTemp$Total.FA_G_AREA - FAdataTemp$C19.0_AREA)

# multiple columns called date because I'm not saving the intermediate file anymore
colnames(FAdataTemp)[29] <- "date.2"
# Sum the filters that need to be combined (i.e. were the same sample split across multiple filters)
# To combine filters that are from the same sample, remove cols that do not have unique fields
POMFAsmall <- FAdataTemp %>% select(Date, volume.analyzed, Sample.wt.g, Sample.WETwt.g, IVOL, C04.0_AREA:AREA.Total.STD)

POMFAsmall$Date <- as.Date(FAdataTemp$Date.new, origin="1970-01-01")
POMFAsmall = aggregate(POMFAsmall[,c(2:66)],
                       by = list(POMFAsmall$Date),
                       FUN = sum, na.rm=TRUE)
colnames(POMFAsmall)[1] <- "Date"
# 138 rows instead off 144: that is correct, six samples were combined (plus one that one filter had bad data; one that is still missing its other filter)


# This is to track the date samples were first run- an indication of when they were extracted.
# Can't include this when I aggregate
# I need this to do a manual correction on C17:0
rundates <- FAdataTemp %>% select(Date, Date.of.first.run)
rundates$Date.of.first.run <- as.Date(rundates$Date.of.first.run, "%m/%d/%y")
rundates$Date <- as.Date(rundates$Date, "%m/%d/%y")

# Run dates have not been aggregated but POMFAsmall has
rundates <- rundates[!duplicated(rundates$Date),]

POMFAsmall <- left_join(POMFAsmall, rundates)

# Replace the IVOL column with Date.of.first.run
POMFAsmall$IVOL <- POMFAsmall$Date.of.first.run
POMFAsmall <- POMFAsmall[,-67]
colnames(POMFAsmall)[5] <- "Date.of.first.run"


# Clean and Standardize Data to Analyze -----------------------------------

# Make a new column with the unstadardized area before standarizing the original one
POMFAsmall$AREA.Total.STD_unstandardized <- POMFAsmall$AREA.Total.STD

# To calculate markers from Areas, I want to replace all NA values with zero in columns 6:66 (the FA peaks)
POMFAsmall[6:66][is.na(POMFAsmall[6:66])] <- 0

# standardize peak areas to C19:0 peak area = 500 ug 
hist(POMFAsmall$C19.0_AREA,50)

# standardize based on C19:0 area in each sample
# 47 to 92 because I need to standardize the total as well if that's what Im going to use to calc percentages later
for(i in c(6:35,37:66)){
  POMFAsmall[,i] <- ((POMFAsmall[i]/POMFAsmall$C19.0_AREA)*500)
}



# C17:0 correction --------------------------------------------------------

# There is a change in the minimum standardized quantity of C17:0 in each sample 
# that suggests a minute contamination in the C19:0 standard that is added, that was more
# apparent during the first time period samples were run. I am applying a manual correction of this. 

# now that the data are standardized I need to do the correction on 17:0 to account for the C19 standard
POMFAsmall$Date.of.first.run <- as.character(POMFAsmall$Date.of.first.run)

# Split the dataset into two based on the date samples were run 
POMFAdata.early <- POMFAsmall %>% filter(Date.of.first.run %in% c("2018-03-24", "2018-05-30", "2018-06-27", "2018-06-12", "2018-06-05"))
POMFAdata.late <- POMFAsmall %>% filter(Date.of.first.run %in% c("2019-10-24", "2019-04-16", "2019-06-10", "2019-06-06"))


# I have decided that the minimum area for 17:0 from the first half = 1.65 and second half = 0.9 
# I am going to subtract 0.75 from each run in that first set before I aggregate the runs. 
POMFAdata.early <- POMFAdata.early %>% mutate(C17.0_AREA.modified = C17.0_AREA-0.75)
POMFAdata.early$timeperiod <- "early"
# Am not modifying 17:0 in samples processed later
POMFAdata.late <- POMFAdata.late %>% mutate(C17.0_AREA.modified = C17.0_AREA)
POMFAdata.late$timeperiod <- "late"

POMFAdata.modified <- full_join(POMFAdata.early, POMFAdata.late)


POMFAdata.modified.sm <- POMFAdata.modified %>% select(Date, C17.0_AREA.modified)

POMFAsmall <- full_join(POMFAsmall, POMFAdata.modified.sm)

# Visual check
plot(POMFAdata.modified$C17.0_AREA ~ POMFAdata.modified$Date, 
     col=as.factor(POMFAdata.modified$timeperiod), pch=16, ylab="C17:0 μg sample", xlab="", yaxs="i")
plot(POMFAdata.modified$C17.0_AREA.modified ~ POMFAdata.modified$Date, 
     col=as.factor(POMFAdata.modified$timeperiod), pch=16, ylab="C17:0 μg sample", xlab="", yaxs="i")
#This brings the older samples in line with the newer ones. It is the best correction I can do. 


png("C17 pre-correction.png", width=5, height=4, units="in", res=150)
plot(POMFAsmall$C17.0_AREA ~ POMFAsmall$Date, 
     col=as.factor(POMFAsmall$timeperiod), pch=16, ylab="C17:0 μg/sample", xlab="", yaxs="i")
dev.off()

png("C17 with correction.png", width=5, height=4, units="in", res=150)
plot(POMFAsmall$C17.0_AREA.modified ~ POMFAsmall$Date, 
     col=as.factor(POMFAsmall$timeperiod), pch=16, ylab="C17:0 μg/sample", xlab="", yaxs="i")
dev.off()

# Replace C17:0 w corrected values 
POMFAsmall$C17.0_AREA <- POMFAsmall$C17.0_AREA.modified

 
 
 

# Peaks that are ID'ed correctly ------------------------------------------

# I don't trust 15:1 or 17:1 because they don't line up right
# C21:0 is a contaminant coming in with the C19:0 standard
# Take out C18:5n-3, it is not ID'ed w the standards (I'm not confident in it) and is rare
# 20:2n-6 looks like a bad ID
# 20:3n-6 is late/bad ID
# More notes in Peaks I trust.xlsx and POM FA data peak notes UPTODATE.xlsx

# Still open to question: C10.0_AREA, C11.0_AREA, C14.1_AREA, C18.1n.7t_AREA, C22.2n.6_AREA, C22.4n.6_AREA, C23.0_AREA

# Of those, the only ones I am actually tempted to throw out are C18.1n.7t_AREA (low and a little strange, higher in 2015) 
# and C23.0_AREA (low and erratic). I am leaving them in for now. 

POMFAsmall <- POMFAsmall %>%  select(Date, volume.analyzed, Sample.wt.g, Sample.WETwt.g, Date.of.first.run,
                                     C10.0_AREA,     C11.0_AREA ,       C12.0_AREA ,
                                     C13.0_AREA ,    C14.0_AREA ,       C14.1_AREA ,    C15.0_AREA , 
                                     C16.0_AREA ,    C16.1n.7_AREA,        
                                     C16.1n.9_AREA,  C17.0_AREA ,       C18.0_AREA , 
                                     C18.1n.12_AREA, C18.1n.7_AREA,     C18.1n.7t_AREA, C18.1n.9c_AREA,               
                                     C18.1n.9t_AREA, C18.2n.6c_AREA,    C18.2n.6t_AREA, C18.3n.3_AREA,            
                                     C18.3n.6_AREA,  C18.4n.3.SDA_AREA, 
                                     C20.0_AREA,     C20.1n.9_AREA,               
                                     C20.3n.3_AREA,  C20.4n.6_AREA,           
                                     C20.5n.3_AREA,  C22.0_AREA,        C22.1n.9_AREA,                
                                     C22.2n.6_AREA,  C22.4n.6_AREA,  C22.5n.3_AREA,                 
                                     C22.5n.6._AREA, C22.6n.3_AREA,     C23.0_AREA ,    C24.0_AREA , 
                                     C24.1n.9_AREA,  Diatom.II_AREA,    Diatom.III_AREA, UNK.II_AREA,            
                                     UNK.I_AREA, AREA.Total.STD)





# *Calculate Percentages --------------------------------------------------

# Calculate sum of IDed FAs 
POMFAsmall <- POMFAsmall %>% mutate(SumFA=rowSums(.[c(6:46)]))

# Calculate percentages, but not for Total.FA_G_AREA 
percentages <- (data.frame(matrix(ncol = 41, nrow = nrow(POMFAsmall))))
columns <- c(6:46)
for(i in 1:length(columns)){
  percentages[i] <- ((POMFAsmall[columns[i]]/POMFAsmall$SumFA))
}

colnames(percentages) <- paste0(str_remove(colnames(POMFAsmall)[columns], "_AREA"), "_PERCENT")

POMFAsmall <- cbind(POMFAsmall, percentages)


# *Calculate DW concentrations --------------------------------------------

# Divide standardized area by 1000 to convert from ug to mg, then divide by sample mass (g) to get mg/g
# Need to do this for total FA as well, USE SumFA!
colnames(POMFAsmall)[c(6:46, 48)]
columns <- c(6:46, 48)

FAweights <- (data.frame(matrix(ncol = 42, nrow = nrow(POMFAsmall))))
for(i in 1:length(columns)){
  FAweights[i] <- ((POMFAsmall[columns[i]]/1000)/POMFAsmall$Sample.wt.g)
}

colnames(FAweights) <- paste(str_remove(colnames(POMFAsmall)[columns], "_AREA"), "mg.g", sep = "_")

POMFAsmall <- cbind(POMFAsmall, FAweights)



# *Calculate WW concentrations -------------------------------------------

# Divide standardized area by 1000 to convert from ug to mg, then divide by sample mass (g) to get mg/g
# Need to do this for total FA as well
FAWETweights <- (data.frame(matrix(ncol = 42, nrow = nrow(POMFAsmall))))
for(i in 1:length(columns)){
  FAWETweights[i] <- ((POMFAsmall[columns[i]]/1000)/POMFAsmall$Sample.WETwt.g)
}

colnames(FAWETweights) <- paste(str_remove(colnames(POMFAsmall)[columns], "_AREA"), "mg.g_WetWt", sep = "_")

POMFAsmall <- cbind(POMFAsmall, FAWETweights)




# *Calculate /L concentrations --------------------------------------------

# Divide standardized area by 1000 to convert from ug to mg, then divide by sample mass (g) to get mg/g
# Need to do this for total FA as well
FAliter <- (data.frame(matrix(ncol = 42, nrow = nrow(POMFAsmall))))
for(i in 1:length(columns)){
  FAliter[i] <- (POMFAsmall[columns[i]]/(POMFAsmall$volume.analyzed/1000))
}

colnames(FAliter) <- paste(str_remove(colnames(POMFAsmall)[columns], "_AREA"), "ug.L", sep = "_")

POMFAsmall <- cbind(POMFAsmall, FAliter)










# Write file for Data Package ---------------------------------------------

# All FA as percentages, total FA in all units
POMpackage <- POMFAsmall %>% select(Date:Sample.WETwt.g, C10.0_PERCENT:C24.1n.9_PERCENT, Diatom.III_PERCENT:UNK.I_PERCENT, SumFA)

# Up to this point I'v left some peak names as what they were historically called in the DFO raw data
# Changing to what they actually are
colnames(POMpackage)[c(2:4,17,26,42:45)] <- c("Volume.filtered_mL", "Sample.dryweight_g", "Sample.wetweight_g",
                                        "C16.3n.4_PERCENT", "C18.4n.3_PERCENT", "C16.2n.4_PERCENT", 
                                           "C20.4n.3_PERCENT", "C22.1n.11_PERCENT", "TotalFA_ug")


# Change three lowest dry weights to NA (two negative and one unrealistically low)
POMpackage$Sample.dryweight_g[POMpackage$Sample.dryweight_g<0.002] <- NA


# Subset data to remove 2019 samples
POMpackage <- POMpackage %>% 
  filter(Date<"2019-01-01")


#write.csv(POMpackage, "POMFA data package 20210728.csv", row.names = F)


#write.csv(POMFAsmall, "POMFA full data 20210728.csv", row.names = F)

