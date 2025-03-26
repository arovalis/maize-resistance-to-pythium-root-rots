################################################################################
## Field Experiment Script, Fall 2024  
## Analysis by Harrison Hall (harrisonpiercehall@gmail.com)
################################################################################

## Libraries used 
library("dplyr") #contains various helpful functions
library("tidyverse")   ## Adds ggplot2 and other data manipulation tools
library("ggplot2") ## Some functions look to see if ggplot2 is actually in the library or not
library("ggthemes") ## More control over ggplots with more options
library("esquisse") ## GUI for ggplots to look at data quickly
library("ggcorrplot")## More ggplot options for correlations
library("emmeans") ## Analysis of estimated means
library("multcompView") ## To add letters to plots for multiple comparison
library("magrittr") ## Needed for quick manipulation of ggplots and for chaining commands
library("lmtest") ## For testing linear regression models
library("car") ## Provides advanced regression tools
library("lme4") ## Fits linear and generalized linear mixed-effects models
library("agricolae") ## Statistical tests and graphics 
library("MASS") ## More functions and datasets 
library("readxl") ## To open excel files easily


##### File path setup ##########################################################
ProjDir <- rstudioapi::getActiveDocumentContext()$path # Gets path of script
setwd(dirname(ProjDir)) # Sets the working directory to the same location as the script location
ResultDir <- file.path(dirname(ProjDir), "Results")

dir.create(ResultDir, recursive = TRUE)
fTime <- format(Sys.time(), "%y-%m-%d %H%M") #Include time on each file created
AnalysisDir <- file.path(dirname(ResultDir), "data_analysis", fTime, sep="") #make analysis folder for each time
dir.create(AnalysisDir, recursive = TRUE)
##### Get files and folders#####################################################
dfDir <- rstudioapi::selectFile(caption="Choose dataset to work with", path=ProjDir) #To manually select data
data <- read_excel(dfDir) # Reads if it is a .xlsx file. If it is a csv, change read_excel to read_csv
dataName <- tools::file_path_sans_ext(basename(dfDir)) # Grabs file name without extension

AnalysisDataDir <- paste0(AnalysisDir, "/", dataName, sep="") # Creates a folder for analysis
experiment_name_here <- dataName # Names all future files after the file name. You can manually change this.
dir.create(AnalysisDataDir, recursive = TRUE) #Guarantees the file is made.
# Create paths for anova, graph, and tables directories with a timestamp
anovaDir <- paste0(AnalysisDataDir, "/analysis", sep="")
graphDir <- paste0(AnalysisDataDir, "/visual", sep="")
tablesDir <- paste0(AnalysisDataDir, "/tables", sep="")
dir.create(anovaDir, recursive = TRUE)
dir.create(graphDir, recursive = TRUE)
dir.create(tablesDir, recursive = TRUE)
##### Graph Universal Variables#################################################
options(max.print=100000000) # Guarantees all datapoints will fit in the graphs
fontSize <- 16 # Can be manually changed later
# Set the default exported chart size to 5000x2000 units, this is print quality
ChartWidth <- (5000)
ChartHeight <- (2000)
# Set the size of any dots 
DotSize <- (15)
# Set default fill colors - try to make these colorblind friendly 
fillOne <- "#F0E442"
fillTwo <- "#56B4E9"
fillThree <- "#D55E00"
fillFour <- "#0072B2"
fillFive <- "#E69F00"
fillSix <- "#009E73"
fillSeven <- "#CC79A7"
fillEight <- "#E8E7A2"
fillNine <- "#86DAB2"
fillTen <- "#555555"
# All of the above colors work for the top three most common colorblind types

################################################################################
##### Beginning of analysis ####################################################
################################################################################

## First let's look at our data:

head(data) # Take a peek at the data and make sure we grabbed the right file
str(data) # This shows us our column names and examples of data. We will want to make sure things are formatted the way we want.

data <- data[!grepl("1", data$escape),]

## If it is a category or qualitative (red/blue) trait, we want it to be a factor
## If it is a quantitative trait (length or mass), we want it to be numeric

data$PlotID <- as.factor(data$PlotID)
data$row <- as.factor(data$row)
data$column <- as.factor(data$column)
data$inoculation <- as.factor(data$inoculation)
data$seed <- as.factor(data$seed)

data$stuntCountJune <- as.numeric(data$stuntCountJune)
data$standCount <- as.numeric(data$standCount)
data$earCount <- as.numeric(data$earCount)

data$treatmentStuntCountAvg <- as.numeric(data$treatmentStuntCountAvg)
data$treatmentStandCountAvg <- as.numeric(data$treatmentStandCountAvg)
data$treatmentEarCountAvg <- as.numeric(data$treatmentEarCountAvg)

data$treatmentStuntCountPercent <- as.numeric(data$treatmentStuntCountPercent)*100
data$treatmentStandCountPercent <- as.numeric(data$treatmentStandCountPercent)*100
data$treatmentEarCountPercent <- as.numeric(data$treatmentEarCountPercent)*100

data$lineStuntCountAvg <- as.numeric(data$lineStuntCountAvg)
data$lineStandCountAvg <- as.numeric(data$lineStandCountAvg)
data$lineEarCountAvg <- as.numeric(data$lineEarCountAvg)

data$lineStuntCountPercent <- as.numeric(data$lineStuntCountPercent)*100
data$lineStandCountPercent <- as.numeric(data$lineStandCountPercent)*100
data$lineEarCountPercent <- as.numeric(data$lineEarCountPercent)*100

### Remove outliers:
# This will remove anything that is more than 3 standard deviations above the mean
stats <- data %>%
  group_by(inoculation) %>%
  summarise(
    meanlineStuntCountPercent = mean(lineStuntCountPercent, na.rm = TRUE),
    sdlineStuntCountPercent = sd(lineStuntCountPercent, na.rm = TRUE),
    meanlineStandCountPercent = mean(lineStandCountPercent, na.rm = TRUE),
    sdlineStandCountPercent = sd(lineStandCountPercent, na.rm = TRUE),
    meanlineEarCountPercent = mean(lineEarCountPercent, na.rm = TRUE),
    sdlineEarCountPercent = sd(lineEarCountPercent, na.rm = TRUE)
  )

# Join this stats back to the original data
filteredData <- data %>%
  left_join(stats, by = "inoculation") %>%
  filter(
    lineStuntCountPercent <= meanlineStuntCountPercent + 2.1 * sdlineStuntCountPercent,
    lineStandCountPercent <= meanlineStandCountPercent + 2.1 * sdlineStandCountPercent,
    lineEarCountPercent <= meanlineEarCountPercent + 2.1 * sdlineEarCountPercent
  )

# Select columns to avoid clutter
data_cleaned <- filteredData %>%
  dplyr::select(-starts_with("mean"), -starts_with("sd"))




ogData <- data_cleaned

## Let's open the data in esquisse and see what kind of graphs we can make, and visualize some of the data.
#esquisser(data)

################################################################################################################################################################
# Experiment Analysis - Did the experimental design work?

exp_data <- data_cleaned

AnovaData <- exp_data[!grepl("Mock Inocs", exp_data$inoculation),]
AnovaData <- exp_data[!grepl("Mock Inoc", exp_data$inoculation),]


summary(aov(lineStuntCountPercent ~ seed * inoculation + row, AnovaData))
summary(aov(lineEarCountPercent ~ seed * inoculation + row, AnovaData))
summary(aov(lineStandCountPercent ~ seed * inoculation + row, AnovaData))



anovaFile <- paste0(anovaDir,"/anova_summaries.txt",sep="")

# Use capture.output to save each summary to the file
capture.output({
  cat("Summary of ANOVA for Stunt Count as Percent of Standards:\n")
  print(summary(aov(lineStuntCountPercent ~ seed * inoculation + row, AnovaData)))
  
  cat("\n\nSummary of ANOVA for Stand Count as Percent of Standards:\n")
  print(summary(aov(lineStandCountPercent ~ seed * inoculation + row, AnovaData)))
  
  cat("\n\nSummary of ANOVA for percentMassPedByIsolate:\n")
  print(summary(aov(lineEarCountPercent ~ seed * inoculation + row, AnovaData)))
}, file = anovaFile)




model_lesion <- aov(lineStuntCountPercent ~ inoculation, exp_data)
model_length <- aov(lineStandCountPercent ~ inoculation, exp_data)
model_mass <- aov(lineEarCountPercent ~ inoculation, exp_data)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$inoculation$Letters)
letters_lesionDf$inoculation <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
exp_data <- merge(exp_data, letters_lesionDf, by = "inoculation")
summaryData <- exp_data %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters))

ggplot(exp_data) +
  aes(
    x = inoculation,
    y = lineStuntCountPercent,
    fill = inoculation
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  labs(
    x = "Inoculation Species",  
    y = "Stunt Count as Percent of Standard" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = inoculation, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Treatment Stunt Count.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Stunt Count.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$inoculation$Letters)
letters_lengthDf$inoculation <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
exp_data <- merge(exp_data, letters_lengthDf, by = "inoculation")
summaryData <- exp_data %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent = median(lineStandCountPercent), lengthLetters = first(lengthLetters))

ggplot(exp_data) +
  aes(
    x = inoculation,
    y = lineStandCountPercent,
    fill = inoculation
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = inoculation, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Treatment Stand Count of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Stand Count of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$inoculation$Letters)
letters_massDf$inoculation <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
exp_data <- merge(exp_data, letters_massDf, by = "inoculation")
summaryData <- exp_data %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent = median(lineEarCountPercent), massLetters = first(massLetters))

ggplot(exp_data) +
  aes(
    x = inoculation,
    y = lineEarCountPercent,
    fill = inoculation
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = inoculation, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Treatment Ear Count of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Ear Count of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


######################### IS THERE A DIFFERENCE IN ROOT ROTS? ##################################################################################################

NoControlData <- data[!grepl("Mock Inocs", data$inoculation),]
NoControlData <- data[!grepl("Mock Inoc", data$inoculation),]

model_lesion <- aov(lineStuntCountPercent ~ seed, NoControlData)
model_length <- aov(lineStandCountPercent ~ seed, NoControlData)
model_mass <- aov(lineEarCountPercent ~ seed, NoControlData)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
NoControlData <- merge(NoControlData, letters_lesionDf, by = "seed")
summaryData <- NoControlData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters))

ggplot(NoControlData) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pedigree Root Rot Stunt Count.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Root Rot Stunt Count.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
NoControlData <- merge(NoControlData, letters_lengthDf, by = "seed")
summaryData <- NoControlData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent = median(lineStandCountPercent), lengthLetters = first(lengthLetters))

ggplot(NoControlData) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pedigree Stunt Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Stunt Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
NoControlData <- merge(NoControlData, letters_massDf, by = "seed")
summaryData <- NoControlData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent = median(lineEarCountPercent), massLetters = first(massLetters))

ggplot(NoControlData) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pedigree Ear Count Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Ear Count Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


############ INDIVIDUAL SPECIES ANALYSIS #######################################################################################################################################
fgData <- data[grepl("Fusarium", data$inoculation),]
miData <- data[grepl("Mock Inoc", data$inoculation),]
puuData <- data[grepl("Pythium", data$inoculation),]

##Fusarium only##########################################################################################################################################################

model_lesion <- aov(lineStuntCountPercent ~ seed, fgData)
model_length <- aov(lineStandCountPercent ~ seed, fgData)
model_mass <- aov(lineEarCountPercent ~ seed, fgData)

summary(aov(lineStuntCountPercent ~ seed + row, fgData))
summary(aov(lineStandCountPercent ~ seed + row, fgData))
summary(aov(lineEarCountPercent ~ seed + row, fgData))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
fgData <- merge(fgData, letters_lesionDf, by = "seed")
summaryData <- fgData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters))

ggplot(fgData) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Fusarium Root Rot Stunt Count.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Fusarium Root Rot Stunt Count.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
fgData <- merge(fgData, letters_lengthDf, by = "seed")
summaryData <- fgData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent = median(lineStandCountPercent), lengthLetters = first(lengthLetters))

ggplot(fgData) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Fusarium Stunt Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Fusarium Stunt Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
fgData <- merge(fgData, letters_massDf, by = "seed")
summaryData <- fgData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent = median(lineEarCountPercent), massLetters = first(massLetters))

ggplot(fgData) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Fusarium Ear Count Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Fusarium Ear Count Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##Mock Inoc only##############################################################################################################################################################

model_lesion <- aov(lineStuntCountPercent ~ seed, miData)
model_length <- aov(lineStandCountPercent ~ seed, miData)
model_mass <- aov(lineEarCountPercent ~ seed, miData)

summary(aov(lineStuntCountPercent ~ seed + row, miData))
summary(aov(lineStandCountPercent ~ seed + row, miData))
summary(aov(lineEarCountPercent ~ seed + row, miData))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
miData <- merge(miData, letters_lesionDf, by = "seed")
summaryData <- miData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters))

ggplot(miData) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Mock Inoc Root Rot Stunt Count.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Mock Inoc Root Rot Stunt Count.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
miData <- merge(miData, letters_lengthDf, by = "seed")
summaryData <- miData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent = median(lineStandCountPercent), lengthLetters = first(lengthLetters))

ggplot(miData) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Mock Inoc Stunt Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Mock Inoc Stunt Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
miData <- merge(miData, letters_massDf, by = "seed")
summaryData <- miData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent = median(lineEarCountPercent), massLetters = first(massLetters))

ggplot(miData) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y =  -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Mock Inoc Ear Count Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Mock Inoc Ear Count Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##Pythium only##############################################################################################################################################################

model_lesion <- aov(lineStuntCountPercent ~ seed, puuData)
model_length <- aov(lineStandCountPercent ~ seed, puuData)
model_mass <- aov(lineEarCountPercent ~ seed, puuData)

summary(aov(lineStuntCountPercent ~ seed + row, puuData))
summary(aov(lineStandCountPercent ~ seed + row, puuData))
summary(aov(lineEarCountPercent ~ seed + row, puuData))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
puuData <- merge(puuData, letters_lesionDf, by = "seed")
summaryData <- puuData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters))

ggplot(puuData) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y =  -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pythium Root Rot Stunt Count.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pythium Root Rot Stunt Count.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
puuData <- merge(puuData, letters_lengthDf, by = "seed")
summaryData <- puuData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent = median(lineStandCountPercent), lengthLetters = first(lengthLetters))

ggplot(puuData) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y =  -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pythium Stunt Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pythium Stunt Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
puuData <- merge(puuData, letters_massDf, by = "seed")
summaryData <- puuData %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent = median(lineEarCountPercent), massLetters = first(massLetters))

ggplot(puuData) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y =  -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pythium Ear Count Percent of Standards.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pythium Ear Count Percent of Standards.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

model_lesion <- aov(lineStuntCountPercent ~ seed, puuData)
model_length <- aov(lineStandCountPercent ~ seed, puuData)
model_mass <- aov(lineEarCountPercent ~ seed, puuData)

####################################################################################################################################################################################

seedList <- unique(NoControlData$seed)

# Loop through each seed and create a plot
for(seed in seedList) {
  
  # Filter data for the current seed
  dataSubset <- NoControlData %>% filter(seed == !!seed)
  Sys.sleep(.1)
  print(paste("Processing seed:", seed))
  Sys.sleep(.1)
  
  summary(aov(lineStuntCountPercent ~ inoculation + row, dataSubset))
  summary(aov(lineStandCountPercent ~ inoculation + row, dataSubset))
  summary(aov(lineEarCountPercent ~ inoculation + row, dataSubset))
  
  model_lesion <- aov(lineStuntCountPercent ~ inoculation, dataSubset)
  model_length <- aov(lineStandCountPercent ~ inoculation, dataSubset)
  model_mass <- aov(lineEarCountPercent ~ inoculation, dataSubset)
  
  tukey_lesion <- TukeyHSD(model_lesion)
  tukey_length <- TukeyHSD(model_length)
  tukey_mass <- TukeyHSD(model_mass)
  
  letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
  letters_length <- multcompLetters4(model_length, tukey_length)
  letters_mass <- multcompLetters4(model_mass, tukey_mass)
  
  
  labelName <- paste0("Disease % Stunt Count by line:", seed, sep = " ")
  ylab = labelName
  
  letters_lesionDf <- as.data.frame(letters_lesion$inoculation$Letters)
  
  letters_lesionDf$inoculation <- rownames(letters_lesionDf)
  
  colnames(letters_lesionDf)[1] <- "LesionLetters"
  
  rownames(letters_lesionDf) <- NULL
  
  dataSubset <- merge(dataSubset, letters_lesionDf, by = "inoculation")
  
  letters_lesionDf <- as.data.frame(letters_lesion$inoculation$Letters)
  
  letters_lesionDf$inoculation <- rownames(letters_lesionDf)
  
  colnames(letters_lesionDf)[1] <- "LesionLetters"
  
  rownames(letters_lesionDf) <- NULL
  
  dataSubset <- merge(dataSubset, letters_lesionDf, by = "inoculation")
  
  summaryData <- dataSubset %>%
    group_by(inoculation, seed) %>%
    dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters), .groups = "drop")
  
  
  ggplot(dataSubset) +
    aes(
      x = inoculation,
      y = lineStuntCountPercent,
      fill = inoculation
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey() +
    geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
    labs(
      x = "Inoculum Species",  
      y = "Stunt Count as Percent of Standards" 
    ) +
    # Use the summarized data to plot the letters
    geom_text(
      data = summaryData, 
      aes(
        x = inoculation, 
        y =  -Inf,  # Adjust to place the letters above the boxplot
        label = LesionLetters
      ), 
      position = position_dodge(width = 0.75), 
      vjust = -0.5
    ) +
    labs(y = ylab) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),
    )
  
  filename <- paste0(seed, "_percent Stunt Count.jpg") 
  ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  filename <- paste0(seed, "_percent Stunt Count.wmf")  
  ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  Sys.sleep(.5)
  
  
  labelName <- paste0("% StandCount by Line:", seed, sep = " ")
  ylab = labelName
  
  letters_lengthDf <- as.data.frame(letters_length$inoculation$Letters)
  
  letters_lengthDf$inoculation <- rownames(letters_lengthDf)
  
  colnames(letters_lengthDf)[1] <- "lengthLetters"
  
  rownames(letters_lengthDf) <- NULL
  
  dataSubset <- merge(dataSubset, letters_lengthDf, by = "inoculation")
  
  letters_lengthDf <- as.data.frame(letters_length$inoculation$Letters)
  
  letters_lengthDf$inoculation <- rownames(letters_lengthDf)
  
  colnames(letters_lengthDf)[1] <- "lengthLetters"
  
  rownames(letters_lengthDf) <- NULL
  
  dataSubset <- merge(dataSubset, letters_lengthDf, by = "inoculation")
  
  summaryData <- dataSubset %>%
    group_by(inoculation, seed) %>%
    dplyr::summarize(medianlineStandCountPercent = median(lineStandCountPercent), lengthLetters = first(lengthLetters), .groups = "drop")
  
  
  ggplot(dataSubset) +
    aes(
      x = inoculation,
      y = lineStandCountPercent,
      fill = inoculation
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey() +
    geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
    labs(
      x = "Inoculum Species",  
      y = "Stand Count as Percent of Standards" 
    ) +
    # Use the summarized data to plot the letters
    geom_text(
      data = summaryData, 
      aes(
        x = inoculation, 
        y =  -Inf,  # Adjust to place the letters above the boxplot
        label = lengthLetters
      ), 
      position = position_dodge(width = 0.75), 
      vjust = -0.5
    ) +
    labs(y = ylab) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),
    )
  
  filename <- paste0(seed, "_percent Stand Count.jpg") 
  ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  filename <- paste0(seed, "_percent Stand Count.wmf")  
  ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  Sys.sleep(.5)
  
  labelName <- paste0("% EarCount by Line:", seed, sep = " ")
  ylab = labelName
  
  letters_massDf <- as.data.frame(letters_mass$inoculation$Letters)
  
  letters_massDf$inoculation <- rownames(letters_massDf)
  
  colnames(letters_massDf)[1] <- "massLetters"
  
  rownames(letters_massDf) <- NULL
  
  dataSubset <- merge(dataSubset, letters_massDf, by = "inoculation")
  
  letters_massDf <- as.data.frame(letters_mass$inoculation$Letters)
  
  letters_massDf$inoculation <- rownames(letters_massDf)
  
  colnames(letters_massDf)[1] <- "massLetters"
  
  rownames(letters_massDf) <- NULL
  
  dataSubset <- merge(dataSubset, letters_massDf, by = "inoculation")
  
  summaryData <- dataSubset %>%
    group_by(inoculation, seed) %>%
    dplyr::summarize(medianlineEarCountPercent = median(lineEarCountPercent), massLetters = first(massLetters), .groups = "drop")
  
  
  ggplot(dataSubset) +
    aes(
      x = inoculation,
      y = lineEarCountPercent,
      fill = inoculation
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey() +
    geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
    labs(
      x = "Inoculum Species",  
      y = "Ear Count as Percent of Standards" 
    ) +
    # Use the summarized data to plot the letters
    geom_text(
      data = summaryData, 
      aes(
        x = inoculation, 
        y =  -Inf,  # Adjust to place the letters above the boxplot
        label = massLetters
      ), 
      position = position_dodge(width = 0.75), 
      vjust = -0.5
    ) +
    labs(y = ylab) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),
    )
  
  filename <- paste0(seed, "_percent Ear Count.jpg") 
  ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  filename <- paste0(seed, "_percent Ear Count.wmf")  
  ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  Sys.sleep(1)
  
}

###########################################################################################################################################

###########################################################################################################################################

###########################################################################################################################################


NoControlDeadRemoved <- ogData

histMax <- max(NoControlDeadRemoved$lineStandCountPercent)+10.5
histBy <- max(NoControlDeadRemoved$lineStandCountPercent)/25
hist(NoControlDeadRemoved$lineStandCountPercent, 
     breaks = seq(-0.1, histMax, by = histBy), #Adjust these until you have the fidelity you want
     main = "Percent Lengths", 
     xlab = "", 
     col = "darkblue", 
     border = "white")
# This can be an early indication of normality!

histMax <- max(NoControlDeadRemoved$lineEarCountPercent)+10.5
histBy <- max(NoControlDeadRemoved$lineEarCountPercent)/25
hist(NoControlDeadRemoved$lineEarCountPercent, 
     breaks = seq(-0.1, histMax, by = histBy), #Adjust these until you have the fidelity you want
     main = "Percent Masses", 
     xlab = "", 
     col = "darkred", 
     border = "white")
# This can be an early indication of normality!

##### ASSUMPTION TESTS #########################################################
## These are standard tests of meta analysis, these test if our assumptions are met
# First, an analysis of variance (ANOVA):



MassAnova <- lmer(lineEarCountPercent ~ inoculation * seed + (1|row), NoControlDeadRemoved)
Anova(MassAnova, Type = "III")
summary(MassAnova)

qqnorm(residuals(MassAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(MassAnova)) # The closer the dots follow the line, the more likely the data is normal

shapiro.test(residuals(MassAnova)) # if the p-value is below 0.05, it's probably not normal
leveneTest(lineEarCountPercent ~ inoculation * seed, NoControlDeadRemoved) # if the p-value is below 0.05, we should probably transform it...
# BUT, if all other normality tests are greater than 0.05, this 
# becomes subjective... we don't have to. 

plot(fitted(MassAnova), residuals(MassAnova), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red") # With this plot we want the dots to seem random. If there is an obvious pattern or a cluster,
# that would indicate there is some factor influencing some of our data.
# Ideally, things will be somewhat clustered around the red "0" line.
# If anything is really 'out there' from everything else, it might be a statistical outlier 
# and maybe it should be discarded (this is subjective!)

## Because of non-normality, we will transform the data
massMin <- NoControlDeadRemoved$lineEarCountPercent+1

boxcox_result <- boxcox(massMin ~ 1, data = NoControlDeadRemoved, plotit = TRUE)
lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]

NoControlDeadRemoved$lineEarCountPercent_transformed <- if (lambda_optimal == 0) log(NoControlDeadRemoved$lineEarCountPercent) else (NoControlDeadRemoved$lineEarCountPercent^lambda_optimal - 1) / lambda_optimal

histMax <- max(NoControlDeadRemoved$lineEarCountPercent_transformed)+1.5
histBy <- max(NoControlDeadRemoved$lineEarCountPercent_transformed)/50
hist(NoControlDeadRemoved$lineEarCountPercent, 
     breaks = seq(-0.5, 500, by = 5.5), #Adjust these until you have the fidelity you want
     main = "Percent Masses", 
     xlab = "", 
     col = "darkred", 
     border = "white")

MassAnova <- lmer(lineEarCountPercent_transformed ~ inoculation * seed + (1|row), NoControlDeadRemoved)
Anova(MassAnova, Type = "III")
summary(Anova(MassAnova, Type = "III"))
qqnorm(residuals(MassAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(MassAnova))

######### NOW LENGTHS

LengthAnova <- lmer(lineStandCountPercent ~ inoculation * seed + (1|row), NoControlDeadRemoved)
summary(LengthAnova)
Anova(LengthAnova, Type = "III")

qqnorm(residuals(LengthAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(LengthAnova)) # The closer the dots follow the line, the more likely the data is normal

shapiro.test(residuals(LengthAnova)) # if the p-value is below 0.05, it's probably not normal
leveneTest(lineStandCountPercent ~ inoculation * seed, NoControlDeadRemoved) # if the p-value is below 0.05, we should probably transform it...
# BUT, if all other normality tests are greater than 0.05, this 
# becomes subjective... we don't have to. 

plot(fitted(LengthAnova), residuals(LengthAnova), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red") # With this plot we want the dots to seem random. If there is an obvious pattern or a cluster,
# that would indicate there is some factor influencing some of our data.
# Ideally, things will be somewhat clustered around the red "0" line.
# If anything is really 'out there' from everything else, it might be a statistical outlier 
# and maybe it should be discarded (this is subjective!)

## Because of non-normality, we will transform the data

LengthMin <- NoControlDeadRemoved$lineStandCountPercent + 1

boxcox_result <- boxcox(LengthMin ~ 1, data = NoControlDeadRemoved, plotit = TRUE)
lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]

NoControlDeadRemoved$lineStandCountPercent_transformed <- if (lambda_optimal == 0) log(NoControlDeadRemoved$lineStandCountPercent) else (NoControlDeadRemoved$lineStandCountPercent^lambda_optimal - 1) / lambda_optimal

histMax <- max(NoControlDeadRemoved$lineStandCountPercent_transformed)+1.5
histBy <- max(NoControlDeadRemoved$lineStandCountPercent_transformed)/50
hist(NoControlDeadRemoved$lineStandCountPercent, 
     breaks = seq(-0.5, 500, by = 5.5), #Adjust these until you have the fidelity you want
     main = "Percent Lengthes", 
     xlab = "", 
     col = "darkred", 
     border = "white")

LengthAnova <- lmer(lineStandCountPercent_transformed ~ inoculation * seed + (1|row), NoControlDeadRemoved)
summary(LengthAnova)
Anova(LengthAnova, Type = "III")

qqnorm(residuals(LengthAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(LengthAnova))

########## NOW DISEASE!


StuntCountAnova <- lmer(lineStuntCountPercent ~ inoculation * seed + (1|row), NoControlDeadRemoved)
summary(StuntCountAnova)
Anova(StuntCountAnova, Type = "III")

qqnorm(residuals(StuntCountAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(StuntCountAnova)) # The closer the dots follow the line, the more likely the data is normal

shapiro.test(residuals(StuntCountAnova)) # if the p-value is below 0.05, it's probably not normal
leveneTest(lineStuntCountPercent ~ inoculation * seed, NoControlDeadRemoved) # if the p-value is below 0.05, we should probably transform it...
# BUT, if all other normality tests are greater than 0.05, this 
# becomes subjective... we don't have to. 

plot(fitted(StuntCountAnova), residuals(StuntCountAnova), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red") # With this plot we want the dots to seem random. If there is an obvious pattern or a cluster,
# that would indicate there is some factor influencing some of our data.
# Ideally, things will be somewhat clustered around the red "0" line.
# If anything is really 'out there' from everything else, it might be a statistical outlier 
# and maybe it should be discarded (this is subjective!)

## Because of non-normality, we will transform the data

shifted_lesion_percent <- NoControlDeadRemoved$lineStuntCountPercent + 1

boxcox_result <- boxcox(shifted_lesion_percent ~ 1, data = NoControlDeadRemoved, plotit = TRUE)
lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]

NoControlDeadRemoved$lineStuntCountPercent_transformed <- if (lambda_optimal == 0) log(NoControlDeadRemoved$lineStuntCountPercent) else (NoControlDeadRemoved$lineStuntCountPercent^lambda_optimal - 1) / lambda_optimal

histMax <- max(NoControlDeadRemoved$lineStuntCountPercent_transformed)+1.5
histBy <- max(NoControlDeadRemoved$lineStuntCountPercent_transformed)/50
hist(NoControlDeadRemoved$lineStuntCountPercent, 
     breaks = seq(-0.5, 500, by = 5.5), #Adjust these until you have the fidelity you want
     main = "Percent Stunt Countes", 
     xlab = "", 
     col = "darkred", 
     border = "white")

StuntCountAnova <- lmer(lineStuntCountPercent_transformed ~ inoculation * seed + (1|row), NoControlDeadRemoved)
summary(StuntCountAnova)
Anova(StuntCountAnova, Type = "III")

shapiro.test(residuals(StuntCountAnova))
qqnorm(residuals(StuntCountAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(StuntCountAnova))


# groupVariances <- tapply(dataPuvu$averageOospores, dataPuvu$inoculation, var) # we will look at this when we have more data
# hist(groupVariances, xlab = "Variance", main = "Histogram of Variances", breaks = 20, col = fillTwo)

##### With our data checked we can do more analysis and move into graphs:

TukeyHSD(LengthAnova) # This will run an ANOVA for every variable that we specified earlier
# It is how we will determine which inoculations are significantly different
# from others. The more data we have, the better we can determine!

########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
################################# BUILDING GRAPHS WITH TRANSFORMED SIGNIFICANCE ############
########################################################################################################################################################################################
########################################################################################################################################################################################
fgDataTransformed <- NoControlDeadRemoved[grepl("Fusarium", NoControlDeadRemoved$inoculation),]
miDataTransformed <- NoControlDeadRemoved[grepl("Mock Inoc", NoControlDeadRemoved$inoculation),]
puuDataTransformed <- NoControlDeadRemoved[grepl("Pythium", NoControlDeadRemoved$inoculation),]
########################################################################################################################################################################################

summary(aov(lineStuntCountPercent_transformed ~ seed * inoculation + Error(seed), NoControlDeadRemoved))
summary(aov(lineStandCountPercent_transformed ~ seed * inoculation + Error(seed), NoControlDeadRemoved))
summary(aov(lineEarCountPercent_transformed ~ seed * inoculation + Error(seed), NoControlDeadRemoved))

model_lesion <- aov(lineStuntCountPercent_transformed ~ inoculation, NoControlDeadRemoved)
model_length <- aov(lineStandCountPercent_transformed ~ inoculation, NoControlDeadRemoved)
model_mass <- aov(lineEarCountPercent_transformed ~ inoculation, NoControlDeadRemoved)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$inoculation$Letters)
letters_lesionDf$inoculation <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
NoControlDeadRemoved <- merge(NoControlDeadRemoved, letters_lesionDf, by = "inoculation")
summaryData <- NoControlDeadRemoved %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(NoControlDeadRemoved) +
  aes(
    x = inoculation,
    y = lineStuntCountPercent,
    fill = inoculation
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculum Species",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = inoculation, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Treatment Stunt Count_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Stunt Count_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$inoculation$Letters)
letters_lengthDf$inoculation <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
NoControlDeadRemoved <- merge(NoControlDeadRemoved, letters_lengthDf, by = "inoculation")
summaryData <- NoControlDeadRemoved %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent_transformed = median(lineStandCountPercent_transformed), lengthLetters = first(lengthLetters))

ggplot(NoControlDeadRemoved) +
  aes(
    x = inoculation,
    y = lineStandCountPercent,
    fill = inoculation
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculum Species",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = inoculation, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Treatment Stand Count of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Stand Count of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$inoculation$Letters)
letters_massDf$inoculation <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
NoControlDeadRemoved <- merge(NoControlDeadRemoved, letters_massDf, by = "inoculation")
summaryData <- NoControlDeadRemoved %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent_transformed = median(lineEarCountPercent_transformed), massLetters = first(massLetters))

ggplot(NoControlDeadRemoved) +
  aes(
    x = inoculation,
    y = lineEarCountPercent,
    fill = inoculation
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculum Species",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = inoculation, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Treatment Ear Count of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Ear Count of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


######################### IS THERE A DIFFERENCE IN ROOT ROTS? ##################################################################################################

model_lesion <- aov(lineStuntCountPercent_transformed ~ seed, NoControlDeadRemoved)
model_length <- aov(lineStandCountPercent_transformed ~ seed, NoControlDeadRemoved)
model_mass <- aov(lineEarCountPercent_transformed ~ seed, NoControlDeadRemoved)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
NoControlDeadRemoved <- merge(NoControlDeadRemoved, letters_lesionDf, by = "seed")
summaryData <- NoControlDeadRemoved %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent), LesionLetters = first(LesionLetters))

ggplot(NoControlDeadRemoved) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculum Species",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pedigree Root Rot Stunt Count_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Root Rot Stunt Count_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
NoControlDeadRemoved <- merge(NoControlDeadRemoved, letters_lengthDf, by = "seed")
summaryData <- NoControlDeadRemoved %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent_transformed = median(lineStandCountPercent_transformed), lengthLetters = first(lengthLetters))

ggplot(NoControlDeadRemoved) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculum Species",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pedigree Stunt Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Stunt Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
NoControlDeadRemoved <- merge(NoControlDeadRemoved, letters_massDf, by = "seed")
summaryData <- NoControlDeadRemoved %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent_transformed = median(lineEarCountPercent_transformed), massLetters = first(massLetters))

ggplot(NoControlDeadRemoved) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculum Species",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pedigree Ear Count Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Ear Count Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


############ INDIVIDUAL SPECIES ANALYSIS #######################################################################################################################################


##Fusarium only##########################################################################################################################################################

model_lesion <- aov(lineStuntCountPercent_transformed ~ seed, fgDataTransformed)
model_length <- aov(lineStandCountPercent_transformed ~ seed, fgDataTransformed)
model_mass <- aov(lineEarCountPercent_transformed ~ seed, fgDataTransformed)

summary(aov(lineStuntCountPercent_transformed ~ seed + row, fgDataTransformed))
summary(aov(lineStandCountPercent_transformed ~ seed + row, fgDataTransformed))
summary(aov(lineEarCountPercent_transformed ~ seed + row, fgDataTransformed))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
fgDataTransformed <- merge(fgDataTransformed, letters_lesionDf, by = "seed")
summaryData <- fgDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(fgDataTransformed) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Fusarium Root Rot Stunt Count_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Fusarium Root Rot Stunt Count_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
fgDataTransformed <- merge(fgDataTransformed, letters_lengthDf, by = "seed")
summaryData <- fgDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent_transformed = median(lineStandCountPercent_transformed), lengthLetters = first(lengthLetters))

ggplot(fgDataTransformed) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Fusarium Stunt Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Fusarium Stunt Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
fgDataTransformed <- merge(fgDataTransformed, letters_massDf, by = "seed")
summaryData <- fgDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent_transformed = median(lineEarCountPercent_transformed), massLetters = first(massLetters))

ggplot(fgDataTransformed) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Fusarium Ear Count Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Fusarium Ear Count Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##Mock Inoc only##############################################################################################################################################################

model_lesion <- aov(lineStuntCountPercent_transformed ~ seed, miDataTransformed)
model_length <- aov(lineStandCountPercent_transformed ~ seed, miDataTransformed)
model_mass <- aov(lineEarCountPercent_transformed ~ seed, miDataTransformed)

summary(aov(lineStuntCountPercent_transformed ~ seed + row, miDataTransformed))
summary(aov(lineStandCountPercent_transformed ~ seed + row, miDataTransformed))
summary(aov(lineEarCountPercent_transformed ~ seed + row, miDataTransformed))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
miDataTransformed <- merge(miDataTransformed, letters_lesionDf, by = "seed")
summaryData <- miDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(miDataTransformed) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Mock Inoc Root Rot Stunt Count_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Mock Inoc Root Rot Stunt Count_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
miDataTransformed <- merge(miDataTransformed, letters_lengthDf, by = "seed")
summaryData <- miDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent_transformed = median(lineStandCountPercent_transformed), lengthLetters = first(lengthLetters))

ggplot(miDataTransformed) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Mock Inoc Stunt Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Mock Inoc Stunt Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
miDataTransformed <- merge(miDataTransformed, letters_massDf, by = "seed")
summaryData <- miDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent_transformed = median(lineEarCountPercent_transformed), massLetters = first(massLetters))

ggplot(miDataTransformed) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Mock Inoc Ear Count Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Mock Inoc Ear Count Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##Pythium only##############################################################################################################################################################

model_lesion <- aov(lineStuntCountPercent_transformed ~ seed, puuDataTransformed)
model_length <- aov(lineStandCountPercent_transformed ~ seed, puuDataTransformed)
model_mass <- aov(lineEarCountPercent_transformed ~ seed, puuDataTransformed)

summary(aov(lineStuntCountPercent_transformed ~ seed + row, puuDataTransformed))
summary(aov(lineStandCountPercent_transformed ~ seed + row, puuDataTransformed))
summary(aov(lineEarCountPercent_transformed ~ seed + row, puuDataTransformed))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$seed$Letters)
letters_lesionDf$seed <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
puuDataTransformed <- merge(puuDataTransformed, letters_lesionDf, by = "seed")
summaryData <- puuDataTransformed %>%
  group_by(seed) %>%
  dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(puuDataTransformed) +
  aes(
    x = seed,
    y = lineStuntCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stunt Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pythium Root Rot Stunt Count_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pythium Root Rot Stunt Count_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$seed$Letters)
letters_lengthDf$seed <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
puuDataTransformed <- merge(puuDataTransformed, letters_lengthDf, by = "seed")
summaryData <- puuDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineStandCountPercent_transformed = median(lineStandCountPercent_transformed), lengthLetters = first(lengthLetters))

ggplot(puuDataTransformed) +
  aes(
    x = seed,
    y = lineStandCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Stand Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pythium Stunt Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pythium Stunt Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$seed$Letters)
letters_massDf$seed <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
puuDataTransformed <- merge(puuDataTransformed, letters_massDf, by = "seed")
summaryData <- puuDataTransformed %>%
  group_by(inoculation, seed) %>%
  dplyr::summarize(medianlineEarCountPercent_transformed = median(lineEarCountPercent_transformed), massLetters = first(massLetters))

ggplot(puuDataTransformed) +
  aes(
    x = seed,
    y = lineEarCountPercent,
    fill = seed
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Ear Count as Percent of Standards" 
  ) +
  # Use the summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = seed, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.5
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )
filename <- paste0("Pythium Ear Count Percent of Standards_transformed.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pythium Ear Count Percent of Standards_transformed.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
model_lesion <- aov(lineStuntCountPercent_transformed ~ seed, puuDataTransformed)
model_length <- aov(lineStandCountPercent_transformed ~ seed, puuDataTransformed)
model_mass <- aov(lineEarCountPercent_transformed ~ seed, puuDataTransformed)

####################################################################################################################################################################################

seedList <- unique(NoControlDeadRemoved$seed)

# Loop through each seed and create a plot
for(seed in seedList) {
  
  # Filter data for the current seed
  dataSubsetTransformed <- NoControlDeadRemoved %>% filter(seed == !!seed)
  Sys.sleep(.1)
  print(paste("Processing seed:", seed))
  Sys.sleep(.1)
  
  summary(aov(lineStuntCountPercent_transformed ~ inoculation + row, dataSubsetTransformed))
  summary(aov(lineStandCountPercent_transformed ~ inoculation + row, dataSubsetTransformed))
  summary(aov(lineEarCountPercent_transformed ~ inoculation + row, dataSubsetTransformed))
  
  model_lesion <- aov(lineStuntCountPercent_transformed ~ inoculation, dataSubsetTransformed)
  model_length <- aov(lineStandCountPercent_transformed ~ inoculation, dataSubsetTransformed)
  model_mass <- aov(lineEarCountPercent_transformed ~ inoculation, dataSubsetTransformed)
  
  tukey_lesion <- TukeyHSD(model_lesion)
  tukey_length <- TukeyHSD(model_length)
  tukey_mass <- TukeyHSD(model_mass)
  
  letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
  letters_length <- multcompLetters4(model_length, tukey_length)
  letters_mass <- multcompLetters4(model_mass, tukey_mass)
  
  
  labelName <- paste0("Disease % Stunt Count by line:", seed, sep = " ")
  ylab = labelName
  
  letters_lesionDf <- as.data.frame(letters_lesion$inoculation$Letters)
  
  letters_lesionDf$inoculation <- rownames(letters_lesionDf)
  
  colnames(letters_lesionDf)[1] <- "LesionLetters_ind"
  
  rownames(letters_lesionDf) <- NULL
  
  dataSubsetTransformed <- merge(dataSubsetTransformed, letters_lesionDf, by = "inoculation")
  
  letters_lesionDf <- as.data.frame(letters_lesion$inoculation$Letters)
  
  letters_lesionDf$inoculation <- rownames(letters_lesionDf)
  
  colnames(letters_lesionDf)[1] <- "LesionLetters_ind"
  
  rownames(letters_lesionDf) <- NULL
  
  dataSubsetTransformed <- merge(dataSubsetTransformed, letters_lesionDf, by = "inoculation")
  
  summaryData <- dataSubsetTransformed %>%
    group_by(inoculation, seed) %>%
    dplyr::summarize(medianlineStuntCountPercent = median(lineStuntCountPercent_transformed), LesionLetters = first(LesionLetters), .groups = "drop")
  
  
  ggplot(dataSubsetTransformed) +
    aes(
      x = inoculation,
      y = lineStuntCountPercent,
      fill = inoculation
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey() +
    geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
    labs(
      x = "Inoculum Species",  
      y = "Stunt Count as Percent of Standards" 
    ) +
    # Use the summarized data to plot the letters
    geom_text(
      data = summaryData, 
      aes(
        x = inoculation, 
        y = -Inf,  # Adjust to place the letters above the boxplot
        label = LesionLetters
      ), 
      position = position_dodge(width = 0.75), 
      vjust = -0.5
    ) +
    labs(y = ylab) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),
    )
  
  filename <- paste0(seed, "_percent Stunt Count_transformed.jpg") 
  ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  filename <- paste0(seed, "_percent Stunt Count_transformed.wmf")  
  ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  Sys.sleep(.5)
  
  
  labelName <- paste0("% StandCount by Line:", seed, sep = " ")
  ylab = labelName
  
  letters_lengthDf <- as.data.frame(letters_length$inoculation$Letters)
  
  letters_lengthDf$inoculation <- rownames(letters_lengthDf)
  
  colnames(letters_lengthDf)[1] <- "lengthLetters_ind"
  
  rownames(letters_lengthDf) <- NULL
  
  dataSubsetTransformed <- merge(dataSubsetTransformed, letters_lengthDf, by = "inoculation")
  
  letters_lengthDf <- as.data.frame(letters_length$inoculation$Letters)
  
  letters_lengthDf$inoculation <- rownames(letters_lengthDf)
  
  colnames(letters_lengthDf)[1] <- "lengthLetters_ind"
  
  rownames(letters_lengthDf) <- NULL
  
  dataSubsetTransformed <- merge(dataSubsetTransformed, letters_lengthDf, by = "inoculation")
  
  summaryData <- dataSubsetTransformed %>%
    group_by(inoculation, seed) %>%
    dplyr::summarize(medianlineStandCountPercent_transformed = median(lineStandCountPercent_transformed), lengthLetters = first(lengthLetters), .groups = "drop")
  
  
  ggplot(dataSubsetTransformed) +
    aes(
      x = inoculation,
      y = lineStandCountPercent,
      fill = inoculation
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey() +
    geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
    labs(
      x = "Inoculum Species",  
      y = "Stand Count as Percent of Standards" 
    ) +
    # Use the summarized data to plot the letters
    geom_text(
      data = summaryData, 
      aes(
        x = inoculation, 
        y = -Inf,  # Adjust to place the letters above the boxplot
        label = lengthLetters
      ), 
      position = position_dodge(width = 0.75), 
      vjust = -0.5
    ) +
    labs(y = ylab) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),
    )
  
  filename <- paste0(seed, "_percent Stand Count_transformed.jpg") 
  ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  filename <- paste0(seed, "_percent Stand Count_transformed.wmf")  
  ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  Sys.sleep(.5)
  
  labelName <- paste0("% EarCount by Line:", seed, sep = " ")
  ylab = labelName
  
  letters_massDf <- as.data.frame(letters_mass$inoculation$Letters)
  
  letters_massDf$inoculation <- rownames(letters_massDf)
  
  colnames(letters_massDf)[1] <- "massLetters_ind"
  
  rownames(letters_massDf) <- NULL
  
  dataSubsetTransformed <- merge(dataSubsetTransformed, letters_massDf, by = "inoculation")
  
  letters_massDf <- as.data.frame(letters_mass$inoculation$Letters)
  
  letters_massDf$inoculation <- rownames(letters_massDf)
  
  colnames(letters_massDf)[1] <- "massLetters_ind"
  
  rownames(letters_massDf) <- NULL
  
  dataSubsetTransformed <- merge(dataSubsetTransformed, letters_massDf, by = "inoculation")
  
  summaryData <- dataSubsetTransformed %>%
    group_by(inoculation, seed) %>%
    dplyr::summarize(medianlineEarCountPercent_transformed = median(lineEarCountPercent_transformed), massLetters = first(massLetters), .groups = "drop")
  
  
  ggplot(dataSubsetTransformed) +
    aes(
      x = inoculation,
      y = lineEarCountPercent,
      fill = inoculation
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_grey() +
    geom_point(shape = 25, color = "black", size = 1.5, alpha = 1, position = position_jitter(width = 0.2)) +
    labs(
      x = "Inoculum Species",  
      y = "Ear Count as Percent of Standards" 
    ) +
    # Use the summarized data to plot the letters
    geom_text(
      data = summaryData, 
      aes(
        x = inoculation, 
        y = -Inf,  # Adjust to place the letters above the boxplot
        label = massLetters
      ), 
      position = position_dodge(width = 0.75), 
      vjust = -0.5
    ) +
    labs(y = ylab) +
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),
    )
  
  filename <- paste0(seed, "_percent Ear Count_transformed.jpg") 
  ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  filename <- paste0(seed, "_percent Ear Count_transformed.wmf")  
  ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
  Sys.sleep(1)
  
}

########################################################################################################################################################################################
######################################################################## CORRELATIONS ##################################################################################################
########################################################################################################################################################################################

############################################ MATRIX OF CORRELATIONS VIA PEARSON CORRELATION ###############################################
# Define the variable names

library("Hmisc")
library("reshape2")

vars <- c("Stunt Count","Percent Stand Count","Percent Ear Count")

matrixCor <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
matrixP <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))

rownames(matrixCor) <- vars
colnames(matrixCor) <- vars
# Defaults
diag(matrixCor) <- 1
diag(matrixP) <- NA
#

################################################################################################################
model <- cor.test(exp_data$lineStuntCountPercent , exp_data$lineStandCountPercent , method = "pearson")

matrixCor["Stunt Count", "Percent Stand Count"] <- model$estimate
matrixCor["Percent Stand Count", "Stunt Count"] <- model$estimate
matrixP["Stunt Count", "Percent Stand Count"] <- model$p.value
matrixP["Percent Stand Count", "Stunt Count"] <- model$p.value
#

model <- cor.test(exp_data$lineEarCountPercent , exp_data$lineStandCountPercent , method = "pearson")

matrixCor["Percent Ear Count", "Percent Stand Count"] <- model$estimate 
matrixCor["Percent Stand Count", "Percent Ear Count"] <- model$estimate 
matrixP["Percent Ear Count", "Percent Stand Count"] <- model$p.value
matrixP["Percent Stand Count", "Percent Ear Count"] <- model$p.value
#

model <- cor.test(exp_data$lineEarCountPercent , exp_data$lineStuntCountPercent , method = "pearson")

matrixCor["Percent Ear Count", "Stunt Count"] <- model$estimate 
matrixCor["Stunt Count", "Percent Ear Count"] <- model$estimate 
matrixP["Percent Ear Count", "Stunt Count"] <- model$p.value
matrixP["Stunt Count", "Percent Ear Count"] <- model$p.value
#

################################################ MATRIX GRAPH ##################################################

# Naming rows and columns
rownames(matrixCor) <- vars
rownames(matrixP) <- vars
colnames(matrixCor) <- vars
colnames(matrixP) <- vars

# Melting the data frames
meltedCor <- melt(matrixCor, varnames = c("Var1", "Var2"), value.name = "Correlation")
meltedP <- melt(matrixP, varnames = c("Var1", "Var2"), value.name = "PValue")

# Combining melted data frames
combinedMelted <- cbind(meltedCor, PValue = meltedP$PValue)

# Ensure Var1 and Var2 are factors or characters, with Var1 sorted if necessary
combinedMelted$Var1 <- factor(combinedMelted$Var1, levels = unique(c(combinedMelted$Var1, combinedMelted$Var2)))
combinedMelted$Var2 <- factor(combinedMelted$Var2, levels = levels(combinedMelted$Var1))

# Filter out the upper triangle, including diagonal
combinedMeltedLower <- combinedMelted[as.numeric(combinedMelted$Var1) >= as.numeric(combinedMelted$Var2), ]

ggplot(combinedMeltedLower, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("r = %.2f\np = %.3f", Correlation, PValue)), size = 5, vjust = 0.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, vjust = 1)) +
  labs(title = "Correlation and P-Values")

filename <- paste0("Correlation and P-Values.jpg", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Correlation and P-Values.wmf", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

