################################################################################
## MultiSpecies Script, Summer 2024  
## Analysis by Harrison Hall (harrisonpiercehall@gmail.com) & Kim Hagemann 
################################################################################

## Libraries used 
library("dplyr") #contains various helpful functions
library("tidyr")
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
library("extrafont")

fonts()
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
fontSize <- 7 # Can be manually changed later
# Set the default exported chart size to 5000x2000 units, this is print quality
ChartWidth <- (5760)
ChartHeight <- (3240)
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
theme_set(theme_gray(base_family = "sans"))
################################################################################
##### Beginning of analysis ####################################################
################################################################################

## First let's look at our data:

head(data) # Take a peek at the data and make sure we grabbed the right file
str(data) # This shows us our column names and examples of data. We will want to make sure things are formatted the way we want.


#esquisser(data)
## If it is a category or qualitative (red/blue) trait, we want it to be a factor
## If it is a quantitative trait (length or mass), we want it to be numeric

data$Seed_ID <- as.factor(data$Seed_ID)
data$Packet_ID <- as.factor(data$Packet_ID)
data$pedigree <- as.factor(data$pedigree) 
data$DPI <- as.factor(data$DPI)
data$repetition <- as.factor(data$repetition)
data$isolate <- as.factor(data$isolate)
data$inocMethod <- as.factor(data$inocMethod)
data$isolateSpecies <- as.factor(data$isolateSpecies)
data$agarType <- as.factor(data$agarType)
data$soilType <- as.factor(data$soilType)
data$seedNum <- as.factor(data$seedNum)
data$block <- as.factor(data$block)
data$germinated <- as.numeric(data$germinated)
data$escape <- as.factor(data$escape)
#Numerics:
data$lesionPercent <- as.numeric(data$lesionPercent)
data$rootLength <- as.numeric(data$rootLength)
data$rootMass <- as.numeric(data$rootMass)
data$averageInocLength <- as.numeric(data$averageInocLength)
data$averageInocMass <- as.numeric(data$averageInocMass)
data$percentLengthPedByIsolate <- as.numeric(data$percentLengthPedByIsolate)*100
data$percentMassPedByIsolate <- as.numeric(data$percentMassPedByIsolate)*100
#Character
data$notes <- as.character(data$notes)

# Now remove escapes:

data <- data[!grepl("1", data$escape),]
#data <- data[!grepl("CML277", data$pedigree),] # removing this since almost all of them died
ogData <- data
germData <- ogData
## Let's open the data in esquisse and see what kind of graphs we can make, and visualize some of the data.
#esquisser(data)

germinated_cols <- grep("germinated", names(germData), value = TRUE)
germData$totalGerminatedRow <- rowSums(germData[, germinated_cols], na.rm = TRUE)

germResult <- germData %>%
  group_by(pedigree, isolateSpecies) %>%
  summarise(totalGerminated = sum(totalGerminatedRow, na.rm = TRUE), .groups = 'drop')

filename <- paste0(tablesDir,"/germination rate table.csv", sep = "") 
write.csv(germResult, filename, row.names = FALSE)

germWide <- germResult %>%
  pivot_wider(names_from = isolateSpecies, 
              values_from = totalGerminated, 
              names_prefix = "",
              names_sep = "_totalGerminated") %>%
  rename_with(~ paste0(., "_totalGerminated"), -pedigree)

filename <- paste0(tablesDir,"/germination rate table wide.csv", sep = "") 
write.csv(germWide, filename, row.names = FALSE)



################################################################################################################################################################
# Experiment Analysis - Did the experimental design work?


stats <- ogData %>%
  group_by(pedigree) %>%
  summarise(
    meanLesionPercent = mean(lesionPercent, na.rm = TRUE),
    sdLesionPercent = sd(lesionPercent, na.rm = TRUE),
    meanPercentLength = mean(percentLengthPedByIsolate, na.rm = TRUE),
    sdPercentLength = sd(percentLengthPedByIsolate, na.rm = TRUE),
    meanPercentMass = mean(percentMassPedByIsolate, na.rm = TRUE),
    sdPercentMass = sd(percentMassPedByIsolate, na.rm = TRUE)
  )

# Join this stats back to the original data
filteredData <- ogData %>%
  left_join(stats, by = "pedigree") %>%
  filter(
    lesionPercent <= meanLesionPercent + 2.1 * sdLesionPercent,
    percentLengthPedByIsolate <= meanPercentLength + 2.1 * sdPercentLength,
    percentMassPedByIsolate <= meanPercentMass + 2.1 * sdPercentMass
  )

# Select columns to avoid clutter
exp_data <- filteredData %>%
  dplyr::select(-starts_with("mean"), -starts_with("sd"))

# View the filtered dataset
print(exp_data)

######################################################   FIRST ANOVAs    ######################################
AnovaData <- exp_data[!grepl("Standard", exp_data$isolate),]

summary(aov(lesionPercent ~ pedigree * isolateSpecies + repetition + block, AnovaData))
summary(aov(percentLengthPedByIsolate ~ pedigree * isolateSpecies + repetition + block, AnovaData))
summary(aov(percentMassPedByIsolate ~ pedigree * isolateSpecies + repetition + block, AnovaData))

anovaFile <- paste0(anovaDir,"/anova_summaries.txt",sep="")

# Use capture.output to save each summary to the file
capture.output({
  cat("Summary of ANOVA for lesionPercent:\n")
  print(summary(aov(lesionPercent ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), AnovaData)))
  
  cat("\n\nSummary of ANOVA for percentLengthPedByIsolate:\n")
  print(summary(aov(percentLengthPedByIsolate ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), AnovaData)))
  
  cat("\n\nSummary of ANOVA for percentMassPedByIsolate:\n")
  print(summary(aov(percentMassPedByIsolate ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), AnovaData)))
}, file = anovaFile)


model_lesion <- aov(lesionPercent ~ isolateSpecies, exp_data)
model_length <- aov(percentLengthPedByIsolate ~ isolateSpecies, exp_data)
model_mass <- aov(percentMassPedByIsolate ~ isolateSpecies, exp_data)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

# Calculate means for each species
mean_lesions <- exp_data %>%
  dplyr::group_by(isolateSpecies) %>%
  dplyr::summarize(meanLesion = mean(lesionPercent, na.rm = TRUE)) %>%
  dplyr::arrange(meanLesion)

# Move 'Standard' to the first position
standard_first <- mean_lesions %>%
  dplyr::filter(isolateSpecies == "Standard") %>%
  dplyr::bind_rows(mean_lesions %>%
              filter(isolateSpecies != "Standard"))

# Reorder the factor levels in the original data
exp_data$isolateSpecies <- factor(exp_data$isolateSpecies, levels = standard_first$isolateSpecies)



letters_lesionDf <- as.data.frame(letters_lesion$isolateSpecies$Letters)
letters_lesionDf$isolateSpecies <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
exp_data <- merge(exp_data, letters_lesionDf, by = "isolateSpecies")
summaryData <- exp_data %>%
  group_by(inocMethod, isolateSpecies) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent), LesionLetters = first(LesionLetters))

ggplot(exp_data) +
  aes(
    x = isolateSpecies,
    y = lesionPercent,
    fill = isolateSpecies
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey()+
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Disease Severity (Percent)" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = isolateSpecies, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1,
    fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "sans", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),
  )

filename <- paste0("Treatment Disease Severity (Percent).pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Disease Severity (Percent).png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)



# Calculate means for each species
mean_lengths <- exp_data %>%
  dplyr::group_by(isolateSpecies) %>%
  dplyr::summarize(meanLesion = mean(percentLengthPedByIsolate, na.rm = TRUE)) %>%
  dplyr::arrange(meanLesion)

# Move 'Standard' to the first position
standard_first <- mean_lengths %>%
  dplyr::filter(isolateSpecies == "Standard") %>%
  dplyr::bind_rows(mean_lengths %>%
                     filter(isolateSpecies != "Standard"))

# Reorder the factor levels in the original data
exp_data$isolateSpecies <- factor(exp_data$isolateSpecies, levels = standard_first$isolateSpecies)

letters_lengthDf <- as.data.frame(letters_length$isolateSpecies$Letters)
letters_lengthDf$isolateSpecies <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
exp_data <- merge(exp_data, letters_lengthDf, by = "isolateSpecies")
summaryData <- exp_data %>%
  group_by(inocMethod, isolateSpecies) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate = median(percentLengthPedByIsolate), lengthLetters = first(lengthLetters))

ggplot(exp_data) +
  aes(
    x = isolateSpecies,
    y = percentLengthPedByIsolate,
    fill = isolateSpecies
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Root Length as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = isolateSpecies, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( 
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )

filename <- paste0("Treatment Length Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Length Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

# Calculate means for each species
mean_mass <- exp_data %>%
  dplyr::group_by(isolateSpecies) %>%
  dplyr::summarize(meanLesion = mean(percentMassPedByIsolate, na.rm = TRUE)) %>%
  dplyr::arrange(meanLesion)

# Move 'Standard' to the first position
standard_first <- mean_mass %>%
  dplyr::filter(isolateSpecies == "Standard") %>%
  dplyr::bind_rows(mean_mass %>%
                     filter(isolateSpecies != "Standard"))

# Reorder the factor levels in the original data
exp_data$isolateSpecies <- factor(exp_data$isolateSpecies, levels = standard_first$isolateSpecies)


letters_massDf <- as.data.frame(letters_mass$isolateSpecies$Letters)
letters_massDf$isolateSpecies <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
exp_data <- merge(exp_data, letters_massDf, by = "isolateSpecies")
summaryData <- exp_data %>%
  group_by(inocMethod, isolateSpecies) %>%
  dplyr::summarize(medianpercentMassPedByIsolate = median(percentMassPedByIsolate), massLetters = first(massLetters))

ggplot(exp_data) +
  aes(
    x = isolateSpecies,
    y = percentMassPedByIsolate,
    fill = isolateSpecies
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = .5, alpha = 1, position = position_jitter(width = 0.22)) +
  labs(
    x = "Inoculation Species",  
    y = "Root Mass as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = isolateSpecies, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -.10
  ) +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Treatment Mass Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Mass Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

if("LesionLetters" %in% colnames(exp_data)) {
  exp_data <- exp_data %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(exp_data)) {
  exp_data <- exp_data %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(exp_data)) {
  exp_data <- exp_data %>%
    dplyr::select(-massLetters)
}


######################### IS THERE A DIFFERENCE IN ROOT ROTS? ##################################################################################################

#NoControlData <- data[!grepl("Standard|G\\.a", data$isolate),]
NoControlData <- exp_data[!grepl("Standard", exp_data$isolate),]

model_lesion <- aov(lesionPercent ~ pedigree+repetition, NoControlData)
model_length <- aov(percentLengthPedByIsolate ~ pedigree+repetition, NoControlData)
model_mass <- aov(percentMassPedByIsolate ~ pedigree+repetition, NoControlData)

summary(aov(lesionPercent ~ pedigree*isolate+repetition+block, NoControlData))
summary(aov(percentLengthPedByIsolate ~ pedigree*isolate+repetition+block, NoControlData))
summary(aov(percentMassPedByIsolate ~ pedigree*isolate+repetition+block, NoControlData))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

# Calculate means for each maize line
mean_lesions <- NoControlData %>%
  dplyr::group_by(pedigree) %>%
  dplyr::summarize(meanLesion = mean(lesionPercent, na.rm = TRUE)) %>%
  dplyr::arrange(meanLesion)

# Reorder the factor levels in the original data
NoControlData$pedigree <- factor(NoControlData$pedigree, levels = mean_lesions$pedigree)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
NoControlData <- merge(NoControlData, letters_lesionDf, by = "pedigree")
summaryData <- NoControlData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent), LesionLetters = first(LesionLetters))

ggplot(NoControlData) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.22)) +
  labs(
    x = "Maize Line",  
    y = "Disease Severity (Percent)" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Pedigree Root Rot Disease Severity (Percent).pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Root Rot Disease Severity (Percent).png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


# Calculate means for each maize line
mean_lengths <- NoControlData %>%
  dplyr::group_by(pedigree) %>%
  dplyr::summarize(meanLesion = mean(percentLengthPedByIsolate, na.rm = TRUE)) %>%
  dplyr::arrange(desc(meanLesion))

# Reorder the factor levels in the original data
NoControlData$pedigree <- factor(NoControlData$pedigree, levels = mean_lengths$pedigree)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
NoControlData <- merge(NoControlData, letters_lengthDf, by = "pedigree")
summaryData <- NoControlData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate = median(percentLengthPedByIsolate), lengthLetters = first(lengthLetters))

ggplot(NoControlData) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Root Length as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Pedigree Root Rot Length Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Root Rot Length Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

# Calculate means for each maize line
mean_mass <- NoControlData %>%
  dplyr::group_by(pedigree) %>%
  dplyr::summarize(meanLesion = mean(percentMassPedByIsolate, na.rm = TRUE)) %>%
  dplyr::arrange(desc(meanLesion))

# Reorder the factor levels in the original data
NoControlData$pedigree <- factor(NoControlData$pedigree, levels = mean_mass$pedigree)


letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
NoControlData <- merge(NoControlData, letters_massDf, by = "pedigree")
summaryData <- NoControlData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate = median(percentMassPedByIsolate), massLetters = first(massLetters))

ggplot(NoControlData) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Line",  
    y = "Root Mass as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Pedigree Root Rot Mass Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree Root Rot Mass Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

if("LesionLetters" %in% colnames(NoControlData)) {
  NoControlData <- NoControlData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(NoControlData)) {
  NoControlData <- NoControlData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(NoControlData)) {
  NoControlData <- NoControlData %>%
    dplyr::select(-massLetters)
}


############ INDIVIDUAL SPECIES ANALYSIS #######################################################################################################################################
paData <- NoControlData[grepl("G.a", NoControlData$isolateSpecies),]
piData <- NoControlData[grepl("G.i", NoControlData$isolateSpecies),]
puuData <- NoControlData[grepl("G.uu", NoControlData$isolateSpecies),]

##aphanidermatum only##########################################################################################################################################################

if("LesionLetters" %in% colnames(paData)) {
  paData <- paData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(paData)) {
  paData <- paData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(paData)) {
  paData <- paData %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(piData)) {
  piData <- piData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(piData)) {
  piData <- piData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(piData)) {
  piData <- piData %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(puuData)) {
  puuData <- puuData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(puuData)) {
  puuData <- puuData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(puuData)) {
  puuData <- puuData %>%
    dplyr::select(-massLetters)
}


model_lesion <- aov(lesionPercent ~ pedigree+repetition, paData)
model_length <- aov(percentLengthPedByIsolate ~ pedigree+repetition, paData)
model_mass <- aov(percentMassPedByIsolate ~ pedigree+repetition, paData)

summary(aov(lesionPercent ~ pedigree+repetition + block, paData))
summary(aov(percentLengthPedByIsolate ~ pedigree+repetition + block, paData))
summary(aov(percentMassPedByIsolate ~ pedigree+repetition + block, paData))

anovaFile <- paste0(anovaDir,"/P.aphanidermatum_anova_summaries.txt",sep="")

# Use capture.output to save each summary to the file
capture.output({
  cat("Summary of ANOVA for lesionPercent:\n")
  print(summary(aov(lesionPercent ~ pedigree + block, paData)))
  
  cat("\n\nSummary of ANOVA for percentLengthPedByIsolate:\n")
  print(summary(aov(percentLengthPedByIsolate ~ pedigree + block, paData)))
  
  cat("\n\nSummary of ANOVA for percentMassPedByIsolate:\n")
  print(summary(aov(percentMassPedByIsolate ~ pedigree + block, paData)))
}, file = anovaFile)



tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
paData <- merge(paData, letters_lesionDf, by = "pedigree")
summaryData <- paData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent), LesionLetters = first(LesionLetters))

ggplot(paData) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. aphanidermatum",  
    y = "Disease Severity (Percent)" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Aphanidermatum Root Rot Disease Severity (Percent).pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Aphanidermatum Root Rot Disease Severity (Percent).png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
paData <- merge(paData, letters_lengthDf, by = "pedigree")
summaryData <- paData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate = median(percentLengthPedByIsolate), lengthLetters = first(lengthLetters))

ggplot(paData) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. aphanidermatum",  
    y = "Root Length as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Aphanidermatum Root Rot Length Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Aphanidermatum Root Rot Length Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
paData <- merge(paData, letters_massDf, by = "pedigree")
summaryData <- paData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate = median(percentMassPedByIsolate), massLetters = first(massLetters))

ggplot(paData) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. aphanidermatum",  
    y = "Root Mass as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Aphanidermatum Root Rot Mass Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Aphanidermatum Root Rot Mass Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##irregulare only##############################################################################################################################################################

model_lesion <- aov(lesionPercent ~ pedigree+repetition, piData)
model_length <- aov(percentLengthPedByIsolate ~ pedigree+repetition, piData)
model_mass <- aov(percentMassPedByIsolate ~ pedigree+repetition, piData)

summary(aov(lesionPercent ~ pedigree+repetition + block, piData))
summary(aov(percentLengthPedByIsolate ~ pedigree+repetition + block, piData))
summary(aov(percentMassPedByIsolate ~ pedigree+repetition + block, piData))


anovaFile <- paste0(anovaDir,"/P.irregulare_anova_summaries.txt",sep="")

# Use capture.output to save each summary to the file
capture.output({
  cat("Summary of ANOVA for lesionPercent:\n")
  print(summary(aov(lesionPercent ~ pedigree+repetition + block, piData)))
  
  cat("\n\nSummary of ANOVA for percentLengthPedByIsolate:\n")
  print(summary(aov(percentLengthPedByIsolate ~ pedigree+repetition + block, piData)))
  
  cat("\n\nSummary of ANOVA for percentMassPedByIsolate:\n")
  print(summary(aov(percentMassPedByIsolate ~ pedigree+repetition + block, piData)))
}, file = anovaFile)


tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
piData <- merge(piData, letters_lesionDf, by = "pedigree")
summaryData <- piData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent), LesionLetters = first(LesionLetters))

ggplot(piData) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. irregulare",  
    y = "Disease Severity (Percent)" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Irregulare Root Rot Disease Severity (Percent).pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Irregulare Root Rot Disease Severity (Percent).png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
piData <- merge(piData, letters_lengthDf, by = "pedigree")
summaryData <- piData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate = median(percentLengthPedByIsolate), lengthLetters = first(lengthLetters))

ggplot(piData) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. irregulare",  
    y = "Root Length as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Irregulare Root Rot Length Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Irregulare Root Rot Length Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
piData <- merge(piData, letters_massDf, by = "pedigree")
summaryData <- piData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate = median(percentMassPedByIsolate), massLetters = first(massLetters))

ggplot(piData) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. irregulare",  
    y = "Root Mass as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Irregulare Root Rot Mass Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Irregulare Root Rot Mass Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##ultimum ultimum only##############################################################################################################################################################

model_lesion <- aov(lesionPercent ~ pedigree+repetition, puuData)
model_length <- aov(percentLengthPedByIsolate ~ pedigree+repetition, puuData)
model_mass <- aov(percentMassPedByIsolate ~ pedigree+repetition, puuData)

summary(aov(lesionPercent ~ pedigree+repetition + block, puuData))
summary(aov(percentLengthPedByIsolate ~ pedigree+repetition + block, puuData))
summary(aov(percentMassPedByIsolate ~ pedigree+repetition + block, puuData))

anovaFile <- paste0(anovaDir,"/P.ultimum_anova_summaries.txt",sep="")

# Use capture.output to save each summary to the file
capture.output({
  cat("Summary of ANOVA for lesionPercent:\n")
  print(summary(aov(lesionPercent ~ pedigree+repetition + block, puuData)))
  
  cat("\n\nSummary of ANOVA for percentLengthPedByIsolate:\n")
  print(summary(aov(percentLengthPedByIsolate ~ pedigree+repetition + block, puuData)))
  
  cat("\n\nSummary of ANOVA for percentMassPedByIsolate:\n")
  print(summary(aov(percentMassPedByIsolate ~ pedigree+repetition + block, puuData)))
}, file = anovaFile)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
puuData <- merge(puuData, letters_lesionDf, by = "pedigree")
summaryData <- puuData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent), LesionLetters = first(LesionLetters))

ggplot(puuData) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. ultimum var. ultimum",  
    y = "Disease Severity (Percent)" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Ultimum Root Rot Disease Severity (Percent).pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Ultimum Root Rot Disease Severity (Percent).png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
puuData <- merge(puuData, letters_lengthDf, by = "pedigree")
summaryData <- puuData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate = median(percentLengthPedByIsolate), lengthLetters = first(lengthLetters))

ggplot(puuData) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. ultimum var. ultimum",  
    y = "Root Length as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Ultimum Root Rot Length Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Ultimum Root Rot Length Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
puuData <- merge(puuData, letters_massDf, by = "pedigree")
summaryData <- puuData %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate = median(percentMassPedByIsolate), massLetters = first(massLetters))

ggplot(puuData) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. ultimum var. ultimum",  
    y = "Root Mass as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Ultimum Root Rot Mass Percent of Standards.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Ultimum Root Rot Mass Percent of Standards.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

model_lesion <- aov(lesionPercent ~ pedigree+repetition, puuData)
model_length <- aov(percentLengthPedByIsolate ~ pedigree+repetition, puuData)
model_mass <- aov(percentMassPedByIsolate ~ pedigree+repetition, puuData)

if("LesionLetters" %in% colnames(paData)) {
  paData <- paData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(paData)) {
  paData <- paData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(paData)) {
  paData <- paData %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(piData)) {
  piData <- piData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(piData)) {
  piData <- piData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(piData)) {
  piData <- piData %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(puuData)) {
  puuData <- puuData %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(puuData)) {
  puuData <- puuData %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(puuData)) {
  puuData <- puuData %>%
    dplyr::select(-massLetters)
}

####################################################################################################################################################################################

pedigreeList <- unique(NoControlData$pedigree)

#dataForTransformation <- ogData[!grepl("Standards", ogData$isolate),]
#dataForTransformation <- ogData[!grepl("Standard", ogData$isolate),]
dataForTransformation <- ogData

histMax <- max(dataForTransformation$percentLengthPedByIsolate)+10.5
histBy <- max(dataForTransformation$percentLengthPedByIsolate)/25
hist(dataForTransformation$percentLengthPedByIsolate, 
     breaks = seq(-0.1, histMax, by = histBy), #Adjust these until you have the fidelity you want
     main = "Percent Lengths", 
     xlab = "", 
     col = "darkblue", 
     border = "white")
  # This can be an early indication of normality!

histMax <- max(dataForTransformation$percentMassPedByIsolate)+10.5
histBy <- max(dataForTransformation$percentMassPedByIsolate)/25
hist(dataForTransformation$percentMassPedByIsolate, 
     breaks = seq(-0.1, histMax, by = histBy), #Adjust these until you have the fidelity you want
     main = "Percent Masses", 
     xlab = "", 
     col = "darkred", 
     border = "white")
# This can be an early indication of normality!

##### ASSUMPTION TESTS #########################################################
## These are standard tests of meta analysis, these test if our assumptions are met
# First, an analysis of variance (ANOVA):



MassAnova <- lmer(percentMassPedByIsolate ~ isolate * pedigree + (1|block), dataForTransformation)
Anova(MassAnova, Type = "III")

summary(aov(percentLengthPedByIsolate ~ pedigree + block, dataForTransformation))

summary(MassAnova)

qqnorm(residuals(MassAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(MassAnova)) # The closer the dots follow the line, the more likely the data is normal

shapiro.test(residuals(MassAnova)) # if the p-value is below 0.05, it's probably not normal
leveneTest(percentMassPedByIsolate ~ isolate * pedigree, dataForTransformation) # if the p-value is below 0.05, we should probably transform it...
                                                # BUT, if all other normality tests are greater than 0.05, this 
                                                # becomes subjective... we don't have to. 

plot(fitted(MassAnova), residuals(MassAnova), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red") # With this plot we want the dots to seem random. If there is an obvious pattern or a cluster,
                           # that would indicate there is some factor influencing some of our data.
                           # Ideally, things will be somewhat clustered around the red "0" line.
                           # If anything is really 'out there' from everything else, it might be a statistical outlier 
                           # and maybe it should be discarded (this is subjective!)

## Because of non-normality, we will transform the data
massMin <- dataForTransformation$percentMassPedByIsolate+1

boxcox_result <- boxcox(massMin ~ 1, data = dataForTransformation, plotit = TRUE)
lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]

dataForTransformation$percentMassPedByIsolate_transformed <- if (lambda_optimal == 0) log(dataForTransformation$percentMassPedByIsolate) else (dataForTransformation$percentMassPedByIsolate^lambda_optimal - 1) / lambda_optimal

histMax <- max(dataForTransformation$percentMassPedByIsolate_transformed)+1.5
histBy <- max(dataForTransformation$percentMassPedByIsolate_transformed)/50
hist(dataForTransformation$percentMassPedByIsolate, 
     breaks = seq(-0.5, 500, by = 5.5), #Adjust these until you have the fidelity you want
     main = "Percent Masses", 
     xlab = "", 
     col = "darkred", 
     border = "white")

MassAnova <- lmer(percentMassPedByIsolate_transformed ~ isolate * pedigree + (1|block), dataForTransformation)
Anova(MassAnova, Type = "III")
summary(Anova(MassAnova, Type = "III"))

summary(aov(percentMassPedByIsolate_transformed ~ isolate * pedigree + block, dataForTransformation))


qqnorm(residuals(MassAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(MassAnova))

######### NOW LENGTHS

LengthAnova <- lmer(percentLengthPedByIsolate ~ isolate * pedigree + (1|block), dataForTransformation)
summary(LengthAnova)
Anova(LengthAnova, Type = "III")

qqnorm(residuals(LengthAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(LengthAnova)) # The closer the dots follow the line, the more likely the data is normal

shapiro.test(residuals(LengthAnova)) # if the p-value is below 0.05, it's probably not normal
leveneTest(percentLengthPedByIsolate ~ isolate * pedigree, dataForTransformation) # if the p-value is below 0.05, we should probably transform it...
# BUT, if all other normality tests are greater than 0.05, this 
# becomes subjective... we don't have to. 

plot(fitted(LengthAnova), residuals(LengthAnova), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red") # With this plot we want the dots to seem random. If there is an obvious pattern or a cluster,
# that would indicate there is some factor influencing some of our data.
# Ideally, things will be somewhat clustered around the red "0" line.
# If anything is really 'out there' from everything else, it might be a statistical outlier 
# and maybe it should be discarded (this is subjective!)

## Because of non-normality, we will transform the data

LengthMin <- dataForTransformation$percentLengthPedByIsolate + 1

boxcox_result <- boxcox(LengthMin ~ 1, data = dataForTransformation, plotit = TRUE)
lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]

dataForTransformation$percentLengthPedByIsolate_transformed <- if (lambda_optimal == 0) log(dataForTransformation$percentLengthPedByIsolate) else (dataForTransformation$percentLengthPedByIsolate^lambda_optimal - 1) / lambda_optimal

histMax <- max(dataForTransformation$percentLengthPedByIsolate_transformed)+1.5
histBy <- max(dataForTransformation$percentLengthPedByIsolate_transformed)/50
hist(dataForTransformation$percentLengthPedByIsolate, 
     breaks = seq(-0.5, 500, by = 5.5), #Adjust these until you have the fidelity you want
     main = "Percent Lengthes", 
     xlab = "", 
     col = "darkred", 
     border = "white")

LengthAnova <- lmer(percentLengthPedByIsolate_transformed ~ isolate * pedigree + (1|block), dataForTransformation)
summary(LengthAnova)
Anova(LengthAnova, Type = "III")

qqnorm(residuals(LengthAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(LengthAnova))

########## NOW DISEASE!


NecrosisAnova <- lmer(lesionPercent ~ isolate * pedigree + (1|block), dataForTransformation)
summary(NecrosisAnova)
Anova(NecrosisAnova, Type = "III")

qqnorm(residuals(NecrosisAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(NecrosisAnova)) # The closer the dots follow the line, the more likely the data is normal

shapiro.test(residuals(NecrosisAnova)) # if the p-value is below 0.05, it's probably not normal
leveneTest(lesionPercent ~ isolate * pedigree, dataForTransformation) # if the p-value is below 0.05, we should probably transform it...
# BUT, if all other normality tests are greater than 0.05, this 
# becomes subjective... we don't have to. 

plot(fitted(NecrosisAnova), residuals(NecrosisAnova), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red") # With this plot we want the dots to seem random. If there is an obvious pattern or a cluster,
# that would indicate there is some factor influencing some of our data.
# Ideally, things will be somewhat clustered around the red "0" line.
# If anything is really 'out there' from everything else, it might be a statistical outlier 
# and maybe it should be discarded (this is subjective!)

## Because of non-normality, we will transform the data

shifted_lesion_percent <- dataForTransformation$lesionPercent + 1

boxcox_result <- boxcox(shifted_lesion_percent ~ 1, data = dataForTransformation, plotit = TRUE)
lambda_optimal <- boxcox_result$x[which.max(boxcox_result$y)]

dataForTransformation$lesionPercent_transformed <- if (lambda_optimal == 0) log(dataForTransformation$lesionPercent) else (dataForTransformation$lesionPercent^lambda_optimal - 1) / lambda_optimal

histMax <- max(dataForTransformation$lesionPercent_transformed)+1.5
histBy <- max(dataForTransformation$lesionPercent_transformed)/50
hist(dataForTransformation$lesionPercent, 
     breaks = seq(-0.5, 500, by = 5.5), #Adjust these until you have the fidelity you want
     main = "Percent Necrosises", 
     xlab = "", 
     col = "darkred", 
     border = "white")

NecrosisAnova <- lmer(lesionPercent_transformed ~ isolate * pedigree + (1|block), dataForTransformation)
summary(NecrosisAnova)
Anova(NecrosisAnova, Type = "III")

shapiro.test(residuals(NecrosisAnova))
qqnorm(residuals(NecrosisAnova)) # This is a QQ plot. It is the most common test of normality
qqline(residuals(NecrosisAnova))

#############################################################################################################
############################################## TRANSFORMED ANOVAS ###########################################
#############################################################################################################
summary(aov(percentMassPedByIsolate_transformed ~ isolate * pedigree + block, dataForTransformation))
summary(aov(percentLengthPedByIsolate_transformed ~ isolate * pedigree + block, dataForTransformation))
summary(aov(lesionPercent_transformed ~ isolate * pedigree + block, dataForTransformation))
#############################################################################################################
#############################################################################################################
#############################################################################################################

# groupVariances <- tapply(dataPuvu$averageOospores, dataPuvu$isolate, var) # we will look at this when we have more data
# hist(groupVariances, xlab = "Variance", main = "Histogram of Variances", breaks = 20, col = fillTwo)

##### With our data checked we can do more analysis and move into graphs:

#TukeyHSD(LengthAnova)

if("LesionLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-massLetters)
}


########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
################################# BUILDING GRAPHS WITH TRANSFORMED SIGNIFICANCE ############
########################################################################################################################################################################################
########################################################################################################################################################################################
paDataTransformed <- dataForTransformation[grepl("G.a", dataForTransformation$isolateSpecies),]
piDataTransformed <- dataForTransformation[grepl("G.i", dataForTransformation$isolateSpecies),]
puuDataTransformed <- dataForTransformation[grepl("G.uu", dataForTransformation$isolateSpecies),]
########################################################################################################################################################################################

if("LesionLetters" %in% colnames(paDataTransformed)) {
  paDataTransformed <- paDataTransformed %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(paDataTransformed)) {
  paDataTransformed <- paDataTransformed %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(paDataTransformed)) {
  paDataTransformed <- paDataTransformed %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(piDataTransformed)) {
  piDataTransformed <- piDataTransformed %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(piDataTransformed)) {
  piDataTransformed <- piDataTransformed %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(piDataTransformed)) {
  piDataTransformed <- piDataTransformed %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(puuDataTransformed)) {
  puuDataTransformed <- puuDataTransformed %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(puuDataTransformed)) {
  puuDataTransformed <- puuDataTransformed %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(puuDataTransformed)) {
  puuDataTransformed <- puuDataTransformed %>%
    dplyr::select(-massLetters)
}



summary(aov(lesionPercent_transformed ~ pedigree + isolateSpecies + repetition + Error(block/pedigree), dataForTransformation))
summary(aov(percentLengthPedByIsolate_transformed ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), dataForTransformation))
summary(aov(percentMassPedByIsolate_transformed ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), dataForTransformation))


anovaFile <- paste0(anovaDir,"/Transformed_anova_summaries_ControlRemoved.txt",sep="")

# Use capture.output to save each summary to the file
capture.output({
  cat("Summary of ANOVA for lesionPercent:\n")
  print(summary(aov(lesionPercent_transformed ~ pedigree + isolateSpecies + repetition + Error(block/pedigree), dataForTransformation)))
  
  cat("\n\nSummary of ANOVA for percentLengthPedByIsolate:\n")
  print(summary(aov(percentLengthPedByIsolate_transformed ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), dataForTransformation)))
  
  cat("\n\nSummary of ANOVA for percentMassPedByIsolate:\n")
  print(summary(aov(percentMassPedByIsolate_transformed ~ pedigree * isolateSpecies + repetition + Error(block/pedigree), dataForTransformation)))
}, file = anovaFile)

if("LesionLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-massLetters)
}



model_lesion <- aov(lesionPercent_transformed ~ isolateSpecies, dataForTransformation)
model_length <- aov(percentLengthPedByIsolate_transformed ~ isolateSpecies, dataForTransformation)
model_mass <- aov(percentMassPedByIsolate_transformed ~ isolateSpecies, dataForTransformation)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$isolateSpecies$Letters)
letters_lesionDf$isolateSpecies <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
dataForTransformation <- merge(dataForTransformation, letters_lesionDf, by = "isolateSpecies")
summaryData <- dataForTransformation %>%
  group_by(inocMethod, isolateSpecies) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(dataForTransformation) +
  aes(
    x = isolateSpecies,
    y = lesionPercent,
    fill = isolateSpecies
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Disease Severity (Percent)" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = isolateSpecies, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Treatment Disease Severity (Percent)_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Disease Severity (Percent)_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$isolateSpecies$Letters)
letters_lengthDf$isolateSpecies <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
dataForTransformation <- merge(dataForTransformation, letters_lengthDf, by = "isolateSpecies")
summaryData <- dataForTransformation %>%
  group_by(inocMethod, isolateSpecies) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate_transformed = median(percentLengthPedByIsolate_transformed), lengthLetters = first(lengthLetters))

ggplot(dataForTransformation) +
  aes(
    x = isolateSpecies,
    y = percentLengthPedByIsolate,
    fill = isolateSpecies
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Root Length as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = isolateSpecies, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Treatment Length Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Length Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$isolateSpecies$Letters)
letters_massDf$isolateSpecies <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
dataForTransformation <- merge(dataForTransformation, letters_massDf, by = "isolateSpecies")
summaryData <- dataForTransformation %>%
  group_by(inocMethod, isolateSpecies) %>%
  dplyr::summarize(medianpercentMassPedByIsolate_transformed = median(percentMassPedByIsolate_transformed), massLetters = first(massLetters))

ggplot(dataForTransformation) +
  aes(
    x = isolateSpecies,
    y = percentMassPedByIsolate,
    fill = isolateSpecies
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Inoculation Species",  
    y = "Root Mass as % of Standards" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = isolateSpecies, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Treatment Mass Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Treatment Mass Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth/2.5, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


######################### IS THERE A DIFFERENCE IN ROOT ROTS? ##################################################################################################

dataForTransformation <- dataForTransformation %>%
  filter(isolate != "Standard")


if("LesionLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-massLetters)
}

model_lesion <- aov(lesionPercent_transformed ~ pedigree+repetition, dataForTransformation)
model_length <- aov(percentLengthPedByIsolate_transformed ~ pedigree+repetition, dataForTransformation)
model_mass <- aov(percentMassPedByIsolate_transformed ~ pedigree+repetition, dataForTransformation)

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
dataForTransformation <- merge(dataForTransformation, letters_lesionDf, by = "pedigree")
summaryData <- dataForTransformation %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent), LesionLetters = first(LesionLetters))

ggplot(dataForTransformation) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Treated Maize Lines",  
    y = "Disease Severity (Percent), Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Pedigree PRR Disease Severity (Percent)_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree PRR Disease Severity (Percent)_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
dataForTransformation <- merge(dataForTransformation, letters_lengthDf, by = "pedigree")
summaryData <- dataForTransformation %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate_transformed = median(percentLengthPedByIsolate_transformed), lengthLetters = first(lengthLetters))

ggplot(dataForTransformation) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Treated Maize Lines",  
    y = "Root Length as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Pedigree PRR Length Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree PRR Length Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
dataForTransformation <- merge(dataForTransformation, letters_massDf, by = "pedigree")
summaryData <- dataForTransformation %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate_transformed = median(percentMassPedByIsolate_transformed), massLetters = first(massLetters))

ggplot(dataForTransformation) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Treated Maize Lines",  
    y = "Root Mass as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Pedigree PRR Mass Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Pedigree PRR Mass Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)


############ INDIVIDUAL SPECIES ANALYSIS #######################################################################################################################################


##aphanidermatum only##########################################################################################################################################################

model_lesion <- aov(lesionPercent_transformed ~ pedigree+repetition, paDataTransformed)
model_length <- aov(percentLengthPedByIsolate_transformed ~ pedigree+repetition, paDataTransformed)
model_mass <- aov(percentMassPedByIsolate_transformed ~ pedigree+repetition, paDataTransformed)

summary(aov(lesionPercent_transformed ~ pedigree + block, paDataTransformed))
summary(aov(percentLengthPedByIsolate_transformed ~ pedigree + block, paDataTransformed))
summary(aov(percentMassPedByIsolate_transformed ~ pedigree + block, paDataTransformed))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
paDataTransformed <- merge(paDataTransformed, letters_lesionDf, by = "pedigree")
summaryData <- paDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(paDataTransformed) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. aphanidermatum",  
    y = "Disease Severity (Percent), Transformed Dataset" 
  ) +  
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Aphanidermatum Root Rot Disease Severity (Percent)_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Aphanidermatum Root Rot Disease Severity (Percent)_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
paDataTransformed <- merge(paDataTransformed, letters_lengthDf, by = "pedigree")
summaryData <- paDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate_transformed = median(percentLengthPedByIsolate_transformed), lengthLetters = first(lengthLetters))

ggplot(paDataTransformed) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. aphanidermatum",  
    y = "Root Length, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Aphanidermatum Root Rot Length Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Aphanidermatum Root Rot Length Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
paDataTransformed <- merge(paDataTransformed, letters_massDf, by = "pedigree")
summaryData <- paDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate_transformed = median(percentMassPedByIsolate_transformed), massLetters = first(massLetters))

ggplot(paDataTransformed) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. aphanidermatum",  
    y = "Root Mass as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Aphanidermatum Root Rot Mass Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Aphanidermatum Root Rot Mass Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##irregulare only##############################################################################################################################################################

model_lesion <- aov(lesionPercent_transformed ~ pedigree+repetition, piDataTransformed)
model_length <- aov(percentLengthPedByIsolate_transformed ~ pedigree+repetition, piDataTransformed)
model_mass <- aov(percentMassPedByIsolate_transformed ~ pedigree+repetition, piDataTransformed)

summary(aov(lesionPercent_transformed ~ pedigree + block, piDataTransformed))
summary(aov(percentLengthPedByIsolate_transformed ~ pedigree + block, piDataTransformed))
summary(aov(percentMassPedByIsolate_transformed ~ pedigree + block, piDataTransformed))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
piDataTransformed <- merge(piDataTransformed, letters_lesionDf, by = "pedigree")
summaryData <- piDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(piDataTransformed) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. irregulare",  
    y = "Disease Severity (Percent), Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Irregulare Root Rot Disease Severity (Percent)_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Irregulare Root Rot Disease Severity (Percent)_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
piDataTransformed <- merge(piDataTransformed, letters_lengthDf, by = "pedigree")
summaryData <- piDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate_transformed = median(percentLengthPedByIsolate_transformed), lengthLetters = first(lengthLetters))

ggplot(piDataTransformed) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. irregulare",  
    y = "Root Length as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Irregulare Root Rot Length Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Irregulare Root Rot Length Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
piDataTransformed <- merge(piDataTransformed, letters_massDf, by = "pedigree")
summaryData <- piDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate_transformed = median(percentMassPedByIsolate_transformed), massLetters = first(massLetters))

ggplot(piDataTransformed) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. irregulare",  
    y = "Root Mass as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Irregulare Root Rot Mass Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Irregulare Root Rot Mass Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

##ultimum ultimum only##############################################################################################################################################################

model_lesion <- aov(lesionPercent_transformed ~ pedigree+repetition, puuDataTransformed)
model_length <- aov(percentLengthPedByIsolate_transformed ~ pedigree+repetition, puuDataTransformed)
model_mass <- aov(percentMassPedByIsolate_transformed ~ pedigree+repetition, puuDataTransformed)

summary(aov(lesionPercent_transformed ~ pedigree + block, puuDataTransformed))
summary(aov(percentLengthPedByIsolate_transformed ~ pedigree + block, puuDataTransformed))
summary(aov(percentMassPedByIsolate_transformed ~ pedigree + block, puuDataTransformed))

tukey_lesion <- TukeyHSD(model_lesion)
tukey_length <- TukeyHSD(model_length)
tukey_mass <- TukeyHSD(model_mass)

letters_lesion <- multcompLetters4(model_lesion, tukey_lesion)
letters_length <- multcompLetters4(model_length, tukey_length)
letters_mass <- multcompLetters4(model_mass, tukey_mass)

letters_lesionDf <- as.data.frame(letters_lesion$pedigree$Letters)
letters_lesionDf$pedigree <- rownames(letters_lesionDf)
colnames(letters_lesionDf)[1] <- "LesionLetters"
rownames(letters_lesionDf) <- NULL
puuDataTransformed <- merge(puuDataTransformed, letters_lesionDf, by = "pedigree")
summaryData <- puuDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianLesionPercent = median(lesionPercent_transformed), LesionLetters = first(LesionLetters))

ggplot(puuDataTransformed) +
  aes(
    x = pedigree,
    y = lesionPercent,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. ultimum var. ultimum",  
    y = "Disease Severity (Percent), Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = LesionLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Ultimum Root Rot Disease Severity (Percent)_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Ultimum Root Rot Disease Severity (Percent)_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_lengthDf <- as.data.frame(letters_length$pedigree$Letters)
letters_lengthDf$pedigree <- rownames(letters_lengthDf)
colnames(letters_lengthDf)[1] <- "lengthLetters"
rownames(letters_lengthDf) <- NULL
puuDataTransformed <- merge(puuDataTransformed, letters_lengthDf, by = "pedigree")
summaryData <- puuDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentLengthPedByIsolate_transformed = median(percentLengthPedByIsolate_transformed), lengthLetters = first(lengthLetters))

ggplot(puuDataTransformed) +
  aes(
    x = pedigree,
    y = percentLengthPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. ultimum var. ultimum",  
    y = "Root Length as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = lengthLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Ultimum Root Rot Length Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Ultimum Root Rot Length Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

letters_massDf <- as.data.frame(letters_mass$pedigree$Letters)
letters_massDf$pedigree <- rownames(letters_massDf)
colnames(letters_massDf)[1] <- "massLetters"
rownames(letters_massDf) <- NULL
puuDataTransformed <- merge(puuDataTransformed, letters_massDf, by = "pedigree")
summaryData <- puuDataTransformed %>%
  group_by(inocMethod, pedigree) %>%
  dplyr::summarize(medianpercentMassPedByIsolate_transformed = median(percentMassPedByIsolate_transformed), massLetters = first(massLetters))

ggplot(puuDataTransformed) +
  aes(
    x = pedigree,
    y = percentMassPedByIsolate,
    fill = pedigree
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey() +
  geom_point(shape = 25, color = "black", fill = NA, size = 1, alpha = 1, position = position_jitter(width = 0.2)) +
  labs(
    x = "Maize Lines Inoculated with G. ultimum var. ultimum",  
    y = "Root Mass as % of Standards, Transformed Dataset" 
  ) +
  # Use the dplyr::summarized data to plot the letters
  geom_text(
    data = summaryData, 
    aes(
      x = pedigree, 
      y = -Inf,  # Adjust to place the letters above the boxplot
      label = massLetters
    ), 
    position = position_dodge(width = 0.75), 
    vjust = -0.1
  ) +
  theme_minimal() +
  theme( text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank() 
  )
filename <- paste0("Ultimum Root Rot Mass Percent of Standards_transformed.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Ultimum Root Rot Mass Percent of Standards_transformed.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth, height = ChartHeight, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
model_lesion <- aov(lesionPercent_transformed ~ pedigree+repetition, puuDataTransformed)
model_length <- aov(percentLengthPedByIsolate_transformed ~ pedigree+repetition, puuDataTransformed)
model_mass <- aov(percentMassPedByIsolate_transformed ~ pedigree+repetition, puuDataTransformed)

if("LesionLetters" %in% colnames(paDataTransformed)) {
  paDataTransformed <- paDataTransformed %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(paDataTransformed)) {
  paDataTransformed <- paDataTransformed %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(paDataTransformed)) {
  paDataTransformed <- paDataTransformed %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(piDataTransformed)) {
  piDataTransformed <- piDataTransformed %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(piDataTransformed)) {
  piDataTransformed <- piDataTransformed %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(piDataTransformed)) {
  piDataTransformed <- piDataTransformed %>%
    dplyr::select(-massLetters)
}

if("LesionLetters" %in% colnames(puuDataTransformed)) {
  puuDataTransformed <- puuDataTransformed %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(puuDataTransformed)) {
  puuDataTransformed <- puuDataTransformed %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(puuDataTransformed)) {
  puuDataTransformed <- puuDataTransformed %>%
    dplyr::select(-massLetters)
}

########################################################################################################################################################################################
## ANOVAs with some factors removed
########################################################################################################################################################################################
if("LesionLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-LesionLetters)
}

if("lengthLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-lengthLetters)
}

if("massLetters" %in% colnames(dataForTransformation)) {
  dataForTransformation <- dataForTransformation %>%
    dplyr::select(-massLetters)
}

########################################################################################################################################################################################
######################################################################## CORRELATIONS ##################################################################################################
########################################################################################################################################################################################

############################################ MATRIX OF CORRELATIONS VIA PEARSON CORRELATION ###############################################
# Define the variable names

library("Hmisc")
library("reshape2")

vars <- c("% Root Necrosis","Percent Root Length","Percent Root Mass")

matrixCor <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
matrixP <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))

rownames(matrixCor) <- vars
colnames(matrixCor) <- vars
# Defaults
diag(matrixCor) <- 1
diag(matrixP) <- NA
#

################################################################################################################
model <- cor.test(
  NoControlData$lesionPercent, 
  NoControlData$percentLengthPedByIsolate, 
  method = "pearson"
)

matrixCor["% Root Necrosis", "Percent Root Length"] <- model$estimate
matrixCor["Percent Root Length", "% Root Necrosis"] <- model$estimate
matrixP["% Root Necrosis", "Percent Root Length"] <- model$p.value
matrixP["Percent Root Length", "% Root Necrosis"] <- model$p.value
#

model <- cor.test(NoControlData$percentMassPedByIsolate, NoControlData$percentLengthPedByIsolate, method = "pearson")

matrixCor["Percent Root Mass", "Percent Root Length"] <- model$estimate 
matrixCor["Percent Root Length", "Percent Root Mass"] <- model$estimate 
matrixP["Percent Root Mass", "Percent Root Length"] <- model$p.value
matrixP["Percent Root Length", "Percent Root Mass"] <- model$p.value
#

model <- cor.test(NoControlData$percentMassPedByIsolate, NoControlData$lesionPercent, method = "pearson")

matrixCor["Percent Root Mass", "% Root Necrosis"] <- model$estimate 
matrixCor["% Root Necrosis", "Percent Root Mass"] <- model$estimate 
matrixP["Percent Root Mass", "% Root Necrosis"] <- model$p.value
matrixP["% Root Necrosis", "Percent Root Mass"] <- model$p.value
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
  geom_text(aes(label = sprintf("r = %.2f\np = %.3f", Correlation, PValue)), size = 5, vjust = .45) +
  labs(
    x = " ",  
    y = " " 
  ) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(angle = 45, vjust = 1, face = "bold")) +
  labs(title = "Correlation and P-Values")

filename <- paste0("Correlation and P-Values.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Correlation and P-Values.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

