################################################################################
## Globisporangium Root Rot Correlations Script, Winter 24-25  
## Analysis by Harrison Hall (harrisonpiercehall@gmail.com or hphall2@illinois.edu)
################################################################################

## Libraries used 
library("dplyr") ## Contains various helpful functions for filtering datasets
library("tidyr") ## More functions for filtering datasets
library("tidyverse")   ## Adds ggplot2 and other data manipulation tools
library("ggplot2") ## Some functions look to see if ggplot2 is actually in the library or not
library("ggthemes") ## More control over ggplots with more options
library("esquisse") ## GUI for ggplots to look at data quickly
library("ggcorrplot")## More ggplot options for correlations
library("lmtest") ## For testing linear regression models
library("agricolae") ## Statistical tests and graphics 
library("readxl") ## To open excel files easily
library("reshape2") ## To reshape datasets
library("corrplot") ## For correlation matrices and confidence intervals
# library("ggpubr")

ProjDir <- rstudioapi::getActiveDocumentContext()$path # Gets path of script
setwd(dirname(ProjDir)) # Sets the working directory to the same location as the script location
ResultDir <- file.path(dirname(ProjDir), "Results") 
dir.create(ResultDir, recursive = TRUE)
fTime <- format(Sys.time(), "%y-%m-%d %H%M")
AnalysisDir <- file.path(dirname(ResultDir), "data_analysis", fTime, sep="") 
dir.create(AnalysisDir, recursive = TRUE)

##### Get files and folders#####################################################
dfDir <- rstudioapi::selectFile(caption="Choose dataset to work with", path=ProjDir) #To manually select data
data <- read_excel(dfDir) # Reads if it is a .xlsx file. If it is a csv, change read_excel to read_csv
dataName <- tools::file_path_sans_ext(basename(dfDir)) 
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
fontSize <- 7 
# Set the default exported chart size to >5000x2000 units, this is print quality
ChartWidth <- (5760)/5
ChartHeight <- (3240)/2

################################################################################
##### Beginning of analysis ####################################################
################################################################################

head(data) # Take a peek at the data and make sure we grabbed the right file
str(data) # This shows us our column names and examples of data. We will want to make sure things are formatted the way we want.

data$pedigree <- as.factor(data$pedigree) # Maize lines
data$experiment <- as.factor(data$experiment) # ID for each experiment and repetition
data$experimentType <- as.factor(data$experimentType) # Whether it was a single species or multiple species
data$sampleID <- as.factor(data$sampleID) # Just numbers to separate each row
data$isolate <- as.factor(data$isolate) # Species used in treatment
data$block <- as.factor(data$block) # Blocking group
#Numerics:
data$lesionPercent <- as.numeric(data$lesionPercent) # Disease rating
data$percentLengthPedByIsolate <- as.numeric(data$percentLengthPedByIsolate)*100 # Percent of Length when divided by average of standards
data$percentMassPedByIsolate <- as.numeric(data$percentMassPedByIsolate)*100 # Percent of Mass when divided by average of standards

originalData <- data #isolating data so it is unaltered

filteredData <- data[!grepl("Standard", data$isolate),]
ogData <- filteredData
## Lines that are consistent: CML277, CML103, H100, CML333, Oh7B

# Filter the data to include only the specified pedigrees

warmExpData <- filteredData[grepl("multiWarm", filteredData$experimentType),] # Filtered to just use warm experiments

stats <- warmExpData %>%
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
filteredData <- warmExpData %>%
  left_join(stats, by = "pedigree") %>%
  filter(
    lesionPercent <= meanLesionPercent + 2.0 * sdLesionPercent,
    percentLengthPedByIsolate <= meanPercentLength + 2.0 * sdPercentLength,
    percentMassPedByIsolate <= meanPercentMass + 2.0 * sdPercentMass
  )

# Select columns to avoid clutter
warmData <- filteredData %>%
  select(-starts_with("mean"), -starts_with("sd"))

# View the filtered dataset
print(warmData)

########################################################################################################################################################################################
############################################################## CORRELATIONS BY REPETITION ##################################################################################################
########################################################################################################################################################################################

########################################################################## WARM

##### Guu v Ga #####
#################### DISEASE
df_Guu_Ga <- warmData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 4, alpha = 1) +
  geom_text_repel(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree",
       y = "% Disease of Ga for each pedigree",
       title = "Warm Experiments % Disease: Guu x Ga") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Disease between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Ga <- warmData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Ga for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Guu and Ga across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Length of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)


#################### MASS
df_Guu_Ga <- warmData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by( pedigree, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Ga for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Guu and Ga across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Mass of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)


##### Guu v Gi #####
#################### DISEASE
df_Guu_Gi <- warmData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by( pedigree, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Warm Experiments Correlation of % Disease between Guu and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Disease between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Gi <- warmData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Guu and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Length of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Guu_Gi <- warmData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree,  isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Guu and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Mass of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

##### Ga v Gi #####
#################### DISEASE
df_Ga_Gi <- warmData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree,  isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Ga for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Warm Experiments Correlation of % Disease between Ga and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Disease between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Ga_Gi <- warmData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree,  isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Ga for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Ga and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Length of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Ga_Gi <- warmData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree,  isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Ga for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Ga and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Warm Experiments Mass of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)



########################################################################## COLD
##### By Maize Line ###########################################################
coldExpData <- data[grepl("multiCold", data$experimentType),] # Filtered out experiments that just used Guu

# Calculate mean and standard deviation by pedigree
stats <- coldExpData %>%
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
coldFilteredData <- coldExpData %>%
  left_join(stats, by = "pedigree") %>%
  filter(
    lesionPercent <= meanLesionPercent + 2.0 * sdLesionPercent,
    percentLengthPedByIsolate <= meanPercentLength + 2.0 * sdPercentLength,
    percentMassPedByIsolate <= meanPercentMass + 2.0 * sdPercentMass
  )

# Select columns to avoid clutter
coldData <- coldFilteredData %>%
  select(-starts_with("mean"), -starts_with("sd"))

# View the filtered dataset
print(coldData)


##### Guu v Ga #####
#################### DISEASE
df_Guu_Ga <- coldData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree",
       y = "% Disease of Ga for each pedigree",
       title = "Cold Experiments Correlation of % Disease between Guu and Ga across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Disease between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Ga <- coldData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 4, alpha = 1) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Ga for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Guu and Ga across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Length of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)


#################### MASS
df_Guu_Ga <- coldData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Ga for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Guu and Ga across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Mass of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)


##### Guu v Gi #####
#################### DISEASE
df_Guu_Gi <- coldData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree,isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Cold Experiments Correlation of % Disease between Guu and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Disease between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Gi <- coldData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Guu and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Length of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Guu_Gi <- coldData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Guu and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Mass of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

##### Ga v Gi #####
#################### DISEASE
df_Ga_Gi <- coldData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Ga for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Cold Experiments Correlation of % Disease between Ga and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Disease between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Ga_Gi <- coldData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Ga for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Ga and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Length of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Ga_Gi <- coldData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 4, alpha = 1 ) + 
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),  
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Ga for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Ga and Gi across Maize Lines") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Cold Experiments Mass of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

########################################################################################################################################################################################
############################################################ Warm Vs Cold Environments ##################################################################################################
########################################################################################################################################################################################

corFiltData <- data[!grepl("singleSpecies", data$experimentType),] # Filtered out experiments that just used Guu

# Calculate mean and standard deviation by pedigree
stats <- corFiltData %>%
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
corFilteredData <- corFiltData %>%
  left_join(stats, by = "pedigree") %>%
  filter(
    lesionPercent <= meanLesionPercent + 2.0 * sdLesionPercent,
    percentLengthPedByIsolate <= meanPercentLength + 2.0 * sdPercentLength,
    percentMassPedByIsolate <= meanPercentMass + 2.0 * sdPercentMass
  )

# Select columns to avoid clutter
corData <- corFilteredData %>%
  select(-starts_with("mean"), -starts_with("sd"))

# View the filtered dataset
print(corData)

##### Guu v Ga #####
#################### DISEASE
df_Guu_Ga <- corData %>%
  group_by(pedigree, experimentType) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experimentType, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$multiCold, df_Guu_Ga$multiWarm)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = multiCold, y = multiWarm, color = pedigree)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of PRR in Cold Experiments",
       y = "% Disease of PRR in Warm Experiments",
       title = "Correlation of % Disease between Pythium root rots across warm and cold experiments") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Disease between Pythium root rots across warm and cold.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Disease between Pythium root rots across warm and cold.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)

df_Guu_Ga <- corData %>%
  filter(isolate %in% c("Guu")) %>%
  group_by( pedigree, experimentType) %>% ##### PICK UP HERE
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experimentType, values_from = lesionPercent)# %>%
#drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$multiCold, df_Guu_Ga$multiWarm)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = multiCold, y = multiWarm, color = pedigree)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu in Cold Experiments",
       y = "% Disease of Guu in Warm Experiments",
       title = "Correlation of % Disease between Guu across warm and cold experiments") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Disease between Guu across warm and cold.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Disease between Guu across warm and cold.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)



df_Guu_Ga <- corData %>%
  filter(isolate %in% c("Ga")) %>%
  group_by( pedigree, experimentType) %>% ##### PICK UP HERE
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experimentType, values_from = lesionPercent)# %>%
#drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$multiCold, df_Guu_Ga$multiWarm)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = multiCold, y = multiWarm, color = pedigree)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Ga in Cold Experiments",
       y = "% Disease of Ga in Warm Experiments",
       title = "Correlation of % Disease between Ga across warm and cold experiments") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Disease between Ga across warm and cold.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Disease between Ga across warm and cold.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)






df_Guu_Ga <- corData %>%
  filter(isolate %in% c("Gi")) %>%
  group_by( pedigree, experimentType) %>% 
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experimentType, values_from = lesionPercent)# %>%
#drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$multiCold, df_Guu_Ga$multiWarm)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = multiCold, y = multiWarm, color = pedigree)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = paste(pedigree)), vjust = 1.15, hjust = 0.5, color = "black", fontface = "bold", size = 4) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Gi in Cold Experiments",
       y = "% Disease of Gi in Warm Experiments",
       title = "Correlation of % Disease between Gi across warm and cold experiments") +
  theme_minimal() +
  theme(legend.position = "none") 

filename <- paste0("Disease between Gi across warm and cold.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Disease between Gi across warm and cold.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2.4, height = ChartHeight*2, units = "px", dpi = 400, bg = "#FFFEFE", limitsize = FALSE)


########################################################################################################################################################################################
######################################################################## EXPERIMENT CORRELATIONS ##################################################################################################
########################################################################################################################################################################################


wide_data <- corData %>%
  filter(!is.na(percentLengthPedByIsolate)) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarize(meanPercentLengthPed = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = meanPercentLengthPed)


correlation_matrix <- cor(wide_data[,-1], use = "pairwise.complete.obs", method = "pearson")
corrplot(correlation_matrix, method = "circle")

vars <- names(wide_data)[-1] 
# Initialize the matrices
matrixCor <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
matrixP <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
# Set default values
diag(matrixCor) <- 1  # Set diagonal values to 1 for correlation
diag(matrixP) <- NA   # NA for diagonal in P-value matrix

for (i in 1:length(vars)) {
  for (j in i:length(vars)) {
    if (i != j) {
      test <- cor.test(wide_data[[vars[i]]], wide_data[[vars[j]]], method = "pearson")
      matrixCor[vars[i], vars[j]] <- test$estimate
      matrixCor[vars[j], vars[i]] <- test$estimate  # Symmetric
      matrixP[vars[i], vars[j]] <- test$p.value
      matrixP[vars[j], vars[i]] <- test$p.value
    }
  }
}

# Melt the matrices
meltedCor <- melt(matrixCor, varnames = c("Var1", "Var2"), value.name = "Correlation")
meltedP <- melt(matrixP, varnames = c("Var1", "Var2"), value.name = "PValue")

# Combine the melted data frames
combinedMelted <- cbind(meltedCor, PValue = meltedP$PValue)

combinedMeltedLower <- combinedMelted[as.numeric(combinedMelted$Var1) >= as.numeric(combinedMelted$Var2), ]

ggplot(combinedMeltedLower, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("r = %.2f\np = %.3f", Correlation, PValue)), size = 5, vjust = 0.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(angle = 45, vjust = 1, face = "bold")) +
  labs(title = "Length Correlation and P-Values between Experiments", x = "", y = "")

filename <- paste0("Experiment Length Correlation and P-Values.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*8, height = ChartHeight*2, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Experiment Length Correlation and P-Values.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*8, height = ChartHeight*2, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)




wide_data <- corData %>%
  filter(!is.na(percentMassPedByIsolate)) %>%
  group_by( pedigree, experiment) %>%
  dplyr::summarize(meanPercentMassPed = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = meanPercentMassPed)


correlation_matrix <- cor(wide_data[,-1], use = "pairwise.complete.obs", method = "pearson")
corrplot(correlation_matrix, method = "circle")

vars <- names(wide_data)[-1] 
# Initialize the matrices
matrixCor <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
matrixP <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
# Set default values
diag(matrixCor) <- 1  # Set diagonal values to 1 for correlation
diag(matrixP) <- NA   # NA for diagonal in P-value matrix

for (i in 1:length(vars)) {
  for (j in i:length(vars)) {
    if (i != j) {
      test <- cor.test(wide_data[[vars[i]]], wide_data[[vars[j]]], method = "pearson")
      matrixCor[vars[i], vars[j]] <- test$estimate
      matrixCor[vars[j], vars[i]] <- test$estimate  # Symmetric
      matrixP[vars[i], vars[j]] <- test$p.value
      matrixP[vars[j], vars[i]] <- test$p.value
    }
  }
}

# Melt the matrices
meltedCor <- melt(matrixCor, varnames = c("Var1", "Var2"), value.name = "Correlation")
meltedP <- melt(matrixP, varnames = c("Var1", "Var2"), value.name = "PValue")

# Combine the melted data frames
combinedMelted <- cbind(meltedCor, PValue = meltedP$PValue)

combinedMeltedLower <- combinedMelted[as.numeric(combinedMelted$Var1) >= as.numeric(combinedMelted$Var2), ]

ggplot(combinedMeltedLower, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("r = %.2f\np = %.3f", Correlation, PValue)), size = 5, vjust = 0.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(angle = 45, vjust = 1, face = "bold")) +
  labs(title = "Mass Correlation and P-Values between Experiments", x = "", y = "")

filename <- paste0("Experiment Mass Correlation and P-Values.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*8, height = ChartHeight*2, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Experiment Mass Correlation and P-Values.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*8, height = ChartHeight*2, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)




wide_data <- corData %>%
  filter(!is.na(lesionPercent)) %>%
  group_by( pedigree, experiment) %>%
  dplyr::summarize(meanPercentDiseasePed = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = meanPercentDiseasePed)


correlation_matrix <- cor(wide_data[,-1], use = "pairwise.complete.obs", method = "pearson")
corrplot(correlation_matrix, method = "circle")

vars <- names(wide_data)[-1] 
# Initialize the matrices
matrixCor <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
matrixP <- matrix(nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
# Set default values
diag(matrixCor) <- 1  # Set diagonal values to 1 for correlation
diag(matrixP) <- NA   # NA for diagonal in P-value matrix

for (i in 1:length(vars)) {
  for (j in i:length(vars)) {
    if (i != j) {
      test <- cor.test(wide_data[[vars[i]]], wide_data[[vars[j]]], method = "pearson")
      matrixCor[vars[i], vars[j]] <- test$estimate
      matrixCor[vars[j], vars[i]] <- test$estimate  # Symmetric
      matrixP[vars[i], vars[j]] <- test$p.value
      matrixP[vars[j], vars[i]] <- test$p.value
    }
  }
}

# Melt the matrices
meltedCor <- melt(matrixCor, varnames = c("Var1", "Var2"), value.name = "Correlation")
meltedP <- melt(matrixP, varnames = c("Var1", "Var2"), value.name = "PValue")

# Combine the melted data frames
combinedMelted <- cbind(meltedCor, PValue = meltedP$PValue)

combinedMeltedLower <- combinedMelted[as.numeric(combinedMelted$Var1) >= as.numeric(combinedMelted$Var2), ]

ggplot(combinedMeltedLower, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("r = %.2f\np = %.3f", Correlation, PValue)), size = 5, vjust = 0.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(angle = 45, vjust = 1, face = "bold")) +
  labs(title = "Disease Correlation and P-Values between Experiments", x = "", y = "")

filename <- paste0("Experiment Disease Correlation and P-Values.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*8, height = ChartHeight*2, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Experiment Disease Correlation and P-Values.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*8, height = ChartHeight*2, units = "px", dpi = 600, bg = "#FFFEFE", limitsize = FALSE)

