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
ChartWidth <- (5760)
ChartHeight <- (3240)

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

originalData <- data #isolating data so it is uncontaminated

ogData <- data[!grepl("Standard", data$isolate),] #Remove the control groups, use this going forward.

#esquisser(data) # To look at data
###### TEST OF REMOVAL OF TROUBLESOME LINE:

########################################################################################################################################################################################
############################################################## CORRELATIONS BY REPETITION ##################################################################################################
########################################################################################################################################################################################

########################################################################## WARM
warmData <- ogData[grepl("multiWarm", ogData$experimentType),] # Filtered out experiments that just used Guu

##### Guu v Ga #####
#################### DISEASE
df_Guu_Ga <- warmData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree",
       y = "% Disease of Ga for each pedigree",
       title = "Warm Experiments Correlation of % Disease between Guu and Ga across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Disease between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Ga <- warmData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Ga for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Guu and Ga across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Length of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


#################### MASS
df_Guu_Ga <- warmData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Ga for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Guu and Ga across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Mass of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


##### Guu v Gi #####
#################### DISEASE
df_Guu_Gi <- warmData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Warm Experiments Correlation of % Disease between Guu and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Disease between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Gi <- warmData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Guu and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Length of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Guu_Gi <- warmData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Guu and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Mass of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

##### Ga v Gi #####
#################### DISEASE
df_Ga_Gi <- warmData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Ga for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Warm Experiments Correlation of % Disease between Ga and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Disease between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Ga_Gi <- warmData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Ga for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Ga and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Length of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Ga_Gi <- warmData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Ga for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Ga and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Mass of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

##### By Repetition ############################################################
##### By Repetition ############################################################

##### Guu v Ga #####
#################### DISEASE
df_Guu_Comp <- warmData %>%
  filter(isolate == "Guu", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Comp)

cor_test_Guu <- cor.test(df_Guu_Comp$MSW1, df_Guu_Comp$MSW2)
r_value <- cor_test_Guu$estimate
p_value <- cor_test_Guu$p.value

ggplot(df_Guu_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of MSW1 for each pedigree",
       y = "% Disease of MSW2 for each pedigree",
       title = "Warm Experiments Correlation of % Disease between Guu Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Disease between Guu Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Guu Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Comp <- warmData %>%
  filter(isolate == "Guu", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Comp)

cor_test_Guu <- cor.test(df_Guu_Comp$MSW1, df_Guu_Comp$MSW2)
r_value <- cor_test_Guu$estimate
p_value <- cor_test_Guu$p.value

ggplot(df_Guu_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of MSW1 for each pedigree", 
       y = "% Length of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Guu Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Length of Standards between Guu Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Guu Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


#################### MASS
df_Guu_Comp <- warmData %>%
  filter(isolate == "Guu", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Comp)

cor_test_Guu <- cor.test(df_Guu_Comp$MSW1, df_Guu_Comp$MSW2)
r_value <- cor_test_Guu$estimate
p_value <- cor_test_Guu$p.value

ggplot(df_Guu_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of MSW1 for each pedigree", 
       y = "% Mass of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Guu Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Mass of Standards between Guu Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Guu Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


#####  Gi #####
#################### DISEASE
df_Gi_Comp <- warmData %>%
  filter(isolate == "Gi", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = lesionPercent) %>%
  drop_na()
print(df_Gi_Comp)

cor_test_Gi <- cor.test(df_Gi_Comp$MSW1, df_Gi_Comp$MSW2)
r_value <- cor_test_Gi$estimate
p_value <- cor_test_Gi$p.value

ggplot(df_Gi_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of MSW1 for each pedigree", 
       y = "% Disease of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of % Disease between Gi Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Disease between Gi Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Gi Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Gi_Comp <- warmData %>%
  filter(isolate == "Gi", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Gi_Comp)

cor_test_Gi <- cor.test(df_Gi_Comp$MSW1, df_Gi_Comp$MSW2)
r_value <- cor_test_Gi$estimate
p_value <- cor_test_Gi$p.value

ggplot(df_Gi_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of MSW1 for each pedigree", 
       y = "% Length of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Gi Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Length of Standards between Gi Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Gi Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Gi_Comp <- warmData %>%
  filter(isolate == "Gi", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Gi_Comp)

cor_test_Gi <- cor.test(df_Gi_Comp$MSW1, df_Gi_Comp$MSW2)
r_value <- cor_test_Gi$estimate
p_value <- cor_test_Gi$p.value

ggplot(df_Gi_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of MSW1 for each pedigree", 
       y = "% Mass of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Gi Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Mass of Standards between Gi Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Gi Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

##### Ga #####
#################### DISEASE
df_Ga_Comp <- warmData %>%
  filter(isolate == "Ga", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = lesionPercent) %>%
  drop_na()
print(df_Ga_Comp)

cor_test_Ga <- cor.test(df_Ga_Comp$MSW1, df_Ga_Comp$MSW2)
r_value <- cor_test_Ga$estimate
p_value <- cor_test_Ga$p.value

ggplot(df_Ga_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of MSW1 for each pedigree", 
       y = "% Disease of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of % Disease between Ga Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Disease between Ga Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Disease between Ga Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Ga_Comp <- warmData %>%
  filter(isolate == "Ga", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Ga_Comp)

cor_test_Ga <- cor.test(df_Ga_Comp$MSW1, df_Ga_Comp$MSW2)
r_value <- cor_test_Ga$estimate
p_value <- cor_test_Ga$p.value

ggplot(df_Ga_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of MSW1 for each pedigree", 
       y = "% Length of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of Percent Length of Standards between Ga Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Length of Standards between Ga Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Length of Standards between Ga Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Ga_Comp <- warmData %>%
  filter(isolate == "Ga", experiment %in% c("MSW1", "MSW2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Ga_Comp)

cor_test_Ga <- cor.test(df_Ga_Comp$MSW1, df_Ga_Comp$MSW2)
r_value <- cor_test_Ga$estimate
p_value <- cor_test_Ga$p.value

ggplot(df_Ga_Comp, aes(x = MSW1, y = MSW2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of MSW1 for each pedigree", 
       y = "% Mass of MSW2 for each pedigree", 
       title = "Warm Experiments Correlation of Percent Mass of Standards between Ga Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Warm Experiments Mass of Standards between Ga Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Warm Experiments Mass of Standards between Ga Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)




########################################################################## COLD
##### By Maize Line ###########################################################
coldData <- ogData[grepl("multiCold", ogData$experimentType),] # Filtered out experiments that just used Guu

##### Guu v Ga #####
#################### DISEASE
df_Guu_Ga <- coldData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree",
       y = "% Disease of Ga for each pedigree",
       title = "Cold Experiments Correlation of % Disease between Guu and Ga across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Disease between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Ga <- coldData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Ga for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Guu and Ga across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Length of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


#################### MASS
df_Guu_Ga <- coldData %>%
  filter(isolate %in% c("Guu", "Ga")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Ga)

cor_test <- cor.test(df_Guu_Ga$Guu, df_Guu_Ga$Ga)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Ga, aes(x = Guu, y = Ga, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Ga for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Guu and Ga across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Mass of Standards between Guu and Ga.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Guu and Ga.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


##### Guu v Gi #####
#################### DISEASE
df_Guu_Gi <- coldData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Guu for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Cold Experiments Correlation of % Disease between Guu and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Disease between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Gi <- coldData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Guu for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Guu and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Length of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Guu_Gi <- coldData %>%
  filter(isolate %in% c("Guu", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Gi)

cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Guu for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Guu and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Mass of Standards between Guu and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Guu and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

##### Ga v Gi #####
#################### DISEASE
df_Ga_Gi <- coldData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of Ga for each pedigree", 
       y = "% Disease of Gi for each pedigree", 
       title = "Cold Experiments Correlation of % Disease between Ga and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Disease between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Ga_Gi <- coldData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of Ga for each pedigree", 
       y = "% Length of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Ga and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Length of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Ga_Gi <- coldData %>%
  filter(isolate %in% c("Ga", "Gi")) %>%
  group_by(pedigree, experiment, isolate) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = isolate, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Ga_Gi)

cor_test <- cor.test(df_Ga_Gi$Ga, df_Ga_Gi$Gi)
r_value <- cor_test$estimate
p_value <- cor_test$p.value

ggplot(df_Ga_Gi, aes(x = Ga, y = Gi, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9, aes(shape = factor(experiment))) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of Ga for each pedigree", 
       y = "% Mass of Gi for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Ga and Gi across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Mass of Standards between Ga and Gi.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Ga and Gi.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

##### By Repetition ############################################################
##### Guu #####
#################### DISEASE
df_Guu_Comp <- coldData %>%
  filter(isolate == "Guu", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = lesionPercent) %>%
  drop_na()
print(df_Guu_Comp)

cor_test_Guu <- cor.test(df_Guu_Comp$MSC1, df_Guu_Comp$MSC2)
r_value <- cor_test_Guu$estimate
p_value <- cor_test_Guu$p.value

ggplot(df_Guu_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) +
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of MSC1 for each pedigree",
       y = "% Disease of MSC2 for each pedigree",
       title = "Cold Experiments Correlation of % Disease between Guu Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Disease between Guu Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Guu Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Guu_Comp <- coldData %>%
  filter(isolate == "Guu", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Guu_Comp)

cor_test_Guu <- cor.test(df_Guu_Comp$MSC1, df_Guu_Comp$MSC2)
r_value <- cor_test_Guu$estimate
p_value <- cor_test_Guu$p.value

ggplot(df_Guu_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of MSC1 for each pedigree", 
       y = "% Length of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Guu Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Length of Standards between Guu Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Guu Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


#################### MASS
df_Guu_Comp <- coldData %>%
  filter(isolate == "Guu", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Guu_Comp)

cor_test_Guu <- cor.test(df_Guu_Comp$MSC1, df_Guu_Comp$MSC2)
r_value <- cor_test_Guu$estimate
p_value <- cor_test_Guu$p.value

ggplot(df_Guu_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of MSC1 for each pedigree", 
       y = "% Mass of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Guu Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Mass of Standards between Guu Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Guu Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)


#####  Gi #####
#################### DISEASE
df_Gi_Comp <- coldData %>%
  filter(isolate == "Gi", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = lesionPercent) %>%
  drop_na()
print(df_Gi_Comp)

cor_test_Gi <- cor.test(df_Gi_Comp$MSC1, df_Gi_Comp$MSC2)
r_value <- cor_test_Gi$estimate
p_value <- cor_test_Gi$p.value

ggplot(df_Gi_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of MSC1 for each pedigree", 
       y = "% Disease of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of % Disease between Gi Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Disease between Gi Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Gi Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Gi_Comp <- coldData %>%
  filter(isolate == "Gi", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Gi_Comp)

cor_test_Gi <- cor.test(df_Gi_Comp$MSC1, df_Gi_Comp$MSC2)
r_value <- cor_test_Gi$estimate
p_value <- cor_test_Gi$p.value

ggplot(df_Gi_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of MSC1 for each pedigree", 
       y = "% Length of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Gi Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Length of Standards between Gi Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Gi Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Gi_Comp <- coldData %>%
  filter(isolate == "Gi", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Gi_Comp)

cor_test_Gi <- cor.test(df_Gi_Comp$MSC1, df_Gi_Comp$MSC2)
r_value <- cor_test_Gi$estimate
p_value <- cor_test_Gi$p.value

ggplot(df_Gi_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of MSC1 for each pedigree", 
       y = "% Mass of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Gi Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Mass of Standards between Gi Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Gi Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

##### Ga #####
#################### DISEASE
df_Ga_Comp <- coldData %>%
  filter(isolate == "Ga", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = lesionPercent) %>%
  drop_na()
print(df_Ga_Comp)

cor_test_Ga <- cor.test(df_Ga_Comp$MSC1, df_Ga_Comp$MSC2)
r_value <- cor_test_Ga$estimate
p_value <- cor_test_Ga$p.value

ggplot(df_Ga_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Disease of MSC1 for each pedigree", 
       y = "% Disease of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of % Disease between Ga Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Disease between Ga Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Disease between Ga Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### LENGTH
df_Ga_Comp <- coldData %>%
  filter(isolate == "Ga", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentLengthPedByIsolate = mean(percentLengthPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentLengthPedByIsolate) %>%
  drop_na()
print(df_Ga_Comp)

cor_test_Ga <- cor.test(df_Ga_Comp$MSC1, df_Ga_Comp$MSC2)
r_value <- cor_test_Ga$estimate
p_value <- cor_test_Ga$p.value

ggplot(df_Ga_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Length of MSC1 for each pedigree", 
       y = "% Length of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of Percent Length of Standards between Ga Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Length of Standards between Ga Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Length of Standards between Ga Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)

#################### MASS
df_Ga_Comp <- coldData %>%
  filter(isolate == "Ga", experiment %in% c("MSC1", "MSC2")) %>%
  group_by(pedigree, experiment) %>%
  dplyr::summarise(percentMassPedByIsolate = mean(percentMassPedByIsolate, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = experiment, values_from = percentMassPedByIsolate) %>%
  drop_na()
print(df_Ga_Comp)

cor_test_Ga <- cor.test(df_Ga_Comp$MSC1, df_Ga_Comp$MSC2)
r_value <- cor_test_Ga$estimate
p_value <- cor_test_Ga$p.value

ggplot(df_Ga_Comp, aes(x = MSC1, y = MSC2, color = pedigree)) +
  geom_point(size = 8, alpha = 0.9 ) + 
  annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
  labs(x = "% Mass of MSC1 for each pedigree", 
       y = "% Mass of MSC2 for each pedigree", 
       title = "Cold Experiments Correlation of Percent Mass of Standards between Ga Repetitions across Maize Lines") +
  theme_minimal()

filename <- paste0("Cold Experiments Mass of Standards between Ga Repetitions.pdf", sep = "") 
ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
filename2 <- paste0("Cold Experiments Mass of Standards between Ga Repetitions.png", sep = "") 
ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)



########################################################################################################################################################################################
############################################################## CORRELATIONS OVERALL ##################################################################################################
########################################################################################################################################################################################

corData <- ogData[!grepl("singleSpecies", ogData$experimentType),] # Filtered out experiments that just used Guu

  ##### Guu v Ga #####
  #################### DISEASE
  df_Guu_Ga <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) +
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value),
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Disease of Guu for each pedigree",
         y = "% Disease of Ga for each pedigree",
         title = "Correlation of % Disease between Guu and Ga across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Disease between Guu and Ga.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Disease between Guu and Ga.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  #################### LENGTH
  df_Guu_Ga <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Length of Guu for each pedigree", 
         y = "% Length of Ga for each pedigree", 
         title = "Correlation of Percent Length of Standards between Guu and Ga across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Length of Standards between Guu and Ga.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Length of Standards between Guu and Ga.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  
  #################### MASS
  df_Guu_Ga <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Mass of Guu for each pedigree", 
         y = "% Mass of Ga for each pedigree", 
         title = "Correlation of Percent Mass of Standards between Guu and Ga across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Mass of Standards between Guu and Ga.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Mass of Standards between Guu and Ga.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  
  ##### Guu v Gi #####
  #################### DISEASE
  df_Guu_Gi <- corData %>%
    filter(isolate %in% c("Guu", "Gi")) %>%
    group_by(pedigree, isolate) %>%
    dplyr::summarise(lesionPercent = mean(lesionPercent, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = isolate, values_from = lesionPercent) %>%
    drop_na()
  print(df_Guu_Gi)
  
  cor_test <- cor.test(df_Guu_Gi$Guu, df_Guu_Gi$Gi)
  r_value <- cor_test$estimate
  p_value <- cor_test$p.value
  
  ggplot(df_Guu_Gi, aes(x = Guu, y = Gi, color = pedigree)) +
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Disease of Guu for each pedigree", 
         y = "% Disease of Gi for each pedigree", 
         title = "Correlation of % Disease between Guu and Gi across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Disease between Guu and Gi.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Disease between Guu and Gi.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  #################### LENGTH
  df_Guu_Gi <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Length of Guu for each pedigree", 
         y = "% Length of Gi for each pedigree", 
         title = "Correlation of Percent Length of Standards between Guu and Gi across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Length of Standards between Guu and Gi.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Length of Standards between Guu and Gi.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  #################### MASS
  df_Guu_Gi <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Mass of Guu for each pedigree", 
         y = "% Mass of Gi for each pedigree", 
         title = "Correlation of Percent Mass of Standards between Guu and Gi across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Mass of Standards between Guu and Gi.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Mass of Standards between Guu and Gi.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  ##### Ga v Gi #####
  #################### DISEASE
  df_Ga_Gi <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Disease of Ga for each pedigree", 
         y = "% Disease of Gi for each pedigree", 
         title = "Correlation of % Disease between Ga and Gi across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Disease between Ga and Gi.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Disease between Ga and Gi.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  #################### LENGTH
  df_Ga_Gi <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Length of Ga for each pedigree", 
         y = "% Length of Gi for each pedigree", 
         title = "Correlation of Percent Length of Standards between Ga and Gi across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Length of Standards between Ga and Gi.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Length of Standards between Ga and Gi.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  
  #################### MASS
  df_Ga_Gi <- corData %>%
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
    geom_point(size = 8, alpha = 0.9) + 
    annotate("text", x = Inf, y = Inf, label = sprintf("R = %.2f, p = %.3f", r_value, p_value), 
             hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
    geom_smooth(method = "lm", color = "blue") +  # Adds a linear regression line
    labs(x = "% Mass of Ga for each pedigree", 
         y = "% Mass of Gi for each pedigree", 
         title = "Correlation of Percent Mass of Standards between Ga and Gi across Maize Lines") +
    theme_minimal()
  
  filename <- paste0("Mass of Standards between Ga and Gi.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Mass of Standards between Ga and Gi.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  

  ########################################################################################################################################################################################
  ######################################################################## EXPERIMENT CORRELATIONS ##################################################################################################
  ########################################################################################################################################################################################
  
  
  wide_data <- ogData %>%
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
    labs(title = "Correlation and P-Values between Experiments", x = "", y = "")
  
  filename <- paste0("Experiment Correlation and P-Values.pdf", sep = "") 
  ggsave(path = graphDir, filename, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  filename2 <- paste0("Experiment Correlation and P-Values.png", sep = "") 
  ggsave(path = graphDir, filename2, width = ChartWidth*2, height = ChartHeight*2, units = "px", dpi = 900, bg = "#FFFEFE", limitsize = FALSE)
  