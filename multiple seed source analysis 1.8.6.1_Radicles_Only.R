## Version 1.8.5 ###############################################################
##                                                                            ##
## New Features in 1.8.5:                                                     ##
##                                                                            ##
## + Output and save all analysis                                             ##
## + Removes all escapes                                                      ##
## + Analysis of controls and blocks                                          ##
##                                                                            ##
## TO DO:                                                                     ##
## - Rename variables so it is clear what they are for                        ##
## - Improve documentation of each section                                    ##
## - Add a glossary of variable names                                         ##
## - Add Correlation analysis                                                 ##
## - Analysis of blocks/controls is Exploratory Analysis. Add another layer   ##
## with non-standardized data. Do this every time, don't need to report on it.##
##                                                                            ##
##                                                                            ##
## #############################################################################
##                                                                            ##
################################################################################
##                           TABLE OF CONTENTS                                ##
## sec0.0: Future suggestions                                                 ##
## sec1.0: Libraries                                                          ##
## sec2.0: File and data setup                                                ##
## sec3.0: Experiment Analysis                                                ##
## sec4.0: Combined Experiment Analysis                                       ##
## sec5.0: Troubleshooting Tips & Thanks                                      ##
################################################################################

######### 0.0 Future Suggestions ###############################################
## Future Analysis Suggestions by Sarah:
##
# model.aug <- lm(ger.avg + 1 ~ plot + check + env/DTS + env + env/block + plot*env, weights = counts, data =  ge.filter)
# 
# aov(model.aug)
# 
# my.anova <- aov(model.aug)
# my.anova$Pr>F
#
# shapiro.test(model.aug$residuals) (edited) 
# 
# histogram(model.aug$residuals, col=rainbow(8))
## - Suggestions from Sarah Lipps:                                            ##
##                                                                            ##
## SARAH: ADD T-Test if only 2 variables of something                         ##
## Consider Weighting Per Plate for Germs/Non-Germs, will need to add function##
## to calculate and add these counts.                                         ##



####################################sec1.0######################################
##                                                                            ##
##                                  Libraries                                 ##
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
fillyesAnthocyanin <- "#CC79A7"
fillEight <- "#E8E7A2"
fillNine <- "#86DAB2"
fillnoAnthocyanin <- "#555555"
# All of the above colors work for the top three most common colorblind types
theme_set(theme_gray(base_family = "sans"))
################################################################################
##### Beginning of analysis ####################################################
################################################################################

## First let's look at our data:

head(data) # Take a peek at the data and make sure we grabbed the right file
str(data) # This shows us our column names and examples of data. We will want to make sure things are formatted the way we want.



germdata <- data[!grepl("n",data$germinated),]
pydata <- germdata[!grepl("y",germdata$escape),]
pydata$rootLength <- as.numeric(pydata$rootLength)
pydata$coleoptileLength <- as.numeric(pydata$coleoptileLength)
pydata_backup <- paste0(tablesDir,"/",dataName,"-Cleaned-Backup.csv",sep="")
write.csv(pydata,pydata_backup,row.names=FALSE)

combinedData <- data.frame() #Creates a blank data frame to cumulatively build onto


####################################sec3.0######################################
##                              Experiment analysis                           ##
##                              Manpulation of data                           ##
##                        Graph and statistical analysis                      ##
################################################################################

# File Loaded
experimentlist <- c(unique(pydata$assay))
experimentlist

# Begin analysis loop
for (n in experimentlist) {
  
  loopData <- pydata #Insulate the main dataframe from what the loop will use.
  #TEST
  #n = "assay5"
  #ENDTEST
  print(paste("Starting on ",n,sep=""))
  Sys.sleep(.1)
  ### Build Graph Info: ----
  assayNum <- n
  assayCheck <- which(loopData == assayNum,arr.ind=TRUE)  ##Check the assay, grabs first row details
  rowCheck <- paste(assayCheck[1,1]) ## Checks first row of this assay, should be control
  assay <- paste(loopData[rowCheck,2])
  isolate <- paste(loopData[rowCheck,3])
  inocMethod <- paste(loopData[rowCheck,4])
  isolateSpecies <- paste(loopData[rowCheck,5])
  agarType <- paste(loopData[rowCheck,6])
  DAI <- paste(loopData[rowCheck,13])
  escape <- paste(loopData[rowCheck,14])
  totalLength <- paste(loopData[rowCheck,15])
  seedSource <- paste(loopData[rowCheck,16])
  maizeLine <- paste(loopData[rowCheck,17])
  
  loopData
  
  subTitle = ""
  #subTitle <- paste("Assay: ",assay," Species: ,",isolateSpecies," Isolate: ",isolate," Inoculation Method: ",inocMethod," Agar Type: ",agarType," Days after inoculation: ",DAI,sep = "")   ## Build subtitle for chart
  cTitle <- paste("",sep="")
  
  
  loopData <- loopData[grepl(n,loopData$assay),] ### Begin building data frame. This removes anything that isn't in the assay
  loopData
  
  ## 3.2 Get length of standards: ----
  RootAvg <- aggregate(rootLength ~ seedType + isolate, loopData, mean)
  ctrlRoots <- RootAvg[grepl("control",RootAvg$isolate),]
  ctrlRoots$tissue <- 'radicle'
  ctrlRoots_Average_backup <- paste0(tablesDir,"/","Ctrl-Root-Avgs-",n,"-Cleaned-Backup.csv",sep="")
  write.csv(ctrlRoots,ctrlRoots_Average_backup,row.names=FALSE)
  #ctrlRootValues <- ctrlRoots[,3] ## Not useful in this version
  #ctrlRootValues
  ## Get Coleoptile length Control:
  LeafAvg <- aggregate(coleoptileLength ~ seedType + isolate, loopData, mean)
  ctrlLeaves <- LeafAvg[grepl("control",LeafAvg$isolate),]
  ctrlLeaves$tissue <- 'coleoptile'
  ctrlLeaves_Average_backup <- paste0(tablesDir,"/","Ctrl-Leaf-Avgs-",n,"-Cleaned-Backup.csv",sep="")
  write.csv(ctrlLeaves,ctrlLeaves_Average_backup,row.names=FALSE)
  #ctrlLAValue <- ControlLeafAverage[,4]
  totalAvg <- aggregate(totalLength ~ seedType + isolate, loopData, mean)
  totalAvg
  ctrlTotal <- totalAvg[grepl("control",totalAvg$isolate),]
  ctrlTotal
  ctrlTotal$tissue <- 'total'
  ctrlTotal_Average_backup <- paste0(tablesDir,"/","Ctrl-Leaf-Avgs-",n,"-Cleaned-Backup.csv",sep="")
  write.csv(ctrlTotal,ctrlTotal_Average_backup,row.names=FALSE)
  
  ############################ Sort by pigmentation ##########################
  RootPigmentAvg <- aggregate(rootLength ~ seedType + isolate + redPigment, loopData, mean)
  PigmentRoots <- RootAvg[grepl("control",RootPigmentAvg$isolate),]
  PigmentRoots$tissue <- 'radicle'
  PigmentRoots_Average_backup <- paste0(tablesDir,"/","Pigment-Root-Avgs-",n,"-Cleaned-Backup.csv",sep="")
  write.csv(PigmentRoots,PigmentRoots_Average_backup,row.names=FALSE)
  
  LeafPigmentAvg <- aggregate(coleoptileLength ~ seedType + isolate, loopData, mean)
  PigmentLeaves <- LeafAvg[grepl("control",LeafPigmentAvg$isolate),]
  PigmentLeaves$tissue <- 'coleoptile'
  PigmentLeaves_Average_backup <- paste0(tablesDir,"/","Pigment-Leaf-Avgs-",n,"-Cleaned-Backup.csv",sep="")
  write.csv(PigmentLeaves,PigmentLeaves_Average_backup,row.names=FALSE)
  ############################################################################
  speakCtrlsAveraged <- paste0("Controls averaged and backed up!")
  print(speakCtrlsAveraged)
  Sys.sleep(.25)
  
  
  
  ## 3.3 Get lengths of maize lines ----
  seedList <- c(unique(loopData$seedType))
  
  unModifiedLengthsAll <- data.frame()
  standardizedLengthsAll <- data.frame()
  
  ## 3.4 Standardize Dataset ----
  for (s in seedList){
    s
    standardizedLengths <- data.frame()
    SLoopDataNew <- data.frame(seedType=loopData$seedType, assay=loopData$assay, isolate=loopData$isolate, inocMethod=loopData$inocMethod, isolateSpecies=loopData$isolateSpecies, agarType=loopData$agarType, seedNum=loopData$seedNum, block=loopData$block,rootLength=loopData$rootLength,coleoptileLength=loopData$coleoptileLength,redPigment=loopData$redPigment,DAI=loopData$DAI,totalLength=loopData$totalLength, seedSource=loopData$packet_ID, maizeLine=loopData$maize_line)
    SLoopDataNew
    SLoopData <- SLoopDataNew[grepl(s,SLoopDataNew$seedType),]
    SLoopData
    sRootAvg <- ctrlRoots[grepl(s,ctrlRoots$seedType),]
    sRootAvg
    sRootAvgVal <- sRootAvg[1,3]
    sRootAvgVal
    sLeafAvg <- ctrlLeaves[grepl(s,ctrlLeaves$seedType),]
    sLeafAvg
    sLeafAvgVal <- sLeafAvg[1,3]
    sLeafAvgVal
    sTotalAvg <- ctrlTotal[grepl(s,ctrlTotal$seedType),]
    sTotalAvg
    sTotalAvgVal <- sTotalAvg[1,3]
    sTotalAvgVal
    ## NOT standardized
    sRootsNormal <- data.frame(seedType=SLoopData[,1],assay=SLoopData[,2],isolate=SLoopData[,3],inocMethod=SLoopData[,4],isolateSpecies=SLoopData[,5],agarType=SLoopData[,6],seedNum=SLoopData[,7],block=SLoopData[,8],length=SLoopData[,9],redPigment=SLoopData[,11],DAI=SLoopData[,12],seedSource=SLoopData[,14],maizeLine=SLoopData[,15])
    sRootsNormal$tissue <- 'radicle'
    sRootsNormal
    sLeafsNormal <- data.frame(seedType=SLoopData[,1],assay=SLoopData[,2],isolate=SLoopData[,3],inocMethod=SLoopData[,4],isolateSpecies=SLoopData[,5],agarType=SLoopData[,6],seedNum=SLoopData[,7],block=SLoopData[,8],length=SLoopData[,10],redPigment=SLoopData[,11],DAI=SLoopData[,12],seedSource=SLoopData[,14],maizeLine=SLoopData[,15])
    sLeafsNormal$tissue <- 'coleoptile'
    sTotalNormal <- data.frame(seedType=SLoopData[,1],assay=SLoopData[,2],isolate=SLoopData[,3],inocMethod=SLoopData[,4],isolateSpecies=SLoopData[,5],agarType=SLoopData[,6],seedNum=SLoopData[,7],block=SLoopData[,8],length=SLoopData[,13],redPigment=SLoopData[,11],DAI=SLoopData[,12],seedSource=SLoopData[,14],maizeLine=SLoopData[,15])
    sTotalNormal$tissue <- 'total'
    
    ## standardized:
    sRootsAllAveraged <- data.frame(seedType=SLoopData[,1],assay=SLoopData[,2],isolate=SLoopData[,3],inocMethod=SLoopData[,4],isolateSpecies=SLoopData[,5],agarType=SLoopData[,6],seedNum=SLoopData[,7],block=SLoopData[,8],length=SLoopData[,9]/sRootAvgVal*100,redPigment=SLoopData[,11],DAI=SLoopData[,12],seedSource=SLoopData[,14],maizeLine=SLoopData[,15])
    sRootsAllAveraged$tissue <- 'radicle'
    sLeafsAllAveraged <- data.frame(seedType=SLoopData[,1],assay=SLoopData[,2],isolate=SLoopData[,3],inocMethod=SLoopData[,4],isolateSpecies=SLoopData[,5],agarType=SLoopData[,6],seedNum=SLoopData[,7],block=SLoopData[,8],length=SLoopData[,10]/sLeafAvgVal*100,redPigment=SLoopData[,11],DAI=SLoopData[,12],seedSource=SLoopData[,14],maizeLine=SLoopData[,15])
    sLeafsAllAveraged$tissue <- 'coleoptile'
    sTotalAllAveraged <- data.frame(seedType=SLoopData[,1],assay=SLoopData[,2],isolate=SLoopData[,3],inocMethod=SLoopData[,4],isolateSpecies=SLoopData[,5],agarType=SLoopData[,6],seedNum=SLoopData[,7],block=SLoopData[,8],length=SLoopData[,13]/sTotalAvgVal*100,redPigment=SLoopData[,11],DAI=SLoopData[,12],seedSource=SLoopData[,14],maizeLine=SLoopData[,15])
    sTotalAllAveraged$tissue <- 'total'
    
    
    unModifiedLengths <- rbind (sRootsNormal,sLeafsNormal,sTotalNormal)
    
    ##### Combine Data #####
    
    standardizedLengths <- rbind(sRootsAllAveraged,sLeafsAllAveraged,sTotalAllAveraged)
    
    ##### PIGMENT CHECK #####
    
    sTreatedPigments <- standardizedLengths[!grepl("control",standardizedLengths$isolate),]
    
    ## If some seeds on a maize line show anthocyanin, assume all of them in that maize line will
    if ((sum(sTreatedPigments$redPigment == 'y'))>1){
      standardizedLengths$redPigment <- 'y'
    } else {
      standardizedLengths$redPigment <- 'n'
    }
    
    #### Add to total ####
    
    standardizedLengthsAll <-rbind(standardizedLengths,standardizedLengthsAll)
    unModifiedLengthsAll <- rbind(unModifiedLengths,unModifiedLengthsAll)
    print(paste("Values standardized to ",s," controls",sep=""))
    Sys.sleep(.25)
    
  }
  standardizedLengthsAll <- unique(standardizedLengthsAll)
  
  # Calculating mean and standard deviation per group
  stats <- standardizedLengthsAll %>%
    group_by(seedType, tissue, isolate) %>%
    mutate(
      meanLength = mean(length, na.rm = TRUE),
      sdLength = sd(length, na.rm = TRUE)
    ) %>%
    ungroup()
  
  filteredData <- stats %>%
    filter(between(length, meanLength - 2.01 * sdLength, meanLength + 2.01 * sdLength)) %>%
    dplyr::select(-meanLength, -sdLength)  # Remove calculation columns after filtering
  
  ## 3.4 Force proper data type ----
  LoopIIData <- filteredData
  LoopIIData
  LoopIIData$isolate <- factor(LoopIIData$isolate, levels = c(unique(LoopIIData$isolate)))
  lyheight <- max(LoopIIData$length)
  lymin <- min(LoopIIData$length)
  
  LoopIIData$seedType <- as.factor(LoopIIData$seedType)
  LoopIIData$assay <- as.factor(LoopIIData$assay)
  LoopIIData$isolate <- as.factor(LoopIIData$isolate)
  LoopIIData$inocMethod <- as.factor(LoopIIData$inocMethod)
  LoopIIData$isolateSpecies <- as.factor(LoopIIData$isolateSpecies)
  LoopIIData$agarType <- as.factor(LoopIIData$agarType)
  LoopIIData$seedNum <- as.factor(LoopIIData$seedNum)
  LoopIIData$block <- as.factor(LoopIIData$block)
  LoopIIData$length <- as.numeric(LoopIIData$length)
  LoopIIData$redPigment <- as.factor(LoopIIData$redPigment)
  LoopIIData$DAI <- as.factor(LoopIIData$DAI)
  LoopIIData$tissue <- as.factor(LoopIIData$tissue)
  LoopIIData$maizeLine <- as.factor(LoopIIData$maizeLine)
  
  ## 3.5 Begin Analysis ----
  
  ## 3.5.1 Radicle length analysis ----
  tukeyRadiclesOnly <- LoopIIData[grepl("radicle",LoopIIData$tissue),]
  tukeyRadiclesOnly <- tukeyRadiclesOnly[!grepl("control",tukeyRadiclesOnly$isolate),]
  tukeyRadiclesOnly$seedType = with(tukeyRadiclesOnly, reorder(seedType, -length, mean))
  lyheight <- max(tukeyRadiclesOnly$length)
  lymin <- min(tukeyRadiclesOnly$length)
  
  RadAnova <- aov(length~seedType+block,tukeyRadiclesOnly)
  summary(RadAnova)
  RadTukey <- TukeyHSD(RadAnova)
  cld <- multcompLetters4(RadAnova,RadTukey)
  RTk <- group_by(tukeyRadiclesOnly,seedType) %>%
      dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) #%>%
      #arrange(desc(sd))
  cld <- as.data.frame.list(cld$seedType)
  RTk$cld <- cld$Letters
  RTk$redPigment <- "y"
  print(RTk)
  
  #esquisser(tukeyRadiclesOnly)
  
  ## 3.5.2 Radicle length graph ----
  
  ylab = "Radicle lengths standardized to control"
  tukeyRadiclesOnly %>%
    filter(tissue %in% "radicle") %>%
    filter(!(isolate %in% "control")) %>%
    ggplot(outlier.shape = NA) +
    geom_hline(yintercept=100, color = "#999999") +
    aes(x = seedType, y = length, fill = maizeLine, drop=TRUE, middle = mean(length)) +
    ylim(lymin,lyheight+50) +
    #geom_segment(aes(x=0,xend=11,y=100,yend=100))+
    #scale_fill_manual(values=c(fillnoAnthocyanin,fillyesAnthocyanin))+
    geom_boxplot(position = "dodge", width = .9, aes(x = seedType, y = length, fill = maizeLine, middle = mean(length))
    )+
    geom_point(shape = 25, color = "black", fill = NA, size = 0.5, alpha = 1, position = position_jitter(width = 0.22)) +
    scale_fill_brewer(palette = "Set2", direction = 1) +
    geom_text(data = RTk, aes(x=seedType, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
    labs(x = "", y = ylab, title = "", 
         subtitle = n, 
         fill = "Red Pigment", color = "") + 
    stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
    theme_minimal() +
    theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
          panel.grid.major.y=element_line(colour = "#f0f0f0"),
          panel.grid.major.x=element_line(colour = "#cecece"),
          #axis.ticks=element_line(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.subtitle = element_text(hjust=0, size = fontSize/2),
          #legend.position = c(0.5, 1.03),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size=fontSize),
          text = element_text(size = fontSize),
          axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  
  
  Sys.sleep(.1)
  filename <- paste0(assayNum,"-radicles-standardized",".jpg",sep="") 
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0(assayNum,"-radicles-standardized",".wmf",sep="") 
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  csvFilename <- paste0(tablesDir,"/",assayNum,"-radicles-standardized",".csv",sep="")
  write.csv(LoopIIData,csvFilename,row.names=FALSE)
  rm(cld)
  
  stage1Announce <- paste0(filename," and corresponding csv finished and saved.",sep="")
  print(stage1Announce)
  Sys.sleep(.25)
  
  ##### Experimental Analysis
  
  #rm(cld)
  #rm(ATk)
  #rm(BTk)
  #rm(CATk)
  #rm(CBTk)
  #rm(CCTk)
  #rm(CRTk)
  #rm(CTk)
  #rm(RTk)
  #rm(TTk)
  
  ## 3.5.11 Radicle length analysis by maizeline ----
  MTR <- LoopIIData[grepl("radicle",LoopIIData$tissue),]
  
  MtukeyRadiclesOnly <- MTR[!grepl("control",MTR$isolate),]
  MtukeyRadiclesOnly$maizeLine = with(MtukeyRadiclesOnly, reorder(maizeLine, -length, mean))
  lyheight <- max(MtukeyRadiclesOnly$length)
  lymin <- min(MtukeyRadiclesOnly$length)
  
  RadAnova <- aov(length~maizeLine,MtukeyRadiclesOnly)
  RadTukey <- TukeyHSD(RadAnova)
  cld <- multcompLetters4(RadAnova,RadTukey)
  MTk <- group_by(MtukeyRadiclesOnly,maizeLine) %>%
    dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) #%>%
    #arrange(desc(sd))
  cld <- as.data.frame.list(cld$maizeLine)
  MTk
  #MTk <- MTk %>% filter(row_number() <= n()-1)
  MTk$cld <- cld$Letters
  MTk$redPigment <- "y"
  print(MTk)
  
  cld

  ## 3.5.21 Radicle length graph by maizeline----
  
  ylab = "Radicle lengths standardized to control tukey by maizeline"
  MtukeyRadiclesOnly %>%
    filter(tissue %in% "radicle") %>%
    filter(!(isolate %in% "control")) %>%
    ggplot(outlier.shape = NA) +
    geom_hline(yintercept=100, color = "#999999") +
    aes(x = maizeLine, y = length, fill = redPigment, drop=TRUE, middle = mean(length)) +
    ylim(lymin,lyheight+50) +
    #geom_segment(aes(x=0,xend=11,y=100,yend=100))+
    scale_fill_manual(values=c(fillnoAnthocyanin,fillyesAnthocyanin))+
    geom_boxplot(position = "dodge", width = .9, aes(x = maizeLine, y = length, fill = redPigment, middle = mean(length))
    )+
    geom_point(shape = 25, color = "black", fill = NA, size = 0.5, alpha = 1, position = position_jitter(width = 0.22)) +
    geom_text(data = MTk, aes(x=maizeLine, y= -Inf, label = cld), position = position_jitter(width=0), size = fontSize, vjust=-4) +
    labs(x = "", y = ylab, title = "", 
         subtitle = n, 
         fill = "Red Pigment", color = "") + 
    stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
    theme_minimal() +
    theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
          panel.grid.major.y=element_line(colour = "#f0f0f0"),
          panel.grid.major.x=element_line(colour = "#cecece"),
          #axis.ticks=element_line(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.subtitle = element_text(hjust=0, size = fontSize/2),
          #legend.position = c(0.5, 1.03),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size=fontSize),
          text = element_text(size = fontSize),
          axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  
  
  Sys.sleep(.1)
  filename <- paste0(assayNum,"-radicles-standardized-tukey-maizeline",".jpg",sep="") 
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0(assayNum,"-radicles-standardized-tukey-maizeline",".wmf",sep="") 
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  csvFilename <- paste0(tablesDir,"/",assayNum,"-radicles-standardized-tukey-maizeline",".csv",sep="")
  write.csv(LoopIIData,csvFilename,row.names=FALSE)
  rm(cld)
  
  stage1Announce <- paste0(filename," and corresponding csv finished and saved.",sep="")
  print(stage1Announce)
  Sys.sleep(.25)
  
  
  ## 3.6 Analysis of Variance ----  
  # stage3Announce <- paste0("Starting on statistical analysis...",sep="")
  # print(stage3Announce)
  # Sys.sleep(.25)
  # 
  # # Factor forcing:
  # ##### 
  # LoopIIData$seedType <- as.factor(LoopIIData$seedType)
  # LoopIIData$assay <- as.factor(LoopIIData$assay)
  # LoopIIData$isolate <- as.factor(LoopIIData$isolate)
  # LoopIIData$inocMethod <- as.factor(LoopIIData$inocMethod)
  # LoopIIData$isolateSpecies <- as.factor(LoopIIData$isolateSpecies)
  # LoopIIData$agarType <- as.factor(LoopIIData$agarType)
  # LoopIIData$seedNum <- as.factor(LoopIIData$seedNum)
  # LoopIIData$block <- as.factor(LoopIIData$block)
  # LoopIIData$length <- as.numeric(LoopIIData$length)
  # LoopIIData$redPigment <- as.factor(LoopIIData$redPigment)
  # LoopIIData$DAI <- as.factor(LoopIIData$DAI)
  # LoopIIData$tissue <- as.factor(LoopIIData$tissue)
  # # For some reason the data types have to be forced again.
  # # Not sure why this happens, but if I don't, R complains.
  combinedData <- rbind(combinedData, LoopIIData)
  
  
  
  
  # #####
  # # Filters:
  # #####
  # 
  # LoopIIData
  # 
  # # FILTER FALSE FLAGS
  # ControlData <- LoopIIData[grepl("control",LoopIIData$isolate),]
  # LoopIIData <- LoopIIData[!grepl("total",LoopIIData$tissue),]
  # LoopIIData <- LoopIIData[!grepl("control",LoopIIData$isolate),]
  # # If these are left in, I will get false significance
  # #####
  # ## 3.6.1 Exploratory Analysis of Blocks ----
  # 
  # NoRadLoop <- LoopIIData
  # NoRadLoop <- NoRadLoop[!grepl("radicle",NoRadLoop$tissue),]
  # NoColeLoop <- LoopIIData
  # NoColeLoop <- NoColeLoop[!grepl("coleoptile",NoColeLoop$tissue),]
  
  ## 3.6.2 Radicle graph ----
  
  # BlockAnova <- aov(length~block,NoColeLoop)
  # BlockTukey <- TukeyHSD(BlockAnova)
  # cld <- multcompLetters4(BlockAnova,BlockTukey)
  # BTk <- group_by(NoColeLoop,block) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # cld <- as.data.frame.list(cld$block)
  # BTk$cld <- cld$Letters
  # BTk$redPigment <- "y"
  # BTk$tissue <- "radicle"
  # 
  # NoColeLoop %>%
  #   filter(tissue %in% "radicle") %>%
  #   ggplot(outlier.shape = NA) +
  #   aes(x = block, y = length, fill = block) +
  #   geom_boxplot() +
  #   scale_fill_hue(direction = 1) +
  #   labs(y = "Radicle lengths by block") +
  #   geom_text(data = BTk, aes(x=, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = c(0.5, 1.03),
  #         legend.box = "horizontal",
  #         legend.direction = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0(assayNum,"-radicle-lengths-block",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0(assayNum,"-radicle-lengths-block",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # rm(cld)
  # ## 3.6.4 Radicle Standards ----
  # NoRadLoopC <- ControlData
  # NoRadLoopC <- NoRadLoopC[!grepl("radicle",NoRadLoopC$tissue),]
  # NoColeLoopC <- ControlData
  # NoColeLoopC <- NoColeLoopC[!grepl("coleoptile",NoColeLoopC$tissue),]
  # 
  # ControlData
  # NoColeLoopC
  # 
  # CtrlBlockAnova <- aov(length~block,NoColeLoopC)
  # CtrlBlockTukey <- TukeyHSD(CtrlBlockAnova)
  # cld <- multcompLetters4(CtrlBlockAnova,CtrlBlockTukey)
  # cld
  # CtrlBTK <- group_by(NoColeLoopC,block) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # cld <- as.data.frame.list(cld$block)
  # CtrlBTK$cld <- cld$Letters
  # CtrlBTK$redPigment <- "y"
  # CtrlBTK$tissue <- "radicle"
  # 
  # NoColeLoopC %>%
  #   filter(tissue %in% "radicle") %>%
  #   ggplot(outlier.shape = NA) +
  #   aes(x = block, y = length, fill = block) +
  #   geom_boxplot() +
  #   scale_fill_hue(direction = 1) +
  #   labs(y = "Radicle lengths by Control block") +
  #   geom_text(data = CtrlBTK, aes(x=, y= mean, label = cld), position = position_jitter(height=5, width=0), size = fontSize/3, vjust=-4) +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = c(0.5, 1.03),
  #         legend.box = "horizontal",
  #         legend.direction = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0(assayNum,"-radicle-lengths-Controlblock",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0(assayNum,"-radicle-lengths-Controlblock",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # rm(cld)
  # 
  # ## 3.8 Analysis of Variance ----
  # # Length Significance w/Blocks ----
  # anovaDoc <- paste0(anovaDir,"/",assayNum,"block+seedType-RadiclesOnly.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType, NoColeLoop)),file=anovaDoc) 
  # 
  # # Pairwise backup for manual double check of data analysis ----
  # 
  # anovaDoc <- paste0(anovaDir,"/",n,"block+seedType-RadiclesOnly.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType, NoColeCom)),file=anovaDoc) 
  # 
  # summary(aov(length ~ maizeLine + seedType %in% maizeLine, LoopIIData))
  # summary(aov(length ~ seedType, LoopIIData))
  # 
  # ## 3.9 Remove some tables that have a tendency to linger ----
  # rm(cld)
  # rm(ATk)
  # rm(BTk)
  # rm(CATk)
  # rm(CBTk)
  # rm(CCTk)
  # rm(CRTk)
  # rm(CTk)
  # rm(RTk)
  # rm(TTk)
  
}

####################################sec4.0############################################
##                                                                            ##
##                                                                            ##
##                         Combined Experiment Analysis                       ##
################################################################################

# 4.1 Loop announce
loopFinish2Announce <- paste0("Analysis of individual experiments finished, moving onto combined data...",sep="")
print(loopFinish2Announce)
Sys.sleep(.1)

####################################sec4.1 if loop##############################

if ((length(unique(experimentlist)))<2){
  NoMoreAnnounce1 <- paste0("...not enough experiments to do combined data",sep="")
  print(NoMoreAnnounce1)
  Sys.sleep(.1)
  
} else {     #Check Experiment Length
  
  # 4.2 Force Data Types ----

  combinedData$seedType <- as.factor(combinedData$seedType)
  combinedData$assay <- as.factor(combinedData$assay)
  combinedData$isolate <- as.factor(combinedData$isolate)
  combinedData$inocMethod <- as.factor(combinedData$inocMethod)
  combinedData$isolateSpecies <- as.factor(combinedData$isolateSpecies)
  combinedData$agarType <- as.factor(combinedData$agarType)
  combinedData$seedNum <- as.factor(combinedData$seedNum)
  combinedData$block <- as.factor(combinedData$block)
  combinedData$length <- as.numeric(combinedData$length)
  combinedData$redPigment <- as.factor(combinedData$redPigment)
  combinedData$DAI <- as.factor(combinedData$DAI)
  combinedData$tissue <- as.factor(combinedData$tissue)
  # 4.3 Combined Radicle Graph ----
  combinedRad <- combinedData[grepl("radicle",combinedData$tissue),]
  combinedRad <- combinedRad[!grepl("control",combinedRad$isolate),]
  combinedRad$seedType = with(combinedRad, reorder(seedType, -length, mean))
  lyheight <- max(combinedRad$length-70)
  lymin <- min(combinedRad$length)
  
  ## Significance:
  combinedRanova <- aov(length~seedType+assay+block,combinedRad)
  summary(combinedRanova)
  combinedRTukey <- TukeyHSD(combinedRanova)
  cld <- multcompLetters4(combinedRanova,combinedRTukey)
  CRTk <- group_by(combinedRad,seedType) %>%
    dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
    arrange(desc(sd))
  cld <- as.data.frame.list(cld$seedType)
  
  CRTk$cld <- cld[match(CRTk$seedType, rownames(cld)), "Letters"]
  CRTk$redPigment <- "y"
  
  
  
  Expanova <- aov(length~maizeLine*seedSource+assay+block,combinedRad)
  summary(Expanova)
  
  
  # 4.3.1 Graph ----
  
  ylab = "Radicle lengths as % of Standards"
  combinedRad %>%
    filter(tissue %in% "radicle") %>%
    filter(!(isolate %in% "control")) %>%
    ggplot(outlier.shape = NA) +     
    geom_hline(yintercept=100, color = "#999999") +
    aes(x = seedType, y = length, drop=TRUE, middle = mean(length)) +
    ylim(lymin-15,lyheight+10) +
    #geom_segment(aes(x=0,xend=11,y=100,yend=100))+
    geom_boxplot(outlier.shape = NA, position = "dodge", width = .9, aes(x = seedType, fill = maizeLine, y = length, middle = mean(length))
    )+
    geom_point(shape = 25, color = "black", fill = NA, size = .1, alpha = 1, position = position_jitter(width = 0.22)) +
    geom_text(data = CRTk, aes(x=seedType, y= -Inf, label = cld), angle=90, position = position_jitter(width=0), hjust=-0.05) +
    labs(x = "", y = ylab, title = "", 
         subtitle = " ", 
         fill = "Red Pigment: ", color = "") + 
    scale_fill_manual(values = c("H100" = "#353535", "Ms71" = "#BDBDBD", "IL14H" = "#555555", "CML 103" = "#E0E0E0", "B73" = "#9B9B9B", "CML 69" = "#F9F9F9")) + 
    #stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
    theme_minimal() +
    theme(
      text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
      legend.position = "none",
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank() 
    )
  
  Sys.sleep(.1)
  filename <- paste0("SeedSource-radicles-standardized",".png",sep="") 
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=600,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0("SeedSource-standardized",".pdf",sep="") 
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=600,bg="#FFFEFE",limitsize=FALSE)
  csvFilename <- paste0(tablesDir,"/All-radicles-standardized",".csv",sep="")
  write.csv(combinedData,csvFilename,row.names=FALSE)
  
  stage5Announce <- paste0(filename," and corresponding csv finished and saved.",sep="")
  print(stage5Announce)
  Sys.sleep(.25)
  
  
  summary(aov(length~maizeLine*seedType,combinedRad))
  
  ##### TukeyHSD BY MAIZE LINE:
  combinedRad
  
  maizeList <-c(unique(combinedRad$maizeLine))
  
  for (m in maizeList){
    m
    
    mLineOnly <- combinedRad[grepl("radicle",combinedRad$tissue),]
    mLineOnly <- mLineOnly[!grepl("control",mLineOnly$isolate),]
    mLineOnly <- mLineOnly[grepl(m,mLineOnly$maizeLine),]
    mLineOnly$seedType = with(mLineOnly, reorder(seedType, -length, mean))
    
    mRadAnova <- aov(length~seedType,mLineOnly)
    print(summary(mRadAnova))
    mRadTukey <- TukeyHSD(mRadAnova)
    print(mRadTukey)
    cld <- multcompLetters4(mRadAnova,mRadTukey)
    RTk <- group_by(mLineOnly,seedType,maizeLine) %>%
      dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) #%>%
    #arrange(desc(sd))
    
    
    cld <- as.data.frame.list(cld$seedType)
    RTk$cld <- cld$Letters
    RTk$redPigment <- "y"
    print(RTk)
    rootTukeyBackup <- paste0(tablesDir,"/",m," combined Tukey Backup.csv",sep="")
    write.csv(RTk,rootTukeyBackup,row.names=FALSE)
    ############################################################################
    print(paste0("Tukey analysis of roots for ",m," finished and backed up!"))
    Sys.sleep(.1)
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  ## 3.5.1 Root length analysis ----
  combinedTRootsOnly <- combinedRad[grepl("radicle",combinedRad$tissue),]
  combinedTRootsOnly <- combinedTRootsOnly[!grepl("control",combinedTRootsOnly$isolate),]
  combinedTRootsOnly$seedType = with(combinedTRootsOnly, reorder(seedType, -length, mean))
  lyheight <- max(combinedTRootsOnly$length)
  lymin <- min(combinedTRootsOnly$length)
  RadAnova <- aov(length~seedType,combinedTRootsOnly)
  RadTukey <- TukeyHSD(RadAnova)
  cld <- multcompLetters4(RadAnova,RadTukey)
  RTk <- group_by(combinedTRootsOnly,seedType,maizeLine) %>%
    dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) #%>%
    #arrange(desc(sd))
  
  
  cld <- as.data.frame.list(cld$seedType)
  RTk$cld <- cld$Letters
  RTk$redPigment <- "y"
  print(RTk)
  rootTukeyBackup <- paste0(tablesDir,"/","rootTukey Combined",n,"Backup.csv",sep="")
  write.csv(RTk,rootTukeyBackup,row.names=FALSE)
  ############################################################################
  print(paste0("Tukey analysis of roots for ",n," finished and backed up!"))
  Sys.sleep(.1)
  
  ## 3.5.2 Radicle length graph ----
  combinedTRootsOnly
  
  #combinedTRootsOnly$seedType = with(combinedTRootsOnly, reorder(maizeLine, -length, mean))
  
  ylab = "Radicle lengths standardized to standards"
  combinedTRootsOnly %>%
    filter(tissue %in% "radicle") %>%
    filter(!(isolate %in% "control")) %>%
    ggplot(outlier.shape = NA) +
    geom_hline(yintercept=100, color = "#999999") +
    aes(x = reorder(combinedTRootsOnly$seedType, -as.numeric(combinedTRootsOnly$maizeLine)), y = combinedTRootsOnly$length, drop=TRUE, middle = mean(combinedTRootsOnly$length)) +
    ylim(lymin,lyheight+50) +
    geom_boxplot(position = "dodge", width = .9, aes(x = reorder(combinedTRootsOnly$seedType, -as.numeric(combinedTRootsOnly$maizeLine)), y = combinedTRootsOnly$length, fill = combinedTRootsOnly$maizeLine, middle = mean(combinedTRootsOnly$length))
    )+
    geom_point(shape = 25, color = "black", fill = NA, size = 0.5, alpha = 1, position = position_jitter(width = 0.22)) +
    scale_fill_brewer(palette = "Set2", direction = 1) +
    geom_text(data = RTk, aes(x=seedType, y= max+5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
    labs(x = "", y = ylab, title = "", 
         subtitle = n, 
         fill = "", color = "") + 
    stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
    theme_minimal() +
    theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
          panel.grid.major.y=element_line(colour = "#f0f0f0"),
          panel.grid.major.x=element_line(colour = "#cecece"),
          #axis.ticks=element_line(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.subtitle = element_text(hjust=0, size = fontSize/2),
          #legend.position = c(0.5, 1.03),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size=fontSize),
          text = element_text(size = fontSize),
          axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  Sys.sleep(.05)
  filename <- paste0("combined-roots-standardized-maizeLine-groups",".jpg",sep="") 
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0("combined-roots-standardized-maizeLine-groups",".wmf",sep="") 
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  csvFilename <- paste0(tablesDir,"/","combined-roots-standardized",".csv",sep="")
  write.csv(combinedRad,csvFilename,row.names=FALSE)

  
  ylab = "Radicle lengths standardized to standards"
  combinedTRootsOnly %>%
    filter(tissue %in% "radicle") %>%
    filter(!(isolate %in% "control")) %>%
    dplyr::mutate(seedType = factor(seedType, levels = c("22NN_H100","18GV_H100","16ND_H100","22NN_IL14H","19NN_IL14H","16ND_IL14H","B73_21NN","22NN_B73","19NN_CML103","22NN_CML103","17NN_CML69","16NN_CML69","22NN_Ms71","19NN_Ms71","16ND_Ms71"))) %>%
    ggplot(outlier.shape = NA) +
    geom_hline(yintercept=100, color = "#999999") +
    aes(x = seedType, y = combinedTRootsOnly$length, drop=TRUE, middle = mean(combinedTRootsOnly$length)) +
    ylim(lymin,lyheight+50) +
    geom_boxplot(position = "dodge", width = .9, aes(x = seedType, y = combinedTRootsOnly$length, fill = combinedTRootsOnly$maizeLine, middle = mean(combinedTRootsOnly$length))
    )+
    geom_point(shape = 25, color = "black", fill = NA, size = 0.5, alpha = 1, position = position_jitter(width = 0.22)) +
    scale_fill_brewer(palette = "Set2", direction = 1) +
    geom_text(data = RTk, aes(x=seedType, y= max+5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
    labs(x = "", y = ylab, title = "", 
         subtitle = n, 
         fill = "", color = "") + 
    stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
    theme_minimal() +
    theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
          panel.grid.major.y=element_line(colour = "#f0f0f0"),
          panel.grid.major.x=element_line(colour = "#cecece"),
          #axis.ticks=element_line(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.subtitle = element_text(hjust=0, size = fontSize/2),
          #legend.position = c(0.5, 1.03),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size=fontSize),
          text = element_text(size = fontSize),
          axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  Sys.sleep(.05)
  filename <- paste0("combined-roots-manual-maizeLine-groups",".jpg",sep="") 
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0("combined-roots-manual-maizeLine-groups",".wmf",sep="") 
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)

  
  ##############################################################################
  # 4.4 Combined Coleoptile Graph ----
  # combinedCol <- combinedData[grepl("coleoptile",combinedData$tissue),]
  # combinedCol <- combinedCol[!grepl("control",combinedCol$isolate),]
  # combinedCol$seedType = with(combinedCol, reorder(seedType, -length, mean))
  # lyheight <- max(combinedCol$length)
  # lymin <- min(combinedCol$length)
  # # Significance:
  # combinedCanova <- aov(length~seedType,combinedCol)
  # combinedCtuke <- TukeyHSD(combinedCanova)
  # cld <- multcompLetters4(combinedCanova,combinedCtuke)
  # cld
  # CCTk <- group_by(combinedCol,seedType) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # cld <- as.data.frame.list(cld$seedType)
  # CCTk$cld <- cld$Letters
  # CCTk$redPigment <- "y"
  # # 4.4.1 Graph ----
  # ylab = "Coleoptile lengths standardized to control"
  # combinedCol %>%
  #   filter(tissue %in% "coleoptile") %>%
  #   filter(!(isolate %in% "control")) %>%
  #   ggplot(outlier.shape = NA) +     
  #   geom_hline(yintercept=100, color = "#999999") +
  #   aes(x = seedType, y = length, fill = redPigment, drop=TRUE, middle = mean(length)) +
  #   ylim(lymin,lyheight+50) +
  #   #geom_segment(aes(x=0,xend=11,y=100,yend=100))+
  #   scale_fill_manual(values=c(fillnoAnthocyanin,fillyesAnthocyanin))+
  #   geom_boxplot(position = "dodge", width = .9, aes(x = seedType, y = length, fill = redPigment, middle = mean(length))
  #   )+
  #   geom_text(data = CCTk, aes(x=seedType, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
  #   labs(x = "", y = ylab, title = "", 
  #        subtitle = "All repetitions", 
  #        fill = "Red Pigment: ", color = "") + 
  #   stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
  #   theme_minimal() +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = "none",
  #         legend.direction = "horizontal",
  #         legend.box = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0("All-coleoptiles-standardized",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0("All-coleoptiles-standardized",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # csvFilename <- paste0(tablesDir,"/All-coleoptiles-standardized",".csv",sep="")
  # write.csv(combinedData,csvFilename,row.names=FALSE)
  # 
  # stage6Announce <- paste0(filename," and corresponding csv finished and saved.",sep="")
  # print(stage6Announce)
  # Sys.sleep(.25)
  # 
  # 
  # # 4.5 Combined Total Length Graph ----
  # tukeyTotal <- combinedData[grepl("total",combinedData$tissue),]
  # tukeyTotal <- tukeyTotal[!grepl("control",tukeyTotal$isolate),]
  # tukeyTotal$seedType = with(tukeyTotal, reorder(seedType, -length, mean))
  # lyheight <- max(tukeyTotal$length)
  # lymin <- min(tukeyTotal$length)
  # # Significance:
  # TotalAnova <- aov(length~seedType,tukeyTotal)
  # TotalTukey <- TukeyHSD(TotalAnova)
  # tukeyTotal
  # cld <- multcompLetters4(TotalAnova,TotalTukey)
  # cld
  # CTTK <- group_by(tukeyTotal,seedType) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # CTTK
  # cld <- as.data.frame.list(cld$seedType)
  # cld
  # 
  # CTTK$cld <- cld$Letters
  # CTTK$redPigment <- "y"
  # # 4.5.1 Graph ----
  # ylab = "Total tissue lengths standardized to control"
  # tukeyTotal %>%
  #   filter(tissue %in% "total") %>%
  #   filter(!(isolate %in% "control")) %>%
  #   ggplot(outlier.shape = NA) +     
  #   geom_hline(yintercept=100, color = "#999999") +
  #   aes(x = seedType, y = length, fill = redPigment, drop=TRUE, middle = mean(length)) +
  #   ylim(lymin,lyheight+50) +
  #   #geom_segment(aes(x=0,xend=11,y=100,yend=100))+
  #   scale_fill_manual(values=c(fillnoAnthocyanin,fillyesAnthocyanin))+
  #   geom_boxplot(position = "dodge", width = .9, aes(x = seedType, y = length, fill = redPigment, middle = mean(length))
  #   )+
  #   geom_text(data = CTTK, aes(x=seedType, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
  #   labs(x = "", y = ylab, title = "", 
  #        subtitle = "All repetitions", 
  #        fill = "Red Pigment: ", color = "") + 
  #   stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
  #   theme_minimal() +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = "none",
  #         legend.direction = "horizontal",
  #         legend.box = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0("Combined-alltissues-standardized",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0("Combined-alltissues-standardized",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # csvFilename <- paste0(tablesDir,"/",assayNum,"-alltissues-standardized",".csv",sep="")
  # write.csv(combinedData,csvFilename,row.names=FALSE)
  # rm(cld)
  # 
  # # 4.5.2 Separated Total Length Graph ----
  # 
  # tukeyAll <- combinedData
  # tukeyAll <- tukeyAll[!grepl("control",tukeyAll$isolate),]
  # tukeyAll$seedType = with(tukeyAll, reorder(seedType, -length, mean))
  # lyheight <- max(tukeyAll$length)
  # lymin <- min(tukeyAll$length)
  # 
  # AllAnova <- aov(length~seedType,tukeyAll)
  # AllTukey <- TukeyHSD(AllAnova)
  # cld <- multcompLetters4(AllAnova,AllTukey)
  # cld
  # CATK <- group_by(tukeyAll,seedType) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # cld <- as.data.frame.list(cld$seedType)
  # CATK$cld <- cld$Letters
  # CATK$redPigment <- "y"
  # CATK$tissue <- "radicle"
  # # 4.5.3 Graph ----
  # ylab = "All tissue lengths standardized to control, by tissue"
  # tukeyAll %>%
  #   filter(!(isolate %in% "control")) %>%
  #   ggplot(outlier.shape = NA) +     
  #   geom_hline(yintercept=100, color = "#999999") +
  #   aes(x = seedType, y = length, fill = redPigment, colour = tissue, drop=TRUE, middle = mean(length)) +
  #   ylim(lymin,lyheight+50) +
  #   #geom_segment(aes(x=0,xend=11,y=100,yend=100))+
  #   scale_fill_manual(values=c(fillnoAnthocyanin,fillyesAnthocyanin))+
  #   scale_color_manual(values=c("#a0a0a0","#000000","#1111f1")) +
  #   geom_boxplot(position = "dodge", width = .9, aes(x = seedType, y = length, fill = redPigment, colour = tissue, middle = mean(length))
  #   )+
  #   geom_text(data = CATK, aes(x=seedType, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
  #   labs(x = "", y = ylab, title = "", 
  #        subtitle = "All repetitions", 
  #        fill = "Red Pigment: ", color = "") + 
  #   stat_summary(fun = mean, size=5, geom = "point", col = "black", pch=4, position = position_dodge(width = 1)) +  # Add points to plot
  #   theme_minimal() +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = c(0.5, 1.03),
  #         legend.box = "horizontal",
  #         legend.direction = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0("Combined-alltissues-divided-standardized",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0("Combined-alltissues-divided-standardized",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # csvFilename <- paste0(tablesDir,"/",assayNum,"-alltissues-divided-standardized",".csv",sep="")
  # write.csv(combinedData,csvFilename,row.names=FALSE)
  # rm(cld)
  # 
  # stage4Announce <- paste0(filename," and corresponding csv finished and saved.",sep="")
  # print(stage4Announce)
  # Sys.sleep(.25)
  # # 4.6 Variance Analysis
  # 
  # # 4.6.1 Filtering out false flags ----
  # ## Including controls could give us false significance, and totals are just 
  # ## combined lengths of coleoptile and radicle and are better analyzed separately
  # ControlData <- combinedData[grepl("control",combinedData$isolate),]
  # combinedData <- combinedData[!grepl("total",combinedData$tissue),]
  # combinedData <- combinedData[!grepl("control",combinedData$isolate),]
  # 
  # 
  # 4.6.2 Radicle Length Analysis ----
  NoRadCom <- combinedData[!grepl("control",combinedData$isolate),]
  NoRadCom <- NoRadCom[!grepl("radicle",NoRadCom$tissue),]
  NoColeCom <- combinedData[!grepl("control",combinedData$isolate),]
  NoColeCom <- NoColeCom[!grepl("coleoptile",NoColeCom$tissue),]

  BlockAnova <- aov(length~block,NoColeCom)
  BlockTukey <- TukeyHSD(BlockAnova)
  cld <- multcompLetters4(BlockAnova,BlockTukey)
  CBTK <- group_by(NoColeCom,block) %>%
    dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
    arrange(desc(sd))
  cld <- as.data.frame.list(cld$block)
  CBTK$cld <- cld$Letters
  CBTK$redPigment <- "y"
  CBTK$tissue <- "radicle"

  NoColeCom %>%
    filter(tissue %in% "radicle") %>%
    ggplot(outlier.shape = NA) +
    aes(x = block, y = length, fill = block) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    labs(y = "Radicle lengths by block") +
    geom_text(data = CBTK, aes(x=, y= max+min+sd-median-15, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
    geom_point(shape = 25, color = "black", fill = NA, size = 0.5, alpha = 1, position = position_jitter(width = 0.22)) +
    theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
          panel.grid.major.y=element_line(colour = "#f0f0f0"),
          panel.grid.major.x=element_line(colour = "#cecece"),
          #axis.ticks=element_line(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.subtitle = element_text(hjust=0, size = fontSize/2),
          #legend.position = c(0.5, 1.03),
          legend.position = c(0.5, 1.03),
          legend.box = "horizontal",
          legend.direction = "horizontal",
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size=fontSize),
          text = element_text(size = fontSize),
          axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))

  Sys.sleep(.1)
  filename <- paste0("Combined","-radicle-lengths-block",".jpg",sep="")
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0("Combined","-radicle-lengths-block",".wmf",sep="")
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  rm(cld)
  # 
  # # 4.6.3 Coleoptile Length Analysis ----
  # CBlockAnova <- aov(length~block,NoRadCom)
  # CBlockTukey <- TukeyHSD(CBlockAnova)
  # cld <- multcompLetters4(CBlockAnova,CBlockTukey)
  # cld
  # CCBTK <- group_by(NoRadCom,block) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # cld <- as.data.frame.list(cld$block)
  # CCBTK$cld <- cld$Letters
  # CCBTK$redPigment <- "y"
  # CCBTK$tissue <- "radicle" # Setting this to radicle will reveal if any data wasn't properly filtered out
  # 
  # NoRadCom %>%
  #   filter(tissue %in% "coleoptile") %>%
  #   ggplot(outlier.shape = NA) +
  #   aes(x = block, y = length, fill = block) +
  #   geom_boxplot() +
  #   scale_fill_hue(direction = 1) +
  #   labs(y = "Coleoptile lengths by block") +
  #   geom_text(data = CCBTK, aes(x=, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = c(0.5, 1.03),
  #         legend.box = "horizontal",
  #         legend.direction = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0("Combined","-coleoptile-lengths-block",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0("Combined","-coleoptile-lengths-block",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # rm(cld)
  # 
  # # 4.6.4 Exploratory Analysis of Standards/Controls ----
  NoRadComC <- combinedControls
  NoRadComC <- NoRadComC[!grepl("radicle",NoRadComC$tissue),]
  NoColeComC <- combinedControls
  NoColeComC <- NoColeComC[!grepl("coleoptile",NoColeComC$tissue),]

  CtrlBlockAnova <- aov(length~block,NoColeComC)
  CtrlBlockTukey <- TukeyHSD(CtrlBlockAnova)
  cld <- multcompLetters4(CtrlBlockAnova,CtrlBlockTukey)
  cld
  CtrlCBTK <- group_by(NoColeComC,block) %>%
    dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
    arrange(desc(sd))
  cld <- as.data.frame.list(cld$block)
  CtrlCBTK$cld <- cld$Letters
  CtrlCBTK$redPigment <- "y"
  CtrlCBTK$tissue <- "radicle"

  NoColeComC %>%
    filter(tissue %in% "radicle") %>%
    ggplot(outlier.shape = NA) +
    aes(x = block, y = length, fill = block) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    labs(y = "Radicle lengths by Control block") +
    geom_text(data = CtrlCBTK, aes(x=, y= mean, label = cld), position = position_jitter(height=5, width=0), size = fontSize/3, vjust=-4) +
    geom_point(shape = 25, color = "black", fill = NA, size = 0.5, alpha = 1, position = position_jitter(width = 0.22)) +
    theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
          panel.grid.major.y=element_line(colour = "#f0f0f0"),
          panel.grid.major.x=element_line(colour = "#cecece"),
          #axis.ticks=element_line(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.subtitle = element_text(hjust=0, size = fontSize/2),
          #legend.position = c(0.5, 1.03),
          legend.position = c(0.5, 1.03),
          legend.box = "horizontal",
          legend.direction = "horizontal",
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size=fontSize),
          text = element_text(size = fontSize),
          axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))

  Sys.sleep(.1)
  filename <- paste0("Combined","-radicle-lengths-Controlblock",".jpg",sep="")
  ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  filename2 <- paste0("Combined","-radicle-lengths-Controlblock",".wmf",sep="")
  ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  rm(cld)
  # 
  # # 4.6.5 Exploratory Analysis of Coleoptiles ----
  # 
  # CCtrlBlockAnova <- aov(length~block,NoRadComC)
  # CCtrlBlockTukey <- TukeyHSD(CCtrlBlockAnova)
  # cld <- multcompLetters4(CCtrlBlockAnova,CCtrlBlockTukey)
  # cld
  # CCtrlCBTK <- group_by(NoRadComC,block) %>%
  #   dplyr::summarise(mean=mean(length), sd = sd(length), max=max(length), median=median(length), min=min(length)) %>%
  #   arrange(desc(sd))
  # cld <- as.data.frame.list(cld$block)
  # CCtrlCBTK$cld <- cld$Letters
  # CCtrlCBTK$redPigment <- "y"
  # CCtrlCBTK$tissue <- "radicle"
  # 
  # NoRadComC %>%
  #   filter(tissue %in% "coleoptile") %>%
  #   ggplot(outlier.shape = NA) +
  #   aes(x = block, y = length, fill = block) +
  #   geom_boxplot() +
  #   scale_fill_hue(direction = 1) +
  #   labs(y = "Coleoptile lengths by Control block") +
  #   geom_text(data = CCtrlCBTK, aes(x=, y= max+min+sd-median-5, label = cld), position = position_jitter(width=0), size = fontSize/3, vjust=-4) +
  #   theme(panel.grid.minor=element_line(colour = "#f9f9f9"),
  #         panel.grid.major.y=element_line(colour = "#f0f0f0"),
  #         panel.grid.major.x=element_line(colour = "#cecece"),
  #         #axis.ticks=element_line(),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1.5),
  #         plot.subtitle = element_text(hjust=0, size = fontSize/2),
  #         #legend.position = c(0.5, 1.03),
  #         legend.position = c(0.5, 1.03),
  #         legend.box = "horizontal",
  #         legend.direction = "horizontal",
  #         legend.key.size = unit(0.75, 'cm'),
  #         legend.text = element_text(size=fontSize),
  #         text = element_text(size = fontSize),
  #         axis.text.x = element_text(face="bold", color="black", angle = 90, vjust=.25, size = fontSize))
  # 
  # Sys.sleep(.1)
  # filename <- paste0("Combined","-coleoptile-lengths-Controlblock",".jpg",sep="") 
  # ggsave(path=graphDir,filename,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # filename2 <- paste0("Combined","-coleoptile-lengths-Controlblock",".wmf",sep="") 
  # ggsave(path=graphDir,filename2,width=ChartWidth,height=ChartHeight,units="px",dpi=324,bg="#FFFEFE",limitsize=FALSE)
  # rm(cld)
  # 
  # 
  # # 4.6.5 ANOVA of Combined w/Blocks ----
  # anovaDoc <- paste0(anovaDir,"/","Combined","block+seedType-RadiclesOnly.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType, NoColeCom)),file=anovaDoc) 
  # 
  # anovaDoc <- paste0(anovaDir,"/","Combined","block+seedType-ColeoptilesOnly.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType, NoRadCom)),file=anovaDoc) 
  # 
  combinedTRootsOnly
  summary(aov(length ~ block+maizeLine:seedType, combinedTRootsOnly))
  
  
  # anovaDoc <- paste0(anovaDir,"/","Combined","block+seedType+tissue.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType+tissue, combinedData)),file=anovaDoc)
  # 
  # # 4.6.51 ANOVA of Combined w/Assays ----
  # anovaDoc <- paste0(anovaDir,"/","Combined","assay+seedType-RadiclesOnly.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType+assay, NoColeCom)),file=anovaDoc) 
  # 
  # anovaDoc <- paste0(anovaDir,"/","Combined","assay+seedType-ColeoptilesOnly.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType+assay, NoRadCom)),file=anovaDoc) 
  # 
  # anovaDoc <- paste0(anovaDir,"/","Combined","assay+seedType+tissue.doc",sep="")
  # capture.output("Radicles+Block+SeedType:", "",summary(aov(length ~ block+seedType+tissue+assay, combinedData)),file=anovaDoc) 
  # 
  # # 4.6.52 Comparing Anthocyanin Significance ----
  # anovaDoc <- paste0(anovaDir,"/","Combined","redPigment+seedType-RadiclesOnly.doc",sep="")
  # capture.output("Radicles+redPigment+SeedType:", "",summary(aov(length ~ assay+redPigment+seedType, NoColeCom)),file=anovaDoc) 
  # 
  # anovaDoc <- paste0(anovaDir,"/","Combined","redPigment+seedType-ColeoptilesOnly.doc",sep="")
  # capture.output("Radicles+redPigment+SeedType:", "",summary(aov(length ~ assay+redPigment+seedType, NoRadCom)),file=anovaDoc) 
  # 
  # anovaDoc <- paste0(anovaDir,"/","Combined","redPigment+seedType+tissue.doc",sep="")
  # capture.output("Radicles+redPigment+SeedType:", "",summary(aov(length ~ assay+redPigment+seedType+tissue+tissue:redPigment, combinedData)),file=anovaDoc) 
  # 
  # # 4.6.6 Pairwise Analysis to double check results manually ----
  # 
  # av8 <- aov(length~seedType*tissue
  #            ,combinedData)
  # 
  # TukeyDoc <- paste0(anovaDir,"/",n," Tukey-Results.doc",sep="")
  # capture.output(n," Tukey:", "",TukeyHSD(av8),file=TukeyDoc)
  # 
  # av9 <- aov(length~seedType*tissue
  #            ,combinedData)
  # 
  # lsm <- emmeans(av9, ~ seedType*tissue)
  # PairwiseDoc <- paste0(anovaDir,"/",n," Pairwise-Results.csv",sep="")
  # ptest3 <- contrast(lsm, method = "pairwise", adjust = "none")
  # ptest4 <- data.frame(ptest3)
  # 
  # write.csv(ptest4,PairwiseDoc,row.names=FALSE)
  # 
  # 
  # Sys.sleep(.1)
  # 
}

####################################sec5.0######################################
##                                                                            ##
## Harrison Hall (he/him)                                                     ##
## hphall2@illinois.edu                                                       ##
##                                                                            ##
## -When in doubt, restart-                                                   ##
##                                                                            ##
## When the data looks weird, first                                           ##
## ask what it is telling you, not                                            ##
## what is wrong. - Sarah Lipps                                               ##
##                                                                            ##
################################################################################
## DONE! ##
## 