#Loading libraries
library("broom", lib.loc="~/R/R-3.3.1/library")
library("tidyr", lib.loc="~/R/R-3.3.1/library")
library("ggplot2", lib.loc="~/R/R-3.3.1/library")
library("gridExtra", lib.loc="~/R/R-3.3.1/library")
library("dplyr", lib.loc="~/R/R-3.3.1/library")
library("tidyr", lib.loc="~/R/R-3.3.2/library")
library("plotly", lib.loc="~/R/R-3.3.2/library")
library(reshape2)
library("lmerTest", lib.loc="~/R/R-3.3.2/library")
library("lme4", lib.loc="~/R/R-3.3.2/library")
library(car)
library(Rmisc)
library(stringr)
#setting a few values for reuse
pd <- position_dodge(0.1) # move them .05 to the left and right
theme_set(theme_bw(base_size = 18))

#where this can be backed up to github (if create a new version linked through git)
#https://github.com/kgeier/Object-Familiarity-Study-fMRI-and-EyeTracking-.git

#CODE START####

# define some paths: removed







#part 1 - main features: fixation count and sampling variability 

# read in the data
df.main = read.csv(filepath1)

# missing data? 
df.main <- subset(df.main, subset = df.main$SubsequentMorph != '.')
df.main$SID = gsub("_", "", df.main$SID)

# Convert conf to categorical
df.main$Confbinary = df.main$test_confidence.y
df.main[df.main$test_confidence.y < 3, ]$Confbinary = 'new'
df.main[df.main$test_confidence.y > 3, ]$Confbinary = 'old'
df.main$Confbinary = factor(df.main$Confbinary)

#Eliminate bins with fewer than x(9) trials in each bin.
library(dplyr)
#Grouping (Note, session label 9b is blocks 10-12 of subject 9. Next bit of code groups by SID not Recording session label)
main.df.table_binsize <- df.main  %>% group_by(Accuracy, SubsequentMorph, SID) %>%  summarise(count = length(SubsequentAccuracy))
main.littleBinsonly <- main.df.table_binsize[main.df.table_binsize$count < 8.5 ,]
#Making a shorter more managable version of df with just the 2 Relevant columns (so cleaner when merging with main df)
main.littleBinsonly <- main.littleBinsonly[, c("SID", "Accuracy", "SubsequentMorph")]
main.littleBinsonly$Remove <- "removeme"
#Merging the count df to main df by SID and subq acc
df.mainRem <- merge(main.littleBinsonly, df.main, by = c('SID', "Accuracy", "SubsequentMorph"), all = T)

dfs.main <- subset(df.mainRem, subset = is.na(df.mainRem$Remove))
rm(df.mainRem)
rm(df.main)
rm(main.littleBinsonly)
rm(main.df.table_binsize)

#analysis of features####
unloadNamespace("plotly")
unloadNamespace("broom")
library(lme4)
library(lmerTest)
library(dplyr)

d <- dfs.main

for (i in 6:(ncol(dfs.main)-3)){
  print(i)
  dfs.main[i] <- as.numeric(as.character(dfs.main[,i]))
}

#Reformat data 
names(dfs.main)[names(dfs.main) == 'SID'] <- 'Subject.ID'
library(dplyr)
dfs.main$Subject.ID = as.numeric(as.character(dfs.main$Subject.ID))
dfs.main$image_filename = ifelse(dfs.main$Subject.ID <= 9, substr(dfs.main$Sample_ID, 11,17), substr(dfs.main$Sample_ID, 12,18))
dfs.main$Subject.ID = gsub("_","", dfs.main$Subject.ID)
dfs.main$Subject.ID = str_pad(dfs.main$Subject.ID, 2, pad="0")
dfs.main$Subject.ID <- as.numeric(as.character(dfs.main$Subject.ID)); head(dfs.main)
dfs.main$SubsequentMorph <-  as.character(dfs.main$SubsequentMorph)
dfs.main$SubsequentMorph[dfs.main$SubsequentMorph == "Old"] <- "Identical" 

#making repnum column 
dfs.main$repnum = ifelse(dfs.main$Subject.ID <= 9, substr(dfs.main$Sample_ID, 25, 25), substr(dfs.main$Sample_ID, 26, 26))

#change filename back to 1,2,3
dfs.main$newimgname = substr(dfs.main$image_filename, 7, 7)

dfs.main$newimgname <- ifelse(dfs.main$SubsequentMorph == "Identical", "1", 
                        ifelse(dfs.main$SubsequentMorph == "Near", "2",
                               ifelse(dfs.main$SubsequentMorph == "Far", "3", "YEAHHH")))

dfs.main$newimg = substr(dfs.main$image_filename, 1, 6)

dfs.main$image_filename = paste0(dfs.main$newimg, dfs.main$newimgname)

dfs.main = dfs.main[,!(colnames(dfs.main) %in% c("newimg", "newimgname"))]

#removing 3 and no responses
dfs.main = dfs.main[dfs.main$SubsequentAccuracy != 'NoResponse'& dfs.main$Confbinary != '3',]


#remove outliers 

FCbytrial = dfs.main %>% 
  group_by (Subject.ID) %>%
  summarise (q1 = quantile(Fix_Count)[2], q3 = quantile(Fix_Count)[4]) %>% 
  group_by (Subject.ID) %>% 
  summarise (ProbOutB_FC = q1 - (1.5*(q3-q1)), ProbOutT_FC = q3 + (1.5*(q3-q1)))

TFDbytrial = dfs.main %>% 
  group_by (Subject.ID) %>%
  summarise (q1 = quantile(Total_Fix_Duration)[2], q3 = quantile(Total_Fix_Duration)[4]) %>% 
  group_by (Subject.ID) %>% 
  summarise (ProbOutB_TFD = q1 - (1.5*(q3-q1)), ProbOutT_TFD = q3 + (1.5*(q3-q1)))

BT = Reduce(function(x, y) merge(x, y, all=TRUE), list(FCbytrial, TFDbytrial)) #changeme: include TSDbytrial in bracket

FC= subset(dfs.main, select=c("Subject.ID", "Sample_ID", "repnum", "Fix_Count"))
TFD = subset(dfs.main, select=c("Subject.ID", "Sample_ID", "repnum", "Total_Fix_Duration"))

features = join_all(list(FC, TFD), by = c("Subject.ID", "Sample_ID", "repnum"), type='full'); head(features) #changeme: include TSD in list()

FBT=merge(features, BT, by="Subject.ID", all=TRUE); head(FBT)

FBT$goodtrial = ifelse(FBT$Fix_Count >= FBT$ProbOutB_FC & FBT$Fix_Count <= FBT$ProbOutT_FC
                       & FBT$Total_Fix_Duration >= FBT$ProbOutB_TFD & FBT$Total_Fix_Duration <= 2600, "good", "bad")
df.gb = FBT[,(colnames(FBT) %in% c("Subject.ID","Sample_ID", "repnum", "goodtrial"))]
dfs.mm = merge(df.gb, dfs.main, by = c("Subject.ID", "Sample_ID", "repnum"), all = TRUE); head(dfs.mm)
dfs.main = subset(dfs.mm, subset = dfs.mm$goodtrial == "good"); head(dfs.main)
dfs.main = dfs.main[, !(colnames(dfs.main) %in% c("Total_Saccade_Duration", "goodtrial"))]; head(dfs.main)


#between-subject analysis 
dfavg = dfs.main %>%
  group_by(Subject.ID) %>%
  summarise(fixm = mean(Fix_Count, na.rm=TRUE), sam = mean(JRobin_areaSampling, na.rm=TRUE))

dfavg1 = merge(dprime_far, dfavg, by = "Subject.ID", all=TRUE)

##bwsubj-Fixation count####
fcbwsubj = lm(formula = fixm ~ dprime_faronly_noguess, data = dfavg1)
summary(fcbwsubj)

##Fixation count plot####
fcbw = ggplot(dfavg1, aes(x=dprime_faronly_noguess, y=fixm))+ 
  geom_point(color = 'black') + 
  geom_smooth(method='lm', formula = y~x, se = TRUE, color = "darkred", fill="grey") + 
  labs(title="Fixation count by subject performance",
       x="d'", y="Fixation count") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
  theme(
    plot.title=element_text(size=20, face="bold", hjust=0.5),
      axis.title.x=element_text(size=18), 
      axis.title.y=element_text(size=18)) 

fcbwR = fcbw + annotate("text", x=1.8, y=7, label = "R² = 0.18")


##bwsubj-sampling variability#### 
svbwsubj = lm(formula = sam ~ dprime_faronly_noguess, data = dfavg1)
summary(svbwsubj)

#create data file for stats - fixatio count
library(plyr) 
unloadNamespace("dplyr")
d_sum1 = ddply(dfs.main,~SubsequentMorph + Confbinary + Accuracy + Subject.ID+image_filename,
              summarise, Fix_Count=mean(Fix_Count))

d_sum1 <- subset(d_sum1, !is.nan(get("Fix_Count")))


#create data file for visualization - fixation count
d_sumFC.violin = summarySE(dfs.main, measurevar="Fix_Count", groupvars=c("SubsequentMorph", "Accuracy", "Subject.ID"), na.rm=TRUE)      

#First steps are making column accuracy
d_sumFC.violin$AccuracyType2 <- paste0(d_sumFC.violin$SubsequentMorph,d_sumFC.violin$Accuracy)
#Making new column with names of accuracy type (eg. correct rejection Far)
d_sumFC.violin$AccuracyType <- ifelse(d_sumFC.violin$AccuracyType2 == "IdenticalCorrect", "Hit", 
                                      ifelse(d_sumFC.violin$AccuracyType2 == "IdenticalIncorrect", "Miss",
                                             ifelse(d_sumFC.violin$AccuracyType2 == "FarIncorrect", " FA ",
                                                    ifelse(d_sumFC.violin$AccuracyType2 == "FarCorrect", " CR ",
                                                           ifelse(d_sumFC.violin$AccuracyType2 == "NearIncorrect", "FA",
                                                                  ifelse(d_sumFC.violin$AccuracyType2 == "NearCorrect", "CR", "AGHHHWAHHDFH"))))))

#check duplicates
d_sum1$subimg = paste0(d_sum1$Subject.ID, sep = "_", d_sum1$image_filename)
unique(d_sum1$subimg)


#create data file for stats - sampling variability 
unloadNamespace("dplyr")
d_sum2 =  summarySEwithin(dfs.main,
                         measurevar="JRobin_areaSampling",
                         withinvars=c("SubsequentMorph", "Accuracy", "Subject.ID", "image_filename"),
                         na.rm = T)


d_sum2 <- subset(d_sum2, !is.nan(get("JRobin_areaSampling")))

##data file for visualization - sampling variability 
d_sumJR.boxplot = summarySEwithin(dfs.main, measurevar = "JRobin_areaSampling", withinvars=c("SubsequentMorph", "Accuracy", "Subject.ID"), na.rm=TRUE) #doesn't matter that its summarySE bc only using this for the "mean" (not CI etc.)

#First steps are making column accuracy
d_sumJR.boxplot$AccuracyType2 <- paste0(d_sumJR.boxplot$SubsequentMorph,d_sumJR.boxplot$Accuracy)
#Making new column with names of accuracy type (eg. correct rejection Far)
d_sumJR.boxplot$AccuracyType <- ifelse(d_sumJR.boxplot$AccuracyType2 == "IdenticalCorrect", "Hit", 
                                       ifelse(d_sumJR.boxplot$AccuracyType2 == "IdenticalIncorrect", "Miss",
                                              ifelse(d_sumJR.boxplot$AccuracyType2 == "FarIncorrect", " FA ",
                                                     ifelse(d_sumJR.boxplot$AccuracyType2 == "FarCorrect", " CR ",
                                                            ifelse(d_sumJR.boxplot$AccuracyType2 == "NearIncorrect", "FA",
                                                                   ifelse(d_sumJR.boxplot$AccuracyType2 == "NearCorrect", "CR", "AGHHHWAHHDFH"))))))


#check duplicates
d_sum2$subimg = paste0(d_sum2$Subject.ID, d_sum2$image_filename)
unique(d_sum2$subimg)











#part 2 - meta features: eye-movement repetition effect 1-3

# read in the data
df.meta = read.csv(filepath2) #1to3

#missing data?
df.meta <- subset(df.meta, subset = df.meta$SubsequentMorph != '.')

# Convert conf to categorical
df.meta$Confbinary = df.meta$test_confidence.y
df.meta[df.meta$test_confidence.y < 3, ]$Confbinary = 'new'
df.meta[df.meta$test_confidence.y > 3, ]$Confbinary = 'old'
df.meta$Confbinary = factor(df.meta$Confbinary)

df.meta$SID = gsub("_", "", df.meta$SID)

#Eliminate bins with fewer than x(9) trials in each bin
library(dplyr)
#Grouping (Note, session label 9b is blocks 10-12 of subject 9. Next bit of code groups by SID not Recording session label)
meta.df.table_binsize <- df.meta  %>% group_by(Accuracy, SubsequentMorph, SID) %>%  summarise(count = length(SubsequentAccuracy))
meta.littleBinsonly <- meta.df.table_binsize[meta.df.table_binsize$count < 2.5,]
#Making a shorter more managable version of df with just the 2 Relevant columns (so cleaner when merging with main df)
meta.littleBinsonly <- meta.littleBinsonly[, c("SID", "SubsequentMorph", "Accuracy")]
meta.littleBinsonly$Remove <- "removeme"
#Merging the count df to main df by SID and subq acc
df.metaRem <- merge(meta.littleBinsonly, df.meta, by = c('SID', "SubsequentMorph", "Accuracy"), all = T)

dfs.meta <- subset(df.metaRem, subset = is.na(df.metaRem$Remove))
rm(df.metaRem)
rm(df.meta)
rm(meta.littleBinsonly)
rm(meta.df.table_binsize)

#analysis of features
unloadNamespace("plotly")
unloadNamespace("broom")
library(lme4)
library(lmerTest)
library(dplyr)

dfs.meta[] <- lapply(dfs.meta, gsub, pattern = "Inf", replacement = NA, fixed = TRUE) #Getting rid of the few infinite values and replacing that cell with NA

for (i in 7:15){
  print(i)
  dfs.meta[i] <- as.numeric(as.character(dfs.meta[,i]))
}

fact <- c(2,3,16,17)
for (i in fact){
  print(i)
  dfs.meta[i] <- as.factor(dfs.meta[,i])
}

#reformat data
names(dfs.meta)[names(dfs.meta) == 'SID'] <- 'Subject.ID'

library(dplyr)
dfs.meta$Subject.ID = as.numeric(as.character(dfs.meta$Subject.ID))
dfs.meta$image_filename = ifelse(dfs.meta$Subject.ID <= 9, substr(dfs.meta$Sample_ID, 11,17), substr(dfs.meta$Sample_ID, 12,18))
dfs.meta$Subject.ID = gsub("_","", dfs.meta$Subject.ID)
dfs.meta$Subject.ID = str_pad(dfs.meta$Subject.ID, 2, pad="0")
dfs.meta$Subject.ID <- as.numeric(as.character(dfs.meta$Subject.ID)); head(dfs.meta)

dfs.meta$SubsequentMorph <-  as.character(dfs.meta$SubsequentMorph)
dfs.meta$SubsequentMorph[dfs.meta$SubsequentMorph == "Old"] <- "Identical"
dfs.meta$SubsequentMorph <-  as.factor(dfs.meta$SubsequentMorph)

head(dfs.meta)


#change filename back to 1,2,3
dfs.meta$newimgname = substr(dfs.meta$image_filename, 7, 7)

dfs.meta$newimgname <- ifelse(dfs.meta$SubsequentMorph == "Identical", "1", 
                              ifelse(dfs.meta$SubsequentMorph == "Near", "2",
                                     ifelse(dfs.meta$SubsequentMorph == "Far", "3", "YEAHHH")))

dfs.meta$newimg = substr(dfs.meta$image_filename, 1, 6)

dfs.meta$image_filename = paste0(dfs.meta$newimg, dfs.meta$newimgname)

dfs.meta = dfs.meta[,!(colnames(dfs.meta) %in% c("newimg", "newimgname"))]

#removing 3 and no response 
dfs.meta = dfs.meta[dfs.meta$SubsequentAccuracy != 'NoResponse'& dfs.meta$test_confidence.y != '3',]


#remove outliers
badtrials13 = df.gb %>% 
  subset(repnum == 1|repnum == 3)

badtrials13$ISID <- ifelse(badtrials13$Subject.ID <= 9, substr(badtrials13$Sample_ID, 1,17),substr(badtrials13$Sample_ID, 1,18))

trimbt <- badtrials13[,c("ISID", "repnum", "goodtrial")] 
library(tidyr)
wide <- spread(trimbt, repnum, goodtrial) ;head(wide)
library(data.table)
setnames(wide, c("1", "3"), c("rep1", "rep3"))
wide$badeither <- ifelse(wide$rep1 == "bad" | wide$rep3 == "bad", "bad", "good"); head(wide)

dfs.meta$ISID = ifelse(dfs.meta$Subject.ID <= 9, substr(dfs.meta$Sample_ID, 1,17),substr(dfs.meta$Sample_ID, 1,18))
dfs.meta = merge(dfs.meta, wide, by = "ISID", all=TRUE); head(dfs.meta) 
dfs.meta = dfs.meta[dfs.meta$badeither %in% 'good',]; sample_n(dfs.meta, 20)
dfs.meta = dfs.meta[,!(colnames(dfs.meta) %in% c("rep1", "rep3", "badeither"))]; head(dfs.meta)

##bwsubj-eye-movement repetition effect 1-3#### 
library(dplyr)
df13avg = dfs.meta %>% 
  group_by(Subject.ID) %>%
  summarise(df13avg = mean(FixCountDiff))

dfavg2 = merge(df13avg, dprime_far, by = "Subject.ID", all=TRUE)

df13bwsubj = lm(formula = df13avg ~ dprime_faronly_noguess, data = dfavg2)
summary(df13bwsubj)

#create data file for stats - eye-movement repetition effect 1-3
unloadNamespace("dplyr")
d_sum3 =  summarySEwithin(dfs.meta,
                         measurevar="FixCountDiff",
                         withinvars=c("SubsequentMorph", "Accuracy", "Subject.ID", "image_filename"),
                         na.rm = T)

d_sum3 <- subset(d_sum3, !is.nan("FixCountDiff"))
d_sum3 -> sc3meta

dfs.meta3 = dfs.meta

#create data for visualization - eye-movement repetition effect 1-3
d_sum4.violin = summarySE(dfs.meta3, measurevar="FixCountDiff", groupvars=c("SubsequentMorph", "Accuracy", "Subject.ID"), na.rm=TRUE)  
#First steps are making column accuracy
d_sum4.violin$AccuracyType2 <- paste0(d_sum4.violin$SubsequentMorph,d_sum4.violin$Accuracy) 
#Making new column with names of accuracy type (eg. correct rejection Far)
d_sum4.violin$AccuracyType <- ifelse(d_sum4.violin$AccuracyType2 == "IdenticalCorrect", "Hit", 
                                     ifelse(d_sum4.violin$AccuracyType2 == "IdenticalIncorrect", "Miss",
                                            ifelse(d_sum4.violin$AccuracyType2 == "FarIncorrect", " FA ",
                                                   ifelse(d_sum4.violin$AccuracyType2 == "FarCorrect", " CR ",
                                                          ifelse(d_sum4.violin$AccuracyType2 == "NearIncorrect", "FA",
                                                                 ifelse(d_sum4.violin$AccuracyType2 == "NearCorrect", "CR", "AGHHHWAHHDFH"))))))

#OPTIONAL IF WANT TO SWITCH TO MAGNITUDE OF REP EFFECT INSTEAD OF "change in fix count" ###
d_sum4.violin$MagRep <- d_sum4.violin$FixCountDiff * -1 


#check duplicates
d_sum3$subimg = paste0(d_sum3$Subject.ID, d_sum3$image_filename)
unique(d_sum3$subimg)










#part 3 - meta features: eye-movement repetition effect 1-2

# read in the data
df.meta = read.csv(filepath3) 

# missing data?
df.meta <- subset(df.meta, subset = df.meta$SubsequentMorph != '.')

# Convert conf to categorical
df.meta$Confbinary = df.meta$test_confidence.y
df.meta[df.meta$test_confidence.y < 3, ]$Confbinary = 'new'
df.meta[df.meta$test_confidence.y > 3, ]$Confbinary = 'old'
df.meta$Confbinary = factor(df.meta$Confbinary)


#Eliminate bins with fewer than x(9) trials in each bin
library(dplyr)
#Grouping (Note, session label 9b is blocks 10-12 of subject 9. Next bit of code groups by SID not Recording session label)
meta.df.table_binsize <- df.meta  %>% group_by(SubsequentMorph, Accuracy, SID) %>%  summarise(count = length(SubsequentAccuracy))

meta.littleBinsonly <- meta.df.table_binsize[meta.df.table_binsize$count < 2.5 ,] #only if meta each... is 5.5 right?
#Making a shorter more managable version of df with just the 2 Relevant columns (so cleaner when merging with main df)
meta.littleBinsonly <- meta.littleBinsonly[, c("SID", "SubsequentMorph", "Accuracy")]
meta.littleBinsonly$Remove <- "removeme"
meta.littleBinsonly <- as.data.frame(meta.littleBinsonly)
as.data.frame(50) -> a ; as.data.frame("NoResponse") -> c ; as.data.frame("removeme") -> d ; as.data.frame("nah") -> b
ThisStupid <- cbind(a,b,c,d)
colnames(ThisStupid) <- c("SID" ,    "SubsequentMorph" ,"Accuracy" ,   "Remove")
meta.littleBinsonly <- rbind(meta.littleBinsonly,ThisStupid)
#Merging the count df to main df by SID and subq acc
df.metaRem <- merge(meta.littleBinsonly, df.meta, by = c('SID', "SubsequentMorph", "Accuracy"), all = T)

dfs.meta <- df.metaRem
rm(df.meta)
rm(meta.littleBinsonly)
rm(meta.df.table_binsize)

#analysis of features####
unloadNamespace("plotly")
unloadNamespace("broom")
library(lme4)
library(lmerTest)
library(dplyr)


dfs.meta[] <- lapply(dfs.meta, gsub, pattern = "Inf", replacement = NA, fixed = TRUE) #Getting rid of the few infinite values and replacing that cell with NA


#reformat data
names(dfs.meta)[names(dfs.meta) == 'SID'] <- 'Subject.ID'

library(dplyr)
dfs.meta$Subject.ID = gsub("_", "", dfs.meta$Subject.ID)
dfs.meta$Subject.ID = as.numeric(as.character(dfs.meta$Subject.ID))
dfs.meta$image_filename = ifelse(dfs.meta$Subject.ID <= 9, substr(dfs.meta$Sample_ID, 11,17), substr(dfs.meta$Sample_ID, 12,18))
dfs.meta$Subject.ID = str_pad(dfs.meta$Subject.ID, 2, pad="0")
dfs.meta$Subject.ID <- as.numeric(as.character(dfs.meta$Subject.ID)); head(dfs.meta)

dfs.meta$SubsequentMorph <-  as.character(dfs.meta$SubsequentMorph)
dfs.meta$SubsequentMorph[dfs.meta$SubsequentMorph == "Old"] <- "Identical" #Changing morph from Old to Identical to distinguish from "Old" responses
dfs.meta$SubsequentMorph <-  as.factor(dfs.meta$SubsequentMorph)

head(dfs.meta)


for (i in 7:15){
  print(i)
  dfs.meta[i] <- as.numeric(as.character(dfs.meta[,i]))
}

fact <- c(2,3,16,17)
for (i in fact){
  print(i)
  dfs.meta[i] <- as.factor(dfs.meta[,i])
}

#change filename back to 1,2,3
dfs.meta$newimgname = substr(dfs.meta$image_filename, 7, 7)

dfs.meta$newimgname <- ifelse(dfs.meta$SubsequentMorph == "Identical", "1", 
                              ifelse(dfs.meta$SubsequentMorph == "Near", "2",
                                     ifelse(dfs.meta$SubsequentMorph == "Far", "3", "YEAHHH")))

dfs.meta$newimg = substr(dfs.meta$image_filename, 1, 6)

dfs.meta$image_filename = paste0(dfs.meta$newimg, dfs.meta$newimgname)

dfs.meta = dfs.meta[,!(colnames(dfs.meta) %in% c("newimg", "newimgname"))]

#removing 3 and no response
dfs.meta = dfs.meta[dfs.meta$SubsequentAccuracy != 'NoResponse'& dfs.meta$test_confidence.y != '3',]

#remove outliers

badtrials12 = df.gb %>% 
  subset(repnum == 1|repnum == 2)

badtrials12$ISID <- ifelse(badtrials12$Subject.ID <= 9, substr(badtrials12$Sample_ID, 1,17),substr(badtrials12$Sample_ID, 1,18))

trimbt <- badtrials12[,c("ISID", "repnum", "goodtrial")]
library(tidyr)
wide <- spread(trimbt, repnum, goodtrial) ;head(wide)
library(data.table)
setnames(wide, c("1", "2"), c("rep1", "rep2"))
wide$badeither <- ifelse(wide$rep1 == "bad" | wide$rep2 == "bad", "bad", "good"); head(wide)

#dfs.meta = dfs.meta[!dfs.meta$image_filename %in% badtrials$image_filename,]
dfs.meta$ISID = ifelse(dfs.meta$Subject.ID <= 9, substr(dfs.meta$Sample_ID, 1,17),substr(dfs.meta$Sample_ID, 1,18))
dfs.meta = merge(dfs.meta, wide, by = "ISID", all=TRUE); head(dfs.meta)
unloadNamespace("plyr")
dfs.meta = dfs.meta[dfs.meta$badeither %in% 'good',]
dfs.meta = dfs.meta[,!(colnames(dfs.meta) %in% c("rep1", "rep2", "badeither"))]; head(dfs.meta)

##bwsubj-eye-movement repetition effect 1-2#### 
df12avg = dfs.meta %>% 
  group_by(Subject.ID) %>%
  summarise(df12avg = mean(FixCountDiff))
dfavg3 = merge(df12avg, dprime_far, by = "Subject.ID", all=TRUE)

df12bwsubj = lm(formula = df12avg ~ dprime_faronly_noguess, data = dfavg3)
summary(df12bwsubj)

unloadNamespace("dplyr")

#create data file for stats - eye-movement repetition effect 1-2 
d_sum4 =  summarySEwithin(dfs.meta,
                         measurevar="FixCountDiff",
                         withinvars=c("SubsequentMorph", "Accuracy", "Subject.ID", "image_filename"),
                         na.rm = T)

d_sum4 <- subset(d_sum4, !is.nan(FixCountDiff))
d_sum4 -> sc2meta


#create data file for visualization - eye-movement repetition effect 1-2

dfs.meta3 <- dfs.meta

unloadNamespace("dplyr")

d_sum4.violin = summarySE(dfs.meta3, measurevar= "FixCountDiff", groupvars=c("SubsequentMorph", "Accuracy", "Subject.ID"), na.rm=TRUE)      

#Making new column with names of accuracy type (eg. correct rejection Far)
d_sum4.violin$AccuracyType2 <- paste0(d_sum4.violin$SubsequentMorph,d_sum4.violin$Accuracy)
d_sum4.violin$AccuracyType <- ifelse(d_sum4.violin$AccuracyType2 == "IdenticalCorrect", "Hit", 
                                     ifelse(d_sum4.violin$AccuracyType2 == "IdenticalIncorrect", "Miss",
                                            ifelse(d_sum4.violin$AccuracyType2 == "FarIncorrect", " FA ",
                                                   ifelse(d_sum4.violin$AccuracyType2 == "FarCorrect", " CR ",
                                                          ifelse(d_sum4.violin$AccuracyType2 == "NearIncorrect", "FA",
                                                                 ifelse(d_sum4.violin$AccuracyType2 == "NearCorrect", "CR", "AGHHHWAHHDFH"))))))



#OPTIONAL IF WANT TO SWITCH TO MAGNITUDE OF REP EFFECT INSTEAD OF "change in fix count" ###
d_sum4.violin$MagRep <- d_sum4.violin$FixCountDiff * -1 


#check duplicates
d_sum4$subimg = paste0(d_sum4$Subject.ID, d_sum4$image_filename)
unique(d_sum4$subimg)









#part 4 - repetition effect overall: data file for stats and visualization

dfs.main -> df_repef
head(df_repef)


d_rp = summarySE(df_repef, measurevar="Fix_Count", groupvars=c( "Subject.ID", "repnum"), na.rm=TRUE)      
d_rp13 <- subset(d_rp, subset = d_rp$repnum == 1)
d_rp13$RepComp <- "1 & 3"
d_rp$RepComp <- ifelse(d_rp$repnum == 3, "1 & 3", "1 & 2")
d_rp <- rbind(d_rp, d_rp13)


