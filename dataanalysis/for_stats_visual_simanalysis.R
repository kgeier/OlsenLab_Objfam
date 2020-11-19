#loading packages####
library("devtools", lib.loc="~/R/R-3.3.2/library")
library(stringr)
install.packages("tibble")
install.packages("imager")
install_github("bbuchsbaum/eyesim")
library(eyesim)
library(dplyr)
library(tibble)
library(rimage)


#Prepping Data for STS####

#load df data generated from 'feature_generation'
df <- read.csv("removed")
df$RepID <- gsub('T', "t", df$RepID) #some RepID were "Test" not "test"
df$Subject = str_pad(df$Subject, 2, pad="0") #add 0 in front of SID just to be consistent with msgrprt

#change filename back to _1 to match msgreport
df1 = subset(df, df$RepID != "test")
df1$imaga1 = substr(df1$image, 7,7)
df1$imaga1 = 1
df1$imaga2 = substr(df1$image, 1, 6)
df1$image = paste0(df1$imaga2, df1$imaga1)

df1 = df1[, !(colnames(df1) %in% c("imaga1", "imaga2", "subimg"))]

df2 = subset(df, df$RepID == "test")
df2 = df2[, !(colnames(df2) %in% c("subimg"))]

df = rbind(df1, df2, fill=TRUE)

#extract info from msgrprt to get fix-start time 
msgrprt <- read.csv("removed")

msgreport <- subset(msgrprt, subset = msgrprt$CURRENT_MSG_TEXT == "Study_display" | msgrprt$CURRENT_MSG_TEXT == "test_display" )
rm(msgrprt)

colnames(msgreport) <- c("Subject", "CURRENT_MSG_TEXT",        "CURRENT_MSG_TIME",        "trial",            
                         "block" ,                  "image",          "task")
msgreport <- msgreport[,c('Subject', 'image', 'trial', "CURRENT_MSG_TIME")]
msgreport$image = gsub(".jpg", "", msgreport$image)
msgreport$Subject = substr(msgreport$Subject, 3,4)

#merge msgrprt and df to get updated data
dfWithStartTime <- merge(df, msgreport, by = c("Subject", "image", "trial"))

dfWithStartTime$FixStart <- dfWithStartTime$Fixtime - dfWithStartTime$CURRENT_MSG_TIME 
dfWithStartTime$Subject <- as.factor(dfWithStartTime$Subject)
df <- dfWithStartTime

#change filename back to 1,2,3
df$newimgname = substr(df$image, 7, 7)

df$newimgname <- ifelse(df$SubsequentMorph == "Old", "1", 
                           ifelse(df$SubsequentMorph == "Near", "2",
                                  ifelse(df$SubsequentMorph == "Far", "3", "YEAHHH")))

df$newimg = substr(df$image, 1, 6)

df$image = paste0(df$newimg, df$newimgname)

df = df[,!(colnames(df) %in% c("newimg", "newimgname"))]

df$imgGroup = substr(df$image, 1,5)

badtrials = read.csv("removed")
df = df[df$image %in% badtrials$image_filename,]

##study-test similarity####

#EyeSim####
#NEW METHOD WITHOUT CENTRE BIAS####
## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_study <- eye_table("FixX", "FixY", "FixDur", onset="FixStart", groupvar=c("Subject", "imgGroup", "RepID"),
                        vars=c("trial", "Accuracy", "acc", "SubsequentMorph", "test_confidence", "RepID"),
                        data=subset(df, RepID != "test"), clip_bounds=c(212,812, 84, 684))

## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_test<- eye_table("FixX", "FixY", "FixDur", onset="FixStart", groupvar=c("Subject", "imgGroup", "RepID"),
                      vars=c("trial", "Accuracy", "acc", "SubsequentMorph", "test_confidence", "RepID"),
                      data=subset(df, RepID == "test"), clip_bounds=c(212,812, 84, 684))

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
edens_study <- density_by(etab_study, groups=c("Subject", "imgGroup", "Accuracy", "acc", "SubsequentMorph"))

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
edens_test <- density_by(etab_test, groups=c("Subject", "imgGroup", "Accuracy", "acc", "SubsequentMorph"))

edens_study$Subj_Image <- paste0(edens_study$Subject, "_", edens_study$imgGroup)
edens_test$Subj_Image <- paste0(edens_test$Subject, "_", edens_test$imgGroup)

tsim <- template_similarity(edens_study, edens_test, match_on="Subj_Image", method="spearman", permutations=50)

              
#remove small bins & create data file for stats - study-test similarity
tSIM <- tsim
umum  <-  tSIM %>% group_by(Accuracy, Subject) %>% summarise(count = length(acc))
umumSmallBins <- umum[umum$count <2.5 ,]
umumSmallBins$Remove <- "removeme"
tsimRem <- merge(umumSmallBins, tSIM, by = c("Subject", "Accuracy"), all =T)
tSIM <- subset(tsimRem, subset = is.na(tsimRem$Remove))
tsim <- tSIM[, -which(names(tSIM) %in% c("Remove", "count"))]

tsim = tsim[,c("Subject", "Accuracy", "imgGroup", "acc", "SubsequentMorph", "Subj_Image", "eye_sim_diff", "eye_sim", "perm_sim")]


#create data file for visualization - study-test similarity
d_sum6 = tsim

unloadNamespace("dplyr")

d_sum6.violin = summarySE(d_sum6, measurevar= "eye_sim_diff", groupvars=c("SubsequentMorph", "Accuracy", "Subject"), na.rm=TRUE)  




##prep data for SSS####

#load df data from 'objfam_mainfeatureandsim_feb2020'
df <- read.csv("removed")
df$Subject = str_pad(df$Subject, 2, pad="0") #add 0 in front of SID just to be consistent with msgrprt

#change filename back to _1 to match msgreport
df$imaga1 = substr(df$image, 7,7)
df$imaga1 = 1
df$image2 = substr(df$image, 1, 6)
df$image = paste0(df$image2, df$imaga1)

df = df[, !(colnames(df) %in% c("imaga1", "image2"))]

#extract info from msgrprt to get fix-start time 
msgrprt <- read.csv("removed") #Reading in msg rprt so can have fixation start time be matched to trial start time. 
msgreport <- subset(msgrprt, subset = msgrprt$CURRENT_MSG_TEXT == "Study_display" )
rm(msgrprt)

colnames(msgreport) <- c("Subject", "CURRENT_MSG_TEXT",        "CURRENT_MSG_TIME",        "trial",            
                         "block" ,                  "image",          "task")
msgreport <- msgreport[,c('Subject', 'image', 'trial', "CURRENT_MSG_TIME")]
msgreport$image = gsub(".jpg", "", msgreport$image)
msgreport$Subject = substr(msgreport$Subject, 3,4); 

#merge msgrprt and df to get updated data
dfWithStartTime <- merge(df, msgreport, by = c("Subject", "image", "trial"))

dfWithStartTime$FixStart <- dfWithStartTime$Fixtime - dfWithStartTime$CURRENT_MSG_TIME 
dfWithStartTime$Subject <- as.factor(dfWithStartTime$Subject)
df <- dfWithStartTime

df$imgGroup = substr(df$image, 1,5)

df = df[df$image %in% badtrials$image_filename,]


library(dplyr)
#Without Centre bias study-study#########
## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_study <- eye_table("FixX", "FixY", "FixDur", onset="FixStart", groupvar=c("Subject", "imgGroup", "RepID"),
                        vars=c("trial", "Accuracy", "acc", "SubsequentMorph", "test_confidence", "RepID"),
                        data=df, clip_bounds=c(212,812, 84, 684))

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
edens_study_eachrep <- density_by(etab_study, groups=c("Subject", "imgGroup", "RepID", "Accuracy", "acc", "SubsequentMorph"))

edens_study_eachrep$Subj_Image <- paste0(edens_study_eachrep$Subject, "_", edens_study_eachrep$imgGroup)
edens_rep1 <- subset(edens_study_eachrep, RepID == 1)
edens_rep2 <- subset(edens_study_eachrep, RepID == 2)
edens_rep3 <- subset(edens_study_eachrep, RepID == 3)

EncSim12 <- template_similarity( edens_rep1,edens_rep2, match_on="Subj_Image", method="spearman", permutations=50)
EncSim23 <- template_similarity( edens_rep2,edens_rep3, match_on="Subj_Image", method="spearman", permutations=50)
EncSim13 <- template_similarity( edens_rep1,edens_rep3, match_on="Subj_Image", method="spearman", permutations=50)

EncSim3reps <- merge(EncSim12, EncSim13, by = c("Subj_Image", "Subject", "imgGroup", "Accuracy", "acc", "SubsequentMorph"))
EncSim3reps <- merge(EncSim3reps, EncSim23, by = c("Subj_Image", "Subject", "imgGroup", "Accuracy", "acc", "SubsequentMorph"))

EncSim3reps$esd_avg <- (EncSim3reps$eye_sim_diff + EncSim3reps$eye_sim_diff.x + EncSim3reps$eye_sim_diff.y) /3

#remove small bins & create data file for stats - study-study similarity
ENCSim3reps <- EncSim3reps
umum  <-  ENCSim3reps %>% group_by(Accuracy, Subject) %>% summarise(count = length(Accuracy))
umumSmallBins <- umum[umum$count <2.5 ,]
umumSmallBins$Remove <- "removeme"
EncSim3repsRem <- merge(umumSmallBins, ENCSim3reps, by = c("Subject", "Accuracy"), all =T)
ENCSim3reps <- subset(EncSim3repsRem, subset = is.na(EncSim3repsRem$Remove))
EncSim3reps <- ENCSim3reps[, -which(names(ENCSim3reps) %in% c("Remove", "count"))]

EncSim3reps = EncSim3reps[,c("Subject", "Accuracy", "imgGroup", "acc", "SubsequentMorph", "Subj_Image", "esd_avg")]


#create data file for visualization - study-study similarity
d_sum5 = EncSim3reps

unloadNamespace("dplyr")

d_sum5.violin = summarySE(d_sum5, measurevar= "esd_avg", groupvars=c("SubsequentMorph", "Accuracy", "Subject"), na.rm=TRUE)      


