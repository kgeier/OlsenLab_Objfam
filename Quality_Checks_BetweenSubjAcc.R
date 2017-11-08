#Loading libraries
library("broom", lib.loc="~/R/R-3.3.1/library")
library("tidyr", lib.loc="~/R/R-3.3.1/library")
library("ggplot2", lib.loc="~/R/R-3.3.1/library")
library("gridExtra", lib.loc="~/R/R-3.3.1/library")
library("dplyr", lib.loc="~/R/R-3.3.1/library")
library("tidyr", lib.loc="~/R/R-3.3.2/library")
library("plotly", lib.loc="~/R/R-3.3.2/library")
library(plyr)
library(reshape2)
#setting a few values for reuse
pd <- position_dodge(0.1) # move them .05 to the left and right
theme_set(theme_bw(base_size = 18))

#where this can be backed up to github (if create a new version linked through git)
#https://github.com/kgeier/Object-Familiarity-Study-fMRI-and-EyeTracking-.git

#Loading in 3 summary functions for SE etc. 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci<- datac$se * ciMult
  
  return(datac)
}
## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

# #Step 1: histograms #### 
# #of d prime based on generous or strict generation of dprimes. Genorous doesn't count 3's and "noresponses" as incorrect. It also only compares Old to Far, no Near.  
#   p3 <- ggplot(df_wide1, aes(x=d.)) +
#   geom_histogram(binwidth=.23, aes(fill=..x..)) +
#   geom_vline(aes(xintercept=mean(d.,na.rm=T)),
#              linetype="dashed",size=.4) +
#   theme_classic() +
#   scale_fill_gradient(name="d prime",low="purple", high = "darkorange", guide = FALSE) +
#   theme(axis.line.x = element_line(color = "black"), 
#         axis.line.y = element_line(color = "black")) +
#   theme(plot.title = element_text(size=36),
#         axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   scale_y_continuous(limits= c(-0.6, 9)) +
#   scale_x_continuous(limits= c(-0.5, 2.2)) +
#   ylab("Frequency") +
#   xlab("d prime") +
#   ggtitle("Distribution of d'\n Far Only") + theme(plot.title = element_text(hjust = 0.5, color ="#330066")) +
#     geom_text(aes(mean(d. ,na.rm = T),-0.5,label = "mean", hjust = -0.15))
# p3
# grid.arrange(p3,p4, p1,p2, ncol=2,nrow=2)
# 
# #graph should do confidence (1-5) by morph (1-3) (15 bars) 
# p5 <- ggplot(df_newFiltered, aes(x=averageconfidence_old, averageconfidence_near, averageconfidence_far)) +
#   geom_boxplot(aes(fill=..x..)) +
#   theme_classic() +
#   scale_fill_gradient(name="wedkkt",low="blue", high = "darkorange") +
#   theme(axis.line.x = element_line(color = "black"), 
#         axis.line.y = element_line(color = "black")) +
#   theme(plot.title = element_text(size=36),
#         axis.text = element_text(size=24), 
#         axis.title = element_text(size=32),
#         strip.text.x = element_text(size=16)) +
#   ylab("Morph") +
#   xlab("Average Confidence") +
#   ggtitle("confidence for morph") +
#   theme(plot.title = element_text(hjust=0.5, color = "#330066"))
# p5
# 
  
  #Redoing with new excel report####
  # define some paths
  data_dir = 'C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/EyeTracking/ObjFamETStudy/Kirk Graphs/Clean_forStephCodeTesting' 
  filename1 =  'AllSubj_studyPeriod.csv' 
  filename2 = 'AllSubj_testPeriod.csv'
  filepath1 = paste(data_dir, '/', filename1, sep = ''); filepath1
  filepath2 = paste(data_dir, '/', filename2, sep = ''); filepath2
  
  # read in the data
  df.rawstudyIP = read.csv(filepath1, stringsAsFactors = F)
  df.rawtestIP = read.csv(filepath2, stringsAsFactors = F)
  
  #working####
  #convert to numeric 
  df.rawtestIP$AVERAGE_BLINK_DURATION <- as.numeric(as.character(df.rawtestIP$AVERAGE_BLINK_DURATION))
  df.rawtestIP$BLINK_COUNT <-as.numeric(as.character(df.rawtestIP$BLINK_COUNT))
  df.rawtestIP$FIXATION_COUNT <-as.numeric(as.character(df.rawtestIP$FIXATION_COUNT))
  df.rawtestIP$AVERAGE_FIXATION_DURATION <-as.numeric(as.character(df.rawtestIP$AVERAGE_FIXATION_DURATION))
  df.rawtestIP$LAST_TRIAL_TIME <-as.numeric(as.character(df.rawtestIP$LAST_TRIAL_TIME))
  df.rawtestIP$X__Reaction_Time__1 <-as.numeric(as.character(df.rawtestIP$X__Reaction_Time__1))
  df.rawtestIP$displayTime <-as.numeric(as.character(df.rawtestIP$displayTime))
  class(df.rawtestIP$FIXATION_COUNT)
  #AHHHHH. There must be a better way...
  #convert to numeric 
  df.rawstudyIP$AVERAGE_BLINK_DURATION <- as.numeric(as.character(df.rawstudyIP$AVERAGE_BLINK_DURATION))
  df.rawstudyIP$BLINK_COUNT <-as.numeric(as.character(df.rawstudyIP$BLINK_COUNT))
  df.rawstudyIP$FIXATION_COUNT <-as.numeric(as.character(df.rawstudyIP$FIXATION_COUNT))
  df.rawstudyIP$AVERAGE_FIXATION_DURATION <-as.numeric(as.character(df.rawstudyIP$AVERAGE_FIXATION_DURATION))
  df.rawstudyIP$LAST_TRIAL_TIME <-as.numeric(as.character(df.rawstudyIP$LAST_TRIAL_TIME))
  df.rawstudyIP$X__Reaction_Time__1 <-as.numeric(as.character(df.rawstudyIP$X__Reaction_Time__1))
  df.rawstudyIP$displayTime <-as.numeric(as.character(df.rawstudyIP$displayTime))
  class(df.rawstudyIP$FIXATION_COUNT)
  #AHHHHH. There must be a better way...
  
  
  #Combine study and test Interest Periods####
  #Also gets rid of null (odd/even) trials

  df.rawstudyIP <- cbind(df.rawstudyIP, seq(1:nrow(df.rawstudyIP)))
  df.rawtestIP <-cbind(df.rawtestIP, seq(1:nrow(df.rawtestIP)))
  
  #making sure the new column is named 
  names(df.rawstudyIP)[names(df.rawstudyIP) == 'seq(1:nrow(df.rawstudyIP))'] <- 'orderer'
  names(df.rawtestIP)[names(df.rawtestIP) == 'seq(1:nrow(df.rawtestIP))'] <- 'orderer'
  
  #and numeric
  df.rawstudyIP$orderer <- as.numeric(as.character(df.rawstudyIP$orderer))
  df.rawtestIP$orderer <- as.numeric(as.character(df.rawtestIP$orderer))
  
  #subsetting #Confusing bc there are both Test and test phases. Also study phase still has null trials but labeled differently...
  subset_study <-  df.rawstudyIP[grepl("study", df.rawstudyIP$task),]
  subset_test <-  df.rawtestIP[grepl("Test", df.rawtestIP$task),]
  subset_Test <- df.rawtestIP[grepl("test", df.rawtestIP$task),]
  
  
  #rbind #This is the line that introduces the NA row. Still unclear why its introduced but it's only the 1 line and happens even when combining just 1 df. 
  df_rebound <- rbind(subset_study, subset_test, subset_Test,fill = F) 
  
  #Getting rid of 'fill' row from last command
  df_rebound <- subset(df_rebound, subset = df_rebound$orderer != 0)
  
  #re-ordering
  df_final <- df_rebound[order(df_rebound$orderer),]
  
  #removing Null trials. Regular expressions to target image name. BC "null" was not consistent between subjects for some reason
  df_final2 <-  df_final[!grepl("^[1-9].jpg", df_final$image_filename),]
  df_final3 <- df_final2[is.na(df_final2$image_filename) == FALSE,]
  
  #Converting recording session label to subject ID####
  #adding Subject ID to both study and test. 
  df_final3$SessionLabelCopy <- df_final3$RECORDING_SESSION_LABEL #makes a copy of session label to get SID from
  df_final3$SessionLabelCopy <- gsub("(?i)[a-z]", "", df_final3$SessionLabelCopy)
  df_final3$SessionLabelCopy <- gsub("_", "", df_final3$SessionLabelCopy)
  df_final3$SessionLabelCopy <- substring(df_final3$SessionLabelCopy,1, 2) #selects first2 charachters to remove numbers included for repeated subjects 
  names(df_final3)[names(df_final3) == 'SessionLabelCopy'] <- "SID"
  
  #Removing excluded Subjects####
  #getting rid of 5,6,20, 22 and 31 (uncomment to remove subjects)
  df_32final <- df_final3#[!(df_final3$SID == "05" | df_final3$SID == "06" | df_final3$SID == "20"| df_final3$SID == "22"| df_final3$SID == "31"),]
  #unique(df_32final$SID)
  
  #Sorting by Img Name and SID####
  df_32final <- df_32final[order(df_32final$image_filename),] 
  df_32final <- df_32final[order(df_32final$SID),] 
  
  #Accuracy####
  df_32final$Accuracy <- "placeholder"
  
  #converting to numeric for accuracy detection
  df_32final$test_confidence <-  as.numeric(as.character(df_32final$test_confidence))
  
  df_32final$Accuracy[(df_32final$test_confidence > 3.5 & df_32final$oldnearorfar != "Old" )] <- "FalseAlarm"
  df_32final$Accuracy[(df_32final$test_confidence < 2.5 & df_32final$oldnearorfar != "Old" )] <- "CorrectRejection"
  df_32final$Accuracy[(df_32final$test_confidence > 3.5 & df_32final$oldnearorfar == "Old" )] <- "Hit"
  df_32final$Accuracy[(df_32final$test_confidence < 2.5 & df_32final$oldnearorfar == "Old" )] <- "Miss"
  df_32final$Accuracy[(df_32final$test_confidence == "0")] <- "NoResponse"
  df_32final$Accuracy[(df_32final$test_confidence == "3" & df_32final$oldnearorfar == "Old")] <- "Miss"
  df_32final$Accuracy[(df_32final$test_confidence == "3" & df_32final$oldnearorfar != "Old")] <- "CorrectRejection"
  df_32final$Accuracy[(df_32final$test_confidence == "-1")] <- "StudyPhase"
  df_32final$Accuracy[is.na(df_32final$test_confidence)] <- "NA"
  df_32final$IsGuess[(df_32final$test_confidence == "3")] <- "Guess"
  df_32final$IsGuess[(df_32final$test_confidence != "3")] <- "NotAGuess"
   #Note one trial with "." was converted to NA
  
  #checking trials/subject (sanity check)
table(df_32final$SID)
table(df_32final$RECORDING_SESSION_LABEL)

#subj 3 crashed block 9 and wasn't able to restart at block 10
#subj 24 participant was falling asleep way too much so stopped at block 10
#subj 26 Unclear why 3 extra trials..
#Conclusion is Sub 9 crashed in block 11. Repeat (9b) had blocks 10 to 12. So 10 in 9b is a repeat (remove), 11 was partially seen before (remove both 9 and 9b) and block 12 is good (keep!)
df_32final <- df_32final[!(df_32final$RECORDING_SESSION_LABEL == "ofy09alb" & df_32final$block == 10),]
df_32final <- df_32final[!(df_32final$RECORDING_SESSION_LABEL == "ofy09alb" & df_32final$block == 11),]
df_32final <- df_32final[!(df_32final$RECORDING_SESSION_LABEL == "ofy09al" & df_32final$block == 11),]

  #Making a column with image family from image name
  df_32final$ImageFamily <- substring(df_32final$image_filename, 1, 2)
  #Making a column with image Morph from image name
  df_32final$ImageMorph <- substring(df_32final$image_filename, 4, 5)
  
  #Matching test accuracy to study accuracy
  df_32final$SID_block_image <- paste0(df_32final$SID, "_", df_32final$block, "_", df_32final$RECORDING_SESSION_LABEL, "_", df_32final$ImageFamily, "_", df_32final$ImageMorph)
  #first subsetting test phase
  df_test <-  df_32final[grepl("test", df_32final$task),]
  df_Test <- df_32final[grepl("Test", df_32final$task),]
  df.test <- rbind(df_test, df_Test, fill = TRUE)   
  
  df.shorttest <- df.test[,c("SID_block_image", "Accuracy", "oldnearorfar", "test_confidence")]
  
  #Merging test with full to get new accuracy column
  df.merged <- merge(x= df_32final, y = df.shorttest, by = "SID_block_image", all = TRUE)
  
  #Changing names to get rid of  .y 
  names(df.merged)[names(df.merged)== 'Accuracy.y'] <- "SubsequentAccuracy"
  #names(df.merged)[names(df.merged)=='test_confidence.y'] <-"SubsequentConfidence"
  names(df.merged)[names(df.merged)=='oldnearorfar.x'] <-"oldnearorfar" #this is now labeling the subsequentmorph as oldnearorfar. Double check for this not introducing bugs
  names(df.merged)[names(df.merged)=='oldnearorfar.y'] <- "SubsequentMorph" 
  
  #Making a column with high and low confidence binned just in case that's useful?  
  df.merged$confidenceBin <- NA 
  df.merged$confidenceBin[(df.merged$test_confidence.y == 1 | df.merged$test_confidence.y == 5)] <- "High Confidence"
  df.merged$confidenceBin[(df.merged$test_confidence.y == 2 | df.merged$test_confidence.y == 4)] <- "Low Confidence"
  df.merged$confidenceBin[(df.merged$test_confidence.y == 3)] <- "Guess"
  df.merged$confidenceBin[(df.merged$test_confidence.y == 0)] <- "NoResponse"
  
  #Calculating dPrime####
  df.merged_both <- df.merged[df.merged$oldnearorfar != ".",]
  df.merged_nearonly <- df.merged[df.merged$oldnearorfar != "Far"& df.merged$oldnearorfar != ".",]
  df.merged_faronly <- df.merged[df.merged$oldnearorfar != "Near" & df.merged$oldnearorfar != ".",]
  
  df.counts_both <- df.merged_both %>% group_by(SubsequentAccuracy, SID) %>% summarise(count=n())
  df.counts_nearonly <- df.merged_nearonly %>% group_by(SID, SubsequentAccuracy) %>% summarise(count=n())
  df.counts_faronly <- df.merged_faronly %>% group_by(SubsequentAccuracy, SID) %>% summarise(count=n())
  
  df.widecount_both <- spread(df.counts_both, SubsequentAccuracy, count)
  df.widecount_nearonly <- spread(df.counts_nearonly, SubsequentAccuracy, count)
  df.widecount_faronly <- spread(df.counts_faronly, SubsequentAccuracy, count)
  
  ####combined near and far
  
  #HitRate
  df.widecount_both$HitRate_hit_miss <- df.widecount_both$Hit / (df.widecount_both$Hit + df.widecount_both$Miss)

  #FalseAlarmRate
  df.widecount_both$FARate_FA_CR <- df.widecount_both$FalseAlarm / (df.widecount_both$FalseAlarm + df.widecount_both$CorrectRejection)

  #d prime not counting guesses
  df.widecount_both$dprime_both_noguess <- (qnorm(df.widecount_both$HitRate_hit_miss) - qnorm(df.widecount_both$FARate_FA_CR))
  
  ####d prime only Far
  
  #HitRate
  df.widecount_faronly$HitRate_hit_miss <- df.widecount_faronly$Hit / (df.widecount_faronly$Hit + df.widecount_faronly$Miss)

  #FalseAlarmRate
  df.widecount_faronly$FARate_FA_CR <- df.widecount_faronly$FalseAlarm / (df.widecount_faronly$FalseAlarm + df.widecount_faronly$CorrectRejection)
 
  #d prime not counting guesses
  df.widecount_faronly$dprime_faronly_noguess <- (qnorm(df.widecount_faronly$HitRate_hit_miss) - qnorm(df.widecount_faronly$FARate_FA_CR))

  ####d prime only near
  
  #HitRate
  df.widecount_nearonly$HitRate_hit_miss <- df.widecount_nearonly$Hit / (df.widecount_nearonly$Hit + df.widecount_nearonly$Miss)

  #FalseAlarmRate
  df.widecount_nearonly$FARate_FA_CR <- df.widecount_nearonly$FalseAlarm / (df.widecount_nearonly$FalseAlarm + df.widecount_nearonly$CorrectRejection)

  #d prime not counting guesses
  df.widecount_nearonly$dprime_nearonly_noguess <- (qnorm(df.widecount_nearonly$HitRate_hit_miss) - qnorm(df.widecount_nearonly$FARate_FA_CR))
 
  #combining dprimes to 1 dataframe
  df.justDprimes <- df.widecount_both[,c("SID","dprime_both_noguess")]
  df.justDprimes1 <- df.widecount_faronly[,c("dprime_faronly_noguess")]
  df.justDprimes2 <-  df.widecount_nearonly[,c("dprime_nearonly_noguess")]
  df.justdp <- cbind(df.justDprimes, df.justDprimes1, df.justDprimes2)
  
  #save to csv on google drive
  write.csv(df.justdp, file = "C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/EyeTracking/ObjFamETStudy/OF_ET_Accuracyv1.csv", row.names = FALSE)
  
  
  #merge dprimes to main file just for fun / in case needed
  #selecting only desired columns of df.justdp
  df.justdpshort <- df.justdp[, c('SID','dprime_nearonly_noguess', 'dprime_faronly_noguess')]
  
  #IF want to merge the other ways of looking at dprime, change following code to merge df.justdp instead of df.justdpshort.
  df.merged_wDp <- merge(df.merged, df.justdpshort, by = 'SID')
  
  #Writing out file to be merged with Machine Learning sample data####
  #writing one to ML directory
  write.csv(x=df.merged_wDp, file='C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/EyeTracking/ObjFamETStudy/MACHINE_LEARNING_codeForPradeep/dp_and_subqAcc.csv', row.names = F)
  #writing one to home folder
  write.csv(x=df.merged_wDp, file='C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/EyeTracking/ObjFamETStudy/Kirk Graphs/Clean_forStephCodeTesting/dp_and_subqAcc.csv', row.names = F)


  #Loading in Steph's data####
  # define some paths
  data_dir = 'C:/Users/kgeier/Google Drive/OlsenLab_cloud/All Experiments/EyeTracking/ObjFamETStudy/Kirk Graphs/Clean_forStephCodeTesting' 
  filename =  'data_test_Steph.csv' 
  filename =  '_onsets_FromStephForDPrime.csv'
  filepath = paste(data_dir, '/', filename, sep = ''); filepath
  
  # read in the data
  df.StephTest = read.csv(filepath)
  
  #sorting to just test phase
  df.StephTest <- df.StephTest[df.StephTest$Segment == 'Test',]
  
  #changing names (built code for Steph's data_test_Steph file and now checking on onset file to see if dprime changes)
  names(df.StephTest)[names(df.StephTest) == 'Subject'] <- "subid"
  names(df.StephTest)[names(df.StephTest) == 'Morph_Resp'] <- "condition"
  
  #names(df.StephTest)[names(df.StephTest) == 'Subject'] <- "subid"
  
  #finding dprime for Steph data#### 
  
  # df.StephTest$Accuracy <- ifelse((df.StephTest$condition == "1_new"), "Miss", 
  #                                  ifelse((df.StephTest$condition == "2_new"), "CorrectRejection",
  #                                     ifelse((df.StephTest$condition == "3_new"), "CorrectRejection",
  #                                         ifelse((df.StephTest$condition == "1_old"), "Hit",
  #                                                ifelse((df.StephTest$condition == "2_old"), "FalseAlarm",
  #                                     ifelse((df.StephTest$condition == "3_old"), "FalseAlarm",
  #                                            ifelse((df.StephTest$Resp_bin == "Guess"), "Guess", "error")))))))
  # ###########################################################################################################################################################
   #version with Guess to check if that matches Stephs code better
   df.StephTest$Accuracy <- ifelse((df.StephTest$condition == "1_new"), "Miss",
                                   ifelse((df.StephTest$condition == "2_new"), "CorrectRejection",
                                          ifelse((df.StephTest$condition == "3_new"), "CorrectRejection",
                                                 ifelse((df.StephTest$condition == "1_old"), "Hit",
                                                        ifelse((df.StephTest$condition == "2_old"), "FalseAlarm",
                                                             ifelse((df.StephTest$condition == "3_old"), "FalseAlarm",
                                                                      ifelse((df.StephTest$Resp_bin == "Guess" & df.StephTest$Morph > 1.5), "CorrectRejection", "Miss")))))))
  ###########################################################################################################################################################

  
  df.StephTest_nearonly <- df.StephTest[!grepl(3,df.StephTest$condition),]
  df.StephTest_faronly <-  df.StephTest[!grepl(2,df.StephTest$condition),]
  
  df.SGcounts_both <- df.StephTest %>% group_by(Accuracy, subid) %>% summarise(count=n())
  df.SGcounts_nearonly <- df.StephTest_nearonly %>% group_by(Accuracy, subid) %>% summarise(count=n())
  df.SGcounts_faronly <- df.StephTest_faronly %>% group_by(Accuracy, subid) %>% summarise(count=n())
  
  df.SGwidecount_both <- spread(df.SGcounts_both, Accuracy, count)
  df.SGwidecount_nearonly <- spread(df.SGcounts_nearonly, Accuracy, count)
  df.SGwidecount_faronly <- spread(df.SGcounts_faronly, Accuracy, count)
  
  ####d prime only Far Steph
  
  #HitRate
  df.SGwidecount_faronly$HitRate_hit_miss <- df.SGwidecount_faronly$Hit / (df.SGwidecount_faronly$Hit + df.SGwidecount_faronly$Miss)
  
  #FalseAlarmRate
  df.SGwidecount_faronly$FARate_FA_CR <- df.SGwidecount_faronly$FalseAlarm / (df.SGwidecount_faronly$FalseAlarm + df.SGwidecount_faronly$CorrectRejection)
  
  #d prime not counting guesses
  df.SGwidecount_faronly$dprime_faronly_noguess <- (qnorm(df.SGwidecount_faronly$HitRate_hit_miss) - qnorm(df.SGwidecount_faronly$FARate_FA_CR))
  
  
  ####d prime only near Steph
  
  #HitRate
  df.SGwidecount_nearonly$HitRate_hit_miss <- df.SGwidecount_nearonly$Hit / (df.SGwidecount_nearonly$Hit + df.SGwidecount_nearonly$Miss)
  
  #FalseAlarmRate
  df.SGwidecount_nearonly$FARate_FA_CR <- df.SGwidecount_nearonly$FalseAlarm / (df.SGwidecount_nearonly$FalseAlarm + df.SGwidecount_nearonly$CorrectRejection)
  
  #d prime not counting guesses
  df.SGwidecount_nearonly$dprime_nearonly_noguess <- (qnorm(df.SGwidecount_nearonly$HitRate_hit_miss) - qnorm(df.SGwidecount_nearonly$FARate_FA_CR))
  
  #combining dprimes to 1 dataframe
  df.SGjustDprimes1 <- df.SGwidecount_faronly[,c("subid", "dprime_faronly_noguess")]
  df.SGjustDprimes2 <-  df.SGwidecount_nearonly[,c("dprime_nearonly_noguess")]
  df.SGjustdp <- cbind(df.SGjustDprimes1, df.SGjustDprimes2)
  df.SGjustdpshort <- df.SGjustdp[, c('subid','dprime_nearonly_noguess', 'dprime_faronly_noguess')]
  
    #Graphs for dprime####
  
  #making long format df with SID and dprimes
  df.longDprime <- gather(df.justdpshort, `dprime_nearonly_noguess`, `dprime_faronly_noguess`, key = "Generosity", value = "dprime_Value")
  df.SGlongDprime <- gather(df.SGjustdpshort, `dprime_nearonly_noguess`, `dprime_faronly_noguess`, key = "Generosity", value = "dprime_Value")
  
  #Combining Steph and Kirk dPrimes
  df.SGlongDprime$Generosity <- paste0(df.SGlongDprime$Generosity, "_S")
  df.longDprime$Generosity <- paste0(df.longDprime$Generosity, "_B")
  names(df.SGlongDprime)[names(df.SGlongDprime) == 'subid'] <- "SID"
  df.BlongDprime <- rbind(df.SGlongDprime, df.longDprime)
  #
  
  
  plotDP <- ggplot(df.BlongDprime, aes(x=Generosity, y=dprime_Value)) +
    geom_point(aes(fill=Generosity), size=4, shape=21, colour="grey20",
               position=position_jitter(width=0.6, height=0)) +
    geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
    labs(title="Subject Performance- All Subjects") +
    xlab("Location and Type")+
    ylab("d'") +
    theme(legend.position="none")+ 
    scale_x_discrete(breaks=c("dprime_nearonly_noguess_B","dprime_faronly_noguess_B","dprime_nearonly_noguess_S","dprime_faronly_noguess_S"), labels=c("Baycrest-Near Only", "Baycrest- Far Only", "Stanford- Near Only", "Stanford- Far Only")) 
  plotDP
  

  
  #Filtering out subjects with too low dprime. GOing to be defined as below 0 near and 0.5 far
  

  
  lowdprimesSubj <- df.justdpshort[df.justdpshort$dprime_nearonly_noguess < 0 | df.justdpshort$dprime_faronly_noguess <0.5,]
  SGlowdprimesSubj <- df.SGjustdpshort[df.SGjustdpshort$dprime_nearonly_noguess < 0 | df.SGjustdpshort$dprime_faronly_noguess <0.5,]
  
 
  df.lowdp_remov <- df.justdpshort[!(df.justdpshort$dprime_nearonly_noguess < 0 | df.justdpshort$dprime_faronly_noguess <0.5),]
  df.SGlowdp_remov <- df.SGjustdpshort[!(df.SGjustdpshort$dprime_nearonly_noguess < 0 | df.SGjustdpshort$dprime_faronly_noguess <0.5),]

  #Adding yes to new column indicating 'yes, sufficient dprime to be kept'
  df.lowdp_remov$SufDP <- "y"
  highdp_SIDs <- df.lowdp_remov[,c("SID","SufDP")]
  
  #adding no for rows without sufficient dprime to be kept. 
  lowdprimesSubj$SufDP <- "n"
  lowdp_SIDs <- lowdprimesSubj[,c("SID","SufDP")]
  
  SIDs_sufdp <- rbind(highdp_SIDs, lowdp_SIDs)
  
  df.merged_wDp <- merge(df.merged_wDp, SIDs_sufdp, by = 'SID', all = T)
  
  # EXCLUDING SUBJECTS WITH LOW DPRIME. COMMENT THIS PART OUT IF YOU WOULD LIKE TO LEAVE THEM IN
  df.merged_wDp <- df.merged_wDp[df.merged_wDp$SufDP=="y",]
  
  #OPTIONAL - Uncomment following lines to match dprimes####
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 32)
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 14)
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 34)
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 10)
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 29)
  # #What's above is required to have same median. Below is to ahve same average (note bc near and far are different can't match both perfectly)
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 15)
  # df.merged_wDp <- subset(df.merged_wDp, subset = df.merged_wDp$SID != 26)
  #####################################################################################
  #repeating graph after excluding subjects
  df.BlongDprime2 <- df.BlongDprime[grep("dprime_(near|far)only_noguess_B",df.BlongDprime$Generosity, value = F),]
  df.BlongDprime2 <- merge(x=df.BlongDprime2, y= highdp_SIDs, all = F)
  df.BlongDprime2 <- within(df.BlongDprime2, rm(SufDP))
  #OPTIONAL - Uncomment following lines to match dprimes####
  # #REMOVING MORE SUBJECTS MANUALLY TO SEE WHEN AVERAGE DPRIMES BECOME EQUAL 
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 32)
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 14)
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 34)
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 10)
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 29)
  # #What's above is required to have same median. Below is to ahve same average (note bc near and far are different can't match both perfectly)
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 15)
  # df.BlongDprime2 <- subset(df.BlongDprime2, subset = df.BlongDprime2$SID != 26)
  # 
  df.BlongDprime2 <- rbind(df.BlongDprime2, (df.BlongDprime[grep("dprime_(near|far)only_noguess_S", df.BlongDprime$Generosity, value = F),]))
 
  
  plotDP_exl <- ggplot(df.BlongDprime2, aes(x=Generosity, y=dprime_Value)) +
    geom_point(aes(fill=Generosity), size=4, shape=21, colour="grey20",
               position=position_jitter(width=0.6, height=0)) +
    geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") +
    labs(title="Subject Performance- All Subjects") +
    xlab("Location and Type")+
    ylab("d'") +
    theme(legend.position="none")+ 
    scale_x_discrete(breaks=c("dprime_nearonly_noguess_B","dprime_faronly_noguess_B","dprime_nearonly_noguess_S","dprime_faronly_noguess_S"), labels=c("Baycrest-Near Only", "Baycrest- Far Only", "Stanford- Near Only", "Stanford- Far Only")) 
  plotDP_exl
  
  #checking if new dprimes are the same or not. 
  #df.BlongDprime2
  df.BlongDprime2 <- subset(df.BlongDprime2, subset = is.na(SID) ==F)
  ddply(df.BlongDprime2,~Generosity, summarise,avg=mean(dprime_Value, na.rm = T), sd=sd(dprime_Value, na.rm = T))
 
  just_B <- df.BlongDprime2[order(df.BlongDprime2$dprime_Value),]
  

   ##############################################################################
  
  #Cleaning step 2. Removing 0 fixations then doing Tukey boxplot method.#### 
  
  #Find and remove all trials with 0 fixations. 
  df.no0fix <- df.merged_wDp[df.merged_wDp$FIXATION_COUNT !=0 ,] 
  
  #checking trial counts for each subject with 0 fixation trials removed
  table(df.no0fix$SID)
  
  #Calculating quartiles and median #REMOVE "IP_LABEL" if want to group by suject only not subject and study/test
  df.sum <- ddply(df.no0fix, .(SID, IP_LABEL), function(df.no0fix) quantile(df.no0fix$FIXATION_COUNT, na.rm =T))
   
  #calculating pos(ible) and prob(able) ranges for outliers
  df.sum$BotWhisk_pos <- (df.sum$`25%`-(1.5* (df.sum$`75%`-df.sum$`25%`))) 
  df.sum$TopWhisk_pos <- (df.sum$`75%`+ (1.5* (df.sum$`75%`-df.sum$`25%`)))
  df.sum$BotWhisk_prob <- (df.sum$`25%`-(3 * (df.sum$`75%`-df.sum$`25%`)))
  df.sum$TopWhisk_prob <- (df.sum$`75%`+(3 * (df.sum$`75%`-df.sum$`25%`)))

  #Making new columns and merging boxplot extremes with core data
  df.no0fix$Outlier_pos <- "placeholder"
  df.no0fix$Outlier_prob <- "ph"
  df.no0fix <- merge(df.no0fix, df.sum, by = c("SID","IP_LABEL")) #REMOVE "IP_LABEL" if want to group by suject only not subject and study/test
  
  #marking all trials that have fixations outside of the whiskers as yes or no. 
  df.no0fix$Outlier_pos <- ifelse((df.no0fix$FIXATION_COUNT < df.no0fix$BotWhisk_pos | df.no0fix$FIXATION_COUNT > df.no0fix$TopWhisk_pos), "yes", "no")
  df.no0fix$Outlier_prob <- ifelse((df.no0fix$FIXATION_COUNT < df.no0fix$BotWhisk_prob | df.no0fix$FIXATION_COUNT > df.no0fix$TopWhisk_prob), "yes", "no")
  
  #counting number of trials as eliminated or kept. 
  table(df.no0fix$Outlier_pos)
  table(df.no0fix$Outlier_prob)
  
  #in this case removing all possible outliers based on the theory that you can witness 
  #eye-tracker issues which lead to weird fixation counts. Feel free to adjust. 
  df.clean <- df.no0fix[df.no0fix$Outlier_prob == "no",]  
  
  #Cleaning step 3. 
  #Verifying total dwell time and blink time isn't weird.#### 
  #Want to calculate how much of 2.5s trial eyes are actually open (ex. so blinks dont skew data) but, times are adding up to much more than length of trial..
  df.clean$timeEyesOpen <- (df.clean$FIXATION_COUNT *df.clean$AVERAGE_FIXATION_DURATION)
  df.clean$timeEyesclosed <- (df.clean$BLINK_COUNT * df.clean$AVERAGE_BLINK_DURATION)
  
 
  #
  #Make tables with bin sizes####
  #To remove subjects with fewer than 9 trials in each bin. I.e. if comparing 50 Hits to 6 False Alarms, might be too noisy.

  df.merged_wDp <- df.clean
  
  df.merged_wDp$NewOldResponse[(df.merged_wDp$test_confidence.y > 3.5)] <- 'OLD_resp'
  df.merged_wDp$NewOldResponse[(df.merged_wDp$test_confidence.y < 2.5)] <- 'NEW_resp'
  df.merged_wDp$NewOldResponse[(df.merged_wDp$test_confidence.y == 3 | df.merged_wDp$test_confidence.y == 0)] <- 'GuessOrNo_resp' 
  
  #Grouping (Note, session label 9b is blocks 10-12 of subject 9. Next bit of code groups by SID not Recording session label)
  df.table_binsize <- df.merged_wDp  %>% group_by(NewOldResponse, oldnearorfar, SID) %>%  summarise(count = length(oldnearorfar))
  df.table_binsize <- df.table_binsize[df.table_binsize$oldnearorfar != ".",]
  df.fewer3 <- df.table_binsize[df.table_binsize$count < 2.5 ,]
  
  #steph bin size
  df.SGtable_binsize <- df.StephTest  %>% group_by(condition, subid) %>%  summarise(count = length(condition))
  df.SGfewer3 <- df.SGtable_binsize[df.SGtable_binsize$count < 2.5 ,]
  
  #SO, realizing that the cutoff Steph used was actually 9 study trials not 9 test trials. So changing following code to exclude fewer subjects 
  
  #Eliminate bins with fewer than x(2.5) trials in each bin.
  
  #Matching count to test trials of main df 
  df.markingSmallBins <- merge(df.merged_wDp, df.table_binsize, by = c('SID', 'NewOldResponse', 'oldnearorfar'), all = F)
  #Making a shorter mor managable version of df with just the 2 Relevant columns (so cleaner when merging with main df)
  df.markingSmallRel <- df.markingSmallBins[, c("SID_block_image", "count")]
#Merging the count df to main df by SID_block_image
  df.withCount <- merge(df.markingSmallRel, df.merged_wDp, by = 'SID_block_image', all = T)
  
  #Marking rows with fewer than 9 trials in that bin  
  #NOTE - excluded trials that are to be removed rather than marking them as NA (bc so many columns to mark)
  #ALSO NOTE - excludes trials with small bin size even if not one of the 6 bins (ex. guess)
  df.noSmallBins <- df.withCount[(df.withCount$count > 2.5),]

  #Confidence by Morph Graph####
  df.confByMorph <- df.noSmallBins[df.noSmallBins$test_confidence.y != 0,]
  df.confByMorph <- df.confByMorph[df.confByMorph$oldnearorfar != '.',]
  df.confByMorph <- df.confByMorph[!is.na(df.confByMorph$oldnearorfar),]
  d_sumCxM = ddply(df.confByMorph,~oldnearorfar + SID,
                summarise,TestConfidence=mean(test_confidence.y)) 
  

  detach("package:dplyr", unload=TRUE)

   d_sumCxM = summarySEwithin(d_sumCxM, 
                          measurevar="TestConfidence", 
                          withinvars=c("oldnearorfar"),
                          idvar="SID") 
  
  
  
  pCxM <- ggplot(d_sumCxM, aes(x=factor(oldnearorfar), y = TestConfidence)) +
    geom_bar(stat = "summary", fun.y = "mean") +
    # facet_wrap(~SID) +
    theme_classic() +
    theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
    theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
          axis.title = element_text(size=32),
          strip.text.x = element_text(size=16)) +
    scale_y_continuous(limits= c(0, 5)) +
    theme(legend.position="none") +
    ggtitle("Average Confidence by Morph") +
    xlab("Morph") +
    ylab("Average Confidence")+ 
    scale_x_discrete(limits=c("Old","Near", "Far")) +
    geom_errorbar(aes(x=oldnearorfar, ymin= (TestConfidence-se), ymax=(TestConfidence+se), width = 0.25))
  #pCxM
  
  #Making a graph that plots confidence by morph but faceted for each subject
  d_sumCxMfacet = ddply(df.confByMorph,~oldnearorfar + SID,
                   summarise,TestConfidence=mean(test_confidence.y)) 
  
  pCxM_facet <- ggplot(d_sumCxMfacet, aes(x=oldnearorfar, y =TestConfidence, group = 1, colour = factor(SID))) +
    geom_point(stat = "identity") +
    geom_line() +
    facet_grid(~SID) +
    theme_classic() +
    theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
    theme(plot.title = element_text(size=36), axis.text = element_text(size=10), 
          axis.title = element_text(size=32),
          strip.text.x = element_text(size=16)) +
    scale_y_continuous(limits= c(0.9, 5)) +
    theme(legend.position="none") +
    ggtitle("Average Confidence by Morph") +
    xlab("Morph") +
    ylab("Average Confidence")+ 
    scale_x_discrete(limits=c("Old","Near", "Far"), labels=c("O","N","F"))
  #pCxM_facet
  
  #Checking number of subjects remaining for OF-ET. Only 27. why? 
  length(unique(d_sumCxMfacet$SID))
  

  
  # #lmer to use mixed effect model rather than linear model
  library(lme4)
  library(lmerTest)
  library(dplyr)

#Start of Steph Stats####
  d <- df.noSmallBins
  
  #OPTIONAL - Uncomment following lines to run on only repnum =1 #####
  #subj 1-4 don't have values in repnum (studyphase) column (so for this first pass are eliminated)
  d <- subset(d, subset = d$studyphase == 1)
    #OPTIONAL PART FINISHED
  
  #convert to numeric 
  d$AVERAGE_BLINK_DURATION <- as.numeric(as.character(d$AVERAGE_BLINK_DURATION))
  d$BLINK_COUNT <-as.numeric(as.character(d$BLINK_COUNT))
  d$FIXATION_COUNT <-as.numeric(as.character(d$FIXATION_COUNT))
  d$AVERAGE_FIXATION_DURATION <-as.numeric(as.character(d$AVERAGE_FIXATION_DURATION))
  d$LAST_TRIAL_TIME <-as.numeric(as.character(d$LAST_TRIAL_TIME))
  d$X__Reaction_Time__1 <-as.numeric(as.character(d$X__Reaction_Time__1))
  d$displayTime <-as.numeric(as.character(d$displayTime))
  d$SubqMorph <- ifelse(d$SubsequentMorph == "Old", 1, ifelse(d$SubsequentMorph == "Near", 2, ifelse(d$SubsequentMorph == "Far",3, "error"))) 
  d$SubqMorph <- as.numeric(as.character(d$SubqMorph))
  class(d$FIXATION_COUNT)
  #AHHHHH. There must be a better way...
  
  # split into just study df
  ds = d[d$task == 'study',]
  # str(ds)
  
  #Optional Subset to use just 1 and 5 responses
  #ds2 = ds[ds$Confidence.Subsequent == c(1,5),]
  
  #Renaming from kirk's name to steph's so can use same code
  names(ds)[names(ds) == 'test_confidence.y'] <- 'Conf'
  names(ds)[names(ds) == 'SID'] <- 'Subject.ID'
  
  # missing data?
  dim(ds[ds$Conf == '.',])
  dim(ds[ds$Conf == '0',])
  
  # 0 SubqMorphs? Commented out by Kirk because doesn't work for factors. ANd no 0 morphs
  # ds = ds[ds$SubqMorph > 0,]
  
  # Convert conf to quantitative variable
  conf_responses = c('1', '2', '3', '4', '5')
  ds = ds[ds$Conf %in% conf_responses,]
  ds$Confq = factor(ds$Conf)
  ds$Confq = as.integer(ds$Confq)
  # str(ds)
  
  # Convert conf to categorical
  ds$Confbinary = ds$Confq
  ds[ds$Confq < 3, ]$Confbinary = 'new'
  #OPTIONAL. Uncomment to include guesses as new####
  ds[ds$Confq < 3.5, ]$Confbinary = 'new'
  ds[ds$Confq > 3, ]$Confbinary = 'old'
  ds$Confbinary = factor(ds$Confbinary)
  # str(ds)
  
  # set up some contrasts
  contrasts(ds$Confbinary) = cbind(guessVSother = c(-1,.5,.5),
                                   oldVSnew = c(0,-1,1)); contrasts(ds$Confbinary) 
  
  library(dplyr)
  #df_ktest <- ds %>% group_by(Subject.ID) %>% summarise(N= length(Morph))
  
  
  ## How does subsequent test response (categorical) and SubqMorph type relate to fixations at encoding?
  
  d_sum = ddply(ds,~SubqMorph + Confbinary + Subject.ID,
                summarise,FIXATION_COUNT=mean(FIXATION_COUNT)) 
  
  detach("package:dplyr", unload=TRUE)
  d_sum = summarySEwithin(d_sum[d_sum$Confbinary != '3',], 
                          measurevar="FIXATION_COUNT", 
                          withinvars=c("SubqMorph", "Confbinary"),
                          idvar="Subject.ID") 
  d_sum
  
  ggplot(d_sum, aes(x=SubqMorph, y=FIXATION_COUNT, group=Confbinary, 
                    colour=Confbinary)) +
    geom_line(position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
    geom_point(size=5, position=pd) + 
    scale_color_brewer(palette='Set1')
  
  
  d_sum2 = ddply(ds,~SubqMorph + Confbinary + Subject.ID,
                 summarise,FIXATION_COUNT=mean(FIXATION_COUNT))
  
  res2 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * Confbinary + (1|Subject.ID),
              data=d_sum2)
  summary(res2) #no res1 
  
  # using all data (including item effects)
  d_sum3 = ddply(ds,~SubqMorph + Confbinary + Subject.ID+image_filename,
                 summarise, FIXATION_COUNT=mean(FIXATION_COUNT))
  
  #seems to be random interecept fixed slope (captures variance in overall fixations but has fixed slope for effect of SubqMorph and old/new)
  res3 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * Confbinary + (1|Subject.ID) +
                (1|image_filename), data=d_sum3)
  summary(res3)
  
  
  ## How do subsequent test confidence (ordinal? -kg) and SubqMorph type relate to fixations at encoding?
  
  d_sum4 = ddply(ds,~SubqMorph + Confq + Subject.ID,
                 summarise,FIXATION_COUNT=mean(FIXATION_COUNT))
  
  #converting to factor SubqMorph and Confq in case it helps with ggplot not finding a certain object. 
  d_sum4$SubqMorph <- as.factor((d_sum4$SubqMorph))
  class(d_sum4$SubqMorph)
  d_sum4$Confq <- as.factor((d_sum4$Confq))
  class(d_sum4$Confq)
  
  d_sum4 = summarySEwithin(d_sum4, 
                           measurevar="FIXATION_COUNT", 
                           withinvars=c("SubqMorph", "Confq"),
                           idvar="Subject.ID")
 #Uncomment these lines to see just 1 and 5 (extremes of confidence)
   # d_sum4 <- d_sum4[d_sum4$Confq != 2,]
  # d_sum4 <- d_sum4[d_sum4$Confq != 3,]
  # d_sum4 <- d_sum4[d_sum4$Confq != 4,]
  
  ggplot(d_sum4, aes(x=SubqMorph, y=FIXATION_COUNT, group=Confq, 
                     colour=Confq)) +
    geom_line(position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
    geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum4, position=pd, size=1.5) +
    geom_point(size=5, position=pd) + 
    scale_color_brewer(palette='RdYlBu')
  
  
  d_sum5 = ddply(ds,~SubqMorph + Confq + Subject.ID,
                 summarise, FIXATION_COUNT=mean(FIXATION_COUNT))
  
  d_sum5$Confq <- as.numeric(as.character(d_sum5$Confq))
  d_sum5$SubqMorph <- as.numeric(as.character(d_sum5$SubqMorph))
  d_sum5$Subject.ID <- as.numeric(as.character(d_sum5$Subject.ID))
  d_sum5$FIXATION_COUNT <- as.numeric(as.character(d_sum5$FIXATION_COUNT))
  res4 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * scale(Confq) + (1|Subject.ID),
              data=d_sum5)
  summary(res4)
  
  # using all data (including item effects)
  d_sum6 = ddply(ds,~SubqMorph + Confq + Subject.ID+image_filename,
                 summarise, FIXATION_COUNT=mean(FIXATION_COUNT))
  d_sum6$Confq <- as.numeric(as.character(d_sum6$Confq))
  d_sum6$SubqMorph <- as.numeric(as.character(d_sum6$SubqMorph))
  d_sum6$Subject.ID <- as.numeric(as.character(d_sum6$Subject.ID))
  d_sum6$FIXATION_COUNT <- as.numeric(as.character(d_sum6$FIXATION_COUNT))
  
#subject and image are random intercepts. Fixation count to morph AND confidence is fixed effect  
  res5 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * scale(Confq) + (1|Subject.ID) +
                (1|image_filename), data=d_sum6)
  summary(res5)
  
  
  
  
  # mixed effects model for between subjects across 3 repnum.
  summary.df_eachviewing <- ddply(ds, ~Subject.ID + dprime_faronly_noguess + studyphase, summarise,MeanFixations=mean(FIXATION_COUNT))
  
   res_lmer <- lmer(MeanFixations ~ dprime_faronly_noguess * studyphase + (1|Subject.ID), summary.df_eachviewing)
   summary(res_lmer)

  #function to make linear model and R^2
  lm_eqn <- function(dfp){
  m <- lmer(MeanFixations ~ dPrime  + (1|Subject.ID) + (1|studyphase), data = summary.df_eachviewing);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
           list(a = format(coef(m)[1], digits = 2),
                b = format(coef(m)[2], digits = 2),
               r2 = format(summary(m)$r.squared, digits = 3)))
      as.character(as.expression(eq));
  }

    plotA <- ggplot(summary.df_eachviewing, aes(x=dprime_faronly_noguess, y=MeanFixations, colour = studyphase)) +
    #geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
   # geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
    geom_point(size=3) +
    scale_color_brewer(palette='Blues') +
    geom_smooth(method='lmer',formula=y~x + (1|Subject.ID) + (1|studyphase), data = summary.df_eachviewing)
    #geom_text(x = 0.4, y = 8, label = lm_eqn(dfp), parse = TRUE)

    plotA

    #summary of fixations averaged accross repnum this time.
  summary.df_average <- ddply(ds, ~Subject.ID + dprime_faronly_noguess, summarise, MeanFixations=mean(FIXATION_COUNT),  dPrime = mean(dprime_faronly_noguess))
   
  #mixed effects model for between subjects across 3 repnum.
  res_lmerB <- lmer(dPrime ~ MeanFixations + (1|Subject.ID), summary.df_average)
  summary(res_lmerB)


  #function with df_average
    lm_eqn <- function(dfp){
  m <- lm(MeanFixations ~ dPrime, summary.df_average);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
           list(a = format(coef(m)[1], digits = 2),
                b = format(coef(m)[2], digits = 2),
               r2 = format(summary(m)$r.squared, digits = 3)))
      as.character(as.expression(eq));
  }
  plotB <- ggplot(summary.df_average, aes(x=dPrime, y=MeanFixations)) +
    #geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
   # geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
    geom_point(size=3) +
    scale_color_brewer(palette='Blues') +
    geom_smooth(method='lm',formula=y~x, se = TRUE)+
    geom_text(x = 0.4, y = 8, label = lm_eqn(dfp), parse = TRUE)

    plotB

    #load in plotly for interactive graphs
    library(plotly)
    ggplotly(plotA)
    
    
    
    
    
    
    
    
    
    
    
    
    
    ###########################################################
    #Redoing with "time eyes open" instead of FIXATION count to see if that has a bigger impact####
    
    #Now Redoing with "Average FIXATION Duration" instead of FIXATION count to see if that has a bigger impact####
    #Now redoing with "Average FIXATION Duration" but for only repnum 1
    #Now redoing with all subjects (past 0 and 0.5 dprime cutoff) for all reps, and FIXATION count, but counting guess as new
    
    #OPTIONAL - cntrl F and replace FIXATION_COUNT with average FIXATION duration or whatever else 
    
    #Confidence by Morph Graph####
    df.confByMorph <- df.noSmallBins[df.noSmallBins$test_confidence.y != 0,]
    df.confByMorph <- df.confByMorph[df.confByMorph$oldnearorfar != '.',]
    df.confByMorph <- df.confByMorph[!is.na(df.confByMorph$oldnearorfar),]
    d_sumCxM = ddply(df.confByMorph,~oldnearorfar + SID,
                     summarise,TestConfidence=mean(test_confidence.y)) 
    
    detach("package:dplyr", unload=TRUE)
    
    d_sumCxM = summarySEwithin(d_sumCxM, 
                               measurevar="TestConfidence", 
                               withinvars=c("oldnearorfar"),
                               idvar="SID") 
    
    
    
    pCxM <- ggplot(d_sumCxM, aes(x=factor(oldnearorfar), y = TestConfidence)) +
      geom_bar(stat = "summary", fun.y = "mean") +
      # facet_wrap(~SID) +
      theme_classic() +
      theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
      theme(plot.title = element_text(size=36), axis.text = element_text(size=24), 
            axis.title = element_text(size=32),
            strip.text.x = element_text(size=16)) +
      scale_y_continuous(limits= c(0, 5)) +
      theme(legend.position="none") +
      ggtitle("Average Confidence by Morph") +
      xlab("Morph") +
      ylab("Average Confidence")+ 
      scale_x_discrete(limits=c("Old","Near", "Far")) +
      geom_errorbar(aes(x=oldnearorfar, ymin= (TestConfidence-se), ymax=(TestConfidence+se), width = 0.25))
    #pCxM
    
    #Making a graph that plots confidence by morph but faceted for each subject
    d_sumCxMfacet = ddply(df.confByMorph,~oldnearorfar + SID,
                          summarise,TestConfidence=mean(test_confidence.y)) 
    
    pCxM_facet <- ggplot(d_sumCxMfacet, aes(x=oldnearorfar, y =TestConfidence, group = 1, colour = factor(SID))) +
      geom_point(stat = "identity") +
      geom_line() +
      facet_grid(~SID) +
      theme_classic() +
      theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
      theme(plot.title = element_text(size=36), axis.text = element_text(size=10), 
            axis.title = element_text(size=32),
            strip.text.x = element_text(size=16)) +
      scale_y_continuous(limits= c(0.9, 5)) +
      theme(legend.position="none") +
      ggtitle("Average Confidence by Morph") +
      xlab("Morph") +
      ylab("Average Confidence")+ 
      scale_x_discrete(limits=c("Old","Near", "Far"), labels=c("O","N","F"))
    #pCxM_facet
    
    #Checking number of subjects remaining for OF-ET. Only 27. why? 
    length(unique(d_sumCxMfacet$SID))
    
    
    
    # #lmer to use mixed effect model rather than linear model
    library(lme4)
    library(lmerTest)
    library(dplyr)
    
    #Start of Steph Stats####
    d <- df.noSmallBins
    
    #OPTIONAL - Uncomment following lines to run on only repnum =1 #####
    #subj 1-4 don't have values in repnum (studyphase) column (so for this first pass are eliminated)
    #d <- subset(d, subset = d$studyphase == 1)
    #OPTIONAL PART FINISHED
    
    #convert to numeric 
    d$AVERAGE_BLINK_DURATION <- as.numeric(as.character(d$AVERAGE_BLINK_DURATION))
    d$BLINK_COUNT <-as.numeric(as.character(d$BLINK_COUNT))
    d$FIXATION_COUNT <-as.numeric(as.character(d$FIXATION_COUNT))
    d$FIXATION_COUNT <-as.numeric(as.character(d$FIXATION_COUNT))
    d$LAST_TRIAL_TIME <-as.numeric(as.character(d$LAST_TRIAL_TIME))
    d$X__Reaction_Time__1 <-as.numeric(as.character(d$X__Reaction_Time__1))
    d$displayTime <-as.numeric(as.character(d$displayTime))
    d$SubqMorph <- ifelse(d$SubsequentMorph == "Old", 1, ifelse(d$SubsequentMorph == "Near", 2, ifelse(d$SubsequentMorph == "Far",3, "error"))) 
    d$SubqMorph <- as.numeric(as.character(d$SubqMorph))
    class(d$FIXATION_COUNT)
    #AHHHHH. There must be a better way...
    
    # split into just study df
    ds = d[d$task == 'study',]
    # str(ds)
    
    #Optional Subset to use just 1 and 5 responses
    #ds2 = ds[ds$Confidence.Subsequent == c(1,5),]
    
    #Renaming from kirk's name to steph's so can use same code
    names(ds)[names(ds) == 'test_confidence.y'] <- 'Conf'
    names(ds)[names(ds) == 'SID'] <- 'Subject.ID'
    
    # missing data?
    dim(ds[ds$Conf == '.',])
    dim(ds[ds$Conf == '0',])
    
    # 0 SubqMorphs? Commented out by Kirk because doesn't work for factors. ANd no 0 morphs
    # ds = ds[ds$SubqMorph > 0,]
    
    # Convert conf to quantitative variable
    conf_responses = c('1', '2', '3', '4', '5')
    ds = ds[ds$Conf %in% conf_responses,]
    ds$Confq = factor(ds$Conf)
    ds$Confq = as.integer(ds$Confq)
    # str(ds)
    
    # Convert conf to categorical
    ds$Confbinary = ds$Confq
    ds[ds$Confq < 3, ]$Confbinary = 'new'
    ds[ds$Confq > 3, ]$Confbinary = 'old'
    ds$Confbinary = factor(ds$Confbinary)
    # str(ds)
    
    # set up some contrasts
    contrasts(ds$Confbinary) = cbind(guessVSother = c(-1,.5,.5),
                                     oldVSnew = c(0,-1,1)); contrasts(ds$Confbinary) 
    
    library(dplyr)
    #df_ktest <- ds %>% group_by(Subject.ID) %>% summarise(N= length(Morph))
    
    
    ## How does subsequent test response (categorical) and SubqMorph type relate to FIXATIONs at encoding?
    
    d_sum = ddply(ds,~SubqMorph + Confbinary + Subject.ID,
                  summarise,FIXATION_COUNT=mean(FIXATION_COUNT)) 
    
    detach("package:dplyr", unload=TRUE)
    d_sum = summarySEwithin(d_sum[d_sum$Confbinary != '3',], 
                            measurevar="FIXATION_COUNT", 
                            withinvars=c("SubqMorph", "Confbinary"),
                            idvar="Subject.ID") 
    d_sum
    
    ggplot(d_sum, aes(x=SubqMorph, y=FIXATION_COUNT, group=Confbinary, 
                      colour=Confbinary)) +
      geom_line(position=pd, size=1.5) +
      geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
      geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
      geom_point(size=5, position=pd) + 
      scale_color_brewer(palette='Set1')
    
    
    d_sum2 = ddply(ds,~SubqMorph + Confbinary + Subject.ID,
                   summarise,FIXATION_COUNT=mean(FIXATION_COUNT))
    
    res2 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * Confbinary + (1|Subject.ID),
                data=d_sum2)
    summary(res2) #no res1 
    
    # using all data (including item effects)
    d_sum3 = ddply(ds,~SubqMorph + Confbinary + Subject.ID+image_filename,
                   summarise, FIXATION_COUNT=mean(FIXATION_COUNT))
    
    #seems to be random interecept fixed slope (captures variance in overall FIXATIONs but has fixed slope for effect of SubqMorph and old/new)
    res3 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * Confbinary + (1|Subject.ID) +
                  (1|image_filename), data=d_sum3)
    summary(res3)
    
    
    ## How do subsequent test confidence (ordinal? -kg) and SubqMorph type relate to FIXATIONs at encoding?
    
    d_sum4 = ddply(ds,~SubqMorph + Confq + Subject.ID,
                   summarise,FIXATION_COUNT=mean(FIXATION_COUNT))
    
    #converting to factor SubqMorph and Confq in case it helps with ggplot not finding a certain object. 
    d_sum4$SubqMorph <- as.factor((d_sum4$SubqMorph))
    class(d_sum4$SubqMorph)
    d_sum4$Confq <- as.factor((d_sum4$Confq))
    class(d_sum4$Confq)
    
    d_sum4 = summarySEwithin(d_sum4, 
                             measurevar="FIXATION_COUNT", 
                             withinvars=c("SubqMorph", "Confq"),
                             idvar="Subject.ID")
    #Uncomment these lines to see just 1 and 5 (extremes of confidence)
    # d_sum4 <- d_sum4[d_sum4$Confq != 2,]
    # d_sum4 <- d_sum4[d_sum4$Confq != 3,]
    # d_sum4 <- d_sum4[d_sum4$Confq != 4,]
    
    ggplot(d_sum4, aes(x=SubqMorph, y=FIXATION_COUNT, group=Confq, 
                       colour=Confq)) +
      geom_line(position=pd, size=1.5) +
      geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
      geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum4, position=pd, size=1.5) +
      geom_point(size=5, position=pd) + 
      scale_color_brewer(palette='RdYlBu')
    
    
    d_sum5 = ddply(ds,~SubqMorph + Confq + Subject.ID,
                   summarise, FIXATION_COUNT=mean(FIXATION_COUNT))
    
    d_sum5$Confq <- as.numeric(as.character(d_sum5$Confq))
    d_sum5$SubqMorph <- as.numeric(as.character(d_sum5$SubqMorph))
    d_sum5$Subject.ID <- as.numeric(as.character(d_sum5$Subject.ID))
    d_sum5$FIXATION_COUNT <- as.numeric(as.character(d_sum5$FIXATION_COUNT))
    res4 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * scale(Confq) + (1|Subject.ID),
                data=d_sum5)
    summary(res4)
    
    # using all data (including item effects)
    d_sum6 = ddply(ds,~SubqMorph + Confq + Subject.ID+image_filename,
                   summarise, FIXATION_COUNT=mean(FIXATION_COUNT))
    d_sum6$Confq <- as.numeric(as.character(d_sum6$Confq))
    d_sum6$SubqMorph <- as.numeric(as.character(d_sum6$SubqMorph))
    d_sum6$Subject.ID <- as.numeric(as.character(d_sum6$Subject.ID))
    d_sum6$FIXATION_COUNT <- as.numeric(as.character(d_sum6$FIXATION_COUNT))
    
    res5 = lmer(FIXATION_COUNT ~ scale(SubqMorph) * scale(Confq) + (1|Subject.ID) +
                  (1|image_filename), data=d_sum6)
    summary(res5)
    
    
    
    
    # mixed effects model for between subjects across 3 repnum.
    summary.df_eachviewing <- ddply(ds, ~Subject.ID + dprime_faronly_noguess + studyphase, summarise,MeanFIXATIONs=mean(FIXATION_COUNT))
    
    res_lmer <- lmer(MeanFIXATIONs ~ dprime_faronly_noguess * studyphase + (1|Subject.ID), summary.df_eachviewing)
    summary(res_lmer)
    
    #function to make linear model and R^2
    lm_eqn <- function(dfp){
      m <- lmer(MeanFIXATIONs ~ dPrime  + (1|Subject.ID) + (1|studyphase), data = summary.df_eachviewing);
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                       list(a = format(coef(m)[1], digits = 2),
                            b = format(coef(m)[2], digits = 2),
                            r2 = format(summary(m)$r.squared, digits = 3)))
      as.character(as.expression(eq));
    }
    
    plotA <- ggplot(summary.df_eachviewing, aes(x=dprime_faronly_noguess, y=MeanFIXATIONs, colour = studyphase)) +
      #geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
      # geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
      geom_point(size=3) +
      scale_color_brewer(palette='Blues') +
      geom_smooth(method='lmer',formula= y~x + (1|Subject.ID) + (1|studyphase), data = summary.df_eachviewing)
    #geom_text(x = 0.4, y = 8, label = lm_eqn(dfp), parse = TRUE)
    
    plotA
    
    #summary of FIXATIONs averaged accross repnum this time.
    summary.df_average <- ddply(ds, ~Subject.ID + dprime_faronly_noguess, summarise, MeanFIXATIONs=mean(FIXATION_COUNT),  dPrime = mean(dprime_faronly_noguess))
    
    #mixed effects model for between subjects across 3 repnum.
    res_lmerB <- lmer(dPrime ~ MeanFIXATIONs + (1|Subject.ID), summary.df_average)
    summary(res_lmerB)
    
    
    #function with df_average
    lm_eqn <- function(dfp){
      m <- lm(MeanFIXATIONs ~ dPrime, summary.df_average);
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                       list(a = format(coef(m)[1], digits = 2),
                            b = format(coef(m)[2], digits = 2),
                            r2 = format(summary(m)$r.squared, digits = 3)))
      as.character(as.expression(eq));
    }
    plotB <- ggplot(summary.df_average, aes(x=dPrime, y=MeanFIXATIONs)) +
      #geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), colour="red",position=pd, size=1.5) +
      # geom_errorbar(width=0, aes(ymin=FIXATION_COUNT-ci, ymax=FIXATION_COUNT+ci), data=d_sum, position=pd, size=1.5) +
      geom_point(size=3) +
      scale_color_brewer(palette='Blues') +
      geom_smooth(method='lm',formula=y~x, se = TRUE)+
      geom_text(x = 0.4, y = 8, label = lm_eqn(dfp), parse = TRUE)
    
    plotB
    
