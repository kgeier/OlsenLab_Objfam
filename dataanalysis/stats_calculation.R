library("ggplot2")
library(lmerTest)
library("lme4")
library(car)
library(emmeans)
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
  detach("package:plyr", unload=TRUE) 
} 
library(plyr)
library(dplyr)
library(readxl)

# define some paths: removed


#Fixation count####
FCdf = read.csv(Filepath1)
tail(FCdf)

as.factor(FCdf$SubsequentMorph) -> FCdf$SubsequentMorph
FCdf$Accuracy = relevel (FCdf$Accuracy, "Incorrect")
FCdf$SubsequentMorph = factor(FCdf$SubsequentMorph, levels = c("Identical", "Near", "Far"))

#try contrast coding
contrasts(FCdf$SubsequentMorph) <- cbind(c(1, -1463/2773, -1463/2773), c(0, 1, -1410/1363))
colnames(contrasts(FCdf$SubsequentMorph)) <- c('oldvsnew', 'NearvsFar')

#regression model
lmfc3 <- glmer(Accuracy ~  Fix_Count * SubsequentMorph + (1 | Subject.ID) + (1|image_filename),data=FCdf, family=binomial(link="logit"))
summary(lmfc3)

emtrends(lmfc3, pairwise ~ factor(SubsequentMorph), var="Fix_Count")
emmip(ref_grid(lmfc3, transform = "response", cov.reduce = range), SubsequentMorph ~ Fix_Count, CIs = TRUE)

#fixed-effect OR
exp(fixef(lmfc3))
exp(confint(lmfc3, method = "Wald"))

#post-hoc####
#OLD
d_old <- subset(FCdf, SubsequentMorph == 'Identical')
old = glmer(Accuracy~ Fix_Count + (1|Subject.ID) + (1|image_filename), data=d_old, family = "binomial")
old_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=d_old, family = "binomial")
sumOldPH <- anova(old, old_simple)
anova(old, old_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(old)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#New
d_new <- subset(FCdf, SubsequentMorph != 'Identical')
new = glmer(Accuracy~ Fix_Count + (1|Subject.ID) + (1|image_filename), data=d_new, family = "binomial")
new_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=d_new, family = "binomial")
sumOldPH <- anova(new, new_simple)
anova(new, new_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(new)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#Near
d_Near <- subset(FCdf, SubsequentMorph == 'Near')
Near = glmer(Accuracy~ Fix_Count + (1|Subject.ID) + (1|image_filename), data=d_Near, family = "binomial")
Near_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=d_Near, family = "binomial")
sumNearPH <- anova(Near, Near_simple)
anova(Near, Near_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(Near)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#Far
d_Far <- subset(FCdf, SubsequentMorph == 'Far')
Far = glmer(Accuracy~ Fix_Count + (1|Subject.ID) + (1|image_filename), data=d_Far, family = "binomial")
Far_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=d_Far, family = "binomial")
sumFarPH <-   anova(Far, Far_simple)
anova(Far, Far_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(Far)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1


#bias scores check 
FCdf$Confbinary = relevel(FCdf$Confbinary, "new")
lmfc1 = glmer(Confbinary ~ Fix_Count + (1 | Subject.ID) + (1|image_filename),data=FCdf, family="binomial")
summary(lmfc1)



#Sampling variability####
JRdf = read.csv(Filepath2)
head(JRdf)

as.factor(JRdf$SubsequentMorph) -> JRdf$SubsequentMorph
JRdf$Accuracy = relevel (JRdf$Accuracy, "Incorrect")
JRdf$SubsequentMorph = factor(JRdf$SubsequentMorph, levels = c("Identical", "Near", "Far"))
JRdf$scaledArea <- (JRdf$JRobin_areaSampling/1000) 

#try contrast coding
contrasts(JRdf$SubsequentMorph) <- cbind(c(1, -1459/2756, -1459/2756), c(0, 1, -1401/1359))
colnames(contrasts(JRdf$SubsequentMorph)) <- c('oldvsnew', 'NearvsFar')

#regression model
lmJR3 <- glmer(Accuracy ~ scaledArea* SubsequentMorph + (1 | Subject.ID) + (1|image_filename), data=JRdf, family="binomial")
summary(lmJR3)

#fixed-effect OR
exp(fixef(lmJR3))
exp(confint(lmJR3, method = "Wald"))

#post-hoc####
#OLD
dJR_old <- subset(JRdf, SubsequentMorph == "Identical")
old = glmer(Accuracy~ scale(JRobin_areaSampling) + (1|Subject.ID) + (1|image_filename), data=dJR_old, family = "binomial")
old_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dJR_old, family = "binomial")
sumOldPH <- anova(old, old_simple)
anova(old, old_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(old)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#New
dJR_new <- subset(JRdf, SubsequentMorph != "Identical")
new = glmer(Accuracy~ scale(JRobin_areaSampling) + (1|Subject.ID) + (1|image_filename), data=dJR_new, family = "binomial")
new_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dJR_new, family = "binomial")
sumOldPH <- anova(new, new_simple)
anova(new, new_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(new)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#bias scores check 

JRdf$Confbinary = ifelse(JRdf$SubsequentMorph == "Far"|JRdf$Accuracy == "Correct", "New",
                         ifelse(JRdf$SubsequentMorph == "Far"|JRdf$Accuracy == "Incorrect", "Old",
                                ifelse (JRdf$SubsequentMorph == "Near"|JRdf$Accuracy == "Incorrect", "Old",
                                        ifelse(JRdf$SubsequentMorph == "Near"|JRdf$Accuracy == "Correct", "New",
                                               ifelse(JRdf$SubsequentMorph == "Identical"|JRdf$Accuracy == "Correct", "Old",
                                                      ifelse(JRdf$SubsequentMorph == "Identical"|JRdf$Accuracy == "Incorrect", "New", "AHHH"))))))
JRdf$Confbinary = as.factor(JRdf$Confbinary)
lmJR1 = glmer(Confbinary ~ scaledArea + (1 | Subject.ID) + (1|image_filename),data=JRdf, family="binomial")
summary(lmJR1)



#Repetition effect 1-2####
Re12df = read.csv(Filepath3)
head(Re12df)

as.factor(Re12df$SubsequentMorph) -> Re12df$SubsequentMorph
Re12df$Accuracy = relevel (Re12df$Accuracy, "Incorrect")
Re12df$SubsequentMorph = factor(Re12df$SubsequentMorph, levels = c("Identical", "Near", "Far"))

#try contrast coding
contrasts(Re12df$SubsequentMorph) <- cbind(c(1, -1219/2302, -1219/2302), c(0, 1, -1173/1129))
colnames(contrasts(Re12df$SubsequentMorph)) <- c('oldvsnew', 'NearvsFar')

#regression model 
lmRe12.3 <- glmer(Accuracy ~ FixCountDiff * SubsequentMorph + (1 | Subject.ID) + (1|image_filename), data=Re12df, family="binomial")
summary(lmRe12.3)

#fixed-effect OR
exp(fixef(lmRe12.3))
exp(confint(lmRe12.3, method = "Wald"))

#post-hocs####

#OLD
dRe12_old <- subset(Re12df, SubsequentMorph == "Identical")
old = glmer(Accuracy~ scale(FixCountDiff) + (1|Subject.ID) + (1|image_filename), data=dRe12_old, family = "binomial",  nAGQ = 0)
old_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dRe12_old, family = "binomial",  nAGQ = 0)
sumOldPH <- anova(old, old_simple)
anova(old, old_simple)


# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(old)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#New
dRe12_new <- subset(Re12df, SubsequentMorph != "Identical")
new = glmer(Accuracy~ scale(FixCountDiff) + (1|Subject.ID) + (1|image_filename), data=dRe12_new, family = "binomial")
new_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dRe12_new, family = "binomial")
sumOldPH <- anova(new, new_simple)
anova(new, new_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(new)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#Near
dRe12_near <- subset(Re12df, SubsequentMorph == "Near")
near = glmer(Accuracy~ scale(FixCountDiff) + (1|Subject.ID) + (1|image_filename), data=dRe12_near, family = "binomial")
near_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dRe12_near, family = "binomial")
sumOldPH <- anova(near, near_simple)
anova(near, near_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(old)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#Far
dRe12_far <- subset(Re12df, SubsequentMorph == "Far")
far = glmer(Accuracy~ scale(FixCountDiff) + (1|Subject.ID) + (1|image_filename), data=dRe12_far, family = "binomial")
far_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dRe12_far, family = "binomial")
sumOldPH <- anova(far, far_simple)
anova(far, far_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(new)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#bias scores check 
Re12df$Confbinary = ifelse(Re12df$SubsequentMorph == "Far"|Re12df$Accuracy == "Correct", "New",
                           ifelse(Re12df$SubsequentMorph == "Far"|Re12df$Accuracy == "Incorrect", "Old",
                                  ifelse (Re12df$SubsequentMorph == "Near"|Re12df$Accuracy == "Incorrect", "Old",
                                          ifelse(Re12df$SubsequentMorph == "Near"|Re12df$Accuracy == "Correct", "New",
                                                 ifelse(Re12df$SubsequentMorph == "Identical"|Re12df$Accuracy == "Correct", "Old",
                                                        ifelse(Re12df$SubsequentMorph == "Identical"|Re12df$Accuracy == "Incorrect", "New", "AHHH"))))))
Re12df$Confbinary = as.factor(Re12df$Confbinary)
Re12df$Confbinary = relevel(Re12df$Confbinary, "New")

lmRe12.1 = glmer(Confbinary ~  FixCountDiff + (1 | Subject.ID) + (1|image_filename),data=Re12df, family="binomial")
summary(lmRe12.1)


#Repetition effect 1-3####
Re13df = read.csv(Filepath4)
head(Re13df)

as.factor(Re13df$SubsequentMorph) -> Re13df$SubsequentMorph
Re13df$Accuracy = relevel (Re13df$Accuracy, "Incorrect")
Re13df$SubsequentMorph = factor(Re13df$SubsequentMorph, levels = c("Identical", "Near", "Far"))

#try contrast coding
contrasts(Re13df$SubsequentMorph) <- cbind(c(1, -1210/2279, -1210/2279), c(0, 1, -1162/1117))
colnames(contrasts(Re13df$SubsequentMorph)) <- c('oldvsnew', 'NearvsFar')

#regression model
lmRe13.3 <- glmer(Accuracy ~ FixCountDiff * SubsequentMorph + (1 | Subject.ID) + (1|image_filename), data=Re13df, family="binomial")
summary(lmRe13.3)

#fixed-effect OR
exp(fixef(lmRe13.3))
exp(confint(lmRe13.3, method = "Wald"))


#post-hocs####

#OLD
dRe13_old <- subset(Re13df, SubsequentMorph == "Identical")
old = glmer(Accuracy~ scale(FixCountDiff) + (1|Subject.ID) + (1|image_filename), data=dRe13_old, family = "binomial")
old_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dRe13_old, family = "binomial")
sumOldPH <- anova(old, old_simple)
anova(old, old_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(old)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#New
dRe13_new <- subset(Re13df, SubsequentMorph != "Identical")
new = glmer(Accuracy~ scale(FixCountDiff) + (1|Subject.ID) + (1|image_filename), data=dRe13_new, family = "binomial")
new_simple = glmer(Accuracy ~  (1|Subject.ID) + (1|image_filename), data=dRe13_new, family = "binomial")
sumOldPH <- anova(new, new_simple)
anova(new, new_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(new)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1


#check if theres sig difference b/w 1-2 and 1-3
sigall = merge(Re12df, Re13df, by = c("SubsequentMorph", "Accuracy", "Subject.ID", "image_filename"), all=T)
t.test(sigall$FixCountDiff.x, sigall$FixCountDiff.y, paired = T)


#bias scores check
Re13df$Confbinary = ifelse(Re13df$SubsequentMorph == "Far"|Re13df$Accuracy == "Correct", "New",
                           ifelse(Re13df$SubsequentMorph == "Far"|Re13df$Accuracy == "Incorrect", "Old",
                                  ifelse (Re13df$SubsequentMorph == "Near"|Re13df$Accuracy == "Incorrect", "Old",
                                          ifelse(Re13df$SubsequentMorph == "Near"|Re13df$Accuracy == "Correct", "New",
                                                 ifelse(Re13df$SubsequentMorph == "Identical"|Re13df$Accuracy == "Correct", "Old",
                                                        ifelse(Re13df$SubsequentMorph == "Identical"|Re13df$Accuracy == "Incorrect", "New", "AHHH"))))))
Re13df$Confbinary = as.factor(Re13df$Confbinary)

lmRe13.1 = glmer(Confbinary ~  FixCountDiff + (1 | Subject.ID) + (1|image_filename),data=Re13df, family="binomial")
summary(lmRe13.1)


#repetition overall####
Reall = read.csv(Filepath5)
head(Reall)

#Add rep stats here
Reall12 <- subset(Reall, subset = Reall$RepComp == "1 & 2")
Reall13 <- subset(Reall, subset = Reall$RepComp == "1 & 3")
Reall1 <- subset(Reall12, subset = Reall12$repnum ==  1)
Reall2 <- subset(Reall12, subset = Reall12$repnum == 2)
Reall3 <- subset(Reall13, subset = Reall13$repnum == 3)

#check to see whether there is repetition effect
t.test(Reall1$Fix_Count, Reall2$Fix_Count, tail= 2, paired = T)
t.test(Reall1$Fix_Count, Reall3$Fix_Count, tail= 2, paired = T)


#encoding similarity####
SSSdf = read.csv(Filepath6)
head(SSSdf)

as.factor(SSSdf$SubsequentMorph) -> SSSdf$SubsequentMorph
SSSdf$acc= relevel (SSSdf$acc, "Incorrect")
SSSdf$SubsequentMorph = factor(SSSdf$SubsequentMorph, levels = c("Old", "Near", "Far"))

#try contrast coding
contrasts(SSSdf$SubsequentMorph) <- cbind(c(1, -1106/2110, -1106/2110), c(0, 1, -1075/1035))
colnames(contrasts(SSSdf$SubsequentMorph)) <- c('oldvsnew', 'NearvsFar')

#regression model
lmSSS3 <- glmer(acc ~ esd_avg * SubsequentMorph + (1 | Subject) + (1|Subj_Image), data=SSSdf, family="binomial")
summary(lmSSS3)


#fixed-effect OR
exp(fixef(lmSSS3))
exp(confint(lmSSS3, method = "Wald"))

#post-hoc####
#Near
dSSSdf_near <- subset(SSSdf, SubsequentMorph == "Near")
near = glmer(acc~ esd_avg + (1 | Subject) + (1|Subj_Image), data=dSSSdf_near, family = "binomial")
near_simple = glmer(acc ~  (1 | Subject) + (1|Subj_Image), data=dSSSdf_near, family = "binomial")
sumOldPH <- anova(near, near_simple)
anova(near, near_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(near)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#Far
dSSSdf_far <- subset(SSSdf, SubsequentMorph == "Far")
far = glmer(acc~ esd_avg + (1 | Subject) + (1|Subj_Image), data=dSSSdf_far, family = "binomial")
far_simple = glmer(acc ~  (1 | Subject) + (1|Subj_Image), data=dSSSdf_far, family = "binomial")
sumOldPH <- anova(far, far_simple)
anova(far, far_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(far)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#bias scores check
SSSdf$Confbinary = ifelse(SSSdf$SubsequentMorph == "Far"|SSSdf$acc == "Correct", "New",
                          ifelse(SSSdf$SubsequentMorph == "Far"|SSSdf$acc == "Incorrect", "Old",
                                 ifelse (SSSdf$SubsequentMorph == "Near"|SSSdf$acc == "Incorrect", "Old",
                                         ifelse(SSSdf$SubsequentMorph == "Near"|SSSdf$acc == "Correct", "New",
                                                ifelse(SSSdf$SubsequentMorph == "Identical"|SSSdf$acc == "Correct", "Old",
                                                       ifelse(SSSdf$SubsequentMorph == "Identical"|SSSdf$acc == "Incorrect", "New", "AHHH"))))))
SSSdf$Confbinary = as.factor(SSSdf$Confbinary)

lmSSS1 = glmer(Confbinary ~  esd_avg + (1 | Subject) + (1|Subj_Image), data=SSSdf, family="binomial")
summary(lmSSS1)


#study-test similarity####
STSdf = read.csv(Filepath7)
head(STSdf)

as.factor(STSdf$SubsequentMorph) -> STSdf$SubsequentMorph
STSdf$acc= relevel (STSdf$acc, "Incorrect")
STSdf$SubsequentMorph = factor(STSdf$SubsequentMorph, levels = c("Old", "Near", "Far"))

#try contrast coding
contrasts(STSdf$SubsequentMorph) <- cbind(c(1, -1363/2569, -1363/2569), c(0, 1, -1311/1258))
colnames(contrasts(STSdf$SubsequentMorph)) <- c('oldvsnew', 'NearvsFar')

#regression model
lmSTS3 <- glmer(acc ~ eye_sim_diff * SubsequentMorph + (1 | Subject) + (1|Subj_Image), data=STSdf, family="binomial")
summary(lmSTS3)

#fixed-effect OR
exp(fixef(lmSTS3))
exp(confint(lmSTS3, method = "Wald"))

#post-hoc####

#Near
dSTSdf_near <- subset(STSdf, SubsequentMorph == "Near")
near = glmer(acc~ eye_sim_diff + (1 | Subject) + (1|Subj_Image), data=dSTSdf_near, family = "binomial")
near_simple = glmer(acc ~  (1 | Subject) + (1|Subj_Image), data=dSTSdf_near, family = "binomial")
sumOldPH <- anova(near, near_simple)
anova(near, near_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(near)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#Far
dSTSdf_far <- subset(STSdf, SubsequentMorph == "Far")
far = glmer(acc~ eye_sim_diff + (1 | Subject) + (1|Subj_Image), data=dSTSdf_far, family = "binomial")
far_simple = glmer(acc ~  (1 | Subject) + (1|Subj_Image), data=dSTSdf_far, family = "binomial")
sumOldPH <- anova(far, far_simple)
anova(far, far_simple)

# Calculate Odds Ratio and 95% CI and output to a table
s1 <- summary(far)$coefficients ; se1 <- 1.96*(s1[,2])
ll <- s1[,1] - se1 ; ul <- s1[,1] + se1
t1 <- rbind(exp(s1[,1]), exp(ll), exp(ul))
row.names(t1) <- c("OR", "LowerLimit", "UpperLimit")
t1

#bias scores check
STSdf$Confbinary = ifelse(STSdf$SubsequentMorph == "Far"|STSdf$acc == "Correct", "New",
                          ifelse(STSdf$SubsequentMorph == "Far"|STSdf$acc == "Incorrect", "Old",
                                 ifelse (STSdf$SubsequentMorph == "Near"|STSdf$acc == "Incorrect", "Old",
                                         ifelse(STSdf$SubsequentMorph == "Near"|STSdf$acc == "Correct", "New",
                                                ifelse(STSdf$SubsequentMorph == "Identical"|STSdf$acc == "Correct", "Old",
                                                       ifelse(STSdf$SubsequentMorph == "Identical"|STSdf$acc == "Incorrect", "New", "AHHH"))))))
STSdf$Confbinary = as.factor(STSdf$Confbinary)

lmSTS1 = glmer(Confbinary ~  eye_sim_diff + (1 | Subject) + (1|Subj_Image), data=STSdf, family="binomial")
summary(lmSTS1)
