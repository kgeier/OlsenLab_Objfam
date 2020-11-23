library("lmerTest")
library("lme4")
library(car)
library(ggplot2)
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
  detach("package:plyr", unload=TRUE) 
} 
library(plyr)
library(dplyr)
library(purrr)
library(ggsignif)
library(tidyr)
library(data.table)
library(ggsci)
library(grid)
library(gridExtra)


wid = 200
heig = 140

# define some paths: removed


#Fixation Count graph####
FCdf = read.csv(filepath1) 
head(FCdf)

a = subset(FCdf, AccuracyType2 == "FarCorrect", select=c(Subject.ID, Fix_Count))
b = subset(FCdf, AccuracyType2 == "FarIncorrect", select=c(Subject.ID, Fix_Count))
c = subset(FCdf, AccuracyType2 == "NearCorrect", select=c(Subject.ID, Fix_Count))
d = subset(FCdf, AccuracyType2 == "NearIncorrect", select=c(Subject.ID, Fix_Count))
e = subset(FCdf, AccuracyType2 == "IdenticalCorrect", select=c(Subject.ID, Fix_Count))
f = subset(FCdf, AccuracyType2 == "IdenticalIncorrect", select=c(Subject.ID, Fix_Count))

#Far   
far =  bind_rows(b,a) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))   

#Near       
near =  bind_rows(d,c) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))        

#Identical
identical =  bind_rows(f,e) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))  


combined = list(far,near,identical) %>% reduce(full_join, by = "Subject.ID")
combined

library(data.table)
setnames(combined, c("Fix_Count.x", "Fix_Count.y", "Fix_Count"), c("Far", "Near", "Identical"))

combined <- combined %>% gather(Morphtype, FCD, c("Far", "Near", "Identical"))
combined %>% arrange(desc(FCD))

####Trying to make new graph type
FCnew <- 
  ggplot(combined, aes(x=Morphtype, fill = Morphtype, y= FCD)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size=3.5, width = 0.15,  alpha = 0.5) +
  xlim("Identical", "Near", "Far") +
  ylim(-1.3,2.05) +
  geom_signif(comparisons = list(c(2,3), c(1,2.5)), annotations = c("***","*"), textsize = 11, vjust=0.5, step_increase=.1, linetype=1, alpha=1, color = "black") + 
  geom_text(x=1, y=1.45, label = "*", size=11, color = "black") + 
  ggtitle("Effect of Fixation Count on Accuracy") +
  xlab("Morph Type") +
  ylab("Fix # Difference (Cor. - Incor.)")  +
  stat_summary(fun.y=mean, geom = "point", shape=23, size=5,fill= "white", color = "black") +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.3, size=20), axis.text = element_text(size=15), 
        axis.title = element_text(size=17),
        strip.text.x = element_text(size=16)) +
  geom_vline(xintercept=2.5, linetype=3, alpha=.5) +
  geom_vline(xintercept=1.5, linetype="dashed") + 
  geom_hline(yintercept=c(0,0), linetype=4) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#66c2a5",   "#8da0cb",  "#fc8d62"))

FCnew

filename <- "removed"  
ggsave(path= data.save, plot = FCnew, filename , width = wid, height = heig, units = "mm"  )


#JRareasampling graph####
JRdf = read.csv(filepath2)
head(JRdf)

a = subset(JRdf, AccuracyType2 == "FarCorrect", select=c(Subject.ID, JRobin_areaSampling))
b = subset(JRdf, AccuracyType2 == "FarIncorrect", select=c(Subject.ID, JRobin_areaSampling))
c = subset(JRdf, AccuracyType2 == "NearCorrect", select=c(Subject.ID, JRobin_areaSampling))
d = subset(JRdf, AccuracyType2 == "NearIncorrect", select=c(Subject.ID, JRobin_areaSampling))
e = subset(JRdf, AccuracyType2 == "IdenticalCorrect", select=c(Subject.ID, JRobin_areaSampling))                                                                           
f = subset(JRdf, AccuracyType2 == "IdenticalIncorrect", select=c(Subject.ID, JRobin_areaSampling))

#Far   
far =  bind_rows(b,a) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))   

#Near       
near =  bind_rows(d,c) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))        

#Identical
identical =  bind_rows(f,e) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))  


combined = list(far,near,identical) %>% reduce(full_join, by = "Subject.ID")
combined

setnames(combined, c("JRobin_areaSampling.x", "JRobin_areaSampling.y", "JRobin_areaSampling"), c("Far", "Near", "Identical"))

combined <- combined %>% gather(Morphtype, SV, c("Far", "Near", "Identical"))

JRnew <- 
  ggplot(combined, aes(x=Morphtype, fill =Morphtype, y= SV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3.5, width = 0.15,  alpha = 0.5) +
  xlim("Identical", "Near", "Far") +
  ylim(-2300, 3200) + 
  ggtitle(label = "Effect of SV on Accuracy",
          subtitle = "") + 
  geom_signif(comparisons = list(c(2,3), c(1,2.5)), annotations = c("???", "*"), vjust=.5, step_increase = .1, textsize = 11, linetype=1, alpha=1, color = "black") + 
  xlab("Morph Type") +
  ylab("SV Difference (Cor. - Incor.)")  +
  stat_summary(fun.y=mean, geom = "point", shape=23, size=5,fill= "white", color = "black") +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20), axis.text = element_text(size=15), 
        axis.title = element_text(size=17),
        strip.text.x = element_text(size=16)) +
  geom_vline(xintercept=2.5, linetype=3, alpha=.5) +
  geom_vline(xintercept=1.5, linetype="dashed") + 
  geom_hline(yintercept=c(0 , 0), linetype=6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#66c2a5",   "#8da0cb",  "#fc8d62"))

JRnew

#combine FC and SV
library(cowplot)

FCandSV.combined <- plot_grid(FCnew, JRnew, 
  ncol = 2, nrow = 1, labels = c("a)", "b)"), scale=0.9
) 

FCandSV.combined


#Repeff 1to2 ####
Re12df = read.csv(filepath3)
head(Re12df)
sample_n(Re12df, 20)
a = subset(Re12df, AccuracyType2 == "FarCorrect", select=c(Subject.ID, MagRep))
b = subset(Re12df, AccuracyType2 == "FarIncorrect", select=c(Subject.ID, MagRep))
c = subset(Re12df, AccuracyType2 == "NearCorrect", select=c(Subject.ID, MagRep))
d = subset(Re12df, AccuracyType2 == "NearIncorrect", select=c(Subject.ID, MagRep))
e = subset(Re12df, AccuracyType2 == "IdenticalCorrect", select=c(Subject.ID, MagRep))
f = subset(Re12df, AccuracyType2 == "IdenticalIncorrect", select=c(Subject.ID, MagRep))

#Far   
far =  bind_rows(b,a) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))   

#Near       
near =  bind_rows(d,c) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))        

#Identical
identical =  bind_rows(f,e) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))  


combined = list(far,near,identical) %>% reduce(full_join, by = "Subject.ID")
combined

setnames(combined, c("MagRep.x", "MagRep.y", "MagRep"), c("Far", "Near", "Identical"))

combined <- combined %>% gather(Morphtype, FCD, c("Far", "Near", "Identical"))

Rep12new = 
  ggplot(combined, aes(x=Morphtype, fill =Morphtype, y= FCD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3.5, width = 0.15,  alpha = 0.5) +
  xlim("Identical", "Near", "Far") +
  ylim(-2.9, 2.8)+
  geom_signif(comparisons = list(c(2,3), c(1,2.5)), annotations = c("", "*"), textsize = 11, vjust=.5, step_increase = .03, linetype=1, alpha=1, color = "black")+
  geom_text(x=1, y=2.14, label = "*", size=11, color = "black") + 
  ggtitle(label = "Repetition Effect on Accuracy",
          subtitle = "Repetition 1 to 2") +
  xlab("Morph Type") +
  ylab("Mag of Rep Difference (Cor. - Incor.)")  +
  stat_summary(fun.y=mean, geom = "point", shape=23, size=5,fill= "white", color = "black") +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20), axis.text = element_text(size=17), 
        plot.subtitle = element_text(hjust = 0.5, size=18), 
        axis.title = element_text(size=17),
        strip.text.x = element_text(size=16)) +
  geom_vline(xintercept=2.5, linetype=3, alpha=.5) +
  geom_vline(xintercept=1.5, linetype="dashed") + 
  geom_hline(yintercept=c(0 , 0), linetype=6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#66c2a5",   "#8da0cb",  "#fc8d62"))

Rep12new


#Rep 1-3####
Re13df = read.csv(filepath4)
head(Re13df)
sample_n(Re13df, 20)
a = subset(Re13df, AccuracyType2 == "FarCorrect", select=c(Subject.ID, MagRep))
b = subset(Re13df, AccuracyType2 == "FarIncorrect", select=c(Subject.ID, MagRep))
c = subset(Re13df, AccuracyType2 == "NearCorrect", select=c(Subject.ID, MagRep))
d = subset(Re13df, AccuracyType2 == "NearIncorrect", select=c(Subject.ID, MagRep))
e = subset(Re13df, AccuracyType2 == "IdenticalCorrect", select=c(Subject.ID, MagRep))
f = subset(Re13df, AccuracyType2 == "IdenticalIncorrect", select=c(Subject.ID, MagRep))

#Far   
far =  bind_rows(b,a) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))   

#Near       
near =  bind_rows(d,c) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))        

#Identical
identical =  bind_rows(f,e) %>%
  group_by(Subject.ID) %>%
  summarise_each(funs(diff(.)))  


combined = list(far,near,identical) %>% reduce(full_join, by = "Subject.ID")
combined

setnames(combined, c("MagRep.x", "MagRep.y", "MagRep"), c("Far", "Near", "Identical"))

combined <- combined %>% gather(Morphtype, FCD, c("Far", "Near", "Identical"))

Rep13new = 
  ggplot(combined, aes(x=Morphtype, fill =Morphtype, y= FCD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3.5, width = 0.15,  alpha = 0.5) +
  xlim("Identical", "Near", "Far") +
  ylim(-4.2, 2.5) +
  ggtitle(label = "Repetition Effect on Accuracy",
    subtitle = "Repetition 1 to 3") +
  xlab("Morph Type") +
  ylab("Mag of Rep Difference (Cor. - Incor.)")  +
  stat_summary(fun.y=mean, geom = "point", shape=23, size=5,fill= "white", color = "black") +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20), axis.text = element_text(size=17), 
        plot.subtitle = element_text(hjust = 0.5, size=18), 
        axis.title = element_text(size=17),
        strip.text.x = element_text(size=16)) +
  geom_vline(xintercept=2.5, linetype=3, alpha=.5) +
  geom_vline(xintercept=1.5, linetype="dashed") + 
  geom_hline(yintercept=c(0 , 0), linetype=6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#66c2a5",   "#8da0cb",  "#fc8d62"))

Rep13new


#combine Rep12 and Rep13
library(cowplot)

Rep12and13.combined <- plot_grid(Rep12new + theme(#axis.text.x = element_blank(),
                                                  axis.ticks.x = element_blank(),
                                                  axis.title = element_blank(),
                                                  plot.title = element_blank()), 
                                 Rep13new + theme(plot.title = element_blank(),
                                                  axis.ticks.x = element_blank(),
                                                  axis.title = element_blank()
                                                  ), 
                                 nrow = 2,labels = c("a)", "b)")
                                 ) 
Rep12and13.combined 



#repoverall####
d_rp = read.csv(filepath5)
head(d_rp)
d_rp$Rep= as.factor(d_rp$Rep)
d_rp$repnum = as.factor(d_rp$repnum)

repall = ggplot(data = d_rp, aes(x = RepComp, y = Fix_Count, fill = repnum)) +
  geom_boxplot(width = 0.4) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.15)) +
  theme_bw() +
  ggtitle("Repetition Effect") +
  xlab("Repetition Comparison")+ 
  ylab("Mean Fixation Count") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_jco() + 
  theme(axis.line = element_line(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank()) + 
  theme(legend.position="right")  +
  stat_summary(fun=mean,  geom="point", aes(group=repnum), position=position_dodge(.4), shape=23, size=4,fill= "white", color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 23, face="bold")) +
  theme(axis.title = element_text(size=18)) +
  theme(axis.text = element_text(size=18)) +
  ylim(2,8.2) +
  geom_segment(aes(x=0.9,y=8,xend=1.1, yend=8), linetype = 1, color = "black") +
  geom_segment(aes(x=1.9,y=8,xend=2.1, yend=8), linetype = 1, color = "black") +
  geom_text(x=1, y = 8.2, label = "*", size = 11, color = "black") +
  geom_text(x=2, y = 8.2, label = "†", size = 11, color = "black") 


repall


#StudyTestSimilarity  Graph####
STSdf = read.csv(filepath7)
head(STSdf)
sample_n(STSdf, 20)
STSdf$AccuracyType = ifelse(STSdf$Accuracy == "CorrectRejection"|STSdf$Accuracy == "Hit", "Correct", "Incorrect")
STSdf$AccuracyType2 = paste0(STSdf$SubsequentMorph, STSdf$AccuracyType)
a = subset(STSdf, AccuracyType2 == "FarCorrect", select=c(Subject, eye_sim_diff))
b = subset(STSdf, AccuracyType2 == "FarIncorrect", select=c(Subject, eye_sim_diff))
c = subset(STSdf, AccuracyType2 == "NearCorrect", select=c(Subject, eye_sim_diff))
d = subset(STSdf, AccuracyType2 == "NearIncorrect", select=c(Subject, eye_sim_diff))
e = subset(STSdf, AccuracyType2 == "OldCorrect", select=c(Subject, eye_sim_diff))
f = subset(STSdf, AccuracyType2 == "OldIncorrect", select=c(Subject, eye_sim_diff))

#Far   
far =  bind_rows(b,a) %>%
  group_by(Subject) %>%
  summarise_each(funs(diff(.)))   

#Near       
near =  bind_rows(d,c) %>%
  group_by(Subject) %>%
  summarise_each(funs(diff(.)))        

#Identical
old = merge(e,f, by="Subject") 
old$eye_sim_diff = old$eye_sim_diff.x - old$eye_sim_diff.y
old = old[, !(colnames(old) %in% c("eye_sim_diff.x","eye_sim_diff.y"))]

combined = list(far,near,old) %>% reduce(full_join, by = "Subject")
combined

setnames(combined, c("eye_sim_diff.x", "eye_sim_diff.y", "eye_sim_diff"), c("Far", "Near", "Identical"))

combined <- combined %>% gather(Morphtype, FCD, c("Far", "Near", "Identical"))

studytestnew = 
  ggplot(combined, aes(x=Morphtype, fill =Morphtype, y= FCD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3.5, width = 0.15,  alpha = 0.5) +
  xlim("Identical", "Near", "Far") +
  ylim(-0.35, .2) + 
  geom_signif(comparisons = list(c("Near","Far")), annotations = "†", y_position = .18,textsize = 11, vjust=.3, linetype=1, alpha=1, color = "black") + 
  ggtitle("Effect of Study-Test Similarity on Accuracy") +
  xlab("Morph Type") +
  ylab("Unique Similarity Difference (Cor. - Incor.)")  +
  stat_summary(fun.y=mean, geom = "point", shape=23, size=5,fill= "white", color = "black") +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20), axis.text = element_text(size=17), 
        axis.title = element_text(size=17),
        strip.text.x = element_text(size=16)) +
  geom_vline(xintercept=2.5, linetype=3, alpha=.5) +
  geom_vline(xintercept=1.5, linetype="dashed") +     
  geom_hline(yintercept=c(0 , 0), linetype=6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#66c2a5",   "#8da0cb",  "#fc8d62"))

studytestnew


#encoding similarity####
SSSdf = read.csv(filepath6)
head(SSSdf)
sample_n(SSSdf, 20)

SSSdf$AccuracyType = ifelse(SSSdf$Accuracy == "CorrectRejection"|SSSdf$Accuracy == "Hit", "Correct", "Incorrect")
SSSdf$AccuracyType2 = paste0(SSSdf$SubsequentMorph, SSSdf$AccuracyType)
a = subset(SSSdf, AccuracyType2 == "FarCorrect", select=c(Subject, esd_avg))
b = subset(SSSdf, AccuracyType2 == "FarIncorrect", select=c(Subject, esd_avg))
c = subset(SSSdf, AccuracyType2 == "NearCorrect", select=c(Subject, esd_avg))
d = subset(SSSdf, AccuracyType2 == "NearIncorrect", select=c(Subject, esd_avg))
e = subset(SSSdf, AccuracyType2 == "OldCorrect", select=c(Subject, esd_avg))
f = subset(SSSdf, AccuracyType2 == "OldIncorrect", select=c(Subject, esd_avg))

#Far   
far =  bind_rows(b,a) %>%
  group_by(Subject) %>%
  summarise_each(funs(diff(.)))   

#Near       
near =  bind_rows(d,c) %>%
  group_by(Subject) %>%
  summarise_each(funs(diff(.)))        

#Identical
old = merge(e,f, by="Subject") 
old$esd_avg = old$esd_avg.x - old$esd_avg.y
old = old[, !(colnames(old) %in% c("esd_avg.x","esd_avg.y"))]

combined = list(far,near,old) %>% reduce(full_join, by = "Subject")
combined

setnames(combined, c("esd_avg.x", "esd_avg.y", "esd_avg"), c("Far", "Near", "Identical"))

combined <- combined %>% gather(Morphtype, FCD, c("Far", "Near", "Identical"))

encodingnew = 
  ggplot(combined, aes(x=Morphtype, fill =Morphtype, y= FCD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3.5, width = 0.15,  alpha = 0.5) +
  xlim("Identical", "Near", "Far") +
  geom_signif(comparisons = list(c("Near","Far")), annotations = "???", textsize = 11, vjust=.5, linetype=1, alpha=1, color = "black")+
  ggtitle("Effect of Encoding Similarity on Accuracy") +
  xlab("Morph Type") +
  ylab("Unique Similarity Difference (Cor. - Incor.)")  +
  stat_summary(fun.y=mean, geom = "point", shape=23, size=5,fill= "white", color = "black") +
  theme_classic() +
  theme(axis.line.x= element_line(colour = "black"), axis.line.y = element_line(color= "black")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size=20), axis.text = element_text(size=17), 
        axis.title = element_text(size=17),
        strip.text.x = element_text(size=16)) +
  geom_vline(xintercept=2.5, linetype=3, alpha=.5) +
  geom_vline(xintercept=1.5, linetype="dashed") +  
  geom_hline(yintercept=c(0 , 0), linetype=6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#66c2a5",   "#8da0cb",  "#fc8d62"))

encodingnew



#combine SSS and STS
encodingstudytestnew.combined <- plot_grid(encodingnew + theme(#axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title = element_blank()), 
  studytestnew + theme(axis.ticks.x = element_blank(),
                   axis.title = element_blank()
  ), 
  nrow = 2,labels = c("a)", "b)")
) 
encodingstudytestnew.combined 




