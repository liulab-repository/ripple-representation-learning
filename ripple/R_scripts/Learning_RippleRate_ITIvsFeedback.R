 

library(lme4)
library(emmeans)
library(afex)
 

fileName <- 'dt_ripple_feedback_interval.csv'
stat <- "Learning_ripple_iti_vs_feedback";


parentfolder <- file.path(dirname(getwd()),'Results' )  
data_file <- file.path(parentfolder,fileName)
D <- read.table(data_file,header=TRUE,sep=",")
#

path <- file.path(parentfolder,stat)
file1 <- paste("mixed_effects_model_",stat,".csv",sep="")
file2 <- paste("tvalues_",stat,".csv",sep="")
file3 <- paste("emmeans_phase_",stat,".csv",sep="")
file4 <- paste("contrast_phase_",stat,".csv",sep="")
dir.create(path,showWarnings = FALSE)
file1 <- file.path(path,file1)
file2 <- file.path(path,file2)
file3 <- file.path(path,file3)
file4 <- file.path(path,file4)

#
model<-mixed(ripple~phase+ACC+phase*ACC+(1|subject/channel),data=D)
coef(summary(model))
write.csv(model$anova_table,file1)
summary(model$full_model)
anova(model)

#
coef(summary(model))
coef(summary(model$full_model))
write.csv(coef(summary(model)),file2)

#
output3 <- emmeans(model,pairwise~phase,adjust="None",mode="kenward-roger")
write.csv(summary(output3$emmeans),file3)
write.csv(summary(output3$contrast),file4)

#==================================
file5 <- paste("contrast_phase_by_ACC_",stat,".csv",sep="")
file5 <- file.path(path,file5)
file6 <- paste("contrast_emmeans_phase_by_ACC_",stat,".csv",sep="")
file6 <- file.path(path,file6)
emm_i1 <- emmeans(model, "phase", by = c("ACC"))
output <- update(pairs(emm_i1), by = NULL, adjust = "holm")
write.csv(summary(output),file5)
write.csv(summary(emm_i1),file6)
file5 <- paste("contrast_ACC_by_phase_",stat,".csv",sep="")
file5 <- file.path(path,file5)
file6 <- paste("contrast_emmeans_ACC_by_phase_",stat,".csv",sep="")
file6 <- file.path(path,file6)
emm_i1 <- emmeans(model, "ACC", by = c("phase"))
output <- update(pairs(emm_i1), by = NULL, adjust = "holm")
write.csv(summary(output),file5)
write.csv(summary(emm_i1),file6)

