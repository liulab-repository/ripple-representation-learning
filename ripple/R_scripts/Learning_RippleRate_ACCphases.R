library(lme4)
library(emmeans)
library(afex)


fileName <- 'dt_trial_level_Ripple.csv'

 
parentfolder <- file.path(dirname(getwd()),'Results' )  
data_file <- file.path(parentfolder,fileName)
DD <- read.table(data_file,header=TRUE,sep=",")
#
for (phase  in c('stimulus','feedback','interval')){
  D <- DD
  D$ripple <- D[[phase]]
  D <- D[!is.nan(D$ripple),]
  
  stat <- "Learning_ripple_ACC_phase";
  path<-file.path(parentfolder,stat)
  file1<-paste("mixed_effects_model_",phase,".csv",sep="")
  file2<-paste("tvalues_",phase,".csv",sep="")
  file3<-paste("emmeans_ACC_",phase,".csv",sep="")
  file4<-paste("contrast_ACC_",phase,".csv",sep="")
  dir.create(path)
  file1=file.path(path,file1)
  file2=file.path(path,file2)
  file3=file.path(path,file3)
  file4=file.path(path,file4)
  #
  model<-mixed(ripple~ACC+(1|subject/channel),data=D,test_intercept=FALSE)
  write.csv(model$anova_table,file1)
  summary(model$full.model)
  anova(model)
  #
  coef(summary(model))
  write.csv(coef(summary(model)),file2)
  #
  output3<-emmeans(model,pairwise~ACC,adjust="none",mode="kenward-roger")
  write.csv(summary(output3$emmeans),file3)
  write.csv(summary(output3$contrast),file4)
}
#==================================
