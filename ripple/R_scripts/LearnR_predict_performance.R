
library(robustbase)


parentfolder <- file.path(dirname(getwd()),'Results' , .Platform$file.sep)  

setwd(parentfolder)

D <- read.table('dt_sublevel_LearnR_behav.csv',header=TRUE,sep=",")

p <- file.path(parentfolder,'LearnR_predict_infer')
dir.create(p,showWarnings = FALSE)
setwd(p)


#===================
model <- lmrob(eleMemory_perf~stimulus,data=D,  setting = "KS2014")
summary(model)$coefficients
cef <- summary(model)$coefficients
 
model <- lmrob(eleInfer_perf~stimulus,data=D,  setting = "KS2014")
summary(model)$coefficients
cef <- rbind(cef,summary(model)$coefficients)
 

#===================
model <- lmrob(eleMemory_perf~feedback,data=D,   setting = "KS2014")
summary(model)$coefficients
cef <- rbind(cef,summary(model)$coefficients)

model <- lmrob(eleInfer_perf~feedback,data=D,  setting = "KS2014")
summary(model)$coefficients
cef <- rbind(cef,summary(model)$coefficients)
 

#===================
model <- lmrob(eleMemory_perf~interval,data=D,   setting = "KS2014")
summary(model)$coefficients
cef <- rbind(cef,summary(model)$coefficients)
 
model <- lmrob(eleInfer_perf~interval,data=D,   setting = "KS2014")
summary(model)$coefficients
cef <- rbind(cef,summary(model)$coefficients)

 
#==================

H <- ifelse(cef[,4]<0.05,"Sig","")
cef <- cbind(cef,H)
out <- subset(cef, rownames(cef) != "(Intercept)")
write.csv(out,'learnR_predict_feature_performance.csv')
