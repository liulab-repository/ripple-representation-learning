
fileName <- 'dt_rest_ripple_feature_behave.csv'
parentfolder <- file.path(dirname(getwd()),'Results')  
path <- file.path(parentfolder, 'Rest_Ripple_feature')
dir.create(path, showWarnings = FALSE)
setwd(path)

data_file <- file.path(parentfolder, fileName)
D <- read.table(data_file, header = TRUE, sep = ",")

 
D$enhance_testMemory <- D$testAfterMemory - D$testBeforeMemory
D$enhance_testInfer <- D$testAfterInfer - D$testBeforeInfer
D$enhance_testAll <- D$testAfter - D$testBefore

library(robustbase)

model <- lmrob(enhance_testAll ~ RR_Rest_enhance, data =D,   setting = "KS2014")
s <- summary(model)
s$coefficients
write.csv(s$coefficients,'lmrob_testAll.csv')

model <- lmrob(enhance_testMemory ~ RR_Rest_enhance, data =D,   setting = "KS2014")
s <- summary(model)
write.csv(s$coefficients,'lmrob_testMemory.csv')

model <- lmrob(enhance_testInfer ~ RR_Rest_enhance, data =D,   setting = "KS2014")
s <- summary(model)
write.csv(s$coefficients,'lmrob_testInfer.csv')


