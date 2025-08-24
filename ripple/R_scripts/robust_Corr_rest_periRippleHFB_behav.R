 
#install.packages("robustbase")
 
library(robustbase)
 
 
WkDir= file.path(dirname(getwd()),'Results/Peri_Ripple_HFB_behav', .Platform$file.sep)  

setwd(WkDir)
WkDir
DR = read.table('dt_periRippleHFB_behave.csv',header=TRUE,sep=",")




x_vars = c( 'Y7_Default','Y7_DorsalAttention',	'Y7_Frontoparietal',	'Y7_Limbic',	'Y7_Somatomotor','Y7_VentralAttention',	'Y7_Visual','DMNTL',	'DMNTPJ',	'DMNDLPFC',	'DMNLOFC',	'DMNmPFC',	'DMNPCC')
######################
 
results <- list()

templateList <- data.frame(
  y_var = character(),
  x_var = character(),
  t_value = numeric(),
  df = numeric(),
  cor = numeric(),  
  p_value = numeric(),
  p_flag = character(),
  stringsAsFactors = FALSE
)

tmp_res_xy  = data.frame(
  res_x = numeric(),
  res_y = numeric(),
  y_var = character(),
  x_var = character()
)

res_xy  =  tmp_res_xy
summary_results = templateList

for (y_var  in c('compNonAdjacent','compAdjacent')){
  D = DR
  D$y = D[,y_var]
  D$z = D$testAfter
  
  
  for (x_var in x_vars) {
    # 
    temp_data <- D[, c(x_var, "y", "z")]
    colnames(temp_data) <- c("x", "y", "z")
    #temp_data = temp_data[temp_data$x>-999,]
    temp_data = na.omit(temp_data)
    
    # 
    robust_model_x <- lmrob(x ~ z, data = temp_data,   setting = "KS2014")
    robust_model_y <- lmrob(y ~ z, data = temp_data,   setting = "KS2014")
    
    # 
    res_x <- residuals(robust_model_x)
    res_y <- residuals(robust_model_y)
    
    # 
    partial_r <- cor.test(res_x, res_y)
    
    # 
    results[[x_var]] <- partial_r
    
    if (partial_r$p.value<0.05){
      p_flag = '*'
    }else{
      p_flag = ' '
    }
    
    
    res_xy  = rbind(res_xy, data.frame(
      res_x = res_x,
      res_y = res_y,
      y_var = rep(y_var, length(res_x)),
      x_var = rep(x_var, length(res_x))       
    ))    
    
    summary_results <- rbind(summary_results, data.frame(
      y_var = y_var,
      x_var = x_var,
      t_value = partial_r$statistic,
      df = partial_r$parameter,
      cor = partial_r$estimate,
      p_value = partial_r$p.value,
      p_flag = p_flag
    ))
    
  }
}


for (x_var in x_vars) {
  cat("Results for", x_var, ":\n")
  print(results[[x_var]])
  cat("\n")
}

summary_results1 = summary_results;
print(summary_results)
write.csv(summary_results, file = "R_robCorr_periHFB_CompInfer.csv", row.names = FALSE)
write.csv(res_xy, file = "R_robCorr_periHFB_CompInfer_xy.csv", row.names = FALSE)


 
 

x_vars = c( 'Y7_Default', 'DMNmPFC')

results <-  list()
summary_results = templateList
res_xy  =  tmp_res_xy

for (y_var  in c('compNonAdjacent','compAdjacent')){
  D = DR
  D$y = D[,y_var]
  D$z1 = D$testAfter
  D$z2 = D$mean_Rest
  D$z3 = D$mean_Overall
  for (x_var in x_vars) {

    zx = paste('g_',x_var,sep="")
    D$zx = D[,zx]
    temp_data <- D[, c(x_var, "y", "z1","z2","z3","zx")]
    colnames(temp_data) <- c("x", "y", "z1","z2","z3","zx")
    temp_data = temp_data[temp_data$x>-999,]
    

    robust_model_x <- lmrob(x ~ z1+z2+zx, data = temp_data,   setting = "KS2014")
    robust_model_y <- lmrob(y ~ z1+z2+zx, data = temp_data,   setting = "KS2014")
    

    res_x <- residuals(robust_model_x)
    res_y <- residuals(robust_model_y)
    

    partial_r <- cor.test(res_x, res_y)
    

    results[[x_var]] <- partial_r

    
    if (partial_r$p.value<0.05){
      p_flag = '*'
    }else{
      p_flag = ' '
    }
    
    res_xy  = rbind(res_xy, data.frame(
      res_x = res_x,
      res_y = res_y,
      y_var = rep(y_var, length(res_x)),
      x_var = rep(x_var, length(res_x))       
    ))        
    
    summary_results <- rbind(summary_results, data.frame(
      y_var = y_var,
      x_var = x_var,
      t_value = partial_r$statistic,
      df = partial_r$parameter,
      cor = partial_r$estimate,
      p_value = partial_r$p.value,
      p_flag = p_flag
    ))
    
    
  }
}


for (x_var in x_vars) {
  cat("Results for", x_var, ":\n")
  print(results[[x_var]])
  cat("\n")
}

print(summary_results)
write.csv(summary_results, file = "R_robCorr_periHFB_CompInfer_additionalControll.csv", row.names = FALSE)
write.csv(res_xy, file = "R_robCorr_periHFB_CompInfer_additionalControll_xy.csv", row.names = FALSE)



print(summary_results1)


