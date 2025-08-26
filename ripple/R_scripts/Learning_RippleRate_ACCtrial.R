library(lme4)
library(interactions)
library(effects)
library(afex)

fileName <- 'dt_trial_level_Ripple.csv'
parentfolder <- file.path(dirname(getwd()),'Results' )  
path <- file.path(parentfolder, 'Learning_Ripple_Trial')
dir.create(path, showWarnings = FALSE)
setwd(path)

data_file <- file.path(parentfolder, fileName)
DD <- read.table(data_file, header = TRUE, sep = ",")

variables_to_analyze <- c('stimulus', 'feedback', 'interval')

for (var in variables_to_analyze) {
  D <- DD
  D$RippleRate <- D[[var]]
  D  <- D[!is.nan(D$RippleRate),]
  print(paste("Analyzing:", var))
  
  # 
  model <- mixed(RippleRate ~ ACC * trial + (1 | subject/channel), D, test_intercept = FALSE)
  write.csv(model$anova_table,paste0('mixed_model_F_', var, '.csv'))
  
  # mixed model
  model <- lmer(RippleRate ~ ACC * trial + (1 | subject/channel), D)
  write.csv(coef(summary(model)), paste0('mixed_model_', var, '.csv'))
  print(paste("Model coefficients for", var, ":"))
  print(coef(summary(model)))
  
  # interaction
  sim_slopes_results <- sim_slopes(model, pred = "trial", modx = "ACC", johnson_neyman = FALSE, cond.int = TRUE)
  write.csv(sim_slopes_results$ints, paste0('interaction_', var, '_ACCtrial_int.csv'))
  write.csv(sim_slopes_results$slopes, paste0('interaction_', var, '_ACCtrial_slope.csv'))
  print(paste("Simulated slopes for", var, ":"))
  print(sim_slopes_results$slopes)
  
  # effect
  inter.sd <- effect(c('ACC*trial'), mod = model, xlevels = list(ACC = c(1, 0), trial = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)))
  inter.sd <- as.data.frame(inter.sd)
  write.csv(inter.sd, paste0('interaction_', var, '_ACCtrial_fit.csv'))
}


