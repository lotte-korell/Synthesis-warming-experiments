##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##              Total abundance at local and gamma scale 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## set seeds where generated in R version 4.4.1, numbers may differ if you use a different version 
## (including the packages)
rm(list = ls())
gc()

## set random number generator (RNG), so it is always the same across different systems
RNGkind("L'Ecuyer-CMRG")

## load packages and install if required
packages <- c("lme4", "MASS", "car", "MuMIn", "arm", "boot", "glmmTMB", "here", "dplyr")

lapply(packages, function(pkg) 
{
  if (!requireNamespace(pkg, quietly = T)) 
    install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##              1. Functions----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to get bootstrap estimates of model and the ci
boot_ci <- function(model_data, seed = NULL) 
{
  if(!is.null(seed)) set.seed(seed) ## set a seed so we always get the same result
  est_boot <- bootMer(model_data, FUN = fixef, nsim = 1000) ## bootstrap the model and calc effect size
  est_ci <- boot.ci(est_boot, index = 1, type = "perc") ## calculate the ci's of the bootstrap
  return(list(est_boot = est_boot, est_ci = est_ci))
}

results_diagnostics <- function(model)
{
  results <- summary(model) ## calculate model results
  rsquare <- r.squaredGLMM(model) ## calculate the rsquare value
  
  ## plots with input to go through
  ## Store current settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))  ## Restore original settings when function exits
  
  ## Ask user before showing each plot
  par(ask = TRUE)
  invisible(qqnorm(resid(model))) ## plot the qqplot, invisibale silences the consol output
  invisible(hist(resid(model))) ## show residuals in a histogram, silence consol output
  
  return(list(results = results, rsquare = rsquare))
}

## read data, using here for cross machine compatibility and it looks in 
## Rproject folder
csv_read <- c(Div.data_local = "Temp_data_final.csv", Div.data_gamma = "Temp_data_gamma_final.csv")

## read data
data_list <- Map(function(data)
  {read.csv(here(data), dec = ".", sep = ";", h = T)}, csv_read)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      2.  Grand mean----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2.1 Calculate models----
## loop to calculate every model needed
mod_list <- c() 
for(i in 1:length(data_list))
{
  name <- names(data_list)[i]
  data <- data_list[[i]]
  
  if(name == "Div.data_local")
  {
    mod_list$grand_mean.S_T <- lmer(LRR_S ~ 1 + (1|study:site:SiteBlock), data = data)
    mod_list$grand_mean.betaSn_T <- lmer(LRR_betaSn ~ 1 + (1|study:site:SiteBlock),  data = data)
    mod_list$grand_mean.SPie_T <- lmer(LRR_SPie ~ 1 + (1|study:site:SiteBlock),  data = data)  
    mod_list$grand_mean.betaSPie_T <- lmer(LRR_betaSPie ~ 1 + (1|study:site:SiteBlock),  data = data) 
    mod_list$grand_mean.beta.sim_T <- lmer(beta.sim ~ 1 + (1|study:site:SiteBlock), data = data)
    mod_list$grand_mean.beta.sne_T <- glmmTMB(beta.sne ~ 1 + (1|study:site:SiteBlock), ziformula=~1,
                                              beta_family(link = "logit"), data = data, na.action = "na.omit")
    mod_list$grand_mean.beta.bray.bal_T <- lmer(beta.bray.bal ~ 1 + (1|study:site:SiteBlock),  data = data) 
    mod_list$grand_mean.beta.bray.gra_T <-glmmTMB(beta.bray.gra ~ 1 + (1|study:site), ziformula=~1, 
            beta_family(link = "logit"), data = data, na.action = "na.omit") 
    
    mod_list$m1.mod <- lmer(LRR_S ~ delta.T1_rescale + duration_rescale + T_warm_rescale + 
                              (1|study:site:SiteBlock), REML = T, data = data, na.action ="na.fail")
    mod_list$m3.mod <- lmer(LRR_betaSn ~ delta.T1_rescale + (1|study:site:SiteBlock),
                            REML = T, data = data, na.action = "na.fail")
    mod_list$m3a.mod <- lmer(LRR_betaSn ~ delta.T1_rescale * T_warm_rescale + (1|study:site:SiteBlock),
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m2.mod <- lmer(LRR_SPie ~ delta.T1_rescale + duration_rescale +
                     (1|study:site:SiteBlock), REML = T, data = data, na.action="na.fail")
    mod_list$m4.mod3 <- lmer(LRR_betaSPie ~ delta.T1_rescale + T_warm_rescale + duration_rescale +
                    (1|study:site:SiteBlock), REML = T, data = data, na.action = "na.fail")
    mod_list$m8.mod <- lmer(beta.sim ~ delta.T1_rescale + T_warm_rescale + (1|study:site:SiteBlock),
                            REML = T, data = data, na.action = "na.fail")
    mod_list$m8a.mod <- lmer(beta.sne ~ delta.T1_rescale + (1|study:site:SiteBlock), 
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m8b.mod <- lmer(beta.bray.bal ~ delta.T1_rescale + (1|study:site:SiteBlock),
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m8c.mod <- lmer(beta.bray.gra ~ delta.T1_rescale + (1|study:site), 
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m1.mod_dur <- lmer(LRR_S ~ duration_rescale + T_warm_rescale + delta.T1_rescale + 
                                  (1|study:site:SiteBlock), REML = T, 
                                data = data, na.action = "na.fail")
    mod_list$m3.mod_dur <- lmer(LRR_betaSn ~ duration_rescale + (1|study:site:SiteBlock),
                                REML = T, data = data, na.action = "na.fail")
    mod_list$m3a.mod_dur <- lmer(LRR_betaSn ~ delta.T1_rescale * T_warm_rescale + duration_rescale +
                                 (1|study:site:SiteBlock), REML = T, data = data, 
                                 na.action = "na.fail")
    mod_list$m2.mod_dur <- lmer(LRR_SPie ~ duration_rescale + (1|study:site:SiteBlock), 
                                REML = T, data = data, na.action = "na.fail")
    mod_list$m4.mod3_dur <- lmer(LRR_betaSPie ~ duration_rescale + T_warm_rescale + (1|study:site:SiteBlock),
                                 REML = T, data = data, na.action = "na.fail")
    mod_list$m9.mod_dur <- lmer(beta.sim ~ duration_rescale + (1|study:site:SiteBlock), REML = T,
                                data = data, na.action = "na.fail")
    mod_list$m9a.mod_dur <- lmer(beta.sne ~ duration_rescale + (1|study:site:SiteBlock),
                                 REML = T, data = data, na.action = "na.fail")
    mod_list$m9b.mod_dur <- lmer(beta.bray.bal ~ duration_rescale + (1|study:site:SiteBlock),
                                 REML = T, data = data, na.action = "na.fail")
    mod_list$m9c.mod_dur <- lmer(beta.bray.gra ~ duration_rescale + (1|study:site),
                                 REML = T, data = data, na.action = "na.fail")
    mod_list$m1.mod1 <- lmer(LRR_S~T_warm_rescale + duration_rescale + delta.T1_rescale + 
                               (1|study:site:SiteBlock), REML = T, data = data, na.action = "na.fail")
    mod_list$m3.mod1 <- lmer(LRR_betaSn ~ T_warm_rescale + (1|study:site:SiteBlock),
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m3a.mod1 <- lmer(LRR_betaSn ~ delta.T1_rescale * T_warm_rescale + (1|study:site:SiteBlock),
                              REML = T, data = data, na.action = "na.fail")
    mod_list$m2.mod1 <- lmer(LRR_SPie ~ T_warm_rescale + delta.T1_rescale + (1|study:site:SiteBlock),
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m4.mod <- lmer(LRR_betaSPie ~ T_warm_rescale + duration_rescale + (1|study:site:SiteBlock),
                            REML = T, data = data, na.action = "na.fail")
    mod_list$m10.mod <- lmer(beta.sim ~ T_warm_rescale + delta.T1_rescale + (1|study:site:SiteBlock),
                             REML = T, data = data, na.action = "na.fail")
    mod_list$m10a.mod <- lmer(beta.sne ~ T_warm_rescale + (1|study:site:SiteBlock), REML = T, 
                     data = data, na.action = "na.fail")
    mod_list$m10b.mod <- lmer(beta.bray.bal ~ T_warm_rescale + (1|study:site:SiteBlock),
                              REML = T, data = data, na.action = "na.fail")
    mod_list$m10c.mod <- lmer(beta.bray.gra ~ T_warm_rescale + (1|study:site),
                              REML = T, data = data, na.action = "na.fail")
  }

  if(name == "Div.data_gamma")
  {
    mod_list$grand_mean.gammaSn_T <- lmer(LRR_gammaSn ~ 1 + (1|study), data = data)
    mod_list$grand_mean.gammaSPie_T <- lmer(LRR_gammaSPie ~ 1 + (1|study), data = data)
    mod_list$m7.mod <- lmer(LRR_gammaSn ~ delta.T1_rescale + (1|study), REML = T, 
                            data = data, na.action = "na.fail")
    mod_list$m6.mod <- lmer(LRR_gammaSPie ~ delta.T1_rescale + T_warm_rescale + (1|study),
                            REML = T, data = data, na.action = "na.fail") 
    mod_list$m5.mod_Dur <- lmer(LRR_gammaSn ~ duration_rescale + delta.T1_rescale + (1|study),
                                REML = T, data = data, na.action = "na.fail")
    mod_list$m7.mod_dur <- lmer(LRR_gammaSn ~ duration_rescale + T_warm_rescale + (1|study),
                                REML = T, data = data, na.action = "na.fail")
    mod_list$m6.mod_dur <- lmer(LRR_gammaSPie ~ duration_rescale + T_warm_rescale + (1|study),
                                REML = T, data = data, na.action = "na.fail")
    mod_list$m7.mod1 <- lmer(LRR_gammaSn ~ T_warm_rescale + (1|study), REML = T, data = data,
                             na.action = "na.fail")
    mod_list$m6.mod <- lmer(LRR_gammaSPie ~ T_warm_rescale + (1|study),
                            REML = T, data = data, na.action = "na.fail")
  }
}

## Model results and diagnostics. Prints two plots (qqplot and histogram), you have to press enter to see both 
results_diagnostics(mod_list$grand_mean.S_T) ## please cycle through the models from the list as you like

## Species richness at local, turnover and gamma scale
m0c <- boot_ci(mod_list$grand_mean.S_T, seed = 1234) %>% print() 
m0e <- boot_ci(mod_list$grand_mean.betaSn_T, seed = 1234) %>% print()
m0g <- boot_ci(mod_list$grand_mean.gammaSn_T, seed = 1234) %>% print()

## Evenness at local, turnover and gamma scale 
m0h <- boot_ci(mod_list$grand_mean.SPie_T, seed = 1234) %>% print()
m0i <- boot_ci(mod_list$grand_mean.betaSPie_T, seed = 1234) %>% print()
m0j <- boot_ci(mod_list$grand_mean.gammaSPie_T, seed = 1234) %>% print()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      Compositional turnover----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Incidence-based pair-wise dissimilarities (Sorensen)
m0i <- boot_ci(mod_list$grand_mean.beta.sim_T, seed = 1234) %>% print()

## Nestedness 
pred <- predict(mod_list$grand_mean.beta.sne_T, type = "response", se.fit = TRUE)

mean_beta.sne <- pred$fit
se_beta.sne <- pred$se.fit

upper <- pred$fit + 1.96 * pred$se.fit
lower <- pred$fit - 1.96 * pred$se.fit

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      Abundance-based pair-wise dissimilarities (Bray-Curtis)----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Turnover
m0k <- boot_ci(mod_list$grand_mean.beta.bray.bal_T, seed = 1234) %>% print()

#nestedness
summary(mod_list$grand_mean.beta.bray.gra_T)
pred_bray.gra <- predict(mod_list$grand_mean.beta.bray.gra_T, type="response", se.fit =TRUE )

mean_beta.sne <- mean(pred_bray.gra$fit)
se_beta.sne <- mean(pred_bray.gra$se.fit)

upper <- mean(pred_bray.gra$fit + 1.96 * pred_bray.gra$se.fit)
lower <- mean(pred_bray.gra$fit - 1.96 * pred_bray.gra$se.fit)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    3. Delta T----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## get a list without the grand means as we did those before
mod_list_sub <- mod_list[-which(names(mod_list) %in% c("grand_mean.S_T","grand_mean.betaSn_T", "grand_mean.SPie_T",         
                                     "grand_mean.betaSPie_T", "grand_mean.beta.sim_T", "grand_mean.beta.sne_T",
                                     "grand_mean.beta.bray.bal_T", "grand_mean.beta.bray.gra_T",
                                     "grand_mean.gammaSn_T", "grand_mean.gammaSPie_T"))]

## create a list with all summaries, coefficients and the bootstrapped CI's 
m_list <- vector(mode = "list", length = length(mod_list_sub)) ## initialize list (needs as many entries
                                                               ## as there will be)

## this will take longer!
for(i in 1:length(mod_list_sub))
{
  ## create names used to name list entries
  level <- sapply(strsplit(names(mod_list_sub)[i], "\\."), '[', 1) ## get the model name
  summary_name <- paste0(level, ".summary") ## create the name for summary entry in the list
  coeffic_name <- paste0(level, ".coeffic") ## create the name for coefficient entry
  ci_name <- paste0(level, ".ci.boot") ## create the name for the 
  
  m_list[[i]][[summary_name]] <- summary(mod_list_sub[[i]]) 
  m_list[[i]][[coeffic_name]] <- summary(mod_list_sub[[i]])$coefficients 
  m_list[[i]][[ci_name]] <- boot_ci(mod_list_sub[[i]], seed = 1234)
  
  names(m_list)[i] <- names(mod_list_sub)[i]
}
## You can just call the results from the list as follow:
m_list$m1.mod$m1.summary ## this calls the summary of the model m1
m_list$m1.mod$m1.coeffic ## this calls the coefficients of model m1 
m_list$m1.mod$m1.ci.boot ## calls the results of the bootstrapped ci's for model m1


