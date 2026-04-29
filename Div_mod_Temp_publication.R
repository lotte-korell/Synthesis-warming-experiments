rm(list=ls())
gc()

## Set random number generator (RNG), so it is always the same across different systems
RNGkind("L'Ecuyer-CMRG")

## Functions
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
  invisible(qqnorm(resid(model)))## plot the qqplot, invisibale silences the consol output
  invisible(qqline(resid(model)))## Show a line 
  invisible(hist(resid(model))) ## show residuals in a histogram, silence consol output
  
  return(list(results = results, rsquare = rsquare))
}

## Load packages
packages <- c("lme4", "MASS", "car", "MuMIn", "arm", "boot", "here")

lapply(packages, function(pkg) 
{
  if (!requireNamespace(pkg, quietly = T)) 
    install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})

## read data from project location
csv_read <- c( Div.data_local = "Temp_data_final.csv", Div.data_gamma = "Temp_data_gamma_final.csv")

## read data
data_list <- Map(function(data)
{read.csv(here(data), dec = ".", sep = ";", h = T)}, csv_read)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                Species richness 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Full model and dredge for selection
d1.beta <- lmer(LRR_S ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
     REML = F, data = data_list$Div.data_local, na.action="na.fail") %>% 
  dredge() 

model.avg(d1.beta) %>% summary() ## get the model average, only one model has delta <2

## Best model
m1.best1 <- lmer(LRR_S ~ T_warm_rescale + delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock), 
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m1.best1, test = "Chisq")

#results_diagnostics(m1.best1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                ENS PIE     
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d2.beta <- lmer(LRR_SPie ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>% dredge()

model.avg(d2.beta, subset = delta < 2) %>% summary()

## Best model
m2.best <- lmer(LRR_SPie ~ duration_rescale + (1|study:site:SiteBlock),
              REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m2.best, test = "Chisq")
r.squaredGLMM(m2.best) 

#results_diagnostics(m2.best)

## Second best model
m2.best2 <- lmer(LRR_SPie ~ delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
               REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m2.best2)
r.squaredGLMM(m2.best2)
#results_diagnostics(m2.best2)

## Third best model
m2.best3 <- lmer(LRR_SPie ~ T_warm_rescale + duration_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m2.best3)
r.squaredGLMM(m2.best3)
#results_diagnostics(m2.best3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                  Turnover 
##                  Species richness
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Rarefied richness
d3a.beta <- lmer(LRR_betaSn ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                 REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>%
  dredge()


model.avg(d3a.beta, subset = delta < 2) %>% summary(d3a.AIC)

## Best Model
m3a.best1 <- lmer(LRR_betaSn~+(1|study:site:SiteBlock),
                  REML = T,data = data_list$Div.data_local, na.action = "na.fail")

r.squaredGLMM(m3a.best1) 
#results_diagnostics(m3a.best1)

## Second best model
m3a.best2 <- lmer(LRR_betaSn ~ T_warm_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m3a.best2)
r.squaredGLMM(m3a.best2) 
#results_diagnostics(m3a.best2)

## Third best model
m3a.best3 <- lmer(LRR_betaSn ~ delta.T1_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m3a.best3)
r.squaredGLMM(m3a.best3) 
#results_diagnostics(m3a.best3)

## Fourth best model
m3a.best4 <- lmer(LRR_betaSn ~ duration_rescale + (1|study:site:SiteBlock),
                  REML = T, data = Div.data_local, na.action = "na.fail")

Anova(m3a.best4)
r.squaredGLMM(m3a.best4) 
#results_diagnostics(m3a.best4)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##               ENS PIE
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d4.beta <- lmer(LRR_betaSPie ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>% dredge()

model.avg(d4.beta, subset = delta < 2) %>% summary(d4.AIC)

## Best model 
m4.best <- lmer(LRR_betaSPie ~ duration_rescale + T_warm_rescale + (1|study:site:SiteBlock),
                REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m4.best)
r.squaredGLMM(m4.best)
#results_diagnostics(m4.best)

## Second best model
m4.best1 <- lmer(LRR_betaSPie ~ T_warm_rescale + delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m4.best1)
r.squaredGLMM(m4.best1) 
#results_diagnostics(m4.best1)

## Third best model
m4.best2 <- lmer(LRR_betaSPie ~ duration_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m4.best2)
r.squaredGLMM(m4.best2) 
#results_diagnostics(m4.best2)

## Fourth best model
m4.best3 <- lmer(LRR_betaSPie ~ delta.T1_rescale + (1|study:site:SiteBlock),
               REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m4.best3)
r.squaredGLMM(m4.best3) 
#results_diagnostics(m4.best3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                 Compositional turnover 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Incidence-based pair-wise dissimilarities (Sorensen)
## Turnover
d5a.beta <- lmer(beta.sim ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                 REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>% dredge()

model.avg(d5a.beta, subset = delta < 2) %>% summary(d5a.AIC)

## Best model 
m5a.best <- lmer(beta.sim ~ duration_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5a.best)
r.squaredGLMM(m5a.best)
#results_diagnostics(m5a.best)

## Second best model
m5a.best1 <- lmer(beta.sim ~ duration_rescale + T_warm_rescale +(1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5a.best1)
r.squaredGLMM(m5a.best1)
#results_diagnostics(m5a.best1)

## Thirds best model
m5a.best2 <- lmer(beta.sim ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5a.best2)
r.squaredGLMM(m5a.best2)
#results_diagnostics(m5a.best2)

## Nestedness 
d5b.beta <- lmer(beta.sne ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                 REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>% 
  dredge()

model.avg(d5b.beta, subset = delta < 2) %>% summary(d5b.AIC)

## Best model
m5b.best <- lmer(beta.sne ~ T_warm_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5b.best)
r.squaredGLMM(m5b.best)
#results_diagnostics(m5b.best)

## Second best model
m5b.best1 <- lmer(beta.sne ~ delta.T1_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5b.best1)
r.squaredGLMM(m5b.best1)
#results_diagnostics(m5b.best1)

## Third best model
m5b.best2 <-lmer(beta.sne ~ duration_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5b.best2)
r.squaredGLMM(m5b.best2)
#results_diagnostics(m5b.best2)

## Abundance-based pair-wise dissimilarities (Bray-Curtis)
## Turnover
d6a.beta <- lmer(beta.bray.bal ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                 REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>%
  dredge()

model.avg(d6a.beta, subset = delta < 2) %>% summary()

## Best model 
m5c.best <- lmer(beta.bray.bal ~ duration_rescale + (1|study:site:SiteBlock),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5c.best)
r.squaredGLMM(m5c.best)
#results_diagnostics(m5c.best)

## Second best model
m5c.best1 <- lmer(beta.bray.bal ~ T_warm_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5c.best1)
r.squaredGLMM(m5c.best1)
#results_diagnostics(m5c.best1)

## Third best model
m5c.best2 <- lmer(beta.bray.bal ~ delta.T1_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5c.best2)
r.squaredGLMM(m5c.best2)
#results_diagnostics(m5c.best2)

## Fourth best model
m5c.best3 <- lmer(beta.bray.bal ~ delta.T1_rescale + duration_rescale + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m5c.best3)
r.squaredGLMM(m5c.best3)
#results_diagnostics(m5c.best3)

## Nestedness
d6b.beta <- lmer(beta.bray.gra ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study:site),
                 REML = F, data = data_list$Div.data_local, na.action = "na.fail") %>% 
  dredge()## Random effect needed to be reduced 

model.avg(d6b.beta, subset = delta < 2) %>% summary()

## Best model
m6b.best <- lmer(beta.bray.gra ~ duration_rescale + (1|study:site),
                 REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m6b.best)
r.squaredGLMM(m6b.best)
#results_diagnostics(m6b.best)

## Second best model
m6b.best1 <- lmer(beta.bray.gra ~ duration_rescale + T_warm_rescale + (1|study:site),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m6b.best1)
r.squaredGLMM(m6b.best1)
#results_diagnostics(m6b.best1)

## Third best model
m6b.best2 <- lmer(beta.bray.gra ~ duration_rescale + delta.T1_rescale + (1|study:site),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m6b.best2)
r.squaredGLMM(m6b.best2)
#results_diagnostics(m6b.best2)

## Fourth best model
m6b.best3 <- lmer(beta.bray.gra ~ delta.T1_rescale +(1|study:site), 
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

Anova(m6b.best3)
r.squaredGLMM(m6b.best3)
#results_diagnostics(m6b.best3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Gamma
##                        Species richness  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Rarefied richness
d7a.beta <- lmer(LRR_gammaSn ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study),
                 REML = F, data = data_list$Div.data_gamma, na.action = "na.fail") %>% 
  dredge()

model.avg(d7a.beta, subset = delta < 2) %>% summary()

## Best model 
m7a.best <- lmer(LRR_gammaSn ~ delta.T1_rescale + (1|study),
                 REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m7a.best)
r.squaredGLMM(m7a.best) 
#results_diagnostics(m7a.best)

## Second best model
m7a.best1 <- lmer(LRR_gammaSn ~ delta.T1_rescale + duration_rescale + (1|study),
                  REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m7a.best1)
r.squaredGLMM(m7a.best1) 
#results_diagnostics(m7a.best1)

## Third best model
m7a.best2 <- lmer(LRR_gammaSn ~ duration_rescale + (1|study),
                  REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m7a.best2)
r.squaredGLMM(m7a.best2)
#results_diagnostics(m7a.best2)

## Fourth best model
m7a.best3 <- lmer(LRR_gammaSn ~ delta.T1_rescale + duration_rescale + T_warm_rescale + (1|study),
                  REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m7a.best3)
r.squaredGLMM(m7a.best3) 
#results_diagnostics(m7a.best3)

## Fifth best model
m7a.best4 <- lmer(LRR_gammaSn ~ 1 + (1|study), 
                  REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

r.squaredGLMM(m7a.best4) 
#results_diagnostics(m7a.best4)

## Sixth best model
m7a.best5 <- lmer(LRR_gammaSn ~ duration_rescale + T_warm_rescale + (1|study),
                  REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m7a.best5)
r.squaredGLMM(m7a.best5) 
#results_diagnostics(m7a.best5)

## Seventh best model
m7a.best6 <- lmer(LRR_gammaSn ~ T_warm_rescale * delta.T1_rescale + (1|study),
                  REML = F, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m7a.best6)
r.squaredGLMM(m7a.best6) 
#results_diagnostics(m7a.best6)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                    ENS PIE 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d6.beta <- lmer(LRR_gammaSPie ~ T_warm_rescale * delta.T1_rescale + duration_rescale + (1|study),
                REML = F, data = data_list$Div.data_gamma, na.action = "na.fail") %>% 
  dredge()

model.avg(d6.beta, subset = delta < 2) %>% summary()

## Best model
m6.best <- lmer(LRR_gammaSPie ~ T_warm_rescale + (1|study),
                REML = T, data = data_list$Div.data_gamma, na.action = "na.fail")

Anova(m6.best)
r.squaredGLMM(m6.best)
#results_diagnostics(m6.best)

## Second best model
m6.best1 <- lmer(LRR_gammaSPie ~ delta.T1_rescale + T_warm_rescale + (1|study),
                 REML = T, data = Div.data_gamma, na.action = "na.fail")

Anova(m6.best1)
r.squaredGLMM(m6.best1)
#results_diagnostics(m6.best1)

## Third best model
m6.best2 <- lmer(LRR_gammaSPie ~ duration_rescale + T_warm_rescale + (1|study),
                 REML = T, data = Div.data_gamma, na.action = "na.fail")

Anova(m6.best2)
r.squaredGLMM(m6.best2)
#results_diagnostics(m6.best2)