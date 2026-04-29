rm(list = ls())
gc()

## Load packages
packages <- c("ggplot2", "ggpubr", "lme4", "gridExtra", "effects", "car", "devtools", 
              "plotbiomes", "dplyr", "ggdist", "gghalves", "here")

lapply(packages, function(pkg) 
{
  if (!requireNamespace(pkg, quietly = T)) 
    install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})

## Load necessary functions from Github (Valentin, https://github.com/valentinitnelav)
devtools::install_github("valentinitnelav/plotbiomes")
devtools::install_github('erocoar/gghalves')

## Load data, it will load the function from a R project, so you either need to create on
## or you replace here(data) by a path leading to where you stored your data
csv_read <- c(Div.data_local = "Temp_data_final.csv", Div.data_gamma = "Temp_data_gamma_final.csv", 
              Comp_data_final = "Comp_data_final.csv", biome = "Biome.csv", 
              grand_mean = "grand_mean.csv")

data_list <- Map(function(data)
{read.csv(here(data), dec = ".", sep = ";", h = T)}, csv_read)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Figure 1 biome plot    -----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
means_stud1 <- data_list$biome %>% dplyr::select(ID, MAT_chelsa, MAP_chelsa) %>% 
  mutate(MAP_chelsa = MAP_chelsa * 0.1)

Fig.1 <- whittaker_base_plot() +
  geom_point(data = means_stud1, aes(x = MAT_chelsa, y = MAP_chelsa), size = 3,
    shape = 16) +
  labs(title = "Figure 1") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0, face = "bold"), 
    axis.text = element_text(size = 14, colour = "black"), 
    axis.title = element_text(size = 14))

Fig.1
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      Figure 2 Species richness -----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fig.2a <- ggplot(data_list$grand_mean, aes(y = ES_S, x = scale)) +
  labs(y = "LRR species richness", x = "scale") +
  geom_errorbar(width = 0, linewidth = 0.3, aes(ymin = lower_CI_S, ymax = upper_CI_S)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 6, color = "mediumpurple") +
  scale_x_discrete(labels = c("alpha", "gamma", "beta"))+
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        aspect.ratio = 1)

###delta T alpha
m1a.delta.T <- lmer(LRR_S ~ delta.T + treat.duration.y + T_warm_chelsa + (1|study:site:SiteBlock),
                  REML = T, data = data_list$Div.data_local, na.action = "na.fail")

newdat.lmer1 <- data_list$Div.data_local %>% 
  dplyr::select(LRR_S, delta.T, treat.duration.y, T_warm_chelsa) %>%
  mutate(treat.duration.y = median(treat.duration.y), 
         T_warm_chelsa = median(T_warm_chelsa))
newdat.lmer1$predlmer <- predict(m1a.delta.T, newdata = newdat.lmer1, re.form = NA)

des <- model.matrix(terms(m1a.delta.T), newdat.lmer1)

predvar <- diag(des %*% vcov(m1a.delta.T) %*% t(des))

## Calculate CI's
newdat.lmer1$lower <- with(newdat.lmer1, predlmer - 1.96 * sqrt(predvar))
newdat.lmer1$upper <- with(newdat.lmer1, predlmer + 1.96 * sqrt(predvar))

Fig.2b <- ggplot(data_list$Div.data_local, aes(x = delta.T, y = LRR_S)) +
  labs(y = "LRR species richness", x = "Magnitude of manipulation (°C)") +
  annotate("text", x = 3, y = 1, label = "alpha", size = 5) +
  scale_x_continuous(limits = c(0, 5), breaks = c(0:5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(size = 3, color = "mediumpurple", stroke = 0.3, alpha = 0.4) + 
  geom_ribbon(data = newdat.lmer1, aes(y = NULL, ymin = lower, ymax = upper, 
                                       color = NULL),alpha = .15) +
  geom_line(data = newdat.lmer1, aes(y = predlmer), size = .75) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        aspect.ratio = 1) 

#delta T gamma
Fig.2c <- ggplot(data_list$Div.data_gamma, aes(x = delta.T, y = LRR_gammaSn)) +
  labs(y = "LRR species richness", x = "Magnitude of manipulation (°C)") +
  annotate("text", x = 3, y = 1, label = "gamma", size = 5) +
  scale_x_continuous(limits = c(0,5), breaks = c(0:5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  geom_point(size = 3, color = "mediumpurple", stroke = 0.3, alpha = 0.4)+ 
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        aspect.ratio = 1) 

Fig.2d <- ggplot(data_list$Div.data_local, aes(x = delta.T, y = LRR_betaSn)) +
  labs(y = "LRR species richness", x = "Magnitude of manipulation (°C)") +
  annotate("text", x = 3, y = 1, label = "beta", size = 5) +
  scale_x_continuous(limits = c(0,5), breaks = c(0:5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  geom_hline(yintercept = 0, linetype =2, size = 0.3) +
  geom_point(size = 3, color = "mediumpurple", stroke = 0.3, alpha = 0.4 ) + 
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        aspect.ratio = 1) 

#T warm alpha
m1a.T_warm <- lmer(LRR_S ~ T_warm_chelsa + treat.duration.y + (1|study:site:SiteBlock), 
                   REML = T, data = data_list$Div.data_local, na.action = "na.fail")

newdat.lmer1 <- data_list$Div.data_local %>% 
  dplyr::select(LRR_S, T_warm_chelsa, treat.duration.y) %>% 
  mutate(treat.duration.y = median(treat.duration.y))

newdat.lmer1$predlmer <- predict(m1a.T_warm, newdata = newdat.lmer1, re.form = NA)

des <- model.matrix(terms(m1a.T_warm), newdat.lmer1)

predvar <- diag(des %*% vcov(m1a.T_warm) %*% t(des))

newdat.lmer1$lower <- with(newdat.lmer1, predlmer - 1.96 * sqrt(predvar))
newdat.lmer1$upper <- with(newdat.lmer1, predlmer + 1.96 * sqrt(predvar))

Fig.2e <- ggplot(data_list$Div.data_local, aes(x = T_warm_chelsa, y = LRR_S)) +
  labs(y = "LRR species richness", x = "Warmest quarter (°C)") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3)+
  geom_point(size = 3, color =  "mediumpurple", stroke = 0.3, alpha = 0.4)+ 
  geom_ribbon(data = newdat.lmer1, aes(y = NULL, ymin = lower, ymax = upper, 
                                       color = NULL),alpha = .15) +
  geom_line(data = newdat.lmer1, aes(y = predlmer), size = .75)+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        aspect.ratio = 1) 

###T warm gamma
Fig.2f <- ggplot(data_list$Div.data_gamma, aes(x = T_warm_chelsa, y = LRR_gammaSn)) +
  labs(y = "LRR species richness", x = "Warmest quarter (°C)") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  geom_point(size = 3, color = "mediumpurple", stroke = 0.3, alpha = 0.4) + 
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        aspect.ratio = 1) 

Fig.2g <- ggplot(data_list$Div.data_local, aes(x = T_warm_chelsa, y = LRR_betaSn)) +
  labs(y = "LRR species richness", x = "Warmest quarter (°C)") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  geom_point(size = 3, color =  "mediumpurple", stroke = 0.3, alpha = 0.4)+ 
  theme_classic() +
  theme(axis.text = element_text(size=12, colour = "black"), 
        axis.title = element_text(size=12),
        aspect.ratio = 1) 

Fig.2 <- ggarrange(Fig.2a, NULL, NULL,
                   Fig.2b, Fig.2c, Fig.2d,
                   Fig.2e, Fig.2f, Fig.2g,
                   labels = c("a","", "", "b", "c", "d", "e", "f", "g"),
                   ncol = 3, nrow = 3)

Fig.2 <- annotate_figure(Fig.2, top = text_grob("Figure 2", face = "bold",
                                                size = 14, x = 0, hjust = 0))

Fig.2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Figure 3 ENS Pie -----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fig.3a <- ggplot(data_list$grand_mean, aes(y = ES_Pie, x = scale)) +
  labs(y = expression(LRR ~ evenness ~ (ENS[Pie])), x = "scale")+
  geom_errorbar(width = 0,linewidth = 0.3, aes(ymin = lower_CI_Pie, ymax = upper_CI_Pie))+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 6, color = "violet") +
  scale_x_discrete(labels=c("alpha", "gamma", "beta"))+
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        aspect.ratio = 1)

## Duration Pie alpha scale 
m3a.duration <- lmer(LRR_SPie ~ T_warm_chelsa + treat.duration.y + (1|study:site:SiteBlock), 
                   REML = T,data = data_list$Div.data_local, na.action="na.fail")

newdat.lmer1 <- data_list$Div.data_local %>% dplyr::select(LRR_SPie, treat.duration.y, T_warm_chelsa) %>% 
  mutate(T_warm_chelsa = median(T_warm_chelsa))

newdat.lmer1$predlmer <- predict(m3a.duration, newdata = newdat.lmer1, re.form = NA)

des <- model.matrix(terms(m3a.duration), newdat.lmer1)

predvar <- diag( des %*% vcov(m3a.duration) %*% t(des))

newdat.lmer1$lower <- with(newdat.lmer1, predlmer - 1.96*sqrt(predvar))
newdat.lmer1$upper <- with(newdat.lmer1, predlmer + 1.96*sqrt(predvar))

Fig.3b <- ggplot(data_list$Div.data_local, aes(x = treat.duration.y, y = LRR_SPie)) +
  labs(y = expression(LRR ~ evenness ~ (ENS[Pie])), x = "Duration (years)") +
  scale_x_continuous(limits = c(0,25), breaks = seq(0,25, 5)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
  annotate("text", x = 12, y = 2, label = "alpha", size = 5) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  geom_point(size = 3, color = "violet", stroke = 0.3, alpha = 0.4) + 
  geom_ribbon(data = newdat.lmer1, aes(y = NULL, ymin = lower, ymax = upper, 
                                      color = NULL),alpha = .15) +
  geom_line(data = newdat.lmer1, aes(y = predlmer), size = .75) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        aspect.ratio = 1) 

## Duration Pie gamma scale 
Fig.3c <- ggplot(data_list$Div.data_gamma, aes(x = treat.duration.y, y = LRR_gammaSPie)) +
  labs(y = expression(LRR ~ evenness ~ (ENS[Pie])), x = "Duration (years)") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
  annotate("text", x = 12, y = 2, label = "gamma", size = 5) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  geom_point(size = 3, color = "violet", stroke = 0.3, alpha = 0.4) + 
   theme_classic()+
  theme(axis.text= element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 12),
        aspect.ratio = 1) 

## Duration Pie beta scale 
m4a.duration <- lmer(LRR_betaSPie ~ treat.duration.y + (1|study:site:SiteBlock), 
                     REML = T, data = data_list$Div.data_local, na.action = "na.fail")

newdat.lmer4 <- data_list$Div.data_local %>% dplyr::select(LRR_betaSPie, treat.duration.y)

newdat.lmer4$predlmer <- predict(m4a.duration, newdata = newdat.lmer4, re.form = NA)

des <- model.matrix(terms(m4a.duration), newdat.lmer4)
predvar <- diag(des %*% vcov(m4a.duration) %*% t(des))

newdat.lmer4$lower = with(newdat.lmer4, predlmer - 1.96*sqrt(predvar) )
newdat.lmer4$upper = with(newdat.lmer4, predlmer + 1.96*sqrt(predvar) )

Fig.3d <- ggplot(data_list$Div.data_local, aes(x = treat.duration.y, y = LRR_betaSPie)) +
  labs(y =  expression(LRR ~ evenness ~ (ENS[Pie])), x = "Duration (years)") +
  annotate("text", x = 12.5, y = 2, label = "beta", size = 5) +
  scale_x_continuous(limits = c(0,25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  geom_point(size = 3, color = "violet", stroke = 0.3, alpha = 0.4) + 
  geom_ribbon(data = newdat.lmer4, aes(y = NULL, ymin = lower, ymax = upper, 
                                       color = NULL), alpha = .15) +
  geom_line(data = newdat.lmer4, aes(y = predlmer), size = .75) +
  theme_classic() +
  theme(axis.text = element_text(size=12, colour = "black"), 
        axis.title = element_text(size=12),
        aspect.ratio = 1) 

Fig.3 <- ggarrange(Fig.3a, NULL, NULL, Fig.3b, Fig.3c, Fig.3d,
         labels = c("a","", "", "b", "c", "d"), ncol = 3, nrow = 2)

Fig.3 <- annotate_figure(Fig.3, top = text_grob("Figure 2", face = "bold",
                                                size = 14, x = 0, hjust = 0))

Fig.3
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##            Figure 4 composition -----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rain <- data_list$Comp_data_final %>%
  select(measure, value) %>%
  filter(measure %in% c("beta.bray.bal", "beta.bray.gra", 
                        "beta.sim", "beta.sne")) %>%
  mutate(measure = factor(measure, levels = c("beta.bray.bal", "beta.sim",
                                              "beta.bray.gra", "beta.sne"),
                          labels = c("Bray–Curtis turnover",
                                     "Sørensen turnover",
                                     "Bray–Curtis nestedness",
                                     "Sørensen nestedness")))

Fig.4 <- ggplot(rain, aes(x = measure, y = value, fill = measure)) +
  geom_half_violin(side = "r", alpha = 0.6, width = 0.8,trim = FALSE,
    scale = "width", colour = NA) +
  geom_point(position = position_jitter(width = 0.06), size = 1, alpha = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("#1B9E77", "#2CA25F", "white", "#1D91C0")) +
  labs(title = "Figure 4", x = NULL, y = "Compositional dissimilarity") +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold", hjust = 0, 
                                  margin = margin(l = 5)),   # moves title further left
    plot.title.position = "plot",
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 14),
    aspect.ratio = 1.1)

Fig.4