# r script for analysis of within seasonal guano stain variation
# creator: Alexandra Strang
# created: 2025

# Load packages
library(ggplot2)
library(ggpubr)
library(dplyr)
library(corrplot)
library(car) # for VIF values
library(nlme) # Use REML for low sample sizes
library(performance) # can't use vif for lmm

# Read in data
Dataset.1.0 <- read.csv("Merged_masterdata.csv")

# Colours for each colony
colours <- c("darkblue","royalblue","skyblue")

##########################################################################
# Create dataset with only within season data and covariates
##########################################################################

# Extract only within season data and needed variables 
Dataset.1.1 <- Dataset.1.0[,c("Colony_name","Colony_code","Analysis","Season","Date",
"GA","MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE")]

# Keep only intra-seaonal data
Dataset.1.1$Analysis[Dataset.1.1$Analysis == ""] <- NA
Dataset.1.1$Analysis[Dataset.1.1$Analysis == "NA"] <- NA
sum(is.na(Dataset.1.1$Analysis))

# Remove NAs
Dataset.1.2 <- na.omit(Dataset.1.1)

# Condense/ remove analysis column
Dataset.1.3 <- Dataset.1.2[,c("Colony_name","Colony_code","Season","Date","GA",
"MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE")]

# View
View(Dataset.1.3)

# Turn image date into days since December 1st (Day_D1)

# convert image dates to date variable in r
Dataset.1.3$r_date <- as.character(Dataset.1.3$Date)
Dataset.1.3$r_date <- as.Date(Dataset.1.3$r_date, format="%Y%m%d")

# Calculate day deviation since December 1st
Dataset.1.3$Day_D1 <- as.numeric(difftime(Dataset.1.3$r_date,
                               as.Date(paste0(
                                 ifelse
                                 (format(Dataset.1.3$r_date, "%m") == "12",
                                   format(Dataset.1.3$r_date, "%Y"),
                                   as.numeric(format(Dataset.1.3$r_date, "%Y")) 
                                   -1), "-12-01")), units = "days"))
# Days since December 1st: numerical but discrete without decimals 

# Condense/ remove r date column
Dataset.1.4 <- Dataset.1.3[,c("Colony_name","Colony_code","Season","Date","GA", "Day_D1",
"MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE")]

# View
View(Dataset.1.4)

# Change colony_name to colony for figures
Dataset.1.4$Colony <- Dataset.1.4$Colony_name

##########################################################################
# Visualise within season data by colony across time
##########################################################################

# Subset the data by the three colonies
ADARdf <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="ADAR")
CROZdf <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="CROZ")
HALLdf <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="HALL")

# Plot guano areas over time from days since December 1st for each colony

# Adare
ADAR_plot <- ggplot(ADARdf, aes(x = Day_D1, y = GA)) +
  geom_point(colour = "royalblue") +
  geom_line(colour = "royalblue") +
  xlab("Days since December 1st") +
  ylab("Guano area (m2)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

ADAR_plot

# Crozier
CROZ_plot <- ggplot(CROZdf, aes(x = Day_D1, y = GA)) +
  geom_point(colour = "darkblue") +
  geom_line(colour = "darkblue") +
  xlab("Days since December 1st") +
  ylab("Guano area (m2)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

CROZ_plot

# Hallett
HALL_plot <- ggplot(HALLdf, aes(x = Day_D1, y = GA)) +
  geom_point(colour = "skyblue") +
  geom_line(colour = "skyblue") +
  xlab("Days since December 1st") +
  ylab("Guano area (m2)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

HALL_plot

# Plot together
All_3_plot <- plot(ggarrange(ADAR_plot,
                              CROZ_plot,
                              HALL_plot,
                              ncol = 1,
                              nrow = 3,
                              labels = c("Cape Adare", "Cape Crozier", "Cape Hallett"),
                              label.x = 0.1))


# Try guano area on log scale
Dataset.1.4$Log_GA <- log(Dataset.1.4$GA)

# View
View(Dataset.1.4)

# Subset the data by the three colonies
ADARdf2 <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="ADAR")
CROZdf2 <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="CROZ")
HALLdf2 <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="HALL")

# Adare
ADAR_plot2 <- ggplot(ADARdf2, aes(x = Day_D1, y = Log_GA)) +
  geom_point(colour = "royalblue") +
  geom_line(colour = "royalblue") +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

ADAR_plot2

# Crozier
CROZ_plot2 <- ggplot(CROZdf2, aes(x = Day_D1, y = Log_GA)) +
  geom_point(colour = "darkblue") +
  geom_line(colour = "darkblue") +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

CROZ_plot2

# Hallett
HALL_plot2 <- ggplot(HALLdf2, aes(x = Day_D1, y = Log_GA)) +
  geom_point(colour = "skyblue") +
  geom_line(colour = "skyblue") +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

HALL_plot2

# Plot together
All_3_plot <- plot(ggarrange(ADAR_plot2,
                             CROZ_plot2,
                             HALL_plot2,
                             ncol = 1,
                             nrow = 3,
                             labels = c("Cape Adare", "Cape Crozier", "Cape Hallett"),
                             label.x = 0.1))
All_3_plot <- annotate_figure(All_3_plot, left = "Log guano area (m²)", bottom = "Days since December 1st")

# Figure 3 of Chap 1
ggsave("Chap1_outputs/Seasonal_GA_by_colony.png", All_3_plot,
       width = 8, height = 5, units = "in",
       dpi = 600)

# View on the same plot
All_together <- ggplot(Dataset.1.4, aes(x = Day_D1, y = Log_GA, group = Colony_code)) +
  geom_point(aes(colour = Colony_code)) +
  geom_line(aes(colour = Colony_code)) +
  xlab("Days since December 1st") +
  ylab("Log guano area (m²)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_color_manual(values = colours) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

All_together

# Differences between and within colonies (means and standard deviations)

# Plot as box plot
Box_plot <- ggplot(Dataset.1.4, aes(x = Colony, y = Log_GA, fill = Colony)) +
  geom_boxplot() +
  xlab("Colony") +
  ylab("Log guano area (m²)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "none"
        ) +
  scale_fill_manual(values = colours)

Box_plot

# Figure 4 of Chap 1
ggsave("Chap1_outputs/Seasonal_GA_boxplot.png", Box_plot,
       width = 8, height = 5, units = "in",
       dpi = 600)

# Calculate means and standard deviations for each colony

# Create table
summary_table <- Dataset.1.4 %>%
  group_by(Colony, Season) %>%
  summarise(n = n(),
            mean_GA = mean(GA),
            median_GA = median(GA),
            sd = sd(GA))

# Create a table plot
table_plot <- ggtexttable(summary_table, rows = NULL)

# Show table in Plots window
plot(table_plot)

# The larger the colony, the larger the variation

##########################################################################
# Test for multi-colinearity between covariates
##########################################################################

# Investigate correlation between covariates using correlogram
covariates <- data.frame(Dataset.1.4$Day_D1, Dataset.1.4$MEANCOLLECTEDGSD,Dataset.1.4$MEANSUNEL, 
                         Dataset.1.4$MEANSUNAZ, Dataset.1.4$MEANOFFNADIRVIEWANGLE)

cor.matrix <- cor(covariates) # default method is pearsons

corrplot(cor.matrix, method = "number", type = "lower", tl.cex = 1)
# plot

# correlated variables of greater than +/- 0.50
# sun elevation and days since December first correlated - expected
# check correlations with vif analysis too

##########################################################################
# Build base model with satellite-related factors
##########################################################################

# non-log base model without random effect of colony
base_model <- lm(Dataset.1.4$GA ~ Dataset.1.4$Day_D1 + Dataset.1.4$MEANCOLLECTEDGSD + 
                   Dataset.1.4$MEANSUNAZ + Dataset.1.4$MEANSUNEL +
                   Dataset.1.4$MEANOFFNADIRVIEWANGLE)

summary(base_model) # too many parameters for the no. of observations

# Day since Dec 1 is an integer = not continuous
vif(base_model) # based on vif remove SUNEL or Day_D1

# check residuals
par(mfrow = c(2,2))
plot(base_model)

Dataset.1.4$fittedbase <- fitted(base_model)
Dataset.1.4$residbase <- resid(base_model)

Resids_base <- ggplot(Dataset.1.4, aes(x=fittedbase, y=residbase, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-400000,400000), breaks = seq(-400000,400000, by=50000))

Resids_base
# resid pattern as it doesn't account for colony effect

# try log guano area due to colonies of different sizes and pattern in resids
Dataset.1.4$Log_GA <- log(Dataset.1.4$GA)

log_base_model <- lm(Dataset.1.4$Log_GA ~ Dataset.1.4$Day_D1 + Dataset.1.4$MEANCOLLECTEDGSD + 
                   Dataset.1.4$MEANSUNAZ + Dataset.1.4$MEANSUNEL +
                   Dataset.1.4$MEANOFFNADIRVIEWANGLE)

summary(log_base_model)
vif(log_base_model) # based on vif remove SUNEL or Day_D1

# check residuals
par(mfrow = c(2,2))
plot(log_base_model)

Dataset.1.4$fittedlog <- fitted(log_base_model)
Dataset.1.4$residlog <- resid(log_base_model)

Resids_log <- ggplot(Dataset.1.4, aes(x=fittedlog, y=residlog, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-2.0,2.0), breaks = seq(-2.0, 2.0, by=0.2))

Resids_log
# also doesn't account for colony effect
# spread of residuals is huge (1.0 to -1.2)

# Include random effect of colony

# base mixed-effects model
log_base_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANCOLLECTEDGSD + MEANSUNAZ + MEANSUNEL + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_base_lmm) #check other correlations
anova(log_base_lmm)

check_collinearity(log_base_lmm)
# remove Days since December 1st (correlated with sun el)
# days since December first may be significant
# colonies of different sizes are significantly different

# check residuals
par(mfrow = c(1,1))
plot(log_base_lmm)

Dataset.1.4$fittedlmm <- fitted(log_base_lmm)
Dataset.1.4$residlmm <- resid(log_base_lmm)

Resids_lmm <- ggplot(Dataset.1.4, aes(x=fittedlmm, y=residlmm, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-1.0,1.0), breaks = seq(-1.0, 1.0, by=0.2))

Resids_lmm
# resids look okay?

# Try AR(1) - auto regressive error structure (points will be correlated to the previous points)
# with random effect of colony
AR1_Model <- lme(
  fixed = Log_GA ~ Day_D1 + MEANCOLLECTEDGSD + MEANSUNAZ + MEANSUNEL + MEANOFFNADIRVIEWANGLE, # or try null 
  random = ~ 1 | Colony_code, # random effect
  data = Dataset.1.4, 
  correlation = corAR1(form = ~ Day_D1 | Colony_code) # AR(1) correlation within colonies over time
) 

summary(AR1_Model) 
# Phi1 = 0 (means there is no autocorrelation)
# phi1 is zero indicating that auto regressive error structure isn't needed
# no auto regressive correlation is detected 

anova(log_base_lmm, AR1_Model)
# model without autoregressive structure is better 
# no evidence that adding the AR(1) structure improves model fit

# check residuals
qqnorm(resid(AR1_Model))  # Q-Q plot for residuals
qqline(resid(AR1_Model))  # reference line

Dataset.1.4$fittedAR <- fitted(AR1_Model)
Dataset.1.4$residAR <- resid(AR1_Model)

ResidsAR <- ggplot(Dataset.1.4, aes(x=fittedAR, y=residAR, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-1.0,1.0), breaks = seq(-1.0, 1.0, by=0.2))

ResidsAR
# look okay

acf(resid(AR1_Model))
# no autocorrelation
# check again with top model

##########################################################################
# Build candidate models with satellite-related factors
##########################################################################

# compare to log_base_lmm
# include only max of two covariates
# need to use subsets of base model excluding collinear variables

# removed days since December 1st first?

# single fixed-effects

# off-nadir model
log_reduced1_lmm <- lme(
  fixed = Log_GA ~ MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced1_lmm)
anova(log_reduced1_lmm)
M1_AIC <- AIC(log_reduced1_lmm)
M1_AICc <- AICc(log_reduced1_lmm)

# gsd model
log_reduced2_lmm <- lme(
  fixed = Log_GA ~ MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced2_lmm)
anova(log_reduced2_lmm)
M2_AIC <- AIC(log_reduced2_lmm)
M2_AICc <- AICc(log_reduced2_lmm)

# sun az model
log_reduced3_lmm <- lme(
  fixed = Log_GA ~ MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced3_lmm)
anova(log_reduced3_lmm)
M3_AIC <- AIC(log_reduced3_lmm)
M3_AICc <- AICc(log_reduced3_lmm)

# sun el model
log_reduced4_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced4_lmm)
anova(log_reduced4_lmm)
M4_AIC <- AIC(log_reduced4_lmm)
M4_AICc <- AICc(log_reduced4_lmm)
# best model so far

# Day_d1 model ?????
log_reduced5_lmm <- lme(
  fixed = Log_GA ~ Day_D1,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced5_lmm)
anova(log_reduced5_lmm)
M5_AIC <- AIC(log_reduced5_lmm)
M5_AICc <- AICc(log_reduced5_lmm)
# sun el model better

# 9 candidate models with two fixed effects

# gsd and off-nadir model
log_reduced6_lmm <- lme(
  fixed = Log_GA ~ MEANCOLLECTEDGSD + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced6_lmm)
anova(log_reduced6_lmm)
M6_AIC <- AIC(log_reduced6_lmm)
M6_AICc <- AICc(log_reduced6_lmm)

# sun el and off-nadir model
log_reduced7_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced7_lmm)
anova(log_reduced7_lmm)
M7_AIC <- AIC(log_reduced7_lmm)
M7_AICc <- AICc(log_reduced7_lmm)

# sun az and off-nadir model
log_reduced8_lmm <- lme(
  fixed = Log_GA ~ MEANSUNAZ + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced8_lmm)
anova(log_reduced8_lmm)
M8_AIC <- AIC(log_reduced8_lmm)
M8_AICc <- AICc(log_reduced8_lmm)

# sun el and gsd model
log_reduced9_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced9_lmm)
anova(log_reduced9_lmm)
M9_AIC <- AIC(log_reduced9_lmm)
M9_AICc <- AICc(log_reduced9_lmm)

# sun el and sun az model
log_reduced10_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced10_lmm)
anova(log_reduced10_lmm)
M10_AIC <- AIC(log_reduced10_lmm)
M10_AICc <- AICc(log_reduced10_lmm)

# gsd and sun az model
log_reduced11_lmm <- lme(
  fixed = Log_GA ~ MEANCOLLECTEDGSD + MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced11_lmm)
anova(log_reduced11_lmm)
M11_AIC <- AIC(log_reduced11_lmm)
M11_AICc <- AICc(log_reduced11_lmm)

# Day_D1 and gsd model
log_reduced12_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced12_lmm)
anova(log_reduced12_lmm)

# Day_D1 and sun az model
log_reduced13_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced13_lmm)
anova(log_reduced13_lmm)

# Day_D1 and off-nadir model
log_reduced14_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced14_lmm)
anova(log_reduced14_lmm)

# null model
null_model <- lme(
  fixed = Log_GA ~ 1,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(null_model)
anova(null_model)
null_AIC <- AIC(null_model)
null_AICc <- AICc(null_model)

anova(null_model, log_reduced4_lmm)
# null model better?

# model selection using mumin (without null)
model_list <- list(log_reduced1_lmm, log_reduced2_lmm, log_reduced3_lmm, log_reduced4_lmm, 
                   log_reduced5_lmm, log_reduced6_lmm, log_reduced7_lmm, log_reduced8_lmm,
                   log_reduced9_lmm, log_reduced10_lmm, log_reduced11_lmm, log_reduced12_lmm,
                   log_reduced13_lmm, log_reduced14_lmm)
model_selection <- model.sel(model_list) # default rank AICc
model_selection
# single fixed effects model rank highest
# null would be ranked higher with highest weight

# model averaging
avg_model <- model.avg(model_selection, subset = delta < 2) # models with 2 delta AICc
summary(avg_model)
# full average for conservative estimates - should variable be in model at all
# conditional average for effect size - if it is present in the model

# do i need to exlcude day d1?

# calculate delta AICc scores and weights

# for candidate models within delta 2 AIC
x <- c(null_AIC, M4_AIC, M1_AIC, M2_AIC)

# AICc
x <- c(null_AICc, M4_AICc, M1_AICc, M2_AICc)

# compute delta AICc
delta <- x - min(x)

# compute relative likelihoods
rel.lik <- exp(-0.5 * delta)
rel.lik

model_names <- c("Null (random intercpets) model", "Sun elevation angle", "Off-nadir angle", "Ground resolution")

AICc_table <- data.frame(
  Model = model_names,
  AICc = x,
  Delta_AICc = delta,
  Akaike_Weight = Weights(x),
  Relative_Likelihood = rel.lik
)

AICc_table[, 2:3] <- round(AICc_table[, 2:3], 3)
print(AICc_table)

# variation in guano area is more influential then a change in size

# check residuals of top model containing a fixed effect
qqnorm(resid(log_reduced4_lmm))  # Q-Q plot for residuals
qqline(resid(log_reduced4_lmm))  # reference line

Dataset.1.4$fittedbest <- fitted(log_reduced4_lmm)
Dataset.1.4$residbest <- resid(log_reduced4_lmm)

Best_resids <- ggplot(Dataset.1.4, aes(x=fittedbest, y=residbest, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4, 0.4, by=0.2))

Best_resids
# still spreading

# null resids
qqnorm(resid(null_model))  # Q-Q plot for residuals
qqline(resid(null_model))  # reference line

Dataset.1.4$fittednull <- fitted(null_model)
Dataset.1.4$residnull <- resid(null_model)

null_resids <- ggplot(Dataset.1.4, aes(x=fittednull, y=residnull, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4, 0.4, by=0.2))

null_resids
# resids of the best model look better than the null

# acf on null model
AR1_null <- lme(
  fixed = Log_GA ~ 1,
  random = ~ 1 | Colony_code,
  correlation = corAR1(form = ~ Day_D1 | Colony_code),
  data = Dataset.1.4
)
summary(AR1_null) # phi 0
anova(AR1_null)

acf(resid(AR1_null))
# no autocorrelation

anova(null_model, AR1_null)
# model without autoregressive structure is better 
# no evidence that adding the AR(1) structure improves model fit

##########################################################################
# ANOVA test for February effect
##########################################################################

# subset within season data into estimates before and after 59 days since December 1st
# 59 days since December 1st represents late January
# the last two estimates for each colony are "Feb" and rest are "Pre_Feb"

# Subset the data into the first three points and the last two (Pre_Feb and Feb)
Dataset.1.4$Feb_effect <- ifelse(Dataset.1.4$Day_D1 < 59, "Pre_Feb", "Feb")
Dataset.1.4$Feb_effect <- as.factor(Dataset.1.4$Feb_effect)

View(Dataset.1.4)

# ANOVA that tests whether guano area estimates differ between Pre_Feb and Feb
# also accounts for differences in colony size and the feb effect dependent on colony
anova_model <- aov(Log_GA ~ Feb_effect * Colony_code, data = Dataset.1.4)
summary(anova_model)
# used log guano
# significantly different between Pre_Feb and Feb, where the effect depends on colony

# test February effect accounting for colony differences with random effect

# Feb effect model
Feb_model <- lme(
  fixed = Log_GA ~ Feb_effect,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(Feb_model)
anova(Feb_model)
Feb_AIC <- AIC(Feb_model)
Feb_AICc <- AICc(Feb_model)

anova(null_model, Feb_model)
# Feb effect model is better
# Suggests that the last two data points are higher than the rest

# AICc
x <- c(Feb_AICc, null_AICc, M4_AICc, M1_AICc, M2_AICc)
Weights(x)

# compute delta AICc
delta <- x - min(x)

# compute relative likelihoods
rel.lik <- exp(-0.5 * delta)
rel.lik

model_names <- c("February effect", "Null (random intercpets) model", "Sun elevation angle", "Off-nadir angle", "Ground resolution")

AICc_table <- data.frame(
  Model = model_names,
  AICc = x,
  Delta_AICc = delta,
  Akaike_Weight = Weights(x),
  Relative_Likelihood = rel.lik
)

AICc_table[, 2:5] <- round(AICc_table[, 2:5], 2)
print(AICc_table)

# Feb model resids
qqnorm(resid(Feb_model))  # Q-Q plot for residuals
qqline(resid(Feb_model))  # reference line
# pattern?

Dataset.1.4$fittedfeb <- fitted(Feb_model)
Dataset.1.4$residfeb <- resid(Feb_model)

Feb_resids <- ggplot(Dataset.1.4, aes(x=fittedfeb, y=residfeb, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
  geom_point(size=3) + 
  geom_hline(yintercept = 0) + 
  xlab("Observed") + 
  ylab("Residuals") + 
  scale_shape_manual(values = c(21,21,21)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size=12),
        axis.text.y = element_text(colour = 'black', size=12),) + 
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4, 0.4, by=0.2))

Feb_resids

# acf on Feb model
AR1_Feb <- lme(
  fixed = Log_GA ~ Feb_effect,
  random = ~ 1 | Colony_code,
  correlation = corAR1(form = ~ Day_D1 | Colony_code),
  data = Dataset.1.4
)
summary(AR1_Feb) # phi 0
anova(AR1_Feb)

acf(resid(AR1_Feb))
# no autocorrelation

anova(Feb_model, AR1_Feb)
# model without autoregressive structure is better 
# no evidence that adding the AR(1) structure improves model fit

##########################################################################
# GA and BP relationship (model created in Strang MSc thesis)
##########################################################################

# update to include all 15 colonies

# Extract only interseasonal data and needed variables 
Dataset.2.0 <- Dataset.1.0[,c("Colony_name","GA","BP","Analysis2","Date")]
View(Dataset.2.0)

# Keep only interseaonal data
# interseasonal analysis for within season images is median estimate date
Dataset.2.0$Analysis2[Dataset.2.0$Analysis2 == ""] <- NA
Dataset.2.0$Analysis2[Dataset.2.0$Analysis2 == "NA"] <- NA
sum(is.na(Dataset.2.0$Analysis2))

# Condense/ remove analysis column

# Remove NAs
Dataset.2.1 <- na.omit(Dataset.2.0)
View(Dataset.2.1)

Dataset.2.1$Log_GA <- log(Dataset.2.1$GA)
Dataset.2.1$Log_BP <- log(Dataset.2.1$BP)

GA_BP <- lm(Dataset.2.1$Log_GA ~ Dataset.2.1$Log_BP)
summary(GA_BP)
r2 <- round(summary(GA_BP)$r.squared, 2)

colours_14 <- c("orchid","navy", "indianred", "royalblue", "skyblue", "red", 
             "pink", "orange", "gold", "slateblue", 
             "magenta", "cyan", "grey", "violetred")

All_plot <- ggplot(Dataset.2.1, aes(x = Log_BP, y = Log_GA, colour = Colony_name)) + 
  geom_point(size=3) + 
  geom_smooth(method="lm", col = "black") +
  annotate("text", x = 8, y = 13, 
           label = paste0("R² = ", r2)) +
  xlab("Log BP") +
  ylab("Log Guano area (m2)") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  labs(color = "Colony") +
  scale_color_manual(values = colours_14)

All_plot

##########################################################################
# Test for influence of Feb estimates on BP predictions
##########################################################################

# convert image dates to date variable in r
Dataset.2.1$r_date <- as.character(Dataset.2.1$Date)
Dataset.2.1$r_date <- as.Date(Dataset.2.1$r_date, format="%Y%m%d")

# Calculate day deviation since December 1st
Dataset.2.1$Day_D1 <- as.numeric(difftime(Dataset.2.1$r_date,
                                          as.Date(paste0(
                                            ifelse
                                            (format(Dataset.2.1$r_date, "%m") == "12",
                                              format(Dataset.2.1$r_date, "%Y"),
                                              as.numeric(format(Dataset.2.1$r_date, "%Y")) 
                                              -1), "-12-01")), units = "days"))
# Days since December 1st: numerical but discrete without decimals 

# Condense/ remove r date column
Dataset.2.2 <- Dataset.2.1[,c("Colony_name","Log_GA","Log_BP","Date","GA", "Day_D1")]

# subset data into estimates before and after 59 days since December 1st
# 59 days since December 1st represents late January
# estimates after 59 days are "Feb" and rest are "Pre_Feb"

# Subset the data into Pre_Feb and Feb
Dataset.2.2$Feb_effect <- ifelse(Dataset.2.2$Day_D1 < 59, "Pre_Feb", "Feb")
Dataset.2.2$Feb_effect <- as.factor(Dataset.2.2$Feb_effect)

# View
View(Dataset.2.2)

# test Feb effect in GA ~ BP model
Feb_GA_BP <- lm(Dataset.2.2$Log_GA ~ Dataset.2.2$Log_BP + Dataset.2.2$Feb_effect)
summary(Feb_GA_BP)
# Feb effect not significant

# mixed Feb effect model
Feb_mixed <- lme(
  fixed = Log_GA ~ Log_BP + Feb_effect,
  random = ~ 1 | Colony_name,
  data = Dataset.2.2
)
summary(Feb_mixed)
anova(Feb_mixed)
# Feb effect not significant

# Make BP predictions from GA~BP model (model created in Strang MSc thesis)

# model: Log GA = a + b * Log BP
# a is intercept
# b is slope

# extract model covariates 
a <- coef(GA_BP)[1]  # intercept
b <- coef(GA_BP)[2]  # slope

# invert by:
# Log_BP = (Log_GA - a) / b
# BP = exp((log(GA) - a) / b)

# GA estimates across Dec - Feb for Adare, Crozier, and Hallett

# Extract only needed variables 
Dataset.1.5 <- Dataset.1.4[,c("Colony_code","Season","Date","GA","Day_D1","Log_GA")]
View(Dataset.1.5)

# predict log BP from intra-seasonal GA using updated GA ~ BP model
Dataset.1.5$Log_BP_Pred <- NA  # initialise column
Dataset.1.5$Log_BP_Pred <- (Dataset.1.5$Log_GA - a) / b

# back-transform to get estimated BP
Dataset.1.5$BP_Pred <- exp(Dataset.1.5$Log_BP_Pred)

View(Dataset.1.5)

# scale to view GA and BP on same graph
# find scale range of Log_GA
GA_range <- range(Dataset.1.5$Log_GA, na.rm = TRUE)
BP_range <- range(Dataset.1.5$BP_Pred, na.rm = TRUE)

# rescale BP_Pred to match the Log_GA range
Dataset.1.5$BP_Pred_scaled <- (Dataset.1.5$BP_Pred - BP_range[1]) / diff(BP_range) * diff(GA_range) + GA_range[1]

# View on the same plot
GA_BP_preds <- ggplot(Dataset.1.5, aes(x = Day_D1)) +
  # Guano Area
  geom_point(aes(y = Log_GA, colour = Colony_code), size = 2) +
  geom_line(aes(y = Log_GA, colour = Colony_code), linewidth = 1) +
  
  # predicted BP
  geom_point(aes(y = BP_Pred_scaled, shape = Colony_code), size = 2, alpha = 0.5) +
  geom_line(aes(y = BP_Pred_scaled, linetype = Colony_code), alpha = 0.5, linewidth = 1) +
  
  # axis labels
  xlab("Days since December 1st") +
  scale_y_continuous(
    name = "Log guano area (m²)",
    sec.axis = sec_axis(
      transform = ~ (.-GA_range[1]) * diff(BP_range)/diff(GA_range) + BP_range[1],
      name = "Predicted breeding pairs"
    )
  ) +
  scale_x_continuous(limits = c(1, 90), breaks = seq(1, 90, by = 10)) +
  scale_color_manual(values = colours) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

GA_BP_preds

# just BP preds
BP_preds <- ggplot(Dataset.1.5, aes(x = Day_D1, y = BP_Pred, group = Colony_code)) +
  geom_point(aes(colour = Colony_code)) +
  geom_line(aes(colour = Colony_code)) +
  xlab("Days since December 1st") +
  ylab("Predicted breeding pairs") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_color_manual(values = colours) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

BP_preds
# need to determine if change in predicted BP is large enough to care?
