# r script for exploratory analysis of within seasonal guano stain variation
# creator: Alexandra Strang
# created: 2025

# Read in data
#Dataset.1.0 <- read.csv("Strang_PhD_masterdata.csv")

Dataset.1.0 <- read.csv("Merged_masterdata.csv")
colours <- c("darkblue","royalblue","skyblue")

# Extract only within season data and needed variables 
Dataset.1.1 <- Dataset.1.0[,c("Colony_code","Analysis","Season","Date",
"GA","MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE")]

# Keep only intra-seaonal data
Dataset.1.1$Analysis[Dataset.1.1$Analysis == ""] <- NA
Dataset.1.1$Analysis[Dataset.1.1$Analysis == "NA"] <- NA
sum(is.na(Dataset.1.1$Analysis))

# Remove NAs
Dataset.1.2 <- na.omit(Dataset.1.1)

# Condense/ remove analysis column
Dataset.1.3 <- Dataset.1.2[,c("Colony_code","Season","Date","GA",
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
Dataset.1.4 <- Dataset.1.3[,c("Colony_code","Season","Date","GA", "Day_D1",
"MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE")]

# View
View(Dataset.1.4)

# Subset the data by the three colonies
ADARdf <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="ADAR")
CROZdf <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="CROZ")
HALLdf <- subset(Dataset.1.4, Dataset.1.4$Colony_code=="HALL")

# Plot guano areas over time from days since December 1st for each colony

library(ggplot2)

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

library(ggpubr)

# Combine
All_3_plot <- plot(ggarrange(ADAR_plot,
                              CROZ_plot,
                              HALL_plot,
                              ncol = 1,
                              nrow = 3,
                              labels = c("Cape Adare", "Cape Crozier", "Cape Hallett"),
                              label.x = 0.1))

# Linear model of guano areas and days since December 1st
LM_1 <- lm(Dataset.1.4$GA ~ Dataset.1.4$Day_D1)
summary(LM_1)

# Plot residuals
qqnorm(resid(LM_1))  # Q-Q plot for residuals
qqline(resid(LM_1))  # reference line

Dataset.1.4$fitted1 <- LM_1$fitted.values
Dataset.1.4$resid1 <- LM_1$residuals

Dataset.1.4$colony <- as.factor(Dataset.1.4$Colony_code)

Resids_1 <- ggplot(Dataset.1.4, aes(x=fitted1, y=resid1, colour = colony, fill = colony, shape = colony)) + 
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
  scale_y_continuous(limits = c(-500000,500000), breaks = seq(-500000,500000, by=100000))

Resids_1
# best residuals but not accounting for colony
par(mfrow = c(2,2))
plot(LM_1)

# Pools data across colonies but colonies vary in size
# Need to include the influence of colony

# Linear model of guano areas with days since December 1st and colony code
LM_6 <- lm(Dataset.1.4$GA ~ Dataset.1.4$Day_D1 + Dataset.1.4$Colony_code)
summary(LM_6)

par(mfrow = c(2,2))
plot(LM_6)


# Account for colony as a random effect

library(lme4)

# Linear mixed effect model of guano areas and days since December 1st
# Including a random effect of colony
LMM_1 <- lme4::lmer(GA ~ Day_D1 + (1|Colony_code), data = Dataset.1.4)
summary(LMM_1)

# Plot residuals
Dataset.1.4$fitted1.5 <- fitted(LMM_1)
Dataset.1.4$resid1.5 <- resid(LMM_1)

qqnorm(Dataset.1.4$resid1.5)  # Q-Q plot for residuals
qqline(Dataset.1.4$resid1.5)  # reference line

Resids_1.5 <- ggplot(Dataset.1.4, aes(x=fitted1.5, y=resid1.5, colour = colony, fill = colony, shape = colony)) + 
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
  scale_y_continuous(limits = c(-200000,200000), breaks = seq(-200000,200000, by=50000))

Resids_1.5

# Heteroskedasticity in the residuals

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

# Combine
All_3_plot <- plot(ggarrange(ADAR_plot2,
                             CROZ_plot2,
                             HALL_plot2,
                             ncol = 1,
                             nrow = 3,
                             labels = c("Cape Adare", "Cape Crozier", "Cape Hallett"),
                             label.x = 0.1))
annotate_figure(All_3_plot, left = "Log guano area (m2)", bottom = "Days since December 1st")

# View on the same plot

All_together <- ggplot(Dataset.1.4, aes(x = Day_D1, y = Log_GA, group = Colony_code)) +
  geom_point(aes(colour = Colony_code)) +
  geom_line(aes(colour = Colony_code)) +
  xlab("Days since December 1st") +
  ylab("Log guano area (m2)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_color_manual(values = colours) +
  scale_x_continuous(limits = c(1,90), breaks = seq(1,90, by=10))

All_together


# Linear model of log guano areas and days since December 1st
LM_2 <- lm(Dataset.1.4$Log_GA ~ Dataset.1.4$Day_D1)
summary(LM_2)

# Plot residuals
qqnorm(resid(LM_2))  # Q-Q plot for residuals
qqline(resid(LM_2))  # reference line
# pattern in residuals

par(mfrow = c(2,2))
plot(LM_2)

Dataset.1.4$fitted2 <- LM_2$fitted.values
Dataset.1.4$resid2 <- LM_2$residuals

Resids_2 <- ggplot(Dataset.1.4, aes(x=fitted2, y=resid2, colour = colony, fill = colony, shape = colony)) + 
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
  scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2, by=0.5))

Resids_2

# Heteroskedasticity in the residuals looks better
# Need to account for colony

# Linear mixed effect model of log guano areas and days since December 1st
# Including a random effect of colony
LMM_2 <- lme4::lmer(Log_GA ~ Day_D1 + (1|Colony_code), data = Dataset.1.4)
summary(LMM_2)
anova(LMM_2) # just looks at fixed effects

# Plot residuals
Dataset.1.4$fitted2.5 <- fitted(LMM_2)
Dataset.1.4$resid2.5 <- resid(LMM_2)

qqnorm(Dataset.1.4$resid2.5)  # Q-Q plot for residuals
qqline(Dataset.1.4$resid2.5)  # reference line

Resids_2.5 <- ggplot(Dataset.1.4, aes(x=fitted2.5, y=resid2.5, colour = colony, fill = colony, shape = colony)) + 
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

Resids_2.5

# Heteroskedasticity still present in the residuals

##########################################################################
# Appropriate models
##########################################################################

# Include colony as fixed effect
LM_3 <- lm(Log_GA ~ Day_D1 + Colony_code, data = Dataset.1.4)
summary(LM_3)
anova(LM_3)

# check residuals
# Plot residuals
qqnorm(resid(LM_3))  # Q-Q plot for residuals
qqline(resid(LM_3))  # reference line

par(mfrow = c(2,2))
plot(LM_3)

Dataset.1.4$fitted3 <- LM_3$fitted.values
Dataset.1.4$resid3 <- LM_3$residuals

Resids_3 <- ggplot(Dataset.1.4, aes(x=fitted3, y=resid3, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
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
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4,0.4,by=0.2))

Resids_3
# Heteroskedasticity present in the residuals (spreads out)

# Include colony as a random effect

library(nlme) # nlme to explore AR structure

LMM_3 <- lme(
  fixed = Log_GA ~ Day_D1,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(LMM_3)
anova(LMM_3)

# check residuals
# Plot residuals
qqnorm(resid(LMM_3))  # Q-Q plot for residuals
qqline(resid(LMM_3))  # reference line

Dataset.1.4$fittedLMM_3 <- fitted(LMM_3)
Dataset.1.4$residLMM_3 <- resid(LMM_3)

LMM_3_resids <- ggplot(Dataset.1.4, aes(x=fittedLMM_3, y=residLMM_3, colour = Colony_code, fill = Colony_code, shape = Colony_code)) + 
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

LMM_3_resids
# Heteroskedasticity

# Try AR(1) - auto regressive error structure (points will be correlated to the previous points)
# with random effect of colony
AR1_Model <- nlme::lme(
  fixed = Log_GA ~ Day_D1, # fixed effect
  random = ~ 1 | Colony_code, # random effect
  data = Dataset.1.4, 
  correlation = corAR1(form = ~ Day_D1 | Colony_code) # AR(1) correlation within colonies over time
  ) 

summary(AR1_Model) # the same as LMM_3
# Phi1 = 0 (means there is no autocorrelation)
anova(AR1_Model) # the same as LMM_3

# phi1 is zero indicating that auto regressive error structure isn't needed
# no auto regressive correlation is detected 

anova(LMM_3, AR1_Model)
# model without autoregressive structure is better 
# no evidence that adding the AR(1) structure improves model fit

# check residuals
# Plot residuals
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
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4, 0.4, by=0.2))

ResidsAR
# Heteroskedasticity

acf(resid(AR1_Model))
# no autocorrelation

##########################################################################

# Two-way anova
# Include interaction between Day_D1 and colony to see if the effect of Day_D1 depends on colony
LM_4 <- aov(Log_GA ~ Day_D1 + Colony_code + Day_D1 * Colony_code, data = Dataset.1.4)
summary(LM_4)

# Plot residuals
qqnorm(resid(LM_4))  # Q-Q plot for residuals
qqline(resid(LM_4))  # reference line

Dataset.1.4$fitted4 <- LM_4$fitted.values
Dataset.1.4$resid4 <- LM_4$residuals

Resids_4 <- ggplot(Dataset.1.4, aes(x=fitted4, y=resid4, colour = colony, fill = colony, shape = colony)) + 
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
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4,0.4, by=0.2))

Resids_4

# Linear regression (equivalent to model above)
# Include interaction between Day_D1 and colony to see if the effect of Day_D1 depends on colony
LM_5 <- lm(Log_GA ~ Day_D1 + Colony_code + Day_D1 * Colony_code, data = Dataset.1.4)
summary(LM_5)
anova(LM_5) # same as model above (LM_4)

# low sample size - 15 observations
# may need to do resampling due to low sample size (bootstrapping or other resampling method?)
# or remove interaction? (LM_3)

# Plot residuals
qqnorm(resid(LM_5))  # Q-Q plot for residuals
qqline(resid(LM_5))  # reference line

Dataset.1.4$fitted5 <- LM_5$fitted.values
Dataset.1.4$resid5 <- LM_5$residuals

Resids_5 <- ggplot(Dataset.1.4, aes(x=fitted5, y=resid5, colour = colony, fill = colony, shape = colony)) + 
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
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4,0.4, by=0.2))

Resids_5

# does the interaction improve the model
anova(LM_3, LM_5)
# the interaction doesn't significantly improve the model
# but it does reduce some unexplained variance
# LM 3 has better Q-Q plot


##########################################################################

# Differences between and within colonies
# Means and standard deviations

# Plot as box plot
Box_plot <- ggplot(Dataset.1.4, aes(x = Colony_code, y = Log_GA, fill = Colony_code)) +
  geom_boxplot() +
  xlab("Colony") +
  ylab("Log guano area (m²)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_manual(values = colours)

Box_plot


# Calculate means and standard deviations for each colony
library(dplyr)

# Table
summary_table <- Dataset.1.4 %>%
  group_by(Colony_code, Season) %>%
  summarise(n = n(),
            mean_GA = mean(GA),
            median_GA = median(GA),
            sd = sd(GA))

# Create a table plot
table_plot <- ggtexttable(summary_table, rows = NULL)

# Show table in Plots window
plot(table_plot)

# The larger the colony, the larger the variation
# Try removing shadowed estimates (Adare 2020 4 and 5)

##########################################################################

# Remove rows where Adare guano stain is over 500,000 m2
# Might be a better way to do this
Dataset.1.5 <- Dataset.1.4 %>% filter(GA <= 500000)
View(Dataset.1.5)

# Table without shadowed 
summary_table2 <- Dataset.1.5 %>%
  group_by(Colony_code, Season) %>%
  summarise(n = n(),
            mean_GA = mean(GA),
            median_GA = median(GA),
            sd = sd(GA))

# Create a table plot
table_plot2 <- ggtexttable(summary_table2, rows = NULL)

# Show table in Plots window
plot(table_plot2)

##########################################################################

# Produce non-linear growth model with nls()

# model log GA or just GA?
# non-linear mixed effect model (nlmm)

# needs to handles curves that vary by groups (colony)
# use logistic3 for a 3 parameter logistic growth model

# Fit a 3-parameter logistic model grouped by Colony_code

# no GA above log 14?
# initial day 1?
# growth rate 1?
# inflection point at day 70?
# growth fastest at day 50?
# or build manually in python 

##########################################################################
# Build model withe satellite-related factors
##########################################################################

# Investigate correlation between covariates using correlogram
library(corrplot)

covariates <- data.frame(Dataset.1.4$Day_D1, Dataset.1.4$MEANCOLLECTEDGSD,Dataset.1.4$MEANSUNEL, 
                         Dataset.1.4$MEANSUNAZ, Dataset.1.4$MEANOFFNADIRVIEWANGLE)
# cant include colony code (not numeric)

cor.matrix <- cor(covariates) # default method is pearsons

corrplot(cor.matrix, method = "number", type = "lower", tl.cex = 1)
# plot

# correlated variables of greater than +/- 0.60
# sun elevation and days since December first correlated - expected
# check variables with vif analysis too
# remove sun elevation?

library(car) # for VIF values
library(MuMIn) # for AICc scores

# base model without colony code
# no random effect of colony
base_model <- lm(Dataset.1.4$GA ~ Dataset.1.4$Day_D1 + Dataset.1.4$MEANCOLLECTEDGSD + 
                   Dataset.1.4$MEANSUNAZ + Dataset.1.4$MEANSUNEL +
                   Dataset.1.4$MEANOFFNADIRVIEWANGLE)
# too many parameters for observations?
summary(base_model)
vif(base_model) # remove sun elevation

# check residuals?
par(mfrow = c(2,2))
plot(base_model)

# Try guano area on log scale
Dataset.1.4$Log_GA <- log(Dataset.1.4$GA)

log_base_model <- lm(Dataset.1.4$Log_GA ~ Dataset.1.4$Day_D1 + Dataset.1.4$MEANCOLLECTEDGSD + 
                   Dataset.1.4$MEANSUNAZ + Dataset.1.4$MEANSUNEL +
                   Dataset.1.4$MEANOFFNADIRVIEWANGLE)

summary(log_base_model)
plot(log_base_model)
vif(log_base_model) # remove sun elevation

# removed sun elevation no random effect
M1_model <- lm(Dataset.1.4$Log_GA ~ Dataset.1.4$Day_D1 + Dataset.1.4$MEANCOLLECTEDGSD + 
                 Dataset.1.4$MEANSUNAZ + Dataset.1.4$MEANOFFNADIRVIEWANGLE)
summary(M1_model) # bad model
vif(M1_model) # all look fine
plot(M1_model)

# include random effect of colony in base_model
library(nlme)

# base model (without sun el)
log_base_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANCOLLECTEDGSD + MEANSUNAZ + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_base_lmm)
# GSD and sun AZ correlated
anova(log_base_lmm)
# days since December first may be significant
# colonies of different sizes are significantly different
vif(log_base_lmm)
# might need to remove sun azimuth too
# check residuals

# remove sun az
log_reduced_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANCOLLECTEDGSD + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced_lmm)
# days since December first significant
anova(log_reduced_lmm)
# colonies are different
vif(log_reduced_lmm) # all look fine
AICc(log_reduced_lmm)

# remove GSD?
log_reduced2_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced2_lmm)
anova(log_reduced2_lmm)
M2_AIC <- AIC(log_reduced2_lmm)
M2_AICc <- AICc(log_reduced2_lmm)

# Day_d1 model remove off-nadir
log_reduced3_lmm <- lme(
  fixed = Log_GA ~ Day_D1,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced3_lmm)
anova(log_reduced3_lmm)
M3_AIC <- AIC(log_reduced3_lmm)
M3_AICc <- AICc(log_reduced3_lmm)

# keep off-nadir instead
log_reduced4_lmm <- lme(
  fixed = Log_GA ~ MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced4_lmm)
anova(log_reduced4_lmm)
# better log lik and AIC
M4_AIC <- AIC(log_reduced4_lmm)
M4_AICc <- AICc(log_reduced4_lmm)

# keep gsd instead
log_reduced5_lmm <- lme(
  fixed = Log_GA ~ MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced5_lmm)
anova(log_reduced5_lmm)
# no different to off-nadir model
M5_AIC <- AIC(log_reduced5_lmm)
M5_AICc <- AICc(log_reduced5_lmm)

# sun az model
log_reduced6_lmm <- lme(
  fixed = Log_GA ~ MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced6_lmm)
anova(log_reduced6_lmm)
# quite bad
M6_AIC <- AIC(log_reduced6_lmm)
M6_AICc <- AICc(log_reduced6_lmm)

# sun el model
log_reduced7_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced7_lmm)
anova(log_reduced7_lmm)
# best model
M7_AIC <- AIC(log_reduced7_lmm)
M7_AICc <- AICc(log_reduced7_lmm)

# gsd and off-nadir model
log_reduced8_lmm <- lme(
  fixed = Log_GA ~ MEANCOLLECTEDGSD + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced8_lmm)
anova(log_reduced8_lmm)
M8_AIC <- AIC(log_reduced8_lmm)

# sun el and off-nadir model
log_reduced9_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced9_lmm)
anova(log_reduced9_lmm)
M9_AIC <- AIC(log_reduced9_lmm)

# sun el and gsd model
log_reduced10_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced10_lmm)
anova(log_reduced10_lmm)
M10_AIC <- AIC(log_reduced10_lmm)

# sun el and sun az model
log_reduced11_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced11_lmm)
anova(log_reduced11_lmm)

# sun az and off-nadir model
log_reduced12_lmm <- lme(
  fixed = Log_GA ~ MEANSUNAZ + MEANOFFNADIRVIEWANGLE,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced12_lmm)
anova(log_reduced12_lmm)

# Day_D1 and gsd
log_reduced13_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced13_lmm)
anova(log_reduced13_lmm)

# Day_D1 and sun az
log_reduced14_lmm <- lme(
  fixed = Log_GA ~ Day_D1 + MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced14_lmm)
anova(log_reduced14_lmm)

# GSD and sun az
log_reduced15_lmm <- lme(
  fixed = Log_GA ~ MEANCOLLECTEDGSD + MEANSUNAZ,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
summary(log_reduced15_lmm)
anova(log_reduced15_lmm)

# *order of best*
# *5 individual models*
# sun el 7
# off-nadir 4
# gsd 5
# day_d1 3
# sun az 6

# sun el and off_nadir9
# sun el and gsd 10
# day_d1 and off-nadir 2
# gsd and off-nadir 8
# sun el and sun az 11
# day_d1 and gsd 13
# gsd and sun az (correlated) 15
# sun az and off-nadir 12
# day_d1 and sun az 14

# need to use subsets of full model excluding collinear variables

# all three competitive covariates
log_reduced16_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL + MEANOFFNADIRVIEWANGLE + MEANCOLLECTEDGSD,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4)
summary(log_reduced16_lmm)
vif(log_reduced16_lmm) # but under 2
# not great - correlations

# calculate delta AICc scores and weights
library(qpcR) # package masks MuMin 

null_model <- lme(
  fixed = Log_GA ~ 1,
  random = ~ 1 | Colony_code,
  data = Dataset.1.4
)
null_AIC <- AIC(null_model)
# null model better?
# confirms random effect is needed only?

# for candidate models AIC
x <- c(M3_AIC,M4_AIC,M5_AIC,M6_AIC,M7_AIC)
akaike.weights(x)

# AICc
x <- c(M3_AICc,M4_AICc,M5_AICc,M6_AICc,M7_AICc)
akaike.weights(x)
# same story

# stepAIC?

# variation in satellite related factors
# is more influential then a change in size

# check residuals
# Plot residuals
qqnorm(resid(log_reduced7_lmm))  # Q-Q plot for residuals
qqline(resid(log_reduced7_lmm))  # reference line

Dataset.1.4$fittedbest <- fitted(log_reduced7_lmm)
Dataset.1.4$residbest <- resid(log_reduced7_lmm)

library(ggplot2)

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

# acf on sun el model
AR1_lmm <- lme(
  fixed = Log_GA ~ MEANSUNEL,
  random = ~ 1 | Colony_code,
  correlation = corAR1(form = ~ Day_D1 | Colony_code),
  data = Dataset.1.4
)
summary(AR1_lmm) # phi 0
anova(AR1_lmm)

acf(resid(AR1_lmm))
# no autocorrelation

anova(log_reduced7_lmm, AR1_lmm)
# check which one is better


##########################################################################
# look at shadows
##########################################################################
# out of 15 there are only two shadowed estimates (all others 0)

# Extract only within season data and needed variables 
Dataset.3.1 <- Dataset.1.0[,c("Colony_code","Analysis","Season","Date",
                              "GA","MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE","Shadow")]

# Keep only intra-seaonal data
Dataset.3.1$Analysis[Dataset.3.1$Analysis == ""] <- NA
Dataset.3.1$Analysis[Dataset.3.1$Analysis == "NA"] <- NA
sum(is.na(Dataset.3.1$Analysis))

# Remove NAs
Dataset.3.2 <- na.omit(Dataset.3.1)

# Condense/ remove analysis column
Dataset.3.3 <- Dataset.3.2[,c("Colony_code","Season","Date","GA",
                              "MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE","Shadow")]

# View
View(Dataset.3.3)

# Turn image date into days since December 1st (Day_D1)

# convert image dates to date variable in r
Dataset.3.3$r_date <- as.character(Dataset.3.3$Date)
Dataset.3.3$r_date <- as.Date(Dataset.3.3$r_date, format="%Y%m%d")

# Calculate day deviation since December 1st
Dataset.3.3$Day_D1 <- as.numeric(difftime(Dataset.3.3$r_date,
                                          as.Date(paste0(
                                            ifelse
                                            (format(Dataset.3.3$r_date, "%m") == "12",
                                              format(Dataset.3.3$r_date, "%Y"),
                                              as.numeric(format(Dataset.3.3$r_date, "%Y")) 
                                              -1), "-12-01")), units = "days"))
# Days since December 1st: numerical but discrete without decimals 

# Condense/ remove r date column
Dataset.3.4 <- Dataset.3.3[,c("Colony_code","Season","Date","GA", "Day_D1",
                              "MEANCOLLECTEDGSD","MEANSUNAZ","MEANSUNEL","MEANOFFNADIRVIEWANGLE","Shadow")]

Dataset.3.4$Log_GA <- log(Dataset.3.4$GA)

# View
View(Dataset.3.4)


covariates2 <- data.frame(Dataset.3.4$Day_D1, Dataset.3.4$MEANCOLLECTEDGSD,Dataset.3.4$MEANSUNEL, 
                         Dataset.3.4$MEANSUNAZ, Dataset.3.4$MEANOFFNADIRVIEWANGLE,Dataset.3.4$Shadow)
# cant include colony code (not numeric)

cor.matrix2 <- cor(covariates2) # default method is pearsons

corrplot(cor.matrix2, method = "number", type = "lower", tl.cex = 1)

# Shadow
Shadow_lmm <- lme(
  fixed = Log_GA ~ Shadow,
  random = ~ 1 | Colony_code,
  data = Dataset.3.4
)
summary(Shadow_lmm)
anova(Shadow_lmm)

##########################################################################
# GA and BP relationship
##########################################################################

# Extract only within season data and needed variables 
Dataset.2 <- Dataset.1.0[,c("Colony_name","GA","BP","Analysis2")]
View(Dataset.2)

# Keep only interseaonal data
# interseasonal analysis for within season images is median estimate date
Dataset.2$Analysis2[Dataset.2$Analysis2 == ""] <- NA
Dataset.2$Analysis2[Dataset.2$Analysis2 == "NA"] <- NA
sum(is.na(Dataset.2$Analysis2))

# Condense/ remove analysis column

# Remove NAs
Dataset.2.1 <- na.omit(Dataset.2)
View(Dataset.2.1)

Dataset.2.1$Log_GA <- log(Dataset.2.1$GA)
Dataset.2.1$Log_BP <- log(Dataset.2.1$BP)

GA_BP <- lm(Dataset.2.1$Log_GA ~ Dataset.2.1$Log_BP)
summary(GA_BP)
r2 <- round(summary(GA_BP)$r.squared, 2)

colours <- c("orchid","navy", "indianred", "royalblue", "skyblue", "red", 
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
  scale_color_manual(values = colours)

All_plot # does this include all the intra-seasonal estimates?
