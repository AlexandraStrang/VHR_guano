# r script for exploratory analysis of within seasonal guano stain variation
# creator: Alexandra Strang
# created: 2025

# Read in data
Dataset.1.0 <- read.csv("Strang_PhD_masterdata.csv")

# Extract only within season data and needed variables 
Dataset.1.1 <- Dataset.1.0[,c("Colony_code","Analysis","Season","Date","GA")]

# Keep only intraseaonal data
Dataset.1.1$Analysis[Dataset.1.1$Analysis == ""] <- NA
Dataset.1.1$Analysis[Dataset.1.1$Analysis == "NA"] <- NA
sum(is.na(Dataset.1.1$Analysis))

# Remove NAs
Dataset.1.2 <- na.omit(Dataset.1.1)

# Condense/ remove analysis column
Dataset.1.3 <- Dataset.1.2[,c("Colony_code","Season","Date","GA")]

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
Dataset.1.4 <- Dataset.1.3[,c("Colony_code","Season","Date","GA", "Day_D1")]

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

# Pools data across colonies but colonies vary in size
# Need to include the influence of colony
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
  xlab("Days since December 1st") +
  ylab("Log guano area (m2)") +
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
  xlab("Days since December 1st") +
  ylab("Log guano area (m2)") +
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
  xlab("Days since December 1st") +
  ylab("Log guano area (m2)") +
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

# View on the same plot
colours <- c("darkblue","royalblue","skyblue")

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

# Include colony as fixed effect
LM_3 <- lm(Log_GA ~ Day_D1 + Colony_code, data = Dataset.1.4)
summary(LM_3)
anova(LM_3)

# Plot residuals
qqnorm(resid(LM_3))  # Q-Q plot for residuals
qqline(resid(LM_3))  # reference line

Dataset.1.4$fitted3 <- LM_3$fitted.values
Dataset.1.4$resid3 <- LM_3$residuals

Resids_3 <- ggplot(Dataset.1.4, aes(x=fitted3, y=resid3, colour = colony, fill = colony, shape = colony)) + 
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

Resids_3

# Heteroskedasticity still present in the residuals

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


# Model by colony individually
# poor

# Adare
LM_A <- lm(ADARdf2$Log_GA ~ ADARdf2$Day_D1)
summary(LM_A)

# Crozier
LM_C <- lm(CROZdf2$Log_GA ~ CROZdf2$Day_D1)
summary(LM_C)

# Hallett
LM_H <- lm(HALLdf2$Log_GA ~ HALLdf2$Day_D1)
summary(LM_H)

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


# Produce non-linear growth model with nls() or nlgm() 

nls()

library(nlgm)

nlgm

# or build manually in python 
