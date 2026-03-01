# r script for plotting GA and BP trends over time for 15 colonies
# creator: Alexandra Strang
# created: 2026

# Load packages
library(ggplot2)
library(ggpubr)

setwd("C:/Users/astra/OneDrive - University of Canterbury/ANTA - PhD/Data/Data sheets")

# Read in data
Dataset.1.0 <- read.csv("Merged_masterdata.csv")

# Extract only GA data first
Dataset.GA <- Dataset.1.0[,c("Colony_name","GA","Analysis2","Season")]
View(Dataset.GA)

# Keep only interseaonal GA data
# interseasonal analysis for within season images is median estimate date
Dataset.GA$Analysis2[Dataset.GA$Analysis2 == ""] <- NA
Dataset.GA$Analysis2[Dataset.GA$Analysis2 == "NA"] <- NA
sum(is.na(Dataset.GA$Analysis2))

# Remove NAs
Dataset.GA.1 <- na.omit(Dataset.GA)
View(Dataset.GA.1) # all 15 colonies should have at least 3 GA estimates

# Remove Analysis column
Dataset.GA.2 <- Dataset.GA.1[,c("Colony_name","GA","Season")]
View(Dataset.GA.2)

# Log GA
Dataset.GA.2$Log_GA <- log(Dataset.GA.2$GA)

# Extract BP only data next
Dataset.BP <- Dataset.1.0[,c("Colony_name","BP","Analysis2","Season")]
View(Dataset.BP)

# Keep only interseaonal BP data
# interseasonal analysis for within season images is median estimate date
Dataset.BP$Analysis2[Dataset.BP$Analysis2 == ""] <- NA
Dataset.BP$Analysis2[Dataset.BP$Analysis2 == "NA"] <- NA
sum(is.na(Dataset.BP$Analysis2))

# Remove NAs
Dataset.BP.1 <- na.omit(Dataset.BP)
View(Dataset.BP.1)

# Remove Analysis column
Dataset.BP.2 <- Dataset.BP.1[,c("Colony_name","BP","Season")]
View(Dataset.BP.2)

# Log BP
Dataset.BP.2$Log_BP <- log(Dataset.BP.2$BP)

# Plot each individually and then combine

# Subset data by colony
BEAU_GA <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Beaufort Island")
BEAU_BP <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Beaufort Island")

BEA_N_GA <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Beaufort Island North")
BEA_N_BP <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Beaufort Island North")
# BEAU_N from datasets

ADAR_GA <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Cape Adare")
ADAR_BP <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Cape Adare")

BIRD_GA <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Cape Bird")
BIRD_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Cape Bird")

CROZ_GA <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Cape Crozier")
CROZ_BP <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Cape Crozier")

HALL_GA <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Cape Hallett")
HALL_BP <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Cape Hallett")

ROYD_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Cape Royds")
ROYD_BP <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Cape Royds")

COU_M_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Coulman Island Middle")
COU_M_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Coulman Island Middle")

COU_N_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Coulman Island North")
COU_N_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Coulman Island North")

COU_S_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Coulman Island South")
COU_S_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Coulman Island South")

DOWN_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Downshire Cliffs")
DOWN_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Downshire Cliffs")

DUKE_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Duke of York Island")
DUKE_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Duke of York Island")

INEX_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Inexpressible Island")
INEX_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Inexpressible Island")

POSS_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Possession Island")
POSS_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Possession Island")

SVEN_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Sven Foyn Island")
SVEN_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Sven Foyn Island")

TERR_GA  <- subset(Dataset.GA.2, Dataset.GA.2$Colony_name=="Terra Nova Bay")
TERR_BP  <- subset(Dataset.BP.2, Dataset.BP.2$Colony_name=="Terra Nova Bay")


# Plot Beaufort Island
BEAU_trends_plot <- ggplot(BEAU_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "orchid") + 
  geom_line(colour = "orchid") +
  xlab("") +
  ylab("") +
  labs(title = "Beaufort Island") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

BEAU_trends_plot

# Plot Beaufort Island North
BEA_N_trends_plot <- ggplot(BEA_N_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "blue") + 
  geom_line(colour = "blue") +
  xlab("") +
  ylab("") +
  labs(title = "Beaufort Island North") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

BEA_N_trends_plot

# Plot Cape Adare
ADAR_trends_plot <- ggplot(ADAR_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "navy") + 
  geom_line(colour = "navy") +
  xlab("") +
  ylab("") +
  labs(title = "Cape Adare") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

ADAR_trends_plot

# Plot Cape Bird
BIRD_trends_plot <- ggplot(BIRD_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "indianred") + 
  geom_line(colour = "indianred") +
  xlab("") +
  ylab("") +
  labs(title = "Cape Bird") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

BIRD_trends_plot
# 2009 seems wrong

# Plot Cape Crozier
CROZ_trends_plot <- ggplot(CROZ_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "royalblue") + 
  geom_line(colour = "royalblue") +
  xlab("") +
  ylab("") +
  labs(title = "Cape Crozier") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

CROZ_trends_plot

# Plot Cape Hallett
HALL_trends_plot <- ggplot(HALL_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "royalblue") + 
  geom_line(colour = "royalblue") +
  xlab("") +
  ylab("") +
  labs(title = "Cape Hallett") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

HALL_trends_plot

# Plot Cape Royds
ROYD_trends_plot <- ggplot(ROYD_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "red") + 
  geom_line(colour = "red") +
  xlab("") +
  ylab("") +
  labs(title = "Cape Royds") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

ROYD_trends_plot

# Plot Coulman Island Middle
COU_M_trends_plot <- ggplot(COU_M_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "pink") + 
  geom_line(colour = "pink") +
  xlab("") +
  ylab("") +
  labs(title = "Coulman Island Middle") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

COU_M_trends_plot

# Plot Coulman Island North
COU_N_trends_plot <- ggplot(COU_N_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "hotpink") + 
  geom_line(colour = "hotpink") +
  xlab("") +
  ylab("") +
  labs(title = "Coulman Island North") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

COU_N_trends_plot

# Plot Coulman Island South
COU_S_trends_plot <- ggplot(COU_S_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "orange") + 
  geom_line(colour = "orange") +
  xlab("") +
  ylab("") +
  labs(title = "Coulman Island South") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

COU_S_trends_plot

# Plot Downshire cliffs
DOWN_trends_plot <- ggplot(DOWN_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "gold") + 
  geom_line(colour = "gold") +
  xlab("") +
  ylab("") +
  labs(title = "Downshire Cliffs") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

DOWN_trends_plot

# Plot Duke of York Island
DUKE_trends_plot <- ggplot(DUKE_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "slateblue") + 
  geom_line(colour = "slateblue") +
  xlab("") +
  ylab("") +
  labs(title = "Duke of York Island") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

DUKE_trends_plot

# Plot Inexpressible Island
INEX_trends_plot <- ggplot(INEX_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "magenta") + 
  geom_line(colour = "magenta") +
  xlab("") +
  ylab("") +
  labs(title = "Inexpressible Island") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

INEX_trends_plot

# Plot Possession Island
POSS_trends_plot <- ggplot(POSS_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "cyan") + 
  geom_line(colour = "cyan") +
  xlab("") +
  ylab("") +
  labs(title = "Possession Island") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

POSS_trends_plot

# Plot Sven Foyn Island
SVEN_trends_plot <- ggplot(SVEN_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "grey") + 
  geom_line(colour = "grey") +
  xlab("") +
  ylab("") +
  labs(title = "Sven Foyn Island") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

SVEN_trends_plot

# Plot Terra Nova Bay
TERR_trends_plot <- ggplot(TERR_GA, aes(x = Season, y = GA)) + 
  geom_point(colour = "violetred") + 
  geom_line(colour = "violetred") +
  xlab("") +
  ylab("") +
  labs(title = "Terra Nova Bay") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_continuous(limits = c(2009,2024), breaks = seq(2009,2024, by=2))

TERR_trends_plot

# 16 colonies
Together <- plot(ggarrange(BEAU_trends_plot,
                           BEA_N_trends_plot,
                           ADAR_trends_plot,
                           BIRD_trends_plot,
                           CROZ_trends_plot,
                           HALL_trends_plot,
                           ROYD_trends_plot,
                           COU_M_trends_plot,
                           COU_N_trends_plot,
                           COU_S_trends_plot,
                           DOWN_trends_plot,
                           DUKE_trends_plot,
                           INEX_trends_plot,
                           POSS_trends_plot, 
                           SVEN_trends_plot,
                           TERR_trends_plot,
                           ncol = 4, nrow = 4))
annotate_figure(Together, left = "Guano area (m²)", bottom = "Season")

# To add BP from aerial census need to get BP from GA