require(pracma)
library(signal)
library(pspline)
library(lme4)
library(emmeans)
library(ggplot2)

# Import the data
filename <- file.choose()
data <- (read.csv(filename, sep = ",", skip=0))
data$B.A <- factor(data$B.A , levels=c("Before", "After"))

# Separate the data by treatment
data20 <- data[which(data$Load=="Iso20"),]
data20$B.A <- factor(data20$B.A , levels=c("Before", "After"))

data60 <- data[which(data$Load=="Iso60"),]
data60$B.A <- factor(data60$B.A , levels=c("Before", "After"))

data60B <- data60[which(data60$B.A=="Before"),]
data60A <- data60[which(data60$B.A=="After"),]

data20B <- data20[which(data20$B.A=="Before"),]
data20A <- data20[which(data20$B.A=="After"),]

# Descriptive Statistics - replace variable of interest
mean_data60B <- mean(data60B$time.to.short)
SEM_data60B <- std_err(data60B$time.to.short)
mean_data60A <- mean(data60A$time.to.short)
SEM_data60A <- std_err(data60A$time.to.short)
mean_data20B <- mean(data20B$time.to.short)
SEM_data20B <- std_err(data20B$time.to.short)
mean_data20A <- mean(data20A$time.to.short)
SEM_data20A <- std_err(data20A$time.to.short)

c(mean_data60B,SEM_data60B,mean_data60A,SEM_data60A,mean_data20B,SEM_data20B,mean_data20A,SEM_data20A)

# Muscle Power
aovPower <- aov(data$mms.p ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.p, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovPower)

# Effective Stiffness of SEE
aovStiffness <- aov(data$tendon.stiffness ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, tendon.stiffness, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovStiffness)

# Energy Stored in SEE
aovSEE <- aov(data$mms.elastic ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.elastic, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovSEE)

# Loading Time
aovTime <- aov(data$time.to.short ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, time.to.short, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovTime)

# Direct Muscle Work
aovMuscle <- aov(data$mms.muscle ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.muscle, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovMuscle)










