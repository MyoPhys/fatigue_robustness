require(pracma)
library(signal)
library(pspline)
library(lme4)
library(emmeans)
library(ggplot2)

# Import the data
filename <- file.choose()
data <- (read.csv(filename, sep = ",", skip=0))
data$B.A <- factor(data$B.A , levels=c("B", "A"))

# Separate data by treatment and before/after
data20 <- data[which(data$Load=="LAMSA20"),]
data20$B.A <- factor(data20$B.A , levels=c("B", "A"))
data60 <- data[which(data$Load=="LAMSA60"),]
data60$B.A <- factor(data60$B.A , levels=c("B", "A"))
data60B <- data60[which(data60$B.A=="B"),]
data60A <- data60[which(data60$B.A=="A"),]
data20B <- data20[which(data20$B.A=="B"),]
data20A <- data20[which(data20$B.A=="A"),]

# Descriptive Statistics - replace variable of interest
mean_data60B <- mean(data60B$tendonp)
SEM_data60B <- std_err(data60B$tendonp)
mean_data60A <- mean(data60A$tendonp)
SEM_data60A <- std_err(data60A$tendonp)
mean_data20B <- mean(data20B$tendonp)
SEM_data20B <- std_err(data20B$tendonp)
mean_data20A <- mean(data20A$tendonp)
SEM_data20A <- std_err(data20A$tendonp)

c(mean_data60B,SEM_data60B,mean_data60A,SEM_data60A,mean_data20B,SEM_data20B,mean_data20A,SEM_data20A)


# Preloaded Elastic Energy
aovPre <- aov( data$mms.pre ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.pre, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovPre)

# Overshoot Elastic Energy
aovOvershoot <- aov( data$mms.over ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.over, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovOvershoot)

# Loading Muscle Power
aovLoading <- aov( data$loading_power ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, loading_power, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovLoading)

# Effective Stiffness
aovStiffness <- aov( data$tendon.stiffness ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, tendon.stiffness, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovStiffness)

# Loading Time
aovLoadTime <- aov( data$loading_time ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, loading_time, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovLoadTime)

# Recoil Time
aovRecoilTime <- aov( data$recoil_time ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, recoil_time, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovRecoilTime)

# Elastic Energy Return
aovReturn <- aov( data$mms.return ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.return, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovReturn)

# Elastic Recoil Power
aovSEEPower <- aov( data$tendonp ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, tendonp, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovSEEPower)

# Elastic Efficiency
aovEfficiency <- aov( data$efficiency ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, efficiency, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovEfficiency)

# Direct Muscle Work
aovDirect <- aov( data$mms.direct ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, mms.direct, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovDirect)

# Direct Muscle Power
aovDirectP <- aov( data$fasciclep ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, fasciclep, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovDirectP)

# MTU Power
aovMTUP <- aov( data$MTUp ~ data$Load * data$B.A + Error(data$Indiv))
elasticbox <- ggplot(data, aes(Load, MTUp, fill = B.A))
elasticbox +  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75))
summary(aovMTUP)










