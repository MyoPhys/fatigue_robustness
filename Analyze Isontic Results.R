require(pracma)
library(signal)
library(pspline)
library(lme4)
library(emmeans)

# Import the data
filename <- file.choose()
data <- (read.csv(filename, sep = ",", skip=0))

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

# Energy Stored in SEE
ttest60 <- t.test(data60B$mms.elastic, data60A$mms.elastic, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.elastic, data20A$mms.elastic, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.elastic ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,2),outline = FALSE)
points(data20$B.A, data20$mms.elastic, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.elastic ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,8),outline = FALSE)
points(data60$B.A, data60$mms.elastic, pch = c(1,2,3,4,5,1,2,3,4,5))

# Direct Muscle Work
ttest60 <- t.test(data60B$mms.muscle, data60A$mms.muscle, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.muscle, data20A$mms.muscle, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.muscle ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,20),outline = FALSE)
points(data20$B.A, data20$mms.muscle, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.muscle ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,45),outline = FALSE)
points(data60$B.A, data60$mms.muscle, pch = c(1,2,3,4,5,1,2,3,4,5))

# Tendon Stiffness
ttest60 <- t.test(data60B$tendon.stiffness, data60A$tendon.stiffness, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$tendon.stiffness, data20A$tendon.stiffness, paired = TRUE, alternative = "two.sided")

boxplot(data20$tendon.stiffness ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,8),outline = FALSE)
points(data20$B.A, data20$tendon.stiffness, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$tendon.stiffness ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,12),outline = FALSE)
points(data60$B.A, data60$tendon.stiffness, pch = c(1,2,3,4,5,1,2,3,4,5))

# Power
ttest60 <- t.test(data60B$mms.p, data60A$mms.p, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.p, data20A$mms.p, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.p ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,8),outline = FALSE)
points(data20$B.A, data20$mms.p, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.p ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,12),outline = FALSE)
points(data60$B.A, data60$mms.p, pch = c(1,2,3,4,5,1,2,3,4,5))

# Loading Time
ttest60 <- t.test(data60B$time.to.short, data60A$time.to.short, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$time.to.short, data20A$time.to.short, paired = TRUE, alternative = "two.sided")

boxplot(data20$time.to.short ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,60),outline = FALSE)
points(data20$B.A, data20$time.to.short, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$time.to.short ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,100),outline = FALSE)
points(data60$B.A, data60$time.to.short, pch = c(1,2,3,4,5,1,2,3,4,5))








