require(pracma)
library(signal)
library(pspline)
library(lme4)
library(emmeans)

# Import the data
filename <- file.choose()
data <- (read.csv(filename, sep = ",", skip=0))

# Separate data by treatment
data20 <- data[which(data$Load=="LAMSA20"),]
data20$B.A <- factor(data20$B.A , levels=c("B", "A"))

data60 <- data[which(data$Load=="LAMSA60"),]
data60$B.A <- factor(data60$B.A , levels=c("B", "A"))

data60B <- data60[which(data60$B.A=="B"),]
data60A <- data60[which(data60$B.A=="A"),]

data20B <- data20[which(data20$B.A=="B"),]
data20A <- data20[which(data20$B.A=="A"),]

# Descriptive Statistics - replace variable of interest
mean_data60B <- mean(data60B$MTUp)
SEM_data60B <- std_err(data60B$MTUp)
mean_data60A <- mean(data60A$MTUp)
SEM_data60A <- std_err(data60A$MTUp)
mean_data20B <- mean(data20B$MTUp)
SEM_data20B <- std_err(data20B$MTUp)
mean_data20A <- mean(data20A$MTUp)
SEM_data20A <- std_err(data20A$MTUp)

c(mean_data60B,SEM_data60B,mean_data60A,SEM_data60A,mean_data20B,SEM_data20B,mean_data20A,SEM_data20A)


# Muscle Work
ttest60 <- t.test(data60B$mms.direct, data60A$mms.direct, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.direct, data20A$mms.direct, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.direct ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,15),outline = FALSE)
points(data20$B.A, data20$mms.direct, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.direct ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,15),outline = FALSE)
points(data60$B.A, data60$mms.direct, pch = c(1,2,3,4,5,1,2,3,4,5))

# Pre Load Elastic
ttest60 <- t.test(data60B$mms.pre, data60A$mms.pre, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.pre, data20A$mms.pre, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.pre ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,2),outline = FALSE)
points(data20$B.A, data20$mms.pre, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.pre ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,10),outline = FALSE)
points(data60$B.A, data60$mms.pre, pch = c(1,2,3,4,5,1,2,3,4,5))

# overshoot
ttest60 <- t.test(data60B$mms.over, data60A$mms.over, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.over, data20A$mms.over, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.over ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,25),outline = FALSE)
points(data20$B.A, data20$mms.over, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.over ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,25),outline = FALSE)
points(data60$B.A, data60$mms.over, pch = c(1,2,3,4,5,1,2,3,4,5))

#efficiency
ttest60 <- t.test(data60B$efficiency, data60A$efficiency, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$efficiency, data20A$efficiency, paired = TRUE, alternative = "two.sided")

boxplot(data20$efficiency ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,150),outline = FALSE)
points(data20$B.A, data20$efficiency, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$efficiency ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,100),outline = FALSE)
points(data60$B.A, data60$efficiency, pch = c(1,2,3,4,5,1,2,3,4,5))

#stiffness
ttest60 <- t.test(data60B$tendon.stiffness, data60A$tendon.stiffness, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$tendon.stiffness, data20A$tendon.stiffness, paired = TRUE, alternative = "two.sided")

boxplot(data20$tendon.stiffness ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,10),outline = FALSE)
points(data20$B.A, data20$tendon.stiffness, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$tendon.stiffness ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,10),outline = FALSE)
points(data60$B.A, data60$tendon.stiffness, pch = c(1,2,3,4,5,1,2,3,4,5))

#tendonp
ttest60 <- t.test(data60B$tendonp, data60A$tendonp, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$tendonp, data20A$tendonp, paired = TRUE, alternative = "two.sided")

boxplot(data20$tendonp ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,2000),outline = FALSE)
points(data20$B.A, data20$tendonp, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$tendonp ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,2000),outline = FALSE)
points(data60$B.A, data60$tendonp, pch = c(1,2,3,4,5,1,2,3,4,5))

#fasciclep
ttest60 <- t.test(data60B$fasciclep, data60A$fasciclep, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$fasciclep, data20A$fasciclep, paired = TRUE, alternative = "two.sided")

boxplot(data20$fasciclep ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,500),outline = FALSE)
points(data20$B.A, data20$fasciclep, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$fasciclep ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,500),outline = FALSE)
points(data60$B.A, data60$fasciclep, pch = c(1,2,3,4,5,1,2,3,4,5))

#MTUp
ttest60 <- t.test(data60B$MTUp, data60A$MTUp, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$MTUp, data20A$MTUp, paired = TRUE, alternative = "two.sided")

boxplot(data20$MTUp ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,2000),outline = FALSE)
points(data20$B.A, data20$MTUp, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$MTUp ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,2000),outline = FALSE)
points(data60$B.A, data60$MTUp, pch = c(1,2,3,4,5,1,2,3,4,5))

#Return
ttest60 <- t.test(data60B$mms.return, data60A$mms.return, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$mms.return, data20A$mms.return, paired = TRUE, alternative = "two.sided")

boxplot(data20$mms.return ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,20),outline = FALSE)
points(data20$B.A, data20$mms.return, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$mms.return ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,20),outline = FALSE)
points(data60$B.A, data60$mms.return, pch = c(1,2,3,4,5,1,2,3,4,5))

#Recoil Time
ttest60 <- t.test(data60B$recoil_time, data60A$recoil_time, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$recoil_time, data20A$recoil_time, paired = TRUE, alternative = "two.sided")

boxplot(data20$recoil_time ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,200),outline = FALSE)
points(data20$B.A, data20$recoil_time, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$recoil_time ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,200),outline = FALSE)
points(data60$B.A, data60$recoil_time, pch = c(1,2,3,4,5,1,2,3,4,5))

#Loading Time
ttest60 <- t.test(data60B$loading_time, data60A$loading_time, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$loading_time, data20A$loading_time, paired = TRUE, alternative = "two.sided")

boxplot(data20$loading_time ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,200),outline = FALSE)
points(data20$B.A, data20$loading_time, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$loading_time ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,200),outline = FALSE)
points(data60$B.A, data60$loading_time, pch = c(1,2,3,4,5,1,2,3,4,5))

#Loading power
ttest60 <- t.test(data60B$loading_power, data60A$loading_power, paired = TRUE, alternative = "two.sided")
ttest20 <- t.test(data20B$loading_power, data20A$loading_power, paired = TRUE, alternative = "two.sided")

boxplot(data20$loading_power ~ data20$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,200),outline = FALSE)
points(data20$B.A, data20$loading_power, pch = c(1,2,3,4,5,1,2,3,4,5))

boxplot(data60$loading_power ~ data60$B.A, boxwex = 0.8, xlim = c(0.5,2.5),ylim = c(0,200),outline = FALSE)
points(data60$B.A, data60$loading_power, pch = c(1,2,3,4,5,1,2,3,4,5))




