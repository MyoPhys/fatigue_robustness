require(pracma)
library(signal)
library(pspline)

# Edit the metadata for the trial being run
MTU_length <- 43

# Import the data. Select a .txt file from Raw Data
filename <- file.choose()
data <- (read.csv(filename, sep = "\t", skip=0))

# For all trials, the stimulation starts 1000 ms after the start of data collection. We start visualizing at 500ms
start <- 500;

# Trim the Force data to the start
Force <- data[start:4999,1] ;               #newtons (N)

# The system is set at L0, so there is passive tension. Adjust for passive tension
Force <- Force - Force[1];

# Calibration in Igor code was set to Azizi motor, this applies the calibration for Olberding motor
Force <- Force/4.84*5.27

# Set the end point to trim sono data based on the size of the force recording
end <- length(Force);

# Trim the Sono data to the start and end points
sono_SS <- data[start:end,4];                    #(mm) 
# Calibration in Igor file was set for Azizi sono; this applies calibration for Olberding sono.
sono_SS <- sono_SS/12.77*13.823-1.206;

# Trim the level position to the start and end points
Lout <- data[start:end,3]; #* 2.61;                #(mm)
# Calibration in Igor file was set for Azizi motor; this applies calibration for Olberding motor.
Lout <- Lout/2.61*2.052;

# Adjust for starting position of the lever so it starts at 0
Lout <- Lout - Lout[1];

# Save a copy of the raw Sono data
sono_RR <- sono_SS

# set a tolerance for outliers
tolerance <- 10*(sono_SS[1] - sono_SS[length(sono_SS)])/length(sono_SS) #the average velocity of the whole thing (m/s)
sono_diff <- NULL
sono_TT <- sono_SS
outliers <- NULL

# create index of outliers in Sono data
for(i in 2:length(sono_SS)){
  sono_diff[i] <- sono_TT[i] - sono_TT[i-1]
  if(abs(sono_diff[i]) > tolerance){
    outliers <- c(outliers, i)
  }
}

# correct each outlier
for(i in outliers){
  sono_TT[i] <- (sono_TT[i-1] + sono_TT[i+1])/2
}


# plot the raw data
plot(sono_RR)
points(sono_SS, col="blue")
# plot data without outliers
points(sono_TT, col="red")

#search for baseline shifts
meandif <- NULL
sono_WW <- sono_TT[1:(length(sono_TT)-1)]
for(i in 11:length(sono_TT)){
  meandif[i] <- mean(sono_WW[(i-10):i]) - mean(sono_WW[i:(i+10)])
}
for(i in 11:(length(meandif)-11)){
  if(meandif[i] < -1.9){
    localmin <- i - (10 - which.min(sono_WW[(i-10):i]))
    localmax <- i + (which.max(sono_WW[i:(i+10)])) 
    localdif <- sono_WW[localmax] - sono_WW[localmin]
    interplength <- length(sono_WW[(localmin+1):(localmax-1)])
    interpoints <- rep(sono_WW[(localmin+1)],interplength)
    sono_WW <- c(sono_WW[1:localmin],interpoints,sono_WW[localmax:length(sono_WW)]-localdif)
  }
}


# plot the raw data
plot(sono_WW)
# plot data without baseline shifts
points(sono_TT, col="red")


# Assigning the filter specification variables
freq <- 1000 # Frequency of data
freqN <- freq/2 # Frequency normalized to Nyquist frequency
passbandFreqN <- 10/freqN  # Passband frequency AKA Wp; normalized to Nyquist frequency
stopbandFreqN <- 500/freqN  # Stopband frequency AKA Ws; normalized to Nyquist frequency
passbandRip <- 2 # Passband Ripple (dB) AKA Rp
stopbandAtt <- 40 # Stopband Attenuation (dB) AKA Rs

# Determing the order and cut-off frequency based on the filter specifications
buttOrderCut <- buttord(passbandFreqN, stopbandFreqN, passbandRip, stopbandAtt)

# Creating low pass Butterworth filter of order n
buttFiltLP <- butter(buttOrderCut)

# Normalized cut-off frequency
cutFreqN <- buttOrderCut$Wc/freqN

# Filter sono data
videodata1 <- sono_WW
videodata1 <- c(rep(videodata1[1],500),videodata1,rep(videodata1[length(videodata1)],500))

# Using forward and reverse filtering to prevent phase shifts (i.e., making this a zero phase filter)
filterVPad <- filtfilt(buttFiltLP$b, buttFiltLP$a, videodata1)
endplot <- length(filterVPad)-500
trimdata <- filterVPad[501:endplot]
sono_VV <- trimdata

# Plot raw sono versus filtered
plot(sono_WW)
points(sono_TT, col = "blue")
points(sono_SS, col = "red")
points(sono_VV, type = "l", lwd = 2)

#Filter motor position
# Assigning the filter specification variables
freq <- 1000 # Frequency of data
freqN <- freq/2 # Frequency normalized to Nyquist frequency
passbandFreqN <- 200/freqN  # Passband frequency AKA Wp; normalized to Nyquist frequency
stopbandFreqN <- 500/freqN  # Stopband frequency AKA Ws; normalized to Nyquist frequency
passbandRip <- 2 # Passband Ripple (dB) AKA Rp
stopbandAtt <- 40 # Stopband Attenuation (dB) AKA Rs

# Determing the order and cut-off frequency based on the filter specifications
buttOrderCut <- buttord(passbandFreqN, stopbandFreqN, passbandRip, stopbandAtt)

# Creating low pass Butterworth filter of order n
buttFiltLP <- butter(buttOrderCut)

# Normalized cut-off frequency
cutFreqN <- buttOrderCut$Wc/freqN

videodata1 <- Lout
videodata1 <- c(rep(videodata1[1],500),videodata1,rep(videodata1[length(videodata1)],500))

# Using forward and reverse filtering to prevent phase shifts (i.e., making this a zero phase filter)
filterVPad <- filtfilt(buttFiltLP$b, buttFiltLP$a, videodata1)
endplot <- length(filterVPad)-500
trimdata <- filterVPad[501:endplot]
plot(Lout)
points(trimdata,col="red",type="l")
Lout <- trimdata

# Calculate change in MTU length
MTU_Lin <- MTU_length + Lout;
MTU_Lin <- MTU_Lin[1:length(MTU_Lin)-1]

# Calculate change in fascicle length assuming sono is 66% fascicle length
sono_L <- mean(sono_SS[1:200]);
fascicle_L <- 1.5*sono_L
fascicle <- fascicle_L * sono_VV/sono_L;

## Calculate fascicle velocity using quintic spline
smoothingcoef = 10^8
Time <- seq(from=1,to=length(fascicle),by=1)
PointModel <- smooth.Pspline(Time,fascicle,norder=5,method=1,spar=smoothingcoef)
fascicle_vel <- signif(predict(PointModel, Time, nderiv=1), digits=3)

# Calculate fascicle power
Force <- Force[1:4000]
fascicle_power <- fascicle_vel*Force

# Shortening starts when Force reaches maximum
shortstart <- which(Force[500:1000] >= mean(Force[1000:500]))[1]
shortstart <- shortstart+500

#Find the time at which we see max force
timemaxforce <- which(Force == max(Force));
timemaxforce <- timemaxforce[length(timemaxforce)]

#Find the time at which force drops below zero after that
timezero <- which(Force[timemaxforce:length(Force)] < 0)[1];
timezero <- timezero + timemaxforce;
shortend <- timezero;

# Calculate change in tendon length
tendon <- MTU_Lin- fascicle;
tendon <- tendon - tendon[1];

shortening <- MTU_Lin - MTU_Lin[1];
contraction <- fascicle - fascicle[1]

# The end of shortening is define as the point where the maximum shortening 
# distance (minimum length of MTU_Lin) is reached
shortend <- which(shortening==min(shortening))[1]
plot(contraction[1:shortend],Force[1:shortend],type="l")
points(shortening[shortstart:shortend],Force[shortstart:shortend],type="l")
points(-tendon[1:shortstart],Force[1:shortstart],type="l")

# Calculate variables of interest
averageforce <- mean(Force[shortstart:shortend])
tendonlength <- tendon[shortstart]
maxshortening <- shortening[shortend]
maxcontraction <- contraction[shortend]
timetoshort <- shortstart-500
loadenergy <- trapz(tendon[1:shortstart],Force[1:shortstart])

timezero <- 1100
peakfasciclev <- -min(fascicle_vel[shortstart:timezero]/fascicle_L*1000)
peakfasciclep <- min(fascicle_power[shortstart:timezero])

quartz(width=12,height=6)
par(mfrow=c(6,1),mar=c(2,5,2,2))
plot(Force[400:timezero])
plot(contraction[400:timezero])
plot(tendon[400:timezero])
plot(shortening[400:timezero])
plot(fascicle_vel[shortstart:timezero]/fascicle_L)
plot(fascicle_power[shortstart:timezero])

c(averageforce,tendonlength,maxshortening,maxcontraction,timetoshort,loadenergy,peakfasciclev,peakfasciclep)
