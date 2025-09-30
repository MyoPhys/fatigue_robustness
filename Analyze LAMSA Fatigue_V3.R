require(pracma)
library(signal)
library(pspline)

# Edit the metadata for the trial being run
MTU_length <- 43 #MTU length (mm)
mass <- 0.0020041 #muscle mass (kg)
threshold <- 4.4156 #force threshold in model (N)

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
plot(shortening[1:timezero])
points(contraction[1:timezero])
points(-tendon[1:timezero])

## Calculate fascicle velocity using a quintic spline
smoothingcoef = 10^8
Time <- seq(from=1,to=length(fascicle),by=1)
PointModel <- smooth.Pspline(Time,fascicle,norder=5,method=1,spar=smoothingcoef)
fascicle_vel <- signif(predict(PointModel, Time, nderiv=1), digits=3)

## Calculate tendon velocity using a quintic spline
smoothingcoef = 10^8
Time <- seq(from=1,to=length(tendon),by=1)
PointModel <- smooth.Pspline(Time,tendon,norder=5,method=1,spar=smoothingcoef)
tendon_vel <- signif(predict(PointModel, Time, nderiv=1), digits=3)

## Calculate MTU velocity using a qunitic spline
smoothingcoef = 10^8
Time <- seq(from=1,to=length(MTU_Lin),by=1)
PointModel <- smooth.Pspline(Time,MTU_Lin,norder=5,method=1,spar=smoothingcoef)
MTU_vel <- signif(predict(PointModel, Time, nderiv=1), digits=3)

Force <- Force[1:4000]
MTU_power <- -MTU_vel*Force/mass
fascicle_power <- -fascicle_vel*Force/mass
tendon_power <- -tendon_vel*Force/mass

thresholdtime <- which(Force>threshold)[1]

#Separate the tendon into loading and unloading
tendon_load_power <- tendon_power
for(i in 1:length(tendon_load_power)){
  if(tendon_power[i] >= 0){
    tendon_load_power[i] <- 0
  }
}

tendon_unload_power <- tendon_power
for(i in 1:length(tendon_unload_power)){
  if(tendon_power[i] <= 0){
    tendon_unload_power[i] <- 0
  }
}

tendontime <- which(tendon[1:timezero]==max(tendon[1:timezero]))
forcetime <- which.max(Force)

# Plot results
quartz(width=12,height=6)
par(mfrow=c(2,3),mar=c(2,5,2,2))
plot(Force[500:timezero],type="l",lwd=2, ylab = "Force (N)",cex.lab=2,cex.axis=1.5)

plot(Lout[500:timezero]-Lout[1],type="l",lwd=2,ylim=c(min((fascicle - fascicle[1])),max(tendon)), ylab = "Shortening (mm)",cex.lab=2,cex.axis=1.5)
points(fascicle[500:timezero]-fascicle[1], type="l",col="red",lwd=2)
points(tendon[500:timezero], type="l",col="blue",lwd=2)
abline(h=0,lty=2)

plot(Time[500:timezero],MTU_vel[500:timezero],type="l",lwd=2,ylim = c(min(fascicle_vel),1.5*max(tendon_vel)), ylab = "Velocity (m/s)",cex.lab=2,cex.axis=1.5)
points(Time[500:timezero],fascicle_vel[500:timezero], type="l",col="red",lwd=2)
points(Time[500:timezero],tendon_vel[500:timezero], type="l",col="blue",lwd=2)
abline(h=0,lty=2)

plot(Time[500:timezero],MTU_power[500:timezero],type="l",lwd=2,ylim = c(min(tendon_power),1.5*max(MTU_power)), ylab = "Power (W/kg)",cex.lab=2,cex.axis=1.5)
points(Time[500:timezero],fascicle_power[500:timezero], type="l",col="red",lwd=2)
points(Time[500:timezero],tendon_power[500:timezero], type="l",col="blue",lwd=2)
abline(h=0,lty=2)


# Calculate variables of interest
# Energy loaded into tendon before threshold
load_pre <- trapz(tendon[1:thresholdtime],Force[1:thresholdtime])

# Energy loaded into tendon during overshoot
overshoot <- trapz(tendon[thresholdtime:tendontime],Force[thresholdtime:tendontime])

# Energy released from tendon
e_return <- trapz(tendon[tendontime:timezero],Force[tendontime:timezero])

# Direct Muscle Work
direct <- trapz(fascicle[tendontime:timezero],Force[tendontime:timezero])

plot(tendon[1:timezero],Force[1:timezero])

max_tendon <- max(tendon[1:1000])

max_force <- max(Force)

Force[tendontime]

peak_tendonp <- max(tendon_unload_power)
peak_fasciclep <- max(fascicle_power[tendontime:timezero])
loading_power <- max(fascicle_power[1:tendontime])
peak_MTUp <- max(MTU_power)
peak_fasciclev <- min(fascicle_vel[1:tendontime])*1000/fascicle_L
loading_time <- tendontime-500
recoiltime <- timezero - tendontime

c(max_tendon,max_force,load_pre,overshoot,e_return,direct, peak_tendonp,peak_fasciclep,peak_MTUp,peak_fasciclev,loading_power,loading_time, recoiltime)

sono_L