#### libraries ----
library(lubridate)
library(ggplot2)
library(dplyr)

#### data directory ----
#sapflow and sensor data parent directory
#source("/Users/ebenson/Desktop")


dirData <- "/Users/ebenson/Desktop/"

sapRaw <- read.csv(paste0(dirData,"/Sapflow_TableDT1.dat"),
                   header=FALSE,skip=4,na.strings=c("NAN"))

datSap <- sapRaw[,1:18]


#parse date
datSap$dateF <- ymd_hms(datSap$V1)
datSap$year <- year(datSap$dateF)
datSap$doy <- yday(datSap$dateF)
datSap$hour <- hour(datSap$dateF)+(minute(datSap$dateF)/60)
datSap$DD <- datSap$doy + (datSap$hour/24)

colnames(datSap ) <- c("date","record",paste0("dT",seq(1,16)), "dateF", "year", "doy", "hour", "DD")

datSap <- datSap[datSap$doy >= 160 & datSap$year == 2022,]

#### read in data ----

weather <- read.csv("/Users/ebenson/Desktop/weatherseptember22.csv")

basswoodmeas <- read.csv("/Users/ebenson/Desktop/basswoodmeas.csv")
hemlockmeas <- read.csv("/Users/ebenson/Desktop/hemlockmeas.csv")


#### initial plots ----

sensors <- read.csv("/Users/ebenson/Desktop/sensordata_061522.csv")

#add sapwood to sensors
#hemlock allometry from "Water use by eastern hemlock (Tsuga canadensis) 
#and black birch (Betula lenta): implications of effects of the hemlock woolly adelgid" by Daley et Al
#basswood allometry by Ewers et. Al: "Tree species effects on stand transpiration"

sapdepth <- lm(basswoodmeas$Sapwooddepth.cm ~ basswoodmeas$DBH.cm)

summary(sapdepth)

sapdepthHem <- lm(hemlockmeas$Sapwooddepth.cm ~ hemlockmeas$DBH.cm)

summary(sapdepthHem)

sensors$sd.cm <- ifelse(sensors$Tree.Type == "Hemlock", #if sensors is hemlock,
                        -0.0133 + (0.1252*sensors$DBH..cm.),
                        -0.7783 + (0.24546*sensors$DBH..cm.))



#organize data for easier calculations
tabledt <- datSap


dtAll <- data.frame(date= rep(tabledt$date, times = 16),
                    doy = rep(tabledt$doy, times = 16),
                    hour = rep(tabledt$hour, times = 16),
                    DD = rep(tabledt$DD, times = 16),
                    sensor = rep(seq(1,16), each = nrow(tabledt)),
                    dT = c(tabledt[,3],
                           tabledt[,4],
                           tabledt[,5],
                           tabledt[,6],
                           tabledt[,7],
                           tabledt[,8],
                           tabledt[,9],
                           tabledt[,10],
                           tabledt[,11],
                           tabledt[,12],
                           tabledt[,13],
                           tabledt[,14],
                           tabledt[,15],
                           tabledt[,16],
                           tabledt[,17],
                           tabledt[,18]))



#################
#check for dt outliers
#quantile(dtAll$dT, prob=seq(0,1,by=0.001))
#definitely few outliers. 99.5% and above are unusually high
#dtAll <- dtAll[dtAll$dT <= quantile(dtAll$dT, prob=0.995),]

#join sensor info into table dt
#make a doy that contains the same night
#so new day actually starts at 5 am not midnight
dtAll$doy5 <- ifelse(dtAll$hour < 5, dtAll$doy-1,dtAll$doy)

night <- dtAll[dtAll$hour < 5|dtAll$hour >= 22,]

#filter night so maximum in day and sensor is provided
maxnight1 <- night %>%
  group_by(sensor, doy5) %>%
  filter(dT == max(dT),na.rm=TRUE)
#remove duplicate maximums that occur for longer than 15 min
#just take earliest measurement
maxnight <- maxnight1  %>%
  group_by(sensor, doy5) %>%
  filter(hour == min(hour),na.rm=TRUE)

#ggplot(maxnight, aes(doy5,dT, color=sensor))+
#geom_point()
#isolate max and join back into table
maxJoin <- data.frame(sensor=maxnight$sensor,
                      doy5=maxnight$doy5,
                      maxDT = maxnight$dT)

ggplot(maxJoin, aes(doy5, maxDT, color = as.factor(sensor)))+
       geom_line()+
    geom_point()


#join backinto tabledt
dtCalct1 <- left_join(dtAll, maxJoin, by=c("sensor","doy5"))
#join sensor info
dtCalc <- left_join(dtCalct1 , sensors, by=c("sensor"="Sensor.Number"))

#from clearwater

#sap velocity m s-1 (v)
#v = 0.000119*k^1.231
#flow is F (L s-1) = v* A (m2, sapwood area)

#K= (dTmax - dT)/dT if sensor is fully within sapwood

#otherwise correction is:
#dt sap = (dT - b* Dtmax)/a

#a = proportion of probe in sapwood and b=1-a

dtCalc$a <- ifelse(dtCalc$sd.cm >= 3,1,
                   dtCalc$sd.cm/3)

dtCalc$b <- 1 - dtCalc$a

dtCalc$dTCor <- (dtCalc$dT - (dtCalc$b * dtCalc$maxDT))/dtCalc$a
dtCalc$K <- (dtCalc$maxDT - dtCalc$dTCor)/dtCalc$dTCor
dtCalc$velo <- 0.000119*(dtCalc$K^1.231)


#separate types
hemlockT <- dtCalc[dtCalc$Tree.Type == "Hemlock",]
basswoodT <- dtCalc[dtCalc$Tree.Type == "Basswood",]

hemlockQ <- quantile(hemlockT$velo,probs = seq(0,1,by=0.01),na.rm=TRUE)

hemlock <- hemlockT[hemlockT$velo<= hemlockQ[100],]

basswoodQ <- quantile(basswoodT$velo,probs = seq(0,1,by=0.01),na.rm=TRUE)


basswood <-basswoodT[basswoodT$velo<= basswoodQ[100],]



#############
#compare N & S sensors for hemlock
sens5 <- na.omit(data.frame(date = hemlock$date[hemlock$sensor == 5],
                    veloN = hemlock$velo[hemlock$sensor == 5]))

sens6 <- na.omit(data.frame(date = hemlock$date[hemlock$sensor == 6],
                    veloS = hemlock$velo[hemlock$sensor == 6]))

treeD1 <- inner_join(sens5,sens6, by="date")

plot(treeD1$veloN, treeD1$veloS)
#compare N & S sensors for hemlock
sens12 <- na.omit(data.frame(date = hemlock$date[hemlock$sensor == 12],
                     veloN = hemlock$velo[hemlock$sensor == 12]))

sens11 <- na.omit(data.frame(date = hemlock$date[hemlock$sensor == 11],
                     veloS = hemlock$velo[hemlock$sensor == 11]))

treeD2 <- inner_join(sens12,sens11, by="date")

sens15 <- na.omit(data.frame(date = hemlock$date[hemlock$sensor == 15],
                     veloN = hemlock$velo[hemlock$sensor == 15]))

sens14 <- na.omit(data.frame(date = hemlock$date[hemlock$sensor == 14],
                     veloS = hemlock$velo[hemlock$sensor == 14]))

treeD3 <- inner_join(sens15,sens14, by="date")

treeDirHem <- na.omit(rbind(treeD1,treeD2,treeD3))
#check relationship
azim.rel <- lm(treeDirHem$veloS ~ treeDirHem$veloN)
summary(azim.rel)

ggplot(treeDirHem, aes(veloN,veloS))+
  geom_point()+
  geom_abline()

# check basswood


sens1 <- na.omit(data.frame(date = basswood$date[basswood$sensor == 1],
                    veloN = basswood$velo[basswood$sensor == 1]))

sens2 <- na.omit(data.frame(date = basswood$date[basswood$sensor == 2],
                    veloS = basswood$velo[basswood$sensor == 2]))

treeB1 <- inner_join(sens1,sens2, by="date")

sens3 <- na.omit(data.frame(date = basswood$date[basswood$sensor == 3],
                    veloN = basswood$velo[basswood$sensor == 3]))

sens4 <- na.omit(data.frame(date = basswood$date[basswood$sensor == 4],
                    veloS = basswood$velo[basswood$sensor == 4]))

treeB2 <- inner_join(sens3,sens4, by="date")

sens8 <- na.omit(data.frame(date = basswood$date[basswood$sensor == 8],
                    veloN = basswood$velo[basswood$sensor == 8]))

sens9 <- na.omit(data.frame(date = basswood$date[basswood$sensor == 9],
                    veloS = basswood$velo[basswood$sensor == 9]))

treeB3 <- inner_join(sens8,sens9, by="date")

treeBDir <- rbind(treeB1,treeB2, treeB3)

azimB.rel <- lm(treeBDir$veloS ~ treeBDir$veloN)
summary(azimB.rel)

ggplot(treeBDir, aes(veloN,veloS))+
  geom_point()+
  geom_abline()



hemlock.tree <- hemlock %>%
  filter(Direction == "N")

basswood.tree <- basswood %>%
  filter(Direction == "N")

# running radial correction
hemlock.cor <- coefficients(azim.rel)
hemlock.tree$velo.cor <- (hemlock.tree$velo*0.5)+(((hemlock.tree$velo*hemlock.cor[2])+hemlock.cor[1])*0.5)

basswood.cor <- coefficients(azimB.rel)
basswood.tree$velo.cor <- (basswood.tree$velo*0.5)+(((basswood.tree$velo*basswood.cor[2])+basswood.cor[1])*0.5)



##############################
#### canopy leaf allometry   ----


#Ash allometry from literature

hemlockmeas$treeArea <- ((hemlockmeas$DBH.cm /2)^2)*pi

#plot(greenwood$dbh.cm,greenwood$sap.area.cm)

#saparea.reg <- lm(greenwood$sap.area.cm ~ greenwood$dbh.cm)
#summary(saparea.reg)


#sap cm2 = -9.6 + 8.854*DBH cm


##############################
#### Canopy calculations   ----

## sapwood area

# hemlock tree
plot(log(hemlockmeas$DBH.cm), log(hemlockmeas$SapwoodArea.cm2))
sapareaHem <- lm(log(hemlockmeas$SapwoodArea.cm2) ~ log(hemlockmeas$DBH.cm))
abline(sapareaHem)
summary(sapareaHem)

hemlock.tree$sap.areacm2 <- exp(-1.192 + 2.010*log(hemlock.tree$DBH..cm.))
#convert sap area to m2
hemlock.tree$sap.aream2 <- 0.0001*hemlock.tree$sap.areacm2

# basswood

#calculate heartwood

basswood.tree$bark <- (basswood.tree$DBH..cm.*0.0326) - 0.1708

basswood.tree$Htwd <- basswood.tree$DBH..cm.  - (basswood.tree$sd.cm*2) - (basswood.tree$bark*2)



#calculate sapwood area

basswood.tree$sap.areacm2 <- (pi*(((basswood.tree$sd.cm/2)+(basswood.tree$Htwd/2))^2))-(pi*((basswood.tree$Htwd/2)^2))
basswood.tree$sap.aream2 <-  0.0001*basswood.tree$sap.areacm2


#check relationship
# lm.log<- lm(log(basswoodLA$LA.m2) ~ log(basswoodLA$DBH.cm))
# summary(lm.log)
# plot(basswoodLA$DBH.cm, basswoodLA$LA.m2)
# plot(log(basswoodLA$DBH.cm), log(basswoodLA$LA.m2))
#regression log(LA (m2)) = -1.058 + 1.828 * log(dbh.cm)

#estimate leaf area in m2

#crown = 1.6961 + 0.4233(DBH)

#basswood leaf area to the best of our ability
#mean basswood weight = 22.1324289 g/m2
#1/slw = 0.04518257

basswoodlm <- read.csv("C:/Users/ebenson/Desktop/DettmannMcfarlane.csv")
blm <- basswoodlm[basswoodlm$SPECIES =="Tilia americana",]
plot(log(blm$DBH_CM), log(blm$LEAF_DRY_MASS_KG))

b.mod <- lm(log(blm$LEAF_DRY_MASS_KG)~log(blm$DBH_CM))
summary(b.mod)

basswood.tree$biomass.kg = exp(-4.25 + 1.79*(log(basswood.tree$DBH..cm.)))

basswood.tree$biomass.g <- basswood.tree$biomass.kg*1000

exp(2.014+(0.94*log(basswood.tree$DBH..cm.)) + (-0.632*(log(basswood.tree$DBH..cm.)^2)))
#*total grams of leaves

#conversion from Jurik 1986
basswood.tree$LA.m2 <- 0.04518257*basswood.tree$biomass.g

#leaf area in m2 for hemlock, from Kenefic et al
hemlock.tree$LA.m2 <- 7.5432 + (0.3659*(hemlock.tree$sap.areacm2))

##############################
#### Flow calculations   ----

#flow rate according to clearwater
#F(L s-1) =  v(m s-1)* A (m2)

hemlock.tree$Flow.m3.s <- hemlock.tree$velo * hemlock.tree$sap.aream2

basswood.tree$Flow.m3.s <- basswood.tree$velo * basswood.tree$sap.aream2

#convert to L per secton

hemlock.tree$Flow.L.s <- hemlock.tree$Flow.m3.s * 1000

basswood.tree$Flow.L.s <- basswood.tree$Flow.m3.s * 1000

#normalize by canopy leaf area
hemlock.tree$Flow.L.m2.s <- hemlock.tree$Flow.L.s /hemlock.tree$LA.m2

basswood.tree$Flow.L.m2.s <- basswood.tree$Flow.L.s /basswood.tree$LA.m2

#summarize total per day for each tree
#remove NA
hemlock.treeNN <- hemlock.tree[is.na(hemlock.tree$Flow.L.s)==FALSE,]
#calculate total water use by each tree in a day
#total liters used in 15 min period
hemlock.treeNN$L.p <- hemlock.treeNN$Flow.L.s* 60 *15
#per canopy area
hemlock.treeNN$L.p.m2  <- hemlock.treeNN$L.p/hemlock.treeNN$LA.m2

#summarize total per day for each tree
#remove NA
basswood.treeNN <- basswood.tree[is.na(basswood.tree$Flow.L.s)==FALSE,]
#calculate total water use by each tree in a day
#total liters used in 15 min period
basswood.treeNN$L.p <- basswood.treeNN$Flow.L.s* 60 *15
# per canopy area
basswood.treeNN$L.p.m2  <- basswood.treeNN$L.p/basswood.treeNN$LA.m2


##############################
#### Summary tables    ----

hemlock.treeNN$hour1 <- floor(hemlock.treeNN$hour)
#summary table
#flow L s every 15 min by treatment
hemlock.Flow <- hemlock.treeNN %>%
  group_by(doy, hour, Tree.Type) %>%
  summarise(mean15 = mean(Flow.L.s),sd=sd(Flow.L.s), n=length(Flow.L.s))

#flow L s every hour
hemlock.Flow.hr <- hemlock.treeNN %>%
  group_by(doy, hour1, sensor, Tree.Type) %>%
  summarise(sd=sd(Flow.L.s), n=length(Flow.L.s), meanHr=mean(Flow.L.s))

#flow L m-2 leaf s-1 by 15min
hemlock.Flow.m2 <- hemlock.treeNN %>%
  group_by(doy, hour, Tree.Type) %>%
  summarise(mean15 = mean(Flow.L.m2.s),
            sd=sd(Flow.L.m2.s),
            n=length(Flow.L.m2.s)) 

#flow L m-2 every hour
hemlock.Flow.m2.hr <- hemlock.treeNN %>%
  group_by(doy, hour1, sensor, Tree.Type) %>%
  summarise(sd=sd(Flow.L.m2.s),
            n=length(Flow.L.m2.s),
            meanHr = mean(Flow.L.m2.s)) 

dayAllhem <- data.frame(doy= rep(rep(seq(166, 251), each=24), times=length(unique(hemlock.treeNN$sensor))), 
                     hour1= rep(rep(seq(0, 23), times=length(seq(166, 251))), times=length(unique(hemlock.treeNN$sensor))), 
                     sensor= rep(unique(hemlock.treeNN$sensor), each= length(rep(seq(0, 23), times=length(seq(166, 251))))))

dayAllhem1 <- data.frame(doy= rep(seq(166, 251), each=24),
                          hour1= rep(seq(0, 23), times=length(seq(166, 251))))

dayAllhemjoin <- full_join(dayAllhem, hemlock.Flow.hr, by =c("doy", "hour1", "sensor"))


dayAllhemjoinNA <- na.omit(dayAllhemjoin)

hemlock.avg.hr <- dayAllhemjoinNA %>%
  group_by(doy, hour1, sensor)%>%
  summarise(Avg=mean(meanHr), count=n())

hemlock.avg.hr1 <- hemlock.avg.hr %>%
  group_by(doy, hour1)%>%
  summarise(mean=mean(Avg), count=n())

#only use time points with at least 2 trees

hemlock.avg.hr1 <- hemlock.avg.hr1[ hemlock.avg.hr1$count >2,]

hemlock.day.count <- hemlock.avg.hr1 %>%
  group_by(doy) %>%
  summarise(count= n())

hemlock.day.all <- hemlock.day.count %>%
  filter(count>=20)

dayallhem.hr1 <- left_join(dayAllhem1, hemlock.avg.hr1, by= c("doy", "hour1"))

dayallhem.hr <- inner_join(dayallhem.hr1, hemlock.day.all, by= ("doy"))

# gapfill hemlock

install.packages("zoo")
library(zoo)

dayallhem.hr$dateT <- paste0(as.Date(dayallhem.hr$doy-1, origin="2022-01-01"), " ", dayallhem.hr$hour1, ":00")
dayallhem.hr$dateP <- ymd_hm(dayallhem.hr$dateT)

hemzoo <- zoo(dayallhem.hr$mean, dayallhem.hr$dateP)
hemgap <- na.approx(hemzoo, na.rm = FALSE)
hemgapD <- data.frame(hemgap)
dayallhem.hr$gap.fill <- hemgap

dayallhem.hr$flow.L.hr <- dayallhem.hr$gap.fill * 3600

hemflow.day <- dayallhem.hr %>%
  group_by(doy)%>%
  summarise(flow.mean = mean(flow.L.hr), flow.day = sum(flow.L.hr))

hemflow.day$type <- rep("Hemlock", nrow(hemflow.day))

#hemlock.Flow <- hemlock.Flow[ hemlock.Flow$n >=2,]
#hemlock.Flow.hr <- hemlock.Flow.hr[ hemlock.Flow.hr$n >=2,]

#hemlock.Flow.m2 <- hemlock.Flow.m2[ hemlock.Flow.m2$n >=2,]
#hemlock.Flow.m2.hr <- hemlock.Flow.m2.hr[ hemlock.Flow.m2.hr$n >=2,]

#hemlock.Flow$se <- hemlock.Flow$sd/sqrt(hemlock.Flow$n)
#hemlock.Flow.m2$se <- hemlock.Flow.m2$sd/sqrt(hemlock.Flow.m2$n)

#hemlock.Flow.hr$se <- hemlock.Flow.hr$sd/sqrt(hemlock.Flow.hr$n)
#hemlock.Flow.m2.hr$se <- hemlock.Flow.m2.hr$sd/sqrt(hemlock.Flow.m2.hr$n)


#BASSWOOD
basswood.treeNN$hour1 <- floor(basswood.treeNN$hour)

#flow L s per 15 mins
basswood.Flow <- basswood.treeNN %>%
  group_by(doy, hour, Tree.Type) %>%
  summarise(mean15 = mean(Flow.L.s),sd=sd(Flow.L.s), n=length(Flow.L.s))

# per hour
basswood.Flow.hr <- basswood.treeNN %>%
  group_by(doy, hour1, sensor, Tree.Type) %>%
  summarise(meanHr = mean(Flow.L.s),sd=sd(Flow.L.s), n=length(Flow.L.s))

#flow L m-2 leaf s-1 by 15min
basswood.Flow.m2 <- basswood.treeNN %>%
  group_by(doy, hour, Tree.Type) %>%
  summarise(mean15 = mean(Flow.L.m2.s),sd=sd(Flow.L.m2.s), n=length(Flow.L.m2.s))

#hour
basswood.Flow.m2.hr <- basswood.treeNN %>%
  group_by(doy, hour1, sensor, Tree.Type) %>%
  summarise(meanHr = mean(Flow.L.m2.s),sd=sd(Flow.L.m2.s), count = n())

dayAllbass <- data.frame(doy= rep(rep(seq(166, 251), each=24), times=length(unique(basswood.treeNN$sensor))), 
                     hour1= rep(rep(seq(0, 23), times=length(seq(166, 251))), times=length(unique(basswood.treeNN$sensor))), 
                     sensor= rep(unique(basswood.treeNN$sensor), each= length(rep(seq(0, 23), times=length(seq(166, 251))))))

dayAllbass1 <- data.frame(doy= rep(seq(166, 251), each=24),
                                   hour1= rep(seq(0, 23), times=length(seq(166, 251))))
  
  
dayAllbassjoin <- full_join(dayAllbass, basswood.Flow.hr, by =c("doy", "hour1", "sensor"))

sum(is.na(dayAllbassjoin$meanHr))

dayAllbassjoinNA <- na.omit(dayAllbassjoin)

basswood.avg.hr <- dayAllbassjoinNA %>%
  group_by(doy, hour1, sensor)%>%
  summarise(Avg=mean(meanHr), count=n())

basswood.avg.hr1 <- basswood.avg.hr %>%
  group_by(doy, hour1)%>%
  summarise(mean=mean(Avg), count=n(), sd= sd(Avg))

#only use time points with at least 2 trees

basswood.avg.hr1 <- basswood.avg.hr1[ basswood.avg.hr1$count >2,]

basswood.day.count <- basswood.avg.hr1 %>%
  group_by(doy) %>%
  summarise(count= n())

basswood.day.all <- basswood.day.count %>%
  filter(count>=20)

dayallbass.hr1 <- left_join(dayAllbass1, basswood.avg.hr1, by= c("doy", "hour1"))

dayallbass.hr <- inner_join(dayallbass.hr1, basswood.day.all, by= ("doy"))

install.packages("zoo")
library(zoo)

dayallbass.hr$dateT <- paste0(as.Date(dayallbass.hr$doy-1, origin="2022-01-01"), " ", dayallbass.hr$hour1, ":00")
dayallbass.hr$dateP <- ymd_hm(dayallbass.hr$dateT)

basszoo <- zoo(dayallbass.hr$mean, dayallbass.hr$dateP)
bassgap <- na.approx(basszoo, na.rm = FALSE)
bassgapD <- data.frame(bassgap)
dayallbass.hr$gap.fill <- bassgap

dayallbass.hr$flow.L.hr <- dayallbass.hr$gap.fill* 3600

bassflow.day <- dayallbass.hr %>%
  group_by(doy)%>%
  summarise(flow.mean = mean(flow.L.hr), flow.day = sum(flow.L.hr))

bassflow.day$type <- rep("Basswood", nrow(bassflow.day))

basshemflow.day <- rbind(bassflow.day, hemflow.day)


############# graphs


# average daily sap flow

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/sapflow_dailyavg.png",
    height=10, width=20, units="in",res=300) 
ggplot(basshemflow.day , aes(doy, flow.mean, color = type))+
  geom_point(size=10)+
  geom_line(size=3, lty=1, alpha=0.5)+
  ggtitle("Average Daily Sap Flow in Basswood and Hemlock")+
  scale_color_manual(values=c("mediumorchid4", "slategrey"))+
  xlab("Day of Year")+
  ylab(expression(paste("Sap flow (liters per day)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()


# daily total sap flow all days

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/sapflow_dailytotal.png",
    height=10, width=20, units="in",res=300) 
ggplot(basshemflow.day , aes(doy, flow.day, color = type))+
  geom_point(size=10)+
  geom_line(size=3, lty=1, alpha=0.5)
  ggtitle("Daily Total Sap Flow in Basswood and Hemlock")+
  scale_color_manual(values=c("mediumorchid4", "slategrey"))+
  xlab("Day of Year")+
  ylab(expression(paste("Sap flow (liters per day)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()

#divide up days into dry and wet periods

sapfluxdry223 <- basshemflow.day %>%
  filter(doy >= 223 & doy <= 233)

sapfluxwet234 <- basshemflow.day %>%
  filter(doy >= 234 & doy <= 244)


#DRY PERIOD

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/sapflow_223.png",
    height=7, width=10, units="in",res=300) 
ggplot(sapfluxdry223 , aes(doy, flow.day, color = type))+
  geom_point(size=10)+
  geom_line(size=3, lty=1, alpha=0.5)+
  ggtitle("Daily Total Sap Flow Dry Period")+
  scale_color_manual(values=c("mediumorchid4", "slategrey"))+
  xlab("Day of Year")+
  ylab(expression(paste("Sap flow (liters per day)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()

#WET PERIOD

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/sapflow_234.png",
    height=7, width=10, units="in",res=300) 
ggplot(sapfluxwet234 , aes(doy, flow.day, color = type))+
  geom_point(size=10)+
  geom_line(size=3, lty=1, alpha=0.5)+
  ggtitle("Daily Total Sap Flow Wet Period")+
  scale_color_manual(values=c("mediumorchid4", "slategrey"))+
  xlab("Day of Year")+
  ylab(expression(paste("Sap flow (liters per day)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()

