library(lubridate)
library(dplyr)

### tomst 1025-------

dirT <- "K:/Environmental_Studies/hkropp/Ecohydro/10_21_2022"

soilclass <- read.csv("C:/Users/ebenson/Desktop/soilclass.csv")

colnames(soilclass) <- c("location", "texture")
tomstF <- list.files(paste0(dirT))
fileSN <- character()
for(i in 1:length(tomstF)){
  fileSN[i] <- as.numeric(strsplit(tomstF, "_")[[i]][2])
}

# read in files
datT <- list()

enddate <- list()

for(i in 1: length(fileSN)){
  datT[[i]] <- read.csv(paste0(dirT,"/",tomstF[i]), sep=";",
                        header=FALSE)[,1:9]
  datT[[i]]$SN <- rep(fileSN[i], nrow(datT[[i]]))
  enddate[[i]] <- datT[[i]][nrow(datT[[i]]),]
}

endDF <- do.call("rbind", enddate)
tomst <- do.call("rbind", datT)
colnames(tomst)[1:9] <- c("record","date","tz","Tm6","T2","T15","SM","shake","errFlag")

tomst$Tm6 <- as.numeric(gsub("\\,","\\.", tomst$Tm6))
tomst$T2 <- as.numeric(gsub("\\,","\\.", tomst$T2))
tomst$T15 <- as.numeric(gsub("\\,","\\.", tomst$T15))

tomst$dateF <- ymd_hm(tomst$date) 
tomst$estD <- with_tz(tomst$dateF, tzone="America/New_York")

TOMSTSensor <- data.frame(SN= c(
  94236481,# deciduous
  94236491,# deciduous
  94236490,# deciduous
  94236492, # ag
  94236484, # ag
  94236487, # ag
  94236489, # lawn
  94236495, # lawn
  94236485, # lawn
  94236483, # sapflux
  94214744, # sapflux
  94214743, # sapflux
  94236482, # reforestation
  94236486, # reforestation
  94236493 # reforestation
),
location=c("deciduous",
           "deciduous",
           "deciduous",
           "ag field",
           "ag field",
           "ag field",
           "mowed lawn",
           "mowed lawn",
           "mowed lawn",
           "sapflux",
           "sapflux",
           "sapflux",
           "reforestation",
           "reforestation",
           "reforestation"
))


tomst$SN <- as.numeric(tomst$SN)
tomst <- left_join(tomst, TOMSTSensor, by="SN")

### tomst 1129-----

dirT2 <- "K:/Environmental_Studies/hkropp/Ecohydro/11_29_2022"

tomstF2 <- list.files(paste0(dirT2))
fileSN2 <- character()
for(i in 1:length(tomstF2)){
  fileSN2[i] <- as.numeric(strsplit(tomstF2, "_")[[i]][2])
}

namesSN2 <- data.frame(SN=fileSN2,
                       fileID = seq(1,length(fileSN2)))
end2read <- left_join(namesSN2, endDF, by="SN")

replace94 <- which(end2read$SN == 94236494)
replace84 <- which(end2read$SN ==94236484)
replace91 <- which(end2read$SN == 94207591)


end2read$V2[replace94] <- "2022.10.22 20:15"
end2read$V2[replace91] <- "2022.10.22 20:30"
end2read$V2[replace84] <- "2022.10.22 18:45"

# read in files
datT2 <- list()

for(i in 1: length(fileSN2)){
  datT2[[i]] <- read.csv(paste0(dirT2,"/",tomstF2[i]), sep=";",
                        header=FALSE)[,1:9]
  datT2[[i]]$SN <- rep(fileSN2[i], nrow(datT2[[i]]))
}

whichrow <- numeric()

for(i in 1: length(fileSN2)){
  
whichrow[i] <- which(datT2[[i]]$V2 == end2read$V2[end2read$fileID==i])
}

datT2crop <- list()
for(i in 1: length(fileSN2)){
  
  datT2crop[[i]] <- datT2[[i]][whichrow[i]+1:nrow(datT2[[i]]),]
}


tomst2 <- do.call("rbind", datT2crop)
colnames(tomst2)[1:9] <- c("record","date","tz","Tm6","T2","T15","SM","shake","errFlag")

tomst2$Tm6 <- as.numeric(gsub("\\,","\\.", tomst2$Tm6))
tomst2$T2 <- as.numeric(gsub("\\,","\\.", tomst2$T2))
tomst2$T15 <- as.numeric(gsub("\\,","\\.", tomst2$T15))

tomst2$dateF <- ymd_hm(tomst2$date) 
tomst2$estD <- with_tz(tomst2$dateF, tzone="America/New_York")

TOMSTSensor2 <- data.frame(SN= c(
  94236481,# deciduous
  94236491,# deciduous
  94236490,# deciduous
  94236492, # ag
  94236484, # ag
  94236487, # ag
  94236489, # lawn
  94236495, # lawn
  94236485, # lawn
  94236483, # sapflux
  94214744, # sapflux
  94236494, # sapflux
  94236482, # reforestation
  94207591, # reforestation
  94236493 # reforestation
),
location=c("deciduous",
           "deciduous",
           "deciduous",
           "ag field",
           "ag field",
           "ag field",
           "mowed lawn",
           "mowed lawn",
           "mowed lawn",
           "sapflux",
           "sapflux",
           "sapflux",
           "reforestation",
           "reforestation",
           "reforestation"
))


tomst2$SN <- as.numeric(tomst2$SN)
tomst2 <- left_join(tomst2, TOMSTSensor2, by="SN")
### combine all-----
tomst.bind <- rbind(tomst, tomst2)
tomst.bind2 <- na.omit(tomst.bind)

tomst.all <- left_join(tomst.bind2, soilclass, by="location")

tomst.all$doy <- yday(tomst.all$estD) 
tomst.all$hour <- hour(tomst.all$estD) 


tomst.all$SMcor <- ifelse(tomst.all$texture == "silty loam",
  (-0.00000001*(tomst.all$SM^2)) + (0.0003*tomst.all$SM) -0.1751,
  (-0.00000001*(tomst.all$SM^2)) + (0.0002*tomst.all$SM) -0.1499) #all others are sandy loam

### data----

hourlySM <- tomst.all %>%
  group_by(location, doy, hour) %>%
  summarize(soil.moist = mean(SMcor, na.rm=TRUE),
            n.count = length(na.omit(SMcor)),
            T15 = mean(T15, na.rm=TRUE))


doyCompleteDF <- data.frame(doy=rep(seq(160,333), each=24),
                            hour=rep(seq(0,23), times=length(seq(160,333))))

hourlySMfull <- full_join(hourlySM,doyCompleteDF, by=c("doy","hour"))

dailySM <- tomst.all %>%
  group_by(location, doy) %>%
  summarize(soil.moist = mean(SMcor, na.rm=TRUE))

# soil moisture is negative before it is installed remove

daySM <- dailySM %>%
  filter(soil.moist > 0.12)

hourSM <- hourlySM %>%
  filter(soil.moist > 0.12)
hourlySMfull$soil.moistC <- ifelse(hourlySMfull$soil.moist <= 0.12,NA,hourlySMfull$soil.moist )
hourlySMfinal <- hourlySMfull %>%
  filter(is.na(location) == FALSE)
library(ggplot2)
ggplot(daySM, aes(x=doy, y= soil.moist, color=location))+
  geom_point()+
  geom_line()

######----  Graphs

# all soil moisture days

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/soil_all1.png",
    height=10, width=20, units="in",res=300) 

ggplot(hourlySMfinal, aes(x=doy+(hour/24), y= soil.moistC, color=location))+
  geom_point()+
  geom_line(alpha=0.5)+
  geom_polygon(aes(x = c(290, 290, 315, 315), y=c(0, 0.7, 0.7, 0)))+
  ggtitle("Soil Moisture in Land Cover Classes")+
  scale_color_manual(values=c("#FFBF00", "#E83F6F", "#1C6087", "#4DA167", "#83469B"))+
  xlab("Day of Year")+
  ylab(expression(paste("Soil Moisture (cm"^"3"," Water cm"^"-3","Soil)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.size = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()

# subsetting wet and dry periods

subSMdry223 <- hourlySMfinal %>%
  filter(doy >= 223 & doy <= 233)
subSMwet234 <- hourlySMfinal %>%
  filter(doy >= 234 & doy <= 244)

# DRY DATA

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/soil_2231.png",
    height=10.5, width=15, units="in",res=300) 


ggplot(subSMdry223, aes(x=doy+(hour/24), y=soil.moistC, color=location))+
  geom_point()+
  geom_line(alpha=0.5)+
  ggtitle("Soil Moisture Dry Period")+
  scale_color_manual(values=c("#FFBF00", "#E83F6F", "#1C6087", "#4DA167", "#83469B"))+
  xlab("Day of Year")+
  ylim(0,0.45)+
  ylab(expression(paste("Soil Moisture (cm"^"3"," Water cm"^"-3","Soil)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.size = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()

# WET DATA

png("K:/Environmental_Studies/hkropp/Ecohydro/graphs/soil_2341.png",
    height=10.5, width=15, units="in",res=300) 


ggplot(subSMwet234, aes(x=doy+(hour/24), y=soil.moistC, color=location))+
  geom_point()+
  geom_line(alpha=0.5)+
  ggtitle("Soil Moisture Wet Period")+
  scale_color_manual(values=c("#FFBF00", "#E83F6F", "#1C6087", "#4DA167", "#83469B"))+
  xlab("Day of Year")+
  ylim(0,0.57)+
  ylab(expression(paste("Soil Moisture (cm"^"3"," Water cm"^"-3","Soil)")))+
  theme(axis.text.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title=element_text(size=25),
        title =element_text(size=30, face='bold'),
        legend.key.size = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=25)) #change legend text font size)
dev.off()
