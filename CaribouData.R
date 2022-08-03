# .libPaths("C:/Users/erneilso/Documents/R/win-library/3.5")
# .libPaths("P:/R/R-3.6.1/3.5")
library(sp)
library(rgdal)
library(raster)
library(igraph)
library(pROC)
# library(MASS)
# library(ks)
library(tidyverse)
# library(rgeos)
# library(lme4)
# library(MuMIn)
# library(sf)
# library(ggspatial)
# library(gridExtra)
# library(rasterVis)
# library(maptools)
(amt)
# library(trajr)
library(ggplot2)
# library(effects)
# library(reshape2)

### changes coming though?


### stole this from https://www.rdocumentation.org/packages/spatialEco/versions/1.3-2/source
raster.vol <- function(x, p = 0.95, sample = FALSE, spct = 0.05) {
  if( sample == FALSE ) {
    den <- raster::getValues(x)
    z <- sort(den[!is.na(den)], decreasing=TRUE)
    y <- cumsum(as.numeric(z))
    i <- sum(y <= p * y[length(y)])
    return(raster::setValues(x, den >= z[i])) 
  } else {
    den <- raster::sampleRegular(x, round(raster::ncell(x) * spct, digits = 0), sp = TRUE) 
    den@data <- data.frame(idx = 1:nrow(den), den = den@data[,1], cls = 0)
    sum.p <- sum(den@data[,"den"], na.rm = TRUE) * p
    den@data <- den@data[order(-den@data$den),] 
    i=0; j=0
    while(i <= sum.p) { 
      j=j+1
      if( !is.na(den@data[,"den"][j])) {
        i = i + den@data[,"den"][j]
        den@data[,"cls"][j] <- 1
      } else {
        den@data[,"cls"][j] <- NA
      }	
    }						   
    den@data <- den@data[order(den@data[,"idx"]),]
    return( den[den$cls == 1 ,] )	
  }
}

localPath = 'F:/LocalProjects/Permafrost/LocalData'



###############
## The analysis starts here. We first get a shapefile of the seasonal home ranges for all the caribou data to make a buffer study area.

####################
 # caribou data ####
####################

C <- readOGR(localPath,"NWT_Caribou")
# C <- read_sf("P:/Projects/Permafrost/Analysis/Data/Vector","NWT_Caribou")
# C <- read_sf("C:/Users/erneilso/Documents/LocalProjects/Permafrost/LocalData","NWT_Caribou")

###
# caribou dates
C$DateTime<-as.POSIXct(strptime(as.character(C$timestamp), format="%Y-%m-%d %H:%M:%S",tz="GMT"))
C$Date = format.Date(C$DateTime, "%Y-%m-%d")
C$Time = format.Date(C$DateTime, "%H:%M:%S")
C$Month = format.Date(C$DateTime, "%m")
C$Year = format.Date(C$DateTime, "%Y")


length(C$Date[which(!is.na(C$Date))])

###
# # Add julian day
# for(i in 1:nrow(C)){
#   C$JDate[i] <- round(difftime(C$Date[i], as.POSIXlt(paste0(format.Date(C$Date[i], "%Y"),'-01-01'))))+1
#   print(i)}


###
# seasons # 2019-20 CIMP boreal caribou habitat and nutrition FIGURES_DRAFT
# Date ranges: winter low Dec 27-March 18;
# pre-calving peak March 19-May 15th; calving low May 16-June 29;
# mid-summer peak June 30-Aug 22; Late summer/rut low Aug 23-Oct 17;
# Late fall peak Oct 18-Dec 26.
# calving is actually mid-May with tails of bell curve extending from beginning to end of May. 
wintermonths = c("12","01","02","03") 
calfmonths = c("04","05","06") 
summermonths = c("07","08","09") 
fallmonths = c("10", "11") 

C$Season = ifelse(C$Month %in% wintermonths, "w", ifelse(C$Month %in% calfmonths, "c", ifelse(C$Month %in% summermonths, "s", "f")))
## only need to check if it is winter of the year before
C$SY = ifelse(as.numeric(format.Date(C$Date, "%m")) ==12,
              paste(as.numeric(format.Date(C$Date, "%Y")) + 1, C$Season, sep=""),
              paste(as.numeric(format.Date(C$Date, "%Y")),     C$Season, sep=""))

table(C$SY,C$Year)

###
# individuals
cfid = gsub("WL-", "20",C$individu_1)
C$ID = gsub("-", "_", cfid)
C$cilabel = paste(C$ID,C$SY, sep="_")
max(table(C$cilabel))
table(C$ID)


# ##########
# # resample by fix rate 
# C_sample <- C %>%
#   group_by(ID,Date) %>%
#   sample_n(1) 
#                          
# length(unique(C$ID))



# ######################################
# # caribou durations and fix rates ####
# ######################################
# # sort
# C <- C[order(C$ID,C$DateTime),]
# Cs <- unique(C$ID)
# car_colors = colorRampPalette(c("blue", "red"))( length(unique(C$ID)) )  
# 
# Cdattab = data.frame(C=NA,N=NA,MaxDate=as.POSIXct(NA),MinDate=as.POSIXct(NA),FRs=NA)
# ci<-C[which(C$ID==Cs[1]),]
# plot(ci,xlim=c(CBB[1],CBB[2]), ylim=c(CBB[3],CBB[4]))
# for (i in 1:length(Cs)){
#   ci<-C[which(C$ID==Cs[i]),]
#   plot(ci,xlim=c(CBB[1],CBB[2]), ylim=c(CBB[3],CBB[4]),add=T,col=sample(car_colors,1))
#   Cdattab[i,1]=as.character(Cs[i])
#   Cdattab$N[i]=nrow(ci)
#   Cdattab$MaxDate[i]=max(ci$Date)
#   Cdattab$MinDate[i]=min(ci$Date)
#   Cdattab$Duration[i]=difftime(Cdattab$MaxDate[i],Cdattab$MinDate[i],units="days")
#   Cdattab$FRs[i]=length(unique(table(ci$Date)))
# }
# plot(Cdattab$Duration,Cdattab$N)
# hist(Cdattab$N/Cdattab$Duration)
# 
# 
# # this takes too long and the FRs are all over the place!!!!  
# for (i in 2:nrow(C)){
#   print(i)
#   if(C$ID[i] == C$ID[i-1] & !is.na(C$ID[i])){
#     C$FR[i] = difftime(C$DateTime[i],C$DateTime[i-1],units="hours")
#   }else{C$FR[i] = NA}
# }
# 
# data.frame(table(round(C$FR)))
# 
# 
# ######################################



################

## body condition data
bcd = read.csv(paste(localPath, "/NonSpatialData/BCDat.csv", sep=""))
bcd$Date<-as.POSIXct(strptime(as.character(bcd$Cap.date), format="%Y-%m-%d"))
bcd = bcd%>%
  mutate(BM = as.numeric(as.character(Body.mass)),
         BF = as.numeric(as.character(Body.fat)),
  )%>%
  select(-Body.mass, -Body.fat)
bcd = bcd[order(bcd$Date),]


bcd %>%
  filter(Date<as.POSIXct(strptime(as.character("2019"), format="%Y")))%>%
  ggplot(aes(Date,BM,group=tid))+
  geom_line()+
  theme_bw(base_size=20)+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30))


# body fat dynamics
bf = bcd %>%
  filter(Date<as.POSIXct(strptime(as.character("2019"), format="%Y")))%>%
  select(tid,BF)%>%
  mutate(count=1)%>%
  group_by(tid) %>%
  mutate(ticker = cumsum(count))%>%
  spread(ticker,BF)%>%
  select(-count)
# calculate dynamics
bf$diff=NA
for (i in 1:nrow(bf)) {
  if(!is.na(bf$"2"[i])){
    if(!is.na(bf$"3"[i])){
      bf$diff[i] = bf$"1"[i] - bf$"3"[i]}
    else{bf$diff[i] = bf$"1"[i] - bf$"2"[i]}}
  else{bf$diff[i] = NA}}

hist(bf$diff,breaks=10)
mean(bf$diff,na.rm=T)
sd(bf$diff,na.rm=T)

# body mass dynamics
bm = bcd %>%
  filter(Date<as.POSIXct(strptime(as.character("2019"), format="%Y")))%>%
  select(tid,BM)%>%
  mutate(count=1)%>%
  group_by(tid) %>%
  mutate(ticker = cumsum(count))%>%
  spread(ticker,BM)%>%
  select(-count)
# calculate dynamics
bm$diff=NA
for (i in 1:nrow(bm)) {
  if(!is.na(bm$"2"[i])){
    if(!is.na(bm$"3"[i])){
      bm$diff[i] = bm$"1"[i] - bm$"3"[i]}
    else{bm$diff[i] = bm$"1"[i] - bm$"2"[i]}}
  else{bm$diff[i] = NA}}

hist(bm$diff,breaks=10)
mean(bm$diff,na.rm=T)
(sd(bm$diff,na.rm=T) / sqrt(length(bm$diff))) * 1.96

bootfun(bm$diff)


hist(replicate(100000, mean(sample(bf$diff, replace=TRUE), na.rm=TRUE)),breaks=100)





############################
# seasonal home ranges! ####
############################
C_sf = C_sample
C = as(C_sample, Class = "Spatial")
# here we want to keep everyihtn outside of Carolyns
# delinateaion area to make the kernels

# add coords in meters
C$X = coordinates(C)[,1]
C$Y = coordinates(C)[,2]



## testing
caribou = "WL-13-01"
i = which(IDS==caribou)
s = which(Seasons=="c")
year=2014
ci<-C[which(C$ID==IDS[i] & C$Season==Seasons[s] & C$Year==year),]
nrow(ci)
plot(ci)
writeOGR(ci,"P:/Projects/Permafrost/Analysis/Data/AnalysisData",caribou,driver="ESRI Shapefile")



# loop over individuals and seasons
# table(C$ID,C$Season,C$Year) 
Seasons=unique(C$Season)
CYears=unique(C$Year)
IDS=unique(C$ID)

# home range and availability lists
CarSA_List=list()
AvailList=list()
CHRs=list()
HRTab = data.frame()
tab_counter = 0

for (i in 1:length(IDS)){
  for (s in 1:length(Seasons)){
    ci<-C[which(C$ID==IDS[i] & C$Season==Seasons[s]),]
    ciyears=unique(ci$Year)
    
    for (y in 1:length(ciyears)){
      ci<-C[which(C$ID==IDS[i] & C$Season==Seasons[s] & C$Year==ciyears[y]),]
      
      if(nrow(ci)!=0){
        ciSA = ci[which(!is.na(over(ci,SA)[,1])),]  # remove points outside the study area
        
        # get the points for this caribou in this season
        if(nrow(ciSA) > (0.5 * nrow(ci)) & nrow(ciSA) > 30){  # if not enough points then forget this individual
          print(paste(IDS[i],Seasons[s],ciyears[y], nrow(ciSA),round(nrow(ciSA)/nrow(ci),2), "in study area",sep=" "))
          tab_counter = tab_counter +1
          ciKern= NULL;ciHR=NULL;ciAvail=NULL
          
          ## ID
          fid = gsub("WL-", "20", IDS[i])
          fid = gsub("-", "_", fid)
          cilabel = paste(fid,Seasons[s],ciyears[y], sep="_")
          
          if (cilabel %in% c("2013_01_c_2014","2013_02_c_2014","2013_05_c_2014")){next}  ## this caribou data is fucked up in the calving season - all on one point
          
          # build kernel raster
          buffdist = 1 * (abs(raster::extent(ci)[1])-abs(raster::extent(ci)[2])) # extent set by caribou
          cibuffer = buffer(ci, width=buffdist)
          r <- raster();extent(r)=extent(cibuffer);res(r)=100  # empty raster to populate and set extent
          ciKern <- ks::kde(ci@data[c("X","Y")],
                            gridsize=c(dim(r)[2],dim(r)[1]),
                            xmin=c(extent(r)@xmin,extent(r)@ymin),
                            xmax=c(extent(r)@xmax,extent(r)@ymax))
          ciKern <- raster(ciKern,crs=proj4string(C))
          
          # hr poly
          ciHR=raster.vol(ciKern,p=0.95)
          ciHR=coordinates(rasterToContour(ciHR))
          ciHR=list(Polygons(lapply(ciHR[[1]],Polygon),as.character(cilabel)))
          ciHR<-SpatialPolygonsDataFrame(SpatialPolygons(ciHR),
                                         data=data.frame(row.names=as.character(cilabel),SID = cilabel))
          proj4string(ciHR)=proj4string(C)
          
          # ## plottting
          plot(ciHR,main=cilabel);plot(ci,add=T,col="red")
          plot(ciSA,add=T,col="green");plot(SA,add=T,border="orange")
          
          # join to landscape stack
          ciSA = ciSA[which(!is.na(over(ciSA,ciHR))),] # get rid of points outside
          ciSA = cbind(ciSA,raster::extract(hstack,ciSA))
          ciSA$fid = fid
          ciSA$cilabel = cilabel
          
          ## availablity
          ciAvail = AV[which(!is.na(over(AV,ciHR))),]  # removed areas outside the study area
          ciAvail$fid = fid
          ciAvail$cilabel = cilabel
          ciAvail$Season = as.factor(Seasons[s])
          ciAvail$Year = ciyears[y]
          
          ## populate home range table
          HRTab[tab_counter,"fid"] = fid
          HRTab$Season[tab_counter] = Seasons[s]
          HRTab$Year[tab_counter] = ciyears[y]
          HRTab$cilabel[tab_counter]  = cilabel
          HRTab$HRS[tab_counter] = raster::area(gIntersection(SA,ciHR)) ##nrow(ciAvail))  ## want to have teh area interesecting the study area
          HRTab$decid[tab_counter] = mean(ciAvail$decid,na.rm = T)
          HRTab$decid_count[tab_counter] = length(ciAvail$decid[ciAvail$decid>50])
          HRTab$con[tab_counter] = mean(ciAvail$con,na.rm = T)
          HRTab$conMed[tab_counter] = median(ciAvail$con,na.rm = T)
          HRTab$con_count[tab_counter] = length(ciAvail$con[ciAvail$con>50])
          #wetlands
          wettab = table(ciAvail$wet)
          HRTab$up[tab_counter] = ifelse(is.na(wettab["0"]),0,wettab["0"])
          HRTab$peat[tab_counter] = ifelse(is.na(wettab["1"]),0,wettab["1"])
          HRTab$fen[tab_counter] = ifelse(is.na(wettab["2"]),0,wettab["2"])
          HRTab$pp[tab_counter] = ifelse(is.na(wettab["3"]),0,wettab["3"])
          HRTab$water[tab_counter] = ifelse(is.na(wettab["4"]),0,wettab["4"])
          
          ## also add to polygons for Q
          ciHR@data$fid = fid
          ciHR@data$Season = Seasons[s]
          ciHR@data$Year = ciyears[y]
          ciHR@data$cilabel = cilabel
          ciHR@data$con = mean(ciAvail$con,na.rm = T)
          ciHR@data$conMed = median(ciAvail$con,na.rm = T)
          ciHR@data$con_count = length(ciAvail$con[ciAvail$con>50])
          ciHR@data$decid = mean(ciAvail$decid,na.rm = T)
          ciHR@data$decid_count = length(ciAvail$decid[ciAvail$decid>50])
          ciHR@data$up = ifelse(is.na(wettab["0"]),0,wettab["0"])
          ciHR@data$peat =ifelse(is.na(wettab["1"]),0,wettab["1"])
          ciHR@data$fen = ifelse(is.na(wettab["2"]),0,wettab["2"])
          ciHR@data$pp = ifelse(is.na(wettab["3"]),0,wettab["3"])
          ciHR@data$water = ifelse(is.na(wettab["4"]),0,wettab["4"])
          
          ## movement
          datm = ci@data
          datm = datm[which(!is.na(datm$Date)),]
          dat_track = amt::steps(amt::make_track(datm,X,Y,DateTime))
          dat_sl = dat_track$sl_
          ## movement outputs
          HRTab$SLmean[tab_counter] = mean(dat_sl)
          HRTab$SLsd[tab_counter] = sd(dat_sl)
          HRTab$SLMax[tab_counter] = max(dat_sl)
          HRTab$SLLow[tab_counter] = quantile(dat_sl, probs=c(0.2))
          HRTab$SLHigh[tab_counter] = quantile(dat_sl, probs=c(0.8))
          
          ## get distribution of step lenghts
          dat_sl = dat_sl[which(dat_sl != 0 )]
          sld = list()
          sld[["exp"]] <- fitdistr(dat_sl, "exponential")
          exp_pv =ks.test(dat_sl, "pexp", sld[["exp"]]$estimate)[["p.value"]]
          sld[["gamma"]] <- fitdistr(dat_sl[dat_sl!=0], "gamma", list(shape = 1, rate = 0.1),lower = 0.1)
          gam_pv = ks.test(dat_sl, "pgamma", sld[["gamma"]]$estimate[1], sld[["gamma"]]$estimate[2])[["p.value"]]
          sld[["lnorm"]] <- fitdistr(dat_sl, "lognormal")
          lnorm_pv = ks.test(dat_sl, "plnorm", sld[["lnorm"]]$estimate[1], sld[["lnorm"]]$estimate[2])[["p.value"]]
          ## movemetn distribution parameters
          HRTab$SLDist[tab_counter] = names(sld)[which.max(c(exp_pv,gam_pv,lnorm_pv))]
          HRTab$SL1[tab_counter]= sld[[which.max(c(exp_pv,gam_pv,lnorm_pv))]]$estimate[1]
          HRTab$SL2[tab_counter] = try(sld[[which.max(c(exp_pv,gam_pv,lnorm_pv))]]$estimate[2])
          HRTab$TAsd[tab_counter] =sd(dat_track$ta_,na.rm=T)
          
          ## winter body condition
          if (Seasons[s] =='w') {
            if( length ( bf$diff[which(bf$tid == fid)] ) != 0 ){
              HRTab$bf[tab_counter] = bf$diff[which(bf$tid == fid)]}
            if( length ( bm$diff[which(bm$tid == fid)] ) != 0 ){
              HRTab$bm[tab_counter] = bm$diff[which(bm$tid == fid)]}
            if( length ( bf$'1'[which(bf$tid == fid)] ) != 0 ){
              HRTab$fat[tab_counter] = bf$'1'[which(bf$tid == fid)]}
            if( length ( bm$'1'[which(bm$tid == fid)] ) != 0 ){
              HRTab$mass[tab_counter] = bm$'1'[which(bm$tid == fid)]}
          }
          
          # add to datalists
          CHRs[[cilabel]] <- ciHR
          AvailList[[cilabel]] = ciAvail
          CarSA_List[[cilabel]] = ciSA
          
          # # #plotting
          # plot(ciHR);plot(ciAvail,add=T);plot(ci,add=T,col="red")
          # plot(ciSA,add=T,col="green");plot(SA,add=T,border="orange")
        }
      }
    }
  }
}

# i = which(IDS=="WL-17-22")
# s = which(Seasons=="w")
ci<-C[which(C$ID==IDS[i] & C$Season==Seasons[s]),]
# writeOGR(ci,"P:/Projects/Permafrost/Analysis/Data/AnalysisData","2017_22_w",driver="ESRI Shapefile")nrow(HRTab[which(HRTab$Season == 'w'),])

FullAvail <- do.call(rbind,AvailList)
FullAvail$Use=as.factor(0)
FullUsed <- do.call(rbind,CarSA_List) 
FullUsed$Use=as.factor(1)


# export 
allHR =  do.call(rbind, CHRs) 
names(allHR@data)[1] = "label"
# plot(allHR)
# View(allHR@data)
writeOGR(allHR,"P:/Projects/Permafrost/Analysis/Data/AnalysisData","allHR",driver="ESRI Shapefile")


#################################

