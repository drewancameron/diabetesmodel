### Model for the count of diabetes per SA1 unit over all Australia

setwd("~/VIRTUAL_WA/DIABETES/")

library(terra)
library(MASS)
library(Matrix)
library(truncnorm)
library(TMB)
library(sparseMVN)
library(readxl)
library(spdep)

### Read in and prepare data

sa1_2021 <- vect("/mnt/Z/ewan/DIABETES/sa1_2021_wcovsx.shp")
sa1_2021 <- sa1_2021[sa1_2021$STE_NAME21!="Other Territories",]

sa1_coords <- crds(centroids(sa1_2021))

# Read in RA code table and impute by nearest neighbours etc
ra_table <- read.csv("SA1_RA_correspondence.csv",skip=11,header = FALSE)

lookups <- match(sa1_2021$SA1_CODE21,ra_table$V1)
racodes <- apply(cbind(ra_table$V2,ra_table$V3,ra_table$V4,ra_table$V5,ra_table$V6,ra_table$V7,ra_table$V8)[lookups,],1,which.max)
racodes[which(apply(cbind(ra_table$V2,ra_table$V3,ra_table$V4,ra_table$V5,ra_table$V6,ra_table$V7,ra_table$V8)[lookups,],1,max)==0)] <- NA

nb <- poly2nb(sf::st_as_sf(sa1_2021))

to.infill <- which(is.na(racodes))
for (i in 1:length(to.infill)) {
  av.code <- round(mean(racodes[nb[[to.infill[i]]]],na.rm=TRUE))
  racodes[to.infill[i]] <- av.code
}

mean.with.na <- function(x) {mean(x,na.rm=TRUE)}
mean.ras <- aggregate(racodes,list(match(sa1_2021$SA2_CODE21,unique(sa1_2021$SA2_CODE21))),mean.with.na)
mean.ras <- mean.ras[match(sa1_2021$SA2_CODE21,unique(sa1_2021$SA2_CODE21)),2]

to.infill <- which(is.na(racodes))
racodes[to.infill] <- round(mean.ras[to.infill])

# Read in IRSAD code table and impute by nearest neighbours etc
irsad_table <- read.csv("SA1_IRSAD_lookup.csv.csv",skip=11,header = FALSE)

lookups <- match(sa1_2021$SA1_CODE21,irsad_table$V1)
irsad_codes <- apply(cbind(irsad_table$V2,irsad_table$V3,irsad_table$V4,irsad_table$V5,irsad_table$V6,irsad_table$V7,irsad_table$V8,irsad_table$V9,irsad_table$V10,irsad_table$V11,irsad_table$V12,irsad_table$V13)[lookups,],1,which.max)

### Read in CENSUS DATA : SA1 x HDIAP

filename_SA1 <- "SA1_HDIAP_final_0_9.csv"
filename_SA2 <- "SA2_HDIAP_final_0_9.csv"
filename_SA3 <- "SA3_HDIAP_final_0_9.csv"
filename_SA4 <- "SA4_HDIAP_final_0_9.csv"
filename_STE <- "STE_HDIAP_final_0_9.csv"

sa2_corespondences <- read.csv("CG_SA2_2016_SA2_2021.csv")
sa2_2016_young <- read.csv("T1D counts at SA2 in 2021.csv")

# SA1 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=12,header = FALSE)
SA1_num_dwellings_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_num_dwellings_by_DATUM <- cbind(SA1_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA1_num_dwellings_by_DATUM) <- xcolnames
SA1_num_dwellings_by_DATUM <- as.data.frame(SA1_num_dwellings_by_DATUM)
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[-which(is.na(SA1_num_dwellings_by_DATUM$SA1)),]
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[match(as.numeric(unique(sa1_2021$SA1_CODE21)),SA1_num_dwellings_by_DATUM$SA1),]
SA1_num_dwellings_by_DATUM <- as.matrix(SA1_num_dwellings_by_DATUM)

# SA2 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=12,header = FALSE)
SA2_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA2_num_dwellings_by_DATUM <- cbind(SA2_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(SA2_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA2_num_dwellings_by_DATUM) <- xcolnames
SA2_num_dwellings_by_DATUM <- as.data.frame(SA2_num_dwellings_by_DATUM)
SA2_num_dwellings_by_DATUM <- SA2_num_dwellings_by_DATUM[match((unique(sa1_2021$SA2_NAME21)),SA2_num_dwellings_by_DATUM$SA2),]
SA2_num_dwellings_by_DATUM[,1] <- 1:length(SA2_num_dwellings_by_DATUM[,1])
SA2_num_dwellings_by_DATUM[,2] <- as.numeric(SA2_num_dwellings_by_DATUM[,2])
SA2_num_dwellings_by_DATUM[,3] <- as.numeric(SA2_num_dwellings_by_DATUM[,3])
SA2_num_dwellings_by_DATUM[,4] <- as.numeric(SA2_num_dwellings_by_DATUM[,4])
SA2_num_dwellings_by_DATUM[,5] <- as.numeric(SA2_num_dwellings_by_DATUM[,5])
SA2_num_dwellings_by_DATUM <- as.matrix(SA2_num_dwellings_by_DATUM)

# SA3 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=12,header = FALSE)
SA3_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA3_num_dwellings_by_DATUM <- cbind(SA3_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA3",ydata[-1][1:(length(SA3_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA3_num_dwellings_by_DATUM) <- xcolnames
SA3_num_dwellings_by_DATUM <- as.data.frame(SA3_num_dwellings_by_DATUM)
SA3_num_dwellings_by_DATUM <- SA3_num_dwellings_by_DATUM[match((unique(sa1_2021$SA3_NAME21)),SA3_num_dwellings_by_DATUM$SA3),]
SA3_num_dwellings_by_DATUM[,1] <- 1:length(SA3_num_dwellings_by_DATUM[,1])
SA3_num_dwellings_by_DATUM[,2] <- as.numeric(SA3_num_dwellings_by_DATUM[,2])
SA3_num_dwellings_by_DATUM[,3] <- as.numeric(SA3_num_dwellings_by_DATUM[,3])
SA3_num_dwellings_by_DATUM[,4] <- as.numeric(SA3_num_dwellings_by_DATUM[,4])
SA3_num_dwellings_by_DATUM[,5] <- as.numeric(SA3_num_dwellings_by_DATUM[,5])
SA3_num_dwellings_by_DATUM <- as.matrix(SA3_num_dwellings_by_DATUM)

# SA4 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=12,header = FALSE)
SA4_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA4_num_dwellings_by_DATUM <- cbind(SA4_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA4",ydata[-1][1:(length(SA4_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA4_num_dwellings_by_DATUM) <- xcolnames
SA4_num_dwellings_by_DATUM <- as.data.frame(SA4_num_dwellings_by_DATUM)
SA4_num_dwellings_by_DATUM <- SA4_num_dwellings_by_DATUM[match((unique(sa1_2021$SA4_NAME21)),SA4_num_dwellings_by_DATUM$SA4),]
SA4_num_dwellings_by_DATUM[,1] <- 1:length(SA4_num_dwellings_by_DATUM[,1])
SA4_num_dwellings_by_DATUM[,2] <- as.numeric(SA4_num_dwellings_by_DATUM[,2])
SA4_num_dwellings_by_DATUM[,3] <- as.numeric(SA4_num_dwellings_by_DATUM[,3])
SA4_num_dwellings_by_DATUM[,4] <- as.numeric(SA4_num_dwellings_by_DATUM[,4])
SA4_num_dwellings_by_DATUM[,5] <- as.numeric(SA4_num_dwellings_by_DATUM[,5])
SA4_num_dwellings_by_DATUM <- as.matrix(SA4_num_dwellings_by_DATUM)

# STE level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=12,header = FALSE)
STE_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("STE_num_dwellings_by_DATUM <- cbind(STE_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("STE",ydata[-1][1:(length(STE_num_dwellings_by_DATUM[1,])-1)]))
colnames(STE_num_dwellings_by_DATUM) <- xcolnames
STE_num_dwellings_by_DATUM <- as.data.frame(STE_num_dwellings_by_DATUM)
STE_num_dwellings_by_DATUM <- STE_num_dwellings_by_DATUM[match((unique(sa1_2021$STE_NAME21)),STE_num_dwellings_by_DATUM$STE),]
STE_num_dwellings_by_DATUM[,1] <- 1:length(STE_num_dwellings_by_DATUM[,1])
STE_num_dwellings_by_DATUM[,2] <- as.numeric(STE_num_dwellings_by_DATUM[,2])
STE_num_dwellings_by_DATUM[,3] <- as.numeric(STE_num_dwellings_by_DATUM[,3])
STE_num_dwellings_by_DATUM[,4] <- as.numeric(STE_num_dwellings_by_DATUM[,4])
STE_num_dwellings_by_DATUM[,5] <- as.numeric(STE_num_dwellings_by_DATUM[,5])
STE_num_dwellings_by_DATUM <- as.matrix(STE_num_dwellings_by_DATUM)

sa1sa2codes <- match(sa1_2021$SA2_NAME21,unique(sa1_2021$SA2_NAME21))
sa1sa3codes <- match(sa1_2021$SA3_NAME21,unique(sa1_2021$SA3_NAME21))
sa1sa4codes <- match(sa1_2021$SA4_NAME21,unique(sa1_2021$SA4_NAME21))
sa1stecodes <- match(sa1_2021$STE_NAME21,unique(sa1_2021$STE_NAME21))

M <- dim(SA1_num_dwellings_by_DATUM)[2]-2

SA1_num_dwellings_by_DATUM_obs <- SA1_num_dwellings_by_DATUM[,-c(1,M+2)]
SA1_num_dwellings_obs <- SA1_num_dwellings_by_DATUM[,M+2]
SA2_num_dwellings_by_DATUM_obs <- SA2_num_dwellings_by_DATUM[,-c(1,M+2)]
SA2_num_dwellings_obs <- SA2_num_dwellings_by_DATUM[,M+2]
SA3_num_dwellings_by_DATUM_obs <- SA3_num_dwellings_by_DATUM[,-c(1,M+2)]
SA3_num_dwellings_obs <- SA3_num_dwellings_by_DATUM[,M+2]
SA4_num_dwellings_by_DATUM_obs <- SA4_num_dwellings_by_DATUM[,-c(1,M+2)]
SA4_num_dwellings_obs <- SA4_num_dwellings_by_DATUM[,M+2]
STE_num_dwellings_by_DATUM_obs <- STE_num_dwellings_by_DATUM[,-c(1,M+2)]
STE_num_dwellings_obs <- STE_num_dwellings_by_DATUM[,M+2]

# Load 10-14 EXTRA

filename_SA1 <- "SA1_HDIAP_final_10_14.csv"
filename_SA2 <- "SA2_HDIAP_final_10_14.csv"
filename_SA3 <- "SA3_HDIAP_final_10_14.csv"
filename_SA4 <- "SA4_HDIAP_final_10_14.csv"
filename_STE <- "STE_HDIAP_final_10_14.csv"

# SA1 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=12,header = FALSE)
SA1_num_dwellings_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_num_dwellings_by_DATUM <- cbind(SA1_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA1_num_dwellings_by_DATUM) <- xcolnames
SA1_num_dwellings_by_DATUM <- as.data.frame(SA1_num_dwellings_by_DATUM)
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[-which(is.na(SA1_num_dwellings_by_DATUM$SA1)),]
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[match(as.numeric(unique(sa1_2021$SA1_CODE21)),SA1_num_dwellings_by_DATUM$SA1),]
SA1_num_dwellings_by_DATUM <- as.matrix(SA1_num_dwellings_by_DATUM)

# SA2 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=12,header = FALSE)
SA2_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA2_num_dwellings_by_DATUM <- cbind(SA2_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(SA2_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA2_num_dwellings_by_DATUM) <- xcolnames
SA2_num_dwellings_by_DATUM <- as.data.frame(SA2_num_dwellings_by_DATUM)
SA2_num_dwellings_by_DATUM <- SA2_num_dwellings_by_DATUM[match((unique(sa1_2021$SA2_NAME21)),SA2_num_dwellings_by_DATUM$SA2),]
SA2_num_dwellings_by_DATUM[,1] <- 1:length(SA2_num_dwellings_by_DATUM[,1])
SA2_num_dwellings_by_DATUM[,2] <- as.numeric(SA2_num_dwellings_by_DATUM[,2])
SA2_num_dwellings_by_DATUM[,3] <- as.numeric(SA2_num_dwellings_by_DATUM[,3])
SA2_num_dwellings_by_DATUM[,4] <- as.numeric(SA2_num_dwellings_by_DATUM[,4])
SA2_num_dwellings_by_DATUM[,5] <- as.numeric(SA2_num_dwellings_by_DATUM[,5])
SA2_num_dwellings_by_DATUM <- as.matrix(SA2_num_dwellings_by_DATUM)

# SA3 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=12,header = FALSE)
SA3_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA3_num_dwellings_by_DATUM <- cbind(SA3_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA3",ydata[-1][1:(length(SA3_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA3_num_dwellings_by_DATUM) <- xcolnames
SA3_num_dwellings_by_DATUM <- as.data.frame(SA3_num_dwellings_by_DATUM)
SA3_num_dwellings_by_DATUM <- SA3_num_dwellings_by_DATUM[match((unique(sa1_2021$SA3_NAME21)),SA3_num_dwellings_by_DATUM$SA3),]
SA3_num_dwellings_by_DATUM[,1] <- 1:length(SA3_num_dwellings_by_DATUM[,1])
SA3_num_dwellings_by_DATUM[,2] <- as.numeric(SA3_num_dwellings_by_DATUM[,2])
SA3_num_dwellings_by_DATUM[,3] <- as.numeric(SA3_num_dwellings_by_DATUM[,3])
SA3_num_dwellings_by_DATUM[,4] <- as.numeric(SA3_num_dwellings_by_DATUM[,4])
SA3_num_dwellings_by_DATUM[,5] <- as.numeric(SA3_num_dwellings_by_DATUM[,5])
SA3_num_dwellings_by_DATUM <- as.matrix(SA3_num_dwellings_by_DATUM)

# SA4 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=12,header = FALSE)
SA4_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA4_num_dwellings_by_DATUM <- cbind(SA4_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA4",ydata[-1][1:(length(SA4_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA4_num_dwellings_by_DATUM) <- xcolnames
SA4_num_dwellings_by_DATUM <- as.data.frame(SA4_num_dwellings_by_DATUM)
SA4_num_dwellings_by_DATUM <- SA4_num_dwellings_by_DATUM[match((unique(sa1_2021$SA4_NAME21)),SA4_num_dwellings_by_DATUM$SA4),]
SA4_num_dwellings_by_DATUM[,1] <- 1:length(SA4_num_dwellings_by_DATUM[,1])
SA4_num_dwellings_by_DATUM[,2] <- as.numeric(SA4_num_dwellings_by_DATUM[,2])
SA4_num_dwellings_by_DATUM[,3] <- as.numeric(SA4_num_dwellings_by_DATUM[,3])
SA4_num_dwellings_by_DATUM[,4] <- as.numeric(SA4_num_dwellings_by_DATUM[,4])
SA4_num_dwellings_by_DATUM[,5] <- as.numeric(SA4_num_dwellings_by_DATUM[,5])
SA4_num_dwellings_by_DATUM <- as.matrix(SA4_num_dwellings_by_DATUM)

# STE level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=12,header = FALSE)
STE_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("STE_num_dwellings_by_DATUM <- cbind(STE_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("STE",ydata[-1][1:(length(STE_num_dwellings_by_DATUM[1,])-1)]))
colnames(STE_num_dwellings_by_DATUM) <- xcolnames
STE_num_dwellings_by_DATUM <- as.data.frame(STE_num_dwellings_by_DATUM)
STE_num_dwellings_by_DATUM <- STE_num_dwellings_by_DATUM[match((unique(sa1_2021$STE_NAME21)),STE_num_dwellings_by_DATUM$STE),]
STE_num_dwellings_by_DATUM[,1] <- 1:length(STE_num_dwellings_by_DATUM[,1])
STE_num_dwellings_by_DATUM[,2] <- as.numeric(STE_num_dwellings_by_DATUM[,2])
STE_num_dwellings_by_DATUM[,3] <- as.numeric(STE_num_dwellings_by_DATUM[,3])
STE_num_dwellings_by_DATUM[,4] <- as.numeric(STE_num_dwellings_by_DATUM[,4])
STE_num_dwellings_by_DATUM[,5] <- as.numeric(STE_num_dwellings_by_DATUM[,5])
STE_num_dwellings_by_DATUM <- as.matrix(STE_num_dwellings_by_DATUM)

SA1_num_dwellings_by_DATUM_obsX <- SA1_num_dwellings_by_DATUM[,-c(1,M+2)]
SA1_num_dwellings_obsX <- SA1_num_dwellings_by_DATUM[,M+2]
SA2_num_dwellings_by_DATUM_obsX <- SA2_num_dwellings_by_DATUM[,-c(1,M+2)]
SA2_num_dwellings_obsX <- SA2_num_dwellings_by_DATUM[,M+2]
SA3_num_dwellings_by_DATUM_obsX <- SA3_num_dwellings_by_DATUM[,-c(1,M+2)]
SA3_num_dwellings_obsX <- SA3_num_dwellings_by_DATUM[,M+2]
SA4_num_dwellings_by_DATUM_obsX <- SA4_num_dwellings_by_DATUM[,-c(1,M+2)]
SA4_num_dwellings_obsX <- SA4_num_dwellings_by_DATUM[,M+2]
STE_num_dwellings_by_DATUM_obsX <- STE_num_dwellings_by_DATUM[,-c(1,M+2)]
STE_num_dwellings_obsX <- STE_num_dwellings_by_DATUM[,M+2]

# Load 15-19 EXTRA

filename_SA1 <- "SA1_HDIAP_final_15_19.csv"
filename_SA2 <- "SA2_HDIAP_final_15_19.csv"
filename_SA3 <- "SA3_HDIAP_final_15_19.csv"
filename_SA4 <- "SA4_HDIAP_final_15_19.csv"
filename_STE <- "STE_HDIAP_final_15_19.csv"

# SA1 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=12,header = FALSE)
SA1_num_dwellings_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_num_dwellings_by_DATUM <- cbind(SA1_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA1_num_dwellings_by_DATUM) <- xcolnames
SA1_num_dwellings_by_DATUM <- as.data.frame(SA1_num_dwellings_by_DATUM)
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[-which(is.na(SA1_num_dwellings_by_DATUM$SA1)),]
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[match(as.numeric(unique(sa1_2021$SA1_CODE21)),SA1_num_dwellings_by_DATUM$SA1),]
SA1_num_dwellings_by_DATUM <- as.matrix(SA1_num_dwellings_by_DATUM)

# SA2 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=12,header = FALSE)
SA2_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA2_num_dwellings_by_DATUM <- cbind(SA2_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(SA2_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA2_num_dwellings_by_DATUM) <- xcolnames
SA2_num_dwellings_by_DATUM <- as.data.frame(SA2_num_dwellings_by_DATUM)
SA2_num_dwellings_by_DATUM <- SA2_num_dwellings_by_DATUM[match((unique(sa1_2021$SA2_NAME21)),SA2_num_dwellings_by_DATUM$SA2),]
SA2_num_dwellings_by_DATUM[,1] <- 1:length(SA2_num_dwellings_by_DATUM[,1])
SA2_num_dwellings_by_DATUM[,2] <- as.numeric(SA2_num_dwellings_by_DATUM[,2])
SA2_num_dwellings_by_DATUM[,3] <- as.numeric(SA2_num_dwellings_by_DATUM[,3])
SA2_num_dwellings_by_DATUM[,4] <- as.numeric(SA2_num_dwellings_by_DATUM[,4])
SA2_num_dwellings_by_DATUM[,5] <- as.numeric(SA2_num_dwellings_by_DATUM[,5])
SA2_num_dwellings_by_DATUM <- as.matrix(SA2_num_dwellings_by_DATUM)

# SA3 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=12,header = FALSE)
SA3_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA3_num_dwellings_by_DATUM <- cbind(SA3_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA3",ydata[-1][1:(length(SA3_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA3_num_dwellings_by_DATUM) <- xcolnames
SA3_num_dwellings_by_DATUM <- as.data.frame(SA3_num_dwellings_by_DATUM)
SA3_num_dwellings_by_DATUM <- SA3_num_dwellings_by_DATUM[match((unique(sa1_2021$SA3_NAME21)),SA3_num_dwellings_by_DATUM$SA3),]
SA3_num_dwellings_by_DATUM[,1] <- 1:length(SA3_num_dwellings_by_DATUM[,1])
SA3_num_dwellings_by_DATUM[,2] <- as.numeric(SA3_num_dwellings_by_DATUM[,2])
SA3_num_dwellings_by_DATUM[,3] <- as.numeric(SA3_num_dwellings_by_DATUM[,3])
SA3_num_dwellings_by_DATUM[,4] <- as.numeric(SA3_num_dwellings_by_DATUM[,4])
SA3_num_dwellings_by_DATUM[,5] <- as.numeric(SA3_num_dwellings_by_DATUM[,5])
SA3_num_dwellings_by_DATUM <- as.matrix(SA3_num_dwellings_by_DATUM)

# SA4 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=12,header = FALSE)
SA4_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA4_num_dwellings_by_DATUM <- cbind(SA4_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA4",ydata[-1][1:(length(SA4_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA4_num_dwellings_by_DATUM) <- xcolnames
SA4_num_dwellings_by_DATUM <- as.data.frame(SA4_num_dwellings_by_DATUM)
SA4_num_dwellings_by_DATUM <- SA4_num_dwellings_by_DATUM[match((unique(sa1_2021$SA4_NAME21)),SA4_num_dwellings_by_DATUM$SA4),]
SA4_num_dwellings_by_DATUM[,1] <- 1:length(SA4_num_dwellings_by_DATUM[,1])
SA4_num_dwellings_by_DATUM[,2] <- as.numeric(SA4_num_dwellings_by_DATUM[,2])
SA4_num_dwellings_by_DATUM[,3] <- as.numeric(SA4_num_dwellings_by_DATUM[,3])
SA4_num_dwellings_by_DATUM[,4] <- as.numeric(SA4_num_dwellings_by_DATUM[,4])
SA4_num_dwellings_by_DATUM[,5] <- as.numeric(SA4_num_dwellings_by_DATUM[,5])
SA4_num_dwellings_by_DATUM <- as.matrix(SA4_num_dwellings_by_DATUM)

# STE level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=12,header = FALSE)
STE_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("STE_num_dwellings_by_DATUM <- cbind(STE_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("STE",ydata[-1][1:(length(STE_num_dwellings_by_DATUM[1,])-1)]))
colnames(STE_num_dwellings_by_DATUM) <- xcolnames
STE_num_dwellings_by_DATUM <- as.data.frame(STE_num_dwellings_by_DATUM)
STE_num_dwellings_by_DATUM <- STE_num_dwellings_by_DATUM[match((unique(sa1_2021$STE_NAME21)),STE_num_dwellings_by_DATUM$STE),]
STE_num_dwellings_by_DATUM[,1] <- 1:length(STE_num_dwellings_by_DATUM[,1])
STE_num_dwellings_by_DATUM[,2] <- as.numeric(STE_num_dwellings_by_DATUM[,2])
STE_num_dwellings_by_DATUM[,3] <- as.numeric(STE_num_dwellings_by_DATUM[,3])
STE_num_dwellings_by_DATUM[,4] <- as.numeric(STE_num_dwellings_by_DATUM[,4])
STE_num_dwellings_by_DATUM[,5] <- as.numeric(STE_num_dwellings_by_DATUM[,5])
STE_num_dwellings_by_DATUM <- as.matrix(STE_num_dwellings_by_DATUM)

SA1_num_dwellings_by_DATUM_obsY <- SA1_num_dwellings_by_DATUM[,-c(1,M+2)]
SA1_num_dwellings_obsY <- SA1_num_dwellings_by_DATUM[,M+2]
SA2_num_dwellings_by_DATUM_obsY <- SA2_num_dwellings_by_DATUM[,-c(1,M+2)]
SA2_num_dwellings_obsY <- SA2_num_dwellings_by_DATUM[,M+2]
SA3_num_dwellings_by_DATUM_obsY <- SA3_num_dwellings_by_DATUM[,-c(1,M+2)]
SA3_num_dwellings_obsY <- SA3_num_dwellings_by_DATUM[,M+2]
SA4_num_dwellings_by_DATUM_obsY <- SA4_num_dwellings_by_DATUM[,-c(1,M+2)]
SA4_num_dwellings_obsY <- SA4_num_dwellings_by_DATUM[,M+2]
STE_num_dwellings_by_DATUM_obsY <- STE_num_dwellings_by_DATUM[,-c(1,M+2)]
STE_num_dwellings_obsY <- STE_num_dwellings_by_DATUM[,M+2]

# IRSAD 0-9
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/IRSAD_HDIAP_09.csv",skip=12,header = FALSE)[1:11,]
IRSAD_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("IRSAD_by_DATUM <- cbind(IRSAD_by_DATUM,as.numeric(xdata$V",i,"))")))
}
IRSAD_by_DATUM <- as.matrix(IRSAD_by_DATUM[,2:4])
IRSAD_by_DATUM_09 <- IRSAD_by_DATUM

# IRSAD 10-14
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/IRSAD_HDIAP_1014.csv",skip=12,header = FALSE)[1:11,]
IRSAD_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("IRSAD_by_DATUM <- cbind(IRSAD_by_DATUM,as.numeric(xdata$V",i,"))")))
}
IRSAD_by_DATUM <- as.matrix(IRSAD_by_DATUM[,2:4])
IRSAD_by_DATUM_1014 <- IRSAD_by_DATUM

# IRSAD 15-19
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/IRSAD_HDIAP_1519.csv",skip=12,header = FALSE)[1:11,]
IRSAD_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("IRSAD_by_DATUM <- cbind(IRSAD_by_DATUM,as.numeric(xdata$V",i,"))")))
}
IRSAD_by_DATUM <- as.matrix(IRSAD_by_DATUM[,2:4])
IRSAD_by_DATUM_1519 <- IRSAD_by_DATUM

# RA 0-9
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/RA_HDIAP_09.csv",skip=12,header = FALSE)[1:5,]
RA_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("RA_by_DATUM <- cbind(RA_by_DATUM,as.numeric(xdata$V",i,"))")))
}
RA_by_DATUM <- as.matrix(RA_by_DATUM[,2:4])
RA_by_DATUM_09 <- RA_by_DATUM

# RA 10-14
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/RA_HDIAP_1014.csv",skip=12,header = FALSE)[1:5,]
RA_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("RA_by_DATUM <- cbind(RA_by_DATUM,as.numeric(xdata$V",i,"))")))
}
RA_by_DATUM <- as.matrix(RA_by_DATUM[,2:4])
RA_by_DATUM_1014 <- RA_by_DATUM

# RA 15-19
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/RA_HDIAP_1519.csv",skip=12,header = FALSE)[1:5,]
RA_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("RA_by_DATUM <- cbind(RA_by_DATUM,as.numeric(xdata$V",i,"))")))
}
RA_by_DATUM <- as.matrix(RA_by_DATUM[,2:4])
RA_by_DATUM_1519 <- RA_by_DATUM

### Construct geospatial mesh model

library(INLA)

sydney.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Sydney",],cutoff = 0.015,max.n=300)
melbourne.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Melbourne",],cutoff = 0.015,max.n=300)
brisbane.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Brisbane",],cutoff = 0.025,max.n=300)
adelaide.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Adelaide",],cutoff = 0.015,max.n=300)
perth.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Perth",],cutoff = 0.015,max.n=300)
darwin.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Darwin",],cutoff = 0.01,max.n=300)
hobart.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Hobart",],cutoff = 0.01,max.n=300)
canberra.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Australian Capital Territory",],cutoff = 0.01,max.n=300)
remainder.mesh <- inla.mesh.2d(sa1_coords[!(sa1_2021$GCC_NAME21 %in% c("Greater Sydney","Greater Melbourne","Greater Brisbane", "Greater Adelaide", "Greater Darwin", "Greater Hobart","Greater Perth","Australian Capital Territory")),],cutoff = 0.05,max.n=300)
xcoords <- rbind(sydney.mesh$loc,
                 melbourne.mesh$loc,
                 brisbane.mesh$loc,
                 adelaide.mesh$loc,
                 perth.mesh$loc,
                 darwin.mesh$loc,
                 hobart.mesh$loc,
                 canberra.mesh$loc,
                 remainder.mesh$loc)
aus_mesh_fine <- inla.mesh.2d(xcoords,cutoff = 0.01,max.n=1000)
aus_spde_fine <- (inla.spde2.matern(aus_mesh_fine,alpha=0.5)$param.inla)[c("M0","M1","M2")]
aus_A_fine <- inla.mesh.projector(aus_mesh_fine,sa1_coords)$proj$A

### Helper functions

library(extraDistr)

error_sd <- 2.0

logdiffexp <- function(y,x) {
  # x > y
  x+log(1-exp(y-x))
}
logsumexp <- function(y,x) {
  if (length(x)==1) {
    max(c(x,y))+log(1+exp(min(c(x,y))-max(c(x,y))))} else {
      mmax <- apply(cbind(x,y),1,max)
      mmin <- apply(cbind(x,y),1,min)
      mmax+log(1+exp(mmin-mmax))
    }
}

log_likelihood_diff_fn <- function(proposed_value,current_value,observed_value) {
  
  if (proposed_value==0 & observed_value>0) {
    proposed_likelihood <- NA
  } else if (proposed_value==0 & observed_value==0) {
    proposed_likelihood <- 0
  } else if (proposed_value>0 & observed_value==0) {
    proposed_likelihood <- pnorm(2.5,proposed_value,error_sd,log.p=TRUE)
  } else if (abs(proposed_value-observed_value) < 5 & abs(current_value-observed_value) < 5) {
    proposed_likelihood <- logdiffexp(pnorm(observed_value-0.5,proposed_value,error_sd,log.p = TRUE),pnorm(observed_value+0.5,proposed_value,error_sd,log.p = TRUE))
  } else {
    proposed_likelihood <- dnorm(observed_value,proposed_value,error_sd,log = TRUE)
  }
  
  if (current_value==0 & observed_value==0) {
    current_likelihood <- 0
  } else if (current_value>0 & observed_value==0) {
    current_likelihood <- pnorm(2.5,current_value,error_sd,log.p=TRUE)
  } else if (abs(proposed_value-observed_value) < 5 & abs(current_value-observed_value) < 5) {
    current_likelihood <- logdiffexp(pnorm(observed_value-0.5,current_value,error_sd,log.p = TRUE),pnorm(observed_value+0.5,current_value,error_sd,log.p = TRUE))
  } else {
    current_likelihood <- dnorm(observed_value,current_value,error_sd,log = TRUE)
  }
  
  return(proposed_likelihood - current_likelihood)
}

### MCMC Chain

N_SA1 <- length(sa1_2021)
true_SA1_HDIAP_AGE <- array(1,dim=c(N_SA1,3,3))
true_SA1_HDIAP_AGE[,,1] <- SA1_num_dwellings_by_DATUM_obs+1
true_SA1_HDIAP_AGE[,,2] <- SA1_num_dwellings_by_DATUM_obsX+1
true_SA1_HDIAP_AGE[,,3] <- SA1_num_dwellings_by_DATUM_obsY+1

true_SA1_HDIAP_09 <- true_SA1_HDIAP_AGE[,,1]
true_SA2_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1sa2codes),sum)[,-1]
true_SA3_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1sa3codes),sum)[,-1]
true_SA4_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1sa4codes),sum)[,-1]
true_STE_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1stecodes),sum)[,-1]

true_SA1_HDIAP_1014 <- true_SA1_HDIAP_AGE[,,2]
true_SA2_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1sa2codes),sum)[,-1]
true_SA3_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1sa3codes),sum)[,-1]
true_SA4_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1sa4codes),sum)[,-1]
true_STE_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1stecodes),sum)[,-1]

true_SA1_HDIAP_1519 <- true_SA1_HDIAP_AGE[,,3]
true_SA2_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1sa2codes),sum)[,-1]
true_SA3_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1sa3codes),sum)[,-1]
true_SA4_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1sa4codes),sum)[,-1]
true_STE_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1stecodes),sum)[,-1]

true_IRSAD_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(irsad_codes),sum)[,-1]
true_IRSAD_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(irsad_codes),sum)[,-1]
true_IRSAD_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(irsad_codes),sum)[,-1]

true_RA_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(racodes),sum)[,-1]
true_RA_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(racodes),sum)[,-1]
true_RA_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(racodes),sum)[,-1]

## Fit with spatial model

compile("inla_diabetes_allages.cpp")
dyn.load(dynlib("inla_diabetes_allages"))
ilogit <- function(x) {1/(1+exp(-x))}

count <- 0
prior_draw <- t(matrix(c(sum(STE_num_dwellings_by_DATUM_obs[,1])/(sum(STE_num_dwellings_by_DATUM_obs[,2])+sum(STE_num_dwellings_by_DATUM_obs[,1])),sum(STE_num_dwellings_by_DATUM_obsX[,1])/(sum(STE_num_dwellings_by_DATUM_obsX[,2])+sum(STE_num_dwellings_by_DATUM_obsX[,1])),sum(STE_num_dwellings_by_DATUM_obsY[,1])/(sum(STE_num_dwellings_by_DATUM_obsY[,2])+sum(STE_num_dwellings_by_DATUM_obsY[,1]))),ncol=N_SA1,nrow=3))

accepted_track <- array(0,dim=c(N_SA1,3,3))

for (z in 1:300000000) {
  
  if ((z%%1000000)==0) {
    count <- count + 1
    input.data <- list('NSA1'=N_SA1,
                         'spde'=aus_spde_fine,
                         'A'=aus_A_fine,
                         'npos' = true_SA1_HDIAP_AGE[,1,],
                         'nneg' = true_SA1_HDIAP_AGE[,2,],
                         'sa1_ses' = irsad_codes-1,
                         'sa1_ra' = racodes-1
    )

    if (count==1) {parameters <- list(
      'log_range_young'=0,
      'log_sd_young'=0,
      'log_range_old'=0,
      'log_sd_old'=0,
      'log_sd_ses_young'=0,
      'log_sd_ses_old'=0,
      'log_sd_ra_young'=0,
      'log_sd_ra_old'=0,
      'intercept_young'=-5,
      'intercept_old'=-2,
      'intercept_middle'=0,
      'logit_middle_age_prop'=0,
      'ses_effects_young'=rep(0,11),
      'ra_effects_young'=rep(0,5),
      'ses_effects_old'=rep(0,11),
      'ra_effects_old'=rep(0,5),
      'field_young'=rep(0,aus_mesh_fine$n),
      'field_old'=rep(0,aus_mesh_fine$n)
    )}

    obj <- MakeADFun(input.data,parameters,DLL = "inla_diabetes_allages",random=c('intercept_young','intercept_old','intercept_middle','logit_middle_age_prop','ses_effects_young','ra_effects_young','ses_effects_old','ra_effects_old','field_young','field_old'))
    obj$fn()
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    rep <- sdreport(obj,getJointPrecision = TRUE)

    xsample <- rmvn.sparse(10, c(rep$par.fixed, rep$par.random), Cholesky(rep$jointPrecision), prec = TRUE)
    colnames(xsample) <- names(c(rep$par.fixed, rep$par.random))

    parameters <- list(
      'log_range_young'=xsample[1,colnames(xsample)=="log_range_young"],
      'log_sd_young'=xsample[1,colnames(xsample)=="log_sd_young"],
      'log_range_old'=xsample[1,colnames(xsample)=="log_range_old"],
      'log_sd_old'=xsample[1,colnames(xsample)=="log_sd_old"],
      'log_sd_ses_young'=xsample[1,colnames(xsample)=="log_sd_ses_young"],
      'log_sd_ses_old'=xsample[1,colnames(xsample)=="log_sd_ses_old"],
      'log_sd_ra_young'=xsample[1,colnames(xsample)=="log_sd_ra_young"],
      'log_sd_ra_old'=xsample[1,colnames(xsample)=="log_sd_ra_old"],
      'intercept_young'=0,
      'intercept_old'=0,
      'intercept_middle'=0,
      'logit_middle_age_prop'=0,
      'ses_effects_young'=rep(0,11),
      'ra_effects_young'=rep(0,5),
      'ses_effects_old'=rep(0,11),
      'ra_effects_old'=rep(0,5),
      'field_young'=rep(0,aus_mesh_fine$n),
      'field_old'=rep(0,aus_mesh_fine$n)
    )

    objx <- MakeADFun(input.data,parameters,DLL = "inla_diabetes_allages")

    save(xsample,file=paste0("/mnt/Z/ewan/DIABETES/xsample_final",sprintf("%03i",count)))

    prior_draw <- objx$report(xsample[1,])$joint_field
    prior_draw[,1] <- ilogit(prior_draw[,1])
    prior_draw[,2] <- ilogit(prior_draw[,2])
    prior_draw[,3] <- ilogit(prior_draw[,3])

    xcol <- prior_draw[,1]
    q1 <- quantile(prior_draw[,1],0.1)
    q2 <- quantile(prior_draw[,1],0.9)
    xcol[xcol<q1] <- q1
    xcol[xcol>q2] <- q2
    xcol <- (xcol-min(xcol))/diff(range(xcol))
    plot(sa1_2021,border=NA,col=hsv((1-xcol)*0.666))

    xcol <- prior_draw[,3]
    q1 <- quantile(prior_draw[,3],0.1)
    q2 <- quantile(prior_draw[,3],0.9)
    xcol[xcol<q1] <- q1
    xcol[xcol>q2] <- q2
    xcol <- (xcol-min(xcol))/diff(range(xcol))
    plot(sa1_2021,border=NA,col=hsv((1-xcol)*0.666))
    
    cat(sum(z/10000000),"xx ",xsample[1,colnames(xsample)=="logit_middle_age_prop"],"\n")
    
    gc()
  }

  proposed_location <- cbind(sample(1:N_SA1,1),sample(1:M,1,prob=c(10,2,1)),sample(1:3,1))
  proposed_move <- sample(c(-1,1),1)
  
  current_likelihood_diff <- 0
  
  if ((true_SA1_HDIAP_AGE[proposed_location[,1],proposed_location[,2],proposed_location[,3]] > 0) | (proposed_move==1)) {
    
    # SA1 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA1_HDIAP_09[proposed_location[,1],proposed_location[,2]]
      observed_value <- SA1_num_dwellings_by_DATUM_obs[proposed_location[,1],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA1_HDIAP_1014[proposed_location[,1],proposed_location[,2]]
      observed_value <- SA1_num_dwellings_by_DATUM_obsX[proposed_location[,1],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA1_HDIAP_1519[proposed_location[,1],proposed_location[,2]]
      observed_value <- SA1_num_dwellings_by_DATUM_obsY[proposed_location[,1],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # SA2 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA2_HDIAP_09[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA2_num_dwellings_by_DATUM_obs[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA2_HDIAP_1014[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA2_num_dwellings_by_DATUM_obsX[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA2_HDIAP_1519[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA2_num_dwellings_by_DATUM_obsY[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # SA3 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA3_HDIAP_09[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA3_num_dwellings_by_DATUM_obs[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA3_HDIAP_1014[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA3_num_dwellings_by_DATUM_obsX[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA3_HDIAP_1519[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA3_num_dwellings_by_DATUM_obsY[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # SA4 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA4_HDIAP_09[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA4_num_dwellings_by_DATUM_obs[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA4_HDIAP_1014[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA4_num_dwellings_by_DATUM_obsX[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA4_HDIAP_1519[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA4_num_dwellings_by_DATUM_obsY[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # STE likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_STE_HDIAP_09[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- STE_num_dwellings_by_DATUM_obs[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_STE_HDIAP_1014[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- STE_num_dwellings_by_DATUM_obsX[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_STE_HDIAP_1519[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- STE_num_dwellings_by_DATUM_obsY[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # IRSAD likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_IRSAD_HDIAP_09[irsad_codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- IRSAD_by_DATUM_09[irsad_codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_IRSAD_HDIAP_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- IRSAD_by_DATUM_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_IRSAD_HDIAP_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- IRSAD_by_DATUM_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move

    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)

    # RA likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_RA_HDIAP_09[racodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- RA_by_DATUM_09[racodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_RA_HDIAP_1014[racodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- RA_by_DATUM_1014[racodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_RA_HDIAP_1519[racodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- RA_by_DATUM_1519[racodes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move

    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # Regression priors
    
    current_pos <- true_SA1_HDIAP_AGE[proposed_location[1],1,proposed_location[3]]
    current_neg <- true_SA1_HDIAP_AGE[proposed_location[1],2,proposed_location[3]]

    proposed_pos <- current_pos
    proposed_neg <- current_neg
    if (proposed_location[2]==1) {proposed_pos <- proposed_pos+proposed_move}
    if (proposed_location[2]==2) {proposed_neg <- proposed_neg+proposed_move}

    proposed_priors <- current_priors <- 0

    if ((proposed_location[2] %in% c(1,2)) & ((current_pos+current_neg)>0 & (proposed_pos+proposed_neg)>0)) {
      current_priors <- current_priors + dbinom(current_pos,(current_pos+current_neg),prior_draw[proposed_location[1],proposed_location[3]],log=TRUE)
      proposed_priors <- proposed_priors + dbinom(proposed_pos,(proposed_pos+proposed_neg),prior_draw[proposed_location[1],proposed_location[3]],log=TRUE)
    }

    current_likelihood_diff <- current_likelihood_diff + proposed_priors - current_priors
    
    # if accepted
    if (is.na(current_likelihood_diff) || current_likelihood_diff < log(runif(1))) {
      proposed_likelihood <- NA
    } else {
      
      accepted_track[proposed_location[,1],proposed_location[,2],proposed_location[,3]] <-     accepted_track[proposed_location[,1],proposed_location[,2],proposed_location[,3]] + 1
      
      true_SA1_HDIAP_AGE[proposed_location[,1],proposed_location[,2],proposed_location[,3]] <- true_SA1_HDIAP_AGE[proposed_location[,1],proposed_location[,2],proposed_location[,3]] + proposed_move
      
      if (proposed_location[,3]==1) {
        true_SA1_HDIAP_09[proposed_location[,1],proposed_location[,2]] <- true_SA1_HDIAP_09[proposed_location[,1],proposed_location[,2]] + proposed_move
        true_SA2_HDIAP_09[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] <- true_SA2_HDIAP_09[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA3_HDIAP_09[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] <- true_SA3_HDIAP_09[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA4_HDIAP_09[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] <- true_SA4_HDIAP_09[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_STE_HDIAP_09[sa1stecodes[proposed_location[,1]],proposed_location[,2]] <- true_STE_HDIAP_09[sa1stecodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_IRSAD_HDIAP_09[irsad_codes[proposed_location[,1]],proposed_location[,2]] <- true_IRSAD_HDIAP_09[irsad_codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_RA_HDIAP_09[racodes[proposed_location[,1]],proposed_location[,2]] <- true_RA_HDIAP_09[racodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
      }
      if (proposed_location[,3]==2) {
        true_SA1_HDIAP_1014[proposed_location[,1],proposed_location[,2]] <- true_SA1_HDIAP_1014[proposed_location[,1],proposed_location[,2]] + proposed_move
        true_SA2_HDIAP_1014[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] <- true_SA2_HDIAP_1014[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA3_HDIAP_1014[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] <- true_SA3_HDIAP_1014[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA4_HDIAP_1014[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] <- true_SA4_HDIAP_1014[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_STE_HDIAP_1014[sa1stecodes[proposed_location[,1]],proposed_location[,2]] <- true_STE_HDIAP_1014[sa1stecodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_IRSAD_HDIAP_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]] <- true_IRSAD_HDIAP_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_RA_HDIAP_1014[racodes[proposed_location[,1]],proposed_location[,2]] <- true_RA_HDIAP_1014[racodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
      }
      if (proposed_location[,3]==3) {
        true_SA1_HDIAP_1519[proposed_location[,1],proposed_location[,2]] <- true_SA1_HDIAP_1519[proposed_location[,1],proposed_location[,2]] + proposed_move
        true_SA2_HDIAP_1519[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] <- true_SA2_HDIAP_1519[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA3_HDIAP_1519[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] <- true_SA3_HDIAP_1519[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA4_HDIAP_1519[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] <- true_SA4_HDIAP_1519[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_STE_HDIAP_1519[sa1stecodes[proposed_location[,1]],proposed_location[,2]] <- true_STE_HDIAP_1519[sa1stecodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_IRSAD_HDIAP_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]] <- true_IRSAD_HDIAP_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_RA_HDIAP_1519[racodes[proposed_location[,1]],proposed_location[,2]] <- true_RA_HDIAP_1519[racodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
      }
      
    }
  }
  
  if ((z %% 10000)==1) {cat(sum(z/10000000),sum(true_SA1_HDIAP_AGE[,1,1])/sum(STE_num_dwellings_by_DATUM_obs[,1])," ",sum(true_SA1_HDIAP_AGE[,1,2])/sum(STE_num_dwellings_by_DATUM_obsX[,1])," ",sum(true_SA1_HDIAP_AGE[,1,3])/sum(STE_num_dwellings_by_DATUM_obsY[,1]),"\n")}
  if ((z %% 10000)==1) {
    xx <- true_IRSAD_HDIAP_09[,1]/(true_IRSAD_HDIAP_09[,1]+true_IRSAD_HDIAP_09[,2])
    yy <- IRSAD_by_DATUM_09[,1]/(IRSAD_by_DATUM_09[,1]+IRSAD_by_DATUM_09[,2])
    cat(sum(z/10000000),(xx-yy)/yy,"\n")}
  if ((z %% 10000)==1) {
    xx <- true_IRSAD_HDIAP_1014[,1]/(true_IRSAD_HDIAP_1014[,1]+true_IRSAD_HDIAP_1014[,2])
    yy <- IRSAD_by_DATUM_1014[,1]/(IRSAD_by_DATUM_1014[,1]+IRSAD_by_DATUM_1014[,2])
    cat(sum(z/10000000),(xx-yy)/yy,"\n")}
  if ((z %% 10000)==1) {
    xx <- true_IRSAD_HDIAP_1519[,1]/(true_IRSAD_HDIAP_1519[,1]+true_IRSAD_HDIAP_1519[,2])
    yy <- IRSAD_by_DATUM_1519[,1]/(IRSAD_by_DATUM_1519[,1]+IRSAD_by_DATUM_1519[,2])
    cat(sum(z/10000000),(xx-yy)/yy,"\n")}
  if ((z%%1000000)==0) {save(true_SA1_HDIAP_AGE,file=paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",count)))}
} # Z

### Process results

lgas <- vect("/mnt/Z/ewan/GEOMETRIES/POLYGONS/LGA_2021_AUST_GDA2020.shp")
N_LGA <- length(lgas$LGA_NAME21)
lga_names <- lgas$LGA_NAME21

library(rjson)
ndss_0_9_lga <- fromJSON(file='NDSS_0_9_lga.json')
alt_lga_names <- list(N_LGA)
for (i in 1:length(ndss_0_9_lga$MultiPolygons)) {
  alt_lga_names[[i]] <- ndss_0_9_lga$MultiPolygons[[i]]$Name
}
alt_lga_names <- unlist(alt_lga_names)
alt_lga_names <- gsub('\ \\(C\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(M\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(T\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(S\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(R\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(A\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(B\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(DC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(AC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RegC\\)','',alt_lga_names)

ndss_stat_value <- numeric(N_LGA)+NA
for (i in 1:length(ndss_0_9_lga$MultiPolygons)) {
  ndss_stat_value[match(alt_lga_names[[i]],lga_names)] <- ndss_0_9_lga$MultiPolygons[[i]]$StatisticValue
}
genuine_na_LGA_0_9 <- which(is.na(ndss_stat_value))
suppressed_LGA_0_9 <- which(!is.na(ndss_stat_value) & ndss_stat_value==-2)
ndss_stat_value[genuine_na_LGA_0_9] <- NA
ndss_stat_value[suppressed_LGA_0_9] <- NA
ndss_LGA_0_9 <- ndss_stat_value

ndss_10_19_lga <- fromJSON(file='NDSS_10_19_lga.json')
alt_lga_names <- list(N_LGA)
for (i in 1:length(ndss_10_19_lga$MultiPolygons)) {
  alt_lga_names[[i]] <- ndss_10_19_lga$MultiPolygons[[i]]$Name
}
alt_lga_names <- unlist(alt_lga_names)
alt_lga_names <- gsub('\ \\(C\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(M\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(T\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(S\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(R\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(A\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(B\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(DC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(AC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RegC\\)','',alt_lga_names)

ndss_stat_value <- numeric(N_LGA)+NA
for (i in 1:length(ndss_10_19_lga$MultiPolygons)) {
  ndss_stat_value[match(alt_lga_names[[i]],lga_names)] <- ndss_10_19_lga$MultiPolygons[[i]]$StatisticValue
}
genuine_na_LGA_10_19 <- which(is.na(ndss_stat_value))
suppressed_LGA_10_19 <- which(!is.na(ndss_stat_value) & ndss_stat_value==-2)
ndss_stat_value[genuine_na_LGA_10_19] <- NA
ndss_stat_value[suppressed_LGA_10_19] <- NA
ndss_LGA_10_19 <- ndss_stat_value

mb2sa1 <- read_xlsx('MB_2021_AUST.xlsx')
mb2sa1 <- mb2sa1[which(mb2sa1$SA1_CODE_2021 %in% unique(sa1_2021$SA1_CODE21)),]
mbsa1codes <- match(mb2sa1$SA1_CODE_2021,sa1_2021$SA1_CODE21)

mb2lga <- read_xlsx('LGA_2021_AUST.xlsx')
mb2lga <- mb2lga[match(mb2sa1$MB_CODE_2021,mb2lga$MB_CODE_2021),]
mblgacodes <- match(mb2lga$LGA_NAME_2021,lga_names)

sa12lga <- aggregate(rep(1,length(mbsa1codes)),list(mbsa1codes*10^5+mblgacodes),sum)
sa12lga$sa1 <- floor(sa12lga$Group.1/10^5)
sa12lga$lga <- sa12lga$Group.1 %% 10^5
sa12lga <- sa12lga[sort.list(sa12lga$sa1),]

sa12lgax <- sparseMatrix(sa12lga$sa1,sa12lga$lga,x=sa12lga$x,dims = c(length(sa1_2021$SA1_CODE21),N_LGA))
normalisation <- as.numeric(sa12lgax%*%rep(1,N_LGA))
sa12lga <- sparseMatrix(sa12lga$lga,sa12lga$sa1,x=sa12lga$x/normalisation[sa12lga$sa1],dims = c(N_LGA,length(sa1_2021$SA1_CODE21)))

ntot <- floor(z/1000000)

quantilelow <- function(x) {quantile(x,0.025)}
quantilehigh <- function(x) {quantile(x,0.975)}
quantilemedian <- function(x) {quantile(x,0.5)}

sa1_2021_outputs <- sa1_2021

posterior_draws <- list()
posterior_draws_old <- list()
posterior_draws_NDSS_LGA <- list()
posterior_draws_NDSS_LGA_old <- list()
posterior_draws_NDSS_LGAw <- list()
posterior_draws_NDSS_LGA_oldw <- list()
tot_pops_pred <- list()
tot_pops_pred_old <- list()

for (i in 2:ntot) {
  
  load(paste0("/mnt/Z/ewan/DIABETES/xsample_final",sprintf("%03i",i)))
  load(paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",i)))
 
  for (j in 1:10) {
    xout <- ilogit(objx$report(xsample[j,])$joint_field)
    posterior_draws[[length(posterior_draws)+1]] <- (xout[,1])
    posterior_draws_old[[length(posterior_draws_old)+1]] <- (xout[,2])/2+(xout[,3])/2
  
    pd <- rowSums(true_SA1_HDIAP_AGE[,,1])*posterior_draws[[length(posterior_draws)]]
    pd_old <- (rowSums(true_SA1_HDIAP_AGE[,,2])+rowSums(true_SA1_HDIAP_AGE[,,3]))*posterior_draws_old[[length(posterior_draws_old)]]
    posterior_draws_NDSS_LGA[[length(posterior_draws_NDSS_LGA)+1]] <- as.numeric(sa12lga%*%pd)
    posterior_draws_NDSS_LGA_old[[length(posterior_draws_NDSS_LGA_old)+1]] <- as.numeric(sa12lga%*%pd_old)
  }

  posterior_draws_NDSS_LGAw[[length(posterior_draws_NDSS_LGAw)+1]] <- as.numeric(sa12lga%*%true_SA1_HDIAP_AGE[,1,1])
  posterior_draws_NDSS_LGA_oldw[[length(posterior_draws_NDSS_LGA_oldw)+1]] <- as.numeric(sa12lga%*%(true_SA1_HDIAP_AGE[,1,2]+true_SA1_HDIAP_AGE[,1,3]))
  tot_pops_pred[[length(tot_pops_pred)+1]] <- rowSums(true_SA1_HDIAP_AGE[,,1])
  tot_pops_pred_old[[length(tot_pops_pred_old)+1]] <- rowSums((true_SA1_HDIAP_AGE[,,1]+true_SA1_HDIAP_AGE[,,2]))
  
  cat(i,"\n")
}

posterior_draws <- do.call(rbind,posterior_draws)
posterior_draws_old <- do.call(rbind,posterior_draws_old)
posterior_draws_NDSS_LGA <- do.call(rbind,posterior_draws_NDSS_LGA)
posterior_draws_NDSS_LGA_old <- do.call(rbind,posterior_draws_NDSS_LGA_old)
posterior_draws_NDSS_LGAw <- do.call(rbind,posterior_draws_NDSS_LGAw)
posterior_draws_NDSS_LGA_oldw <- do.call(rbind,posterior_draws_NDSS_LGA_oldw)

posterior_risk_median_young <- apply((posterior_draws),2,median)
posterior_risk_median_old <- apply((posterior_draws_old),2,median)
posterior_risk_median_NDSS_LGA <- apply((posterior_draws_NDSS_LGA),2,median)
posterior_risk_median_NDSS_LGA_old <- apply((posterior_draws_NDSS_LGA_old),2,median)
posterior_risk_median_NDSS_LGAw <- apply((posterior_draws_NDSS_LGAw),2,median)
posterior_risk_median_NDSS_LGA_oldw <- apply((posterior_draws_NDSS_LGA_oldw),2,median)
posterior_risk_low_NDSS_LGAw <- apply((posterior_draws_NDSS_LGAw),2,quantilelow)
posterior_risk_low_NDSS_LGA_oldw <- apply((posterior_draws_NDSS_LGA_oldw),2,quantilelow)
posterior_risk_high_NDSS_LGAw <- apply((posterior_draws_NDSS_LGAw),2,quantilehigh)
posterior_risk_high_NDSS_LGA_oldw <- apply((posterior_draws_NDSS_LGA_oldw),2,quantilehigh)

posterior_risk_sd_young <- apply((posterior_draws),2,sd)
posterior_risk_sd_old <- apply((posterior_draws_old),2,sd)

sa1_2021_outputs$RMEAY <- posterior_risk_median_young
sa1_2021_outputs$RMEAO <- posterior_risk_median_old
sa1_2021_outputs$SSMEY <- posterior_risk_sd_young
sa1_2021_outputs$SSMEO <- posterior_risk_sd_old

exceedence_draws <- posterior_draws>median(sa1_2021_outputs$RMEAY)
exceedence <- apply(exceedence_draws,2,mean)
sa1_2021_outputs$XXXY <- exceedence
exceedence_draws <- posterior_draws_old>median(sa1_2021_outputs$RMEAO)
exceedence <- apply(exceedence_draws,2,mean)
sa1_2021_outputs$XXXO <- exceedence

surprised_lga <- numeric(N_LGA)
for (i in 1:N_LGA) {
  xsample <- posterior_draws_NDSS_LGAw[,i]
  xsample[xsample<20] <- 0
  xsample <- round(xsample/10)*10
  xsample[xsample==0] <- -99
  surprised_lga[i] <- as.integer(mean(ndss_LGA_0_9[i]<=xsample) < 0.05 | mean(ndss_LGA_0_9[i]>=xsample) < 0.05)
}
surprised_lga_old <- numeric(N_LGA)
for (i in 1:N_LGA) {
  xsample <- posterior_draws_NDSS_LGA_oldw[,i]
  xsample[xsample<20] <- 0
  xsample <- round(xsample/10)*10
  xsample[xsample==0] <- -99
  surprised_lga_old[i] <- as.integer(mean(ndss_LGA_10_19[i]<=xsample) < 0.05 | mean(ndss_LGA_10_19[i]>=xsample) < 0.05)
}

lga_2021_outputs <- lgas
lga_2021_outputs$RMEAY <- posterior_risk_median_NDSS_LGAw
lga_2021_outputs$RMEAO <- posterior_risk_median_NDSS_LGA_oldw
lga_2021_outputs$LMEAY <- posterior_risk_low_NDSS_LGAw
lga_2021_outputs$LMEAO <- posterior_risk_low_NDSS_LGA_oldw
lga_2021_outputs$UMEAY <- posterior_risk_high_NDSS_LGAw
lga_2021_outputs$UMEAO <- posterior_risk_high_NDSS_LGA_oldw
lga_2021_outputs$OMEAY <- ndss_LGA_0_9
lga_2021_outputs$OMEAO <- ndss_LGA_10_19

lga_2021_outputs$surp <- surprised_lga
lga_2021_outputs$surpo <- surprised_lga_old

writeVector(lga_2021_outputs,"diabetesLGA",filetype="ESRI Shapefile",overwrite=TRUE)

tot_pops_pred <- do.call(rbind,tot_pops_pred)
tot_pops_pred_old <- do.call(rbind,tot_pops_pred_old)
tot_pops_pred_mean <- apply(tot_pops_pred,2,mean)
tot_pops_pred_mean_old <- apply(tot_pops_pred_old,2,mean)
sa1_2021_outputs$NPOPO <- tot_pops_pred_mean_old
sa1_2021_outputs$NPOPY <- tot_pops_pred_mean

writeVector(sa1_2021_outputs,"diabetes",filetype="ESRI Shapefile",overwrite=TRUE)

ses_effect <- list()
ses_effect_old <- list()
ra_effect <- list()
ra_effect_old <- list()
prop <- list()
prop_old <- list()

for (i in 2:ntot) {

    load(paste0("/mnt/Z/ewan/DIABETES/xsample_final",sprintf("%03i",i)))

    for (j in 1:10) {
      ses_effect[[length(ses_effect)+1]] <- xsample[j,colnames(xsample)=="ses_effects_young"]
      ses_effect_old[[length(ses_effect_old)+1]] <- xsample[j,colnames(xsample)=="ses_effects_old"]
      ra_effect[[length(ra_effect)+1]] <- xsample[j,colnames(xsample)=="ra_effects_young"]
      ra_effect_old[[length(ra_effect_old)+1]] <- xsample[j,colnames(xsample)=="ra_effects_old"]
      prop[[length(prop)+1]] <- ilogit(xsample[j,colnames(xsample)=="logit_middle_age_prop"])
    }
}

ses_effect <- do.call(rbind,ses_effect)
ses_effect_old <- do.call(rbind,ses_effect_old)
ra_effect <- do.call(rbind,ra_effect)
ra_effect_old <- do.call(rbind,ra_effect_old)
prop <- do.call(rbind,prop)

ses_effect <- ses_effect-ses_effect[,1]
ses_effect_old <- ses_effect_old-ses_effect_old[,1]
ra_effect <- ra_effect-ra_effect[,1]
ra_effect_old <- ra_effect_old-ra_effect_old[,1]

quantilelow <- function(x) {quantile(x,0.025)}
quantilehigh <- function(x) {quantile(x,0.975)}

ses_effect_median <- apply((ses_effect),2,median)
ses_effect_low <- apply((ses_effect),2,quantilelow)
ses_effect_high <- apply((ses_effect),2,quantilehigh)
cat(sprintf("%3.2f",ses_effect_median),"\n")
cat(sprintf("%3.2f",ses_effect_low),"\n")
cat(sprintf("%3.2f",ses_effect_high),"\n")

ses_effect_old_median <- apply((ses_effect_old),2,median)
ses_effect_old_low <- apply((ses_effect_old),2,quantilelow)
ses_effect_old_high <- apply((ses_effect_old),2,quantilehigh)
cat(sprintf("%3.2f",ses_effect_old_median),"\n")
cat(sprintf("%3.2f",ses_effect_old_low),"\n")
cat(sprintf("%3.2f",ses_effect_old_high),"\n")

ra_effect_median <- apply((ra_effect),2,median)
ra_effect_low <- apply((ra_effect),2,quantilelow)
ra_effect_high <- apply((ra_effect),2,quantilehigh)
cat(sprintf("%3.2f",ra_effect_median),"\n")
cat(sprintf("%3.2f",ra_effect_low),"\n")
cat(sprintf("%3.2f",ra_effect_high),"\n")

ra_effect_old_median <- apply((ra_effect_old),2,median)
ra_effect_old_low <- apply((ra_effect_old),2,quantilelow)
ra_effect_old_high <- apply((ra_effect_old),2,quantilehigh)
cat(sprintf("%3.2f",ra_effect_old_median),"\n")
cat(sprintf("%3.2f",ra_effect_old_low),"\n")
cat(sprintf("%3.2f",ra_effect_old_high),"\n")

prop_median <- apply((prop),2,median)
prop_low <- apply((prop),2,quantilelow)
prop_high <- apply((prop),2,quantilehigh)
cat(sprintf("%3.2f",prop_median),"\n")
cat(sprintf("%3.2f",prop_low),"\n")
cat(sprintf("%3.2f",prop_high),"\n")

young2counts <- read.csv("T1D_T2D counts at -10 yrs at SA2 in 2021.csv")

sa2_code_matches <- match(young2counts$SA2_CODE_2021,unique(sa1_2021$SA2_CODE21))

sa2_tots <- young2counts$T1D_cases+young2counts$T2D_cases

sa2sa3 <- match(sa1_2021$SA3_CODE21[!duplicated(sa1_2021$SA2_CODE21)],unique(sa1_2021$SA3_CODE21))
xx <- aggregate(sa2_tots,list(sa2sa3[sa2_code_matches]),sum)
sa3_tots <- xx$x

tot_sa2s_pred <- list()
for (i in 2:ntot) {
  load(paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",i)))
  tot_sa2s_pred[[length(tot_sa2s_pred)+1]] <- aggregate(true_SA1_HDIAP_AGE[,1,1],list(sa1sa2codes),sum)[,-1][sa2_code_matches]
}
tot_sa2s_pred <- do.call(rbind,tot_sa2s_pred)

tot_sa2s_pred_lower <- apply(tot_sa2s_pred,2,quantilelow)
tot_sa2s_pred_upper <- apply(tot_sa2s_pred,2,quantilehigh)
tot_sa2s_pred_median <- apply(tot_sa2s_pred,2,quantilemedian)
tot_sa2s_pred_mean <- apply(tot_sa2s_pred,2,mean)
tot_sa2s_pred_sd <- apply(tot_sa2s_pred,2,sd)

sa2_tots_jittered <- sa2_tots+runif(length(sa2_tots),-0.25,0.25)
tot_sa2s_pred_median_jittered <- tot_sa2s_pred_median+runif(length(sa2_tots),-0.1,0.1)

tot_sa3s_pred <- list()
for (i in 2:ntot) {
  load(paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",i)))
  tot_sa3s_pred[[length(tot_sa3s_pred)+1]] <- aggregate(true_SA1_HDIAP_AGE[,1,1],list(sa1sa3codes),sum)[,-1][xx$Group.1]
}
tot_sa3s_pred <- do.call(rbind,tot_sa3s_pred)

tot_sa3s_pred_lower <- apply(tot_sa3s_pred,2,quantilelow)
tot_sa3s_pred_upper <- apply(tot_sa3s_pred,2,quantilehigh)
tot_sa3s_pred_median <- apply(tot_sa3s_pred,2,quantilemedian)
tot_sa3s_pred_mean <- apply(tot_sa3s_pred,2,mean)
tot_sa3s_pred_sd <- apply(tot_sa3s_pred,2,sd)

sa3_tots_jittered <- sa3_tots+runif(length(sa3_tots),-0.25,0.25)
tot_sa3s_pred_median_jittered <- tot_sa3s_pred_mean+runif(length(sa3_tots),-0.1,0.1)

ABS_sa3s_pred_median_jittered <- SA3_num_dwellings_by_DATUM_obs[xx$Group.1,1]+runif(length(sa3_tots),-0.1,0.1)

cor(SA3_num_dwellings_by_DATUM_obs[xx$Group.1,1],sa3_tots)
cor(tot_sa3s_pred_mean,sa3_tots)


pdf("figa.pdf",width=5,height=5)
layout(1)
par(mai=c(1,1,0.05,0.05))

plot(tot_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.4,col="white",xlab="",ylab="",xlim=c(-0.5,23.5),ylim=c(-0.5,23.5),xaxt='n',yaxt='n',xaxs='i',yaxs='i')
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(i-0.5,-0.5,i+0.5,23.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(-0.5,i-0.5,23.5,i+0.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
lines(c(-0.5,23.5),c(-0.5,23.5),col="grey50")
points(tot_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.7,col=hsv(0.82,v=1,alpha=0.25))
for (i in 1:length(tot_sa3s_pred_median_jittered)) {
  lines(c(tot_sa3s_pred_lower[i]+tot_sa3s_pred_median_jittered[i]-tot_sa3s_pred_median[i],tot_sa3s_pred_upper[i]+tot_sa3s_pred_median_jittered[i]-tot_sa3s_pred_median[i]),c(sa3_tots_jittered[i],sa3_tots_jittered[i]),col=hsv(0.82,v=1,alpha=0.25))
}
dev.off()

pdf("figb.pdf",width=5,height=5)
layout(1)
par(mai=c(1,1,0.05,0.05))

plot(ABS_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.4,col="white",xlab="",ylab="",xlim=c(-0.5,23.5),ylim=c(-0.5,23.5),xaxt='n',yaxt='n',xaxs='i',yaxs='i')
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(i-0.5,-0.5,i+0.5,23.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(-0.5,i-0.5,23.5,i+0.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
lines(c(-0.5,23.5),c(-0.5,23.5),col="grey50")
points(ABS_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.7,col=hsv(0.82,v=1,alpha=0.25))
dev.off()

### Model for the count of diabetes per SA1 unit over all Australia

setwd("~/VIRTUAL_WA/DIABETES/")

library(terra)
library(MASS)
library(Matrix)
library(truncnorm)
library(TMB)
library(sparseMVN)
library(readxl)
library(spdep)

### Read in and prepare data

sa1_2021 <- vect("/mnt/Z/ewan/DIABETES/sa1_2021_wcovsx.shp")
sa1_2021 <- sa1_2021[sa1_2021$STE_NAME21!="Other Territories",]

sa1_coords <- crds(centroids(sa1_2021))

# Read in RA code table and impute by nearest neighbours etc
ra_table <- read.csv("SA1_RA_correspondence.csv",skip=11,header = FALSE)

lookups <- match(sa1_2021$SA1_CODE21,ra_table$V1)
racodes <- apply(cbind(ra_table$V2,ra_table$V3,ra_table$V4,ra_table$V5,ra_table$V6,ra_table$V7,ra_table$V8)[lookups,],1,which.max)
racodes[which(apply(cbind(ra_table$V2,ra_table$V3,ra_table$V4,ra_table$V5,ra_table$V6,ra_table$V7,ra_table$V8)[lookups,],1,max)==0)] <- NA

nb <- poly2nb(sf::st_as_sf(sa1_2021))

to.infill <- which(is.na(racodes))
for (i in 1:length(to.infill)) {
  av.code <- round(mean(racodes[nb[[to.infill[i]]]],na.rm=TRUE))
  racodes[to.infill[i]] <- av.code
}

mean.with.na <- function(x) {mean(x,na.rm=TRUE)}
mean.ras <- aggregate(racodes,list(match(sa1_2021$SA2_CODE21,unique(sa1_2021$SA2_CODE21))),mean.with.na)
mean.ras <- mean.ras[match(sa1_2021$SA2_CODE21,unique(sa1_2021$SA2_CODE21)),2]

to.infill <- which(is.na(racodes))
racodes[to.infill] <- round(mean.ras[to.infill])

# Read in IRSAD code table and impute by nearest neighbours etc
irsad_table <- read.csv("SA1_IRSAD_lookup.csv.csv",skip=11,header = FALSE)

lookups <- match(sa1_2021$SA1_CODE21,irsad_table$V1)
irsad_codes <- apply(cbind(irsad_table$V2,irsad_table$V3,irsad_table$V4,irsad_table$V5,irsad_table$V6,irsad_table$V7,irsad_table$V8,irsad_table$V9,irsad_table$V10,irsad_table$V11,irsad_table$V12,irsad_table$V13)[lookups,],1,which.max)

### Read in CENSUS DATA : SA1 x HDIAP

filename_SA1 <- "SA1_HDIAP_final_0_9.csv"
filename_SA2 <- "SA2_HDIAP_final_0_9.csv"
filename_SA3 <- "SA3_HDIAP_final_0_9.csv"
filename_SA4 <- "SA4_HDIAP_final_0_9.csv"
filename_STE <- "STE_HDIAP_final_0_9.csv"

sa2_corespondences <- read.csv("CG_SA2_2016_SA2_2021.csv")
sa2_2016_young <- read.csv("T1D counts at SA2 in 2021.csv")

# SA1 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=12,header = FALSE)
SA1_num_dwellings_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_num_dwellings_by_DATUM <- cbind(SA1_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA1_num_dwellings_by_DATUM) <- xcolnames
SA1_num_dwellings_by_DATUM <- as.data.frame(SA1_num_dwellings_by_DATUM)
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[-which(is.na(SA1_num_dwellings_by_DATUM$SA1)),]
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[match(as.numeric(unique(sa1_2021$SA1_CODE21)),SA1_num_dwellings_by_DATUM$SA1),]
SA1_num_dwellings_by_DATUM <- as.matrix(SA1_num_dwellings_by_DATUM)

# SA2 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=12,header = FALSE)
SA2_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA2_num_dwellings_by_DATUM <- cbind(SA2_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(SA2_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA2_num_dwellings_by_DATUM) <- xcolnames
SA2_num_dwellings_by_DATUM <- as.data.frame(SA2_num_dwellings_by_DATUM)
SA2_num_dwellings_by_DATUM <- SA2_num_dwellings_by_DATUM[match((unique(sa1_2021$SA2_NAME21)),SA2_num_dwellings_by_DATUM$SA2),]
SA2_num_dwellings_by_DATUM[,1] <- 1:length(SA2_num_dwellings_by_DATUM[,1])
SA2_num_dwellings_by_DATUM[,2] <- as.numeric(SA2_num_dwellings_by_DATUM[,2])
SA2_num_dwellings_by_DATUM[,3] <- as.numeric(SA2_num_dwellings_by_DATUM[,3])
SA2_num_dwellings_by_DATUM[,4] <- as.numeric(SA2_num_dwellings_by_DATUM[,4])
SA2_num_dwellings_by_DATUM[,5] <- as.numeric(SA2_num_dwellings_by_DATUM[,5])
SA2_num_dwellings_by_DATUM <- as.matrix(SA2_num_dwellings_by_DATUM)

# SA3 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=12,header = FALSE)
SA3_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA3_num_dwellings_by_DATUM <- cbind(SA3_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA3",ydata[-1][1:(length(SA3_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA3_num_dwellings_by_DATUM) <- xcolnames
SA3_num_dwellings_by_DATUM <- as.data.frame(SA3_num_dwellings_by_DATUM)
SA3_num_dwellings_by_DATUM <- SA3_num_dwellings_by_DATUM[match((unique(sa1_2021$SA3_NAME21)),SA3_num_dwellings_by_DATUM$SA3),]
SA3_num_dwellings_by_DATUM[,1] <- 1:length(SA3_num_dwellings_by_DATUM[,1])
SA3_num_dwellings_by_DATUM[,2] <- as.numeric(SA3_num_dwellings_by_DATUM[,2])
SA3_num_dwellings_by_DATUM[,3] <- as.numeric(SA3_num_dwellings_by_DATUM[,3])
SA3_num_dwellings_by_DATUM[,4] <- as.numeric(SA3_num_dwellings_by_DATUM[,4])
SA3_num_dwellings_by_DATUM[,5] <- as.numeric(SA3_num_dwellings_by_DATUM[,5])
SA3_num_dwellings_by_DATUM <- as.matrix(SA3_num_dwellings_by_DATUM)

# SA4 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=12,header = FALSE)
SA4_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA4_num_dwellings_by_DATUM <- cbind(SA4_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA4",ydata[-1][1:(length(SA4_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA4_num_dwellings_by_DATUM) <- xcolnames
SA4_num_dwellings_by_DATUM <- as.data.frame(SA4_num_dwellings_by_DATUM)
SA4_num_dwellings_by_DATUM <- SA4_num_dwellings_by_DATUM[match((unique(sa1_2021$SA4_NAME21)),SA4_num_dwellings_by_DATUM$SA4),]
SA4_num_dwellings_by_DATUM[,1] <- 1:length(SA4_num_dwellings_by_DATUM[,1])
SA4_num_dwellings_by_DATUM[,2] <- as.numeric(SA4_num_dwellings_by_DATUM[,2])
SA4_num_dwellings_by_DATUM[,3] <- as.numeric(SA4_num_dwellings_by_DATUM[,3])
SA4_num_dwellings_by_DATUM[,4] <- as.numeric(SA4_num_dwellings_by_DATUM[,4])
SA4_num_dwellings_by_DATUM[,5] <- as.numeric(SA4_num_dwellings_by_DATUM[,5])
SA4_num_dwellings_by_DATUM <- as.matrix(SA4_num_dwellings_by_DATUM)

# STE level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=12,header = FALSE)
STE_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("STE_num_dwellings_by_DATUM <- cbind(STE_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("STE",ydata[-1][1:(length(STE_num_dwellings_by_DATUM[1,])-1)]))
colnames(STE_num_dwellings_by_DATUM) <- xcolnames
STE_num_dwellings_by_DATUM <- as.data.frame(STE_num_dwellings_by_DATUM)
STE_num_dwellings_by_DATUM <- STE_num_dwellings_by_DATUM[match((unique(sa1_2021$STE_NAME21)),STE_num_dwellings_by_DATUM$STE),]
STE_num_dwellings_by_DATUM[,1] <- 1:length(STE_num_dwellings_by_DATUM[,1])
STE_num_dwellings_by_DATUM[,2] <- as.numeric(STE_num_dwellings_by_DATUM[,2])
STE_num_dwellings_by_DATUM[,3] <- as.numeric(STE_num_dwellings_by_DATUM[,3])
STE_num_dwellings_by_DATUM[,4] <- as.numeric(STE_num_dwellings_by_DATUM[,4])
STE_num_dwellings_by_DATUM[,5] <- as.numeric(STE_num_dwellings_by_DATUM[,5])
STE_num_dwellings_by_DATUM <- as.matrix(STE_num_dwellings_by_DATUM)

sa1sa2codes <- match(sa1_2021$SA2_NAME21,unique(sa1_2021$SA2_NAME21))
sa1sa3codes <- match(sa1_2021$SA3_NAME21,unique(sa1_2021$SA3_NAME21))
sa1sa4codes <- match(sa1_2021$SA4_NAME21,unique(sa1_2021$SA4_NAME21))
sa1stecodes <- match(sa1_2021$STE_NAME21,unique(sa1_2021$STE_NAME21))

M <- dim(SA1_num_dwellings_by_DATUM)[2]-2

SA1_num_dwellings_by_DATUM_obs <- SA1_num_dwellings_by_DATUM[,-c(1,M+2)]
SA1_num_dwellings_obs <- SA1_num_dwellings_by_DATUM[,M+2]
SA2_num_dwellings_by_DATUM_obs <- SA2_num_dwellings_by_DATUM[,-c(1,M+2)]
SA2_num_dwellings_obs <- SA2_num_dwellings_by_DATUM[,M+2]
SA3_num_dwellings_by_DATUM_obs <- SA3_num_dwellings_by_DATUM[,-c(1,M+2)]
SA3_num_dwellings_obs <- SA3_num_dwellings_by_DATUM[,M+2]
SA4_num_dwellings_by_DATUM_obs <- SA4_num_dwellings_by_DATUM[,-c(1,M+2)]
SA4_num_dwellings_obs <- SA4_num_dwellings_by_DATUM[,M+2]
STE_num_dwellings_by_DATUM_obs <- STE_num_dwellings_by_DATUM[,-c(1,M+2)]
STE_num_dwellings_obs <- STE_num_dwellings_by_DATUM[,M+2]

# Load 10-14 EXTRA

filename_SA1 <- "SA1_HDIAP_final_10_14.csv"
filename_SA2 <- "SA2_HDIAP_final_10_14.csv"
filename_SA3 <- "SA3_HDIAP_final_10_14.csv"
filename_SA4 <- "SA4_HDIAP_final_10_14.csv"
filename_STE <- "STE_HDIAP_final_10_14.csv"

# SA1 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=12,header = FALSE)
SA1_num_dwellings_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_num_dwellings_by_DATUM <- cbind(SA1_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA1_num_dwellings_by_DATUM) <- xcolnames
SA1_num_dwellings_by_DATUM <- as.data.frame(SA1_num_dwellings_by_DATUM)
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[-which(is.na(SA1_num_dwellings_by_DATUM$SA1)),]
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[match(as.numeric(unique(sa1_2021$SA1_CODE21)),SA1_num_dwellings_by_DATUM$SA1),]
SA1_num_dwellings_by_DATUM <- as.matrix(SA1_num_dwellings_by_DATUM)

# SA2 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=12,header = FALSE)
SA2_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA2_num_dwellings_by_DATUM <- cbind(SA2_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(SA2_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA2_num_dwellings_by_DATUM) <- xcolnames
SA2_num_dwellings_by_DATUM <- as.data.frame(SA2_num_dwellings_by_DATUM)
SA2_num_dwellings_by_DATUM <- SA2_num_dwellings_by_DATUM[match((unique(sa1_2021$SA2_NAME21)),SA2_num_dwellings_by_DATUM$SA2),]
SA2_num_dwellings_by_DATUM[,1] <- 1:length(SA2_num_dwellings_by_DATUM[,1])
SA2_num_dwellings_by_DATUM[,2] <- as.numeric(SA2_num_dwellings_by_DATUM[,2])
SA2_num_dwellings_by_DATUM[,3] <- as.numeric(SA2_num_dwellings_by_DATUM[,3])
SA2_num_dwellings_by_DATUM[,4] <- as.numeric(SA2_num_dwellings_by_DATUM[,4])
SA2_num_dwellings_by_DATUM[,5] <- as.numeric(SA2_num_dwellings_by_DATUM[,5])
SA2_num_dwellings_by_DATUM <- as.matrix(SA2_num_dwellings_by_DATUM)

# SA3 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=12,header = FALSE)
SA3_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA3_num_dwellings_by_DATUM <- cbind(SA3_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA3",ydata[-1][1:(length(SA3_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA3_num_dwellings_by_DATUM) <- xcolnames
SA3_num_dwellings_by_DATUM <- as.data.frame(SA3_num_dwellings_by_DATUM)
SA3_num_dwellings_by_DATUM <- SA3_num_dwellings_by_DATUM[match((unique(sa1_2021$SA3_NAME21)),SA3_num_dwellings_by_DATUM$SA3),]
SA3_num_dwellings_by_DATUM[,1] <- 1:length(SA3_num_dwellings_by_DATUM[,1])
SA3_num_dwellings_by_DATUM[,2] <- as.numeric(SA3_num_dwellings_by_DATUM[,2])
SA3_num_dwellings_by_DATUM[,3] <- as.numeric(SA3_num_dwellings_by_DATUM[,3])
SA3_num_dwellings_by_DATUM[,4] <- as.numeric(SA3_num_dwellings_by_DATUM[,4])
SA3_num_dwellings_by_DATUM[,5] <- as.numeric(SA3_num_dwellings_by_DATUM[,5])
SA3_num_dwellings_by_DATUM <- as.matrix(SA3_num_dwellings_by_DATUM)

# SA4 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=12,header = FALSE)
SA4_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA4_num_dwellings_by_DATUM <- cbind(SA4_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA4",ydata[-1][1:(length(SA4_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA4_num_dwellings_by_DATUM) <- xcolnames
SA4_num_dwellings_by_DATUM <- as.data.frame(SA4_num_dwellings_by_DATUM)
SA4_num_dwellings_by_DATUM <- SA4_num_dwellings_by_DATUM[match((unique(sa1_2021$SA4_NAME21)),SA4_num_dwellings_by_DATUM$SA4),]
SA4_num_dwellings_by_DATUM[,1] <- 1:length(SA4_num_dwellings_by_DATUM[,1])
SA4_num_dwellings_by_DATUM[,2] <- as.numeric(SA4_num_dwellings_by_DATUM[,2])
SA4_num_dwellings_by_DATUM[,3] <- as.numeric(SA4_num_dwellings_by_DATUM[,3])
SA4_num_dwellings_by_DATUM[,4] <- as.numeric(SA4_num_dwellings_by_DATUM[,4])
SA4_num_dwellings_by_DATUM[,5] <- as.numeric(SA4_num_dwellings_by_DATUM[,5])
SA4_num_dwellings_by_DATUM <- as.matrix(SA4_num_dwellings_by_DATUM)

# STE level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=12,header = FALSE)
STE_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("STE_num_dwellings_by_DATUM <- cbind(STE_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("STE",ydata[-1][1:(length(STE_num_dwellings_by_DATUM[1,])-1)]))
colnames(STE_num_dwellings_by_DATUM) <- xcolnames
STE_num_dwellings_by_DATUM <- as.data.frame(STE_num_dwellings_by_DATUM)
STE_num_dwellings_by_DATUM <- STE_num_dwellings_by_DATUM[match((unique(sa1_2021$STE_NAME21)),STE_num_dwellings_by_DATUM$STE),]
STE_num_dwellings_by_DATUM[,1] <- 1:length(STE_num_dwellings_by_DATUM[,1])
STE_num_dwellings_by_DATUM[,2] <- as.numeric(STE_num_dwellings_by_DATUM[,2])
STE_num_dwellings_by_DATUM[,3] <- as.numeric(STE_num_dwellings_by_DATUM[,3])
STE_num_dwellings_by_DATUM[,4] <- as.numeric(STE_num_dwellings_by_DATUM[,4])
STE_num_dwellings_by_DATUM[,5] <- as.numeric(STE_num_dwellings_by_DATUM[,5])
STE_num_dwellings_by_DATUM <- as.matrix(STE_num_dwellings_by_DATUM)

SA1_num_dwellings_by_DATUM_obsX <- SA1_num_dwellings_by_DATUM[,-c(1,M+2)]
SA1_num_dwellings_obsX <- SA1_num_dwellings_by_DATUM[,M+2]
SA2_num_dwellings_by_DATUM_obsX <- SA2_num_dwellings_by_DATUM[,-c(1,M+2)]
SA2_num_dwellings_obsX <- SA2_num_dwellings_by_DATUM[,M+2]
SA3_num_dwellings_by_DATUM_obsX <- SA3_num_dwellings_by_DATUM[,-c(1,M+2)]
SA3_num_dwellings_obsX <- SA3_num_dwellings_by_DATUM[,M+2]
SA4_num_dwellings_by_DATUM_obsX <- SA4_num_dwellings_by_DATUM[,-c(1,M+2)]
SA4_num_dwellings_obsX <- SA4_num_dwellings_by_DATUM[,M+2]
STE_num_dwellings_by_DATUM_obsX <- STE_num_dwellings_by_DATUM[,-c(1,M+2)]
STE_num_dwellings_obsX <- STE_num_dwellings_by_DATUM[,M+2]

# Load 15-19 EXTRA

filename_SA1 <- "SA1_HDIAP_final_15_19.csv"
filename_SA2 <- "SA2_HDIAP_final_15_19.csv"
filename_SA3 <- "SA3_HDIAP_final_15_19.csv"
filename_SA4 <- "SA4_HDIAP_final_15_19.csv"
filename_STE <- "STE_HDIAP_final_15_19.csv"

# SA1 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=12,header = FALSE)
SA1_num_dwellings_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_num_dwellings_by_DATUM <- cbind(SA1_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA1),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA1_num_dwellings_by_DATUM) <- xcolnames
SA1_num_dwellings_by_DATUM <- as.data.frame(SA1_num_dwellings_by_DATUM)
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[-which(is.na(SA1_num_dwellings_by_DATUM$SA1)),]
SA1_num_dwellings_by_DATUM <- SA1_num_dwellings_by_DATUM[match(as.numeric(unique(sa1_2021$SA1_CODE21)),SA1_num_dwellings_by_DATUM$SA1),]
SA1_num_dwellings_by_DATUM <- as.matrix(SA1_num_dwellings_by_DATUM)

# SA2 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=12,header = FALSE)
SA2_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA2_num_dwellings_by_DATUM <- cbind(SA2_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA2),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(SA2_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA2_num_dwellings_by_DATUM) <- xcolnames
SA2_num_dwellings_by_DATUM <- as.data.frame(SA2_num_dwellings_by_DATUM)
SA2_num_dwellings_by_DATUM <- SA2_num_dwellings_by_DATUM[match((unique(sa1_2021$SA2_NAME21)),SA2_num_dwellings_by_DATUM$SA2),]
SA2_num_dwellings_by_DATUM[,1] <- 1:length(SA2_num_dwellings_by_DATUM[,1])
SA2_num_dwellings_by_DATUM[,2] <- as.numeric(SA2_num_dwellings_by_DATUM[,2])
SA2_num_dwellings_by_DATUM[,3] <- as.numeric(SA2_num_dwellings_by_DATUM[,3])
SA2_num_dwellings_by_DATUM[,4] <- as.numeric(SA2_num_dwellings_by_DATUM[,4])
SA2_num_dwellings_by_DATUM[,5] <- as.numeric(SA2_num_dwellings_by_DATUM[,5])
SA2_num_dwellings_by_DATUM <- as.matrix(SA2_num_dwellings_by_DATUM)

# SA3 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=12,header = FALSE)
SA3_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA3_num_dwellings_by_DATUM <- cbind(SA3_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA3),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA3",ydata[-1][1:(length(SA3_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA3_num_dwellings_by_DATUM) <- xcolnames
SA3_num_dwellings_by_DATUM <- as.data.frame(SA3_num_dwellings_by_DATUM)
SA3_num_dwellings_by_DATUM <- SA3_num_dwellings_by_DATUM[match((unique(sa1_2021$SA3_NAME21)),SA3_num_dwellings_by_DATUM$SA3),]
SA3_num_dwellings_by_DATUM[,1] <- 1:length(SA3_num_dwellings_by_DATUM[,1])
SA3_num_dwellings_by_DATUM[,2] <- as.numeric(SA3_num_dwellings_by_DATUM[,2])
SA3_num_dwellings_by_DATUM[,3] <- as.numeric(SA3_num_dwellings_by_DATUM[,3])
SA3_num_dwellings_by_DATUM[,4] <- as.numeric(SA3_num_dwellings_by_DATUM[,4])
SA3_num_dwellings_by_DATUM[,5] <- as.numeric(SA3_num_dwellings_by_DATUM[,5])
SA3_num_dwellings_by_DATUM <- as.matrix(SA3_num_dwellings_by_DATUM)

# SA4 level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=12,header = FALSE)
SA4_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA4_num_dwellings_by_DATUM <- cbind(SA4_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_SA4),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA4",ydata[-1][1:(length(SA4_num_dwellings_by_DATUM[1,])-1)]))
colnames(SA4_num_dwellings_by_DATUM) <- xcolnames
SA4_num_dwellings_by_DATUM <- as.data.frame(SA4_num_dwellings_by_DATUM)
SA4_num_dwellings_by_DATUM <- SA4_num_dwellings_by_DATUM[match((unique(sa1_2021$SA4_NAME21)),SA4_num_dwellings_by_DATUM$SA4),]
SA4_num_dwellings_by_DATUM[,1] <- 1:length(SA4_num_dwellings_by_DATUM[,1])
SA4_num_dwellings_by_DATUM[,2] <- as.numeric(SA4_num_dwellings_by_DATUM[,2])
SA4_num_dwellings_by_DATUM[,3] <- as.numeric(SA4_num_dwellings_by_DATUM[,3])
SA4_num_dwellings_by_DATUM[,4] <- as.numeric(SA4_num_dwellings_by_DATUM[,4])
SA4_num_dwellings_by_DATUM[,5] <- as.numeric(SA4_num_dwellings_by_DATUM[,5])
SA4_num_dwellings_by_DATUM <- as.matrix(SA4_num_dwellings_by_DATUM)

# STE level Tablebuilder Outputs
xdata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=12,header = FALSE)
STE_num_dwellings_by_DATUM <- xdata$V1
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("STE_num_dwellings_by_DATUM <- cbind(STE_num_dwellings_by_DATUM,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(paste0("~/VIRTUAL_WA/DIABETES/",filename_STE),skip=10,header = FALSE,nrows = 1)
xcolnames <- unlist(c("STE",ydata[-1][1:(length(STE_num_dwellings_by_DATUM[1,])-1)]))
colnames(STE_num_dwellings_by_DATUM) <- xcolnames
STE_num_dwellings_by_DATUM <- as.data.frame(STE_num_dwellings_by_DATUM)
STE_num_dwellings_by_DATUM <- STE_num_dwellings_by_DATUM[match((unique(sa1_2021$STE_NAME21)),STE_num_dwellings_by_DATUM$STE),]
STE_num_dwellings_by_DATUM[,1] <- 1:length(STE_num_dwellings_by_DATUM[,1])
STE_num_dwellings_by_DATUM[,2] <- as.numeric(STE_num_dwellings_by_DATUM[,2])
STE_num_dwellings_by_DATUM[,3] <- as.numeric(STE_num_dwellings_by_DATUM[,3])
STE_num_dwellings_by_DATUM[,4] <- as.numeric(STE_num_dwellings_by_DATUM[,4])
STE_num_dwellings_by_DATUM[,5] <- as.numeric(STE_num_dwellings_by_DATUM[,5])
STE_num_dwellings_by_DATUM <- as.matrix(STE_num_dwellings_by_DATUM)

SA1_num_dwellings_by_DATUM_obsY <- SA1_num_dwellings_by_DATUM[,-c(1,M+2)]
SA1_num_dwellings_obsY <- SA1_num_dwellings_by_DATUM[,M+2]
SA2_num_dwellings_by_DATUM_obsY <- SA2_num_dwellings_by_DATUM[,-c(1,M+2)]
SA2_num_dwellings_obsY <- SA2_num_dwellings_by_DATUM[,M+2]
SA3_num_dwellings_by_DATUM_obsY <- SA3_num_dwellings_by_DATUM[,-c(1,M+2)]
SA3_num_dwellings_obsY <- SA3_num_dwellings_by_DATUM[,M+2]
SA4_num_dwellings_by_DATUM_obsY <- SA4_num_dwellings_by_DATUM[,-c(1,M+2)]
SA4_num_dwellings_obsY <- SA4_num_dwellings_by_DATUM[,M+2]
STE_num_dwellings_by_DATUM_obsY <- STE_num_dwellings_by_DATUM[,-c(1,M+2)]
STE_num_dwellings_obsY <- STE_num_dwellings_by_DATUM[,M+2]

# IRSAD 0-9
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/IRSAD_HDIAP_09.csv",skip=12,header = FALSE)[1:11,]
IRSAD_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("IRSAD_by_DATUM <- cbind(IRSAD_by_DATUM,as.numeric(xdata$V",i,"))")))
}
IRSAD_by_DATUM <- as.matrix(IRSAD_by_DATUM[,2:4])
IRSAD_by_DATUM_09 <- IRSAD_by_DATUM

# IRSAD 10-14
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/IRSAD_HDIAP_1014.csv",skip=12,header = FALSE)[1:11,]
IRSAD_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("IRSAD_by_DATUM <- cbind(IRSAD_by_DATUM,as.numeric(xdata$V",i,"))")))
}
IRSAD_by_DATUM <- as.matrix(IRSAD_by_DATUM[,2:4])
IRSAD_by_DATUM_1014 <- IRSAD_by_DATUM

# IRSAD 15-19
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/IRSAD_HDIAP_1519.csv",skip=12,header = FALSE)[1:11,]
IRSAD_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("IRSAD_by_DATUM <- cbind(IRSAD_by_DATUM,as.numeric(xdata$V",i,"))")))
}
IRSAD_by_DATUM <- as.matrix(IRSAD_by_DATUM[,2:4])
IRSAD_by_DATUM_1519 <- IRSAD_by_DATUM

# RA 0-9
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/RA_HDIAP_09.csv",skip=12,header = FALSE)[1:5,]
RA_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("RA_by_DATUM <- cbind(RA_by_DATUM,as.numeric(xdata$V",i,"))")))
}
RA_by_DATUM <- as.matrix(RA_by_DATUM[,2:4])
RA_by_DATUM_09 <- RA_by_DATUM

# RA 10-14
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/RA_HDIAP_1014.csv",skip=12,header = FALSE)[1:5,]
RA_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("RA_by_DATUM <- cbind(RA_by_DATUM,as.numeric(xdata$V",i,"))")))
}
RA_by_DATUM <- as.matrix(RA_by_DATUM[,2:4])
RA_by_DATUM_1014 <- RA_by_DATUM

# RA 15-19
xdata <- read.csv("~/VIRTUAL_WA/DIABETES/RA_HDIAP_1519.csv",skip=12,header = FALSE)[1:5,]
RA_by_DATUM <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("RA_by_DATUM <- cbind(RA_by_DATUM,as.numeric(xdata$V",i,"))")))
}
RA_by_DATUM <- as.matrix(RA_by_DATUM[,2:4])
RA_by_DATUM_1519 <- RA_by_DATUM

### Construct geospatial mesh model

library(INLA)

sydney.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Sydney",],cutoff = 0.015,max.n=300)
melbourne.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Melbourne",],cutoff = 0.015,max.n=300)
brisbane.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Brisbane",],cutoff = 0.025,max.n=300)
adelaide.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Adelaide",],cutoff = 0.015,max.n=300)
perth.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Perth",],cutoff = 0.015,max.n=300)
darwin.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Darwin",],cutoff = 0.01,max.n=300)
hobart.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Greater Hobart",],cutoff = 0.01,max.n=300)
canberra.mesh <-  inla.mesh.2d(sa1_coords[sa1_2021$GCC_NAME21=="Australian Capital Territory",],cutoff = 0.01,max.n=300)
remainder.mesh <- inla.mesh.2d(sa1_coords[!(sa1_2021$GCC_NAME21 %in% c("Greater Sydney","Greater Melbourne","Greater Brisbane", "Greater Adelaide", "Greater Darwin", "Greater Hobart","Greater Perth","Australian Capital Territory")),],cutoff = 0.05,max.n=300)
xcoords <- rbind(sydney.mesh$loc,
                 melbourne.mesh$loc,
                 brisbane.mesh$loc,
                 adelaide.mesh$loc,
                 perth.mesh$loc,
                 darwin.mesh$loc,
                 hobart.mesh$loc,
                 canberra.mesh$loc,
                 remainder.mesh$loc)
aus_mesh_fine <- inla.mesh.2d(xcoords,cutoff = 0.01,max.n=1000)
aus_spde_fine <- (inla.spde2.matern(aus_mesh_fine,alpha=0.5)$param.inla)[c("M0","M1","M2")]
aus_A_fine <- inla.mesh.projector(aus_mesh_fine,sa1_coords)$proj$A

### Helper functions

library(extraDistr)

error_sd <- 2.0

logdiffexp <- function(y,x) {
  # x > y
  x+log(1-exp(y-x))
}
logsumexp <- function(y,x) {
  if (length(x)==1) {
    max(c(x,y))+log(1+exp(min(c(x,y))-max(c(x,y))))} else {
      mmax <- apply(cbind(x,y),1,max)
      mmin <- apply(cbind(x,y),1,min)
      mmax+log(1+exp(mmin-mmax))
    }
}

log_likelihood_diff_fn <- function(proposed_value,current_value,observed_value) {
  
  if (proposed_value==0 & observed_value>0) {
    proposed_likelihood <- NA
  } else if (proposed_value==0 & observed_value==0) {
    proposed_likelihood <- 0
  } else if (proposed_value>0 & observed_value==0) {
    proposed_likelihood <- pnorm(2.5,proposed_value,error_sd,log.p=TRUE)
  } else if (abs(proposed_value-observed_value) < 5 & abs(current_value-observed_value) < 5) {
    proposed_likelihood <- logdiffexp(pnorm(observed_value-0.5,proposed_value,error_sd,log.p = TRUE),pnorm(observed_value+0.5,proposed_value,error_sd,log.p = TRUE))
  } else {
    proposed_likelihood <- dnorm(observed_value,proposed_value,error_sd,log = TRUE)
  }
  
  if (current_value==0 & observed_value==0) {
    current_likelihood <- 0
  } else if (current_value>0 & observed_value==0) {
    current_likelihood <- pnorm(2.5,current_value,error_sd,log.p=TRUE)
  } else if (abs(proposed_value-observed_value) < 5 & abs(current_value-observed_value) < 5) {
    current_likelihood <- logdiffexp(pnorm(observed_value-0.5,current_value,error_sd,log.p = TRUE),pnorm(observed_value+0.5,current_value,error_sd,log.p = TRUE))
  } else {
    current_likelihood <- dnorm(observed_value,current_value,error_sd,log = TRUE)
  }
  
  return(proposed_likelihood - current_likelihood)
}

### MCMC Chain

N_SA1 <- length(sa1_2021)
true_SA1_HDIAP_AGE <- array(1,dim=c(N_SA1,3,3))
true_SA1_HDIAP_AGE[,,1] <- SA1_num_dwellings_by_DATUM_obs+1
true_SA1_HDIAP_AGE[,,2] <- SA1_num_dwellings_by_DATUM_obsX+1
true_SA1_HDIAP_AGE[,,3] <- SA1_num_dwellings_by_DATUM_obsY+1

true_SA1_HDIAP_09 <- true_SA1_HDIAP_AGE[,,1]
true_SA2_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1sa2codes),sum)[,-1]
true_SA3_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1sa3codes),sum)[,-1]
true_SA4_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1sa4codes),sum)[,-1]
true_STE_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(sa1stecodes),sum)[,-1]

true_SA1_HDIAP_1014 <- true_SA1_HDIAP_AGE[,,2]
true_SA2_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1sa2codes),sum)[,-1]
true_SA3_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1sa3codes),sum)[,-1]
true_SA4_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1sa4codes),sum)[,-1]
true_STE_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(sa1stecodes),sum)[,-1]

true_SA1_HDIAP_1519 <- true_SA1_HDIAP_AGE[,,3]
true_SA2_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1sa2codes),sum)[,-1]
true_SA3_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1sa3codes),sum)[,-1]
true_SA4_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1sa4codes),sum)[,-1]
true_STE_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(sa1stecodes),sum)[,-1]

true_IRSAD_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(irsad_codes),sum)[,-1]
true_IRSAD_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(irsad_codes),sum)[,-1]
true_IRSAD_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(irsad_codes),sum)[,-1]

true_RA_HDIAP_09 <- aggregate(true_SA1_HDIAP_09,list(racodes),sum)[,-1]
true_RA_HDIAP_1014 <- aggregate(true_SA1_HDIAP_1014,list(racodes),sum)[,-1]
true_RA_HDIAP_1519 <- aggregate(true_SA1_HDIAP_1519,list(racodes),sum)[,-1]

## Fit with spatial model

compile("inla_diabetes_allages.cpp")
dyn.load(dynlib("inla_diabetes_allages"))
ilogit <- function(x) {1/(1+exp(-x))}

count <- 0
prior_draw <- t(matrix(c(sum(STE_num_dwellings_by_DATUM_obs[,1])/(sum(STE_num_dwellings_by_DATUM_obs[,2])+sum(STE_num_dwellings_by_DATUM_obs[,1])),sum(STE_num_dwellings_by_DATUM_obsX[,1])/(sum(STE_num_dwellings_by_DATUM_obsX[,2])+sum(STE_num_dwellings_by_DATUM_obsX[,1])),sum(STE_num_dwellings_by_DATUM_obsY[,1])/(sum(STE_num_dwellings_by_DATUM_obsY[,2])+sum(STE_num_dwellings_by_DATUM_obsY[,1]))),ncol=N_SA1,nrow=3))

accepted_track <- array(0,dim=c(N_SA1,3,3))

for (z in 1:300000000) {
  
  if ((z%%1000000)==0) {
    count <- count + 1
    input.data <- list('NSA1'=N_SA1,
                         'spde'=aus_spde_fine,
                         'A'=aus_A_fine,
                         'npos' = true_SA1_HDIAP_AGE[,1,],
                         'nneg' = true_SA1_HDIAP_AGE[,2,],
                         'sa1_ses' = irsad_codes-1,
                         'sa1_ra' = racodes-1
    )

    if (count==1) {parameters <- list(
      'log_range_young'=0,
      'log_sd_young'=0,
      'log_range_old'=0,
      'log_sd_old'=0,
      'log_sd_ses_young'=0,
      'log_sd_ses_old'=0,
      'log_sd_ra_young'=0,
      'log_sd_ra_old'=0,
      'intercept_young'=-5,
      'intercept_old'=-2,
      'intercept_middle'=0,
      'logit_middle_age_prop'=0,
      'ses_effects_young'=rep(0,11),
      'ra_effects_young'=rep(0,5),
      'ses_effects_old'=rep(0,11),
      'ra_effects_old'=rep(0,5),
      'field_young'=rep(0,aus_mesh_fine$n),
      'field_old'=rep(0,aus_mesh_fine$n)
    )}

    obj <- MakeADFun(input.data,parameters,DLL = "inla_diabetes_allages",random=c('intercept_young','intercept_old','intercept_middle','logit_middle_age_prop','ses_effects_young','ra_effects_young','ses_effects_old','ra_effects_old','field_young','field_old'))
    obj$fn()
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    rep <- sdreport(obj,getJointPrecision = TRUE)

    xsample <- rmvn.sparse(10, c(rep$par.fixed, rep$par.random), Cholesky(rep$jointPrecision), prec = TRUE)
    colnames(xsample) <- names(c(rep$par.fixed, rep$par.random))

    parameters <- list(
      'log_range_young'=xsample[1,colnames(xsample)=="log_range_young"],
      'log_sd_young'=xsample[1,colnames(xsample)=="log_sd_young"],
      'log_range_old'=xsample[1,colnames(xsample)=="log_range_old"],
      'log_sd_old'=xsample[1,colnames(xsample)=="log_sd_old"],
      'log_sd_ses_young'=xsample[1,colnames(xsample)=="log_sd_ses_young"],
      'log_sd_ses_old'=xsample[1,colnames(xsample)=="log_sd_ses_old"],
      'log_sd_ra_young'=xsample[1,colnames(xsample)=="log_sd_ra_young"],
      'log_sd_ra_old'=xsample[1,colnames(xsample)=="log_sd_ra_old"],
      'intercept_young'=0,
      'intercept_old'=0,
      'intercept_middle'=0,
      'logit_middle_age_prop'=0,
      'ses_effects_young'=rep(0,11),
      'ra_effects_young'=rep(0,5),
      'ses_effects_old'=rep(0,11),
      'ra_effects_old'=rep(0,5),
      'field_young'=rep(0,aus_mesh_fine$n),
      'field_old'=rep(0,aus_mesh_fine$n)
    )

    objx <- MakeADFun(input.data,parameters,DLL = "inla_diabetes_allages")

    save(xsample,file=paste0("/mnt/Z/ewan/DIABETES/xsample_final",sprintf("%03i",count)))

    prior_draw <- objx$report(xsample[1,])$joint_field
    prior_draw[,1] <- ilogit(prior_draw[,1])
    prior_draw[,2] <- ilogit(prior_draw[,2])
    prior_draw[,3] <- ilogit(prior_draw[,3])

    xcol <- prior_draw[,1]
    q1 <- quantile(prior_draw[,1],0.1)
    q2 <- quantile(prior_draw[,1],0.9)
    xcol[xcol<q1] <- q1
    xcol[xcol>q2] <- q2
    xcol <- (xcol-min(xcol))/diff(range(xcol))
    plot(sa1_2021,border=NA,col=hsv((1-xcol)*0.666))

    xcol <- prior_draw[,3]
    q1 <- quantile(prior_draw[,3],0.1)
    q2 <- quantile(prior_draw[,3],0.9)
    xcol[xcol<q1] <- q1
    xcol[xcol>q2] <- q2
    xcol <- (xcol-min(xcol))/diff(range(xcol))
    plot(sa1_2021,border=NA,col=hsv((1-xcol)*0.666))
    
    cat(sum(z/10000000),"xx ",xsample[1,colnames(xsample)=="logit_middle_age_prop"],"\n")
    
    gc()
  }

  proposed_location <- cbind(sample(1:N_SA1,1),sample(1:M,1,prob=c(10,2,1)),sample(1:3,1))
  proposed_move <- sample(c(-1,1),1)
  
  current_likelihood_diff <- 0
  
  if ((true_SA1_HDIAP_AGE[proposed_location[,1],proposed_location[,2],proposed_location[,3]] > 0) | (proposed_move==1)) {
    
    # SA1 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA1_HDIAP_09[proposed_location[,1],proposed_location[,2]]
      observed_value <- SA1_num_dwellings_by_DATUM_obs[proposed_location[,1],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA1_HDIAP_1014[proposed_location[,1],proposed_location[,2]]
      observed_value <- SA1_num_dwellings_by_DATUM_obsX[proposed_location[,1],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA1_HDIAP_1519[proposed_location[,1],proposed_location[,2]]
      observed_value <- SA1_num_dwellings_by_DATUM_obsY[proposed_location[,1],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # SA2 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA2_HDIAP_09[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA2_num_dwellings_by_DATUM_obs[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA2_HDIAP_1014[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA2_num_dwellings_by_DATUM_obsX[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA2_HDIAP_1519[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA2_num_dwellings_by_DATUM_obsY[sa1sa2codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # SA3 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA3_HDIAP_09[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA3_num_dwellings_by_DATUM_obs[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA3_HDIAP_1014[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA3_num_dwellings_by_DATUM_obsX[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA3_HDIAP_1519[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA3_num_dwellings_by_DATUM_obsY[sa1sa3codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # SA4 likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_SA4_HDIAP_09[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA4_num_dwellings_by_DATUM_obs[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_SA4_HDIAP_1014[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA4_num_dwellings_by_DATUM_obsX[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_SA4_HDIAP_1519[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- SA4_num_dwellings_by_DATUM_obsY[sa1sa4codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # STE likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_STE_HDIAP_09[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- STE_num_dwellings_by_DATUM_obs[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_STE_HDIAP_1014[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- STE_num_dwellings_by_DATUM_obsX[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_STE_HDIAP_1519[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- STE_num_dwellings_by_DATUM_obsY[sa1stecodes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move
    
    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # IRSAD likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_IRSAD_HDIAP_09[irsad_codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- IRSAD_by_DATUM_09[irsad_codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_IRSAD_HDIAP_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- IRSAD_by_DATUM_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_IRSAD_HDIAP_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- IRSAD_by_DATUM_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move

    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)

    # RA likelihood
    if (proposed_location[,3]==1) {
      current_value <- true_RA_HDIAP_09[racodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- RA_by_DATUM_09[racodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==2) {
      current_value <- true_RA_HDIAP_1014[racodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- RA_by_DATUM_1014[racodes[proposed_location[,1]],proposed_location[,2]]
    }
    if (proposed_location[,3]==3) {
      current_value <- true_RA_HDIAP_1519[racodes[proposed_location[,1]],proposed_location[,2]]
      observed_value <- RA_by_DATUM_1519[racodes[proposed_location[,1]],proposed_location[,2]]
    }
    proposed_value <- current_value+proposed_move

    current_likelihood_diff <- current_likelihood_diff+log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    
    # Regression priors
    
    current_pos <- true_SA1_HDIAP_AGE[proposed_location[1],1,proposed_location[3]]
    current_neg <- true_SA1_HDIAP_AGE[proposed_location[1],2,proposed_location[3]]

    proposed_pos <- current_pos
    proposed_neg <- current_neg
    if (proposed_location[2]==1) {proposed_pos <- proposed_pos+proposed_move}
    if (proposed_location[2]==2) {proposed_neg <- proposed_neg+proposed_move}

    proposed_priors <- current_priors <- 0

    if ((proposed_location[2] %in% c(1,2)) & ((current_pos+current_neg)>0 & (proposed_pos+proposed_neg)>0)) {
      current_priors <- current_priors + dbinom(current_pos,(current_pos+current_neg),prior_draw[proposed_location[1],proposed_location[3]],log=TRUE)
      proposed_priors <- proposed_priors + dbinom(proposed_pos,(proposed_pos+proposed_neg),prior_draw[proposed_location[1],proposed_location[3]],log=TRUE)
    }

    current_likelihood_diff <- current_likelihood_diff + proposed_priors - current_priors
    
    # if accepted
    if (is.na(current_likelihood_diff) || current_likelihood_diff < log(runif(1))) {
      proposed_likelihood <- NA
    } else {
      
      accepted_track[proposed_location[,1],proposed_location[,2],proposed_location[,3]] <-     accepted_track[proposed_location[,1],proposed_location[,2],proposed_location[,3]] + 1
      
      true_SA1_HDIAP_AGE[proposed_location[,1],proposed_location[,2],proposed_location[,3]] <- true_SA1_HDIAP_AGE[proposed_location[,1],proposed_location[,2],proposed_location[,3]] + proposed_move
      
      if (proposed_location[,3]==1) {
        true_SA1_HDIAP_09[proposed_location[,1],proposed_location[,2]] <- true_SA1_HDIAP_09[proposed_location[,1],proposed_location[,2]] + proposed_move
        true_SA2_HDIAP_09[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] <- true_SA2_HDIAP_09[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA3_HDIAP_09[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] <- true_SA3_HDIAP_09[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA4_HDIAP_09[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] <- true_SA4_HDIAP_09[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_STE_HDIAP_09[sa1stecodes[proposed_location[,1]],proposed_location[,2]] <- true_STE_HDIAP_09[sa1stecodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_IRSAD_HDIAP_09[irsad_codes[proposed_location[,1]],proposed_location[,2]] <- true_IRSAD_HDIAP_09[irsad_codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_RA_HDIAP_09[racodes[proposed_location[,1]],proposed_location[,2]] <- true_RA_HDIAP_09[racodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
      }
      if (proposed_location[,3]==2) {
        true_SA1_HDIAP_1014[proposed_location[,1],proposed_location[,2]] <- true_SA1_HDIAP_1014[proposed_location[,1],proposed_location[,2]] + proposed_move
        true_SA2_HDIAP_1014[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] <- true_SA2_HDIAP_1014[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA3_HDIAP_1014[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] <- true_SA3_HDIAP_1014[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA4_HDIAP_1014[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] <- true_SA4_HDIAP_1014[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_STE_HDIAP_1014[sa1stecodes[proposed_location[,1]],proposed_location[,2]] <- true_STE_HDIAP_1014[sa1stecodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_IRSAD_HDIAP_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]] <- true_IRSAD_HDIAP_1014[irsad_codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_RA_HDIAP_1014[racodes[proposed_location[,1]],proposed_location[,2]] <- true_RA_HDIAP_1014[racodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
      }
      if (proposed_location[,3]==3) {
        true_SA1_HDIAP_1519[proposed_location[,1],proposed_location[,2]] <- true_SA1_HDIAP_1519[proposed_location[,1],proposed_location[,2]] + proposed_move
        true_SA2_HDIAP_1519[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] <- true_SA2_HDIAP_1519[sa1sa2codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA3_HDIAP_1519[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] <- true_SA3_HDIAP_1519[sa1sa3codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_SA4_HDIAP_1519[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] <- true_SA4_HDIAP_1519[sa1sa4codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_STE_HDIAP_1519[sa1stecodes[proposed_location[,1]],proposed_location[,2]] <- true_STE_HDIAP_1519[sa1stecodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_IRSAD_HDIAP_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]] <- true_IRSAD_HDIAP_1519[irsad_codes[proposed_location[,1]],proposed_location[,2]] + proposed_move
        true_RA_HDIAP_1519[racodes[proposed_location[,1]],proposed_location[,2]] <- true_RA_HDIAP_1519[racodes[proposed_location[,1]],proposed_location[,2]] + proposed_move
      }
      
    }
  }
  
  if ((z %% 10000)==1) {cat(sum(z/10000000),sum(true_SA1_HDIAP_AGE[,1,1])/sum(STE_num_dwellings_by_DATUM_obs[,1])," ",sum(true_SA1_HDIAP_AGE[,1,2])/sum(STE_num_dwellings_by_DATUM_obsX[,1])," ",sum(true_SA1_HDIAP_AGE[,1,3])/sum(STE_num_dwellings_by_DATUM_obsY[,1]),"\n")}
  if ((z %% 10000)==1) {
    xx <- true_IRSAD_HDIAP_09[,1]/(true_IRSAD_HDIAP_09[,1]+true_IRSAD_HDIAP_09[,2])
    yy <- IRSAD_by_DATUM_09[,1]/(IRSAD_by_DATUM_09[,1]+IRSAD_by_DATUM_09[,2])
    cat(sum(z/10000000),(xx-yy)/yy,"\n")}
  if ((z %% 10000)==1) {
    xx <- true_IRSAD_HDIAP_1014[,1]/(true_IRSAD_HDIAP_1014[,1]+true_IRSAD_HDIAP_1014[,2])
    yy <- IRSAD_by_DATUM_1014[,1]/(IRSAD_by_DATUM_1014[,1]+IRSAD_by_DATUM_1014[,2])
    cat(sum(z/10000000),(xx-yy)/yy,"\n")}
  if ((z %% 10000)==1) {
    xx <- true_IRSAD_HDIAP_1519[,1]/(true_IRSAD_HDIAP_1519[,1]+true_IRSAD_HDIAP_1519[,2])
    yy <- IRSAD_by_DATUM_1519[,1]/(IRSAD_by_DATUM_1519[,1]+IRSAD_by_DATUM_1519[,2])
    cat(sum(z/10000000),(xx-yy)/yy,"\n")}
  if ((z%%1000000)==0) {save(true_SA1_HDIAP_AGE,file=paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",count)))}
} # Z

### Process results

lgas <- vect("/mnt/Z/ewan/GEOMETRIES/POLYGONS/LGA_2021_AUST_GDA2020.shp")
N_LGA <- length(lgas$LGA_NAME21)
lga_names <- lgas$LGA_NAME21

library(rjson)
ndss_0_9_lga <- fromJSON(file='NDSS_0_9_lga.json')
alt_lga_names <- list(N_LGA)
for (i in 1:length(ndss_0_9_lga$MultiPolygons)) {
  alt_lga_names[[i]] <- ndss_0_9_lga$MultiPolygons[[i]]$Name
}
alt_lga_names <- unlist(alt_lga_names)
alt_lga_names <- gsub('\ \\(C\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(M\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(T\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(S\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(R\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(A\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(B\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(DC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(AC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RegC\\)','',alt_lga_names)

ndss_stat_value <- numeric(N_LGA)+NA
for (i in 1:length(ndss_0_9_lga$MultiPolygons)) {
  ndss_stat_value[match(alt_lga_names[[i]],lga_names)] <- ndss_0_9_lga$MultiPolygons[[i]]$StatisticValue
}
genuine_na_LGA_0_9 <- which(is.na(ndss_stat_value))
suppressed_LGA_0_9 <- which(!is.na(ndss_stat_value) & ndss_stat_value==-2)
ndss_stat_value[genuine_na_LGA_0_9] <- NA
ndss_stat_value[suppressed_LGA_0_9] <- NA
ndss_LGA_0_9 <- ndss_stat_value

ndss_10_19_lga <- fromJSON(file='NDSS_10_19_lga.json')
alt_lga_names <- list(N_LGA)
for (i in 1:length(ndss_10_19_lga$MultiPolygons)) {
  alt_lga_names[[i]] <- ndss_10_19_lga$MultiPolygons[[i]]$Name
}
alt_lga_names <- unlist(alt_lga_names)
alt_lga_names <- gsub('\ \\(C\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(M\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(T\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(S\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(R\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(A\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(B\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(DC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(AC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RC\\)','',alt_lga_names)
alt_lga_names <- gsub('\ \\(RegC\\)','',alt_lga_names)

ndss_stat_value <- numeric(N_LGA)+NA
for (i in 1:length(ndss_10_19_lga$MultiPolygons)) {
  ndss_stat_value[match(alt_lga_names[[i]],lga_names)] <- ndss_10_19_lga$MultiPolygons[[i]]$StatisticValue
}
genuine_na_LGA_10_19 <- which(is.na(ndss_stat_value))
suppressed_LGA_10_19 <- which(!is.na(ndss_stat_value) & ndss_stat_value==-2)
ndss_stat_value[genuine_na_LGA_10_19] <- NA
ndss_stat_value[suppressed_LGA_10_19] <- NA
ndss_LGA_10_19 <- ndss_stat_value

mb2sa1 <- read_xlsx('MB_2021_AUST.xlsx')
mb2sa1 <- mb2sa1[which(mb2sa1$SA1_CODE_2021 %in% unique(sa1_2021$SA1_CODE21)),]
mbsa1codes <- match(mb2sa1$SA1_CODE_2021,sa1_2021$SA1_CODE21)

mb2lga <- read_xlsx('LGA_2021_AUST.xlsx')
mb2lga <- mb2lga[match(mb2sa1$MB_CODE_2021,mb2lga$MB_CODE_2021),]
mblgacodes <- match(mb2lga$LGA_NAME_2021,lga_names)

sa12lga <- aggregate(rep(1,length(mbsa1codes)),list(mbsa1codes*10^5+mblgacodes),sum)
sa12lga$sa1 <- floor(sa12lga$Group.1/10^5)
sa12lga$lga <- sa12lga$Group.1 %% 10^5
sa12lga <- sa12lga[sort.list(sa12lga$sa1),]

sa12lgax <- sparseMatrix(sa12lga$sa1,sa12lga$lga,x=sa12lga$x,dims = c(length(sa1_2021$SA1_CODE21),N_LGA))
normalisation <- as.numeric(sa12lgax%*%rep(1,N_LGA))
sa12lga <- sparseMatrix(sa12lga$lga,sa12lga$sa1,x=sa12lga$x/normalisation[sa12lga$sa1],dims = c(N_LGA,length(sa1_2021$SA1_CODE21)))

ntot <- floor(z/1000000)

quantilelow <- function(x) {quantile(x,0.025)}
quantilehigh <- function(x) {quantile(x,0.975)}
quantilemedian <- function(x) {quantile(x,0.5)}

sa1_2021_outputs <- sa1_2021

posterior_draws <- list()
posterior_draws_old <- list()
posterior_draws_NDSS_LGA <- list()
posterior_draws_NDSS_LGA_old <- list()
posterior_draws_NDSS_LGAw <- list()
posterior_draws_NDSS_LGA_oldw <- list()
tot_pops_pred <- list()
tot_pops_pred_old <- list()

for (i in 2:ntot) {
  
  load(paste0("/mnt/Z/ewan/DIABETES/xsample_final",sprintf("%03i",i)))
  load(paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",i)))
 
  for (j in 1:10) {
    xout <- ilogit(objx$report(xsample[j,])$joint_field)
    posterior_draws[[length(posterior_draws)+1]] <- (xout[,1]*2/3 + xout[,2]/3)
    posterior_draws_old[[length(posterior_draws_old)+1]] <- xout[,3]
  
    pd <- rowSums(true_SA1_HDIAP_AGE[,,1])*posterior_draws[[length(posterior_draws)]]
    pd_old <- (rowSums(true_SA1_HDIAP_AGE[,,2])+rowSums(true_SA1_HDIAP_AGE[,,3]))*posterior_draws_old[[length(posterior_draws_old)]]
    posterior_draws_NDSS_LGA[[length(posterior_draws_NDSS_LGA)+1]] <- as.numeric(sa12lga%*%pd)
    posterior_draws_NDSS_LGA_old[[length(posterior_draws_NDSS_LGA_old)+1]] <- as.numeric(sa12lga%*%pd_old)
  }

  posterior_draws_NDSS_LGAw[[length(posterior_draws_NDSS_LGAw)+1]] <- as.numeric(sa12lga%*%true_SA1_HDIAP_AGE[,1,1])
  posterior_draws_NDSS_LGA_oldw[[length(posterior_draws_NDSS_LGA_oldw)+1]] <- as.numeric(sa12lga%*%(true_SA1_HDIAP_AGE[,1,2]+true_SA1_HDIAP_AGE[,1,3]))
  tot_pops_pred[[length(tot_pops_pred)+1]] <- rowSums(true_SA1_HDIAP_AGE[,,1])
  tot_pops_pred_old[[length(tot_pops_pred_old)+1]] <- rowSums((true_SA1_HDIAP_AGE[,,1]+true_SA1_HDIAP_AGE[,,2]))
  
  cat(i,"\n")
}

posterior_draws <- do.call(rbind,posterior_draws)
posterior_draws_old <- do.call(rbind,posterior_draws_old)
posterior_draws_NDSS_LGA <- do.call(rbind,posterior_draws_NDSS_LGA)
posterior_draws_NDSS_LGA_old <- do.call(rbind,posterior_draws_NDSS_LGA_old)
posterior_draws_NDSS_LGAw <- do.call(rbind,posterior_draws_NDSS_LGAw)
posterior_draws_NDSS_LGA_oldw <- do.call(rbind,posterior_draws_NDSS_LGA_oldw)

posterior_risk_median_young <- apply((posterior_draws),2,median)
posterior_risk_median_old <- apply((posterior_draws_old),2,median)
posterior_risk_median_NDSS_LGA <- apply((posterior_draws_NDSS_LGA),2,median)
posterior_risk_median_NDSS_LGA_old <- apply((posterior_draws_NDSS_LGA_old),2,median)
posterior_risk_median_NDSS_LGAw <- apply((posterior_draws_NDSS_LGAw),2,median)
posterior_risk_median_NDSS_LGA_oldw <- apply((posterior_draws_NDSS_LGA_oldw),2,median)
posterior_risk_low_NDSS_LGAw <- apply((posterior_draws_NDSS_LGAw),2,quantilelow)
posterior_risk_low_NDSS_LGA_oldw <- apply((posterior_draws_NDSS_LGA_oldw),2,quantilelow)
posterior_risk_high_NDSS_LGAw <- apply((posterior_draws_NDSS_LGAw),2,quantilehigh)
posterior_risk_high_NDSS_LGA_oldw <- apply((posterior_draws_NDSS_LGA_oldw),2,quantilehigh)

posterior_risk_sd_young <- apply((posterior_draws),2,sd)
posterior_risk_sd_old <- apply((posterior_draws_old),2,sd)

sa1_2021_outputs$RMEAY <- posterior_risk_median_young
sa1_2021_outputs$RMEAO <- posterior_risk_median_old
sa1_2021_outputs$SSMEY <- posterior_risk_sd_young
sa1_2021_outputs$SSMEO <- posterior_risk_sd_old

exceedence_draws <- posterior_draws>median(sa1_2021_outputs$RMEAY)
exceedence <- apply(exceedence_draws,2,mean)
sa1_2021_outputs$XXXY <- exceedence
exceedence_draws <- posterior_draws_old>median(sa1_2021_outputs$RMEAO)
exceedence <- apply(exceedence_draws,2,mean)
sa1_2021_outputs$XXXO <- exceedence

surprised_lga <- numeric(N_LGA)
for (i in 1:N_LGA) {
  xsample <- posterior_draws_NDSS_LGAw[,i]
  xsample[xsample<20] <- 0
  xsample <- round(xsample/10)*10
  xsample[xsample==0] <- -99
  surprised_lga[i] <- as.integer(mean(ndss_LGA_0_9[i]<=xsample) < 0.05 | mean(ndss_LGA_0_9[i]>=xsample) < 0.05)
}
surprised_lga_old <- numeric(N_LGA)
for (i in 1:N_LGA) {
  xsample <- posterior_draws_NDSS_LGA_oldw[,i]
  xsample[xsample<20] <- 0
  xsample <- round(xsample/10)*10
  xsample[xsample==0] <- -99
  surprised_lga_old[i] <- as.integer(mean(ndss_LGA_10_19[i]<=xsample) < 0.05 | mean(ndss_LGA_10_19[i]>=xsample) < 0.05)
}

lga_2021_outputs <- lgas
lga_2021_outputs$RMEAY <- posterior_risk_median_NDSS_LGAw
lga_2021_outputs$RMEAO <- posterior_risk_median_NDSS_LGA_oldw
lga_2021_outputs$LMEAY <- posterior_risk_low_NDSS_LGAw
lga_2021_outputs$LMEAO <- posterior_risk_low_NDSS_LGA_oldw
lga_2021_outputs$UMEAY <- posterior_risk_high_NDSS_LGAw
lga_2021_outputs$UMEAO <- posterior_risk_high_NDSS_LGA_oldw
lga_2021_outputs$OMEAY <- ndss_LGA_0_9
lga_2021_outputs$OMEAO <- ndss_LGA_10_19

lga_2021_outputs$surp <- surprised_lga
lga_2021_outputs$surpo <- surprised_lga_old

writeVector(lga_2021_outputs,"diabetesLGA",filetype="ESRI Shapefile",overwrite=TRUE)

tot_pops_pred <- do.call(rbind,tot_pops_pred)
tot_pops_pred_old <- do.call(rbind,tot_pops_pred_old)
tot_pops_pred_mean <- apply(tot_pops_pred,2,mean)
tot_pops_pred_mean_old <- apply(tot_pops_pred_old,2,mean)
sa1_2021_outputs$NPOPO <- tot_pops_pred_mean_old
sa1_2021_outputs$NPOPY <- tot_pops_pred_mean

writeVector(sa1_2021_outputs,"diabetes",filetype="ESRI Shapefile",overwrite=TRUE)

ses_effect <- list()
ses_effect_old <- list()
ra_effect <- list()
ra_effect_old <- list()
prop <- list()
prop_old <- list()

for (i in 2:ntot) {

    load(paste0("/mnt/Z/ewan/DIABETES/xsample_final",sprintf("%03i",i)))

    for (j in 1:10) {
      ses_effect[[length(ses_effect)+1]] <- xsample[j,colnames(xsample)=="ses_effects_young"]
      ses_effect_old[[length(ses_effect_old)+1]] <- xsample[j,colnames(xsample)=="ses_effects_old"]
      ra_effect[[length(ra_effect)+1]] <- xsample[j,colnames(xsample)=="ra_effects_young"]
      ra_effect_old[[length(ra_effect_old)+1]] <- xsample[j,colnames(xsample)=="ra_effects_old"]
      prop[[length(prop)+1]] <- ilogit(xsample[j,colnames(xsample)=="logit_middle_age_prop"])
    }
}

ses_effect <- do.call(rbind,ses_effect)
ses_effect_old <- do.call(rbind,ses_effect_old)
ra_effect <- do.call(rbind,ra_effect)
ra_effect_old <- do.call(rbind,ra_effect_old)
prop <- do.call(rbind,prop)

ses_effect <- ses_effect-ses_effect[,1]
ses_effect_old <- ses_effect_old-ses_effect_old[,1]
ra_effect <- ra_effect-ra_effect[,1]
ra_effect_old <- ra_effect_old-ra_effect_old[,1]

quantilelow <- function(x) {quantile(x,0.025)}
quantilehigh <- function(x) {quantile(x,0.975)}

ses_effect_median <- apply((ses_effect),2,median)
ses_effect_low <- apply((ses_effect),2,quantilelow)
ses_effect_high <- apply((ses_effect),2,quantilehigh)
cat(sprintf("%3.2f",ses_effect_median),"\n")
cat(sprintf("%3.2f",ses_effect_low),"\n")
cat(sprintf("%3.2f",ses_effect_high),"\n")

ses_effect_old_median <- apply((ses_effect_old),2,median)
ses_effect_old_low <- apply((ses_effect_old),2,quantilelow)
ses_effect_old_high <- apply((ses_effect_old),2,quantilehigh)
cat(sprintf("%3.2f",ses_effect_old_median),"\n")
cat(sprintf("%3.2f",ses_effect_old_low),"\n")
cat(sprintf("%3.2f",ses_effect_old_high),"\n")

ra_effect_median <- apply((ra_effect),2,median)
ra_effect_low <- apply((ra_effect),2,quantilelow)
ra_effect_high <- apply((ra_effect),2,quantilehigh)
cat(sprintf("%3.2f",ra_effect_median),"\n")
cat(sprintf("%3.2f",ra_effect_low),"\n")
cat(sprintf("%3.2f",ra_effect_high),"\n")

ra_effect_old_median <- apply((ra_effect_old),2,median)
ra_effect_old_low <- apply((ra_effect_old),2,quantilelow)
ra_effect_old_high <- apply((ra_effect_old),2,quantilehigh)
cat(sprintf("%3.2f",ra_effect_old_median),"\n")
cat(sprintf("%3.2f",ra_effect_old_low),"\n")
cat(sprintf("%3.2f",ra_effect_old_high),"\n")

prop_median <- apply((prop),2,median)
prop_low <- apply((prop),2,quantilelow)
prop_high <- apply((prop),2,quantilehigh)
cat(sprintf("%3.2f",prop_median),"\n")
cat(sprintf("%3.2f",prop_low),"\n")
cat(sprintf("%3.2f",prop_high),"\n")

young2counts <- read.csv("T1D_T2D counts at -10 yrs at SA2 in 2021.csv")

sa2_code_matches <- match(young2counts$SA2_CODE_2021,unique(sa1_2021$SA2_CODE21))

sa2_tots <- young2counts$T1D_cases+young2counts$T2D_cases

sa2sa3 <- match(sa1_2021$SA3_CODE21[!duplicated(sa1_2021$SA2_CODE21)],unique(sa1_2021$SA3_CODE21))
xx <- aggregate(sa2_tots,list(sa2sa3[sa2_code_matches]),sum)
sa3_tots <- xx$x

tot_sa2s_pred <- list()
for (i in 2:ntot) {
  load(paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",i)))
  tot_sa2s_pred[[length(tot_sa2s_pred)+1]] <- aggregate(true_SA1_HDIAP_AGE[,1,1],list(sa1sa2codes),sum)[,-1][sa2_code_matches]
}
tot_sa2s_pred <- do.call(rbind,tot_sa2s_pred)

tot_sa2s_pred_lower <- apply(tot_sa2s_pred,2,quantilelow)
tot_sa2s_pred_upper <- apply(tot_sa2s_pred,2,quantilehigh)
tot_sa2s_pred_median <- apply(tot_sa2s_pred,2,quantilemedian)
tot_sa2s_pred_mean <- apply(tot_sa2s_pred,2,mean)
tot_sa2s_pred_sd <- apply(tot_sa2s_pred,2,sd)

sa2_tots_jittered <- sa2_tots+runif(length(sa2_tots),-0.25,0.25)
tot_sa2s_pred_median_jittered <- tot_sa2s_pred_median+runif(length(sa2_tots),-0.1,0.1)

tot_sa3s_pred <- list()
for (i in 2:ntot) {
  load(paste0("/mnt/Z/ewan/DIABETES/true_SA1_HDIAP_AGE_saved",sprintf("%03i",i)))
  tot_sa3s_pred[[length(tot_sa3s_pred)+1]] <- aggregate(true_SA1_HDIAP_AGE[,1,1],list(sa1sa3codes),sum)[,-1][xx$Group.1]
}
tot_sa3s_pred <- do.call(rbind,tot_sa3s_pred)

tot_sa3s_pred_lower <- apply(tot_sa3s_pred,2,quantilelow)
tot_sa3s_pred_upper <- apply(tot_sa3s_pred,2,quantilehigh)
tot_sa3s_pred_median <- apply(tot_sa3s_pred,2,quantilemedian)
tot_sa3s_pred_mean <- apply(tot_sa3s_pred,2,mean)
tot_sa3s_pred_sd <- apply(tot_sa3s_pred,2,sd)

sa3_tots_jittered <- sa3_tots+runif(length(sa3_tots),-0.25,0.25)
tot_sa3s_pred_median_jittered <- tot_sa3s_pred_mean+runif(length(sa3_tots),-0.1,0.1)

ABS_sa3s_pred_median_jittered <- SA3_num_dwellings_by_DATUM_obs[xx$Group.1,1]+runif(length(sa3_tots),-0.1,0.1)

cor(SA3_num_dwellings_by_DATUM_obs[xx$Group.1,1],sa3_tots)
cor(tot_sa3s_pred_mean,sa3_tots)


pdf("figa.pdf",width=5,height=5)
layout(1)
par(mai=c(1,1,0.05,0.05))

plot(tot_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.4,col="white",xlab="",ylab="",xlim=c(-0.5,23.5),ylim=c(-0.5,23.5),xaxt='n',yaxt='n',xaxs='i',yaxs='i')
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(i-0.5,-0.5,i+0.5,23.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(-0.5,i-0.5,23.5,i+0.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
lines(c(-0.5,23.5),c(-0.5,23.5),col="grey50")
points(tot_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.7,col=hsv(0.82,v=1,alpha=0.25))
for (i in 1:length(tot_sa3s_pred_median_jittered)) {
  lines(c(tot_sa3s_pred_lower[i]+tot_sa3s_pred_median_jittered[i]-tot_sa3s_pred_median[i],tot_sa3s_pred_upper[i]+tot_sa3s_pred_median_jittered[i]-tot_sa3s_pred_median[i]),c(sa3_tots_jittered[i],sa3_tots_jittered[i]),col=hsv(0.82,v=1,alpha=0.25))
}
dev.off()

pdf("figb.pdf",width=5,height=5)
layout(1)
par(mai=c(1,1,0.05,0.05))

plot(ABS_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.4,col="white",xlab="",ylab="",xlim=c(-0.5,23.5),ylim=c(-0.5,23.5),xaxt='n',yaxt='n',xaxs='i',yaxs='i')
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(i-0.5,-0.5,i+0.5,23.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25)) {
  rect(-0.5,i-0.5,23.5,i+0.5,col = hsv(0,s=0,v=0.9,alpha=0.2),border=NA)
}
lines(c(-0.5,23.5),c(-0.5,23.5),col="grey50")
points(ABS_sa3s_pred_median_jittered,sa3_tots_jittered,pch=19,cex=0.7,col=hsv(0.82,v=1,alpha=0.25))
dev.off()



