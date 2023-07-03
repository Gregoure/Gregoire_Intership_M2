#### SCRIPT FOR SIM-BY-SIM BRANCHING OBSERVEATION #####

#Try for analysis of multiple output files (densities & traits) for a particular range of parameters
#in order to produce figure of evolutionnary trajectories as function of parameter
library(readr)
library(dplyr)
library(plot3D)
library(plotly)
library(reshape2)
library(mixtools)
################# Get Files #######
dir <- choose.dir() # Has to go to global folder in which all sims are stored in subfolders

#### Build a list of folder DENSITIES & TRAITS ####
#Read densities files
cwd <- setwd(dir)
listfiles <- list.files(cwd)
listfiles


### FUNCTIONS USED TO BUILD MEAN SIMULATION ###
MeanSimValue <- function(lsdf, var){ # return mean value for each time of a given variable
  
  nbsims <- length(lsdf)
  values <- vector('list', nbsims)
  
  for (i in 1: nbsims){
    idcol <- which(colnames(lsdf[[i]])== var)
    values[[i]] <- lsdf[[i]][,idcol]
  }
  as.data.frame(values)
  meanvalues <- rowSums(as.data.frame(values), na.rm = T)/nbsims
  meanvalues
  
} 



BuildMeanSim <- function(lsdf){ # Build the mean dataframe from the whole set of simulations
  namecols<- names(lsdf[[1]])
  MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(lsdf[[1]]$Time)))
  colnames(MeanSim) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(lsdf,namecols[i])
    MeanSim[,i] <- meanvar
  }
}


#####
ls_files_Distribtions <- list.files(dir) # list files
ls_files_Distribtions[2]
dir
#currentdir <- paste(dir, ls_files_Distribtions[10], sep = "\\") #Browse subfolders in main folders
#DistribDir <- paste(currentdir, "Distributions", sep = "\\") # Get the subsubfolder with traits files
#cwd <- setwd(DistribDir)
DistribDir <- dir
ls_files_Distrib <- list.files(DistribDir)
ls_files_Distrib
list_tibbles_distrib <- lapply(ls_files_Distrib, read_csv) # Read CSV returns tibble object but i don't know how to deal with thoses guys so...
list_df_Distributions <- lapply(list_tibbles_distrib, as.data.frame) #... we turn them to good old dataframes

DistribFiles <- c()

# I. Shaping dataframes and getting common rows

for (i in 1:length(list_df_Distributions)){
  
  list_df_Distributions[[i]]<- subset(list_df_Distributions[[i]][c(-1)]) # suppress first column
  list_df_Distributions[[i]]$Time <- round(list_df_Distributions[[i]]$Time, 1) #round t values
  list_df_Distributions[[i]] <- distinct(list_df_Distributions[[i]], Time, .keep_all = TRUE) # remove duplicates
  
}
#Get only commons rows between all dataframes 
#Get all t values
t_values <- c(0)
for (i in 1: length(list_df_Distributions)){
  times <- list_df_Distributions[[i]]$Time
  t_values <- c(t_values, times)
}
t_values <- t_values[-1] # Suppress the value given for initialization
counts <- as.data.frame(table(t_values)) # Count occurrences into a table
counts
#Get values that appear nbsim times (common to every df)
nbsim <- length(list_df_Distributions)
tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df

#loop to keep only those t in the whole set of dfs
for (i in 1:length(list_df_Distributions)){
  indexrows <- which(list_df_Distributions[[i]]$Time %in% tcommon)
  list_df_Distributions[[i]] <- subset(list_df_Distributions[[i]][indexrows,])
}

# II. Building the mean simulation ( returns a single DataFrame)
namecols<- names(list_df_Distributions[[1]])
MeanSim3 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_Distributions[[1]]$Time)))
colnames(MeanSim3) <- namecols

for (i in 1:length(namecols)){
  
  meanvar <- MeanSimValue(list_df_Distributions,namecols[i])
  MeanSim3[,i] <- meanvar
}


DistribFiles <- c(DistribFiles, list(MeanSim3))

Mydf <- list_df_Distributions[[8]]
Mydf <- DistribFiles[[1]]
alphas <- colnames(Mydf[2:length(Mydf)])
b <- as.numeric(gsub("Alpha", "", alphas))
b
MySubDf <- Mydf[-c(1:40),-c(1)]
Finaldata <- data.matrix(MySubDf, rownames.force = NA) 
fig <- plot_ly(
  
  x = b, y = Mydf$Time,
  
  z = Finaldata, type = "heatmap"
  
)


fig

# Try for Multinomial mixture model using Mixtools
#Select row
Mydistrib <- Finaldata[1000,]
Mydistrib
#Create the distribution
x<-rep(c(0.18,0.19,0.20,0.21,0.22, 0.33,0.34,0.35,0.36,0.37),times=c(13,106,237,160,1,13,106,237,160,1))
hist(x)
wait1 <- normalmixEM(x)
wait1$loglik
plot(wait1, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, xlab2="Minutes")

for (i in length(alphas)){
  str_remove(alphas[i], "Alpha")
}
b <- as.numeric(gsub("Alpha", "", alphas))
b
MultDatas <- makemultdata(Finaldata, cuts = b)
Mod <-multmixEM(MultDatas)
cdf3 <- compCDF(Finaldata, Mod$posterior, lwd=2, lab=c(7, 5, 7),xlab="Phenotypes", ylab="Component CDFs")
Mod$y
nrow(MultDatas)
