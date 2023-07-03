rm(list = ls())
graphics.off()
dev.off()
#Try for analysis of multiple output files (densities & traits) for a particular range of parameters
#in order to produce figure of evolutionnary trajectories as function of parameter
library(readr)
library(dplyr)
library(plot3D)
library(plotly)
library(mixtools)
library(factoextra)
library(NbClust)
library(tidyr)
library(cluster)

### FUNCTIONS USED TO BUILD MEAN SIMULATION ###
MeanSimValue <- function(lsdf, var){ # return mean value for each time of a given variable
  
  nbsims <- length(lsdf)
  values <- vector('list', nbsims)
  
  for (i in 1: nbsims){
    idcol <- which(colnames(lsdf[[i]])== var)
    values[[i]] <- lsdf[[i]][,idcol]
  }
  as.data.frame(values)
  meanvalues <- rowMeans(as.data.frame(values), na.rm = T)
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


#################Get Files #######
dir <- choose.dir() # Has to go to global folder in which all sims are stored in subfolders

#### Build a list of folder DENSITIES & TRAITS ####
#Read densities files
cwd <- setwd(dir)
cwd
listfiles <- list.files(cwd)
listfiles

TraitsFiles <- list()

for (i in 1:length(listfiles)){
  i=9
  currentdir <- paste(dir, listfiles[i], sep = "\\") #Browse subfolders in main folders
  TraitsDir <- paste(currentdir, "Traits", sep = "\\") # Get the subsubfolder with traits files

  cwd <- setwd(TraitsDir)
  cwd
  
  ls_files_Traits <- list.files(TraitsDir) # list files
  list_tibbles_traits <- lapply(ls_files_Traits, read_csv) # Read CSV returns tibble object but i don't know how to deal with thoses guys so...
  list_df_traits <- lapply(list_tibbles_traits, as.data.frame) #... we turn them to good old dataframes
  
  
  ### This part contruct the mean simulation from the Traits list of file we got earlier ###
  #Can be improved using Xapply mastery
  
  # I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_traits)){
    
    list_df_traits[[i]]<- subset(list_df_traits[[i]][c(-1)]) # suppress first column
    list_df_traits[[i]]$Time <- round(list_df_traits[[i]]$Time, 1) #round t values
    list_df_traits[[i]] <- distinct(list_df_traits[[i]], Time, .keep_all = TRUE)
    
  }
  #Get only commons rows between all dataframes 
  #Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_traits)){
    times <- list_df_traits[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  counts
  #Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_traits)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  tcommon
  #loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_traits)){
    indexrows <- which(list_df_traits[[i]]$Time %in% tcommon)
    list_df_traits[[i]] <- subset(list_df_traits[[i]][indexrows,])
  }
  
  # II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_traits[[1]])
  MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_traits[[1]]$Time)))
  colnames(MeanSim) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_traits,namecols[i])
    MeanSim[,i] <- meanvar
  }
  
  
  #And add to the list
  TraitsFiles <- c(TraitsFiles, list(MeanSim))
  
  
} # GIANT LOOP THAT EXTRACT, COMPILE AND BUILD A MEAN SIM WITH ALL FILES

for (i in 1:length(listfiles)){ 
  i=8
  currentdir <- paste(dir, listfiles[i], sep = "\\") #Browse subfolders in main folders
  TraitsDir <- paste(currentdir, "Traits", sep = "\\") # Get the subsubfolder with traits files
  DensitiesDir <- paste(currentdir, "Densities", sep = "\\") # Get the subsubfolder with densities files
  DistribDir <- paste(currentdir, "Distributions", sep = "\\") # Get the subfolder for distribution of traits when discrete trait values are used
  
  cwd <- setwd(TraitsDir)
  cwd
  
  ### Using a single run
  j = 5 # Files desired for the one run simulation
  ls_files_Traits <- list.files(TraitsDir) 
  list_tibbles_traits <- lapply(ls_files_Traits[j], read_csv) 
  list_df_traits <- lapply(list_tibbles_traits, as.data.frame)
  
  ### This part contruct the mean simulation from the Traits list of file we got earlier ###
  #Can be improved using Xapply mastery
  
  # I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_traits)){
    
    list_df_traits[[i]]<- subset(list_df_traits[[i]][c(-1)]) # suppress first column
    list_df_traits[[i]]$Time <- round(list_df_traits[[i]]$Time, 1) #round t values
    list_df_traits[[i]] <- distinct(list_df_traits[[i]], Time, .keep_all = TRUE)
    
  }
  #Get only commons rows between all dataframes 
  #Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_traits)){
    times <- list_df_traits[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  counts
  #Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_traits)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  tcommon
  #loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_traits)){
    indexrows <- which(list_df_traits[[i]]$Time %in% tcommon)
    list_df_traits[[i]] <- subset(list_df_traits[[i]][indexrows,])
  }
  
  # II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_traits[[1]])
  MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_traits[[1]]$Time)))
  colnames(MeanSim) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_traits,namecols[i])
    MeanSim[,i] <- meanvar
  }
  
  
  #And add to the list
  TraitsFiles <- c(TraitsFiles, list(MeanSim))
  
  
  ### Now we repeat the same process for Densities
  
  cwd <- setwd(DensitiesDir)
  
  ### Using a single run
  ls_files_Densities <- list.files(DensitiesDir) 
  list_tibbles_densities <- lapply(ls_files_Densities[j], read_csv) 
  list_df_densities <- lapply(list_tibbles_densities, as.data.frame)
  
  # I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_densities)){
    list_df_densities[[i]]<- subset(list_df_densities[[i]][c(-1)]) # suppress first column
    list_df_densities[[i]]$Time <- round(list_df_densities[[i]]$Time, 1) #round t values
    list_df_densities[[i]] <- distinct(list_df_densities[[i]], Time, .keep_all = TRUE)#remove duplicates
    
  }
  #Get only commons rows between all dataframes 
  #Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_densities)){
    times <- list_df_densities[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  counts
  #Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_densities)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  #loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_densities)){
    indexrows <- which(list_df_densities[[i]]$Time %in% tcommon)
    list_df_densities[[i]] <- subset(list_df_densities[[i]][indexrows,])
  }
  
  # II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_densities[[1]])
  MeanSim2 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_densities[[1]]$Time)))
  colnames(MeanSim2) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_densities,namecols[i])
    MeanSim2[,i] <- meanvar
  }
  
  
  DensitiesFiles <- c(DensitiesFiles, list(MeanSim2))
  
  
  # Now we re-repeat the same process for distributions
  ### Now we repeat the same process for Densities
  
  cwd <- setwd(DistribDir)
  
  ### Using a single run
  ls_files_Distributions <- list.files(DistribDir) 
  list_tibbles_Distributions <- lapply(ls_files_Distributions[j], read_csv) 
  list_df_Distributions <- lapply(list_tibbles_Distributions, as.data.frame)
  
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
  
  
} # GIANT LOOP THAT SERVES TO DO SINGLE RUN ET SERT A ETRE UTILISER PARTIELLEMENT


########################
i = 1 # Files desired for analysis in the repository
########################

################# Analysis of Traits Files #####################################

ActualDf = list_df_traits[1]
ActualDf = TraitsFiles[[1]]
FinalDf <- subset(ActualDf[c(-1)])
#ActualDf = TraitsFiles[1]
plot(ActualDf$site1, type = 'l', col = 1, xlim=c(1500,2500))
for (i in 2:(length(ActualDf)/8)){
  lines(ActualDf[,i],type = 'l', col = i)
}



