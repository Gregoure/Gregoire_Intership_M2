# Petit nettoyage
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
library(BAMMtools)

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

DensitiesFiles <- list()
TraitsFiles <- list()
DistribFiles <- list()

for (i in 1:length(listfiles)){
  currentdir <- paste(dir, listfiles[i], sep = "\\") #Browse subfolders in main folders
  TraitsDir <- paste(currentdir, "Traits", sep = "\\") # Get the subsubfolder with traits files
  DensitiesDir <- paste(currentdir, "Densities", sep = "\\") # Get the subsubfolder with densities files
  DistribDir <- paste(currentdir, "Distributions", sep = "\\") # Get the subfolder for distribution of traits when discrete trait values are used
  
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
  
  
  ### Now we repeat the same process for Densities
  
  cwd <- setwd(DensitiesDir)
  
  ls_files_Densities <- list.files(DensitiesDir) #list files
  list_tibbles_densities <- lapply(ls_files_Densities, read_csv) 
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
  
  ls_files_Distribtions <- list.files(DistribDir) #list files
  list_tibbles_Distributions <- lapply(ls_files_Distribtions, read_csv) 
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


### DEMOGRAPHY ANALYSIS ###

######### Percentage of infected sites with au moins 1 infectés in the metapopulation at each time ############ --------------------------------------------
i = 2
ActualDf <- DensitiesFiles[[i]]
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])
plot(Ipops$I27,type = "l")
plot(Spops$S27, type = "l")
Ipourcents  = data.frame(matrix(NA, nrow = nrow(Ipops), ncol = length(Ipops)))
for (i in 1:length(Ipops)){
  for (j in 1:nrow(Ipops)){
    Ipourcents[j,i] = (Ipops[j,i]*100)/(Spops[j,i]+Ipops[j,i])
  }
}

Ipops[Ipops < 1] <- 0
Ipops[Ipops >= 1 ] <- 1
Imean <- rowSums(Ipops)
plot(ActualDf$Time, Imean, type='l', xlim = c(0,1500), ylim = c(50,80))

### Pourcentage d'individus infectés dans un site, moyenne de tous les sites
Imeanprc <- rowSums(Ipourcents)/length(Ipops)
plot(ActualDf$Time, Imeanprc, type='l', xlim = c(0,1500), ylim = c(0,100))





####### Average Dynamics of each metapopulation ##### -----------------------------------------------------------------------------
i=1
ActualDf <- DensitiesFiles[[i]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean

Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean, I=Imean, N = Nmean)

plot(Metapop_Data$t, Metapop_Data$N, type = 'l', col='black', main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,7000), xlim = c(0,500))
lines(Metapop_Data$t,Metapop_Data$S, col = 'blue')
lines(Metapop_Data$t, Metapop_Data$I, col = 'red')
legend("topright",legend = c("N","S", "I"), col = c("black", "blue", "red"), lty = 1 )

# Converge values for susceptible density
taillemax <- length(Smean)
ConvergeValueS <- Smean[(taillemax-1000):taillemax]
ConvergeValuesS <- sum(ConvergeValueS)/length(ConvergeValueS)
PopsiteS = ConvergeValuesS/length(Spops)

# Converge values for infected density
taillemax <- length(Imean)
ConvergeValueI <- Imean[(taillemax-1000):taillemax]
ConvergeValuesI <- sum(ConvergeValueI)/length(ConvergeValueI)
PopsiteI = ConvergeValuesI/length(Ipops)




##### Pourcentage de sites où on a une certaine densité d'infectés avec une densité de la pop conséquente

### PAS FONCTIONNEL #######

i=1
ActualDf <- DensitiesFiles[[i]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean

## Converge values for all sites
# Converge values for susceptible density
taillemax <- length(Smean)
ConvergeValueS <- Smean[(taillemax-1000):taillemax]
ConvergeValuesS <- sum(ConvergeValueS)/length(ConvergeValueS)
PopsiteS <- ConvergeValuesS/length(Spops)
# Converge values for infected density
ConvergeValueI <- Imean[(taillemax-1000):taillemax]
ConvergeValuesI <- sum(ConvergeValueI)/length(ConvergeValueI)
PopsiteI <- ConvergeValuesI/length(Ipops)


Ipourcents  = data.frame(matrix(NA, nrow = nrow(Ipops), ncol = length(Ipops)))
for (i in 1:length(Ipops)){      # Pourcentage d'infectés
  for (j in 1:nrow(Ipops)){
    Ipourcents[j,i] = (Ipops[j,i]*100)/(Spops[j,i]+Ipops[j,i])
  }
}

Pr = 20
Ipourcents[Ipourcents < Pr] <- 0
Ipourcents[Ipourcents >= Pr ] <- 1
Imeanprc <- rowSums(Ipourcents)



## Dynamique de la métapop d'un seul site
Onesite <- subset(ActualDf[c(2,3)])
S_Onesite <- Onesite$S0
I_Onesite <- Onesite$I0

# Converge values for susceptible density
taillemax <- length(S_Onesite)
ConvergeValueS <- S_Onesite[(taillemax-100):taillemax]
ConvergeValuesS <- sum(ConvergeValueS)/length(ConvergeValueS)

# Converge values for infected density
ConvergeValueI <- I_Onesite[(taillemax-100):taillemax]
ConvergeValuesI <- sum(ConvergeValueI)/length(ConvergeValueI)




###### Phase Space ##### ----------------------------------------------------------------------------------------------------------
#Select only equilibrium values
plot(Metapop_Data$S, Metapop_Data$I, type = 'l', xlab = 'Susceptibles', ylab = 'Infected')


########## -----------------------------------------------------------------------------------------------------------------------

### FIGURE THAT GIVES THE TOTAL POPULATION SIZE AS A FUNCTION OF PARAMETER ####
ActualDf <- DensitiesFiles[[1]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean

Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean, I=Imean, N = Nmean)
plot(Metapop_Data$t, Metapop_Data$N, type = 'l', col=1, main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,10000))


for (i in 2:length(DensitiesFiles)){
  ActualDf <- DensitiesFiles[[i]]
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  
  Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean, I=Imean, N = Nmean)
  lines(Metapop_Data$t, Metapop_Data$N, col = i)

}
legend(1000,1000, legend = c('dI = 0.1','dI =0.3','dI =0.7'), col = c(1,2,3), lty = 1:2, cex = 0.8)




### FIGURE THAT GIVE THE SUSCEPTIBLE (AND INFECTED) DENSITY AS A FUNCTION OF PARAMETERS ####

ActualDf <- DensitiesFiles[[1]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)

Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean)
plot(Metapop_Data$t, Metapop_Data$S, type = 'l', col=1, main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,7000), xlim = c(0,1500))



for (i in 2:length(DensitiesFiles)){
  
  ActualDf <- DensitiesFiles[[i]]
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  Metapop_Data
  Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean)
  lines(Metapop_Data$t, Metapop_Data$S, col = i)
  
}
legend(1000,1000, legend = c('dI = 0.1','dI =0.3','dI =0.7'), col = c(1,2,3), lty = 1:2, cex = 0.8)


### FIGURE THAT GIVE THE INFECTED FRACTION OF POPULATION AS A FUNCTION OF PARAMETER ###

ActualDf <- DensitiesFiles[[2]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean
Ifrac <- Imean / Nmean

Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean, I=Ifrac, N = Nmean)
plot(Metapop_Data$t, Metapop_Data$I, type = 'l', col=1, main = "Metapopulation dynamics", xlab = "t", ylab = "Fraction I", ylim = c(0,1), xlim = c(0,500))


for (i in 2:length(DensitiesFiles)){
  ActualDf <- DensitiesFiles[[i]]
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  Ifrac <- Imean / Nmean
  
  Metapop_Data <- data.frame(t = ActualDf$Time, I=Ifrac)
  lines(Metapop_Data$t, Metapop_Data$I, col = i)
  
  
}
legend(400,0.2, legend = c('dI = 0.1','dI =0.3','dI =0.7'), col = c(1,2,3), lty = 1:2, cex = 0.8)





#### TRAITS VALUES ANALYSIS #### -----------------------------------------------------------------------------------------------------

#Mean Trait dynamics for average trait value

i=1
#ActualDf = list_df_Distributions[[1]]
ActualDf = DistribFiles[[i]]
Alp = seq(0.01, 0.49, 0.01) 
FinalDf <- subset(ActualDf[c(-1)])
Alfinal = NULL
Moyalpha = NULL


for (i in 1:nrow(FinalDf)){
  Alfinal = NULL
  for (j in 1:length(Alp)){
    Alfi = (FinalDf[i,j]*Alp[j])
    Alfinal = c(Alfinal, Alfi)
  }
  Moyalf = sum(Alfinal)/sum(FinalDf[i,])
  Moyalpha = c(Moyalpha,Moyalf)
}
i=1
plot(ActualDf$Time, Moyalpha, type = 'l',xlab = 'time', ylab = 'Mean Trait Value Alpha', ylim = c(0,0.25))
#plot(list_df_Distributions[[4]]$Time, Moyalpha, type = 'l',xlab = 'time', ylab = 'Mean Trait Value', ylim = c(0,0.5))
legend("bottomleft",legend = c("0.1","0.2", "0.3", "0.4","0.5", "0.6", "0.7", "0.8", "0.9", "1.0"), col = c(1,2,3,4,5,6,7,8,9,10), lty = 1)


for (i in 2:length(DistribFiles)){
  ActualDf = DistribFiles[[i]]
  Alp = seq(0.01, 0.49, 0.01) 
  FinalDf <- subset(ActualDf[c(-1)])
  Alfinal = NULL
  Moyalpha = NULL
  
  for (j in 1:nrow(FinalDf)){
    Alfinal = NULL
    for (k in 1:length(Alp)){
      Alfi = (FinalDf[j,k]*Alp[k])
      Alfinal = c(Alfinal, Alfi)
    }
    Moyalf = sum(Alfinal)/sum(FinalDf[j,])
    Moyalpha = c(Moyalpha,Moyalf)
  }
  lines(ActualDf$Time, Moyalpha, type = 'l', col = i)
}



# Converge values for evolution of the trait virulence (alpha)
taillemax <- length(Moyalpha)
ConvergeValue <- Moyalpha[(taillemax-1000):taillemax]
ConvergeValues <- sum(ConvergeValue)/length(ConvergeValue)
ConvergeValues





###### ANALYSIS OF THE TRAIT DISTRIBUTIONS FOR DISCRETE ALPHA VALUES ########### -------------------------------------------------------------------


#ActualDf = DistribFiles[[1]]
ActualDf = list_df_Distributions[[6]]
ActualDf = ActualDf[ActualDf$Time > 50,]
todrop <- c('Time','t')
df2 <- ActualDf[,!(names(ActualDf) %in% todrop)]
FinalDf <- subset(df2[c(-1)])
time <- as.vector(ActualDf$Time)
alpha <- seq(0.0,1.45,0.05)
Densities <- as.matrix(subset(FinalDf[,-1]))
colnames(Densities) <- NULL
  
fig <- plot_ly( x = alpha, y = time,z = ~Densities, type = "heatmap")%>%layout(title = 'Phenotypic Distribution dynamics, dS=0.5, dI=0.5', plot_bgcolor = "#e5ecf6", xaxis = list(title = 'Alpha'), 
                                                                                 yaxis = list( title = 'Time'),zaxis = list(range = c(0,6000),title = "Density"), legend = list(title=list(text='<b> What is that ? </b>')))

plotly_build(fig)



################# d = dS = dI VS  dS fixé et dI varie ############

dSI = c(0.08334129,0.09719584,0.11136409,0.12640946,0.14388560,0.14368753,0.14651693,0.15550300,0.21072297,0.28627977)

dSfixI <- c(0.09236266,0.10008508,0.10987805,0.12717639,0.14388560,0.14166545,0.15709891,0.17053767,0.16904752,0.18559184)
dSIfix = c(0.14283853,0.13503640,0.14067641,0.14135839,0.14388560,0.14129463,0.13744057,0.14167341,0.14023871,0.18072271)
Alpvalues = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

plot(dSI,dSIfix, col = "red",pch = 16, ylim = c(0.07,0.3),xlim=c(0.07,0.3), xlab = "alpha_dS=dI", ylab = "alpha_dIfix_S", main = "Comparaison entre valeur de virulence 2")
abline(coef = c(0,1))






######################################### Partie un peu en bordel où plein de packages ont été testés pour tenter d'étudier les branchements évolutifs ###########

######## Gaussian Mixture Models - Test ##

ActualDf1 = DistribFiles[[1]]
ActualDf1 = ActualDf1[,-1] # data (excluding the response variable)
ActualDf2 = list_df_Distributions[[6]]
ActualDf2 = ActualDf2[,-1]
Alp = seq(0.01, 0.49, 0.01) 

Tr = ActualDf2[8400,]
Tr <- round(Tr,0)
FrqAlp = NULL
for (i in 1:length(Alp)){
  Frq = rep(Alp[i],Tr[i])
  FrqAlp = c(FrqAlp,Frq)
}
# for (i in 1:length(Tr)){
#   FrqAlp = c(FrqAlp,Tr[i])
# }
# FrqAlp = log(FrqAlp)
# hist(FrqAlp, breaks = 5)

FrqAlp <- as.data.frame(FrqAlp)
FrqAlp$FrqAlp2 = FrqAlp
mydi1 = dist(FrqAlp)
FrqAlp <- log(FrqAlp)
FrqAlp <- t(FrqAlp)

Tr <- as.numeric(Tr[1,])
Tr <- as.data.frame(Tr)
Tr$Tr2 = Tr
mydi2 = dist(Tr)


myclust1 <- hclust(mydi1, method="ward.D2")
plot(myclust1)
myclust2 <- hclust(mydi2, method="ward.D2")
plot(myclust2)


inertie <- sort(myclust1$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie",lwd=2);grid()
k <- 2
abline(v=k,col="red",lty=3)
points(k,inertie[k],pch=16,cex=2,col="red")

inertie <- sort(myclust2$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie",lwd=2);grid()
k <- 2
abline(v=k,col="red",lty=3)
points(k,inertie[k],pch=16,cex=2,col="red")



### Package ClusterR
library(ClusterR)
dat = center_scale(FrqAlp, mean_center = T, sd_scale = T)  # centering and scaling the data
gmm = GMM(dat, 1, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
          em_iter = 10, verbose = F)    
pr = predict(gmm, newdata = dat)

opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 10 , criterion = "AIC",
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               km_iter = 0, em_iter = 100, var_floor = 1,
                               plot_data = T)

#### Normal law vs Gaussian - Example
library(ClusterR)
x = rnorm(240, mean = 0.14, sd = 0.01)
x<-round(x, 2)
Normlaw = as.data.frame(x)
nrow(unique(Normlaw))
dat = center_scale(Normlaw, mean_center = T, sd_scale = T) 
opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 40, criterion = "BIC", 
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               km_iter = 150, em_iter = 150, var_floor = 1e-10, 
                               plot_data = T)
y = rnorm(150, mean = 0, sd = 1)
hist(x, breaks = 5)
hist(y)

Normlaw = c(x,FrqAlp)
hist(Normlaw)
Normlaw = as.data.frame(Normlaw)
dat = center_scale(Normlaw, mean_center = T, sd_scale = T)  # centering and scaling the data
gmm = GMM(dat, 2, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
          em_iter = 10, verbose = F) 

pr = predict(gmm, newdata = dat)
opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 5, criterion = "BIC", 
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               km_iter = 10, em_iter = 10, var_floor = 1e-10, 
                               plot_data = T)





######################### Mixtools #########################################################

Mydf <- list_df_Distributions[[1]]
#Mydf <- DistribFiles[[1]]
alphas <- colnames(Mydf[2:length(Mydf)])
b <- as.numeric(gsub("Alpha", "", alphas))
b
MySubDf <- Mydf[-c(1:40),-c(1)]
Finaldata <- data.matrix(MySubDf, rownames.force = NA) 
fig <- plot_ly(x = b, y = Mydf$Time,z = Finaldata, type = "heatmap")
fig


# Try for Multinomial mixture model using Mixtools
#Select row
Mydistrib <- Finaldata[666,]
Mydistrib
Alp = seq(0.01, 0.49, 0.01) 
FrqAlp = NULL
for (i in 1:length(Alp)){
  Frq = rep(Alp[i],Mydistrib[i])
  FrqAlp = c(FrqAlp,Frq)
}
hist(FrqAlp, breaks=6)
FrqAlp <- as.data.frame(FrqAlp)

wait1 <- normalmixEM(FrqAlp)
wait1$loglik
plot(wait1, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8, xlab2="Minutes")


#Create the distribution
x<-rep(c(0.18,0.19,0.20,0.21,0.22, 0.33,0.34,0.35,0.36,0.37),times=c(13,106,237,160,1,13,106,237,160,1))
hist(x)
wait1 <- normalmixEM(x)
wait1$lambda
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



############################# NbCLust ############################################


set.seed(123)


Mydistrib <- Finaldata[4000,]
Mydistrib
Alp = seq(0.001, 0.499, 0.001) 
FrqAlp = NULL
for (i in 1:length(Alp)){
  Frq = rep(Alp[i],Mydistrib[i])
  FrqAlp = c(FrqAlp,Frq)
}
#hist(FrqAlp)
FrqAlpdf <- as.data.frame(FrqAlp)

### Optimal number of clusters in the data
res = fviz_nbclust(FrqAlpdf, kmeans, method = "gap_stat", k.max=10, nboot = 100, nstart = 25)
fviz_nbclust(FrqAlpdf, kmeans, method = "gap_stat", k.max=10, nboot = 100, nstart = 25) +
  geom_vline(xintercept = 3, linetype = 2)


# Remove species column (5) and scale the data
FrqAlp.scaled <- scale(FrqAlp)

### Gap statistic

# Compute gap statistic for kmeans
# we used B = 10 for demo. Recommended value is ~500
gap_stat <- clusGap(FrqAlp.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 100)
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)
Tab_gap = gap_stat$Tab[,3]
Max_clust = which.max(Tab_gap)


GapData = Finaldata[4000:5000,]
GapData_Clust = GapData[seq(1, nrow(GapData), 10), ]

Max_cluster = NULL
for (i in 1:nrow(GapData_Clust)){
  cat(i,"/", nrow(GapData_Clust))
  Mydistrib <- GapData_Clust[i,]
  Alp = seq(0.001, 0.499, 0.001) 
  FrqAlp = NULL
  for (i in 1:length(Alp)){
    Frq = rep(Alp[i],Mydistrib[i])
    FrqAlp = c(FrqAlp,Frq)
  }
  FrqAlpdf <- as.data.frame(FrqAlp)
  FrqAlp.scaled <- scale(FrqAlp)
  gap_stat <- clusGap(FrqAlp.scaled, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 100)
  Tab_gap = gap_stat$Tab[,3]
  Max_clust = which.max(Tab_gap)
  Max_cluster = c(Max_cluster, Max_clust)
}
Max_cluster
plot(Max_cluster,type='l')




######## DBSCAN ################

library(dbscan)

ActualDf1 = DistribFiles[[1]]
ActualDf1 = ActualDf1[,-1] # data (excluding the response variable)
ActualDf2 = list_df_Distributions[[1]]
ActualDf2 = ActualDf2[,-1]
Alp = seq(0.01, 0.49, 0.01) 

Tr = ActualDf2[3333,]
Tr <- round(Tr,0)
FrqAlp = NULL
for (i in 1:length(Alp)){
  Frq = rep(Alp[i],Tr[i])
  FrqAlp = c(FrqAlp,Frq)
}
# for (i in 1:length(Tr)){
#   FrqAlp = c(FrqAlp,Tr[i])
# }
# FrqAlp = log(FrqAlp)
# hist(FrqAlp, breaks = 5)
FrqAlp <- as.data.frame(FrqAlp)
FrqAlp <- log(FrqAlp)
FrqAlp <- t(FrqAlp)


res <- dbscan(FrqAlp, eps = .3, minPts = 5)
res
?dbscan
## plot clusters and add noise (cluster 0) as crosses.
plot(FrqAlp, col = res$cluster)
points(FrqAlp[res$cluster == 1, ], pch = 3, col = "grey")
hullplot(FrqAlp, res)

hullplot(FrqAlp, res)
points(FrqAlp, pch = 3 , col = "red", lwd = 3)
text(newdata, pos = 1)

pred_label <- predict(res, FrqAlp, data = FrqAlp)
pred_label
points(newdata, col = pred_label + 1L,  cex = 2, lwd = 2)


