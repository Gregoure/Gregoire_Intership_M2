##### Rscript for contour plot and ESS 
##### but this time, Evolution of Virulence and infected dispersal 

# Clean and clear
rm(list = ls())
graphics.off()

# Library used
library(readr)
library(dplyr)
library(plot3D)
library(plotly)
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(tibble)
library(ggpubr)


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


### GET FILES ############

dir <- choose.dir() 

#### Build a list of folder DENSITIES, TRAITS alpha and dispersal, DISTRIBUTION alpha and dispersal ####

cwd <- setwd(dir)
cwd
listfiles <- list.files(cwd)
listfiles

DensitiesFiles <- list()
TraitsAlphaFiles <- list()
TraitsDIFiles <- list()
DistribAlphaFiles <- list()
DistribDIFiles <- list()
DistribAlphaDIFiles <- list()


for (j in 1:length(listfiles)){
  currentdir <- paste(dir, listfiles[j], sep = "\\") #Browse subfolders in main folders
  DensitiesDir <- paste(currentdir, "Densities", sep = "\\") # Get the subsubfolder with densities files
  DistribAlphaDir <- paste(currentdir, "Distributions_alpha", sep = "\\") # Get the subfolder for distribution of Alpha traits
  DistribDIDir <- paste(currentdir, "Distributions_dI", sep = "\\") # Get the subfolder for distribution of Dispersal I traits
  TraitsAlphaDir <- paste(currentdir, "Traits_alpha", sep = "\\") # Get the subsubfolder with Alpha traits files
  TraitsDIDir <- paste(currentdir, "Traits_dI", sep = "\\") # Get the subsubfolder with Dispersal I traits files
  
  
  ########## Alpha Traits #############################################################################################
  
  cwd <- setwd(TraitsAlphaDir)
  cwd
  
  ls_files_Traits_alpha <- list.files(TraitsAlphaDir) # list files
  list_tibbles_traits_alpha <- lapply(ls_files_Traits_alpha, read_csv) # Read CSV returns tibble object
  list_df_traits_alpha <- lapply(list_tibbles_traits_alpha, as.data.frame) # Turn them to good old dataframes
  
  ### I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_traits_alpha)){
    list_df_traits_alpha[[i]]<- subset(list_df_traits_alpha[[i]][c(-1)])
    list_df_traits_alpha[[i]]$Time <- round(list_df_traits_alpha[[i]]$Time, 1)
    list_df_traits_alpha[[i]] <- distinct(list_df_traits_alpha[[i]], Time, .keep_all = TRUE)
  }
  
  # Get only commons rows between all dataframes 
  # Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_traits_alpha)){
    times <- list_df_traits_alpha[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  
  # Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_traits_alpha)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  # loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_traits_alpha)){
    indexrows <- which(list_df_traits_alpha[[i]]$Time %in% tcommon)
    list_df_traits_alpha[[i]] <- subset(list_df_traits_alpha[[i]][indexrows,])
  }
  
  
  ### II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_traits_alpha[[1]])
  MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_traits_alpha[[1]]$Time)))
  colnames(MeanSim) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_traits_alpha,namecols[i])
    MeanSim[,i] <- meanvar
  }
  write.csv2(MeanSim, paste0("C:\\Users\\gaze\\Downloads\\TraitsAlphaFiles",j,".csv"),row.names = FALSE)
  # And add to the list
  TraitsAlphaFiles <- c(TraitsAlphaFiles, list(MeanSim))
  
  
  
  
  ########## Dispersal I Traits #############################################################################################
  
  cwd <- setwd(TraitsDIDir)
  cwd
  
  ls_files_Traits_dI <- list.files(TraitsDIDir) # list files
  list_tibbles_traits_dI <- lapply(ls_files_Traits_dI, read_csv) # Read CSV returns tibble object
  list_df_traits_dI <- lapply(list_tibbles_traits_dI, as.data.frame) # Turn them to good old dataframes
  
  ### I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_traits_dI)){
    list_df_traits_dI[[i]]<- subset(list_df_traits_dI[[i]][c(-1)])
    list_df_traits_dI[[i]]$Time <- round(list_df_traits_dI[[i]]$Time, 1)
    list_df_traits_dI[[i]] <- distinct(list_df_traits_dI[[i]], Time, .keep_all = TRUE)
  }
  
  # Get only commons rows between all dataframes 
  # Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_traits_dI)){
    times <- list_df_traits_dI[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  
  # Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_traits_dI)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  # loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_traits_dI)){
    indexrows <- which(list_df_traits_dI[[i]]$Time %in% tcommon)
    list_df_traits_dI[[i]] <- subset(list_df_traits_dI[[i]][indexrows,])
  }
  
  
  ### II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_traits_dI[[1]])
  MeanSim2 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_traits_dI[[1]]$Time)))
  colnames(MeanSim2) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_traits_dI,namecols[i])
    MeanSim2[,i] <- meanvar
  }
  write.csv2(MeanSim2, paste0("C:\\Users\\gaze\\Downloads\\TraitsDIFiles",j,".csv"),row.names = FALSE)
  TraitsDIFiles <- c(TraitsDIFiles, list(MeanSim2))
  
  
  
  ########## Densities #############################################################################################
  
  cwd <- setwd(DensitiesDir)
  
  ls_files_Densities <- list.files(DensitiesDir) #list files
  list_tibbles_densities <- lapply(ls_files_Densities, read_csv) 
  list_df_densities <- lapply(list_tibbles_densities, as.data.frame) 
  
  ### I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_densities)){
    list_df_densities[[i]]<- subset(list_df_densities[[i]][c(-1)]) # suppress first column
    list_df_densities[[i]]$Time <- round(list_df_densities[[i]]$Time, 1) #round t values
    list_df_densities[[i]] <- distinct(list_df_densities[[i]], Time, .keep_all = TRUE)#remove duplicates
  }
  
  # Get only commons rows between all dataframes 
  # Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_densities)){
    times <- list_df_densities[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1]
  counts <- as.data.frame(table(t_values))

  # Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_densities)
  tcommon<- counts$t_values[counts$Freq==nbsim]
  
  # loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_densities)){
    indexrows <- which(list_df_densities[[i]]$Time %in% tcommon)
    list_df_densities[[i]] <- subset(list_df_densities[[i]][indexrows,])
  }
  
  
  ### II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_densities[[1]])
  MeanSim3 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_densities[[1]]$Time)))
  colnames(MeanSim3) <- namecols
  
  for (i in 1:length(namecols)){
    meanvar <- MeanSimValue(list_df_densities,namecols[i])
    MeanSim3[,i] <- meanvar
  }
  write.csv2(MeanSim3, paste0("C:\\Users\\gaze\\Downloads\\DensitiesFiles",j,".csv"),row.names = FALSE)
  DensitiesFiles <- c(DensitiesFiles, list(MeanSim3))
  
  
  
  ########## Alpha Distribution #############################################################################################
  
  cwd <- setwd(DistribAlphaDir)
  
  ls_files_Distribtions_alpha <- list.files(DistribAlphaDir) #list files
  list_tibbles_Distributions_alpha <- lapply(ls_files_Distribtions_alpha, read_csv) 
  list_df_Distributions_alpha <- lapply(list_tibbles_Distributions_alpha, as.data.frame) 
  
  ### I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_Distributions_alpha)){
    
    list_df_Distributions_alpha[[i]]<- subset(list_df_Distributions_alpha[[i]][c(-1)]) # suppress first column
    list_df_Distributions_alpha[[i]]$Time <- round(list_df_Distributions_alpha[[i]]$Time, 1) #round t values
    list_df_Distributions_alpha[[i]] <- distinct(list_df_Distributions_alpha[[i]], Time, .keep_all = TRUE) # remove duplicates
  }
  
  # Get only commons rows between all dataframes 
  # Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_Distributions_alpha)){
    times <- list_df_Distributions_alpha[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1]
  counts <- as.data.frame(table(t_values))

  # Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_Distributions_alpha)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  # loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_Distributions_alpha)){
    indexrows <- which(list_df_Distributions_alpha[[i]]$Time %in% tcommon)
    list_df_Distributions_alpha[[i]] <- subset(list_df_Distributions_alpha[[i]][indexrows,])
  }
  
  ### II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_Distributions_alpha[[1]])
  MeanSim4 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_Distributions_alpha[[1]]$Time)))
  colnames(MeanSim4) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_Distributions_alpha,namecols[i])
    MeanSim4[,i] <- meanvar
  }
  
  write.csv2(MeanSim4, paste0("C:\\Users\\gaze\\Downloads\\DistribAlphaFiles",j,".csv"),row.names = FALSE)
  
  DistribAlphaFiles <- c(DistribAlphaFiles, list(MeanSim4))
  
  
  ########## Dispersal I Distribution #############################################################################################
  
  cwd <- setwd(DistribDIDir)
  
  ls_files_Distribtions_dI <- list.files(DistribDIDir) #list files
  list_tibbles_Distributions_dI <- lapply(ls_files_Distribtions_dI, read_csv) 
  list_df_Distributions_dI <- lapply(list_tibbles_Distributions_dI, as.data.frame) 
  
  ### I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_Distributions_dI)){
    
    list_df_Distributions_dI[[i]]<- subset(list_df_Distributions_dI[[i]][c(-1)])
    list_df_Distributions_dI[[i]]$Time <- round(list_df_Distributions_dI[[i]]$Time, 1)
    list_df_Distributions_dI[[i]] <- distinct(list_df_Distributions_dI[[i]], Time, .keep_all = TRUE)
  }
  
  # Get only commons rows between all dataframes 
  # Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_Distributions_dI)){
    times <- list_df_Distributions_dI[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1]
  counts <- as.data.frame(table(t_values))
  
  # Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_Distributions_dI)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  # loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_Distributions_dI)){
    indexrows <- which(list_df_Distributions_dI[[i]]$Time %in% tcommon)
    list_df_Distributions_dI[[i]] <- subset(list_df_Distributions_dI[[i]][indexrows,])
  }
  
  ### II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_Distributions_dI[[1]])
  MeanSim5 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_Distributions_dI[[1]]$Time)))
  colnames(MeanSim5) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_Distributions_dI,namecols[i])
    MeanSim5[,i] <- meanvar
  }
  
  write.csv2(MeanSim5, paste0("C:\\Users\\gaze\\Downloads\\DistribDIFiles",j,".csv"),row.names = FALSE)
  
  DistribDIFiles <- c(DistribDIFiles, list(MeanSim5))
  
  
  ########## Two Traits Distribution #############################################################################################
  
  DistribAlphaDIDir <- paste(currentdir, "Distributions_AlpdI", sep = "\\")
  
  cwd <- setwd(DistribAlphaDIDir)
  
  ls_files_Distribtions_AlphadI <- list.files(DistribAlphaDIDir) #list files
  list_tibbles_Distributions_AlphadI <- lapply(ls_files_Distribtions_AlphadI, read_csv) 
  list_df_Distributions_AlphadI <- lapply(list_tibbles_Distributions_AlphadI, as.data.frame) 
  
  ### I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_Distributions_AlphadI)){
    
    list_df_Distributions_AlphadI[[i]]<- subset(list_df_Distributions_AlphadI[[i]][c(-1)])
    list_df_Distributions_AlphadI[[i]]$Time <- round(list_df_Distributions_AlphadI[[i]]$Time, 1)
    list_df_Distributions_AlphadI[[i]] <- distinct(list_df_Distributions_AlphadI[[i]], Time, .keep_all = TRUE)
  }
  
  # Get only commons rows between all dataframes 
  # Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_Distributions_AlphadI)){
    times <- list_df_Distributions_AlphadI[[i]]$Time
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1]
  counts <- as.data.frame(table(t_values))
  
  # Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_Distributions_AlphadI)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  # loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_Distributions_AlphadI)){
    indexrows <- which(list_df_Distributions_AlphadI[[i]]$Time %in% tcommon)
    list_df_Distributions_AlphadI[[i]] <- subset(list_df_Distributions_AlphadI[[i]][indexrows,])
  }
  
  ### II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_Distributions_AlphadI[[1]])
  MeanSim6 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_Distributions_AlphadI[[1]]$Time)))
  colnames(MeanSim6) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_Distributions_AlphadI,namecols[i])
    MeanSim6[,i] <- meanvar
  }
  
  write.csv2(MeanSim6, paste0("C:\\Users\\gaze\\Downloads\\DistribAlphaDIFiles",j,".csv"),row.names = FALSE)
  
  DistribAlphaDIFiles <- c(DistribAlphaDIFiles, list(MeanSim6))

  
} # GIANT LOOP THAT EXTRACT, COMPILE AND BUILD A MEAN SIM WITH ALL FILES


### DEMOGRAPHY ANALYSIS ### ----------------------------------------------------------------------------------------------------------
ActualDf <- DensitiesFiles[[i]]
taille_col = (length(ActualDf)-1)
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2)
generate_odd_indexes = seq(3, taille_col+1, 2)

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean

Metapop_Data <- data.frame(t = ActualDf$Time, S = Smean, I=Imean, N = Nmean)

###### Phase Space #####
plot(Metapop_Data$S, Metapop_Data$I, type = 'l', xlab = 'Susceptibles', ylab = 'Infected')






#### TRAITS VALUES ANALYSIS #### -----------------------------------------------------------------------------------------------------

#Mean Trait dynamics for average trait value

i=1 # Choice of file, here Value of dS and/or epsilon fixed

ActualDf1 = DistribAlphaFiles[[i]] # For Virulence Value
ActualDf2 = DistribDIFiles[[i]] # For Lambda Value

Alp = seq(0.01, 0.49, 0.01) # Possible values of Virulence
dI = seq(0.1,5.0,0.1) # Possible values of dispersal
FinalDf1 <- subset(ActualDf1[c(-1)]) # Virulence
FinalDf2 <- subset(ActualDf2[c(-1)]) # Lambda
Alfinal = NULL
Moyalpha = NULL
dIfinal = NULL
MoydI = NULL

# Virulence
for (i in 1:nrow(FinalDf1)){
  Alfinal = NULL
  for (j in 1:length(Alp)){
    Alfi = (FinalDf1[i,j]*Alp[j])
    Alfinal = c(Alfinal, Alfi)
  }
  Moyalf = sum(Alfinal)/sum(FinalDf1[i,])
  Moyalpha = c(Moyalpha,Moyalf)
}
# Lambda
for (i in 1:nrow(FinalDf2)){
  dIfinal = NULL
  for (j in 1:length(dI)){
    dIfi = (FinalDf2[i,j]*dI[j])
    dIfinal = c(dIfinal, dIfi)
  }
  MoyI = sum(dIfinal)/sum(FinalDf2[i,])
  MoydI = c(MoydI,MoyI)
}


plot(ActualDf1$Time, Moyalpha, type = 'l',xlab = 'time', ylab = 'Mean Trait Value alpha', ylim = c(0,0.25))
plot(ActualDf2$Time, MoydI, type = 'l',xlab = 'time', ylab = 'Mean Trait Value dI')



### Heatmap Part #### -------------------------------------------------------------------------------------------------------------
Alp = seq(0.01, 0.49, 0.01)
dI = seq(0.1,4.9,0.1)
tmax = 1500
Possible_Timesteps = c(tmax-50,tmax-100,tmax-200)


### Virulence Data
NBROW = 11
HeatFinalVir = data.frame(matrix(NA, nrow = NBROW, ncol = 0))
for (j in 1:length(Possible_Timesteps)){
  cat("TimeSteps", j, "\n")
  HeatValues = NULL
  Step = 1
  for (i in 1:length(DistribAlphaFiles)){
    cat(Step, "/", length(DistribAlphaFiles), "\n")
    ActualDf <- DistribAlphaFiles[[i]]
    ActualDf <- ActualDf[ActualDf$Time > Possible_Timesteps[j],]
    FinalDf <- subset(ActualDf[c(-1)])
    ConvergeValue <- NULL
    
    for (k in 1:nrow(FinalDf)){
      Alfinal = NULL
      for (l in 1:length(Alp)){
        Alfi <- (FinalDf[k,l]*Alp[l])
        Alfinal <- c(Alfinal, Alfi)
      }
      Moyalf <- sum(Alfinal)/sum(FinalDf[k,])
      ConvergeValue <- c(ConvergeValue,Moyalf)
    }
    ConvergeValues <- sum(ConvergeValue)/length(ConvergeValue)
    HeatValues <- c(HeatValues,ConvergeValues)
    Step = Step + 1
  }
  HeatFinalVir[j]=HeatValues
}


### Infected dispersal data
NBROW = 11
HeatFinalDis = data.frame(matrix(NA, nrow = NBROW, ncol = 0))
for (j in 1:length(Possible_Timesteps)){
  cat("TimeSteps", j, "\n")
  HeatValues = NULL
  Step = 1
  for (i in 1:length(DistribDIFiles)){
    cat(Step, "/", length(DistribDIFiles), "\n")
    ActualDf <- DistribDIFiles[[i]]
    ActualDf <- ActualDf[ActualDf$Time > Possible_Timesteps[j],]
    FinalDf <- subset(ActualDf[c(-1)])
    ConvergeValue <- NULL
    
    for (k in 1:nrow(FinalDf)){
      DIfinal = NULL
      for (l in 1:length(dI)){
        DIfi <- (FinalDf[k,l]*dI[l])
        DIfinal <- c(DIfinal, DIfi)
      }
      MoyDI <- sum(DIfinal)/sum(FinalDf[k,])
      ConvergeValue <- c(ConvergeValue,MoyDI)
    }
    ConvergeValues <- sum(ConvergeValue)/length(ConvergeValue)
    HeatValues <- c(HeatValues,ConvergeValues)
    Step = Step + 1
  }
  HeatFinalDis[j]=HeatValues
}

# Values of dS and epsilon fixed for sumulation 
dS = c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0,1.1,1.1,1.1)
eps = c(0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2)
HeatFinalDis$dS = dS
HeatFinalVir$dS = dS
HeatFinalDis$eps = eps
HeatFinalVir$eps = eps



## Plot Ds_lambda and dS_Virulence

HeatAnalysis = HeatFinalVir # Or HeatFinalDis
permutation <- order(HeatAnalysis$dS)
HeatAnalysis[permutation,]

plot(HeatAnalysis[permutation,]$dS,HeatAnalysis[permutation,]$V3, xlab = "", ylab = "",ylim = c(0,0.25),
     pch = 16, col = "red")
mtext(expression("Taux de dispersion des individus susceptibles (t"^-1*")"), side=1, line=2.5, cex=1.5)
mtext(expression("Virulence du parasite (t"^-1*")"), side=2, line=2.3, cex=1.5)




## Plot Heatmap_Dis and Heatmap_Virulence

Matrixdf = data.frame(matrix(nrow = length(HeatFinalVir$dS) , ncol = 0))
Matrixdf$dS = as.factor(HeatFinalVir$dS)
Matrixdf$eps = as.factor(HeatFinalVir$eps)

## Convergence Values for Virulence
Matrixdf$HeatValues = as.numeric(HeatFinalVir[,3])
Matrix = xtabs(HeatValues ~dS + eps, Matrixdf)

## Convergence Values for Infected Dispersal
Matrixdf$HeatValues = as.numeric(HeatFinalDis[,3])
Matrix = xtabs(HeatValues ~dS + eps, Matrixdf)

Matrix[Matrix == 0] <- NA
Matridf = as.data.frame(Matrix)


p <- ggplot(Matridf, aes(eps, dS, fill= Freq)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  theme_ipsum() +  
  ggtitle("") + 
  xlab("Epsilon") + ylab("mS") + labs(fill = "Dis") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 15))


ggplotly(p, tooltip="text")





### Analysis of the two traits, Virulence and dispersal ### ------------------------------------------------------------------------------------------------------------------------------------

# Two vectors of value of parameters for Virulence and Dispersal
dS <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1)
Files <- seq(2,11,1)
Alp <- seq(0.01, 0.49, 0.01)
Time = c(0,1,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400)
dI <- seq(0.1,4.9,0.1)
AlpD <-rep(Alp,each=length(dI))
AlpD2 <- rep(AlpD, times = length(Time))
dID <- rep(dI, times=length(Alp)) 
dID2 <- rep(dID, times=length(Time))
TimeD2 = rep(Time, each = length(AlpD))

DistribAlphaDI = NULL
for (i in 1:nrow(DistribAlphaDIFiles[[1]])){
  Distrib <- c(t(DistribAlphaDIFiles[[1]][i,]))
  Distrib <- Distrib[-1]
  DistribAlphaDI = c(DistribAlphaDI,Distrib)
}
DistribAlphaDI <- as.data.frame(DistribAlphaDI)
DistribAlphaDI$Alpha = AlpD2
DistribAlphaDI$mI = dID2
DistribAlphaDI$Time = TimeD2

AlphaDI <- DistribAlphaDI[DistribAlphaDI$DistribAlphaDI != 0, ] 
AlphaDI <- AlphaDI[AlphaDI$Time == 1400,]
dSAlphaDI <- rep(0.1,times = nrow(AlphaDI))
AlphaDI$dS <- dSAlphaDI
# AlphaDI <- AlphaDI[AlphaDI$DistribAlphaDI > 0.01*sum(AlphaDI$DistribAlphaDI),]
AlphaDI <- AlphaDI[AlphaDI$DistribAlphaDI > 5,]
FinalDf = AlphaDI

for (j in Files){
  cat(j,"\n")
  DistribAlphaDI = NULL
  for (i in 1:nrow(DistribAlphaDIFiles[[j]])){
    Distrib <- c(t(DistribAlphaDIFiles[[j]][i,]))
    Distrib <- Distrib[-1]
    DistribAlphaDI = c(DistribAlphaDI,Distrib)
  }
  DistribAlphaDI <- as.data.frame(DistribAlphaDI)
  DistribAlphaDI$Alpha = AlpD2
  DistribAlphaDI$mI = dID2
  DistribAlphaDI$Time = TimeD2
  
  AlphaDI <- DistribAlphaDI[DistribAlphaDI$DistribAlphaDI != 0, ] 
  AlphaDI <- AlphaDI[AlphaDI$Time == 1400,]
  dSAlphaDI <- rep(dS[j],times = nrow(AlphaDI))
  AlphaDI$dS <- dSAlphaDI
  sum(AlphaDI$DistribAlphaDI)
  #AlphaDI <- AlphaDI[AlphaDI$DistribAlphaDI > 0.05*sum(AlphaDI$DistribAlphaDI),]
  AlphaDI <- AlphaDI[AlphaDI$DistribAlphaDI > 5,]
  FinalDf <- rbind(FinalDf,AlphaDI)
}
FinalDf





GraphData <- FinalDf[FinalDf$dS == 0.5,]
write.csv(GraphData, "C:\\Users\\gaze\\Downloads\\GraphData.csv", row.names=FALSE)
GraphData = read.csv("GraphData.csv",header = TRUE)


GraphData$Nb = GraphData$DistribAlphaDI
p <- ggplot(GraphData, aes(Alpha, mI)) + geom_point(aes(size = Nb)) +
  xlim(0,1.3) + ylim(0,1.3) + xlab("Virulence du parasite") + ylab("Influence parasitaire sur la dispersion des infectés") + labs(fill = "Nb") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 20), legend.text = element_text(size = 15))
ggplotly(p, tooltip="text")





# Correlation between dispersal and virulence

dS <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1)
Correl = NULL
for (i in dS){
  GraphData <- FinalDf[FinalDf$dS == i,]
  Shapiro1 = shapiro.test(GraphData$Alpha)
  Shapiro2 = shapiro.test(GraphData$mI)
  res <- cor.test(GraphData$mI, GraphData$Alpha, 
                  method = "pearson")
  cat("Valeur de dS = ",i,"\n","Test de shapiro", Shapiro1$p.value, Shapiro2$p.value,"\n", "Résultat corrélation =",res$p.value,"\n")
  Resultats <- res$estimate
  Correl <- c(Correl,Resultats)
}

plot(dS,Correl, xlim = c(0,1.1), ylim = c(-0.2, 0.45), pch = 16, col = "red", xlab = "", ylab = "")
mtext(expression("Taux de dispersion des susceptibles (t"^-1*")"), side=1, line=2.5, cex=1.5)
mtext(expression("Corrélation entre Virulence et dispersion des infectés"), side=2, line=2.3, cex=1.5)
abline(h=0)



cor(GraphData[,c('DistribAlphaDI', 'Alpha', 'mI')])
cor(GraphData$Alpha,GraphData$mI, method = "spearman")

shapiro.test(GraphData$Alpha)
shapiro.test(GraphData$mI)

res <- cor.test(GraphData$mI, GraphData$Alpha, 
                method = "pearson")
res

ggscatter(GraphData, x = "Alpha", y = "mI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Alpha", ylab = "mI")




