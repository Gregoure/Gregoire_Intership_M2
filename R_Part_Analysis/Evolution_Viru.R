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
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(tibble)
library(tidyverse)


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

#######################################################
##### Script for analysis evolution of Virulence ######
#######################################################


### GET FILES - Results of simulation ###
dir <- choose.dir() # Has to go to global folder in which all sims are stored in subfolders

cwd <- setwd(dir) 
cwd
listfiles <- list.files(cwd)
listfiles

DistribFiles <- list()
DensitiesFiles <- list()
TraitsFiles <- list()

tmax = 1000 # tmax utilisé dans les simu à compiler
Messages= NULL # Stockage des messages d'erreur, pour avoir toutes les simu qui se sont cassés la gueule
list_paramsS = NULL # Stockage des valeurs de dispersion susceptibles (mS) de chaque fichier
list_paramsI = NULL # Stockage des valeurs de dispersion infectés (mI) de chaque fichier

for (i in 1:length(listfiles)){
  cat("Dossier traité N°", i,"/", length(listfiles), "\n") # Check de où on est dans les dossiers traités
  currentdir <- paste(dir, listfiles[i], sep = "\\") # Browse subfolders in main folders
  
  DistribDir <- paste(currentdir, "Distributions", sep = "\\") # Get the subfolder for distribution of traits when discrete trait values are used
  DensitiesDir <- paste(currentdir, "Densities", sep = "\\") # Get the subsubfolder with densities files
  
  
  ##### DISTRIBUTIONS FILES #############
  
  cwd <- setwd(DistribDir)
  
  ls_files_Distributions <- list.files(DistribDir) #list files
  list_tibbles_Distributions <- lapply(ls_files_Distributions, read_csv) 
  list_df_Distributions <- lapply(list_tibbles_Distributions, as.data.frame) 
  
  ## I. Shaping dataframes and getting common rows
  remov = NULL  # Stockage des seed qui n'ont pas atteint tmax
  for (i in 1:length(list_df_Distributions)){
    if (max(list_df_Distributions[[i]][,2]) > tmax - 1){ # Garde les seeds qui ont un temps final > tmax - 1
      list_df_Distributions[[i]]<- subset(list_df_Distributions[[i]][c(-1)]) # suppress first column
      list_df_Distributions[[i]]$Time <- round(list_df_Distributions[[i]]$Time, 1) #round t values
      list_df_Distributions[[i]] <- distinct(list_df_Distributions[[i]], Time, .keep_all = TRUE) # remove duplicates
    }else{
      remov <- c(remov,i) # Seed qui ont un temps final < tmax-1
    }
  }
  # Recuperation des valeurs de mS et mI correspondant au dossier traite, directement dans le nom du fichier
  # A modifier si changement dans le nom du fichier
  Alpha1dS <- as.numeric(substring(ls_files_Distributions[1], 37, 39))  
  Alpha2dI <- as.numeric(substring(ls_files_Distributions[1], 41, 43))
  
  if (length(remov) > 0){  # Si on a des seeds a enlever
    list_df_Distributions <- list_df_Distributions[-c(remov)]  # Sites retires de list_df_Distributions
    sms <- c("Simulation", Alpha1dS,"/",Alpha2dI, ", Seeds retires", remov,"\n")  # Message d'information sur les seeds manquants
    Messages <- c(Messages,sms) # Messages stockés
    
  }
  
  if (length(remov) < length(ls_files_Distributions)/2){ # On garde que les simus où au moins la moitié des réplicats s'en sont sortis
    list_paramsS = c(list_paramsS,Alpha1dS)
    list_paramsI = c(list_paramsI,Alpha2dI)
    
    # Get only commons rows between all dataframes 
    # Get all t values
    t_values <- c(0)
    for (i in 1: length(list_df_Distributions)){
      times <- list_df_Distributions[[i]]$Time
      t_values <- c(t_values, times)
    }
    t_values <- t_values[-1] # Suppress the value given for initialization
    counts <- as.data.frame(table(t_values)) # Count occurrences into a table
    counts
    # Get values that appear nbsim times (common to every df)
    nbsim <- length(list_df_Distributions)
    tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
    
    # Loop to keep only those t in the whole set of dfs
    for (i in 1:length(list_df_Distributions)){
      indexrows <- which(list_df_Distributions[[i]]$Time %in% tcommon)
      list_df_Distributions[[i]] <- subset(list_df_Distributions[[i]][indexrows,])
    }
    
    ## II. Building the mean simulation ( returns a single DataFrame)
    namecols<- names(list_df_Distributions[[1]])
    MeanSim3 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_Distributions[[1]]$Time)))
    colnames(MeanSim3) <- namecols
    
    for (i in 1:length(namecols)){
      meanvar <- MeanSimValue(list_df_Distributions,namecols[i])
      MeanSim3[,i] <- meanvar
    }
    
    DistribFiles <- c(DistribFiles, list(MeanSim3))
  }
  
  ##### DENSITIES FILES #############
  
  cwd <- setwd(DensitiesDir)
  cwd
  
  
  ls_files_Densities <- list.files(DensitiesDir) #list files
  list_tibbles_densities <- lapply(ls_files_Densities, read_csv) 
  list_df_densities <- lapply(list_tibbles_densities, as.data.frame) 
  
  
  # I. Shaping dataframes and getting common rows
  for (i in 1:length(list_df_densities)){
    list_df_densities[[i]]<- subset(list_df_densities[[i]][c(-1)]) # suppress first column
    list_df_densities[[i]]$Time <- round(list_df_densities[[i]]$Time, 1) #round t values
    list_df_densities[[i]] <- distinct(list_df_densities[[i]], Time, .keep_all = TRUE)#remove duplicates
  }
  
  
  if (length(remov) > 0){ # Retire les seeds comme avec les Distributions files
    list_df_densities <- list_df_densities[-c(remov)]
  }
  
  if (length(remov) < length(ls_files_Distributions)/2){
    
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
  
    ## II. Building the mean simulation ( returns a single DataFrame)
    namecols<- names(list_df_densities[[1]])
     MeanSim2 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_densities[[1]]$Time)))
    colnames(MeanSim2) <- namecols
  
    for (i in 1:length(namecols)){
    
      meanvar <- MeanSimValue(list_df_densities,namecols[i])
      MeanSim2[,i] <- meanvar
    }
  
    DensitiesFiles <- c(DensitiesFiles, list(MeanSim2))
    
  }

} # LOOP 1 : THAT EXTRACT, COMPILE AND BUILD A MEAN SIM WITH DISTRIBUTION AND DENSITIES FILES
cat(Messages) # Messages avec tous les seeds "problematique"




### If use of loop 1, Launch all the part between the ------, Few details can be change in this part according to your simulation, like tmax, NBROW, ...
### ---------------------------------------------------------------------------------------------------------------------------------------------------
### Extract values of virulence at convergence for the heatmap
### Heatmap / Contour plot output

Alp = seq(0.01, 0.49, 0.01) # Valeur de virulence possible présente dans les Distributions files
tmax = 1000 # tmax des simulations (1500 pour départ pop monomorphe et 1000 pour depart avec valeurs differentes)
Possible_Timesteps = c(tmax-50,tmax-100,tmax-200) # Differents timesteps pour regarder la moyenne à la fin de la simu 


######## Distribution Files Manipulation for extraction of virulence values #######

NBROW = 35   # NB de simus etudies (Nb de couples dS/dI differents) - a modifier selon le nombre de couples
HeatFinalPrep = data.frame(matrix(NA, nrow = NBROW, ncol = 0)) # Preparation du df
for (j in 1:length(Possible_Timesteps)){ # Moyenne de la virulence pour chaque timesteps
  cat("TimeSteps", j, "\n")
  HeatValues = NULL
  Step = 1
  for (i in 1:length(DistribFiles)){  # Moyenne de la virulence pour chaque simus
    cat(Step, "/", length(DistribFiles), "\n")
    ActualDf <- DistribFiles[[i]]
    ActualDf <- ActualDf[ActualDf$Time > Possible_Timesteps[j],]
    FinalDf <- subset(ActualDf[c(-1)])
    ConvergeValue <- NULL
    
    for (k in 1:nrow(FinalDf)){  # Calcul de la valeur de viru moyenne pour chaque pas de temps
      Alfinal = NULL
      for (l in 1:length(Alp)){
        Alfi <- (FinalDf[k,l]*Alp[l]) 
        Alfinal <- c(Alfinal, Alfi)
      }
      Moyalf <- sum(Alfinal)/sum(FinalDf[k,])
      ConvergeValue <- c(ConvergeValue,Moyalf)
    }
    ConvergeValues <- sum(ConvergeValue)/length(ConvergeValue) # Moy de tous les pas
    HeatValues <- c(HeatValues,ConvergeValues)
    Step = Step + 1
  }
  HeatFinalPrep[j]=HeatValues # Stockage des valeurs de Viru dans le df
}


######## Densities Files Manipulation #######

## Taille totale de la métapop
SommeTotMeta = NULL
for (i in 1:length(DensitiesFiles)) {
  somme <- rowSums(DensitiesFiles[[i]]) # Somme de tous les sites
  Convergeso <- somme[(length(somme)-1000):length(somme)] # On garde que les derniers pas de temps 
  # (Le -1000 correspond aux 1000 dernières lignes des fichiers de données, ce qui correspond aux 100 derniers pas de temps car on garde tous les 0.1)
  Convergesomme <- sum(Convergeso)/length(Convergeso) # Moyenne
  SommeTotMeta <- c(SommeTotMeta,Convergesomme)
}
HeatFinalPrep$SommeTotMeta <- SommeTotMeta


## Taille totale par patch moyenne
SommePatchMean = NULL
Nbpatch = 80  # NB de patch dans la metapop
for (i in 1:length(SommeTotMeta)) {
  SommePatch <- SommeTotMeta[i]/Nbpatch
  SommePatchMean <- c(SommePatchMean,SommePatch)
}
HeatFinalPrep$SommePatchMean <- SommePatchMean


## Prévalence local du parasite et Taille de la population des infectés à "l'équilibre"
PrevLocalMean = NULL
PopulSiteI = NULL
for (i in 1:length(DensitiesFiles)) {
  cat("Infect N°",i,"/",length(DensitiesFiles),"\n")
  ActualDf = DensitiesFiles[[i]]
  taille_col = (length(ActualDf)-1) 
  nbsites <- taille_col/2
  generate_odd_indexes = seq(3, taille_col+1, 2)
  Ipops <- subset(ActualDf[generate_odd_indexes]) # Isolement des populations infectes dans chaque site
  taillemax <- nrow(Ipops)
  Ipops <-Ipops[(taillemax-1000):taillemax,] # Derniers pas de temps pour les moyennes
  
  ConvergeValuesI = NULL
  for (j in 1:nrow(Ipops)){
    ActualLine <- Ipops[j,]
    ActualLine <- ActualLine[cumsum(ActualLine) > 5] 
    ConvergeValue <- sum(ActualLine)/length(ActualLine)
    ConvergeValuesI <- c(ConvergeValuesI,ConvergeValue)
  }
  PopsiteI <- sum(ConvergeValuesI)/length(ConvergeValuesI)
  PopulSiteI <- c(PopulSiteI,PopsiteI)
  Prev <- PopsiteI/SommePatchMean[i]
  PrevLocalMean <- c(PrevLocalMean,Prev)
}
HeatFinalPrep$PrevLocalMean <- PrevLocalMean
HeatFinalPrep$PopulSiteI <- PopulSiteI


## Taille de la population des susceptibles à "l'équilibre"
PopulSiteS = NULL
for (i in 1:length(DensitiesFiles)) {
  cat("Susceptibles N°",i,"/",length(DensitiesFiles),"\n")
  ActualDf = DensitiesFiles[[i]]
  taille_col = (length(ActualDf)-1)
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2)
  Spops <- subset(ActualDf[generate_even_indexes])
  taillemax <- nrow(Spops)
  Spops <- Spops[(taillemax-1000):taillemax,]
  
  ConvergeValuesS = NULL
  for (j in 1:nrow(Spops)){
    ActualLine <- Spops[j,]
    ActualLine <- ActualLine[cumsum(ActualLine) > 5]
    ConvergeValue <- sum(ActualLine)/length(ActualLine)
    ConvergeValuesS <- c(ConvergeValuesS,ConvergeValue)
  }
  PopsiteS <- sum(ConvergeValuesS)/length(ConvergeValuesS)
  PopulSiteS <- c(PopulSiteS,PopsiteS)
}
HeatFinalPrep$PopulSiteS <- PopulSiteS


## Densité de sites avec des infectés et de sites avec des suseptibles à l'équilibre

InfectDens = NULL  # Save Infected Density
SuscepDens = NULL  # Save susceptible Density
for (i in 1:length(DensitiesFiles)){
  cat("Fichier N°", i, "/",length(DensitiesFiles),"\n")
  ActualDf <- DensitiesFiles[[i]]
  taille_col = (length(ActualDf)-1) 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) 
  generate_odd_indexes = seq(3, taille_col+1, 2)
  
  # Separation des pops S et I et en gardant que les derniers pas de temps
  Spops <- subset(ActualDf[generate_even_indexes])
  Spops <- Spops[(nrow(Spops)-500):nrow(Spops),]
  Ipops <- subset(ActualDf[generate_odd_indexes])
  Ipops <- Ipops[(nrow(Ipops)-500):nrow(Ipops),]
  
  ## Densite d'Infectes
  # En pourcentage, non utilise
  Ipourcents  = data.frame(matrix(NA, nrow = nrow(Ipops), ncol = length(Ipops)))
  for (i in 1:length(Ipops)){
    for (j in 1:nrow(Ipops)){
     Ipourcents[j,i] = (Ipops[j,i]*100)/(Spops[j,i]+Ipops[j,i])
    }
  }
  
  # Presence of infected in a patch (Can be modify), if presence -> 1, if not -> 0
  Ipops[Ipops < 1] <- 0
  Ipops[Ipops >= 1 ] <- 1
  Imean <- rowSums(Ipops)
  
  InfDens <- sum(Imean)/length(Imean)
  InfectDens <- c(InfectDens,InfDens)
  
  
  # Densité de Susceptibles
  Spourcents  = data.frame(matrix(NA, nrow = nrow(Spops), ncol = length(Spops)))
  for (i in 1:length(Spops)){
    for (j in 1:nrow(Spops)){
      Spourcents[j,i] = (Spops[j,i]*100)/(Spops[j,i]+Ipops[j,i])
    }
  }
  
  Spops[Spops < 1] <- 0
  Spops[Spops >= 1 ] <- 1
  Smean <- rowSums(Spops)

  SusDens <- sum(Smean)/length(Smean)
  SuscepDens <- c(SuscepDens,SusDens)
}

HeatFinalPrep$InfectDens = InfectDens
HeatFinalPrep$SuscepDens = SuscepDens


### Liste des parametres dS et dI recuperes dans Loop 1
length(list_paramsS)
length(list_paramsI)

HeatFinalPrep$list_paramsS = list_paramsS
HeatFinalPrep$list_paramsI = list_paramsI


### ----------------------------------------------------------------------------------------------------------------------------------------------


### 2 diferents ways at this point

## If You made Loop 1 and the previous part
HeatFinal = HeatFinalPrep
HeatFinalSimu2 = HeatFinalPrep
write.csv(HeatFinalPrep, "C:\\Users\\gaze\\Downloads\\HeatFinal.csv", row.names=FALSE) # Save data


## If you already have the csv files
# GET FILES ALREADY SIMULATE
dir <- choose.dir()
cwd <- setwd(dir)
cwd

HeatFinal_2 = read.csv("HeatFinal_2.csv",header = TRUE)
HeatFinal_added = read.csv("HeatFinal_added.csv",header = TRUE) 
#### Pour HeatFinal_added, changement dans la liste des paramètres dS et dI car Loop 1 déconne au niveau de l'extraction des valeurs, changement manuel
list_paramsS = c(0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15,0.15,0.15,
                 0.2,0.2,0.2,0.2,0.2,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
                 0.3,0.3,0.3,0.3,0.3,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,
                 0.4,0.4,0.4,0.4,0.4,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,
                 0.5,0.5,0.5,0.5,0.5,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,
                 0.6,0.6,0.6,0.6,0.6,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,
                 0.7,0.7,0.7,0.7,0.7,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,
                 0.8,0.8,0.8,0.8,0.8,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,
                 0.9,0.9,0.9,0.9,0.9,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,
                 1.0,1.0,1.0,1.0,1.0,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,
                 1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.15,1.15,1.15,1.15,1.15,1.15,
                 1.2,1.2,1.2,1.2,1.25,1.25)
list_paramsI = c(1.1,1.2,1.3,1.4,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,
                 1.1,1.2,1.3,1.4,1.5,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.35,
                 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,0.05,0.15,0.25,0.35,0.45,0.55,
                 0.1,0.2,0.3,0.4,0.15,0.25)
HeatFinal_added$list_paramsS = list_paramsS
HeatFinal_added$list_paramsI = list_paramsI

HeatFinalNEWHORIZONS = read.csv("HeatFinalNEWHORIZONS.csv",header = TRUE)
HeatFinalFINAL = read.csv("HeatFinalFINAL.csv",header = TRUE)
HeatFinal = read.csv("HeatFinal.csv",header = TRUE)

## Bind of all df
HeatFinal = HeatFinal_2
HeatFinal = rbind(HeatFinal,HeatFinal_added)
HeatFinal = rbind(HeatFinal,HeatFinalNEWHORIZONS)
HeatFinal = rbind(HeatFinal,HeatFinalFINAL)


## Add Column dS-dI and mean of dS and dI for the heatmap
HeatFinal$dS_dI <-(HeatFinal$list_paramsS-HeatFinal$list_paramsI) # dS-dI
HeatFinal$dSdIby2 <- (HeatFinal$list_paramsS + HeatFinal$list_paramsI)/2 # Moyenne des deux
HeatFinal$dS_dI[281] = -0.1 ### -0.1 / 1.1 (dS-dI/dSdIby2)



################################# Calcul du R0 ###########################################################################################
# Parameters value
b = 2
mu = 0.5
gamma = 2.5
Beta0 = 1
K=500

BetaValue = NULL  # Beta value for the parasite
R0Value = NULL    # R0 value with the values of virulence simulated
RatioR0Value = NULL # Ratio R0 (1-1/R0)
for (i in 1:nrow(HeatFinal)){
  Beta = Beta0*HeatFinal[i,3]/(HeatFinal[i,3] + 1)  # Calcul de Beta
  BetaValue = c(BetaValue,Beta)
  R0 = (Beta*(K-K*((mu+HeatFinal$list_paramsS[i])/b)))/(HeatFinal$list_paramsI[i] + HeatFinal[i,3] + gamma + mu) # R0 calculation
  R0Value = c(R0Value,R0)
  RatioR0 = 1-(1/R0)
  RatioR0Value = c(RatioR0Value,RatioR0) # Ratio R0 Calculation
}
HeatFinal$BetaValue = BetaValue
HeatFinal$R0Value = R0Value
HeatFinal$RatioR0Value = RatioR0Value

# Graphique du R0
valuesdS = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)
dataFin = NULL
for (i in 1:length(valuesdS)){
  HeatAnalysis = HeatFinal[HeatFinal$list_paramsS == valuesdS[i],] # Select part of df with values of dS
  data <- as_tibble(HeatAnalysis)
  dataDF <- data %>% arrange(list_paramsI)   # Manipulation to have a good df where values of virulence are in order according to dI            
  ActualDF <- dataDF[,c("list_paramsS", "list_paramsI","R0Value")] # Extraction of few columns for the graphic
  dataFin <- rbind(dataFin, ActualDF)
}
dataFin$mS = dataFin$list_paramsS
ggplot(dataFin, aes(x = list_paramsI, y = R0Value, color = mS)) + geom_line(aes(group = list_paramsS), linewidth=1, linetype = "dashed") +
  scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 0.4) + geom_point(size = 3) +
  xlab(expression("Taux de dispersion des infectés(t"^-1*")")) + ylab(expression("Valeur du R"[0])) + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))


######################################################################################################################################################



################################### Maximisation #####################################################################################################
###### Maximisation R0, Maximisation de la capacité de compétition 

R0 = function(alpha) {
  (k - k*((dS+mu)/b))*(Beta0*alpha/(alpha + 1)) / (dI+alpha+gamma+mu) # R0 formula
}

alphaopt<-c()
for (i in 1:nrow(HeatFinal)){
  dS <- HeatFinal$list_paramsS[i]
  dI <- HeatFinal$list_paramsI[i]
  res <- optimize(R0, interval = c(0,5),maximum = TRUE)
  alphaopt<-c(alphaopt, res$maximum)
}

HeatFinal$R0Max = alphaopt


### Maximisation I*(1-1/R0), Infected colonisation

r=1.5 # (b-mu)

# I*(1-1/R0) Formula
IeqR0 = function(alpha){
  (1-1/(k- dS*k/r)*alpha*Beta0 / ((1+alpha)*(dI+alpha+gamma+mu))) * (-(((1+alpha)*(dI+alpha+gamma+mu)*(dI*r*(1+alpha)+k*alpha*Beta0*(dS+mu)+r*(gamma+mu+alpha*(1+alpha-k*Beta0+gamma+mu))))/(alpha*Beta0*(dI*(r+r*alpha+k*alpha*Beta0)+k*alpha*Beta0*(alpha+mu)+r*(1+alpha)*(alpha+gamma+mu)))))
}

alphaopt<-c()
for (i in 1:nrow(HeatFinal)){
  dS <- HeatFinal$list_paramsS[i]
  dI <- HeatFinal$list_paramsI[i]
  res <- optimize(IeqR0, interval = c(0,50),maximum = TRUE)
  alphaopt<-c(alphaopt, res$maximum)

}

HeatFinal$RatioR0Max = alphaopt


### Graphiques

HeatAnalysis = HeatFinal[HeatFinal$list_paramsS == 0.5,] # Select of a part of data with a values of dS fixed
HeatAnalysis = HeatAnalysis[c(3,11,12,20)] # Sélection des colonnes V3, list_paramsI, RatioR0Max et R0Max

# To change according to the maximisation we want to see on the graphic
AlphaMaximum = HeatAnalysis$RatioR0Max 
AlphaMaximum = HeatAnalysis$R0Max 

V3Maximum = HeatAnalysis$V3
mI = HeatAnalysis$list_paramsI

ggplot(HeatAnalysis, aes(AlphaMaximum,V3Maximum, color = mI)) + geom_point() + geom_abline(intercept = 0, slope = 1) + 
  xlab(expression("Virulence du parasite optimise colonisation (t"^-1*")")) + ylab(expression("Virulence du parasite simulé (t"^-1*")")) +
  xlim(0,0.6) + ylim(0,0.6) + scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 3) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))


######################################################################################################################################################


################################### Colonisation #####################################################################################################
rho = 0.9
epsilon = 0.1
NBpatch = 80

# Colonisation par Infectés (parasite)
ColoniInf = NULL
for (i in 1:nrow(HeatFinal)){
  ColInf <- HeatFinal$list_paramsI[i]*(1-rho)*HeatFinal$PopulSiteI[i]*HeatFinal$RatioR0Value[i] 
  ColoniInf <- c(ColoniInf,ColInf)
}
HeatFinal$ColoniInf = ColoniInf

# Colonisation par Hôtes 
gamma = 2.5
mu = 0.5
ColoniHost = NULL
for (i in 1:nrow(HeatFinal)){
  ud <- (1-((mu+HeatFinal$list_paramsS[i])/b))
  Param <- gamma/(gamma + mu + HeatFinal$list_paramsI[i] + HeatFinal$V3[i])
  ColHost <- (((HeatFinal$list_paramsS[i]*(1-rho)*HeatFinal$PopulSiteS[i])+
                ((HeatFinal$list_paramsI[i]*(1-rho)*HeatFinal$PopulSiteI[i])*Param)) * ud)
  
  ColoniHost <- c(ColoniHost,ColHost)
}
HeatFinal$ColoniHost = ColoniHost


## Graphique Colonisation des hôtes
valuesdS = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)
valuesdI = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2)
dataFin = NULL
for (i in 1:length(valuesdI)){
  HeatAnalysis = HeatFinal[HeatFinal$list_paramsI == valuesdI[i],]
  data <- as_tibble(HeatAnalysis)
  dataDF <- data %>% arrange(list_paramsI)
  ActualDF <- dataDF[,c("list_paramsS", "list_paramsI","ColoniHost")]
  dataFin <- rbind(dataFin, ActualDF)
}
dataFin$mI = dataFin$list_paramsI
ggplot(dataFin, aes(x = list_paramsS, y = ColoniHost, color = mI)) + geom_line(aes(group = list_paramsI), linewidth=1, linetype = "dashed") +
  scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 1.6) + geom_point(size = 3) +
  xlab(expression("Taux de dispersion des individus susceptibles (t"^-1*")")) + ylab("Nb de colonisateurs susceptibles par site et pour sa durée de vie")+ 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))

## Graphique Colonisation des hôtes infectés
valuesdS = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)
dataFin = NULL
for (i in 1:length(valuesdS)){
  HeatAnalysis = HeatFinal[HeatFinal$list_paramsS == valuesdS[i],]
  data <- as_tibble(HeatAnalysis)
  dataDF <- data %>% arrange(list_paramsS)
  ActualDF <- dataDF[,c("list_paramsS", "list_paramsI","ColoniInf")]
  dataFin <- rbind(dataFin, ActualDF)
}
dataFin$mS = dataFin$list_paramsS
ggplot(dataFin, aes(x = list_paramsI, y = ColoniInf, color = mS)) + geom_line(aes(group = list_paramsS), linewidth=1, linetype = "dashed") +
  scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 0.6) + geom_point(size = 3) +
  xlab(expression("Taux de dispersion des individus infectés (t"^-1*")")) + ylab("Nb de colonisateurs infectés par site et pour sa durée de vie")+ 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))


######################################################################################################################################################


################################### Other graphics #####################################################################################################


### Virulence values Graphic (not the big heatmap)
valuesdS = c(0.1,0.3,0.5,0.7,0.9,1.1) # Selection of few values of dS
dataFin = NULL
for (i in 1:length(valuesdS)){
  HeatAnalysis = HeatFinal[HeatFinal$list_paramsS == valuesdS[i],]
  data <- as_tibble(HeatAnalysis)
  dataDF <- data %>% arrange(list_paramsI)
  ActualDF <- dataDF[,c("list_paramsS", "list_paramsI","V3")]
  dataFin <- rbind(dataFin, ActualDF)
}
dataFin$mS = dataFin$list_paramsS
ggplot(dataFin, aes(x = list_paramsI, y = V3, color = mS)) + geom_line(aes(group = list_paramsS), linewidth=1, linetype = "dashed") +
  scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 0.5) + geom_point(size = 3) +
  xlab(expression("Taux de dispersion des infectés(t"^-1*")")) + ylab(expression("Virulence du parasite (t"^-1*")")) + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))



## Graphique Prévalence local du parasite
valuesdS = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)
dataFin = NULL
for (i in 1:length(valuesdS)){
  HeatAnalysis = HeatFinal[HeatFinal$list_paramsS == valuesdS[i],]
  data <- as_tibble(HeatAnalysis)
  dataDF <- data %>% arrange(list_paramsI)
  ActualDF <- dataDF[,c("list_paramsS", "list_paramsI","PrevLocalMean")]
  dataFin <- rbind(dataFin, ActualDF)
}
dataFin$mS = dataFin$list_paramsS
ggplot(dataFin, aes(x = list_paramsI, y = PrevLocalMean, color = mS)) + geom_line(aes(group = list_paramsS), linewidth=1, linetype = "dashed") +
  scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 0.6) + geom_point(size = 3) +
  xlab(expression("Taux de dispersion des individus infectés (t"^-1*")")) + ylab("Prévalence locale du parasite")+ 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))


## Graphique Nb de Sites avec des infectés
valuesdS = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)
dataFin = NULL
for (i in 1:length(valuesdS)){
  HeatAnalysis = HeatFinal[HeatFinal$list_paramsS == valuesdS[i],]
  data <- as_tibble(HeatAnalysis)
  dataDF <- data %>% arrange(list_paramsI)
  ActualDF <- dataDF[,c("list_paramsS", "list_paramsI","InfectDens")]
  dataFin <- rbind(dataFin, ActualDF)
}
dataFin$mS = dataFin$list_paramsS
ggplot(dataFin, aes(x = list_paramsI, y = InfectDens, color = mS)) + geom_line(aes(group = list_paramsS), linewidth=1, linetype = "dashed") +
  scale_colour_gradient2(low = "blue",mid ="green", high = "yellow",midpoint = 0.6) + geom_point(size = 3) +
  xlab(expression("Taux de dispersion des individus infectés (t"^-1*")")) + ylab("Nombre de sites avec individus infectés")+ 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12),
        legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'))



######################################################################################################################################################


################################### Heatmap part #####################################################################################################

### Manipulation of the dataframe
Matrixdf = data.frame(matrix(nrow = length(HeatFinal$dS_dI) , ncol = 0))
Matrixdf$dS_dI = as.factor(HeatFinal$dS_dI)
Matrixdf$dSdIby2 = as.factor(HeatFinal$dSdIby2)


### Choice of the value we want to study in fonction of dS-dI and dSdI/2

## Taille totale de la métapop
Matrixdf$SommeTotMeta = as.numeric(HeatFinal$SommeTotMeta)
Matrix = xtabs(SommeTotMeta ~dS_dI + dSdIby2, Matrixdf)

## Taille totale par patch moyenne
Matrixdf$SommePatchMean = as.numeric(HeatFinal$SommePatchMean)
Matrix = xtabs(SommePatchMean ~dS_dI + dSdIby2, Matrixdf)

## Taille totale par patch moyenne
Matrixdf$PopulSiteI = as.numeric(HeatFinal$PopulSiteI)
Matrix = xtabs(PopulSiteI ~dS_dI + dSdIby2, Matrixdf)

## Prévalence local du parasite
Matrixdf$PrevLocalMean = as.numeric(HeatFinal$PrevLocalMean)
Matrix = xtabs(PrevLocalMean ~dS_dI + dSdIby2, Matrixdf)

## Convergence Values
Matrixdf$HeatValues = as.numeric(HeatFinal[,3])  # Valeur dépendant du Timesteps désirés
Matrix = xtabs(HeatValues ~dS_dI + dSdIby2, Matrixdf)

## RO
Matrixdf$R0Value = as.numeric(HeatFinal$R0Value)
Matrix = xtabs(R0Value ~dS_dI + dSdIby2, Matrixdf)

## Ratio R0
Matrixdf$RatioR0Value = as.numeric(HeatFinal$RatioR0Value)
Matrix = xtabs(RatioR0Value ~dS_dI + dSdIby2, Matrixdf)

## Colonisation par infectés (Parasites)
Matrixdf$ColoniInf = as.numeric(HeatFinal$ColoniInf)
Matrix = xtabs(ColoniInf ~dS_dI + dSdIby2, Matrixdf)

## Colonisation par les hôtes
Matrixdf$ColoniHost = as.numeric(HeatFinal$ColoniHost)
Matrix = xtabs(ColoniHost ~dS_dI + dSdIby2, Matrixdf)


### And a final manipulation to replace every 0 by NA
Matrix[Matrix == 0] <- NA
Matridf = as.data.frame(Matrix)


### ggplot with lines showing viability limits of the model
p <- ggplot(Matridf, aes(dSdIby2, dS_dI, fill= Freq)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  theme_ipsum() +  
  ggtitle("") + 
  geom_abline(intercept = 54, slope = -1) +
  geom_abline(intercept = 64, slope = -1, color = "red") +
  geom_abline(intercept = 56, slope = 1) +
  geom_abline(intercept = 84, slope = -1) +
  geom_text(x=10, y=20, label="a", size = 15) +
  geom_text(x=4, y=63, label="a", size = 15) +
  geom_text(x=45, y=30, label="c", size = 15) +
  geom_text(x=45, y=55, label="b", size = 15) +
  xlab("(mS + mI)/2") + ylab("mS - mI") + labs(fill = "Vir")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 15))

### Version without lines
p <- ggplot(Matridf, aes(dSdIby2, dS_dI, fill= Freq)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  theme_ipsum() +  
  ggtitle("") + 
  xlab("(mS + mI)/2") + ylab("mS - mI") + labs(fill = "Vir")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 15))

ggplotly(p, tooltip="text")


######################################################################################################################################################






######--------------------------------------------------------------------------------------------------------------------

## Traits dynamics of simulation, evolution of the virulence values through time
i = 1
ActualDf = DistribFiles_4[[i]]
# ActualDf[ActualDf$Time > 80 & ActualDf$Time < 200,]
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
plot(ActualDf$Time, Moyalpha, type = 'l',xlab = 'time', ylab = 'Mean Trait Value', xlim = c(0,1500),ylim = c(0,0.5))

#####----------------------------------------------------------------------------------------------------------------------



##### RMSE Part ###############################################################################################

DataBased = HeatFinal[HeatFinal$list_paramsS == 0.1,]
DataBased = DataBased[DataBased$list_paramsI < 1.9,]
DataBased = DataBased[DataBased$list_paramsI > 1.0,]
DataMulti = HeatFinalSimu1[HeatFinalSimu1$list_paramsS == 0.1,]


RMSEValue <- sqrt(mean((DataBased$V3 - DataMulti2$V3)^2))
MSEValue <- mean((DataBased$V3 - DataMulti$V3)^2)

DataMulti$Percentage = (RMSEValue*100)/DataMulti$V3
mean(DataMulti$Percentage)

plot(DataBased$SommeTotMeta, DataMulti$SommeTotMeta)
plot(DataBased$V3, DataMulti2$V3, xlim = c(0.1,0.3), ylim = c(0.1,0.3))
plot(DataBased$PrevLocalMean, DataMulti$PrevLocalMean)
abline(0,1)

