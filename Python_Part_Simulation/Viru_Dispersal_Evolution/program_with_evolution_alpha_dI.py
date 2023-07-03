# STOCHASTIC SIMULATION ALGORITHM FOR EVOLUTION VIRULENCE AND INFECTED DISPERSAL IN SIR METAPOPULATION
import os
from copy import deepcopy
from datetime import datetime
from collections import defaultdict
import NetworkFunctions
import classes_double_evolution
import fonctions_double_evolution
import numpy as np
import pandas as pd
import multiprocessing
import Params_double_evolution

# A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
# Adapted for metapopulations modelling and networks

#POPULATION DYNAMICS PARAMETERS
beta0 = Params_double_evolution.beta0 # Infectious contact rate
bi = Params_double_evolution.bi  # Per capita net growth rate
mu = Params_double_evolution.mu # Per capita natural death rate
dI = Params_double_evolution.dI  # Dispersal infected
omega = Params_double_evolution.omega # Strength of density dependance on births
k = Params_double_evolution.k  # Carrying capacity
#d = 0.5# Dispersal propensity
gamma = Params_double_evolution.gamma  # Parasite Clearance
alpha = Params_double_evolution.alpha  # Parasite Virulence
rho = Params_double_evolution.rho  # Dispersal Cost
epsilon = Params_double_evolution.epsilon  # Extinction rate

def RunModel(seed, paramSU, epsilon, paramINF) :

    #SIMULATION PARAMETERS
    Sp_config = "Island"
    # This part changes parameters (here d) values for multisim runs where numerous are tested
    dS = paramSU
    eps = epsilon
    classes_double_evolution.epsilon = epsilon
    classes_double_evolution.dS = paramSU

    fonctions_double_evolution.epsilon = epsilon
    fonctions_double_evolution.dS = paramSU


    print('Seed', seed)
    np.random.seed(seed) #Set seed for reproducibility
    nb_iterations = 0 # Store the number of interations, used to define times that are saved (later)
    sim_time = 0 # Simulation time (model time, not an iteration number)
    vectime = [0] # to keep track of t variable
    vectime_double = [0] # To keep track of t variable in case of the double evolution files
    tmax = 1500 # Ending time
    Nexactsteps = 20  # Number of steps to do if/when performing direct method (USELESS IF nbsite > 20~30)
    nbsite = 80 # Number of sites
    n = 7 #Number of rows  for lattices configurations
    p= 7 # Number of columns for lattices configurations
    Taillepop = Params_double_evolution.k # Initial local population sizes
    ## Define the trait that will evolve in the simulation
    EvoltraitAlp = classes_double_evolution.EvolvingTraitAlp('alpha', True) # Values of Alpha virulence
    EvoltraitdINF = classes_double_evolution.EvolvingTraitdINF('dI', True) # Values of infected dispersal


    # FOLDER STORAGE
    # Create new folder for storage simulations with subfolder Densities, Distributions and Traits
    dirfile = "dS=" + str(paramSU) +'_' + "e=" + str(eps)
    rep_actuel = "D:/gaze/donnees/bureau/Stage M2/Dispersal-2-modelisation_dev/" + dirfile
    dirtrue = os.path.exists(rep_actuel)

    if dirtrue == False:
        os.mkdir("D:/gaze/donnees/bureau/Stage M2/Dispersal-2-modelisation_dev/" + dirfile)
        os.mkdir(rep_actuel + "/Densities")

        os.mkdir(rep_actuel + "/Distributions_alpha")
        os.mkdir(rep_actuel + "/Distributions_dI")
        os.mkdir(rep_actuel + "/Distributions_AlpdI")

        os.mkdir(rep_actuel + "/Traits_alpha")
        os.mkdir(rep_actuel + "/Traits_dI")
    # Creation of several folder for each output
    rep_densities = rep_actuel + "/Densities"
    rep_distrib_alp = rep_actuel + "/Distributions_alpha"
    rep_distrib_dI = rep_actuel + "/Distributions_dI"
    rep_distrib_AlpdI = rep_actuel + "/Distributions_AlpdI"
    rep_traits_alp = rep_actuel + "/Traits_alpha"
    rep_traits_dI = rep_actuel + "/Traits_dI"

    #ENABLE OUTPUT STORAGE
    dico_densities_df = {} # Track population densities
    dico_distrib_alp_df = {} # Distribution of traits values of alpha virulence
    dico_distrib_dI_df = {}  # Distribution of traits values of infected dispersal
    dico_distrib_Alp_dI_df = {} # Distribution of couple of values for Alpha and virulence
    dico_traits_alp_df = {} # Track mean trait value per site, For alpha virulence
    dico_traits_dI_df = {}  # Track mean trait value per site, For infected dispersal

    #CREATE THE METAPOPULATION
    ListSites = fonctions_double_evolution.SetMetapop(nbsite, Taillepop)

    #SPATIAL CONFIGURATIONS OPTIONS AND ENABLING
    if Sp_config == "Square lattice" :
        Hastable_adjacency = NetworkFunctions.Build_Square_Lattice(n,p,ListSites)
        nb_neighbors = 4
    if Sp_config == "Accordion":
        nb_neighbors = 4
        Hastable_adjacency = NetworkFunctions.Build_Accordion_Neighboring(ListSites, nb_neighbors)
    if Sp_config == "Island" : pass # The programm was initially designed for island (complete graph) pop so we just go with the flow
    if Sp_config == "Hexagonal lattice" :
        Hastable_adjacency = NetworkFunctions.Build_Hexagonal_Lattice(n,p,ListSites) # n the number or row desired (int), p the number of columns (int) , ListSite a list of site objects

     #Gives adjacency lists
    #print(Hastable_adjacency)

    #EVENT DEFINITIONS
    ReproductionS = classes_double_evolution.Event(name='Reproduction S',propensity='(bi - omega * (self.S+self.I) ) * self.S', Schange='1', Ichange='0', order=1,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    DeathS = classes_double_evolution.Event(name='Death S',propensity='mu*self.S', Schange='-1', Ichange='0', order=3,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    DispersalS = classes_double_evolution.Event(name='Dispersal S',propensity='dS*self.S', Schange='-1', Ichange='0', order=1,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    DispersalI = classes_double_evolution.Event(name='Dispersal I',propensity='(dS*dI)*self.I', Schange='0', Ichange='-1', order=1,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    Extinction = classes_double_evolution.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    Infection = classes_double_evolution.Event(name='Infection',propensity='beta0 *(alpha / (1+alpha) )*self.S*self.I', Schange='-1', Ichange='1', order=2,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    Recovery = classes_double_evolution.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    DeathI = classes_double_evolution.Event(name='Death I',propensity='(mu + alpha) *self.I', Schange='0', Ichange='-1', order=1,EvolvingTraitAlp=EvoltraitAlp, EvolvingTraitdINF=EvoltraitdINF)
    Events = [ReproductionS,Infection, DispersalI, DispersalS, DeathS, DeathI,Recovery, Extinction]

    #Initializing outputs storage
    Propensities_out =[] # Collect propensities outputs
    IsTrackPropensites = False #Set false/ true to not track/track propensities
    for i in range(len(Events)):
        Propensities_out.append([0]) #Store times series of propensities

    # INITIALIZE OUTPUTS FOT T=0
    # TRACK OF DENSITIES (one list of S and I per site) and MEAN TRAIT VALUES PER SITE
    for index, i in enumerate(ListSites):
        dico_densities_df[f"S{index}"]= [i.effectifS]
        dico_densities_df[f"I{index}"]= [i.effectifI]
        if i.effectifI > 0:
            dico_traits_alp_df[f"site{index}"] = [sum(i.traitvaluesAlp) / i.effectifI]
            dico_traits_dI_df[f"site{index}"] = [sum(i.traitvaluesdINF) / i.effectifI]
        else:
            dico_traits_alp_df[f"site{index}"] = ["NA"]
            dico_traits_dI_df[f"site{index}"] = ["NA"]
    #TRACK TRAIT DISTRIBUTION (more exaclty the densitiy of each trait value trough time)
    for i in fonctions_double_evolution.Rounded_alpha_values:  # For each value of alpha virulence defined
        count = 0
        for j in ListSites:  # We browse the different sites
            # We count the given value and sum it
            count += j.traitvaluesAlp.count(i)
        dico_distrib_alp_df[f"Alpha{i}"] = [count]

    for i in fonctions_double_evolution.Rounded_dI_values:  # For each value of infected dispersal defined
        count = 0
        for j in ListSites:
            count += j.traitvaluesdINF.count(i)
        dico_distrib_dI_df[f"dI={i}"] = [count]

    # TRACK TRAIT DISTRIBUTION 2 (Couple of values)
    from operator import itemgetter
    import array as arr
    for i in fonctions_double_evolution.Rounded_alpha_values:
        for j in fonctions_double_evolution.Rounded_dI_values:
            count = 0
            for l in ListSites:
                Pos = []
                length = len(l.traitvaluesAlp)
                for v in range(length):
                    if l.traitvaluesAlp[v] == i:
                        Pos.append(v)

                if Pos:
                    res_list = np.take(l.traitvaluesdINF, Pos)
                    res_list = res_list.tolist()
                    count += res_list.count(j)
                else:
                    pass
            dico_distrib_Alp_dI_df[f"Alp{i},dI{j}"] = [count]

    ############################# MODEL MAIN LOOP ########################################
    while sim_time < tmax :
        #DEFINE HOW MUCH YOU LIKE TO SAVE DATAS
        if nb_iterations % 15 == 0: # % x => we save each x times
            sim_time_double = sim_time
            vectime.append(sim_time) # Update time vector

        #COMPUTE PROPENSITIES
        Propensities, Sum_propensities = fonctions_double_evolution.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites
        #print('Propensities', Propensities)
        #print('SumPropensities', Sum_propensities)
        #print('Propensions', Propensities)
        SumS, SumI = fonctions_double_evolution.SumDensities(ListSites) # Get total densities of species
        #print('les effectifs de pop ( I puis S)', SumI, SumS)

        # Break the main loop if there are no infected remaining ( This happens essentially at start if the 1st infected dies)
        if SumI == 0:
            print('WARNING : ABORTED SIMULATION, No infected remaining')
            break
        print(sim_time, 'time')  # Kind of a loading bar but much uglier

        ################################# TAU-LEAP PART #################################

        # Most of this section is about getting Tau or die trying
        # All explicit computations are available in cao et. al. (2006)

        #Get Critical Reactions (NOT USE : NEGATIVE POPULATION ARE AVOIDED IN THE DIRTY WAY)
        #Criticals = fonctions.GetCriticals(Propensities, ListSites, Events)

        # Preliminaries
        Mus = fonctions_double_evolution.ComputeMuNSigma(Sum_propensities, Events) # As each statechange is 1 , 0, or -1 we have sigma = mu
        Epsis = fonctions_double_evolution.GetEpsilonI(SumS, SumI)
        TauPrime = fonctions_double_evolution.GetTauPrime(Epsis, Mus)

        #Main algorithm Decision tree
        aox = sum(Sum_propensities)

        if TauPrime < 10/aox : # Take 10/aox 1 is left for ignoring this part
            print('Direct Method performed')
            Tau = fonctions_double_evolution.DoDirectMethodV2(Propensities,Sum_propensities,Nexactsteps, Events,ListSites)
        else:
            #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
            #We expect that random sample in the place of reactions will be kind equivalent
            #As critical reactions in critical population will have low ocurrences
            #And if not : we let populations having -1 individuals after an event, and set it back to zero after
            Tau=TauPrime


            ################## DETERMINE HOW MANY EVENTS WILL TRIGGER, AND HOW MANY TIMES PER SITE #####################

            #Sample the mean number of realisation of a given event during tau from a poisson distribution
            Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
            #Sample the kjs in poisson law, aka the number of trigger of each event
            triggers = []
            for i in Poisson_means :
                kj = np.random.poisson(i,1)
                triggers.append(kj[0]) # The [0] is due to array structure of kj
            #Creating a vector with sites indexes, used later and put here to not be computed at each subloop iteration
            sites_indexes = []
            for i in range(len(ListSites)):
                sites_indexes.append(i)
            #Now we sample the sites where events will occur from multinomial law
            for index,event in enumerate(Events) : # For each event
                Noccur = triggers[index] #We get the number of times it should trigger during tau
                props = Propensities[index] # We get the propensity per sites
                SumProp = sum(props) # We get the total propensity
                Probas = [float(i /SumProp) for i in props] # We get probability of occurence in each site
                if Noccur == 0 : #Case where the event can't happen
                    trigger_persite=[0 for i in range(nbsite)]
                else : # Normal cases
                    trigger_persite = np.random.multinomial(Noccur, Probas)


                ########## APPLY CHANGES DUE TO EVENTS IN POPULATIONS ##################

                for index, Site in enumerate(ListSites) :

                    if 'Dispersal' in event.name :
                        # Multiply the state change in population by the number of triggers
                        Site.effectifS += trigger_persite[index] * event.Schange
                        Site.effectifI += trigger_persite[index] * event.Ichange
                        nbmigrants = max(abs(trigger_persite[index] * event.Schange), abs(trigger_persite[index] * event.Ichange))

                        # Here we delete the trait value corresponding to dispersing individual
                        # For virulence evolution, hold only for I indviduals
                        dispersers_traitvaluesAlp = []
                        dispersers_traitvaluesdINF = []
                        dispersers_beta = []
                        if Site.traitvaluesAlp and abs(
                                event.Ichange) > 0:  # If there are dispersers AND that those dispersers are infected
                            for i in range(trigger_persite[index]):  # for each disperser
                                if Site.traitvaluesAlp:  # In case we remove all trait values during the loop (induces error message
                                    #disperser = float(np.random.choice(
                                    #    Site.traitvaluesAlp))  # Get the value that is to be depleted and added to the receiving site
                                    #Index_Disperser = Site.traitvaluesAlp.index(disperser)  # Get his corresponding beta
                                    #dispersers_traitvaluesAlp.append(disperser)
                                    #dispersers_traitvaluesdINF.append(Site.traitvaluesdINF[Index_Disperser])
                                    #dispersers_beta.append(Site.betaI[Index_Disperser])
                                    #Site.traitvaluesAlp.remove(disperser)  # Remove disperser from actual site
                                    #Site.traitvaluesdINF.pop(Index_Disperser)
                                    #Site.betaI.pop(Index_Disperser)  # Remove corresponding beta

                                    SumdI = sum(Site.traitvaluesdINF)
                                    IndividualsProbasDisp = [float(i / SumdI) for i in
                                                             Site.traitvaluesdINF]  # Get Individual probability to reproduce
                                    disperser = list(
                                        np.random.multinomial(1, IndividualsProbasDisp))
                                    Index_Disperser = disperser.index(1)
                                    dispersers_traitvaluesAlp.append(Site.traitvaluesAlp[Index_Disperser])
                                    dispersers_traitvaluesdINF.append(Site.traitvaluesdINF[Index_Disperser])
                                    dispersers_beta.append(Site.betaI[Index_Disperser])
                                    Site.traitvaluesAlp.pop(Index_Disperser)  # Remove disperser from actual site
                                    Site.traitvaluesdINF.pop(Index_Disperser)
                                    Site.betaI.pop(Index_Disperser)  # Remove corresponding beta
                                else:  break

                        #Aplpy dispersal Cost
                        SuccessfulMigrants = 0
                        for i in range(nbmigrants):
                            roll4urlife = np.random.uniform(0,1,1)
                            if roll4urlife > rho : SuccessfulMigrants += 1

                        #DISTRIBUTE SUCCESSFUL MIGRANTS TO NEIGHBORS
                        for i in range(SuccessfulMigrants):
                            if Site.traitvaluesAlp and abs(event.Ichange) > 0:
                                # Get trait values of surviving individuals
                                SurvivorTraitAlp = np.random.choice(
                                    dispersers_traitvaluesAlp)  # All values are chosen with same probability
                                Index_survivor = dispersers_traitvaluesAlp.index(SurvivorTraitAlp) # Get index
                                SurvivorTraitdINF = dispersers_traitvaluesdINF[Index_survivor]
                                SurvivorBeta = dispersers_beta[Index_survivor]  # Get corresponding beta

                                dispersers_traitvaluesAlp.remove(SurvivorTraitAlp)  # Chosen value is removed
                                dispersers_traitvaluesdINF.remove(SurvivorTraitdINF)
                                dispersers_beta.remove(SurvivorBeta)  # Remove it from the list

                            #MIGRANT DISTRIBUTION FOR COMPLETE GRAPH "ISLAND"
                            if Sp_config == "Island":
                                index_sites = deepcopy(sites_indexes)  # working copy of site indexes vector
                                del index_sites[
                                    index]  # Drop the current site from the list cause you can't emigrate to the place from which you departed
                                Index_destination = np.random.choice(index_sites)  # Get index of destination site

                            #MIGRANT DISTRIBUTION FOR ANY KIND OF THING THAT RETURNED AND ADJACENCY LIST IN THE BEGINNING
                            elif Sp_config != "Island" :
                                Current_site = index # Get the current site index
                                Nearest_neighbors = Hastable_adjacency[f"{Current_site}"] # We get in the topology dictionnary neighbors corresponding to actual site
                                Index_destination = np.random.choice(Nearest_neighbors) # And choose one of them at random
                            #Add individual to destination
                            if abs(event.Schange) > 0 : #if S are dispersers
                                ListSites[Index_destination].effectifS += 1
                            elif abs(event.Ichange) > 0 : # if I are dispersers
                                ListSites[Index_destination].effectifI += 1
                                # Add trait value of individual
                                ListSites[Index_destination].traitvaluesAlp.append(SurvivorTraitAlp)
                                ListSites[Index_destination].traitvaluesdINF.append(SurvivorTraitdINF)
                                ListSites[Index_destination].betaI.append(SurvivorBeta)
                            else : print('ERROR : disperser is neither S nor I and that is very curious !') #This is useless, the error never raises
                    else:
                        if event.name == 'Extinction' :
                            #When extinction occur we need to retrieve densities values cause they're initialized to zero otherwise in class definition
                            event.Schange=-Site.effectifS
                            event.Ichange=-Site.effectifI
                            Site.effectifS += trigger_persite[index]*event.Schange
                            Site.effectifI += trigger_persite[index]*event.Ichange

                            if trigger_persite[index] > 0:  # If the extinction has really occured
                                Site.traitvaluesAlp = []
                                Site.traitvaluesdINF = []
                                Site.betaI = []
                        else :
                            if abs(event.Ichange) > 0:  # If the event has an effect on 'evolving population'
                                #print('trigger_persite[index]',trigger_persite[index])
                                #print('event.SandIchange',event.Schange,event.Ichange)
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange
                                #print('length traitvalueAlp Before', len(Site.traitvaluesAlp))
                                #print('length traitvaluedINF Before', len(Site.traitvaluesdINF))
                                #print('length BetaI Before', len(Site.betaI))
                                Site.traitvaluesAlp, Site.traitvaluesdINF, Site.betaI = fonctions_double_evolution.ChooseMuteTraitValue(EvoltraitAlp, EvoltraitdINF, trigger_persite[index],
                                                                                          event.Ichange,Site.traitvaluesAlp, Site.traitvaluesdINF, Site.betaI)
                                #print('traitvalueAlp', Site.traitvaluesAlp)
                                #print('traitvaluedINF', Site.traitvaluesdINF)
                                #print('BetaI', Site.betaI)
                            else :
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange

        ############# UPDATE THINGS AT THE END OF ITERATION #############
        #Update time
        sim_time += Tau
        #Update the output tracking
        # 1. Densities
        indexlist = 0
        for index, i in enumerate(ListSites):
            if i.effectifS < 0:  # Avoid negative population in the "big fat brute" way
                i.effectifS = 0
            if nb_iterations % 15 == 0 :
                dico_densities_df[f"S{index}"].append(i.effectifS)
            indexlist += 1
            if i.effectifI < 0:
                i.effectifI = 0
            if nb_iterations % 15 == 0:
                dico_densities_df[f"I{index}"].append(i.effectifI)
            indexlist += 1
        #2. Propensities
        if IsTrackPropensites == True :
             # Propensities of each event in a list sorted by event
            for index,propensitiy in enumerate(Sum_propensities) :
                Propensities_out[index].append(propensitiy)
        # 3. Trait Values
        indexlist2 = 0
        list_value = []
        # For alpha virulence
        for index, i in enumerate(ListSites):
            if i.effectifI > 0:
                SumTraitValuesAlp = sum(i.traitvaluesAlp)
                #print('length traitvaluesAlp AND Effectif I', len(i.traitvaluesAlp),i.effectifI)
                MeanTraitValueAlp = SumTraitValuesAlp / i.effectifI
                if nb_iterations % 15 == 0:
                    dico_traits_alp_df[f"site{index}"].append(MeanTraitValueAlp)
                indexlist2 += 1
            else:
                MeanTraitValueAlp = 'NA'
                list_value.append(MeanTraitValueAlp)
                if nb_iterations % 15 == 0:
                    dico_traits_alp_df[f"site{index}"].append(MeanTraitValueAlp)

        # For infected dispersal
        for index, i in enumerate(ListSites):
            if i.effectifI > 0:
                SumTraitValuesdINF = sum(i.traitvaluesdINF)
                #print('length traitvaluesdINF AND Effectif I',len(i.traitvaluesdINF), i.effectifI)
                MeanTraitValuedINF = SumTraitValuesdINF / i.effectifI
                if nb_iterations % 15 == 0:
                    dico_traits_dI_df[f"site{index}"].append(MeanTraitValuedINF)
                indexlist2 += 1
            else:
                MeanTraitValuedINF = 'NA'
                list_value.append(MeanTraitValuedINF)
                if nb_iterations % 15 == 0:
                    dico_traits_dI_df[f"site{index}"].append(MeanTraitValuedINF)

        # 4 Count the different phenotypes in the metapopulation, in order to follow their distribution over time
        Possible_valuesAlp = fonctions_double_evolution.Rounded_alpha_values  # For alpha virulence
        for i in Possible_valuesAlp:  # For each value defined
            count = 0
            if nb_iterations % 15 == 0:
                for j in ListSites:  # We browse the different sites
                    # We count the given value and sum it
                    count += j.traitvaluesAlp.count(i)
                dico_distrib_alp_df[f"Alpha{i}"].append(count)

        Possible_valuesdINF = fonctions_double_evolution.Rounded_dI_values  # For infected dispersal
        for i in Possible_valuesdINF:  # For each value defined
            count = 0
            if nb_iterations % 15 == 0:
                for j in ListSites:  # We browse the different sites
                    # We count the given value and sum it
                    count += j.traitvaluesdINF.count(i)
                dico_distrib_dI_df[f"dI={i}"].append(count)


        # TRACK TRAIT DISTRIBUTION 2 (Couple of values)
        if nb_iterations % 15 == 0:
            simforresults = round(sim_time, 1)
            # timetosave = [1,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400]
            timetosave = [1,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400] # Save of the situation in the sim at several timesteps
            if simforresults in timetosave:
                vectime_double.append(sim_time)
                for i in Possible_valuesAlp:       # Search for each couple of virulence/dispersal value the number of infected individuals with those values
                    for j in Possible_valuesdINF:
                        count = 0
                        # if nb_iterations % 15 == 0:
                        for l in ListSites:
                            Pos = []
                            length = len(l.traitvaluesAlp)
                            for v in range(length):
                                if l.traitvaluesAlp[v] == i:  # Extract each individuals who have a certain value of virulence
                                    Pos.append(v)
                            if Pos:
                                res_list = np.take(l.traitvaluesdINF, Pos)  # For those indivduals, extract those who a have a certain value of dispersal
                                res_list = res_list.tolist()
                                count += res_list.count(j) # COunt of individuals
                            else:
                                pass
                        dico_distrib_Alp_dI_df[f"Alp{i},dI{j}"].append(count)


        #Update the number of iterations (for times where datas are saved)
        nb_iterations += 1


    ################ AFTER THE MAIN LOOP IS ACHIEVED ################
    if IsTrackPropensites == True : #if we track propensities
        #Creating propensities dataframe
        dataprop = pd.DataFrame(columns=['t'])
        for event in Events :
            EventName = event.name
            dataprop[EventName]= []
        #Filling the dataframe
        for index, colname in enumerate(dataprop):
            if index == 0 : dataprop[colname] = vectime
            else : dataprop[colname] = Propensities_out[index-1]
        #Saving into .csv file
        dataprop.to_csv('Propensities_outputs_LongRun_Step005'+str(seed)+'.csv')


    #DENSITIES TIME SERIES
    datadensity = pd.DataFrame.from_dict(data=dico_densities_df)
    VectimeDf = pd.DataFrame(data=vectime)
    datadensity.insert(0, "Time", VectimeDf, allow_duplicates=False)
    datadensity.to_csv(os.path.join(rep_densities, 'Metapop_outputs_ESSfigure_0811' + '_' + str(dS) + "e=" + str(eps) + '_' + str(seed) + '.csv'))

    #MEAN ALPHA TRAITS TIME SERIES
    datatraitAlp = pd.DataFrame.from_dict(data=dico_traits_alp_df)
    datatraitAlp.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datatraitAlp.to_csv(os.path.join(rep_traits_alp,'Traits_Alpha_outputs_ESSfigure_0811' + '_' + str(dS) + "e=" + str(eps) + '_' + str(seed) + '.csv'))

    # MEAN dI TRAITS TIME SERIES
    datatraitdINF = pd.DataFrame.from_dict(data=dico_traits_dI_df)
    datatraitdINF.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datatraitdINF.to_csv(os.path.join(rep_traits_dI,'Traits_dI_outputs_ESSfigure_0811' + '_' + str(dS) + "e=" + str(eps) + '_' + str(seed) + '.csv'))

    #JUST THE TRAITS ALPHA TIME SERIES
    datadistribAlp = pd.DataFrame.from_dict(data=dico_distrib_alp_df)
    datadistribAlp.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datadistribAlp.to_csv(os.path.join(rep_distrib_alp,'Distribution_Alpha_outputs_ESSfigure_0811' + '_' + str(dS) + "e=" + str(eps) + '_' + str(seed) + '.csv'))

    #JUST THE TRAITS dI TIME SERIES
    datadistribdINF = pd.DataFrame.from_dict(data=dico_distrib_dI_df)
    datadistribdINF.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datadistribdINF.to_csv(os.path.join(rep_distrib_dI,'Distribution_dI_outputs_ESSfigure_0811' + '_' + str(dS) + "e=" + str(eps) + '_' + str(seed) + '.csv'))

    # COMBINATION OF TWO TRAITS
    datadistribdIAlp = pd.DataFrame.from_dict(data=dico_distrib_Alp_dI_df)
    VectimeDfDouble = pd.DataFrame(data=vectime_double)
    datadistribdIAlp.insert(0, 'Time', VectimeDfDouble, allow_duplicates=False)
    datadistribdIAlp.to_csv(os.path.join(rep_distrib_AlpdI,'Distribution_AlpdI_outputs_ESSfigure_0811' + '_' + str(dS) + "e=" + str(eps) + '_' + str(seed) + '.csv'))

################## MULTIPROCESSING PART ################

# Multiprocessing parameters
list_seeds = [7,8,9,10,11,12] # The list of seed you want to test
list_paramsS =[0.1] # The list of params values you want to test (has to be changed also at the beginning)
list_epsilon =[0.1] # Values of epsilon who can be fixed
list_paramsI = 0
#list_config =["Island"]
nbsims = len(list_seeds)

#Launch a batch of nbsims simulations
# WARNING : MULTISIM MAKES ERROR MESSAGES VANISH
if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count()-2 #Number of CPU, minus 2 by precaution. And to be able to do things meanwhile
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb) #I don't know what that is doing exactly, but it is necessary.
    BeginHour = datetime.now()
    for j in range(len(list_paramsS)):
        for i in range(nbsims):
            # Switch the two following lines to enable/disable multisim
            pool.apply_async(RunModel, args=(list_seeds[i], list_paramsS[j],list_epsilon[j],list_paramsI))  # Launch multisim
            #RunModel(list_seeds[i], list_paramsS[j],list_epsilon[j],list_paramsI)  # Launch without multisim (restores errors messages)
    pool.close() # Mandatory
    pool.join() # Idem
    LastHour = datetime.now()
    print("Départ =", BeginHour)
    print("Fin =", LastHour)