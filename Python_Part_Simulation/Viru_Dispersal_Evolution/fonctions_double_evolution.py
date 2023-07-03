# FUNCTION FOR THE BASIC METAPOPULATION MODEL with the coevolution of virulence and infected dispersal
import classes_double_evolution
import numpy as np
from copy import deepcopy
import Params_double_evolution

#Parameters extracted from param files
beta0 = Params_double_evolution.beta0 # Infectious contact rate
bi = Params_double_evolution.bi  # Per capita net growth rate
mu = Params_double_evolution.mu # Per capita natural death rate
dI = Params_double_evolution.dI  # Dispersal infected
omega = Params_double_evolution.omega # Strength of density dependance on births
k = Params_double_evolution.k  # Carrying capacity
# d = 0.5 # Dispersal propensity
gamma = Params_double_evolution.gamma  # Parasite Clearance
alpha = Params_double_evolution.alpha  # Parasite Virulence
rho = Params_double_evolution.rho  # Dispersal Cost
epsilon = Params_double_evolution.epsilon  # Extinction rate

### FOR VIRULENCE EVOLUTION, WE STATE HERE THE RANGE OF VIRULENCE VALUES THAT CAN BE TAKEN
Possible_alpha_values = list(np.arange(0.01, 0.5, 0.01)) # Here we get values from 1st to 2nd by 3rd
Rounded_alpha_values = []
for i in Possible_alpha_values :
    Rounded_alpha_values.append(round(i,2))
print('Possible alphas',Rounded_alpha_values)

### FOR INFECTED DISPERSAL EVOLUTION, WE STATE HERE THE RANGE OF DISPERSAL VALUES THAT CAN BE TAKEN
Possible_dI_values = list(np.arange(0.1,5.0,0.1))
Rounded_dI_values = []
for i in Possible_dI_values :
    Rounded_dI_values.append(round(i,2))
print('Possible dI', Rounded_dI_values)


##### THIS FUNCTION ENABLES A METAPOP : ENABLE INITIAL CONDITION FOR BOTH DEMOGRAPHY AND EVOLUTION

def SetMetapop(nbsite, taillepop): #Creates sites objects containing populations
    ListSites=[] # List that will contain all sites
    for i in range(nbsite): # Creates sites, the 1st will always contain one infected and the other 0
        if i == 0:
            newsite = classes_double_evolution.Site(effectifS=400, effectifI=100) # Anciennement S = 30 et I = 55
            newsite.Index = i
            for j in range(newsite.effectifI):
                # Start of simulation with a monomorph population (can be a change in a "multi start" of values)
                newsite.traitvaluesAlp.append(0.2)  # Alpha virulence
                newsite.traitvaluesdINF.append(1.0) # Infected dispersal
            newsite.betaI = GetBetaI(newsite.traitvaluesAlp)

        else:
            newsite = classes_double_evolution.Site(effectifS=400, effectifI=100) # Anciennement S = 30 et I = 55
            newsite.Index = i
            for j in range(newsite.effectifI):
                newsite.traitvaluesAlp.append(0.2)  # Alpha virulence
                newsite.traitvaluesdINF.append(1.0)
            newsite.betaI = GetBetaI(newsite.traitvaluesAlp)
        ListSites.append(newsite)
    return ListSites



#### FUNCTIONS FOR TAU-LEAP ALGORITHM :

def GetPropensites(Sites, Events): # Compute the propensities
    Propensities = []
    for i in Events: # For each event
        PropEvent =[]
        for j in Sites : # And each site
            S, I = j.effectifS, j.effectifI # Get the xi
            traitvalsAlp = j.traitvaluesAlp
            traitvalsdINF = j.traitvaluesdINF # Second trait add for the trait of dispersal
            betaI = j.betaI
            Prop = i.UpdatePropensity(S,I, traitvalsAlp,traitvalsdINF, betaI) #Compute propensity
            PropEvent.append(Prop)
        Propensities.append(PropEvent) # List by event AND by site
    sumpropensities = []
    for i in Propensities :
        sumpropensities.append(sum(i))
    return Propensities, sumpropensities


def GetPropensites_Per_sites(Sites, Events): # Compute the propensities
    Propensities_per_sites = []
    for i in Sites: # For each event
        PropSite =[]
        S, I = i.effectifS, i.effectifI # get the xi
        traitvalsAlp = i.traitvaluesAlp
        traitvalsdINF = i.traitvaluesdINF
        betaI = i.betaI
        for j in Events : # And each site
            Prop = j.UpdatePropensity(S,I,traitvalsAlp,traitvalsdINF, betaI) #Compute propensity
            PropSite.append(Prop)
        Propensities_per_sites.append(PropSite)
    sumpropensities = []
    for i in Propensities_per_sites :
        sumpropensities.append(sum(i))
    return Propensities_per_sites, sumpropensities

def SumDensities(Sites) : # Count all S and I individual
    SumS = 0
    SumI = 0
    for i in Sites:
        SumS += i.effectifS
        SumI += i.effectifI
    return SumS, SumI

def GetCriticals(Propensities, Sites, Events): # Not used anymore but stored in case
    Crit_treshold = 11 #Critical number of individuals, should be between 2 and 20 (Cao et.al. 2006)
    Criticals = []
    for indexi ,i in enumerate(Events) :
        CriticalEvent = []
        for indexj,j in enumerate(Sites) :
            S, I = j.effectifS, j.effectifI  # Get the xi
            if 'epsilon' in i.formula : # Case of extinction which is always critic
                CriticalEvent.append(1)
                continue
            if i.Schange < 0 : # Case where an S individual is depleted
                if S < Crit_treshold : # If S subpop is low
                    if Propensities[indexi][indexj] > 0 : # But not zero
                        CriticalEvent.append(1) # Its critical
                    else: CriticalEvent.append(0) # Its not critical
                else : CriticalEvent.append(0)
            if i.Ichange < 0 : # Case where an I individual is depleted
                if I < Crit_treshold : # If I subpop is low
                    if Propensities[indexi][indexj] > 0 : # But not zero
                        CriticalEvent.append(1) # Its critical
                    else : CriticalEvent.append(0) # Its not critical
                else: CriticalEvent.append(0)
            elif i.Schange >0 and i.Ichange == 0 : CriticalEvent.append(0) # Case reproduction S
        Criticals.append(CriticalEvent)
    return Criticals

def ComputeMuNSigma(SumPropensities , Events):
    MuS_vector = []
    MuI_vector = []
    Output_vector = []

    for index, i in enumerate(Events) :
        if i.name == 'Extinction':
            continue
        elif i.Schange < 0 : # Case where an S individual is depleted
            vij = i.Schange
            PropS = SumPropensities[index]
            MuS_vector.append(vij*PropS)
        elif i.Ichange < 0 :
            vij = i.Schange
            PropI = SumPropensities[index]
            MuI_vector.append(vij * PropI)
        elif i.Schange > 0 and i.Ichange == 0 :
            vij = 1
            PropS = SumPropensities[index]
            MuS_vector.append(vij*PropS)
    Output_vector.append(sum(MuS_vector))
    Output_vector.append(sum(MuI_vector))
    out = np.array(Output_vector)

    return out

def GetEpsilonI(SS, SI):
    #Hor_S = 3 # Highest order reaction consuming S
    #Hor_I = 2 #Idem for I
    epsilon2 = 0.03
    Output_vector = []
    gi_s = 3/2*(2+1/(SS-1)) # Given in the paper
    gi_i = 2 # same
    Output_vector.append(epsilon2/gi_s * SS)
    Output_vector.append(epsilon2/gi_i * SI)

    out = np.array(Output_vector)
    return out

def GetTauPrime(Epsis, Mus):
    Tau_candidates = []
    Tau_finalists =[]
    for i in range(len(Mus)): # Get all Tau candidates (2 per specie)
        candidates = []
        upperterm = max(Epsis[i], 1) # Numerator
        upperterm_squared = upperterm**2 # Squared

        # Cases where all infected died in the early sim, we use a very small value for mu
        # (which leads to a huge value for tau which will never be selected but allows program to continue)
        #if Mus[i] == 0 :
            #Mus[i]= 1 / pow(10, 6)
            #print('WARNING : Avoiding division by zero by doing Trickster shit')
        #Abandonned because simulation is cancelled if all infected die


        candidates.append(upperterm/ abs(Mus[i]))
        candidates.append(upperterm_squared/abs(Mus[i])) # No need to square because mu = sigma here
        Tau_candidates.append(candidates)
    #print('Tau candidates', Tau_candidates)
    for i in range(len(Tau_candidates)): # get the tau 'finalists' (1 per specie)
        Tau_finalists.append(min(Tau_candidates[i]))
    TauPrime = min(Tau_finalists) # TauPrime is the minimum of all
    #print('Tau candidats', Tau_candidates)
    #print('Tau finalistes', Tau_finalists)
    #print('and the winner is', TauPrime)

    return TauPrime

#### FUNCTIONS FOR THE GILLESPIE ALGORITHM
#### MAYBE NOT FUNCTIONAL ANYMORE, BUT NEVER USED UNLESS THE NUMBER OF SITES IS < 30

def DoDirectMethod(Propensities, Sum_propensities, exactsteps, events, sites):
    rho = Params_double_evolution.rho  # Not very convenient to put it here but it has to....

    for i in range(exactsteps):
        r1 = np.random.uniform(0,1,1) #Draw random numbers
        a0 = sum(Sum_propensities) # Overall propensity

        Tau = (1/a0) * np.log(1/r1) #Time increment
        Probas = [i / a0 for i in Sum_propensities] #Individual propensity
        NextReaction = list(np.random.multinomial(1, Probas)) # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        #print('The next reaction Index', NextReactionIndex)

        #Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        sumprop = sum(props)
        proba_site = [i/sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1)
        #print('Next places', NextPlaceIndex)


        # This part apply the effect of events in site populations
        event = events[NextReactionIndex]
        site = sites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
            nbmigrants = 1
            # Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
            SuccessfulMigrants = 0

            roll4urlife = np.random.uniform(0, 1, 1)
            if roll4urlife > Params_double_evolution.rho: SuccessfulMigrants += 1
            # Here we distribute successful migrants among neighboring sites
            # This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
            if SuccessfulMigrants == 1:
                # Determine which site will receive the dispersing individual
                receiving_sites = deepcopy(sites)  # Create a working copy of sites
                # print('receivers', receiving_sites)
                del receiving_sites[NextPlaceIndex]  # removing departure site from the copy
                # print('receivers post suppression', receiving_sites)
                site_destination = np.random.choice(receiving_sites)  # destination is a site object
                # print('The destination is', site_destination)

                # add individual to destination
                if abs(event.Schange) > 0:  # if S are dispersers
                    site_destination.effectifS += 1
                elif abs(event.Ichange) > 0:
                    site_destination.effectifI += 1
                else:
                    pass
                    #print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
        return Tau

def DoDirectMethodV2(Propensities, Sum_propensities, exactsteps, events, sites, index_sites):
    rho = Params_double_evolution.rho  # Not very convenient to put it here but it has to....
    for i in range(exactsteps):
        r1 = np.random.uniform(0, 1, 1)  # Draw random numbers
        a0 = sum(Sum_propensities)  # Overall propensity

        Tau = (1 / a0) * np.log(1 / r1)  # Time increment
        Probas = [i / a0 for i in Sum_propensities]  # Individual propensity
        NextReaction = list(np.random.multinomial(1, Probas))  # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        # print('The next reaction Index', NextReactionIndex)

        # Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        sumprop = sum(props)
        proba_site = [i / sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1)
        # print('Next places', NextPlaceIndex)

        # This part apply the effect of events in site populations
        event = events[NextReactionIndex]
        site = sites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS += event.Schange
            site.effectifI += event.Ichange
            nbmigrants = 1
            # Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
            SuccessfulMigrants = 0

            roll4urlife = np.random.uniform(0, 1, 1)
            if roll4urlife > Params_double_evolution.rho: SuccessfulMigrants += 1
            # Here we distribute successful migrants among neighboring sites
            # This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
            if SuccessfulMigrants == 1:

                sites_index = deepcopy(index_sites)  # working copy of site indexes vector
                del index_sites[
                    NextPlaceIndex]  # Drop the current site from the list cause you can't emigrate to the place from which you departed

                Index_destination = np.random.choice(index_sites)  # Get index of destination site

                # add individual to destination
                if abs(event.Schange) > 0:  # if S are dispersers
                    sites[Index_destination].effectifS += 1
                elif abs(event.Ichange) > 0:
                    sites[Index_destination].effectifI += 1
                else:
                    pass
                    # print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            # Multiply the state change in population by the number of triggers
            if event.name == 'Extinction':
                # When extinction occur we need to retrieve densities values cause they're initialized to zero otherwise in class definition
                event.Schange = -site.effectifS  # I dont really understand why attribute is not protected but it's good news
                event.Ichange = -site.effectifI
                # print('Event Schange', trigger_persite[index])
                site.effectifS += event.Schange
                # print('Effectif après extinction', Site.effectifS)
                site.effectifI += event.Ichange
            else:
                site.effectifS += event.Schange
                site.effectifI += event.Ichange
        return Tau

#### FUNCTIONS FOR THE EVOLUTIONNARY DYNAMICS OF VIRULENCE

def GetBetaI(traitvaluesAlp): # At initialization ONLY , When evolution with trade-off transmission virulence, get beta values from alpha
    IndividualValues = [] #Outlist
    if traitvaluesAlp : # If Not empty
        for i in range(len(traitvaluesAlp)): #For each value of alpha
            value = beta0 * traitvaluesAlp[i] / (1+ (traitvaluesAlp[i])) # Get an individual value of beta
            IndividualValues.append(value) # store them
    return IndividualValues



def ChooseMuteTraitValue(EvolvingTraitAlp,EvolvingTraitdINF,NbTrigger,Statechange, TraitvaluesAlp,TraitvaluesdINF, BetaI) : #Function for Value of alpha virulence
    newtraitsvaluesAlp = deepcopy(TraitvaluesAlp) # This line is probably useless, but the safer the wiser
    newtraitsvaluesdINF = deepcopy(TraitvaluesdINF)
    #print('length newtraitsvaluesAlp Before', len(newtraitsvaluesAlp))
    #print('length newtraitsvaluesdINF Before', len(newtraitsvaluesdINF))
    newBetaI = deepcopy(BetaI)
    traittomutate = ('alpha','dI') # Two traits that can mutate
    for i in range(NbTrigger):
        if newtraitsvaluesAlp and newtraitsvaluesdINF: # If the trait values vector of each vector is not empty
            if Statechange > 0:  # If it's a birth
                if EvolvingTraitAlp.TraitMutation == True and EvolvingTraitdINF.TraitMutation == True: # If Mutation is allowed in simulation for each Traits
                    probamut = 0.01 # Proba of alpha + Proba of lambda
                    roll2mutate = np.random.uniform(0,1,1)
                    if roll2mutate < probamut : #Here there is mutation
                        traitchoose = np.random.choice(traittomutate) # Choice of the mutate traits
                        #traitchoose = traittomutate
                        if traitchoose == 'alpha' :
                            SumbetaI = sum(newBetaI)  # Get total reproduction propensity
                            IndividualsProbas = [float(i / SumbetaI) for i in
                                                 newBetaI]  # Get Individual probability to reproduce
                            Reproducer = list(
                                np.random.multinomial(1, IndividualsProbas))  # Choose one ( this return a beta_i)
                            Index_reproducer = Reproducer.index(1)  # Get the index of the reproducer

                            Parent_Value = newtraitsvaluesAlp[Index_reproducer]  # We get the trait value of the parent

                            Value_postmut = float(np.random.normal(Parent_Value, 0.02, 1))
                            NewValue = round(Value_postmut, 2)

                            #ChangeMut = [-0.01,0.01]
                            #NewValue = round(float(Parent_Value + np.random.choice(ChangeMut, 1)),2)
                            if NewValue <= 0:
                                NewValue = 0.01

                            print('Nouvelle valeur alpha après mutation', Parent_Value, NewValue)
                            newtraitsvaluesAlp.append(NewValue)
                            newtraitsvaluesdINF.append(newtraitsvaluesdINF[Index_reproducer])
                            NewbetaValue = beta0 * NewValue / (NewValue + 1)
                            newBetaI.append(NewbetaValue)


                        elif traitchoose == 'dI' :
                            SumbetaI = sum(newBetaI)  # Get total reproduction propensity
                            IndividualsProbas = [float(i / SumbetaI) for i in
                                                 newBetaI]  # Get Individual probability to reproduce
                            Reproducer = list(
                                np.random.multinomial(1, IndividualsProbas))  # Choose one ( this return a beta_i)
                            Index_reproducer = Reproducer.index(1)  # Get the index of the reproducer

                            Parent_Value = newtraitsvaluesdINF[Index_reproducer]  # We get the trait value of the parent

                            Value_postmut = float(np.random.normal(Parent_Value, 0.1, 1))
                            NewValue = round(Value_postmut, 1)

                            #ChangeMut = [-0.1,0.1]
                            #NewValue = round(float(Parent_Value + np.random.choice(ChangeMut, 1)),2)
                            #if NewValue <= -0.5:
                            #    NewValue = -0.49

                            if NewValue <= 0:
                                NewValue=0.1

                            print('Nouvelle valeur de dI après mutation', Parent_Value, NewValue)
                            newtraitsvaluesAlp.append(newtraitsvaluesAlp[Index_reproducer])
                            newtraitsvaluesdINF.append(NewValue)
                            newBetaI.append(newBetaI[Index_reproducer])

                    else:
                        # sample the individual that gives birth
                        SumbetaI = sum(newBetaI)  # Get total reproduction propensity
                        IndividualsProbas = [float(i / SumbetaI) for i in
                                             newBetaI]  # Get Individual probability to reproduce

                        Reproducer = list(
                            np.random.multinomial(1, IndividualsProbas))  # Choose one ( this return a beta_i)
                        Index_reproducer = Reproducer.index(1)  # Get the index of the reproducer

                        newtraitsvaluesAlp.append(newtraitsvaluesAlp[Index_reproducer])
                        newtraitsvaluesdINF.append(newtraitsvaluesdINF[Index_reproducer])
                        newBetaI.append(newBetaI[Index_reproducer])

                else : # If mutation not allowed
                    # sample the individual that gives birth
                    #print('coucou',Traitvalues)
                    SumbetaI = sum(newBetaI)  # Get total reproduction propensity

                    IndividualsProbas = [float(i / SumbetaI) for i in newBetaI]  # Get Individual probability to reproduce
                    Reproducer = list(np.random.multinomial(1, IndividualsProbas))  # Choose one ( this return a beta_i)

                    Index_reproducer = Reproducer.index(1)  # Get the index of the reproducer
                    newtraitsvaluesAlp.append(newtraitsvaluesAlp[Index_reproducer])
                    newtraitsvaluesdINF.append(newtraitsvaluesdINF[Index_reproducer])
                    newBetaI.append(newBetaI[Index_reproducer])


            elif Statechange < 0:  # If it's a death
                SumTrait = sum(newtraitsvaluesAlp)  # We get the total of traits values
                IndividualsProbas = [float(i / SumTrait) for i in newtraitsvaluesAlp]  # We get individuals probabilities
                DepletedGuy = list(np.random.multinomial(1, IndividualsProbas))  # Sample the removed individual
                Index_depleted = DepletedGuy.index(1) #Which value was chosen
                #print('length newtraitsvaluesAlp',len(newtraitsvaluesAlp))
                #print('length newtraitsvaluesdINF',len(newtraitsvaluesdINF))
                #print('length newtraitsvaluesBetaI',len(newBetaI))
                #print('Index_depleted',Index_depleted)
                newtraitsvaluesAlp.pop(Index_depleted)  # And remove it
                newtraitsvaluesdINF.pop(Index_depleted) # Also for dI
                newBetaI.pop(Index_depleted) # and his corresponding beta value

        else :
            print("Manque une valeur pour Alpha ou pour dI", len(newtraitsvaluesAlp) - len(newtraitsvaluesdINF))
            break
    return newtraitsvaluesAlp, newtraitsvaluesdINF, newBetaI