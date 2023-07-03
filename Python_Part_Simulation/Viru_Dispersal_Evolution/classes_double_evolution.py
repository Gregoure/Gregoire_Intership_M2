### Classes file for a coevolution of trait between alpha virulence and infected dispersal

import Params_double_evolution

# Parameters extracted from param files
beta0 = Params_double_evolution.beta0  # Infectious contact rate
bi = Params_double_evolution.bi  # Per capita net growth rate
mu = Params_double_evolution.mu  # Per capita natural death rate
dS = 0.5 # Susceptible Dispersal
dI = Params_double_evolution.dI #Dispersal infected
omega = Params_double_evolution.omega  # Strength of density dependance on births
k = Params_double_evolution.k  # Carrying capacity
#d = 0.5 # Dispersal propensity
gamma = Params_double_evolution.gamma  # Parasite Clearance
alpha = Params_double_evolution.alpha  # Parasite Virulence
rho = Params_double_evolution.rho  # Dispersal Cost
epsilon = Params_double_evolution.epsilon  # Extinction rate

### SPECIFIC PYTHON OBJECTS ###

# A SITE CONTAINS INDIVIDUAL DENSITIES AND TRAITS VALUES FOR A LOCAL POPULATION

class Site(): #Site object containing (non explicit) individuals
    def __init__(self,effectifS,effectifI, *args): #First try with arg way to implement feature unsure
        self.effectifS = effectifS #S density
        self.effectifI = effectifI #I density
        self.traitvaluesAlp = [] #traitvalue Affect a trait to each individual (vector) without being individual-based, Trait = alpha
        self.traitvaluesdINF = [] # Trait = dI or lambda
        self.betaI = []
        self.Index = 0

# AN EVOLVING TRAIT IS THE TRAIT IN EVOLUTIONARY MODEL

## Evolving Trait = Alpha
class EvolvingTraitAlp():
    def __init__(self, name, IsMutation): #str, bool
        self.Traitname = name
        self.TraitMutation = IsMutation

## Evolving Trait = Infected dispersal
class EvolvingTraitdINF():
    def __init__(self, name, IsMutation): #str, bool
        self.Traitname = name
        self.TraitMutation = IsMutation



# EVENTS THAT CAN HAPPEN IN THE MODEL, WITH THEIR EFFECT ON POPULATION, THEIR PROPENSITY AND STUFF
class Event():
    def __init__(self,name, propensity, Schange, Ichange, order,EvolvingTraitAlp, EvolvingTraitdINF):
        self.name = name # Event name in letter and not in memory address, handful to identify what's happening
        self.S = 0 #Has to take density values to understand the maths
        self.I = 0
        self.formula = propensity # The unique formule (str) given by model construction
        self.propensity = eval(self.formula)#Convert string in maths instruction, very useful to externalise model building
        self.Ichange = eval(Ichange) # State Change due to event, Typically -1, 0, 1 except for extinctions
        self.Schange = eval(Schange)
        self.order = order # Reaction order, not really useful but we never know
        self.EvolvingTraitAlp = EvolvingTraitAlp
        self.EvolvingTraitdINF = EvolvingTraitdINF


    def UpdatePropensity(self, S, I, TraitValuesAlp, TraitValuesdINF, BetaI):  # Class method to compute propensities without creating new objects
        if self.EvolvingTraitAlp.Traitname in self.formula:
            self.S = S
            self.I = I
            if self.EvolvingTraitAlp.Traitname == 'alpha':  # if the evolving trait is virulence
                if self.name == 'Death I':
                    SumtraitValuesAlp = sum(TraitValuesAlp)  # Do some shit to take into account the sum of trait values instead of parameter value
                    self.propensity = SumtraitValuesAlp

                if self.name == 'Infection':
                    SumBetaI = sum(BetaI)
                    self.propensity = SumBetaI * self.S  # Proper calculation of propensity with individual beta

        elif self.EvolvingTraitdINF.Traitname in self.formula:
            self.S = S
            self.I = I
            if self.EvolvingTraitdINF.Traitname == 'dI':  # if the evolving trait is infected dispersal
                if self.name == 'Dispersal I':
                    SumtraitValuesdINF = sum(TraitValuesdINF)*dS
                    self.propensity = SumtraitValuesdINF

        else:  # Otherwise take the fixed value given
            self.S = S
            self.I = I
            self.propensity = eval(self.formula)  # Changes propensity values while keeping formula untouched

        return self.propensity






