#####   MODEL PARAMETERS FOR POPULATION DYNAMICS ####
beta0 = 1  # Infectious contact rate
bi = 2  # Per capita birth rate
mu = 0.5 # Per capita natural death rate
dI = 1.0 # Dispersal infected
k = 500  # Carrying capacity
omega = (bi -mu) / k # Strength of density dependance on natality
gamma = 2.5  # Parasite Clearance
alpha = 1.5  # Parasite Virulence
rho = 0.9 # Dispersal Cost
epsilon = 0.1 # Extinction rate