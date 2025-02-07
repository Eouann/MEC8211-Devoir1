#############################
#Importation de bibliothèques
#############################
import sympy as sp
import numpy as np
import config


############################
# Importation des constantes
############################
Ntot = config.Ntot      # nombre de point pour la méthodes des différences finies
Deff = config.Deff      # en m^2/s
Ce = config.Ce          # en mol/m^3
D = config.D            # en m
R = config.R            # en m
S = config.S            # en mol/m^3/s


######################################
# Définition de la fonction analytique
######################################
r = sp.symbols('r')
C = sp.Function('C')
C_r = 0.25*S/Deff*R**2*(r**2/R**2 - 1) + Ce

# Transformation de la fonction analytique sympy en fonction traçable par matplotlib
C_analytique = sp.lambdify(r, C_r, modules=['numpy'])
x_values = np.linspace(0, R, 500)
y_values = C_analytique(x_values)