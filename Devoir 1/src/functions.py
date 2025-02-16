"""
Fichier de calcul des fonction analytique et des concentrations avec la methode des differences finies
"""


#############################
#Importation de bibliothèques
#############################
import sympy as sp
import numpy as np
import config


############################
# Importation des constantes
############################
N_TOT = config.N_TOT    # nombre de point pour la méthodes des différences finies
D_EFF = config.D_EFF    # en m^2/s
C_E = config.C_E        # en mol/m^3
D = config.D            # en m
R = config.R            # en m
S = config.S            # en mol/m^3/s


######################################
# Définition de la fonction analytique
######################################
r = sp.symbols('r')
C = sp.Function('C')
C_r = 0.25*S/D_EFF*R**2*(r**2/R**2 - 1) + C_E

# Transformation de la fonction analytique sympy en fonction traçable par matplotlib
C_analytique = sp.lambdify(r, C_r, modules=['numpy'])
x_values = np.linspace(0, R, 500)
y_values = C_analytique(x_values)

def C_analytique_N(N):
    """Fonction de transformation de la fonction analytique sympy en fonction traçable par matplotlib avec N points"""
    C_analytique = sp.lambdify(r, C_r, modules=['numpy'])
    x_Nvalues = np.linspace(0, R, N)
    y_Nvalues = C_analytique(x_Nvalues)
    return x_Nvalues,y_Nvalues


#########################################
# Calcul des concentrations pour N points
#########################################
def Concentrations(N,numCas):
    """Fonction de calcul des N concentrations en différences finies"""
    C_i = np.zeros(N)           # Vecteur des N concentrations numériques calculées C_i
    matA = np.zeros((N,N))      # Matrice A pour la résolution du système matriciel
    vectB = np.zeros(N)         # Vecteur B pour la résolution du système matriciel
    delta_r=R/(N-1)             # Pas de discrétisation
    r_i=np.linspace(0, R, N)    # Vecteur des N points r_i également espacées

    # Condition de Neumann à i = 0
    if numCas == 1:     # Cas 1
        matA[0,0] = 1
        matA[0,1] = -1
        vectB[0] = 0
    elif numCas == 2:   # Cas 2 avec la méthode de Gear
        matA[0,0] = -3
        matA[0,1] = 4
        matA[0,2] = -1
        vectB[0] = 0

    # Condition de Dirichlet à i = N
    matA[-1,-1] = 1
    vectB[-1] = C_E

    # Algorithmes differences finies pour les deux cas
    if numCas == 1:
        for i in range(1,len(C_i)-1):
            matA[i,i-1] = 1/(delta_r**2)                          # Coeff devant C_i-1
            matA[i,i] = -(2/(delta_r**2)+1/(r_i[i]*delta_r))      # Coeff devant C_i
            matA[i,i+1] = 1/(delta_r**2)+1/(r_i[i]*delta_r)       # Coeff devant C_i+1
            vectB[i] = S/D_EFF

    elif numCas == 2:
        for i in range(1,len(C_i)-1):
            matA[i,i-1] = 1/(delta_r**2)-1/(r_i[i]*2*delta_r)     # Coeff devant C_i-1
            matA[i,i] = -2/(delta_r**2)                           # Coeff devant C_i
            matA[i,i+1] = 1/(delta_r**2)+1/(r_i[i]*2*delta_r)     # Coeff devant C_i+1
            vectB[i] = S/D_EFF

    C_i = np.linalg.solve(matA,vectB)
    return C_i,delta_r,r_i
