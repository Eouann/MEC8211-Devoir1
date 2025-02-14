"""
Fichier de calcul des fonction analytique, des coefficients et des concentrations pour les deux cas
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
    y_Nvalues = C_analytique(x_values)
    return x_Nvalues,y_Nvalues


##################################################
# Calcul des coefficients pour N points CAS 1 et 2
##################################################
def Coefficients(N):
    """Fonction de calcul des coefficients devant les Ci+1, Ci et Ci-1"""
    delta_r=R/(N-1)
    r_i=np.zeros(N)          # Vecteur des N points r_i également espacées
    for i in range(len(r_i)):
        r_i[i] = delta_r*i

    # Calcul des coefficients devant les Ci+1, Ci et Ci-1
    # CAS 1
    a1=np.zeros(N)       # Coefficient devant Ci+1
    b1=np.zeros(N)       # Coefficient devant Ci
    c1=np.zeros(N)       # Coefficient devant Ci-1
    for i in range(1,len(r_i)):
        a1[i]=1/delta_r**2+1/(r_i[i]*delta_r)
        b1[i]=-(2/delta_r**2+1/(r_i[i]*delta_r))
        c1[i]=1/delta_r**2

    # CAS 2
    a2=np.zeros(N)       # Coefficient devant Ci+1
    b2=np.zeros(N)       # Coefficient devant Ci
    c2=np.zeros(N)       # Coefficient devant Ci-1
    for i in range(1,len(r_i)):
        a2[i]=1/delta_r**2+1/(r_i[i]*2*delta_r)
        b2[i]=-2/delta_r**2
        c2[i]=1/delta_r**2-1/(r_i[i]*2*delta_r)
    return a1,b1,c1,a2,b2,c2,r_i,delta_r


###########################################
# Calcul des concentrations pour Ntot points
###########################################
def Concentrations(a,b,c,N,numCas):
    """Fonction de calcul des N concentrations en différences finies"""
    C_i = np.zeros(N)
    matA = np.zeros((N,N))
    vectB = np.zeros(N)

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

    # Algorithmes differences finies
    for i in range(1,len(C_i)-1):
        matA[i,i-1] = c[i]
        matA[i,i] = b[i]
        matA[i,i+1] = a[i]

        vectB[i] = S/D_EFF

    C_i = np.linalg.solve(matA,vectB)
    return C_i
