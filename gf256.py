from finite_field import *

## Construction de F256


# On va prendre F256 = F2[X]/(P) où P = X8 + X4 + X3 + X2 + 1 (P facteur irréductible de phi_255 dans F2[X]), alors alpha = Cl(X) convient
p = 2
q = 256
F2 = ["Corps",[2]]
P = [1,0,1,1,1,0,0,0,1]
n = 255 # n = q-1 et donc r = 1 (pas besoin de passer par un surcorps de K : paticularité de Reed Solomon !) Un générateur de F256 sera aussi une racine primitive n-ième de l'unité comme n = q-1
alpha = [0,1,0,0,0,0,0,0]   #au vu de sa forme partiulière, il sera noté T dans les fonctions d'affichage
g = alpha
F2560 = ["Corps",[2],[P]]

table_p = table_puissance(F2560,g)

table_l = table_log(F2560,g)

nb = cardinal(F2560)

F256 = ["Corps",[2],[P,alpha,table_p,table_l,nb]]
K = F256


X_n = monome_poly(F256,un(F256),n)
Q = somme_poly(F256,[moins_un(F256)],X_n)       # X**n - 1

F256_n = ["Anneau",[2],[P,alpha],[Q]]                       # F256/(X**n - 1)
K_n = F256_n

