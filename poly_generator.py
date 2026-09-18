from finite_field import *
from gf256 import n, P
## Polynôme générateur



def poly_gene_RS(Fq,alpha,a,b):
    """
    à utiliser avec F256, alpha un générateur de Fq*, renvoie g, élément de Fq[X]/(X**n - 1) qui est le produit des X - alpha**i avec i variant de a à b inclus
    """

    g = [un(Fq)]   #initialisation

    for i in range(a,b+1):
        # on fait les calculs dans Fq[X] parce c'est plus agréable plus on repasse dans Fq[X]/(X**n -1) (ne change pas g car deg(g) < n mais il faut le bon typage)
        alpha_i = puissance(Fq,alpha,i)
        P_deg1 = [oppose(Fq,alpha_i),un(Fq)]
        g = produit_poly(Fq,g,P_deg1)

    X_n = monome_poly(Fq,un(Fq),n)
    Q = somme_poly(Fq,[moins_un(Fq)],X_n)     # X**n - 1
    Fq_n = ["Anneau",[2],[P, alpha],[Q]]

    return norme_quotient(Fq_n,g)