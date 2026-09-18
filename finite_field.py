
## Implémentation des corps finis

def zero(Fq):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return 0
    if m == 2:
        return 0
    corps_de_base = Fq[:(m-1)]
    d = degre_poly(corps_de_base,Fq[m-1][0])
    return [zero(corps_de_base) for k in range(d)]


def un(Fq):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return 1
    if m == 2:
        return 1
    corps_de_base = Fq[:(m-1)]
    Unite = [un(corps_de_base)]
    d = degre_poly(corps_de_base,Fq[m-1][0])
    for k in range(d-1):
        Unite.append(zero(corps_de_base))
    return Unite


def plus_petit_corps(Fq,a): #renvoie le couple du plus petit sous-corps auquel a appartient, et a exprimé dans ce sous-corps
    m = len(Fq)
    if m <= 2:
        return Fq,a
    else:
        corps_de_base = Fq[:(m-1)]
        a1 = norme_quotient(Fq,a)
        d = len(a1)
        k = 1
        while k < d and a1[k] == zero(corps_de_base):
            k += 1
        if k == d:
            return plus_petit_corps(corps_de_base,a1[0])
        else:
            return Fq,a1



def somme(Fq,a,b):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return a + b
    if m == 2:
        p = Fq[1][0]
        return (a + b)%p
    corps_de_base = Fq[:(m-1)]
    Q = Fq[m-1][0] #Fq = Fq1/(Q) où Q irréductible ou non et Fq1 = corps_de_base
    d = degre_poly(corps_de_base,Q)
    return [somme(corps_de_base,a[k],b[k]) for k in range(d)] #on ne peut pas utiliser somme_poly car ce ne sont pas les mêmes conventions de représentation

def oppose(Fq,a):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return -a
    if m == 2:
        p = Fq[1][0]
        return (-a)%p
    corps_de_base = Fq[:(m-1)]
    opps = mult_scalaire_poly(corps_de_base,moins_un(corps_de_base),a)
    return norme_quotient(Fq,opps)

def moins_un(Fq):
    return oppose(Fq,un(Fq))

def soustraction(Fq,a,b):
    return somme(Fq,a,oppose(Fq,b))

def ancien_produit(Fq,a,b): #à ne pas supprimer, utile si on a pas de générateur

    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return a*b

    if m == 2:
            p = Fq[1][0]
            return (a*b)%p

    corps_de_base = Fq[:(m-1)]
    # on sort temporairement de la représentation avec vecteurs de taille fixé
    c = produit_poly(corps_de_base,a,b)
    # print(c)
    # on revient dans la bonne convention
    return norme_quotient(Fq,c)


def produit(Fq,a,b):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return a*b

    if len(Fq[m-1]) == 1:   # Fq n'est pas munis de générateur
        return ancien_produit(Fq,a,b)

    else:               # Fq est munis d'un générateur
        if a == zero(Fq) or b == zero(Fq):
            return zero(Fq)
        i = log(Fq,a)
        j = log(Fq,b)
        nb = Fq[m-1][4]
        e = (i+j) % (nb-1)
        table_p = Fq[m-1][2]
        return table_p[e]


def puissance(Fq,a,n): #exponention rapide, renvoie a**n, n est un entier relatif
    m = len(Fq)
    if len(Fq[-1]) != 1:     # si Fq* est munis d'un générateur
        if a == Fq[-1][1]:
            return puissance_gene(Fq,n)
            # si on calcule une puissance du générateur (la nouvelle fonction puissance_gene fait le même travail, mais puisque mes camarades ont déjà avancé sur leurs codes je préfère que mes modifications ne les dérangent en rien...)
        else:
            if a == zero(Fq):
                assert n > 0
                if n == 0:
                    return un(Fq)
                else:
                    return zero(Fq)

            nb = Fq[m-1][4]
            n1 = n % (nb-1)
            i = log(Fq,a)
            e = (i*n1) % (nb-1)
            table_p = Fq[m-1][2]
            return table_p[e]
    else:
        return ancien_puissance(Fq,a,n)

def ancien_puissance(Fq,a,n):
    if n == 0:
        return un(Fq)
    if n == 1:
        return a
    if n < 0:
        return puissance(Fq,inverse(Fq,a),-n)
    k = n//2
    b = puissance(Fq,a,k)
    if n%2 == 0:
        return produit(Fq,b,b)
    else:
        c = produit(Fq,b,b)
        return produit(Fq,a,c)

def puissance_gene(Fq,n):
    assert len(Fq[-1]) != 1  # Fq est bien muni d'un géné
    nb = Fq[-1][4]
    n1 = n % (nb-1)
    table_p = Fq[-1][2]

    return table_p[n1]

def ancien_inverse(Fq,a):  # A ne pas supprimer, je m'en sers si il n'y a pas de générateur
    #je dois trouver un couple de Bézout a priori, que ce soit dans Fp ou dans Fq
    if Fq[0] == "Anneau":
        raise ValueError("Pas de division dans un anneau")
    if a == zero(Fq):
        raise ValueError("Division par zéro interdite")
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return 1/a
    if m == 2 and Fq[1][0] == 2:
        return a
    if m == 2:
        p=Fq[1][0]
        return bezout(a,p)[0]%p
    else:
        corps_de_base = Fq[:(m-1)]
        Q = Fq[m-1][0]
        U = bezout_poly(corps_de_base,a,Q)[0]
        return norme_quotient(Fq,U)

def inverse(Fq,a):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return 1/a

    if len(Fq[m-1]) == 1:
        return ancien_inverse(Fq,a)

    else:               # Fq est munis d'un générateur
        i = log(Fq,a)
        nb = Fq[m-1][4]
        e = (-i) % (nb - 1)
        table_p = Fq[m-1][2]
        return table_p[e]


def quotient(Fq,a,b):
    return produit(Fq,a,inverse(Fq,b))

def exposant(chaine): # s'utilise ainsi print("x" + exposant(str(12))), affiche x**12
    sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return chaine.translate(sup)

def affichage(Fq,a,premier_appel):    # a n'a pas besoin d'être normalisé, il sera normalisé dans le code d'affichage
    # Le troisième est un booléen (il me permet d'afficher le zéro qu'au premier appel) en pratique il faut toujours prendre True
    #je prévois pour l'inconnu d'aller dans l'ordre T,Y,Z,W avec X reservé pour les polynômes (explication de ce que je veux dire par là à détailler) - je pourrais soigner le tout avec plusieurs couleurs et plusieurs tailles d'éccriture - je propose un affichage qui ne s'appuie sur des générateurs, mais il faudrait alors envisager une deuxième version (explication de ce que je veux dire par là à détailler)
    m = len(Fq)

    if premier_appel:
        print() #saute une ligne

    if m == 1 and Fq[0] == "R":
        print(a, end =" ")
    if m == 2:
        p = Fq[1][0]
        b = a%p
        print(b,end=" ")

    if m > 2:

        a1 = deep_norme(Fq,a)
        corps_de_base = Fq[:(m-1)]
        d = degre_poly(corps_de_base,a1)
        x = aux_affichage_numerotation(Fq,m,premier_appel) # désigne X, T, ou encore Y

        if a1 == zero(corps_de_base) and premier_appel:
            print("0")

        for i in range(d+1):

            if a1[i] == un(corps_de_base):

                if i == 0 and d == 0:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                elif i == 0:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                    print("+ ", end="")
                elif i == d and d==1:
                    print(x, end="")
                elif i == d:
                    print(x + exposant(str(i)), end="")
                elif i == 1:
                    print(x, "+ ", end="")
                else:
                    print(x + exposant(str(i)), "+ ", end="")

            elif a1[i] != zero(corps_de_base):

                if i == 0:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                    print("+ ", end="")
                elif i == d and d==1:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                    print(x, end="")
                elif i == d:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                    print(x + exposant(str(i)), end="")
                elif i == 1:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                    print(x, "+ ", end="")
                else:
                    aux_affichage_entier_ou_corps(corps_de_base,a1[i],i)
                    print(x + exposant(str(i)), "+ ", end="")


def aux_affichage_numerotation(Fq,m,premier_appel): # renvoie la lettre associée au niveau de profondeur du sous-corps
    if Fq[0] == "Anneau" and premier_appel:   #affiche X comme plus grosse lettre pour un quotient et W (ou Z ou Y ou T) pour un corps
        return "X"
    if m == 2:
        return "T"
    if m == 3:
        return "T"
    if m == 4:
        return "Y"
    if m == 5:
        return "Z"
    if m == 6:
        return "W"


def deep_norme(Fq,a):
    m = len(Fq)
    if m == 1 and Fq[0] == "R":
        return a
    if m == 2:
        p = Fq[1][0]
        return a%p
    else:
        a1 = a.copy()
        corps_de_base = Fq[:(m-1)]
        for i in range(len(a)):
            a1[i] = deep_norme(corps_de_base,a1[i])
        return norme_quotient(Fq,a1)


def cardinal(Fq):
    """si on envoie F256, il renvoie 256, si on envoie F2, il renvoie F2, fonctionne aussi"""
    assert len(Fq) > 1

    if len(Fq[-1]) != 1:    #muni d'un géné
        return Fq[-1][4]

    if len(Fq) == 2:
        return Fq[1][0]
    else:
        m = len(Fq)
        corps_de_base = Fq[:(m-1)]
        Q = Fq[m-1][0]
        return cardinal(corps_de_base)**degre_poly(corps_de_base,Q)


def suivant(Fq,x):
    """Donne le plus petit élément plus grand que x (on munit Fq d'un ordre lexicographique), pour Fp c'est de 0 à p-1 inclus"""
    m = len(Fq)
    if m == 1:
        raise ValueError("Pas de suivant dans des ensembles indénombrables !")
    if m == 2:
        p = Fq[1][0]
        return (x+1)%p

    corps_de_base = Fq[:(m-1)]
    Q = Fq[m-1][0]
    d = degre_poly(corps_de_base,Q)
    M_base = plus_grand(corps_de_base)
    i = d-1
    while x[i] == M_base:
        i -= 1
        if i == -1:
            return zero(Fq)

    x1 = x.copy()
    x1[i] = suivant(corps_de_base,x1[i])
    for k in range(i+1,d):
        x1[k] = zero(corps_de_base)

    return x1


def plus_grand(Fq):

    m = len(Fq)
    if m == 1:
        raise ValueError("Pas de plus grand dans des ensembles infinis !")
    if m == 2:
        return Fq[1][0]-1

    corps_de_base = Fq[:(m-1)]
    Q = Fq[m-1][0]
    d = degre_poly(corps_de_base,Q)
    M_base = plus_grand(corps_de_base)
    return [M_base for i in range(d)]

def table_puissance(Fq0,g):
    """g est un générateur du corps Fq* ; renvoie une liste tel que liste_puissance[i] = g**i ; ATTENTION Fq0 est le corps sans encore le générateur ni les tables !"""
    #Cela revient à faire l'ancien produit (mais pas besoin de passer par ancien_produit car c'est déjà pris en charge par le nouveau produit
    assert est_generateur(Fq0,g)
    return [puissance(Fq0,g,i) for i in range(cardinal(Fq0)-1)]

def numerotation(Fq,a):
    """donne une numérotation efficace de Fq"""
    assert len(Fq) != 1

    if len(Fq) == 2:
        p = Fq[1][0]
        return a%p

    m = len(Fq)
    corps_de_base = Fq[:(m-1)]
    Q = Fq[m-1][0]
    k = cardinal(corps_de_base)
    d = degre_poly(corps_de_base,Q)
    M_base = plus_grand(corps_de_base)
    Num = 0
    for i in range(d):
        Num += (k**(d-1-i)) * numerotation(corps_de_base,a[i])
    return Num

def reciproque_numerotation(Fq,i):
    """bijection réciproque de numérotation, renvoie l'élément a de Fq telle que numerotation(Fq,a) = i ; n'est utilisé que dans la création de table_log donc n'a pas besoin d'être efficace"""
    m = len(Fq)
    assert m != 1
    if m == 2:
        p = Fq[1][0]
        assert 0 <= i and i <= p-1
        return i

    corps_de_base = Fq[:(m-1)]
    Q = Fq[m-1][0]
    k = cardinal(corps_de_base)
    d = degre_poly(corps_de_base,Q)
    a = [False for _ in range(d)]
    for j in range(d):
        r = i // (k**(d-1-j))
        a[j] = r
        i = i % (k**(d-1-j))
    return a


def table_log(Fq0,g):
    """prend i un entier, g un générateur et renvoie une liste d'exposants telle liste[numerotation(Fq,P)] = exposant associé à i """
    assert est_generateur(Fq0,g)
    table_l = [False for _ in range(cardinal(Fq0))]

    for i in range(cardinal(Fq0)):
        num = numerotation(Fq0,puissance(Fq0,g,i))
        table_l[num] = i

    # pas de log en 0, lèvera une erreur avec log : table_l[0] = false
    return table_l

def log(Fq,a):
    """a est un élément de Fq, table est la table précaculée des logarithmes ; Fq est ici le corps complet"""
    # log(F256,1) rend 255 mais on passera outre

    table_l = Fq[-1][3]
    t = table_l[int(numerotation(Fq,a))]

    if t == False:
        raise ValueError("Log en 0 non défini")

    return table_l[int(numerotation(Fq,a))]


## Trouver un/les générateur de K*

# On teste un peu naïvement en faisant le calcul d'une classe d'un élément et on regarde sa taille. Pour optimiser le tout, on peut supprimer tous les éléments d'une classe d'un élément non générateur.
# Pour trouver tous les générateurs, soit on procède pareil, soit on en trouve un puis on détermine les entiers premiers avec q - 1 dans [1,q-1]

def calcul_classe(Fq,a):
    """a est un élément de Fq, cette fonction renvoie une liste d'éléments de Fq qui correspond à la classe de a, on commence par 1, a et on s'arrête juste avant a de sorte que sa taille soit toujours l'ordre de a"""
    if a == zero(Fq):
        raise ValueError("La classe de 0 n'a pas de sens - on travaille dans Fq*")

    x = un(Fq)
    classe = [un(Fq)]
    x = a

    while x != un(Fq):
        classe.append(x)
        x = produit(Fq,a,x)

    return classe

def est_generateur_rec(Fq,q,a,lst_non_gene):
    """a est un élément de Fq, q est le cardinal de Fq mis en argument pour optimiser et pas le recalculer à chaque fois ; lst_non_gene doit être triée ; cet élément renvoie un couple booléen,liste où la liste est une liste triée et égale à non_gene + la classe de a si non géné et on en a découvert des nouveaux, où juste non_gene si a est gene"""
    if a == un(Fq):
        raise ValueError("1 n'est pas générateur (sauf dans F2...)")

    if appartient(lst_non_gene,a):  #on en a pas découvert de nouveaux
        return False,lst_non_gene

    classe = calcul_classe(Fq,a)

    if len(classe) == q-1:
        return True, lst_non_gene
    else:
        classe.sort()
        return False, fusion_listes_triees(lst_non_gene,classe)

def est_generateur(Fq,a):
    """fait directement le test naïvement"""
    if a == un(Fq):
        raise ValueError("1 n'est pas générateur (sauf dans F2...)")

    classe = calcul_classe(Fq,a)

    return len(classe) == cardinal(Fq) - 1

def trouver_un_generateur(Fq):
    m = len(Fq)
    if m == 1:
        raise ValueError("Pas de générateur d'un corps infini !")

    if m == 2 and Fq[1][0] == 2:
            return 1

    x = un(Fq)
    est_gene,lst_non_gene = False,[]
    q = cardinal(Fq)

    while not est_gene:
        x = suivant(Fq,x)
        est_gene,lst_non_gene = est_generateur_rec(Fq,q,x,lst_non_gene)
    return x

def trouver_les_generateurs_vieille_methode(Fq):
    """renvoie une liste de générateurs de Fq"""

    m = len(Fq)

    if m == 1:
        raise ValueError("Pas de générateur d'un corps infini !")
    if m == 2 and Fq[1][0] == 2:
            return [1]
    if m == 2:
        x = 1
    else:
        x = zero(Fq)

    est_gene,lst_non_gene = False,[]
    lst_gene = []
    q = cardinal(Fq)
    M = plus_grand(Fq)

    while x != M:
        x = suivant(Fq,x)
        if x != un(Fq):
            est_gene,lst_non_gene = est_generateur_rec(Fq,q,x,lst_non_gene)
            if est_gene:
                lst_gene.append(x)

    return lst_gene


def fusion_listes_triees(A, B):
    """ c'est comme pour la fusion du tri fusion"""
    i, j = 0, 0
    result = []

    while i < len(A) and j < len(B):
        if A[i] < B[j]:
            result.append(A[i])
            i += 1
        elif A[i] > B[j]:
            result.append(B[j])
            j += 1
        else:  # A[i] == B[j], éviter les doublons
            result.append(A[i])
            i += 1
            j += 1

    # Ajouter les éléments restants
    result.extend(A[i:])
    result.extend(B[j:])

    return result

import bisect

def appartient(liste_triee, x):
    """Teste si x appartient à liste_triee en utilisant bisect."""
    index = bisect.bisect_left(liste_triee, x)
    return index < len(liste_triee) and liste_triee[index] == x

# Après réfléxion, le nombre de générateurs peut-être assez grand et il faut à chaque fois calculer 256 puissances successives (on pourrait améliorer par mémoïsation mais bon) (et d'ailleurs phi(n)/n est une suite bornée) d'où il vaut mieux en trouver une puis calculer des pcgd d'entiers c'est plus efficace

def trouver_les_generateurs(Fq):
    g = trouver_un_generateur(Fq)
    lst_gene = [g]
    q = cardinal(Fq)
    for i in range(2,q-1):
        if pgcd(i,q-1) == 1:
            lst_gene.append(puissance(Fq,g,i))
    return lst_gene



## Construction des polynômes


def norme_poly(Fq,P): # on renvoie un nouveau polynome L sans modifier P

    if P == []:
        return []

    k = len(P)-1
    while P[k] == zero(Fq):
        k -= 1
        if k == -1:
            return []
    return P[0:k+1]

def norme_quotient(Fq,P):   #P est un élément de Fq = Fq1[X]/(Q1) représenté initaliement par une liste de taille quelconque qui représente donc un polynôme de Fq1[X], il s'agit de réaliser une DE dns Fq1[X] puis de remplir avec des zéros ; cette fonction sert essentiellement pour le produit (peut-être pour l'inverse ? - je n'ai pas encore codé l'inverse).
    m = len(Fq)

    corps_de_base = Fq[:(m-1)] #Fq1
    Q1 = norme_poly(corps_de_base,Fq[m-1][0])  #Polynôme de Fq1[X] (irréductible ou non)
    d = degre_poly(corps_de_base,Q1)
    k = len(P)

# j'ai modifié un peu au pif la condition suivante
    # if k == 0 or (m >= 3 and (isinstance(P[0], int) or isinstance(P[0], float))) or len(P[0]) == degre_poly(Fq[:(m-2)],norme_poly(Fq[:(m-2)],Fq[m-2][0])):
    #     raise ValueError("Problème de typage pour norme_quotient !")     #test élémentaire de bon typage

    if d == k :    # d informations si on quotiente par un polynôme de degré n, et le nombre d'informations correspond à la taille de la liste
        return P

    P1 = norme_poly(corps_de_base,P)
    if k < d:
        for i in range(d-k):
            P1.append(zero(corps_de_base))
        return P1

    else:   # k > d
        R = division_euclidienne_poly(corps_de_base,P1,Q1)[1] #reste dans la DE de P par Q1
        return norme_quotient(Fq,R)    #on remplit par des zéros, ne va pas poser de problème car division_euclidienne renvoie déjà bien norme_poly(Q),norme_poly(R)


def deep_norme_poly(Fq,P):
    P1 = P.copy()
    for i in range(len(P)):
        P1[i] = deep_norme(Fq,P1[i])
    return norme_poly(Fq,P1)


def est_nul_poly(Fq,P):
    return norme_poly(Fq,P) == []

def degre_poly(Fq,P):
    if est_nul_poly(Fq,P):
        return - 1
    else:
        return len(norme_poly(Fq,P))-1


def somme_poly(Fq,P,Q):
    L = []
    n = degre_poly(Fq,P)
    m = degre_poly(Fq,Q)
    if n<m:
        return somme_poly(Fq,Q,P)
    else:
        k = 0
        while k <= m:
            L.append(somme(Fq,P[k],Q[k]))
            k += 1
        while k <= n:
            L.append(P[k])
            k += 1
    return norme_poly(Fq,L)

def mult_scalaire_poly(Fq,a,P): #renvoie a fois P où a scalaire de Fq et P dans Fq[X]
    if a == zero(Fq):
        return []
    else:
        L=[]
        for k in range(len(P)):
            L.append(produit(Fq,a,P[k]))
        return norme_poly(Fq,L)

def unitariser_poly(Fq,P): #renvoie le polynôme multiplié par un scalaire pour obtenir un polynôme unitaire
    x = P[degre_poly(Fq,P)]
    if x != un[Fq]:
        a = inverse(Fq,P[degre_poly(Fq,P)])
        return mult_scalaire_poly(Fq,a,P)
    else:
        return norme_poly(Fq,P)

def oppose_poly(Fq,P):
    a = oppose(Fq,un(Fq))
    return mult_scalaire_poly(Fq,a,P)

def soustraction_poly(Fq,P,Q):   #Renvoie P-Q
    moinsQ = oppose_poly(Fq,Q)
    return somme_poly(Fq,P,moinsQ)


# Optimisation possible en cherchant algo + performant
def produit_poly(Fq,P,Q):
    if est_nul_poly(Fq,P) or est_nul_poly(Fq,Q):
        return []

    L = []

    for k in range(1+degre_poly(Fq,P)+degre_poly(Fq,Q)):
        S = zero(Fq)
        for j in range(k+1):
            if j <= degre_poly(Fq,P) and k-j <= degre_poly(Fq,Q):
                S = somme(Fq,S,produit(Fq,P[j],Q[k-j]))
        L.append(S)

    return norme_poly(Fq,L)

def puissance_poly(Fq,a,n): #exponention rapide, renvoie a**n, n est un entier relatif
    if n == 0:
        return [un(Fq)]
    if n == 1:
        return a
    if n < 0:
        raise ValueError("Puissance de polynôme négative interdit")
    k = n//2
    b = puissance_poly(Fq,a,k)
    if n%2 == 0:
        return produit_poly(Fq,b,b)
    else:
        c = produit_poly(Fq,b,b)
        return produit_poly(Fq,a,c)

def monome_poly(Fq,a,n):    #Renvoie a.X**n où a dans Fq
    assert n >= 0
    P = [zero(Fq) for i in range(n)]
    P.append(a)
    return P

# Optimisation possible en cherchant algo + performant
def division_euclidienne_poly(Fq, A, B):    # A et B sont deux élémets de Fq[X]
    if est_nul_poly(Fq, B):
        raise ValueError("Division par le polynôme nul interdite")

    A1 = norme_poly(Fq, A)
    B1 = norme_poly(Fq, B)

    Q = []  # Quotient
    R = A1.copy() # Reste (copie de A1)

    while degre_poly(Fq, R) >= degre_poly(Fq, B1):
        # Calcul du terme de division
        coeff_div = quotient(Fq, R[degre_poly(Fq,R)], B1[degre_poly(Fq,B1)])
        deg_div = degre_poly(Fq, R) - degre_poly(Fq, B1)
        monome_div_Q = monome_poly(Fq,coeff_div, deg_div)
        monome_div_R = monome_poly(Fq, oppose(Fq,coeff_div), deg_div)

        # Mise à jour du quotient
        Q = somme_poly(Fq, Q, monome_div_Q)

        # Soustraction du terme correspondant de B à R
        R = somme_poly(Fq, R, produit_poly(Fq,monome_div_R,B))
    return norme_poly(Fq, Q), norme_poly(Fq, R)

def affichage_poly(Fq,P,premier_appel):   #on s'inspire grandement de l'affichage des éléments d'un corps
    if premier_appel:
        print() #saute une ligne
    P1 = deep_norme_poly(Fq,P)
    if P1 == []:
        print("0 - (polynôme nul)")
    else:
        d = degre_poly(Fq,P1)

        for i in range(d+1):
            if P1[i] != zero(Fq):
                if i == 0:
                    aux_affichage_entier_ou_corps(Fq,P1[i],i)
                    print("+ ", end="")
                elif i == 1 and d ==1:
                    print("X", end="")
                    aux_affichage_entier_ou_corps(Fq,P1[i],i)
                elif i == d:
                    aux_affichage_entier_ou_corps(Fq,P1[i],i)
                    print("X" + exposant(str(i)), end="")
                elif i == 1:
                    aux_affichage_entier_ou_corps(Fq,P1[i],i)
                    print("X", "+ ", end="")
                else:
                    aux_affichage_entier_ou_corps(Fq,P1[i],i)
                    print("X" + exposant(str(i)), "+ ", end="")

def aux_affichage_entier_ou_corps(Fq,a,i):
    """a est un élément de Fq car coefficient d'un polynôme de Fq[X] : cette fonction est une fonction auxiliaire de affichage_poly qui renvoie un entier ou un réel si Fq est un corps premier ou R, et sinon renvoie avec des parenthèses"""
    L,x = plus_petit_corps(Fq,a)
    if len(L) <= 2:
        if x != un(L) or i == 0:
            affichage(L,x,False)
    else:
        print("(",end="")
        affichage(L,x,False)
        print(") ",end="")

def bezout_poly(Fq,A,B): #renvoie un triplet (U,V,D) tq AU + BV = D et D pgcd de A et B (pas forcément unitaire)
    if est_nul_poly(Fq, B):
        return ([un(Fq)], [], A)
    Q, R = division_euclidienne_poly(Fq, A, B)
    X1, Y1, d = bezout_poly(Fq,B,R)
    return (Y1, soustraction_poly(Fq, X1, produit_poly(Fq, Q, Y1)), d)

def bezout(a,b): #renvoie un triplet (u,v,d) tq au + bv = d et d pgcd de a et b
    if b == 0:
        return (1, 0, a)
    x1, y1, d = bezout(b, (a % b))
    return (y1, (x1 - (a // b) * y1), d)

def pgcd(a,b):
    """classique pgcd récursif, a et b sont des entiers et on renvoie leur pgcd"""
    if b == 0:
        return(a)
    else:
        return pgcd(b,a%b)

