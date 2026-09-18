from finite_field import *
from fourier import *

# méthode de décodage via la transformée de Fourier


def gauss_elimination_mod(Fq, A, B):
    """
    Résout AX = B sur F_q, y compris lorsque A est singulière
    A : matrice carrée n x n
    B : vecteur colonne n x 1

    Met le système sous forme échelonnée réduite. Si une ligne nulle de A fait face à un second membre non nul, le système est incompatible et une ValueError est levée. Sinon on renvoie la solution particulière obtenue en fixant à zéro les inconnues libres

    Cette tolérance aux systèmes singuliers est ce qui permet au décodeur d'accepter un nombre d'erreurs strictement inférieur à t
    """
    n = len(A)
    M = [ligne[:] for ligne in A]
    C = list(B)

    colonnes_pivots = []
    ligne = 0
    for col in range(n):
        pivot = None
        for i in range(ligne, n):
            if deep_norme(Fq, M[i][col]) != deep_norme(Fq, zero(Fq)):
                pivot = i
                break
        if pivot is None:
            continue  # colonne libre

        M[ligne], M[pivot] = M[pivot], M[ligne]
        C[ligne], C[pivot] = C[pivot], C[ligne]

        inv = inverse(Fq, M[ligne][col])
        M[ligne] = [produit(Fq, x, inv) for x in M[ligne]]
        C[ligne] = produit(Fq, C[ligne], inv)

        for i in range(n):
            if i != ligne and deep_norme(Fq, M[i][col]) != deep_norme(Fq, zero(Fq)):
                f = oppose(Fq, M[i][col])
                M[i] = [somme(Fq, M[i][k], produit(Fq, f, M[ligne][k])) for k in range(n)]
                C[i] = somme(Fq, C[i], produit(Fq, f, C[ligne]))

        colonnes_pivots.append(col)
        ligne += 1
        if ligne == n:
            break

    for i in range(ligne, n):  # lignes devenues nulles
        if deep_norme(Fq, C[i]) != deep_norme(Fq, zero(Fq)):
            raise ValueError("Système incompatible")

    X = [zero(Fq) for _ in range(n)]  # inconnues libres à zéro
    for r, col in enumerate(colonnes_pivots):
        X[col] = C[r]
    return X


def decodage_via_fourier(Fq, x_prime, alpha, t):
    """
    Décode un mot reçu via la méthode de la transformée de Fourier

    Le mot reçu est un élément de Fq[X]/(X^n-1). Doit avoir été normalisé par norme_quotient sur l'anneau quotient
    Les composantes d'indice 1 à 2t de sa transformée de Fourier ne dépendent que de l'erreur et pas du mot de code émis,
    ce qui permet de déterminer les positions et les valeurs des erreurs sans connaitre le message d'origine

    Paramètres:
    Fq: le corps de base (F256 pour les QR codes)
    x_prime: le mot reçu, élément de Fq[X]/(X^n-1)
    alpha: racine primitive n-ième de l'unité dans Fq
    t: nombre d'erreurs

    Renvoie:
    Le mot de code corrigé

    Limites:
    t doit majorer le nombre d'erreurs et ne pas dépasser la capacité du code; un dépassement n'est pas détecté

    la complexité est en O(n^2)
    """
    n = len(x_prime)
    dft = dft_finite_field(Fq, x_prime, alpha)
    matrice = []
    for i in range(t):
        matrice.append([])
        for j in range(1, t + 1):
            matrice[i].append(dft[t + j - i])
    second_membre = []
    for i in range(t):
        second_membre.append(oppose(Fq, dft[t - i]))

    # Résolution du système
    sigma = [un(Fq)] + gauss_elimination_mod(Fq, matrice, second_membre)
    # sigma est déterminé

    # Trouvons maintenant les epsilon_chapeau_i
    epsilon_chapeau = [zero(Fq) for i in range(n)]
    for i in range(1, 2 * t + 1):
        epsilon_chapeau[i] = dft[i]
    som = zero(Fq)
    for k in range(1, t + 1):
        som = somme(Fq, som, produit(Fq, epsilon_chapeau[k], sigma[k]))
    epsilon_chapeau[0] = som
    for i in range(1, n - t):
        som = zero(Fq)
        for k in range(1, t + 1):
            index = (n - i + k) % n
            som = somme(Fq, produit(Fq, epsilon_chapeau[index], sigma[k]), som)
        epsilon_chapeau[(n - i) % n] = som
    epsilon = idft_finite_field(Fq, epsilon_chapeau, alpha)
    x_corrige = []

    for i in range(len(epsilon)):
        if epsilon[i] == zero(Fq):
            x_corrige.append(x_prime[i])
        else:
            x_corrige.append(somme(Fq, x_prime[i], epsilon[i]))
    return x_corrige
