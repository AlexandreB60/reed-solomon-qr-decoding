from finite_field import *

# Transformée de Fourier et transformée de Fourier inverse


def dft_finite_field(Fq, f, w):
    """
    Transformée de Fourier discrète sur un corps fini
    Renvoie F[j]=somme_k f[k]*w^(jk), pour j de 0 à n-1

    Fq: le corps
    f: vecteur de longueur n à coefficients dans Fq
    w: racine primitive n-ième de l'unité dans Fq

    complexité en O(n^2)
    """

    n = len(f)
    result = [zero(Fq) for _ in range(n)]
    puissances = [un(Fq)]
    for i in range(1, n):
        puissances.append(produit(Fq, puissances[-1], w))

    for j in range(n):
        sum_value = zero(Fq)
        for k in range(n):
            wj = puissances[(j * k) % n]
            sum_value = somme(Fq, sum_value, produit(Fq, f[k], wj))
        result[j] = sum_value

    return result


def idft_finite_field(Fq, f, w):
    """
    Transformée de Fourier inverse discrète sur un corps fini
    Renvoie F[j]=somme_k f[k]*w^(-jk), pour j de 0 à n-1

    Fq: le corps
    f: vecteur de longueur n à coefficients dans Fq
    w: racine primitive n-ième de l'unité dans Fq

    complexité en O(n^2)

    Remarque:
    La transformée de Fourier inverse devrait diviser par n; ici n=255 est impair et comme Fq de caractéristique 2, n=1 dans Fq d'où la division est inutile. Ce n'est pas vrai en général
    """
    n = len(f)
    result = [zero(Fq) for _ in range(n)]
    w_inv = inverse(Fq, w)
    puissances_neg = [un(Fq)]
    for i in range(1, n):
        puissances_neg.append(produit(Fq, puissances_neg[-1], w_inv))

    for j in range(n):
        sum_value = zero(Fq)
        for k in range(n):
            wj = puissances_neg[(k * j) % n]
            sum_value = somme(Fq, sum_value, produit(Fq, f[k], wj))
        result[j] = sum_value

    return result
