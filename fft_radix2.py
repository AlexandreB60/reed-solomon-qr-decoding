"""
FFT radix-2 sur un corps fini (Cooley-Tukey)

Ce module est indépendant du décodeur : le radix-2 exige une racine
primitive 2^k-ième de l'unité, or GF(256)* est d'ordre 255, impair, donc
aucun élément n'y a un ordre pair. La FFT ne s'applique donc pas au
corps utilisé par les QR codes, et le décodeur passe par la DFT directe
de fourier.py


"""

from finite_field import *


def fft_finite_field(Fq, a, w, n=None):
    """
    Transformée de Fourier rapide sur un corps fini

    w: racine primitive len(a)-ième de l'unité; n reservé à la récursion

    Complexité en O(n log(n)), contre O(n^2) pour la DFT directe
    """
    N = len(a)
    if N==0 or N & (N - 1):                       
        raise ValueError("la taille doit être une puissance de 2")

    if n is None:
        n = N
    if N == 1:
        return a[:]

    w_N = puissance(Fq, w, n // N)  # racine adaptée à la taille N

    pair = fft_finite_field(Fq, a[::2], w, n)
    impair = fft_finite_field(Fq, a[1::2], w, n)

    result = [zero(Fq) for _ in range(N)]
    wk = un(Fq)
    for k in range(N // 2):
        t = produit(Fq, wk, impair[k])
        result[k] = somme(Fq, pair[k], t)
        result[k + N // 2] = soustraction(Fq, pair[k], t)
        wk = produit(Fq, wk, w_N)

    return result


def ifft_finite_field(Fq, a, w):
    """
    Transformée de Fourier inverse rapide sur un corps fini

    Complexité en O(n log(n))
    """
    result = fft_finite_field(Fq, a, inverse(Fq, w))
    inv_n = inverse(Fq, len(a))
    result = [produit(Fq, r, inv_n) for r in result]
    return result
