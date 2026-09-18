"""
tests de la partie transformée de Fourier et décodage

Le script prend environ 30 secondes
"""

import random
from finite_field import produit_poly, norme_quotient, deep_norme
from fourier import dft_finite_field, idft_finite_field
from decoding import decodage_via_fourier
from gf256 import K, K_n, alpha, n
from qr_code import polynome_générateur_correction_data
from fft_radix2 import fft_finite_field, ifft_finite_field

T_max = 8  # capacité de correction du code (générateur de degré 17)

# --------------------------------------------------------------------------------------------------


def vecteur_aleatoire(longueur):
    """vecteur aléatoire à coefficients dans F256"""
    return [[random.randrange(2) for _ in range(8)] for _ in range(longueur)]


def egaux(u, v):
    """egalité de deux vecteurs de F256"""
    if len(u) != len(v):
        return False
    return all(deep_norme(K, a) == deep_norme(K, b) for a, b in zip(u, v))


def mot_aleatoire():
    """encode un message aléatoire et renvoie le mot de code, de longueur n"""
    g = polynome_générateur_correction_data("H")
    message = vecteur_aleatoire(56)
    return norme_quotient(K_n, produit_poly(K, message, g))


def corrompre(mot, nb_erreurs):
    """renvoie une copie de mot avec nb_erreurs symboles modifiés, et leur positions"""
    corrompu = [c[:] for c in mot]
    positions = random.sample(range(len(mot)), nb_erreurs)
    for i in positions:
        new = [random.randrange(2) for _ in range(8)]
        while deep_norme(K, new) == deep_norme(K, corrompu[i]):
            new = [random.randrange(2) for _ in range(8)]
        corrompu[i] = new
    return corrompu, positions


# ----tests-----------------------------------------------------------------------------------------


def test_transformees():    
    """idft(dft(f)) doit redonner f"""
    for _ in range(3):
        f = vecteur_aleatoire(n)
        dft = dft_finite_field(K, f, alpha)
        idft = idft_finite_field(K, dft, alpha)
        assert egaux(idft, f), "idft(dft(f))!=f"


def test_decodage_exact():
    """avec t égal au nombre réel d'erreurs, le mot de code doit être retrouvé"""
    for nb_erreurs in range(1, T_max + 1):
        mot = mot_aleatoire()
        recu, _ = corrompre(mot, nb_erreurs)
        corrige = decodage_via_fourier(K, recu, alpha, nb_erreurs)
        assert egaux(corrige, mot), "échec du décodage"


def test_sans_erreur():
    """un mot de code intact doit rester inchangé"""
    mot = mot_aleatoire()
    corrige = decodage_via_fourier(K, [c[:] for c in mot], alpha, 0)
    assert egaux(corrige, mot), "un mot sans erreur a été modifié par le décodage"


def test_t_majorant():
    """t peut être un majorant du nombre d'erreurs"""
    for nb_erreurs in range(0, T_max + 1):
        mot = mot_aleatoire()
        recu, _ = corrompre(mot, nb_erreurs)
        corrige = decodage_via_fourier(K, recu, alpha, T_max)
        assert egaux(corrige, mot), "échec du décodage"


F17 = ["Corps", [17]]  # corps où le radix-2 est possible

# racines primitives N-ièmes de l'unité dans F17 (F17* est cyclique d'ordre 16)
# les clés sont les tailles des vecteurs N, et les valeurs des racines primitives N-ièmes de l'unité dans F17
racines_F17 = {4: 13, 8: 9, 16: 3}


def test_fft_egale_dft():
    """la FFT radix-2 doit coincider avec la DFT"""
    for taille, w in racines_F17.items():
        f = [random.randrange(17) for _ in range(taille)]
        assert fft_finite_field(F17, f, w) == dft_finite_field(F17, f, w), "FFT!=DFT"


def test_aller_retour_fft():
    """ifft(fft(f)) doit donner f, quelle que soit la taille des vecteurs"""
    for taille, w in racines_F17.items():
        f = [random.randrange(17) for _ in range(taille)]
        assert ifft_finite_field(F17, fft_finite_field(F17, f, w), w) == f, ("l'aller-retour FFT/IFFT échoue")


# --------------------------------------------------------------------------------------------------
random.seed(0)
for fonction in [test_transformees, test_decodage_exact, test_sans_erreur, test_t_majorant, test_fft_egale_dft, test_aller_retour_fft]:
    fonction()
    print("OK", fonction.__name__)
