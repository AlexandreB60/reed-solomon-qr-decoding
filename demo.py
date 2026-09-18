"""
Démonstration : correction d'erreurs sur un QR code.

On génère un QR code, on l'encode en Reed-Solomon sur GF(256), on corrompt
quelques octets, puis on décode par la méthode avec la transformée de Fourier.
Trois images sont affichées : le QR code d'origine, le QR code corrompu, et celui reconstruit.
"""

import matplotlib.pyplot as plt

from finite_field import produit_poly, division_euclidienne_poly, norme_quotient
from gf256 import K, K_n, alpha
from qr_code import générer_QR, polynome_générateur_correction_data
from decoding import decodage_via_fourier


positions_erreurs = [0, 17]  # octets du mot de code que l'on corrompt


def afficher(octets, titre):
    """Réassemble une liste d'octets en grille 21 x 21 et l'affiche."""
    vecteur = []
    for i in range(len(octets)):
        for j in range(8):
            vecteur.append(octets[i][j])
    while len(vecteur) < 21 * 21:  # si les derniers octets nuls ont été retirés
        vecteur.append(0)

    bit = [[0 for i in range(21)] for j in range(21)]
    for i in range(21):
        for j in range(21):
            bit[i][j] = vecteur[21 * i + j]

    plt.imshow(bit, cmap="binary")
    plt.title(titre)
    plt.axis("off")
    plt.show()


# ---Génération du QR code et mise en octets---------------------------------------------------

QR = générer_QR("Hello!", "H", 1)

QRvecteur = []
for i in range(21):
    for j in range(21):
        QRvecteur.append(QR[i][j])

QRbit = [[0, 0, 0, 0, 0, 0, 0, 0] for i in range(56)]
for i in range(55):
    for j in range(8):
        QRbit[i][j] = QRvecteur[8 * i + j]
QRbit[-1] = [QRvecteur[-1], 0, 0, 0, 0, 0, 0, 0]

afficher(QRbit, "QR code d'origine")
# ---Encodage Reed-Solomon---------------------------------------------------------------------

g = polynome_générateur_correction_data("H")
message = QRbit
codeword = produit_poly(K, message, g)

# ---Introduction de l'erreur------------------------------------------------------------------

x_corrompu = [c[:] for c in codeword]
for i in positions_erreurs:
    x_corrompu[i] = [1, 1, 1, 1, 1, 1, 1, 1]

afficher(division_euclidienne_poly(K, x_corrompu, g)[0], "QR code corrompu")

# ---Décodage----------------------------------------------------------------------------------
t = len(positions_erreurs)  # t majore le nombre d'erreurs ; ici on lui donne le nombre exact

res_decodage = decodage_via_fourier(K, norme_quotient(K_n, x_corrompu), alpha, t)
message_retrouve = division_euclidienne_poly(K, res_decodage, g)[0]
afficher(message_retrouve, "QR code reconstruit")

# ---Vérification------------------------------------------------------------------------------

print("decodage correct :", res_decodage[: len(codeword)] == codeword)
