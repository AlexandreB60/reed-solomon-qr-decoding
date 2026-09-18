
import numpy as np
from finite_field import *
from gf256 import F2, F256, alpha
from poly_generator import poly_gene_RS
#génération de QR code

def deci_vers_binaire(n):
    l = bin(n)
    if "b" in l and len(l)==9:
        l = l.replace("b","")
    elif "b" in l :
        l = l.replace("b", "0")
    return l

def text_vers_bytes(text):
    """Transfome un texte en binaire"""
    r = len(text)
    out = "" #sortie en bytes
    for i in range(r):
        l= ord(text[i])
        out = out + deci_vers_binaire(l)
    return out

def xor_liste(s,t):
    out = ""
    S = len(s)
    for i in range(S):
        a = s[i]
        b = t[i]
        if (a=="1" and not b=="1") or (not a=="1" and b=="1"):
            out = out+"1"
        else:
            out = out+"0"
    return out

def text_est_conforme(text,niv_correction):
    """Prend le texte et renvoie true s'il a une taille adéquate, false sinon"""
    r = len(text)
    return not(niv_correction =="L" and r>17) and not(niv_correction == "M" and r>14) and not(niv_correction == "Q" and r>11) and not(niv_correction == "H" and r>7)

def nombre_bits_data(niv_correction):
    """Retourne le nombre de bits de data sans correction requis pour le QR code"""
    if niv_correction == "L":
      return 19*8
    if niv_correction == "M":
        return 16 * 8
    if niv_correction == "Q":
        return 13*8
    if niv_correction == "H":
        return 9*8

def text_vers_sans_correction(text,niv_correction):
    data = ""  #chaine finale à renvoyer
    r = len(text)
    if not text_est_conforme(text,niv_correction) :
        raise TypeError("Le texte est trop long ou la correction trop élevée")
    mode = "0100"  # mot pour signifier que c'est en mode d'encodage binaire
    longueur = deci_vers_binaire(r)
    while len(longueur) != 8 :
        longueur = "0"+longueur
    data = data+mode+longueur

    text_encodé = text_vers_bytes(text)
    data = data+text_encodé
    K = len(data)
    K_voulu = nombre_bits_data(niv_correction)
    if K_voulu-K <= 4 :
        for _ in range(K_voulu-K):
            data = data + "0"
    else :
        for _ in range (4):
            data = data+"0"
    K = len(data)
    while K%8 != 0:
        data = data+"0"
        K = len(data)
    nombre_pad_bytes = (K_voulu-K)//8
    for _ in range(nombre_pad_bytes // 2):
        data = data + "1110110000010001"
    if nombre_pad_bytes % 2 != 0:
        data = data+"11101100"
    return data

def nombre_bytes_correction(niv_correction):
    if niv_correction == "L":
        return 7
    if niv_correction == "M":
        return 10
    if niv_correction == "Q":
        return 13
    if niv_correction == "H":
        return 17

def polynome_générateur_correction_data(niv_correction):
    T = nombre_bytes_correction(niv_correction)
    return deep_norme_poly(F256,poly_gene_RS(F256,alpha,0,T-1))

def data_vers_polynôme(data):
    """Renvoie le polynôme message correspondant à la data"""
    k = len(data)
    message = []  #polynôme correspondant à la data initialisé au polynôme nul
    for i in range(k//8):
        coeff = []
        for p in range(8):
            bit = int(data[8*(i+1)-(p+1)])
            coeff.append(bit)
        message.append(deep_norme(F256,coeff))
    message.reverse()
    return deep_norme_poly(F256,message)

def correction_data(data,niv_correction):
    """Renvoie les bits de correction correspondant à la data et au niveau de correction"""
    correction = ""
    t = nombre_bytes_correction(niv_correction)
    G = polynome_générateur_correction_data(niv_correction)
    M = data_vers_polynôme(data)
    X_t = monome_poly(F256,un(F256),t)
    D = produit_poly(F256,M,X_t)
    R = division_euclidienne_poly(F256,D,G)[1]
    R = deep_norme_poly(F256,R)
    R.reverse()
    r = len(R)
    for i in range(r):
        for p in range(8):
            correction = correction + str(R[i][7-p])
    return correction


def text_vers_data(text,niv_correction):
    """Prend le texte et renvoie la data complète(avec correction, longueur, mode d'encodage) à mettre dans le QR code"""
    sans_correction = text_vers_sans_correction(text,niv_correction)
    correction = correction_data(sans_correction,niv_correction)
    return sans_correction + correction


def symbole_masque(num_masque):
    """Renvoie le symbole du masque à mettre dans format info"""
    out = deci_vers_binaire(num_masque)
    r = len(out)
    while r>3:
        out = out[1:]
        r = len(out)
    return out
def symbole_niveau_correction(niv_correction):
    """Renvoie le symbole du niveau de correction à mettre dans format info"""
    if niv_correction == "L" :
        return"01"
    if niv_correction == "M":
        return"00"
    if niv_correction == "Q":
        return"11"
    if niv_correction == "H":
        return"10"

def générer_format_string_sans_correction(niv_correction,num_masque):
    sym_correction = symbole_niveau_correction(niv_correction)
    sym_masque = symbole_masque(num_masque)
    return sym_correction+sym_masque

def polynome_générateur_correction_format():
    G = monome_poly(F2,1,0)
    for i in [1,2,4,5,8,10]:
        X_i = monome_poly(F2,un(F2),i)
        G = somme_poly(F2,G,X_i)
    return G

def format_vers_polynome(format):
    """Renvoie le polynôme message correspondant à la chaine de format sans correction"""
    k = len(format)
    message = []  #polynôme correspondant à la data initialisé au polynôme nul
    for i in range(k):
        bit = int(format[i])
        message.append(bit)
    message.reverse()
    return deep_norme_poly(F2,message)

def correction_format(format):
    """Renvoie les bits de correction correspondant à chaine de format"""
    correction = ""
    G = polynome_générateur_correction_format()
    M = format_vers_polynome(format)
    X_10 = monome_poly(F2,un(F2),10)
    D = produit_poly(F2,M,X_10)
    R = division_euclidienne_poly(F2,D,G)[1]
    R = deep_norme_poly(F2,R)
    R.reverse()
    r = len(R)
    for i in range(r):
        correction = correction + str(R[i])
    while len(correction) < 10:
        correction = "0" + correction
    return correction

def générer_format_string(niv_correction,num_masque):
    format = générer_format_string_sans_correction(niv_correction,num_masque)
    correction = correction_format(format)
    totale = format+correction
    chaine_masque = "101010000010010"
    return xor_liste(totale,chaine_masque)

def ajouter_finder_pattern(QR,l,k):
    """Ajoute un finder pattern avec comme coin gauche l'indice i,j """
    for j in range(k,k+7):
        QR[l,j] = 1
        QR[l+6,j] = 1
    for i in range(l+1,l+6):
        QR[i,k]= 1
        QR[i,k+6]=1
    for j in range(k+1,k+6):
        QR[l+1,j] = 0
        QR[l+5,j]= 0
    for i in range(l + 2, l + 5):
        QR[i, k+1] = 0
        QR[i, k +5] = 0
    for j in range(k+2,k+5):
        for i in range(l+2,l+5):
            QR[i,j]=1

def ajouter_tous_finder(QR):
    ajouter_finder_pattern(QR,0,0)
    ajouter_finder_pattern(QR, 0,14)
    ajouter_finder_pattern(QR, 14,0)

def ajouter_lignes_separatrices(QR):
    for j in range(0,8):
        QR[7,j] = 0
        QR[13,j] = 0
        QR[7, 20-j] = 0
    for i in range(0,7):
        QR[i,7] = 0
        QR[20-i,7] = 0
        QR[i,13] = 0

def ajouter_timing_patterns(QR):
    for k in range(5):
        l = 1+(-1)**(k)//2
        QR[12-k,6] = l
        QR[6,8+k] = l

def retirer_timing_patterns(QR):
    for k in range(5):
        QR[12-k,6] = 0.5
        QR[6,8+k] = 0.5

def ajouter_carré_noir(QR):
    QR[13,8]=1

def générer_chemin_byte(QR,sens,k,l):
    """Renvoie les indices des modules dans l'ordre pour un byte commençant en (k,l) dans le sens emprunté (1 pour montant,-1 pour descendant)"""
    L = [] #liste des indices dans l'odre
    saut = 0
    for p in range(0,8):
        (i,j) = (k-sens*(p//2),l-p%2)
        (i_suiv,j_suiv) = (i-sens*saut,j)
        if QR[i_suiv,j_suiv] != 0.5 : #s'il y a un timing pattern sur celui qu'on s'apprête à remplir
            saut = saut+1
        L.append((i-sens*saut,j))
    return L

def générer_chemin_data(QR):
    """Renvoie les indices des modules dans l'ordre pour la data"""
    L = []
    for q in range(2):
        for p in range(3):
            L = L + générer_chemin_byte(QR,1,20-4*p,20-4*q)
        for p in range(3):
            L = L+générer_chemin_byte(QR,-1,9+4*p,18-4*q)
    for p in range(4):
        L = L + générer_chemin_byte(QR, 1, 20 - 4 * p, 12)
    L = L+générer_chemin_byte(QR,1,3,12)
    L = L+générer_chemin_byte(QR,-1,0,10)
    L = L + générer_chemin_byte(QR, -1, 4, 10)
    for p in range(3):
        L = L+générer_chemin_byte(QR,-1,9+4*p,10)
    L = L+générer_chemin_byte(QR,1,12,8)
    L = L+générer_chemin_byte(QR,-1,9,5)
    L = L+générer_chemin_byte(QR,1,12,3)
    L = L+générer_chemin_byte(QR,-1,9,1)
    return L

def générer_chemin_format():
    L = [[(8,0),(20,8)],[(8,1),(19,8)],[(8,2),(18,8)],[(8,3),(17,8)],[(8,4),(16,8)],[(8,5),(15,8)],[(8,7),(14,8)],[(8,8),(8,13)],[(7,8),(8,14)],
    [(5,8),(8,15)],[(4,8),(8,16)],[(3,8),(8,17)],[(2,8),(8,18)],[(1,8),(8,19)],[(0,8),(8,20)]]
    return L

def générer_masque(num_masque):
    M = 0.5*np.ones(shape=(21,21))
    if num_masque == 0 :
        for i in range(21):
            for j in range(21):
                if (i+j) % 2==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 1 :
        for i in range(21):
            for j in range(21):
                if i % 2==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 2:
        for i in range(21):
            for j in range(21):
                if j % 3==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 3 :
        for i in range(21):
            for j in range(21):
                if (i+j) % 3==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 4 :
        for i in range(21):
            for j in range(21):
                if (i//2+j//3) % 2==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 5 :
        for i in range(21):
            for j in range(21):
                if (i*j) % 2 + (i*j) % 3  ==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 6 :
        for i in range(21):
            for j in range(21):
                if ((i*j) % 2 + (i*j) % 3)%2  ==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    if num_masque == 7 :
        for i in range(21):
            for j in range(21):
                if ((i+j) % 2 + (i*j) % 3)%2 ==0:
                    M[i,j]=0
                else :
                    M[i,j]=1
    return M
def appliquer_masque(QR,num_masque):
    """Applique le masque au QR code, attention le masque ne s'applique que au données et donc doit être appliquer sur un QR code sans patterns"""
    M = générer_masque(num_masque)
    for i in range(21):
        for j in range(21):
            if M[i,j]==0 and QR[i,j] != 0.5 :
                QR[i,j] = 1-QR[i,j]

def remplir_data(QR,text,niv_correction):
    data_str = text_vers_data(text,niv_correction)
    data = []
    for b in data_str:
        data.append(int(b))
    chemin = générer_chemin_data(QR)
    n = len(chemin)
    for p in range(n):
        (i,j) = chemin[p]
        QR[i,j] = data[p]


def remplir_format(QR,niv_correction,num_masque):
    format = générer_format_string(niv_correction,num_masque)
    chemin = générer_chemin_format()
    r = len(chemin)
    for p in range(r):
        (i0,j0) = chemin[p][0]
        (i1,j1) = chemin[p][1]
        bit = format[p]
        QR[i0,j0] = bit
        QR[i1,j1] = bit

def générer_QR(text,niv_correction,num_masque) :
    QR=0.5*np.ones(shape=(21,21)) #vierge
    ajouter_timing_patterns(QR)
    remplir_data(QR,text,niv_correction)
    retirer_timing_patterns(QR)
    appliquer_masque(QR,num_masque)
    ajouter_timing_patterns(QR)
    ajouter_tous_finder(QR)
    ajouter_lignes_separatrices(QR)
    ajouter_carré_noir(QR)
    remplir_format(QR,niv_correction,num_masque)
    return QR