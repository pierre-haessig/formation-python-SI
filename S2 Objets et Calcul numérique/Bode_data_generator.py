#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Génère les données expérimentales (fictives)
du relevé d'un diagramme de Bode.
(Système passe-bas du 2ème ordre)

Pierre Haessig — Mai 2013
"""

import numpy as np
import matplotlib.pyplot as plt

# Paramètres du passe-bas:
H0 = 2. # gain statique
f0 = 1e3 # [Hz] fréquence canonique
Q = 5. # facteur de qualité
print((H0, f0, Q))

filename = 'bode_data.csv'

# Points où évaluer la fonction de transfert:
N = 100
f = np.logspace(2, 4, 100)# de 100 à 10 kHz

# Tranfert to 2è ordre:
H = H0/(1 + 2j*f/(f0*Q) - (f/f0)**2)

Habs = np.abs(H)
gain = 20*np.log10(Habs) # -> dB
phase = np.angle(H)*180/np.pi

# un peu de bruit pour "faire illusion"
# bruit d'autant plus grand que le gain est petit
noise_scale = (Habs.min()/Habs)**.6 # [entre 0 et 1]
gain = gain   + np.random.normal(scale=5, size=N)*noise_scale
phase = phase + np.random.normal(scale=5, size=N)*noise_scale


# Tracé rapide pour vérifier:
fig, (ax1, ax2) = plt.subplots(2,1)
ax1.semilogx(f, gain)
ax2.semilogx(f, phase)

plt.show()


### Enregistrement des données:
print('ecriture dans le fichier "%s"' % filename)
data = np.vstack([f, gain, phase]).T
print('taille des données : %s' % str(data.shape))
# Formattage des nombre avec notation scientifique ('%e'),
# avec 4 chiffres après la virgule ('%.4e')
np.savetxt(filename, data, fmt='%.6e', delimiter=',')

