#!/usr/bin/python
# -*- coding: UTF-8 -*-
""" Analyse des efforts statiques dans un treillis
(assemblage de barres et de pivots)
Le problème est traité en 2D (dans un plan)

Pierre Haessig — Mars 2013
"""

from __future__ import division, print_function, unicode_literals
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


### 1) Construction du treillis ###

# Points du treillis
N_piv = 5
print('treillis à {} pivots (dont 2 pivots de fixation)'.format(N_piv))

# Choix des points d'accroche:
iP_A = 0 # 1er point du treillis
iP_B = 1 # 2ème point
#iP_B = N_piv-1 # dernier point

# Dimensions d'une cellule élémentaire du treillis
dx = 1.
dy = 2.

h = 0. # amplitude des variations de hauteur
def courbe_y(x, haut):
    '''courbe de "mise en forme" du treillis'''
    L = (N_piv-1)*dx # Longueur totale du treillis
    y = h*x*(L-x)/(L/2)**2 # parabole
    #y = h*np.sin(2*np.pi/L*x)  # sinus (1 période)
    #y = h*np.sqrt((L/2)**2 - (x-L/2)**2) # demi-cercle
    if haut: # Partie haute du treillis
        return y + dy
    else: # Partie basse
        return y
# end courbe_y()

# Calcul des points:
iy = 0
ix = 0

# Liste des points du treillis:
pivots = []

for i in range(N_piv):
    # Calcul des coordonnées:
    x = ix*dx
    y = courbe_y(x, iy==1)
    pivots.append((x,y))
    
    # Incrémentation des compteurs:
    iy = (iy+1) % 2
    ix = ix + 1

# Conversion en positions entières
#pivots = [(int(x),int(y)) for x,y in pivots]

# Coordonnées des points d'accroche:
P_A = pivots[iP_A]
P_B = pivots[iP_B]
print('points d\'accroche: {!s} et {!s}'.format(P_A, P_B))


# Construction des barres:
barres = []

# Accroches à l'intérieur du treillis
for i, Pi in enumerate(pivots[:-2]):
    # Chaque point s'accroche aux 2 suivants:
    barres.append((Pi,pivots[i+1]))
    barres.append((Pi,pivots[i+2]))

# Barre entre les deux derniers points:
barres.append((pivots[-2], pivots[-1]))

## Problèmes d'hyperstaticité:
# Si on détecte une barre entre P_A et P_B on l'enlève (car elle est en hyperstaticité)
barres = [(P1,P2) for (P1,P2) in barres
          if (P1,P2) != (P_A,P_B) and (P1,P2) != (P_B,P_A)]

if len(barres) > (N_piv-2)*2:
    print('Il y a une barre en trop:')
    guess = 1
    ind = raw_input('barre à enlever ({} par défaut) > '.format(guess))
    if ind.strip() == '':
        ind = int(guess)
    else:
        ind = int(ind)
    print('barre enlevée : {!s}'.format(barres.pop(ind)))

N_bar = len(barres)
print('treillis à {} barres'.format(N_bar))

### Calcul du vecteurs directeur des barres:
barres_arr = np.array(barres) # shape is (2N, 2, 2)
# "AB = OB - OA":
barre_dir = barres_arr[:,1,:] - barres_arr[:,0,:]
# Longueur des barres:
barre_l = np.sqrt((barre_dir**2).sum(axis = 1))
# normalisation:
barre_dir = barre_dir/barre_l.reshape(-1,1)


### Matrice d'incidence:
Inc_mat = np.zeros((N_piv, len(barres)), dtype=int)

for j, bj in enumerate(barres):
    P1, P2 = bj
    # Remarquons la convention de signe:
    i1 = pivots.index(P1)
    Inc_mat[i1,j] = -1 # la barre bj "quitte" P1
    i2 = pivots.index(P2)
    Inc_mat[i2,j] = +1 # la barre bj "arrive à" P2

print('Matrice d\'incidence:')
print(str(Inc_mat).replace('0','.'))

# Enlever les lignes correspondant au pivot d'accroche
piv_ind = range(N_piv)
piv_ind.remove(iP_A)
piv_ind.remove(iP_B)

# Matrice d'incidence réduite:
Inc_mat_red = Inc_mat[piv_ind, :]

### Construction du système d'équation à inverser:

# 1) Matrice A
Ax = Inc_mat_red*barre_dir[:,0]
Ay = Inc_mat_red*barre_dir[:,1]
# ou bien: Ax = np.dot(Inc_mat_red,np.diag(barre_dir[:,0]))
A = np.vstack((Ax, Ay))
# Image de la matrice:
# plt.imshow(A, interpolation='nearest')

# 2) Vecteur b : force extérieure appliquée à chaque pivot:
F_ext = np.zeros((N_piv, 2))

# Appui sur le dernier pivot:
F_ext[-1] = (0., -1) # effort vers le bas

## Pesanteur sur tous les pivots:
#F_ext[:,1] = -1/N_piv

## Appui sur la clé de voute:
#F_ext[N_piv//2] = (0,-1)

print('Force extérieure (en cartésien) appliquée à chaque pivot:')
print(F_ext.T)

# Empilement des composantes selon x et y:
b_ext = np.hstack((F_ext[piv_ind,0], F_ext[piv_ind,1]))

# 3) Résolution:  Force de traction sur les barres
trac_barres = np.linalg.solve(A,b_ext)

print('Effort de traction sur chaque barre:')
print(trac_barres.round(2))
trac_max = np.max(np.abs(trac_barres))
print(' -> effort max (en val. absolue) : {:.1f}'.format(trac_max))

# Efforts sur les points d'accroche (de la part du treillis et de l'extérieur):
resul_A = -np.inner(Inc_mat[iP_A,:]*trac_barres, barre_dir.T) + F_ext[iP_A]
resul_B = -np.inner(Inc_mat[iP_B,:]*trac_barres, barre_dir.T) + F_ext[iP_B]
print('action du treillis sur pt A: {!s} ({:.2f})'.format(resul_A, np.linalg.norm(resul_A,2)))
print('action du treillis sur pt B: {!s} ({:.2f})'.format(resul_B, np.linalg.norm(resul_B,2)))

### Tracé du treillis #########################################################

fig = plt.figure('efforts treillis', figsize=(12,5))

ax = fig.add_subplot(111, title='efforts sur le treillis '
                                '({:d} pivots, {:d} barres)'.format(N_piv, N_bar))

# Échelle pour le tracé des forces
F_scale = 0.4*barre_l.mean()/trac_max
# Couleur des efforts:
F_color = (0,0.8,0) # vert


# Colormap pour colorer les barres selon l'effort subit: rouge-bleu
col_list = [(0.9,0,0.0), (0.7,0.7,0.7), (0,0,0.9)]
rb = mpl.colors.LinearSegmentedColormap.from_list('red-blue', col_list)
cm = cm_rb
#cm = plt.cm.coolwarm_r

# TODO: essayer de tracer les forces avec des FancyArrowPatch:
# a = mpl.patches.FancyArrowPatch((0,0), (1,1), arrowstyle='->, head_width=5,head_length=10')
# ax.add_patch(a)

# Tracé des barres et des efforts
for j, bj in enumerate(barres):
    # Coordonnées des 2 pivots d'accroche:
    (x1,y1), (x2, y2) = bj
    # direction :
    uj = barre_dir[j]
    # effort de traction sur la barre bj:
    trac = trac_barres[j]
    color = cm(trac/(2*trac_max)+0.5)
    plt.plot((x1, x2), (y1, y2), '-', color = color, lw=4, zorder=1)
    
    # Tracé des efforts barre -> pivot1 et pivot 2
    plt.arrow(x1, y1, +trac*uj[0]*F_scale, +trac*uj[1]*F_scale,
              zorder=2, head_width=0.05*dx, lw=0, width=0.02*dx, color=F_color)
    plt.arrow(x2, y2, -trac*uj[0]*F_scale, -trac*uj[1]*F_scale,
              zorder=2, head_width=0.05*dx, lw=0, width=0.02*dx, color=F_color)
# end for each barre

# Couleur des pivots
piv_color = (1.,1.,1.) # blanc
piv_color_AB = (1.,1.,0.5) # jaune clair
piv_alpha = 1 # opaque


# Tracé des pivots
for i, Pi in enumerate(pivots):
    # Tracé du pivot:
    marker = 'D' if Pi in (P_A,P_B) else 'o' # marqueur Diamond 'D' ou disque 'o'
    color = piv_color_AB if Pi in (P_A,P_B) else piv_color
    plt.plot(Pi[0], Pi[1], marker, ms=8, c=color, alpha = piv_alpha, zorder=3)
    # Force extérieur s'appliquant sur le pivot
    Fi = F_ext[i]
    if Fi.any():
        plt.arrow(Pi[0], Pi[1], Fi[0]*F_scale, Fi[1]*F_scale,
                  zorder=2, head_width=0.05*dx, lw=0, width=0.02*dx, color=(1,0,0))
    if Pi in (P_A,P_B):
        F_soutien = -resul_A if Pi==P_A else -resul_B
        plt.arrow(Pi[0], Pi[1], F_soutien[0]*F_scale, F_soutien[1]*F_scale,
                  zorder=2, head_width=0.05*dx, lw=0, width=0.02*dx, color=(1,1,0))
# end for each pivot

# Limites du tracé
plt.xlim(min([x for (x,y) in pivots]) - dx*1,
         max([x for (x,y) in pivots]) + dx*1)
plt.ylim(min([y for (x,y) in pivots]) - dy*.3,
         max([y for (x,y) in pivots]) + dy*.3)

# Couleur de fond:
ax.patch.set_fc((0.9,)*3)
ax.set_aspect('equal')
plt.grid(False)


fig.tight_layout()

plt.show()
