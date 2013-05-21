#!/usr/bin/python
# -*- coding: UTF-8 -*-
""" Classe Torseur (dans un plan) pour pouvoir calculer avec
Pierre Haessig — Mars 2013
"""

from __future__ import division, print_function, unicode_literals
import numpy as np
import matplotlib.pyplot as plt

from sympy import symbols

class Torseur2D(object):
    def __init__(self, resul, mom_ref, ref, ref_name=''):
        '''torseur de résultant `resul` et de moment `mom_ref`
        exprimé au point de coordonnées `ref`.
        l'argument optionnel `ref_name` peut indiquer le nom du point `ref`.
        '''
        # Résultante
        self.resul = resul
        # Moment exprimé en ref
        self.mom_ref = mom_ref
        # Point de référence
        self.ref = ref
        self.ref_name = ref_name
        
        print('Création d\'un ', end='')
        print(self)
    # end __init__()
    
    def calc_mom(self, pt):
        '''calcule le moment du torseur en un point `pt` quelconque'''
        # Formule du changement de point: M(X) = M(O) + R^OX
        # 1) Vecteur déplacement OX = X-O:
        delta_ref = (pt[0] - self.ref[0],
                     pt[1] - self.ref[1])
        # 2) Produit vectoriel R^OX
        delta_mom = self.resul[0]*delta_ref[1] - self.resul[1]*delta_ref[0]
        # 3) Nouveau moment:
        mom = self.mom_ref + delta_mom
        return mom
    # end calc_mom()
    
    def set_ref(self, new_ref, ref_name=''):
        '''change le point de référence'''
        # 1) Calcul du moment au nouveau point de référence
        self.mom_ref = self.calc_mom(new_ref)
        # 2) Mise à jour du point de référence
        self.ref = new_ref
        self.ref_name = ref_name
    # end set_ref()
    
    def __str__(self):
        res = 'Torseur 2D '
        res += 'de Résultante {0} et de Moment {1}'.format(str(self.resul),
                                                            str(self.mom_ref))
        if self.ref_name:
            res += ' exprimé au point {name}={0}'.format(str(self.ref), name=self.ref_name)
        else:
            res += ' exprimé au point {0}'.format(str(self.ref))
        
        return res.encode('utf-8')
    # end __str__
    
    ### Opération arithmétiques
    def __add__(self, other):
        '''addition de 2 torseurs: T1 + T2
        Le pt de référence utilisé est celui de T1 (opérande gauche)
        '''
        # 1) addition des Résultantes:
        resul = (self.resul[0] + other.resul[0],
                 self.resul[1] + other.resul[1])
        # 2) addition des moments au point de référence de self:
        mom_ref = self.mom_ref
        other_mom = other.calc_mom(self.ref)
        mom_ref = mom_ref + other_mom
        # 3) Création de l'objet Torseur:
        T = Torseur2D(resul, mom_ref, self.ref, self.ref_name)
        return T
    # end __add__
    
    def __neg__(self):
        '''opposé du torseur'''
        resul = (-self.resul[0], -self.resul[1])
        mom_ref = -self.mom_ref
        T = Torseur2D(resul, mom_ref, self.ref, self.ref_name)
        return T
    # end __neg__
    
    def __sub__(self, other):
        '''soustraction de 2 torseurs: T1-T2'''
        return self + (-other)
    
# end Torseur2D

if __name__=='__main__':
    fxA,fyA = symbols(b'fx_A,fy_A')

    F1 = Torseur2D((1,0), 1., (0,0), ref_name='O')
    FA = Torseur2D((fxA, fyA), 0., (1,1), ref_name='A')
    print(F1+FA)

