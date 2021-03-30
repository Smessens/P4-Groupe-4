# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:17:00 2020

Examples of simple t-tests

@author: fschiltz
"""
import numpy as np
from scipy import stats

a= [9.5,11,11,10.5,11]
b= [9.5,9,10,9,8.5]

#%% t-test 1: Les moyennes des deux distributions sont elles égales?
(tvalue,pvalue)=stats.ttest_ind(a,b,equal_var=False)

#%% t-test 2: La moyenne des différences pairées des éléments est-elle égale à 0?
(tvalue2,pvalue2)=stats.ttest_rel(a,b)

#%% t-test 3: La moyenne de la première distribution est-elle égale à 10?
(tvalue3,pvalue3)=stats.ttest_1samp(a,10)
#Ici, p-value=0.1: on ne peut pas rejeter l'hypothèse nulle, donc on ne peut
#rien dire.

#%% Anova: que fait on si on a une 3eme distribution?
c= [10, 10.5,11,9,12]
(tvalue4,pvalue4)=stats.f_oneway(a,b,c)