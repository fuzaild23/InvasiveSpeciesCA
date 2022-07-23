#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 20:27:58 2022

@author: fuzail
"""
#%%
### THIS SCRIPT CONTAINS THE FULL CA. Part A: No control. Part B: With control
#%% 


#PART A: No control strategy
#%%
import pandas as pd
import numpy as np 
import sys
from PyQt5 import QtCore, QtGui, QtWidgets 
#pip install hexutil
import hexutil
from collections import namedtuple
from heapq import heappush, heappop
import operator
import math
import random


#Parameter initialisation
#%%
delta_t = 1 #time step increment
r = 1.18 #species growth rate
K = 0.05 # species carrying sapacity
    
i = 0 #loop counter 1
j=0 #loop counter 2

#Import data set containing suitability scores and initial species densities
#%%
speciesData = pd.read_csv("/Users/Desktop/SpeciesData.csv")


#Note!
#%%
# The notation: 'S0','S1',... refers to the list that stores the state variable (species density) at time = 0,1,... respectively



#PART A: No control strategy. 
#3 time steps are explicitly computed as an example. Repeat for as many time steps as required.
#%%

S1=[] 
for i in (range((len(speciesData)))):
    Cell_MLscore = speciesData.loc[i]['Likelihoods_bestfeatures']
    Cell_prevState = speciesData.loc[i][-1]
    Cell_Sdiff = 0
    
    for j in speciesData.loc[i]['Neighbours'].strip('][').split(', '):
        Neighb_MLscore = speciesData.loc[int(j)]['Likelihoods_bestfeatures']
        Neighb_prevState = speciesData.loc[int(j)][-1]
        Cell_Sdiff =  Cell_Sdiff + max((Cell_MLscore*max(Neighb_prevState-Cell_prevState , 0) - 
                          Neighb_MLscore*max(Cell_prevState - Neighb_prevState , 0)) ,0) 

    Cell_Sdiff = Cell_Sdiff*K
    Cell_Sgrowth = (speciesData.loc[i]['S0'])/((speciesData.loc[i]['S0'])+(1-(speciesData.loc[i]['S0']))*math.exp(-r*speciesData.loc[i]['Likelihoods_bestfeatures']*delta_t))
    
    S1.append(Cell_Sdiff+Cell_Sgrowth)
        
speciesData['S1'] = S1



S2=[]
for i in (range((len(speciesData)))):
    Cell_MLscore = speciesData.loc[i]['Likelihoods_bestfeatures']
    Cell_prevState = speciesData.loc[i][-1]
    Cell_Sdiff = 0
    
    for j in speciesData.loc[i]['Neighbours'].strip('][').split(', '):
        Neighb_MLscore = speciesData.loc[int(j)]['Likelihoods_bestfeatures']
        Neighb_prevState = speciesData.loc[int(j)][-1]
        Cell_Sdiff = Cell_Sdiff + max((Cell_MLscore*max(Neighb_prevState-Cell_prevState , 0) - 
                          Neighb_MLscore*max(Cell_prevState - Neighb_prevState , 0)) , 0) 
            
    Cell_Sdiff = Cell_Sdiff*K
    
    Cell_Sgrowth = (speciesData.loc[i]['S1'])/((speciesData.loc[i]['S1'])+(1-(speciesData.loc[i]['S1']))*math.exp(-r*speciesData.loc[i]['Likelihoods_bestfeatures']*delta_t))

    S2.append(Cell_Sdiff+Cell_Sgrowth)
       
speciesData['S2'] = S2
   


S3=[]
for i in (range((len(speciesData)))):
    Cell_MLscore = speciesData.loc[i]['Likelihoods_bestfeatures']
    Cell_prevState = speciesData.loc[i][-1]
    Cell_Sdiff = 0
    
    for j in speciesData.loc[i]['Neighbours'].strip('][').split(', '):
        Neighb_MLscore = speciesData.loc[int(j)]['Likelihoods_bestfeatures']
        Neighb_prevState = speciesData.loc[int(j)][-1]
        CCell_Sdiff =  Cell_Sdiff + max((Cell_MLscore*max(Neighb_prevState-Cell_prevState , 0) - 
                          Neighb_MLscore*max(Cell_prevState - Neighb_prevState , 0)) , 0) 
            
    Cell_Sdiff = Cell_Sdiff*K
    
    Cell_Sgrowth = (speciesData.loc[i]['S2'])/((speciesData.loc[i]['S2'])+(1-(speciesData.loc[i]['S2']))*math.exp(-r*speciesData.loc[i]['Likelihoods_bestfeatures']*delta_t))

    S3.append(Cell_Sdiff+Cell_Sgrowth)
       
speciesData['S3'] = S3





#PART B: With control strategy
#3 time steps are explicitly computed as an example. Repeat for as many time steps as required.
#%%

#Control threshold parameter initialisation
beta = 1.5


S1_Control=[]

for i in (range((len(speciesData)))):
    Cell_MLscore = speciesData.loc[i]['Likelihoods_bestfeatures']
    Cell_prevState = speciesData.loc[i]['S0']
    Cell_Sdiff = 0
    
    for j in speciesData.loc[i]['Neighbours'].strip('][').split(', '):
        Neighb_MLscore = speciesData.loc[int(j)]['Likelihoods_bestfeatures']
        Neighb_prevState = speciesData.loc[int(j)]['S0']
        Cell_Sdiff = Cell_Sdiff + max((Cell_MLscore*max(Neighb_prevState-Cell_prevState , 0) - 
                          Neighb_MLscore*max(Cell_prevState - Neighb_prevState , 0)) , 0) 
            
    Cell_Sdiff = Cell_Sdiff*K
    Cell_Sgrowth = (speciesData.loc[i]['S0'])/((speciesData.loc[i]['S0'])+(1-(speciesData.loc[i]['S0']))*math.exp(-r*speciesData.loc[i]['Likelihoods_bestfeatures']*delta_t))
    
    if ((speciesData.loc[i]['S0'])>0.001) and ( ((Cell_Sdiff+Cell_Sgrowth) - speciesData.loc[i]['S0'])/(speciesData.loc[i]['S0']) )>= beta:
        alpha = random.randint(40 , 60)/100
        S1_Control.append((1-alpha)*(Cell_Sdiff+Cell_Sgrowth))
    else: S1_Control.append((Cell_Sdiff+Cell_Sgrowth))
    
speciesData['S1_Control'] = S1_Control
    


S2_Control=[]

for i in (range((len(speciesData)))):
    Cell_MLscore = speciesData.loc[i]['Likelihoods_bestfeatures']
    Cell_prevState = speciesData.loc[i]['S1_Control']
    Cell_Sdiff = 0
    
    for j in speciesData.loc[i]['Neighbours'].strip('][').split(', '):
        Neighb_MLscore = speciesData.loc[int(j)]['Likelihoods_bestfeatures']
        Neighb_prevState = speciesData.loc[int(j)]['S1_Control']
        Cell_Sdiff =  Cell_Sdiff + max((Cell_MLscore*max(Neighb_prevState-Cell_prevState , 0) - 
                          Neighb_MLscore*max(Cell_prevState - Neighb_prevState , 0)) , 0) 
            
    Cell_Sdiff = Cell_Sdiff*K
    Cell_Sgrowth = (speciesData.loc[i]['S1_Control'])/((speciesData.loc[i]['S1_Control'])+(1-(speciesData.loc[i]['S1_Control']))*math.exp(-r*speciesData.loc[i]['Likelihoods_bestfeatures']*delta_t))
    
    if ((speciesData.loc[i]['S1_Control'])>0.001) and ( ((Cell_Sdiff+Cell_Sgrowth) - speciesData.loc[i]['S1_Control'])/(speciesData.loc[i]['S1_Control']) )>= beta:
        alpha = random.randint(40 , 60)/100
        S2_Control.append((1-alpha)*(Cell_Sdiff+Cell_Sgrowth))
    else: S2_Control.append((Cell_Sdiff+Cell_Sgrowth))
      
speciesData['S2_Control'] = S2_Control
    
   
    
S3_Control=[]

for i in (range((len(speciesData)))):
    Cell_MLscore = speciesData.loc[i]['Likelihoods_bestfeatures']
    Cell_prevState = speciesData.loc[i]['S2_Control']
    Cell_Sdiff = 0
    
    for j in speciesData.loc[i]['Neighbours'].strip('][').split(', '):
        Neighb_MLscore = speciesData.loc[int(j)]['Likelihoods_bestfeatures']
        Neighb_prevState = speciesData.loc[int(j)]['S2_Control']
        Cell_Sdiff = Cell_Sdiff + max((Cell_MLscore*max(Neighb_prevState-Cell_prevState , 0) - 
                          Neighb_MLscore*max(Cell_prevState - Neighb_prevState , 0)) , 0) 
            
    Cell_Sdiff = Cell_Sdiff*K
    Cell_Sgrowth = (speciesData.loc[i]['S2_Control'])/((speciesData.loc[i]['S2_Control'])+(1-(speciesData.loc[i]['S2_Control']))*math.exp(-r*speciesData.loc[i]['Likelihoods_bestfeatures']*delta_t))
    
    if ((speciesData.loc[i]['S2_Control'])>0.001) and ( ((Cell_Sdiff+Cell_Sgrowth) - speciesData.loc[i]['S2_Control'])/(speciesData.loc[i]['S2_Control']) )>= beta:
        alpha = random.randint(40 , 60)/100
        S3_Control.append((1-alpha)*(Cell_Sdiff+Cell_Sgrowth))
    else: S3_Control.append((Cell_Sdiff+Cell_Sgrowth))
    
speciesData['S3_Control'] = S3_Control
































