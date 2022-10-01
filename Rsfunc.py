# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 14:42:40 2022

@author: Jungle
"""

import math

pi = 3.14
def Rs(t,mi):#t是s
    a = 4*t
    b = math.log(4*math.sqrt(pi**3*t**3)/(10**mi))
    Rs = math.sqrt(-a*b)
    return Rs

def Rs_t_mi(t_min,t_max,timestep,mi):
    
    listt = [i*timestep for i in range(t_min,t_max+1)]
    listR  = []
    for i in range(len(listt)):
        t = listt[i]
        Rsi = Rs(t,mi)
        listR.append(Rsi)
        print('t = {}s,Rs = {}m'.format(t,Rs))
    return listt,listR
# timestep = 3600
# listt,listR = Rs_t_mi(, timestep,)
# for i in range(len(listt)):    
#     print(listt[i]//timestep,end ='\t') 
# for i in range(len(listR)):
#     print(listR[i],end = '\t')