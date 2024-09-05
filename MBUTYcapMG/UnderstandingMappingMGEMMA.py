#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:21:00 2024

@author: francescopiscitelli
"""

import numpy as np



table = np.zeros((128,7))

table[0:64,0] = np.arange(0,64)

table[64:128,0] = np.arange(0,64)

table[64:128,1] = 1 

numWiresRow = 20 

#############################

table[:,2] = 63 - table[:,0] + 64*(1-table[:,1])   


table[:,3] = np.mod(table[:,2],numWiresRow)


table[:,4] = (numWiresRow-1) - table[:,3]


table[:,5] = np.floor_divide(table[:,2],numWiresRow)


table[:,6] = table[:,4] + table[:,5]*numWiresRow