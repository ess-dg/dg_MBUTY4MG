#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:30:05 2024

@author: francescopiscitelli
"""

import numpy as np

data = np.zeros((20,1)) 


data[0,0] = 12000
data[1,0] = 12300
data[2,0] = 12900
data[3,0] = 13800
data[4,0] = 15100
data[5,0] = 16000
data[6,0] = 19000
data[7,0] = 23000
data[8,0] = 24000
data[9,0] = 25000
data[10,0] = 26000
data[11,0] = 27000
data[12,0] = 33000
data[13,0] = 33000
data[14,0] = 33000
data[15,0] = 34000
data[16,0] = 55000
data[17,0] = 56000
data[18,0] = 57000
data[19,0] = 89000


timeWindow = 1e-6

timeWindow_ns = int(round(timeWindow*1e9))
   
TimeWindowRecursive  = int(round(timeWindow_ns*1.01))
TimeWindowMax        = int(round(TimeWindowRecursive*1.5))


data1 = np.concatenate( ( np.zeros((1,np.shape(data)[1]), dtype = 'int64'), data ), axis=0)  #add a line at top not to lose the 1st event
data1[0,0] = -2*TimeWindowMax

 # data[:,0] = np.around(data[:,0],decimals=self.resolution) #time rounded at 1us precision is 6 decimals, 7 is 100ns, etc...

deltaTime = np.diff(data1[:,0])                     #1st derivative of time 
deltaTime = np.concatenate(([0],deltaTime),axis=0) #add a zero at top to restore length of vector

clusterlogic = (np.absolute(deltaTime) <= TimeWindowRecursive) #is zero when a new cluster starts 
 
data2 = np.concatenate((data1,clusterlogic[:,None]),axis=1) #this is for debugging 

index = np.argwhere(clusterlogic == False) #find the index where a new cluster may start

NumClusters = np.shape(index)[0]

data2 = np.concatenate((data2,np.zeros((len(data2),2))),axis=1) #this is for debugging 


for kk in range(0,NumClusters-1,1):
    
    
    
    clusterq = (data1[index[kk,0]:index[kk+1,0],:]).astype(int)
    
    acceptWindow = ((clusterq[-1,0] - clusterq[0,0]) <= TimeWindowMax)  #max difference in time between first and last in cluster 
    
    # print(acceptWindow)
    
    data2[index[kk,0],2] = acceptWindow
    
    
data2[:,3] =     np.logical_and(data2[:,2],np.logical_not(data2[:,1]))