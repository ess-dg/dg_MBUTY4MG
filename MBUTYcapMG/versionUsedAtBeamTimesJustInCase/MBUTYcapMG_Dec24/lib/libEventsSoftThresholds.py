#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:04:54 2021

@author: francescopiscitelli
"""


import numpy as np
import time
import os
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
import pandas as pd

# from lib import libSampleData as sdat
# from lib import libMapping as maps
# from lib import libCluster as clu
# from lib import libParameters as para
# from lib import libReadPcapngVMM as pcapr
# from lib import libHistograms as hh
# from lib import libFileManagmentUtil as fd
# from lib import libParameters as para
# from lib import libTerminal as ta
# from lib import libPlotting as plo
# from lib import libAbsUnitsAndLambda as absu


# import libSampleData as sdat
# import libMappingMG as maps
# import libCluster as clu
# import libParameters as para
# import libReadPcapngVMM as pcapr
# import libHistograms as hh
# import libFileManagmentUtil as fd
# import libParameters as para
# import libTerminal as ta
# import libPlotting as plo
# import libAbsUnitsAndLambda as absu


###############################################################################
###############################################################################

class softThresholds():
    
    def __init__(self, cassettes, parameters):
        
        self.parameters = parameters
        self.ThW = np.zeros((self.parameters.config.DETparameters.numOfWires,len(cassettes)))
        self.ThS = np.zeros((self.parameters.config.DETparameters.numOfStrips,len(cassettes)))
        self.CassetteIDs  = cassettes

###############################

class thresholdDef():
    
    def __init__(self, parameters):
        
        self.parameters = parameters

    def load(self, cassettes):
        
        if self.parameters.dataReduction.softThresholdType == 'fromFile' : 
            
            self.softThresholds = softThresholds(cassettes, self.parameters)
             
            sthfullpath = self.parameters.fileManagement.thresholdFilePath+self.parameters.fileManagement.thresholdFileName
            print('\t loading thresholds from file ...')
            
            if os.path.exists(sthfullpath) is False:
                print('\033[1;33m\t ---> WARNING ... Threshold File: '+sthfullpath+' NOT FOUND\033[1;37m')
                # print("\t ... software thresholds switched OFF ... ")
                self.parameters.dataReduction.softThresholdType = 'off'
                time.sleep(2)
                
            else:
                cassInFile  = pd.read_excel(sthfullpath).columns
                temp        = pd.read_excel(sthfullpath).values
                     
                cassNotFound = []
                cassFound    = []
                
                for k, cc in enumerate(cassettes):
               
                    if not(cc in cassInFile):
                         cassNotFound.append(cc)
                    else:
                         cassFound.append(cc)
                         index = np.where(cc == cassInFile)[0]
                         self.softThresholds.ThW[:,k] = temp[:self.parameters.config.DETparameters.numOfWires,index][:,0]
                         self.softThresholds.ThS[:,k] = temp[self.parameters.config.DETparameters.numOfWires:(self.parameters.config.DETparameters.numOfWires+self.parameters.config.DETparameters.numOfStrips),index][:,0]
                    
                if len(cassNotFound) > 0:
                    print('\t \033[1;33mWARNING ... Threshold File does NOT contain all the cassettes IDs')
                    print("\t ... software thresholds switched OFF for cassette IDs:  ",end='')
                    for cass in cassNotFound:
                        print(cass,end=' ')
                    print('\033[1;37m') 
                
        elif self.parameters.dataReduction.softThresholdType == 'userDefined' :
            
            print('\t loading user defined thresholds ...')
            self.softThresholds = self.parameters.dataReduction.softThArray
            
        elif self.parameters.dataReduction.softThresholdType == 'off' :
                
            # print('\t software thresholds OFF ... ')
            self.softThresholds = softThresholds(cassettes,self.parameters)
                
        return self.parameters.dataReduction.softThresholdType
    
    def get_thresholds1Cass(self, cassette1ID):
        
        if cassette1ID in self.softThresholds.CassetteIDs:
        
            index = self.softThresholds.CassetteIDs.index(cassette1ID)
            
            ThW = self.softThresholds.ThW[:,index]
            ThS = self.softThresholds.ThS[:,index]
            
        else:
            print("\t ... software thresholds switched OFF for cassette ID:  "+str(cassette1ID))
            ThW = np.zeros((self.parameters.config.DETparameters.numOfWires))
            ThS = np.zeros((self.parameters.config.DETparameters.numOfStrips))
            
        return ThW, ThS    
                
##############################################################

class applyThresholdsToEvents():
    
    def __init__(self, events, cassettes, parameters, showStat):
        
        clu.checkCassIDs.checkIfRepeatedIDs(cassettes)
        
        self.eventsBT   = events
        self.parameters = parameters 
        self.showStat   = showStat
        self.cassettes  = cassettes
        
        self.events = clu.events()
        
        self.thresholds = thresholdDef(self.parameters) 
        self.parameters.dataReduction.softThresholdType = self.thresholds.load(self.cassettes)

    def thresholdize1Cass(self,cassette1ID):
        
        if self.parameters.dataReduction.softThresholdType == 'userDefined' or self.parameters.dataReduction.softThresholdType == 'fromFile' : 
            
            self.events1Cass = clu.events()
            
            selCass = self.eventsBT.Cassette == cassette1ID
            
            if np.any(selCass):

                ThW, ThS = self.thresholds.get_thresholds1Cass(cassette1ID)

                wireCh0to31Round = np.round(np.mod(self.eventsBT.positionW,self.parameters.config.DETparameters.numOfWires))
                    
                stripChRound     = np.round(self.eventsBT.positionS)         
                
                for ww in range(self.parameters.config.DETparameters.numOfWires):
                    
                    belowThW = np.logical_and(np.logical_and(selCass, wireCh0to31Round == ww) , self.eventsBT.PHW <= ThW[ww])
                    
                    self.eventsBT.multW[belowThW] = -1
                    
                for ss in range(self.parameters.config.DETparameters.numOfStrips):
                        
                    belowThS  =  np.logical_and(np.logical_and(selCass, stripChRound == ss) , self.eventsBT.PHS <= ThS[ss])
                    
                    self.eventsBT.multS[belowThS] = -1
                    
                selection = np.logical_and(np.logical_and( self.eventsBT.multW >= 0 ,  self.eventsBT.multS >= 0  ), selCass)
                
                self.events1Cass.appendSelection(self.eventsBT,selection)
                
                ################################
                if self.showStat == 'individualStat':
                    print('\n\t thresholdizing ... Cassette ID ',str(cassette1ID),' ',end='')
                    self.someIndStat(self.eventsBT,selection)
                ################################
        
        elif self.parameters.dataReduction.softThresholdType == 'off' :
            
            pass
            
        else:
            
           print('\n \033[1;31msoftware thresholds -> no valid method selected, check spelling ... \033[1;37m')
            
    def thresholdizeAllCassettes(self):
        
        if self.parameters.dataReduction.softThresholdType == 'userDefined' or self.parameters.dataReduction.softThresholdType == 'fromFile' : 
            
            print(" \t software thresholds applied ... ",end='')
            
            for cc in self.cassettes:
                
                  present = clu.checkCassIDs.checkIfPresentInEvents(self.eventsBT,cc)
                  
                  if present:
                      self.thresholdize1Cass(cc)
                      self.events.append(self.events1Cass)
                      # del self.events1Cass 

            del self.events1Cass   
             
            if self.showStat == 'globalStat':
                  self.someGlobStat(self.events)
                  
        elif self.parameters.dataReduction.softThresholdType == 'off' :
           
            print(" \t software thresholds OFF ... ",end='')
           
            self.events = self.eventsBT
 
        else:
            
            print('\n \033[1;31msoftware thresholds -> no valid method selected, check spelling ... \033[1;37m')
            
            
    def someGlobStat(self, events):
        
        NeventsNotRejAfterTh  = np.sum(events.NeventsNotRejAfterTh)
        
        selc = events.positionS >= 0 
        
        NeventsNotRej2D = np.sum(selc)
        NeventsNotRej1D = np.sum(~selc)
        
        print(" --> N of events after thresholds: %d (2D: %d, 1D: %d)" % (NeventsNotRejAfterTh,NeventsNotRej2D,NeventsNotRej1D))
        
    def someIndStat(self, events, selection):
        
        NeventsNotRejAfterTh  = np.sum(selection)
        
        selc = events.positionS[selection] >= 0 
        
        NeventsNotRej2D = np.sum(selc)
        NeventsNotRej1D = np.sum(~selc)
        
        print(' --> N of events after thresholds: ',str(NeventsNotRejAfterTh),' (2D: ',str(NeventsNotRej2D),', 1D: ',str(NeventsNotRej1D),end='  ' )
        
        # print(" --> N of events after thresholds: %d (2D: %d, 1D: %d)" % (NeventsNotRejAfterTh,NeventsNotRej2D,NeventsNotRej1D))
        # print(end=' ')
      
  


###############################################################################
###############################################################################

if __name__ == '__main__':

    plt.close("all")

    configFilePath  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/config/'+"MB300_AMOR_config.json"
    filePathD       = './'+"VMM3a_Freia.pcapng"
   
    tProfilingStart = time.time()
    parameters  = para.parameters('/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/')
    config = maps.read_json_config(configFilePath)
   
    parameters.loadConfigParameters(config)
   
    parameters.cassettes.cassettes = [1,2,3]
    
  
    parameters.fileManagement.thresholdFilePath = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/config/'
    parameters.fileManagement.thresholdFileName = 'MB300L_thresholds.xlsx'
    
    # parameters.configJsonFile.offset1stWires = 10
    
    parameters.dataReduction.softThresholdType = 'off'
    
    parameters.dataReduction.softThresholdType = 'fromFile'
    
    # parameters.dataReduction.softThresholdType = 'userDefined'
    
    parameters.dataReduction.createThArrays(parameters.cassettes.cassettes, parameters)    
    parameters.dataReduction.softThArray.ThW[:,0] = 5000
    parameters.dataReduction.softThArray.ThS[:,0] = 2000 
    parameters.dataReduction.softThArray.ThW[:,1] = 5000
    
    # parameters.configJsonFile.numOfStrips = 32   
    
    # thresholds = thresholdDef(parameters)
    # thresholds.load([45,1,2,89,67,8])
    
    # ThW, ThS = thresholds.get_thresholds(1)

    #  generate sample hits 
    Nhits = 1e4
    cassettes1 = [1,2]
      
    bb = sdat.sampleHitsMultipleCassettes(cassettes1)
    bb.generateGlob(Nhits)
    hits = bb.hits

    cc = clu.clusterHits(hits,'globalStat')
    cc.clusterizeManyCassettes(parameters.cassettes.cassettes, parameters.dataReduction.timeWindow)
    events = cc.events
    
    # vv = absu.calculateAbsUnits(events, parameters)
    # vv.calculatePositionAbsUnit()
    # vv.calculateToF()
    # vv.calculateWavelength()
    
    vv = absu.calculateAbsUnits(events, parameters)
    vv.calculatePositionAbsUnit()
    vv.calculateToFandWavelength()
    events = vv.events
    eventsArray = events.concatenateEventsInArrayForDebug()
    
    # print(parameters.dataReduction.softThresholdType)
    
    thre = applyThresholdsToEvents(events, parameters.cassettes.cassettes, parameters, 'globalStat')
    
    # thre.thresholdize1Cass(1)
    
    # print(parameters.dataReduction.softThresholdType)
    
    thre.thresholdizeAllCassettes()
    
    # print(parameters.dataReduction.softThresholdType)
    
    eventsAT  = thre.events
    
    eventsArrayAT = eventsAT.concatenateEventsInArrayForDebug()
    
    # allAxis = hh.allAxis()
    # allAxis.createAllAxis(parameters)
    
    # pp = plo.plottingEvents(eventsAT,allAxis,True)
    # # pp.plotXYToF(logScale=True, absUnits = True)
    
    # # pp.plotXLambda(logScale=False, absUnits = False)
    
    # # aa = pp.plotMultiplicity([1,2,3])
    
    # pp.plotPHS(parameters.cassettes.cassettes,parameters)
    
    # pp.plotPHScorrelation([1,2,3])