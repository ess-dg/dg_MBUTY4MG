#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 09:59:00 2024

@author: francescopiscitelli
"""

import numpy as np
import json
import os
import sys

# import pandas as pd
# from lib import libReadPcapngVMM as pcapr
# from lib import libSampleData as sdat
# import libSampleData as sdat
# 
# import libReadPcapngVMM as pcapr


# NOTE: THIS SUPPORTS ONLY 1 MONITOR

###############################################################################
###############################################################################

class hits():            
    def __init__(self): 
        
        datype = 'int64'
        
        self.Cassette    = np.zeros((0), dtype = datype)
        self.ADC         = np.zeros((0), dtype = datype)
        self.timeStamp   = np.zeros((0), dtype = datype)
        self.WorS        = -1*np.ones((0), dtype = datype)
        self.WiresStrips = -1*np.ones((0), dtype = datype)
        self.PulseT      = np.zeros((0), dtype = datype)
        self.PrevPT      = np.zeros((0), dtype = datype)
        self.Durations   = np.zeros((0), dtype = datype)
        self.Duration    = np.zeros((1), dtype = datype)
        self.WiresStrips1= -1*np.ones((0), dtype = datype)
        self.ADC1        = np.zeros((0), dtype = datype)
        self.mult0       = np.zeros((0), dtype = datype)
        self.mult1       = np.zeros((0), dtype = datype)
        
    def importReadouts(self,readouts):    
        
        self.ADC       = readouts.ADC
        self.ADC1      = readouts.ADC1
        self.timeStamp = readouts.timeStamp
        self.PulseT    = readouts.PulseT
        self.PrevPT    = readouts.PrevPT
        self.mult0     = readouts.mult0
        self.mult1     = readouts.mult1
        
        self.Durations = readouts.Durations
        self.Duration  = np.ones((1), dtype = 'int64')*(np.sum(readouts.Durations))
        
        leng = len(readouts.Channel)
        
        # 0 if wire, 1 if strip once mapped 
        # self.WorS = np.copy(readouts.ASIC)
        self.WorS = -1*np.ones((leng), dtype = 'int64') 
        
        self.Cassette     = -1*np.ones((leng), dtype = 'int64')        

        self.WiresStrips  = -1*np.ones((leng), dtype = 'int64')
        self.WiresStrips1 = -1*np.ones((leng), dtype = 'int64')
        
    def append(self, hit):
        
        self.Cassette    = np.concatenate((self.Cassette, hit.Cassette), axis=0)
        self.ADC         = np.concatenate((self.ADC, hit.ADC), axis=0)
        self.WorS        = np.concatenate((self.WorS, hit.WorS), axis=0)
        self.timeStamp   = np.concatenate((self.timeStamp, hit.timeStamp), axis=0)
        self.WiresStrips = np.concatenate((self.WiresStrips, hit.WiresStrips), axis=0)
        self.PulseT      = np.concatenate((self.PulseT, hit.PulseT), axis=0)
        self.PrevPT      = np.concatenate((self.PrevPT, hit.PrevPT), axis=0)
        self.ADC1         = np.concatenate((self.ADC1, hit.ADC1), axis=0)
        self.WiresStrips1 = np.concatenate((self.WiresStrips1, hit.WiresStrips1), axis=0)
        self.mult0        = np.concatenate((self.mult0, hit.mult0), axis=0)
        self.mult1        = np.concatenate((self.mult1, hit.mult1), axis=0)
        
        self.Durations   = np.append(self.Durations, hit.Durations)
        self.Duration    = self.Duration+hit.Duration
             
    def concatenateHitsInArrayForDebug(self):
        
        leng = len(self.WiresStrips)
        
        hitsArray = np.zeros((leng,11),dtype = 'int64')
        
        hitsArray[:,0] = self.timeStamp
        hitsArray[:,1] = self.Cassette
        hitsArray[:,2] = self.WorS
        hitsArray[:,3] = self.WiresStrips
        hitsArray[:,4] = self.ADC
        hitsArray[:,5] = self.PulseT
        hitsArray[:,6] = self.PrevPT
        hitsArray[:,7] = self.WiresStrips1
        hitsArray[:,8] = self.ADC1
        hitsArray[:,9] = self.mult0
        hitsArray[:,10] = self.mult1

        return hitsArray
   
    
class extractHitsPortion():
    
    def extract(hits,start,stop):
        
        hitsEx = hits
    
        hitsEx.Cassette  = hits.Cassette[start:stop]
        hitsEx.ADC       = hits.ADC[start:stop]
        hitsEx.WorS      = hits.WorS[start:stop]
        hitsEx.timeStamp = hits.timeStamp[start:stop]
        hitsEx.WiresStrips = hits.WiresStrips[start:stop]
        hitsEx.PulseT    = hits.PulseT[start:stop]
        hitsEx.PrevPT    = hits.PrevPT[start:stop]
        hitsEx.ADC1         = hits.ADC1[start:stop]
        hitsEx.WiresStrips1 = hits.WiresStrips1[start:stop]
        hitsEx.mult0        = hits.mult0[start:stop]
        hitsEx.mult1        = hits.mult1[start:stop]
        
        return hitsEx
    

###############################################################################
###############################################################################   

class MONmap():
    def __init__(self):
  
        self.ID       = None
        self.RingID   = None
        self.FenID    = None
        self.hybridID = None
        self.hybridSerial = None
        self.ASICID   = None
        self.channel  = None
        
class DETmap():
    def __init__(self):
  
        self.cassettesMap = []
        
class cassMap():
    def __init__(self):
        
        self.RingID   = None
        self.FenID    = None
        self.hybridWID = None
        self.hybridSID = None
        self.hybridSerial = None
        
class channelMap():
    def __init__(self): 
         
        self.AdapterType = None
        self.WireASIC    = 0
        self.StripASIC   = 1
        
class DETparameters():
    def __init__(self):
        
        self.name     = None
        
        self.operationMode = 'empty'
        
        self.orientation = 'vertical'
        
        self.numOfWires  = 0
        self.numOfStrips = 0
        
        self.numOfCassettes = 0
        
        self.cassInConfig = []
        
        # in degrees 
        self.angularOffset  = float(0)
        
        # in mm 
        self.linearOffset1stWires  = float(0)
        
        # in mm
        self.wirePitchX  = float(0)
        self.wirePitchZ  = float(0)
        self.stripPitchY = float(0)
        
        self.wiresPerRow = 0
      
################################################        

class read_json_config():
    def __init__(self, configFile_PathAndFileName):
        
        temp =  os.path.split(configFile_PathAndFileName)
            
        self.configFilePath = temp[0]+'/'
        self.configFileName = temp[1]
        
        self.MONmap = MONmap()
        self.DETmap = DETmap()
        self.DETparameters = DETparameters()
        self.cassMap = cassMap()
        # self.channelMap = channelMap()
        
        self.debug    = False
         
        try:
            self.ff   = open(configFile_PathAndFileName,'r') 
        except:
            print('\n \033[1;31m---> Config File: ' + self.configFileName + ' not found \033[1;37m')
            print('\n ---> in folder: ' + self.configFilePath + ' \n -> exiting.')
            sys.exit()
            
        try:    
            self.conf = json.load(self.ff)
        except:
            print('\n \033[1;31m---> Error in config File: ' + self.configFileName + ' \033[1;37m',end='')
            print(' ---> common mistake: last line of Cassette2ElectronicsConfig: {"ID" : X, "Ring" : X, "Fen" : X, "Hybrid" : X} must not have comma! \n -> exiting.')
            sys.exit()
        
        self.get_allParameters()
        
        self.print_DETname()
        self.print_check_operationMode()
        
        self.check_cassetteLabelling()
              
    def __del__(self):
        try:
            self.ff.close()
        except:
            pass
        
    def dprint(self, msg):
        if self.debug:
            print("{}".format(msg))
            
    def get_allParameters(self):
        self.get_DETname()
        self.get_DETparameters()
        self.get_DETmap()
        self.get_DETcassettesInConfig()
        # self.get_channelMap()
        self.get_MONmap()
            
    def print_DETname(self):
        print('\033[1;36mConfiguration for Detector: {}\033[1;37m'.format(self.DETparameters.name))
        
    def print_check_operationMode(self):
    
        if self.DETparameters.operationMode == "normal" or self.DETparameters.operationMode == "clustered":
            print('\033[1;36mOperation Mode: {}\033[1;37m'.format(self.DETparameters.operationMode))
        else:
            print('\n\t\033[1;31mERROR: Operation mode (found {}) can only be either normal or clustered -> check config file! ---> Exiting ... \n\033[1;37m'.format(self.DETparameters.operationMode),end='') 
            sys.exit()
 
            
        
    def get_DETname(self):
        self.DETparameters.name = self.conf.get('Detector')
        return self.DETparameters.name
    
    def get_DETparameters(self):
        self.DETparameters.numOfWires  = self.conf.get('wires')
        self.DETparameters.numOfStrips = self.conf.get('strips')
        self.DETparameters.wirePitchX  = float(self.conf.get('wirePitchX_mm'))
        self.DETparameters.wirePitchZ  = float(self.conf.get('wirePitchZ_mm'))
        self.DETparameters.stripPitchY = float(self.conf.get('stripPitchY_mm'))
        self.DETparameters.angularOffset         = float(self.conf.get('angularOffset_deg'))
        self.DETparameters.linearOffset1stWires  = float(self.conf.get('linearOffset1stWires_mm'))
        self.DETparameters.numOfCassettes  = self.conf.get('cassettes')
        self.DETparameters.orientation     = self.conf.get('orientation')
        self.DETparameters.operationMode   = self.conf.get('operationMode')
        self.DETparameters.wiresPerRow  = self.conf.get('wiresPerRow')
        
    def get_DETmap(self):  
        self.DETmap.cassettesMap = self.conf.get('Cassette2ElectronicsConfig')

    def get_DETcassettesInConfig(self):
                
        for cc in self.DETmap.cassettesMap:
            ID = cc.get("ID")
            self.DETparameters.cassInConfig.append(ID)
                
        self.dprint('cassettes in config file: {}'.format(self.DETparameters.cassInConfig))   
        
    def get_cassID2RingFenHybrid(self,cassetteID):

        if not cassetteID in self.DETparameters.cassInConfig:
                self.cassMap.RingID   = None
                self.cassMap.FenID    = None
                self.cassMap.hybridWID = None
                self.cassMap.hybridSID = None
                self.cassMap.hybridSerial = None
                print('\t \033[1;33mWARNING: Cassette ID ',str(cassetteID),' not found! Skipped! This CONFIG file only contains cassettes:', end=' ')
                for cass in self.DETparameters.cassInConfig:
                    print(cass,end=' ')
                print('\033[1;37m')
                
        else:
             for cc in self.DETmap.cassettesMap:
                 if cc.get("ID") == cassetteID:
                     self.cassMap.RingID       = cc.get("Ring")
                     self.cassMap.FenID        = cc.get("Fen")
                     self.cassMap.hybridWID    = cc.get("HybridW")
                     self.cassMap.hybridSID    = cc.get("HybridS")
                     self.cassMap.hybridSerial = cc.get("HybridSerial")
                     
    def check_cassetteLabelling(self):

        numOfCass               = np.shape(self.DETparameters.cassInConfig)[0]
        numOfCassFromConfigFile = self.DETparameters.numOfCassettes
        
        if numOfCass != numOfCassFromConfigFile:
            print('\033[1;31m CONFIG FILE JSON ERROR: Num of cassettes ({}) not matching num of cassettes in list ({}) in Config file\033[1;31m'.format(numOfCass,numOfCassFromConfigFile))
            print(' \n -> exiting.')
            sys.exit()
        
        # print(self.DETparameters.cassInConfig)
      
        # sortedCass = np.sort(self.DETparameters.cassInConfig)
        
        # sortedCass = self.DETparameters.cassInConfig
        
        # print(self.DETparameters.cassInConfig)
        
        # print(sortedCass)
        
        # if sortedCass[0] != 0:
        #     print('\033[1;31m CONFIG FILE JSON ERROR: Cassettes MUST be labelled from 0 -> please correct in JSON file\033[1;31m')
        #     print(' \n -> exiting.')
        #     sys.exit()
            
                     
    # def get_channelMap(self):  
        
    #     try:
    #         temp = self.conf.get('ChannelMapping')[0]
    #     except:
    #         temp = None
        
    #     # temp = temp[0]
        
    #     if temp is not None:
    #         # self.channelMap.AdapterType = temp.get("AdapterType")
    #         self.channelMap.WireASIC    = temp.get("WireASIC")
    #         self.channelMap.StripASIC   = temp.get("StripASIC")
 
    def get_MONmap(self):
             
        mapm = self.conf.get('Monitor')
        
        if mapm is None:
            print('\t \033[1;33mWARNING: No monitor config found in json file \033[1;37m')
        
        else:
            
             # for cc in self.mapm:
                 
                 # here supports only one monitor !!!
                 
                 cc = mapm[0]
            
                 self.MONmap.ID      = cc.get("ID")
   
                 TTLtypeTemp = (cc.get("TTLtype"))
                 TTLtypeBool = (TTLtypeTemp == "True") or (TTLtypeTemp == "true") or (TTLtypeTemp == "t") or (TTLtypeTemp == "T") or (TTLtypeTemp == "yes") or (TTLtypeTemp == "1")
                 
                 # if TTLtypeBool is True:
                 #        print("is true")
                 # else:
                 #        print("is false")


                 self.MONmap.TTLtype = TTLtypeBool
                 self.MONmap.RingID  = cc.get("Ring")
                 self.MONmap.FenID   = cc.get("Fen")
                 self.MONmap.hybridID      = cc.get("Hybrid")
                 self.MONmap.hybridSerial  = cc.get("HybridSerial")
                 self.MONmap.ASICID     = cc.get("ASIC")
                 self.MONmap.channel    = cc.get("Channel")
        

###############################################################################
###############################################################################

                
class mapDetector():
    def __init__(self, readouts, config):
        
        print('\033[1;36m\nMapping Detector channels ... \033[1;37m')
    
        self.debug    = False
        
        self.readouts = readouts 
        self.config   = config

        self.hits = hits()
        self.hits.importReadouts(self.readouts)
   
           
    def initCatData(self):    # debug

        if self.debug:
            self.catData = np.zeros((len(self.readouts.Ring),10), dtype = 'int64')
            self.catData[:,0] = self.hits.Cassette
            self.catData[:,1] = self.readouts.Ring
            self.catData[:,2] = self.readouts.Fen
            self.catData[:,3] = self.readouts.hybrid
            self.catData[:,4] = self.readouts.ASIC
            self.catData[:,5] = self.readouts.Channel
            self.catData[:,6] = self.hits.WiresStrips
            self.catData[:,7] = self.hits.WorS
            self.catData[:,8] = self.readouts.Channel1
            self.catData[:,9] = self.hits.WiresStrips1
            # self.catData[:,8] = self.hits.WiresStripsGlob
            
    def dprint(self, msg):
        if self.debug:
            print("{}".format(msg))     
        
    def mapp1cass(self, cassette1ID):
        
        if self.debug:
            print(cassette1ID)
        
        self.config.get_cassID2RingFenHybrid(cassette1ID)
        
        if self.debug:
            print(self.config.cassMap.RingID,self.config.cassMap.FenID,self.config.cassMap.hybridWID,self.config.cassMap.hybridSID)
        
        RingLoc = self.readouts.Ring   == self.config.cassMap.RingID
        FenLoc  = self.readouts.Fen    == self.config.cassMap.FenID
        HyWLoc  = self.readouts.hybrid == self.config.cassMap.hybridWID
        HySLoc  = self.readouts.hybrid == self.config.cassMap.hybridSID
        
        selectionCol = np.logical_and(RingLoc,FenLoc)

        selection = np.logical_and(selectionCol , np.logical_or(HyWLoc,HySLoc))
        
        if np.any(selection):
            self.hits.Cassette[selection] = cassette1ID
            flag = True
        else:
            flag = False
            # print('\t \033[1;33m No data found for Cassette ID ',str(cassette1ID),'\033[1;37m')

        #########################################
        # MAPPING CHANNELS HERE : 
        
        ########################################## 
        # WIRES
        selectionWires  =  np.logical_and(selectionCol , HyWLoc)
        
        tempW_1 = (63 - self.readouts.Channel[selectionWires] + 64*(1-self.readouts.ASIC[selectionWires]))
        
        tempW_2 = (self.config.DETparameters.wiresPerRow-1) - np.mod( tempW_1 , self.config.DETparameters.wiresPerRow )
        
        tempW_3 = np.floor_divide( tempW_1 , self.config.DETparameters.wiresPerRow )
        
        tempW_4 = tempW_2 + tempW_3*self.config.DETparameters.wiresPerRow
        
        is_even = np.mod(tempW_4 , 2 ) == 0
        
        tempW  = np.copy(tempW_4)
        
        tempW[is_even]  = tempW_4[is_even]  + 1
        tempW[~is_even] = tempW_4[~is_even] - 1

        # tempW    = (self.readouts.Channel[selectionWires] + 64*self.readouts.ASIC[selectionWires]) 
        # after mapping wires are 0 and strips 1 
        self.hits.WorS[selectionWires] = 0 
        temp1 = self.hits.WorS[selectionWires] 
        
        acceptW  =  tempW < self.config.DETparameters.numOfWires

        if np.any(~acceptW):
             print('Warning: found Wires above {} in MG column {}'.format(self.config.DETparameters.numOfWires-1,cassette1ID))
       
             
        tempW[~acceptW] = -1
        temp1[~acceptW] = -1
        
        self.hits.WiresStrips[selectionWires] = tempW
        self.hits.WorS[selectionWires]        = temp1
        
        ############################################################### 
        # GRIDS
 
        selectionStrips = np.logical_and(selectionCol , HySLoc)
        
        tempS_1 = self.readouts.Channel[selectionStrips] + 64*self.readouts.ASIC[selectionStrips]
        
        tempS_2 = 11 - tempS_1
        
        is_evenS = np.mod(tempS_2 , 2 ) == 0
        
        tempS  = np.copy(tempS_2)
        
        tempS[is_evenS]  = tempS_2[is_evenS]  + 1
        tempS[~is_evenS] = tempS_2[~is_evenS] - 1
        
        
        # after mapping wires are 0 and strips 1 
        self.hits.WorS[selectionStrips] = 1 
        temp2 = self.hits.WorS[selectionStrips] 
        
        acceptS  =  tempS < self.config.DETparameters.numOfStrips

        if np.any(~acceptS):
              print('Warning: found Grids above {} in MG column {}'.format(self.config.DETparameters.numOfStrips-1,cassette1ID))
        
        tempS[~acceptS] = -1
        temp2[~acceptS] = -1

        self.hits.WiresStrips[selectionStrips] = tempS
        self.hits.WorS[selectionStrips]        = temp2
        
        ########################################## 
        
        return flag    
            
    
    def mappAllCass(self):
        
        # mapped according to the json file
        cassettesIDs = self.config.DETparameters.cassInConfig
        
        cassNotFound = []
        
        for cass in cassettesIDs:
            flag = self.mapp1cass(cass)
            
            if flag is False:
                cassNotFound.append(cass)
        
        if len(cassNotFound) > 0:
            print('\t \033[1;33mWARNING: Cassette IDs: ', end=' ')
            for cass in cassNotFound:
                  print(cass,end=' ')
            print(' not found in DATA file. Mapping skipped for those. \033[1;37m')        
                

    # def mapChannels(self):
        
        
    #     # if  self.config.channelMap.AdapterType is None:
    #     #     print('\t \033[1;33mWARNING: No adapter type in json config file, used reverse as default. Wires on ASIC 1 and Strips on ASIC 0', end=' ')
    #     #     self.config.channelMap.AdapterType = 'reverse'
    #     # else:    
    #     #     if self.config.channelMap.AdapterType != 'straight' and self.config.channelMap.AdapterType != 'reverse':
    #     #         print('\t \033[1;33mWARNING: Adapter type can only be straight or reverse, check json config file, used reverse as default.', end=' ')
    #     #         self.config.channelMap.AdapterType = 'reverse'
       
    #     # if self.config.channelMap.WireASIC == self.config.channelMap.StripASIC:
    #     #     print('\t \033[1;31mWARNING: Wires and strips on same ASIC in config file ---> exiting!', end=' ')
    #     #     sys.exit()
            
    #     # if self.config.channelMap.WireASIC == 1 and self.config.channelMap.StripASIC == 0 and self.config.channelMap.AdapterType == 'straight' :
    #     #     print('\t \033[1;33mWARNING: Usually with straight adapter, wires are on ASIC 0 and strips on ASIC 1', end=' ')
    #     #     sys.exit()
    #     # if self.config.channelMap.WireASIC == 0 and self.config.channelMap.StripASIC == 1 and self.config.channelMap.AdapterType == 'reverse' :
    #     #    print('\t \033[1;33mWARNING: Usually with reverse adapter, wires are on ASIC 1 and strips on ASIC 0', end=' ')
    #     #    sys.exit()
        
    #     #  here we define the adapter type if VMM1 is on wires the oprder o channels is reversed for strips 
    #     # if self.config.channelMap.WireASIC == 0 and self.config.channelMap.StripASIC == 1 :
    #     #     AdapterType = 'straight'
    #     # if self.config.channelMap.WireASIC == 1 and self.config.channelMap.StripASIC == 0 :
    #     #     AdapterType = 'reverse'
        
    #     if self.config.DETparameters.operationMode == 'normal':
            
    #             #############
    #             #  wires
    #             # sel1 = self.readouts.ASIC == self.config.channelMap.WireASIC 
    #             # sel2 = self.readouts.Channel >= 16
    #             # sel3 = self.readouts.Channel <= 47
                
    #             # selectionWires  = sel1 & sel2 & sel3
                
    #             selectionWires  = self.readouts.hybrid == self.config.cassMap.hybridWID
        
                
    #             # if AdapterType == 'straight':
    #             #     self.hits.WiresStrips[selectionWires] = 31 - (self.readouts.Channel[selectionWires] - 16)
    #             # elif AdapterType == 'reverse':   
    #             #     self.hits.WiresStrips[selectionWires] = self.readouts.Channel[selectionWires] - 16
                    
    #             self.hits.WiresStrips[selectionWires] = self.readouts.Channel[selectionWires]
        
    #             # self.hits.WiresStrips[~selectionWires] = np.ma.masked # same as np.nan for int64 instead of floats
                
    #             # after mapping wires are 0 and strips 1 
    #             self.hits.WorS[selectionWires] = 0 
                
    #             #############
    #             # strips
    #             sel1s = self.readouts.ASIC == self.config.channelMap.StripASIC 
    #             sel2s = self.readouts.Channel >= 0
    #             sel3s = self.readouts.Channel <= 63
                
    #             selectionStrips  = sel1s & sel2s & sel3s
                
    #             if AdapterType == 'straight':
    #                 self.hits.WiresStrips[selectionStrips] = self.readouts.Channel[selectionStrips]
    #             elif AdapterType == 'reverse':  
    #                 self.hits.WiresStrips[selectionStrips] = 63 - self.readouts.Channel[selectionStrips]
                    
    #             # after mapping wires are 0 and strips 1 
    #             self.hits.WorS[selectionStrips] = 1 
        
    #     #############
        
    #     # if some ch of vmm are not used then nan, it is not a wire or strip 
    #     # self.hits.WorS[~selectionWires & ~selectionStrips] = np.ma.masked # same as np.nan for int64 instead of floats
        
    #     # self.hits.WorS[~selectionWires & ~selectionStrips] = -1
        
    #     elif self.config.DETparameters.operationMode == 'clustered':
            
    #         #############
    #         #  wires
    #         if self.config.channelMap.WireASIC == 0:
    #             sel2 = self.readouts.Channel >= 16
    #             sel3 = self.readouts.Channel <= 47
    #         elif self.config.channelMap.WireASIC == 1:
    #             sel2 = self.readouts.Channel1 >= 16
    #             sel3 = self.readouts.Channel1 <= 47
            
    #         selectionWires  = sel2 & sel3
            
    #         # always for any adapter the wirestips1 is for wires !!! 
    #         if self.config.channelMap.WireASIC == 0:
    #             self.hits.WiresStrips1[selectionWires] = 31 - (self.readouts.Channel[selectionWires] - 16)
    #         elif self.config.channelMap.WireASIC == 1: 
    #             self.hits.WiresStrips1[selectionWires] = self.readouts.Channel1[selectionWires] - 16
  
    #         #############
    #         # strips
    #         if self.config.channelMap.StripASIC == 0:  
    #             sel2s = self.readouts.Channel >= 0
    #             sel3s = self.readouts.Channel <= 63
    #         elif self.config.channelMap.StripASIC == 1:
    #             sel2s = self.readouts.Channel1 >= 0
    #             sel3s = self.readouts.Channel1 <= 63
            
    #         selectionStrips  = sel2s & sel3s
            
    #         # always for any adapter the wirestips is for strips !!! 
    #         if self.config.channelMap.StripASIC == 1:
    #             self.hits.WiresStrips[selectionStrips] = self.readouts.Channel1[selectionStrips]
    #         elif self.config.channelMap.StripASIC == 0:  
    #             self.hits.WiresStrips[selectionStrips] = 63 - self.readouts.Channel[selectionStrips]
      
            
    def mapChannelsGlob(self):
        
    #     self.mapChannels()
        
    #     # # global coord
    #     # self.hits.WiresStripsGlob  = np.copy(self.hits.WiresStrips)
    #     # #  if one cassette is absent in eiither the config file or the data there will be a blank space for the wanted cassette
        
    #     # # for k, cass in enumerate(cassettesIDs):
    #     #     #  does not matter the cassette ID the cassettes are arrangted one after the other according to the order in cassettesIDs,   
        
    #     #  # in this case the order is the arrangement in the json file 
    #     # for k, cass in enumerate(self.config.presentCass): 
    #     #     selection = np.logical_and( self.hits.Cassette == cass , self.hits.WorS == 0 ) #  wires is WorS = 0
    #     #     self.hits.WiresStripsGlob[selection] = self.hits.WiresStrips[selection] + k*self.config.numOfWires
        
                  # in this case the order is the arrangement in the json file 
        for k, cass in enumerate(self.config.DETparameters.cassInConfig): 
            
            index = k
            
            # if self.config.DETparameters.operationMode == 'normal':
            selection = np.logical_and( self.hits.Cassette == cass , self.hits.WorS == 0 ) #  wires is WorS = 0
            self.hits.WiresStrips[selection] = self.hits.WiresStrips[selection] + index*self.config.DETparameters.numOfWires
            # elif self.config.DETparameters.operationMode == 'clustered':
            #     selection = self.hits.Cassette == cass
            #     if self.config.channelMap.WireASIC == 0:
            #        self.hits.WiresStrips[selection] = self.hits.WiresStrips[selection] + index*self.config.DETparameters.numOfWires
            #     elif self.config.channelMap.WireASIC == 1: 
            #        self.hits.WiresStrips1[selection] = self.hits.WiresStrips1[selection] + index*self.config.DETparameters.numOfWires
            
            #  IMPORTANT NOTE 
            #  if just add +32 every cassette in json config does not matter the ID
            # index = k
            #  if the cassette ID drives the position in the space, 1 is the bottom cassette or the most left 
            # index = cass-1
            #  if the cassette ID drives the position in the space, 0 is the bottom cassette or the most left 
            # index = cass
                
            
        
    
    # def mappAllCassAndChannels(self):
        
    #     self.mappAllCass()
    #     self.mapChannels()  
            
    def mappAllCassAndChannelsGlob(self):
        
        self.mappAllCass()
        self.mapChannelsGlob()


class mapMonitor():
    def __init__(self, readouts, config):
        
        self.readouts = readouts 
        
        self.hits = hits()
           
        self.config = config
        self.config.get_MONmap()
        
        self.flagMONfound = False
        
        RingLoc = self.readouts.Ring    == self.config.MONmap.RingID
        FenLoc  = self.readouts.Fen     == self.config.MONmap.FenID
        HyLoc   = self.readouts.hybrid  == self.config.MONmap.hybridID
        ASICLoc = self.readouts.ASIC    == self.config.MONmap.ASICID
        ChLoc   = self.readouts.Channel == self.config.MONmap.channel
    
        selection = RingLoc & FenLoc & HyLoc & ASICLoc & ChLoc
 
        if np.any(selection):
            print('\033[1;36mMapping Monitor ... \033[1;37m')
            self.hits.timeStamp = self.readouts.timeStamp[selection]
            self.hits.ADC       = self.readouts.ADC[selection]
            self.hits.PrevPT    = self.readouts.PrevPT[selection]
            self.hits.PulseT    = self.readouts.PulseT[selection]
            self.hits.Durations = self.readouts.Durations
            self.hits.Duration  = np.ones((1), dtype = 'int64')*(np.sum(self.readouts.Durations))
            self.flagMONfound = True
        
        else:
            print('\t \033[1;33mNo MONITOR data found in DATA file\033[1;37m')
            self.flagMONfound = False
            



###############################################################################
###############################################################################

if __name__ == '__main__':

   filePath  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcapMG/config/'+"MG.json"
   # filePathD = './'+"VMM3a_Freia.pcapng"

   config = read_json_config(filePath)
   
   filePath = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcapMG/data/'
   file = 'MG_2col_2clusters.pcapng'
     # file = 'sampleData_ClusteredMode.pcapng'
     
   # file = 'MG_2col_1cluster.pcapng'

   filePathAndFileName = filePath+file
   
   NSperClockTick = 11.356860963629653  #ns per tick ESS for 88.0525 MHz
   
   pcapng = pcapr.pcapng_reader(filePathAndFileName, NSperClockTick, config.MONmap.TTLtype, config.MONmap.RingID,  timeResolutionType='fine', sortByTimeStampsONOFF = False, operationMode=config.DETparameters.operationMode)

   readouts = pcapng.readouts
   readoutsArray = readouts.concatenateReadoutsInArrayForDebug()
   

   md  = mapDetector(readouts, config)
   md.mappAllCassAndChannelsGlob()
   hits = md.hits
   hitsArray  = hits.concatenateHitsInArrayForDebug()
   
   # parameters.loadConfigParameters(config)
   
   # config.check_cassetteLabelling()
   
   # name = config.get_DETname()
   # config.get_DETcassettesInConfig()
   
   # config.get_cassID2RingFenHybrid(5)
   # print(config.cassMap.RingID, config.cassMap.FenID, config.cassMap.hybridID, config.cassMap.hybridSerial)
   # # config.get_cassID2RingFenHybrid(55)
   # # print(config.RingID, config.FenID, config.hybridID)
   # # config.get_cassID2RingFenHybrid(1)
   # # print(config.RingID, config.FenID, config.hybridID)
   
   # # config.get_cassID2RingFenHybrid([1,41])
   
   # config.get_bladesInclination()
   # config.get_pitches()
   # config.get_offset1stWires()
   
   # pr = pcapr.pcap_reader(filePathD)
   #   # pr.debug = True
   #  pr.read()
   #  vmm1 = pr.readouts
   
   # vmm2 = sdat.sampleReadouts_2()
   # vmm2.fill()
   # readouts = vmm2.readouts

   
   
   # re  = mapDetector(readouts, config)
   # re.debug = True
   # re.initCatData()
   # data1 = re.catData
   

   # re.mapp1cass(1)
   # # # bb = re.mapp1cass(45)
   
   # aa  = re.mappAllCass([1,2,45,3,56])
   
   # # # cassettes = np.arange(0,34,1)
    
   # re.mappAllCassAndChannels()
   
   # re.mappAllCassAndChannelsGlob()
   
   # re.initCatData()
   # data = re.catData
   
   # re.mappAllCassAndChannelsGlob()
   # hits = re.hits
   
   # hits = re.hits
   # hitsArray = hits.concatenateHitsInArrayForDebug()
    
   # # # bb = re.mappAllCassAndChannels([345,6,25])
   
   # # aa = config.get_monitor()
   
   # bb = mapMonitor(vmm2, config)
   
   # ori = getALl