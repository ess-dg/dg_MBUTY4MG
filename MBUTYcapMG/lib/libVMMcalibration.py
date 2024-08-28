#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 11:05:36 2021

@author: francescopiscitelli
"""
import numpy as np
import json
import os

import sys

from lib import libMapping as maps
from lib import libReadPcapngVMM as pcapr

# import libMapping as maps
# import libReadPcapngVMM as pcapr

###############################################################################
###############################################################################

# NOTES for devel

# TDC calib not finished! 

#  add if calib file corrupted add slope 1 and offset 0 to missing values ! 


###############################################################################
###############################################################################

class calibrationMapVMM():
    def __init__(self):
        
        self.ADCslope  = np.ones((64,),dtype=float)
        self.ADCoffset = np.zeros((64,),dtype=float)

class calibrationMapHyb():
    def __init__(self):
        
        self.ID      = -1
        self.Ring    = -1
        self.Fen     = 0
        self.Hybrid  = -1
        self.vmm0    = calibrationMapVMM()
        self.vmm1    = calibrationMapVMM()
        
class calibrationMap():
     def __init__(self,config):
         
         self.config = config
         
         self.numOfHyb = len(config.DETmap.cassettesMap)
         
         self.map = []
         
         for nn in range(self.numOfHyb):
             
             temp = calibrationMapHyb()             
             self.map.append(temp)
             

###############################################################################
###############################################################################

class read_json_calib():
    def __init__(self, calibFile_PathAndFileName, config):
        
        temp =  os.path.split(calibFile_PathAndFileName)
            
        self.calibFilePath = temp[0]+'/'
        self.calibFileName = temp[1]
        
        # self.calibFlag = True
        
        self.config = config
        
        self.calibrationMap = calibrationMap(self.config)
        
        self.numOfCalibs = 0
        
        try:
            self.ff    = open(calibFile_PathAndFileName,'r') 
            self.calib = json.load(self.ff)
            self.calibFlag = True
        except:
            print('\n \033[1;33m---> WARNING: Calibration File ' + self.calibFileName + ' not found \033[1;37m')
            print(' ---> in folder: ' + self.calibFilePath + ' \n ---> calibration skipped.')
            self.calibFlag = False
        else:
            self.print_calib()
            self.calibrations = self.calib.get('Calibrations')
            self.numOfCalibs  = len(self.calibrations)
            self.load()   
            # self.calibFlag = True
        

    def __del__(self):
        try:
            self.ff.close()
        except:
            pass
        
    def print_calib(self):
        print('\033[1;36m\nLoading VMM calibration file: {}\033[1;37m'.format(self.calibFileName))
        
    def loadCalibrationMetadata(self):
        
        self.calib.get('Date')
        
        # etc etc etc 

        
    def get_calibration1hybrid(self,Ring=0,Fen=0,Hybrid=0):
        
        # Fen not used for now!!!! 
        
        foundFlag = False
        
        for nn in range(self.numOfCalibs):
            
            temp = self.calibrations[nn].get('VMMHybridCalibration')
            
            hybridIDtext = temp.get('HybridId')
            
            temp2 = hybridIDtext.split('_',1)
            
            RingCal   = int(temp2[0].split('FEN')[1])
            FenCal    = 0
            HybridCal = int(temp2[1])
            
            # print('index {} -> r {} == {}, f {} == {}, h {} == {}'.format(nn,RingCal,Ring,FenCal,Fen,HybridCal,Hybrid))
            
            if (Ring == RingCal) and (Fen == FenCal) and (Hybrid == HybridCal) :
                
                # print(hybridIDtext)
                
                # print('ciao')
                
                foundFlag = True
                
                vmm0 = temp.get('vmm0')
                vmm1 = temp.get('vmm1')
                
                self.offset_vmm0 = np.array(vmm0.get('adc_offset'))
                self.offset_vmm1 = np.array(vmm1.get('adc_offset'))
                self.slope_vmm0  = np.array(vmm0.get('adc_slope'))
                self.slope_vmm1  = np.array(vmm1.get('adc_slope'))
                
                if foundFlag: 
                    break
                
            else:
 
                foundFlag = False
                
        return foundFlag
    
      
    def load(self):
        
        self.numOfHyb = self.calibrationMap.numOfHyb
        
        for jj in range(self.numOfHyb):
            
            # print(nn)
            
            temp = self.config.DETmap.cassettesMap[jj]
            
            ID       = temp.get('ID')
            RingID   = temp.get('Ring')
            FenID    = temp.get('Fen')
            HybridID = temp.get('Hybrid')
            
            # print('index {} -> r {}, f {}, h {}'.format(jj,RingID,FenID,HybridID))
            
            foundFlag = self.get_calibration1hybrid(RingID,FenID,HybridID)
            
            # print(foundFlag)
            # foundFlag = True
            
            self.calibrationMap.map[jj].ID     = ID
            self.calibrationMap.map[jj].Ring   = RingID
            self.calibrationMap.map[jj].Fen    = FenID
            self.calibrationMap.map[jj].Hybrid = HybridID
            
            if foundFlag is True:
             
                self.calibrationMap.map[jj].vmm0.ADCoffset = self.offset_vmm0
                self.calibrationMap.map[jj].vmm0.ADCslope  = self.slope_vmm0
                self.calibrationMap.map[jj].vmm1.ADCoffset = self.offset_vmm1
                self.calibrationMap.map[jj].vmm1.ADCslope  = self.slope_vmm1
                
            else:
                
                print("No calib found in calib file for Ring {}, Fen {}, Hybrid {} -> using defaults: slope 1, offset 0".format(RingID,FenID,HybridID))
            
            
        
        
###############################################################################
###############################################################################     

class calibrateVMMasic():
    def __init__(self,ADCin,channel,ADCoffset,ADCslope,vmm=0):
    
        rangeSpan = [0,64]
        
        #  if you want to make specific for MB and for speed only calib 32 ch in VMM1
        # if   vmm == 0:
        #     rangeSpan = [0,64]
        # elif vmm == 1:
        #     rangeSpan = [16,48]
            
            
        self.ADCout = np.copy(ADCin)  
        
        # rangeSpan = [0,1]
        
        # loop over all channels 
        for kk in range(rangeSpan[0],rangeSpan[1]):
            
            # print('-->'+str(kk))
            
            sel = channel == kk 
            
            self.ADCout[sel] = pcapr.VMM3A_calibrate_ADC(ADCin[sel],ADCoffset[kk],ADCslope[kk]).ADC_calibrated  
      

class calibrate():
    def __init__(self,readouts,config,calib):
        
        self.readouts = readouts
        self.config   = config
        self.calib    = calib
        
  
    def calibrateTimeStamp(self,NSperClockTick):
        
        # not finished! 
        
        time_offset = 1
        
        time_slope  = 1
        
        TDC_calib = pcapr.VMM3A_convertCalibrate_TDCinSec(self.readouts.TDC,NSperClockTick,time_offset,time_slope).TDC_ns
          
        self.readouts.timeStamp  = self.readouts.timeCoarse + TDC_calib
       
        
    def calibrateADC(self):  
        
        print('\033[1;36mCalibrating ADCs ... \033[1;37m')
        
        self.numOfHyb = len(self.config.DETmap.cassettesMap)
        
        # self.numOfHyb = 1
        
        for nn in range(self.numOfHyb):
            
            temp = self.calib.calibrationMap.map[nn]
            
            # ID     = temp.ID
            Ring   = temp.Ring
            Fen    = temp.Fen
            Hybrid = temp.Hybrid
            
            selRF  = np.logical_and(self.readouts.Ring == Ring , self.readouts.Fen == Fen)
            selHyb = np.logical_and(selRF, self.readouts.hybrid == Hybrid )
            
            if self.config.DETparameters.operationMode == 'normal':
                selVmm0 = np.logical_and( self.readouts.ASIC == 0, selHyb)
                
                # test 
                # if Ring == 0 and Hybrid == 3:
                    
                #     temp.vmm0.ADCoffset = np.zeros((64,))
                #     temp.vmm0.ADCslope  = 10000*np.ones((64,))
    
                cal0 = calibrateVMMasic(self.readouts.ADC[selVmm0],self.readouts.Channel[selVmm0],temp.vmm0.ADCoffset,temp.vmm0.ADCslope,vmm=0)
                
                self.readouts.ADC[selVmm0] = cal0.ADCout
                
            elif self.config.DETparameters.operationMode == 'clustered':
                cal0 = calibrateVMMasic(self.readouts.ADC,self.readouts.Channel,temp.vmm0.ADCoffset,temp.vmm0.ADCslope,vmm=0)
                self.readouts.ADC = cal0.ADCout
            
            # print(self.readouts.ADC[selVmm0])
            # print(cal0.ADCout)
            # # self.deb = np.concatenate(self.readouts.ADC[selVmm0],cal0.ADCout)
            # # print(cal0.ADCout)
            
            if self.config.DETparameters.operationMode == 'normal':
                selVmm1 = np.logical_and( self.readouts.ASIC == 1, selHyb)
                cal1 = calibrateVMMasic(self.readouts.ADC[selVmm1],self.readouts.Channel[selVmm1],temp.vmm1.ADCoffset,temp.vmm1.ADCslope,vmm=1)
                self.readouts.ADC[selVmm1] = cal1.ADCout
            elif self.config.DETparameters.operationMode == 'clustered':
                cal1 = calibrateVMMasic(self.readouts.ADC1,self.readouts.Channel1,temp.vmm1.ADCoffset,temp.vmm1.ADCslope,vmm=1)
                self.readouts.ADC1 = cal1.ADCout
            
      
  
###############################################################################
###############################################################################

if __name__ == '__main__':

   filePathCalib   = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/calib/'+'AMOR_calib_20231111002842.json'
   
   
   filePathCalib   = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/calib/'+'testCalib.json'

   
   filePathConfig  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/config/'+'AMOR.json'
   
   
   filePathData    = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/data/'+'freiatest.pcapng'


   filePathData    = '/Users/francescopiscitelli/Documents/DOC/DATA/202311_PSI_AMOR_MBnewAMOR_VMM_neutrons/SamplesAndMasks/20231106_142811_duration_s_5_YESneutrons1240K1070Rth280_maskESS_00000.pcapng'
   
   config = maps.read_json_config(filePathConfig)
   
   
   # filePath = path+'pcap_for_fra.pcapng'
   # filePath = path+'pcap_for_fra_ch2test.pcapng'
   # filePath = path+'pcap_for_fra_ch2test_take2.pcapng'
   # filePath = path+'pcap_for_fra_coinc.pcapng'
   # filePath = path+'freiatest.pcapng'
   NSperClockTick = 11.356860963629653  #ns per tick ESS for 88.0525 MHz
   # pcap = pcapr.pcapng_reader(filePathData, timeResolution)
   # readouts = pcap.readouts 
   # readoutsArrayIn = readouts.concatenateReadoutsInArrayForDebug()
   
   pcap = pcapr.pcapng_reader(filePathData, NSperClockTick, sortByTimeStampsONOFF = False)
   readouts  = pcap.readouts
   
   readouts_array = readouts.concatenateReadoutsInArrayForDebug()
  

   calib = read_json_calib(filePathCalib,config)
   
   
   # print(calib.calibFlag)
   
   # aa = calib.calibrationMap.map[0]
   
   # ad = calibrateVMMasic(100*np.ones((64,)),34*np.ones((64,)),aa.vmm0.ADCslope,aa.vmm0.ADCoffset,vmm=0)
   
   # adc = ad.ADCout
   
   
   
   
   cal = calibrate(readouts,config,calib)
   cal.calibrateADC()
   
   readoutsOut = cal.readouts
   
   readouts_arrayOut = readoutsOut.concatenateReadoutsInArrayForDebug()
   
   readall = np.concatenate((readouts_array[:,0:6],readouts_arrayOut[:,0:6]),axis = 1)
   
 
