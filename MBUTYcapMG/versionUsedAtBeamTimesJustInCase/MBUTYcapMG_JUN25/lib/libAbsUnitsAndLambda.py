#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 14:18:41 2021

@author: francescopiscitelli
"""

import numpy as np
import time
import matplotlib.pyplot as plt


from lib import libSampleData as sdat
from lib import libMapping as maps
from lib import libCluster as clu
from lib import libParameters as para
from lib import libHistograms as hh
from lib import libPlotting as plo
# from lib import libReadPcapngVMM as pcapr
# from lib import libFileManagmentUtil as fd
# from lib import libParameters as para
# from lib import libTerminal as ta



# import libSampleData as sdat
# import libMapping as maps
# import libCluster as clu
# import libParameters as para
# import libHistograms as hh
# import libPlotting as plo
# import libReadPcapngVMM as pcapr
# import libFileManagmentUtil as fd
# import libParameters as para
# import libTerminal as ta


###############################################################################
###############################################################################

class Tof2LambdaConverter():
    def __init__(self):
        
        ### constants ###
        self.ht        = 1.054e-34  #J*s
        self.mneutr    = 1.67e-27   #Kg
        self.constant  = ( ( (self.ht**2)*4*(3.14159**2) )/(2*self.mneutr) )*1e20*(6.24e18)   #A^2 * eV
        ##############################
    
    def ToF2lambda (self,  distance, ToF):
        
        velocity = distance/ToF               #m/s
        energy   = (1/2)*(self.mneutr/1.6e-19)*(velocity**2)     #eV
        
        selGoods =  energy > 0
        
        self.lamb = np.zeros(len(energy))
         
        self.lamb[selGoods]     = np.sqrt(self.constant/energy[selGoods])                  #A
        self.lamb[~selGoods]    = np.nan                  #A
        # lamb   = np.round(lamb,decimals=2)
        
        return self.lamb
    
    def lambda2ToF (self, distance, lamb):
        
        energy   = self.constant/(lamb**2)
        velocity = np.sqrt(2*energy/((self.mneutr/1.6e-19)))
        self.ToF      = distance/velocity
        
        return self.ToF

###############################################################################

class calculateAbsUnits():
     def __init__(self, events, parameters, text=''): 
         
         if text == '':
             print('\033[1;36m\nCalculating absolute units ... \033[1;37m',end='')
         else: 
             print('\033[1;36m\nCalculating absolute units '+text+'  ... \033[1;37m',end='')
         
         self.events = events
         self.events.createAbsUnitsArrays()
         
         self.parameters = parameters
                 
     def calculatePositionAbsUnit(self):
         
         sine = np.sin(np.deg2rad(self.parameters.config.DETparameters.bladesInclination)) 
         cosi = np.cos(np.deg2rad(self.parameters.config.DETparameters.bladesInclination)) 
         
         #  wires are in global coord so mod 32 to bring it back from 0 to 31
         wireCh0to31 = np.mod(self.events.positionW,self.parameters.config.DETparameters.numOfWires)
         
         for k, cass in enumerate(self.parameters.config.DETparameters.cassInConfig):
             
             selectW = self.events.Cassette == cass
             self.events.positionWmm[selectW]  = np.round((wireCh0to31[selectW] * (self.parameters.config.DETparameters.wirePitch*sine) + self.parameters.config.DETparameters.offset1stWires*k), decimals=2 )  #mm
         
            
         selectS = self.events.positionS >= 0
         self.events.positionSmm[selectS]  = np.round((self.events.positionS[selectS] * self.parameters.config.DETparameters.stripPitch ), decimals = 2) #mm
         self.events.positionSmm[~selectS] = -1                             
           
         #mm Z is 0 at first wire here 
         self.events.positionZmm = np.round((wireCh0to31 * (self.parameters.config.DETparameters.wirePitch*cosi)), decimals = 2) #mm 
             
     def calculateToF(self, removeInvalidToFs = False):
         
          allToFsCounter = np.shape(self.events.ToF)[0]

          self.events.ToF = self.events.timeStamp - self.events.PulseT
          
          invalidToFs = self.events.ToF < 0
          
          invalidToFsCounter1 = np.sum(invalidToFs)
                    
          if invalidToFsCounter1 > 0:
              
              print('\n \033[1;33m\t WARNING ---> %d ToFs (out of %d) are invalid, but corrected with Prev. Pulse Time \033[1;37m' % (invalidToFsCounter1,allToFsCounter))
              
              self.events.ToF[invalidToFs] = self.events.timeStamp[invalidToFs] - self.events.PrevPT[invalidToFs]
              
              
              # self.events.ToF[invalidToFs] = 1000
          
              invalidToFsAgain = self.events.ToF < 0
              
              invalidToFsCounter2 = np.sum(invalidToFsAgain)
          
              if invalidToFsCounter2 > 0:  
              
                  print('\n \033[1;33m\t WARNING ---> %d corrected -> %d ToFs (out of %d) are still invalid after correction (with both PulseT and PrevPT), set as Nan! \033[1;37m' % (invalidToFsCounter1-invalidToFsCounter2,invalidToFsCounter2,allToFsCounter))
                  self.events.ToF[invalidToFsAgain] = np.ma.masked # same as np.nan for int64 instead of floats
                  
                  # self.events.ToF[invalidToFsAgain] = 1000
          
          # self.events.ToF = np.mod(self.events.timeStamp - T0, self.parameters.plotting.ToFduration)
          
          # somehting wiht np.mod the tof range and a shift would work I suspect 
          
          if removeInvalidToFs is True:
              
              self.cleanInvalidToFs()
              
          else:
              
              print('\n \033[1;33m\t Invalid ToFs not removed, kept in the plots! \033[1;37m')

     def calculateWavelength(self):
         
         # self.calculateToF(T0)
         
         # Distance is from chopper to first wire in mm
         ZfromChopper = self.events.positionZmm + self.parameters.wavelength.distance
         
         #  from ns to s, from int ot float ! 
         ToF_s = self.events.ToF/1e9
         
         tcl = Tof2LambdaConverter()

         if self.parameters.wavelength.multipleFramePerReset is True:
               #ToF shifted and corrected by number of bunches
               ToFstart = tcl.lambda2ToF(ZfromChopper*1e-3,self.parameters.wavelength.lambdaMIN)
               tempToF  = ( (ToF_s-ToFstart-self.parameters.wavelength.chopperPickUpDelay) % (self.parameters.wavelength.chopperPeriod/self.parameters.wavelength.numOfBunchesPerPulse) ) + ToFstart
         else:
           
               tempToF  = ToF_s
          
         wavel  = tcl.ToF2lambda(ZfromChopper*1e-3, tempToF) #input m and s, output in A
            
         # append to POPH col 7 of POPH is depth in detector - z (mm) and col 8 is lambda
         
         # tempwavelength = np.round(wavel, decimals=2)
         
         # zeroWave = tempwavelength < 0 
         
         # infWave  = np.isinf(tempwavelength)
         
         # tempwavelength[np.logical_and(zeroWave,infWave)] = np.nan
         
         self.events.wavelength = np.round(wavel, decimals=2)
           
     def calculateToFandWavelength(self,removeInvalidToFs = False):
            
            self.calculateToF(removeInvalidToFs)
            self.calculateWavelength()
            
            
     def cleanInvalidToFs(self):
            
            invalidTofs = self.events.ToF <= 0
            
            condi = np.shape(self.events.positionS)[0] > 0
            
            if condi is True:
            
                sel2D = self.events.positionS >= 0 
                sel1D = self.events.positionS == -1 
                
                invalid2D = np.logical_and(invalidTofs,sel2D)
                invalid1D = np.logical_and(invalidTofs,sel1D)
            
                NumInvalid2D = np.sum(invalid2D)
                NumInvalid1D = np.sum(invalid1D)
            
            NumInvalid = np.sum(invalidTofs)
            
            if NumInvalid > 0:
                
                # print(np.shape(invalidTofs))
                
                if np.shape(self.events.Cassette)[0] != 0:
                    
                    self.events.Cassette     = self.events.Cassette[~invalidTofs] 
                    self.events.multS        = self.events.multS[~invalidTofs] 
                    self.events.multW        = self.events.multW[~invalidTofs] 
                    self.events.positionS    = self.events.positionS[~invalidTofs] 
                    self.events.positionSmm  = self.events.positionSmm[~invalidTofs] 
                    self.events.positionW    = self.events.positionW[~invalidTofs]  
                    self.events.positionWmm  = self.events.positionWmm[~invalidTofs] 
                    self.events.positionZmm  = self.events.positionZmm[~invalidTofs] 
                    
                if np.shape(self.events.wavelength)[0] != 0:
                    self.events.wavelength   = self.events.wavelength[~invalidTofs] 
                
                
                self.events.PHS          = self.events.PHS[~invalidTofs] 
                self.events.PHW          = self.events.PHW[~invalidTofs] 
                self.events.PrevPT       = self.events.PrevPT[~invalidTofs] 
                self.events.PulseT       = self.events.PulseT[~invalidTofs] 
                self.events.ToF          = self.events.ToF[~invalidTofs] 
                self.events.timeStamp    = self.events.timeStamp[~invalidTofs] 
 
                   
                if condi is True:
                    print('\n \033[1;33m\t %d events removed because of Invalid Tofs (%d 2D and %d 1D) \033[1;37m' % (NumInvalid,NumInvalid2D,NumInvalid1D))
                else:
                    print('\n \033[1;33m\t %d events removed because of Invalid Tofs \033[1;37m' % (NumInvalid))
            
                
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
   
   parameters.cassettes.cassettes = [1]
   
       ## distance in mm from chopper and wires 0 of detector
   parameters.wavelength.distance = 8038
    
    ##ON/OFF
   parameters.wavelength.calculateLambda = True
    
    ##ON/OFF
   parameters.wavelength.plotLambda = True
    
   parameters.wavelength.lambdaBins  = 127
   parameters.wavelength.lambdaRange = [1, 25]   #A
    
   parameters.wavelength.chopperPeriod = 0.06 #s
    
    #if chopper has two openings or more per reset of ToF
   parameters.wavelength.multipleFramePerReset = True  #ON/OFF (this only affects the lambda calculation)
   parameters.wavelength.numOfBunchesPerPulse  = 2
   parameters.wavelength.lambdaMIN             = 2.7     #A
    
    # #s on chopper, time shift betweeen pickup and chopper edge 
   parameters.wavelength.chopperPickUpDelay =  13.5/(2.*180.) * parameters.wavelength.chopperPeriod/parameters.wavelength.numOfBunchesPerPulse 
    
    ###############
   parameters.plotting.ToFrange        = 0.1   # s
   parameters.plotting.ToFbinning      = 100e-6 # s
   
   # name = config.get_DetectorName()
   # config.get_cassettesInConfig()
   
   # config.get_cassID2RingFenHybrid(5)
   # print(config.RingID, config.FenID, config.hybridID)
   # config.get_cassID2RingFenHybrid(55)
   # print(config.RingID, config.FenID, config.hybridID)
   # config.get_cassID2RingFenHybrid(1)
   # print(config.RingID, config.FenID, config.hybridID)
   
   # config.get_cassID2RingFenHybrid([1,41])
   

   
   # pr = pcapr.pcap_reader(filePathD)
   #   #  # pr.debug = True
   # pr.read()
   # vmm1 = pr.readouts
   
   # vmm2 = sdat.sampleData()
   # vmm2.fill()

   # cassettes = [1,2,3,4,5]

   # re  = maps.mapDetector(vmm1, filePath)
   #    # re.debug = False
   #    # re.initCatData()

   # re.mappAllCassAndChannels_GlobalCoord(cassettes)
   # hits = re.hits
   
   # Nhits = 30
   # cassettes = [1,2,4,7]
     
   # hits = sdat.sampleHitsMultipleCassettes(cassettes)
   # hits.generate(Nhits)
   
   
    
   # cc = clusterHits(hits) 
   
   # allEvents= events()
   # aa, bb = cc.clusterize1cassette(1,2e-6)
   
   # allEvents.append(cc.events)
   
   # dd = np.loadtxt('dataset1_large_Sorting=True_Filtering=False_Clustered.txt',dtype='float64',delimiter=' ')
   
   
   # aa1, bb1 = cc.clusterize1cassette(3,2e-6)
   
   # allEvents.append(cc.events)
   
   cassette = parameters.cassettes.cassettes
   
   # cc = clu.clusterHits(hits)
   # cc.clusterizeManyCassettes(cassette, 2e-6)
   
   # events =  cc.events
   
    ####################    
# for debug
# generate sample events 
   Nevents = 1e4
   dd = sdat.sampleEventsMultipleCassettes(cassette,'/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/data/')
   dd.generateGlob(Nevents)
   events  = dd.events
   eventsArray2 = events.concatenateEventsInArrayForDebug()
# eventsArray = eventsArray[72:100,:]
#################### 

   events2 = clu.events()
   
   events2.positionW = np.array([0,10])
   events2.positionS = np.array([5,15])
   events2.timeStamp   = np.array([0.01,0.045])
    

   ab = calculateAbsUnits(events, parameters)
   ab.calculatePositionAbsUnit()

   ab.calculateToF()
   
   ab.calculateWavelength()
   
   

   events = ab.events 
   
   eventsArray = events.concatenateEventsInArrayForDebug()
   
   allAxis = hh.allAxis()
   allAxis.createAllAxis(parameters)
   
   #  XY and XToF
   plev = plo.plottingEvents(events,allAxis,True)
   plev.plotXYToF(logScale = parameters.plotting.plotIMGlog, absUnits = parameters.plotting.plotABSunits)
    
    #  lambda
  
   plev.plotXLambda(logScale = parameters.plotting.plotIMGlog, absUnits = parameters.plotting.plotABSunits)
   
   
   
   tElapsedProfiling = time.time() - tProfilingStart
   print('\n Completed in %.2f s' % tElapsedProfiling) 