#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 10:46:16 2021

@author: francescopiscitelli
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# import time

# from lib import libSampleData as sdat
# from lib import libMapping as maps
# from lib import libCluster as clu
# from lib import libAbsUnitsAndLambda as absu
from lib import libHistograms as hh
# from lib import libParameters as para
# from lib import libReadPcapngVMM as pcapr


# import libSampleData as sdat
# import libMappingMG as maps
# import libCluster as clu
# import libAbsUnitsAndLambda as absu
# import libHistograms as hh
# import libParameters as para
# import libReadPcapngVMM as pcapr

###############################################################################
###############################################################################

class preparePlotMatrix():
    def __init__(self, figNum, Nrows,  Ncols=1, figSize=(12,12)):
        
        self.Ncols = Ncols
        self.Nrows = Nrows
        
        self.figSize = figSize
        
        self.sharex='col'
        self.sharey='row' 
        
        self.figHandle, self.axHandle = plt.subplots(num=figNum, figsize=self.figSize, nrows=self.Nrows, ncols=self.Ncols, sharex=self.sharex, sharey=self.sharey)
        self.axHandle.shape      = (self.Nrows,self.Ncols)
        self.axHandle            = np.atleast_2d(self.axHandle)
  
################################

########################################

class checkReadoutsClass():
     def __init__(self, readouts):
         
         self.readouts = readouts
         
         self.flag   =  True
         
         if len(self.readouts.Channel) == 0 :
             
             self.flag  = False
             
             print('\t \033[1;33mWARNING: Readouts array is empty -> skipping plots')
                   
                   
        
class plottingReadouts():
   
    def __init__(self, readouts, config):
    
        # self.readouts = readouts
        
        self.config   = config
        
        checkke = checkReadoutsClass(readouts)
        self.readouts = checkke.readouts
        self.flag = checkke.flag 
        
        # self.parameters = parameters
        if self.flag is True:
            self.xbins = np.linspace(0,63,64)
        
    def histChRaw1hybrid(self,cassette1ID):
        
        if self.flag is True:
            self.config.get_cassID2RingFenHybrid(cassette1ID)
    
            # cassette correpsonds to a particluar triplet of ring fen and hybrid so Icanuse cassette to loop over hybrids
            
            # print(self.config.cassMap.RingID)
            # print(self.config.cassMap.FenID)
            # print(self.config.cassMap.hybridID)
            
            sel1 = self.readouts.Ring   == self.config.cassMap.RingID
            sel2 = self.readouts.Fen    == self.config.cassMap.FenID
            sel3 = self.readouts.hybrid == self.config.cassMap.hybridWID
            sel4 = self.readouts.hybrid == self.config.cassMap.hybridSID
            
            # sel = sel1 & sel2 & sel3
            
            selW = sel1 & sel2 & sel3
            selS = sel1 & sel2 & sel4
            
            # wireCh0to31 = np.mod(self.hits.WiresStrips[cass & wires],parameters.config.DETparameters.numOfWires)
            
            if self.config.DETparameters.operationMode == 'normal':
                asic0  = self.readouts.ASIC == 0
                asic1  = self.readouts.ASIC == 1
                self.histo0 = hh.histog().hist1D(self.xbins, self.readouts.Channel[selW & asic0])
                self.histo1 = hh.histog().hist1D(self.xbins, self.readouts.Channel[selW & asic1])
                self.histo2 = hh.histog().hist1D(self.xbins, self.readouts.Channel[selS & asic0])
                self.histo3 = hh.histog().hist1D(self.xbins, self.readouts.Channel[selS & asic1])
                
            # elif self.config.DETparameters.operationMode == 'clustered':
            #     # self.histo0 = hh.histog().hist1D(self.xbins, self.readouts.Channel[sel])
            #     # self.histo1 = hh.histog().hist1D(self.xbins, self.readouts.Channel1[sel])
            #     pass
        
    def plotChRaw(self,hybrids): 
        
        if self.flag is True:
        
            self.ploth = preparePlotMatrix(1001, 4, len(hybrids))
            
            self.ploth.figHandle.suptitle('Readouts - raw channels')
    
            for k, cc in enumerate(hybrids):
                
                self.histChRaw1hybrid(cc)  
    
                self.ploth.axHandle[0][k].bar(self.xbins,self.histo0,0.8,color='r') 
                self.ploth.axHandle[1][k].bar(self.xbins,self.histo1,0.8,color='m')
                self.ploth.axHandle[2][k].bar(self.xbins,self.histo2,0.8,color='b') 
                self.ploth.axHandle[3][k].bar(self.xbins,self.histo3,0.8,color='c')
                
                self.ploth.axHandle[0][k].set_xlabel('WIRE H ASIC 0 ch no.')
                self.ploth.axHandle[1][k].set_xlabel('WIRE H ASIC 1 ch no.')
                self.ploth.axHandle[2][k].set_xlabel('GRID H ASIC 0 ch no.')
                self.ploth.axHandle[3][k].set_xlabel('GRID H ASIC 1 ch no.')
                
                self.ploth.axHandle[0][k].set_title('MG.'+str(cc)) 
            
    def extractTimeStamp1hybrid(self,cassette1ID):
        
        if self.flag is True:
        
            self.config.get_cassID2RingFenHybrid(cassette1ID)
            
            sel1 = self.readouts.Ring   == self.config.cassMap.RingID
            sel2 = self.readouts.Fen    == self.config.cassMap.FenID
            sel3 = self.readouts.hybrid == self.config.cassMap.hybridWID
            sel4 = self.readouts.hybrid == self.config.cassMap.hybridSID
            
            selW = sel1 & sel2 & sel3
            selS = sel1 & sel2 & sel4
            
            if self.config.DETparameters.operationMode == 'normal':
                asic0  = self.readouts.ASIC == 0
                asic1  = self.readouts.ASIC == 1

                self.timeStamp0 = self.readouts.timeStamp[selW & asic0]
                self.timeStamp1 = self.readouts.timeStamp[selW & asic1]
                self.timeStamp2 = self.readouts.timeStamp[selS & asic0]
                self.timeStamp3 = self.readouts.timeStamp[selS & asic1]
                
            # elif self.config.DETparameters.operationMode == 'clustered':
            #     self.timeStamp0 = self.readouts.timeStamp[sel]
            #     self.timeStamp1 = self.timeStamp0
            
            
    def plotTimeStamps(self,hybrids):
        
        if self.flag is True:
        
            self.plotht = preparePlotMatrix(1002, 4, len(hybrids))
            
            self.plotht.figHandle.suptitle('Readouts - raw channels time stamps')
            
            for k, cc in enumerate(hybrids):
                
                self.extractTimeStamp1hybrid(cc)  
                
                xx0 = np.arange(0,len(self.timeStamp0),1)
                xx1 = np.arange(0,len(self.timeStamp1),1)
                xx2 = np.arange(0,len(self.timeStamp2),1)
                xx3 = np.arange(0,len(self.timeStamp3),1)
                
                self.plotht.axHandle[0][k].scatter(xx0,self.timeStamp0,0.8,color='r',marker='+') 
                self.plotht.axHandle[1][k].scatter(xx1,self.timeStamp1,0.8,color='m',marker='+')
                self.plotht.axHandle[2][k].scatter(xx2,self.timeStamp2,0.8,color='b',marker='+') 
                self.plotht.axHandle[3][k].scatter(xx3,self.timeStamp3,0.8,color='c',marker='+')
                self.plotht.axHandle[0][k].set_xlabel('WIRE H ASIC 0 trigger no.')
                self.plotht.axHandle[1][k].set_xlabel('WIRE H ASIC 1 trigger no.')
                self.plotht.axHandle[2][k].set_xlabel('GRID H ASIC 0 trigger no.')
                self.plotht.axHandle[3][k].set_xlabel('GRID H ASIC 1 trigger no.')
                self.plotht.axHandle[0][k].set_ylabel('time (ns)')
                self.plotht.axHandle[1][k].set_ylabel('time (ns)')
                self.plotht.axHandle[2][k].set_ylabel('time (ns)')
                self.plotht.axHandle[3][k].set_ylabel('time (ns)')
                self.plotht.axHandle[0][k].set_title('MG.'+str(cc)) 
                self.plotht.axHandle[0][k].grid(axis='x', alpha=0.75)
                self.plotht.axHandle[1][k].grid(axis='x', alpha=0.75)
                self.plotht.axHandle[2][k].grid(axis='x', alpha=0.75)
                self.plotht.axHandle[3][k].grid(axis='x', alpha=0.75)
                self.plotht.axHandle[0][k].grid(axis='y', alpha=0.75)
                self.plotht.axHandle[1][k].grid(axis='y', alpha=0.75)
                self.plotht.axHandle[2][k].grid(axis='y', alpha=0.75)
                self.plotht.axHandle[3][k].grid(axis='y', alpha=0.75)
            
    def plotChoppResets(self):
            
            if self.flag is True:
            
                self.plothc = preparePlotMatrix(9999, 1, 2)
                
                self.plothc.figHandle.suptitle('Chopper Resets - all channels')
                
                xax = np.arange(0,len(self.readouts.PulseT))
                
                resetsTime = (self.readouts.PulseT - self.readouts.PulseT[0]) / 1e9    # from ns to s 
                
                self.plothc.axHandle[0][0].scatter(xax,resetsTime,0.8,color='k',marker='+') 
                self.plothc.axHandle[0][0].set_xlabel('trigger no.')
                self.plothc.axHandle[0][0].set_ylabel('time (s)')
                
                deltaTime = np.concatenate((np.diff(resetsTime),[0]),axis=0)
                
                self.plothc.axHandle[0][1].scatter(xax,deltaTime,0.8,color='k',marker='+') 
                self.plothc.axHandle[0][1].set_xlabel('trigger no.')
                self.plothc.axHandle[0][1].set_ylabel('delta time betweeen resets (s)')
            
################################
        
class plottingHits():
   
    def __init__(self, hits, parameters):
        
        self.hits = hits
        
        self.parameters = parameters
        
        self.xbins = np.linspace(0,127,128)
        
    def histChRaw1Cass(self,cassette1ID):
        
        cass = self.hits.Cassette == cassette1ID
        
        if self.parameters.config.DETparameters.operationMode == 'normal':
            wires  = self.hits.WorS == 0
            strips = self.hits.WorS == 1
            wireCh0to31 = np.mod(self.hits.WiresStrips[cass & wires],self.parameters.config.DETparameters.numOfWires)
            self.histow = hh.histog().hist1D(self.xbins, wireCh0to31)
            self.histos = hh.histog().hist1D(self.xbins, self.hits.WiresStrips[cass & strips])
            
        # elif self.parameters.config.DETparameters.operationMode == 'clustered':
        #     wireCh0to31 = np.mod(self.hits.WiresStrips1[cass],self.parameters.config.DETparameters.numOfWires)
        #     self.histow = hh.histog().hist1D(self.xbins, wireCh0to31)
        #     self.histos = hh.histog().hist1D(self.xbins, self.hits.WiresStrips[cass])
        

    def plotChRaw(self,cassettes): 
        
        self.ploth = preparePlotMatrix(1003, 2, len(cassettes))
        
        self.ploth.figHandle.suptitle('Hits - raw channels')

        for k, cc in enumerate(cassettes):
            
            self.histChRaw1Cass(cc)  

            self.ploth.axHandle[0][k].bar(self.xbins,self.histow,0.8,color='r') 
            self.ploth.axHandle[1][k].bar(self.xbins,self.histos,0.8,color='b')
            self.ploth.axHandle[0][k].set_xlabel('hit wire ch no.')
            self.ploth.axHandle[1][k].set_xlabel('hit strip ch no.')
            self.ploth.axHandle[0][k].set_title('cass ID '+str(cc))                       
        
            # self.ploth.axHandle[0][k]
            
    def extractTimeStamp1Cass(self,cassette1ID):
             
        sel = self.hits.Cassette == cassette1ID
        
        if self.parameters.config.DETparameters.operationMode == 'normal':
            isWire   = self.hits.WorS == 0
            isStrip  = self.hits.WorS == 1
            
            # wireCh0to31 = np.mod(self.hits.WiresStrips[sel & isWire],self.parameters.config.DETparameters.numOfWires)
            
            self.timeStampW = self.hits.timeStamp[sel] * isWire[sel]
            self.timeStampS = self.hits.timeStamp[sel] * isStrip[sel]
            
            # self.timeStampW[self.timeStampW == 0] = np.ma.masked # same as np.nan for int64 instead of floats   
            # self.timeStampS[self.timeStampS == 0] = np.ma.masked # same as np.nan for int64 instead of floats
 
        # elif self.parameters.config.DETparameters.operationMode == 'clustered':
            
        #     self.timeStampW = self.hits.timeStamp[sel] 
        #     self.timeStampS = self.hits.timeStamp[sel] 
            
        # self.timeStampW = np.ma.masked_where(self.timeStampW == 0,    self.timeStampW) 
        # self.timeStampS = np.ma.masked_where(self.timeStampS == 0,    self.timeStampS) 
            
        
    def extractTimeStampAndCh1Cass(self,cassette1ID):
           
        self.extractTimeStamp1Cass(cassette1ID) 
        
        sel = self.hits.Cassette == cassette1ID
        
        if self.parameters.config.DETparameters.operationMode == 'normal':
            isWire   = self.hits.WorS == 0
            isStrip  = self.hits.WorS == 1
            
            wireCh0to31 = np.mod(self.hits.WiresStrips,self.parameters.config.DETparameters.numOfWires) 
            
            self.WireCh  = np.round((wireCh0to31[sel]+10) * isWire[sel])
            self.StripCh = np.round((self.hits.WiresStrips[sel]+20) * isStrip[sel])
   
        elif self.parameters.config.DETparameters.operationMode == 'clustered':
            wireCh0to31 = np.mod(self.hits.WiresStrips1,self.parameters.config.DETparameters.numOfWires) 
            
            self.WireCh  = np.round((wireCh0to31[sel]+10))
            self.StripCh = np.round((self.hits.WiresStrips[sel]+20))
            
        self.WireCh[self.WireCh == 0]   = np.ma.masked # same as np.nan for int64 instead of floats   
        self.StripCh[self.StripCh == 0] = np.ma.masked # same as np.nan for int64 instead of floats
        
        self.WireCh  = self.WireCh   - 10
        self.StripCh = self.StripCh  - 20 + self.parameters.config.DETparameters.numOfWires
        
    def plotTimeStamps(self,cassettes):
        
        self.plotht = preparePlotMatrix(1004, 1, len(cassettes))
        
        self.plotht.figHandle.suptitle('Hits - W and S time stamps')
        
        for k, cc in enumerate(cassettes):
            
            self.extractTimeStamp1Cass(cc)  
            
            xx0 = np.arange(0,len(self.timeStampW),1)
            xx1 = np.arange(0,len(self.timeStampS),1)
            
            self.plotht.axHandle[0][k].scatter(xx0,self.timeStampW,0.8,color='r',marker='+') 
            self.plotht.axHandle[0][k].scatter(xx1,self.timeStampS,0.8,color='b',marker='+')
            self.plotht.axHandle[0][k].set_xlabel('trigger no.')   
            self.plotht.axHandle[0][k].set_ylabel('time (ns)')
            self.plotht.axHandle[0][k].set_title('cass ID '+str(cc)) 
            self.plotht.axHandle[0][k].grid(axis='x', alpha=0.75)
            self.plotht.axHandle[0][k].grid(axis='y', alpha=0.75)
            
    def plotTimeStampsVSCh(self,cassettes):
        
        self.plothtvs = preparePlotMatrix(1005, 1, len(cassettes))
        
        self.plothtvs.figHandle.suptitle('Hits - W and S VS time stamps')
        
        for k, cc in enumerate(cassettes):
            
            self.extractTimeStampAndCh1Cass(cc)  
            
            self.plothtvs.axHandle[0][k].scatter(self.WireCh,self.timeStampW,0.8,color='r',marker='+') 
            self.plothtvs.axHandle[0][k].scatter(self.StripCh,self.timeStampS,0.8,color='b',marker='+')
            self.plothtvs.axHandle[0][k].set_ylabel('time (ns)')   
            self.plothtvs.axHandle[0][k].set_xlabel('W or S channel (after mapping)')
            self.plothtvs.axHandle[0][k].set_title('cass ID '+str(cc)) 
            self.plothtvs.axHandle[0][k].grid(axis='x', alpha=0.75)
            self.plothtvs.axHandle[0][k].grid(axis='y', alpha=0.75)
        

################################

class logScaleMap():
    def __init__(self,logScale):
        
        if logScale is True:
            self.normColors = LogNorm()
        elif logScale is False:
            self.normColors = None

########################################

class checkEventsClass():
     def __init__(self, events):
         
         self.events = events
         
         self.flag   =  True
         
         if len(self.events.positionW) == 0 :
             
             self.flag  = False
             
             print('\t \033[1;33mWARNING: Events array is empty -> skipping plots')
             
         else:
        
             if  len(self.events.ToF) == 0 :
                   print('\t \033[1;33mWARNING: ToF array is empty ')
                   self.events.ToF = np.zeros((len(self.events.positionW)),dtype='int64')
                     
        

#######################################
        
class plottingEvents():
    
    def __init__(self, events, allAxis, coincidenceWS_ONOFF, parameters):
        
        # self.Ncass = Ncass
        
        # self.events  = events
        self.allAxis = allAxis
        self.parameters = parameters
        
        # self.selectCoinc = events.positionS >= -2
        
        self.coincidenceWS_ONOFF = coincidenceWS_ONOFF
        
        checkke = checkEventsClass(events)
        self.events = checkke.events
        self.flag   = checkke.flag
    
        if self.flag is True:
            if self.coincidenceWS_ONOFF is True:
                print('\t building histograms ... coincidence W/S ON for ToF and Lambda ...')
                self.selc = events.positionS >= 0 
            else:
                print('\t building histograms ... coincidence W/S OFF for ToF and Lambda ...')
                self.selc = events.positionS >= - np.inf
        
        # self.sharex='col'
        # self.sharey='row' 
            
    def plotXYToF(self, logScale = False, absUnits = False, orientation = 'vertical'):
        
        if self.flag is True:
            normColors = logScaleMap(logScale).normColors
         
            if absUnits == False:
                
                # if  len(self.events.ToF) == 0 :
                #    print('\t \033[1;33mWARNING: ToF array is empty ')
                #    self.events.ToF = np.zeros((len(self.events.positionW)),dtype='int64')
                
                h2D, _, hToF = hh.histog().histXYZ(self.allAxis.axWires.axis, self.events.positionW[self.selc], self.allAxis.axStrips.axis, self.events.positionS[self.selc], self.allAxis.axToF.axis, self.events.ToF[self.selc]/1e9)
        
                hProjAll = hh.histog().hist1D(self.allAxis.axWires.axis, self.events.positionW)
                
                hProj2D  = np.sum(h2D,axis=0)
                
                if orientation == 'vertical':
                
                    fig2D, (ax1, ax2) = plt.subplots(num=101,figsize=(6,12), nrows=2, ncols=1)    
                    # #fig.add_axes([0,0,1,1]) #if you want to position absolute coordinate
                    pos1  = ax1.imshow(h2D,aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axWires.start,self.allAxis.axWires.stop,self.allAxis.axStrips.stop,self.allAxis.axStrips.start], origin='upper',cmap='viridis')
                    
                    for k in range(1,self.parameters.config.DETparameters.numOfCassettes):
                        ax1.plot([k*self.parameters.config.DETparameters.numOfWires, k*self.parameters.config.DETparameters.numOfWires], [0, self.parameters.config.DETparameters.numOfStrips-1], color='r', linewidth = 1)
                    
                    #  temporary fix because LogNorm crashes tihe imShow when Log 
                    try:
                        fig2D.colorbar(pos1, ax=ax1, orientation="horizontal",fraction=0.07,anchor=(1.0,0.0))
                    except:
                        print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                        
                    # cbar1 =fig2D.colorbar(pos1,ax=ax1)
                    # cbar2.minorticks_on()
                    # ax1.set_aspect('tight')
                    ax1.set_xlabel('Wire ch.')
                    ax1.set_ylabel('Strip ch.')
                    fig2D.suptitle('DET image')
                    
                    
                # elif orientation == 'horizontal':   
                           
                #     fig2D, (ax1, ax2) = plt.subplots(num=101,figsize=(12,6), nrows=1, ncols=2) 
                    
                    
                
                #     pos1  = ax1.imshow(np.rot90(h2D,1),aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axStrips.start,self.allAxis.axStrips.stop,self.allAxis.axWires.start,self.allAxis.axWires.stop], origin='upper',cmap='viridis')
                    
                #     #  temporary fix because LogNorm crashes tihe imShow when Log 
                #     try:
                #         fig2D.colorbar(pos1, ax=ax1, orientation="horizontal",fraction=0.07,anchor=(1.0,0.0))
                #     except:
                #         print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                        
                #     # # add red lines to plot IMG
                #     # for k in np.arange(0,445,32):
                #     #     ax1.plot([0,63],[k,k],'r',linewidth=1)
                        
                        
                #     ax1.set_ylabel('Wire ch.')
                #     ax1.set_xlabel('Strip ch.')
                #     fig2D.suptitle('DET image')
                    
    
                
                pos2 = ax2.step(self.allAxis.axWires.axis,hProjAll,'r',where='mid',label='1D')
                ax2.step(self.allAxis.axWires.axis,hProj2D,'b',where='mid',label='2D')
                if logScale is True:
                   ax2.set_yscale('log')
                ax2.set_xlabel('Wire ch.')
                ax2.set_ylabel('counts')
                ax2.set_xlim(self.allAxis.axWires.start,self.allAxis.axWires.stop)
                legend = ax2.legend(loc='upper right', shadow=False, fontsize='large')
    
                ########
                # 2D image of detector ToF vs Wires 
                # ToFxgms = ToFxg*1e3 # in ms 
            
                fig2, ax2 = plt.subplots(num=102,figsize=(6,6), nrows=1, ncols=1) 
                pos2  = ax2.imshow(hToF,aspect='auto',norm=normColors,interpolation='nearest',extent=[self.allAxis.axToF.start*1e3,self.allAxis.axToF.stop*1e3,self.allAxis.axWires.start,self.allAxis.axWires.stop], origin='lower',cmap='viridis')
                
                #  temporary fix because LogNorm crashes tihe imShow when Log 
                try:
                    fig2.colorbar(pos2, ax=ax2)
                except:
                    print('\n --> \033[1;33mWARNING: Cannot plot YToF in Log scale, changed to linear\033[1;37m',end='')
                
                ax2.set_ylabel('Wire ch.')
                ax2.set_xlabel('ToF (ms)')
                fig2.suptitle('DET ToF')
            
            elif absUnits == True:
                
                # if  len(self.events.ToF) == 0 :
                #     print('\t \033[1;33mWARNING: ToF arrasy is empty ')
                #     self.events.ToF = np.zeros(len(self.events.positionWmm),dtype='int64')
                
                
                h2D, hProj, hToF = hh.histog().histXYZ(self.allAxis.axWires_mm.axis, self.events.positionWmm[self.selc], self.allAxis.axStrips_mm.axis, self.events.positionSmm[self.selc], self.allAxis.axToF.axis, self.events.ToF[self.selc]/1e9)    
        
                hProjAll = hh.histog().hist1D(self.allAxis.axWires_mm.axis, self.events.positionWmm)
                
                hProj2D  = np.sum(h2D,axis=0)
                
                if orientation == 'vertical':
                
                    fig2D, (ax1, ax2) = plt.subplots(num=101,figsize=(6,12), nrows=2, ncols=1)    
                    # #fig.add_axes([0,0,1,1]) #if you want to position absolute coordinate
                    pos1  = ax1.imshow(h2D,aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axWires_mm.start,self.allAxis.axWires_mm.stop,self.allAxis.axStrips_mm.stop,self.allAxis.axStrips_mm.start], origin='upper',cmap='viridis')
                    
                    #  temporary fix because LogNorm crashes tihe imShow when Log 
                    try:
                        fig2D.colorbar(pos1, ax=ax1, orientation="horizontal",fraction=0.07,anchor=(1.0,0.0))
                    except:
                        print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                        
                    
                    # cbar1 =fig2D.colorbar(pos1,ax=ax1)
                    # cbar2.minorticks_on()
                    # ax1.set_aspect('tight')
                    ax1.set_xlabel('Wire coord. (mm)')
                    ax1.set_ylabel('Strip (mm)')
                    fig2D.suptitle('DET image')
     
                   
                elif orientation == 'horizontal':  
                    
                    fig2D, (ax1, ax2) = plt.subplots(num=101,figsize=(6,12), nrows=2, ncols=1)    
                    # #fig.add_axes([0,0,1,1]) #if you want to position absolute coordinate
                    pos1  = ax1.imshow(np.rot90(h2D,1),aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axStrips_mm.start,self.allAxis.axStrips_mm.stop,self.allAxis.axWires_mm.start,self.allAxis.axWires_mm.stop], origin='upper',cmap='viridis')
                    
                    #  temporary fix because LogNorm crashes tihe imShow when Log 
                    try:
                        fig2D.colorbar(pos1, ax=ax1, orientation="horizontal",fraction=0.07,anchor=(1.0,0.0))
                    except:
                        print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                    
                    
                    ax1.set_ylabel('Wire coord. (mm)')
                    ax1.set_xlabel('Strip (mm)')
                    fig2D.suptitle('DET image')
        
        
                pos2 = ax2.step(self.allAxis.axWires_mm.axis,hProjAll,'r',where='mid',label='1D')
                ax2.step(self.allAxis.axWires_mm.axis,hProj2D,'b',where='mid',label='2D')
                if logScale is True:
                   ax2.set_yscale('log')
                ax2.set_xlabel('Wire coord. (mm)')
                ax2.set_ylabel('counts')
                ax2.set_xlim(self.allAxis.axWires_mm.start,self.allAxis.axWires_mm.stop)
                legend = ax2.legend(loc='upper right', shadow=False, fontsize='large')
            
            
                fig2, ax2 = plt.subplots(num=102,figsize=(6,6), nrows=1, ncols=1) 
                pos2  = ax2.imshow(hToF,aspect='auto',norm=normColors,interpolation='nearest',extent=[self.allAxis.axToF.start*1e3,self.allAxis.axToF.stop*1e3,self.allAxis.axWires_mm.start,self.allAxis.axWires_mm.stop], origin='lower',cmap='viridis')
                
                #  temporary fix because LogNorm crashes tihe imShow when Log 
                try:
                     fig2.colorbar(pos2, ax=ax2)
                except:
                     print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                
                
                ax2.set_ylabel('Wire coord. (mm)')
                ax2.set_xlabel('ToF (ms)')
                fig2.suptitle('DET ToF')
                
            # return h2D
            

    def plotXLambda(self, logScale=False, absUnits = False):
        
        if self.flag is True:
            
            normColors = logScaleMap(logScale).normColors
            
            if absUnits is False:
                h = hh.histog().hist2D(self.allAxis.axLambda.axis, self.events.wavelength[self.selc], self.allAxis.axWires.axis , self.events.positionW[self.selc])
                
                figl, axl = plt.subplots(num=103,figsize=(6,6), nrows=1, ncols=1) 
                posl1  = axl.imshow(h,aspect='auto',norm=normColors,interpolation='nearest',extent=[self.allAxis.axLambda.start,self.allAxis.axLambda.stop,self.allAxis.axWires.start,self.allAxis.axWires.stop], origin='lower',cmap='viridis')
                
                #  temporary fix because LogNorm crashes tihe imShow when Log 
                try:
                     figl.colorbar(posl1, ax=axl)
                except:
                     print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                
                
                axl.set_ylabel('Wire ch.')
                axl.set_xlabel('wavelength (A)')
                figl.suptitle('DET wavelength')
                
            elif absUnits == True:
                
                h = hh.histog().hist2D(self.allAxis.axLambda.axis, self.events.wavelength[self.selc], self.allAxis.axWires_mm.axis , self.events.positionWmm[self.selc])
                
                figl, axl = plt.subplots(num=103,figsize=(6,6), nrows=1, ncols=1) 
                posl1  = axl.imshow(h,aspect='auto',norm=normColors,interpolation='nearest',extent=[self.allAxis.axLambda.start,self.allAxis.axLambda.stop,self.allAxis.axWires_mm.start,self.allAxis.axWires_mm.stop], origin='lower',cmap='viridis')
                
                #  temporary fix because LogNorm crashes tihe imShow when Log 
                try:
                     figl.colorbar(posl1, ax=axl)
                except:
                     print('\n --> \033[1;33mWARNING: Cannot plot XY in Log scale, changed to linear\033[1;37m',end='')
                                
                axl.set_ylabel('Wire coord. (mm)')
                axl.set_xlabel('wavelength (A)')
                figl.suptitle('DET wavelength')

    def plotMultiplicity(self, cassettes):
        
        if self.flag is True:
        
            self.width      = 0.2
            self.extentplot = 7

            ########

            self.plotMult = preparePlotMatrix(401, 2, len(cassettes))
            
            self.plotMult.figHandle.suptitle('Events - multiplicity')
            
            xx =  self.allAxis.axMult.axis

            for k, cass in enumerate(cassettes):
   
                selc  = self.events.Cassette  == cass
                sel2D = self.events.positionS >= 0

                myw  = hh.histog().hist1D(xx,self.events.multW[selc]) # wires all
                mys  = hh.histog().hist1D(xx,self.events.multS[selc]) # strips all
                mywc = hh.histog().hist1D(xx,self.events.multW[selc & sel2D]) # wires coinc
           
                my2Dwc = hh.histog().hist2D(xx,self.events.multW[selc & sel2D],xx,self.events.multS[selc & sel2D]) # wires coinc with strips 2D
                
                if np.any(selc):
                    mywnorm    = myw/np.sum(myw[1:])
                    mysnorm    = mys/np.sum(mys[1:])
                    mysnormall = mys/np.sum(mys)
                    mywcnorm   = mywc/np.sum(mywc[1:])
                    my2Dwcnorm = my2Dwc/np.sum(my2Dwc)
                else:
                    mywnorm    = np.zeros((len(xx)))
                    mysnorm    = np.zeros((len(xx)))
                    mysnormall = np.zeros((len(xx)))
                    mywcnorm   = np.zeros((len(xx)))
                    my2Dwcnorm = np.zeros((len(xx),len(xx)))

             ########
       
                self.plotMult.axHandle[0][k].bar(xx[:self.extentplot]-self.width,mywnorm[:self.extentplot],self.width,color='m',label='w') 
                self.plotMult.axHandle[0][k].bar(xx[1:self.extentplot]+self.width,mysnorm[1:self.extentplot],self.width,color='b',label='s')
                self.plotMult.axHandle[0][k].bar(xx[0]+self.width,mysnormall[0],self.width,color='c',label='no s')
                self.plotMult.axHandle[0][k].bar(xx[:self.extentplot],mywcnorm[:self.extentplot],self.width,color='r',label='w/s')
                self.plotMult.axHandle[0][k].set_xlabel('multiplicity')
                self.plotMult.axHandle[0][k].set_title('cass ID '+str(cass))
                legend = self.plotMult.axHandle[0][k].legend(loc='upper right', shadow=False, fontsize='large')
                if k == 0:
                    self.plotMult.axHandle[0][k].set_ylabel('probability')
                
                
                pos1 = self.plotMult.axHandle[1][k].imshow(my2Dwcnorm[:self.extentplot,:self.extentplot],aspect='auto',norm=None,interpolation='none',extent=[xx[0]-0.5,xx[self.extentplot]-0.5,xx[0]-0.5,xx[self.extentplot]-0.5], origin='lower',cmap='jet')
                self.plotMult.axHandle[1][k].set_xlabel('multiplicity wires')
                if k == 0:
                    self.plotMult.axHandle[1][k].set_ylabel('multiplicity strips')
                    
                plt.colorbar(pos1,ax=self.plotMult.axHandle[1][k])
       

    def plotPHS(self, cassettes, parameters, logScale = False):
        
        if self.flag is True:
        
            normColors = logScaleMap(logScale).normColors
            
            self.plotPHS = preparePlotMatrix(601, 4, len(cassettes))
            
            self.plotPHS.figHandle.suptitle('Pulse Heigth Spectra')
            
            wireCh0to31Round = np.round(np.mod(self.events.positionW,parameters.config.DETparameters.numOfWires))
                    
            stripChRound     = np.round(self.events.positionS)
            
            wireAx  = np.linspace(0,parameters.config.DETparameters.numOfWires-1, parameters.config.DETparameters.numOfWires)
            
            stripAx = np.linspace(0,parameters.config.DETparameters.numOfStrips-1, parameters.config.DETparameters.numOfStrips)
            
            for k, cass in enumerate(cassettes):
       
                    selc  = self.events.Cassette  == cass
                    sel2D = self.events.positionS >= 0
                    
                    PHSw  = hh.histog().hist2D(self.allAxis.axEnergy.axis,self.events.PHW[selc],wireAx,wireCh0to31Round[selc]) # wires 
                    PHSs  = hh.histog().hist2D(self.allAxis.axEnergy.axis,self.events.PHS[selc & sel2D],stripAx,stripChRound[selc & sel2D]) # strips
                    PHSwc = hh.histog().hist2D(self.allAxis.axEnergy.axis,self.events.PHW[selc & sel2D],wireAx,wireCh0to31Round[selc & sel2D]) # wires coinc with strips 2D
                    
                    self.plotPHS.axHandle[0][k].imshow(PHSw,aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axEnergy.start,self.allAxis.axEnergy.stop,wireAx[0],wireAx[-1]], origin='lower',cmap='jet')
                    self.plotPHS.axHandle[1][k].imshow(PHSs,aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axEnergy.start,self.allAxis.axEnergy.stop,stripAx[0],stripAx[-1]], origin='lower',cmap='jet')
                    self.plotPHS.axHandle[2][k].imshow(PHSwc,aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axEnergy.start,self.allAxis.axEnergy.stop,wireAx[0],wireAx[-1]], origin='lower',cmap='jet')
                    
                    self.plotPHS.axHandle[0][k].set_title('cass ID '+str(cass))
                    if k == 0:
                        self.plotPHS.axHandle[0][k].set_ylabel('wires ch. no.')
                        self.plotPHS.axHandle[1][k].set_ylabel('strips ch. no.')
                        self.plotPHS.axHandle[2][k].set_ylabel('wires coinc. ch. no.')
               
                       #global PHS
                    PHSGw  = np.sum(PHSw,axis=0)
                    PHSGs  = np.sum(PHSs,axis=0)
                    PHSGwc = np.sum(PHSwc,axis=0)
           
                    # global PHS plot
                    self.plotPHS.axHandle[3][k].step(self.allAxis.axEnergy.axis,PHSGw,'r',where='mid',label='w')
                    self.plotPHS.axHandle[3][k].step(self.allAxis.axEnergy.axis,PHSGs,'b',where='mid',label='s')
                    self.plotPHS.axHandle[3][k].step(self.allAxis.axEnergy.axis,PHSGwc,'k',where='mid',label='w/s')
                    self.plotPHS.axHandle[3][k].set_xlabel('pulse height (a.u.)')
                    self.plotPHS.axHandle[3][k].legend(loc='upper right', shadow=False, fontsize='large')
                    if k == 0:
                       self.plotPHS.axHandle[3][k].set_ylabel('counts')
                   
    def plotPHScorrelation(self, cassettes, logScale = False):
        
        if self.flag is True:
              
           normColors = logScaleMap(logScale).normColors
        
           self.plotPHScorr = preparePlotMatrix(602, 1, len(cassettes), figSize=(12,6))
           
           self.plotPHScorr.figHandle.suptitle('Pulse Heigth Spectrum - Correlation W/S')
        
           for k, cass in enumerate(cassettes):
   
                selc  = self.events.Cassette  == cass
                sel2D = self.events.positionS >= 0
                
                PHScorr  = hh.histog().hist2D(self.allAxis.axEnergy.axis,self.events.PHW[selc & sel2D],self.allAxis.axEnergy.axis,self.events.PHS[selc & sel2D]) 
               
                self.plotPHScorr.axHandle[0][k].imshow(PHScorr,aspect='auto',norm=normColors,interpolation='none',extent=[self.allAxis.axEnergy.start,self.allAxis.axEnergy.stop,self.allAxis.axEnergy.start,self.allAxis.axEnergy.stop], origin='lower',cmap='jet')
                
                self.plotPHScorr.axHandle[0][k].set_title('cass ID '+str(cass))
                self.plotPHScorr.axHandle[0][k].set_xlabel('pulse height wires (a.u.)')
                if k == 0:
                    self.plotPHScorr.axHandle[0][k].set_ylabel('pulse height strips (a.u.)')
            
    def plotInstantaneousRate(self, cassettes):
        
        if self.flag is True:
           
           self.plotInst = preparePlotMatrix(209, 1, len(cassettes))
           
           self.plotInst.figHandle.suptitle('Instantaneous Rate')
           
           for k, cass in enumerate(cassettes):
               
               selc  = self.events.Cassette  == cass
               sel2D = self.events.positionS >= 0
               diffeTime = np.diff(self.events.timeStamp[selc & sel2D])
               
               histRate = hh.histog().hist1D(self.allAxis.axInstRate.axis,diffeTime) 
               
               self.plotInst.axHandle[0][k].step(self.allAxis.axInstRate.axis*1e6,histRate,'k',where='mid',label='w')
               self.plotInst.axHandle[0][k].set_xlabel('delta time between events (us)')
               self.plotInst.axHandle[0][k].set_title('cass ID '+str(cass))
               if k == 0:
                   self.plotInst.axHandle[0][k].set_ylabel('num of events')
                   
                   
    def plotToF(self, cassettes):
        
        if self.flag is True:
           
          self.plotTT = preparePlotMatrix(333, 1, len(cassettes))
          
          self.plotTT.figHandle.suptitle('ToF distr per cassette')
          
          for k, cass in enumerate(cassettes):
               
               selc  = self.events.Cassette  == cass
               sel2D = self.events.positionS >= 0
               
               histTT  = hh.histog().hist1D(self.allAxis.axToF.axis,self.events.ToF[selc & sel2D]/1e9) 
               
               histTT1 = hh.histog().hist1D(self.allAxis.axToF.axis,self.events.ToF[selc]/1e9)
               
               self.plotTT.axHandle[0][k].step(self.allAxis.axToF.axis*1e3,histTT1,'r',where='mid',label='all')
               self.plotTT.axHandle[0][k].step(self.allAxis.axToF.axis*1e3,histTT,'b',where='mid',label='2D')
               self.plotTT.axHandle[0][k].set_xlabel('ToF (ms)')
               self.plotTT.axHandle[0][k].set_title('cass ID '+str(cass))
               if k == 0:
                   self.plotTT.axHandle[0][k].set_ylabel('counts')
                   
              
               legend = self.plotTT.axHandle[0][k].legend(loc='upper right', shadow=False, fontsize='large')
               
    def plotLambda(self, cassettes):
        
        if self.flag is True:
           
          self.plotWA = preparePlotMatrix(339, 1, len(cassettes))
          
          self.plotWA.figHandle.suptitle('Wavelength distr per cassette')
          
          for k, cass in enumerate(cassettes):
               
               selc  = self.events.Cassette  == cass
               sel2D = self.events.positionS >= 0
               
               histWA  = hh.histog().hist1D(self.allAxis.axLambda.axis,self.events.wavelength[selc & sel2D]) 
               
               histWA1 = hh.histog().hist1D(self.allAxis.axLambda.axis,self.events.wavelength[selc])
               
               self.plotWA.axHandle[0][k].step(self.allAxis.axLambda.axis,histWA1,'r',where='mid',label='all')
               self.plotWA.axHandle[0][k].step(self.allAxis.axLambda.axis,histWA,'b',where='mid',label='2D')
               self.plotWA.axHandle[0][k].set_xlabel('wavelength (A)')
               self.plotWA.axHandle[0][k].set_title('cass ID '+str(cass))
               if k == 0:
                   self.plotWA.axHandle[0][k].set_ylabel('counts')
                   
              
               legend = self.plotWA.axHandle[0][k].legend(loc='upper right', shadow=False, fontsize='large')
  

  
###############################################################################
  
class plottingMON():
    
    def __init__(self, eventsMON, allAxis):
        
        self.eventsMON = eventsMON
        self.allAxis   = allAxis
        
    def plot_ToF_PHS_MON(self):
        
        histTM  = hh.histog().hist1D(self.allAxis.axToF.axis,self.eventsMON.ToF/1e9) 
        
        histPM  = hh.histog().hist1D(self.allAxis.axEnergy.axis,self.eventsMON.PHW) 

        figMONTOF, (ax1, ax2) = plt.subplots(num=999,figsize=(9,6), nrows=1, ncols=2) 
        figMONTOF.suptitle('MONITOR')
        
        ax1.step(self.allAxis.axToF.axis*1e3,histTM,'k',where='mid',label='MON')
        ax1.set_xlabel('ToF (ms)')
        ax1.set_ylabel('counts')
        ax1.set_title('ToF')
        
        ax2.step(self.allAxis.axEnergy.axis,histPM,'k',where='mid',label='MON')
        ax2.set_xlabel('Pulse Heigth (a.u.)')
        ax2.set_ylabel('counts')
        ax2.set_title('PHS')
    
    
         
    def plotLambda_MON(self):
            
         print('not implemented yet!')
            
###############################################################################
###############################################################################

if __name__ == '__main__' :
    
    plt.close("all")
    
    parameters  = para.parameters('./')
    
    filePath  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap_MGdevel/config/'+"MG.json"

    config = maps.read_json_config(filePath)
    
    parameters.loadConfigAndSetParameters(config)
    
    filePath = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap_MGdevel/data/'
    file = 'MG_2col_2clusters.pcapng'
    
    
    filePathAndFileName = filePath+file
    
    NSperClockTick = 11.356860963629653  #ns per tick ESS for 88.0525 MHz
    
    pcapng = pcapr.pcapng_reader(filePathAndFileName, NSperClockTick, config.MONmap.TTLtype, config.MONmap.RingID,  timeResolutionType='fine', sortByTimeStampsONOFF = False, operationMode=config.DETparameters.operationMode)

    readouts = pcapng.readouts
    readoutsArray = readouts.concatenateReadoutsInArrayForDebug()
 
    plread = plottingReadouts(readouts, config)
    plread.plotChRaw(parameters.config.DETparameters.cassInConfig)
    
    # allAxis = hh.allAxis()
    # allAxis.createAllAxis(parameters)
    

        
    ##################


    
    
    
    