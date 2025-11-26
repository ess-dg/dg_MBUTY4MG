#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 16:51:12 2021

@author: francescopiscitelli
"""

import numpy as np



# from lib import libSampleData as sdat
# from lib import libMapping as maps
# from lib import libCluster as clu
# from lib import libParameters as para
# from lib import libHistograms as hh
# from lib import libAbsUnitsAndLambda as absu

# from lib import libReadPcapngVMM as pcapr
# from lib import libFileManagmentUtil as fd
# from lib import libTerminal as ta
# from lib import libPlotting as plo



# import libSampleData as sdat
# import libMapping as maps
# import libCluster as clu
# import libParameters as para
# import libHistograms as hh
# import libAbsUnitsAndLambda as absu
# import libReadPcapngVMM as pcapr
# import libFileManagmentUtil as fd
# import libTerminal as ta
# import libPlotting as plo

 
# makes 1D or 2D histograms

###############################################################################
############################################################################### 

class createAx():
    def __init__(self,start,stop,steps):
        
        self.start    = start
        self.stop     = stop
        self.steps    = steps
        
        self.axis =  np.linspace(self.start,self.stop,self.steps)
        
    def update(self):
        
        self.axis =  np.linspace(self.start,self.stop,self.steps)
        
        
class allAxis():        
    def __init__(self):
        
        start = 0
        stop  = 1
        steps = 1

        self.axEnergy    = createAx(start, stop, steps)
        self.axLambda    = createAx(start, stop, steps)
        self.axWires     = createAx(start, stop, steps)
        self.axStrips    = createAx(start, stop, steps)
        self.axWires_mm  = createAx(start, stop, steps)
        self.axStrips_mm = createAx(start, stop, steps)
        self.axToF       = createAx(start, stop, steps)
        self.axMult      = createAx(start, stop, steps)
        self.axInstRate  = createAx(start, stop, steps)
        
    def createAllAxis(self,parameters): 
        
        # param.update()
        
        # sine = np.sin(np.deg2rad(parameters.config.DETparameters.bladesInclination))
        
        sinne = 1 
        
        self.axEnergy = createAx(0, parameters.pulseHeigthSpect.maxEnerg, parameters.pulseHeigthSpect.energyBins)
        self.axToF    = createAx(0, parameters.plotting.ToFrange, parameters.plotting.ToFbins)
        self.axLambda = createAx(parameters.wavelength.lambdaRange[0], parameters.wavelength.lambdaRange[1], parameters.wavelength.lambdaBins)
        
        start  = 0
        stop   = len(parameters.config.DETparameters.cassInConfig)*parameters.config.DETparameters.numOfWires-1
        steps  = len(parameters.config.DETparameters.cassInConfig)*parameters.plotting.posWbins - int(parameters.plotting.posWbins/parameters.config.DETparameters.numOfWires - 1)

        self.axWires  = createAx(start, stop, steps)
        
        start  = 0
        stop   = parameters.config.DETparameters.numOfStrips-1
        steps  = parameters.plotting.posSbins - int(parameters.plotting.posSbins/parameters.config.DETparameters.numOfStrips - 1)

        self.axStrips = createAx(start, stop, steps)
        
        start  = 0
        
        sine = 1

        stop   = (len(parameters.config.DETparameters.cassInConfig)*parameters.config.DETparameters.numOfWires-1)*4*sine
        steps  = len(parameters.config.DETparameters.cassInConfig)*parameters.plotting.posWbins - int(parameters.plotting.posWbins/parameters.config.DETparameters.numOfWires - 1)

        self.axWires_mm  = createAx(start, stop, steps)
        
        start  = 0
        stop   = (parameters.config.DETparameters.numOfStrips-1)*4
        steps  = parameters.plotting.posSbins - int(parameters.plotting.posSbins/parameters.config.DETparameters.numOfStrips - 1)
        
        self.axStrips_mm = createAx(start, stop, steps)
        
        self.axMult = createAx(0, parameters.config.DETparameters.numOfStrips-1, parameters.config.DETparameters.numOfStrips)
        
        start = -parameters.plotting.ToFrange
        stop  = parameters.plotting.ToFrange
        steps = round((stop-start)/parameters.plotting.instRateBin)
        self.axInstRate = createAx(start, stop, steps)
        
    def updateAllAxis(self):
        
        self.axEnergy.update()
        self.axLambda.update()
        self.axWires.update()
        self.axStrips.update()
        self.axWires_mm.update()
        self.axStrips_mm.update()
        self.axToF.update()
        self.axMult.update()
        self.axInstRate.update()
        
###############################################################################
###############################################################################

class histog():
    def __init__(self,outBounds=True):
        
        # print('\t histogram outBounds param set as True as default (Events out of bounds stored in first and last bin)')
        self.outBounds=outBounds
    
    def hist1D(self,xbins,xvar):
    
        binX   = len(xbins) 
            
        Xmin   = np.min(xbins) 
        Xmax   = np.max(xbins) 
    
        self.hist   = np.zeros(binX) 
        
        index = np.int_(np.around(((binX-1)*((xvar-Xmin)/(Xmax-Xmin)))))
        
        if self.outBounds == False:
            if not(np.all(index >= 0) and np.all(index <= binX-1)):
                print('\033[1;33mWARNING: hist out of bounds, change limits!\033[1;37m') 
                return self.hist
        
        for k in range(binX):    
            self.hist[k] = np.sum(index == k) 
           
            if self.outBounds == True:
                # fill overflow last bin and first bin
                self.hist[0]  += np.sum(index<0)
                self.hist[-1] += np.sum(index>binX-1)
            
        return self.hist
    
    ###############################################################################
    ############################################################################### 
        
    def hist2D(self, xbins, xvar, ybins, yvar):
        
        # if np.size(self.ybins) == 1 and np.size(self.yvar) == 1:
        #     hist = 0
        #     return hist
        
        binX   = len(xbins) 
        binY   = len(ybins) 
            
        Xmin   = np.min(xbins) 
        Xmax   = np.max(xbins) 
        
        Ymin   = np.min(ybins) 
        Ymax   = np.max(ybins) 
    
        cont = 0
        
        self.hist = np.zeros((binY,binX)) 
        
        if not( (len(xvar) == len(yvar))):
            print('\n \t \033[1;31m----> ABORTED: X and Y not same length! \033[1;37m\n')
            return self.hist
        
        xxtemp =  np.int_(np.around(((binX-1)*((xvar-Xmin)/(Xmax-Xmin)))))
        yytemp =  np.int_(np.around(((binY-1)*((yvar-Ymin)/(Ymax-Ymin)))))
             
        for k in range(len(xvar)):
         
            xx =  xxtemp[k]
            yy =  yytemp[k]
        
            if self.outBounds == True:
                
               if ( (xx >= 0) and (xx <= binX-1) and (yy >= 0) and (yy <= binY-1) ):
                   self.hist[yy,xx] += 1
               elif ( (xx >= 0) and (xx > binX-1) and (yy >= 0) and (yy <= binY-1) ):
                   self.hist[yy,-1] += 1
               elif ( (xx < 0) and (xx <= binX-1) and (yy >= 0) and (yy <= binY-1) ):
                    self.hist[yy,0] += 1
               elif ( (xx >= 0) and (xx <= binX-1) and (yy < 0) and (yy <= binY-1) ):
                   self.hist[0,xx]  += 1
               elif ( (xx >= 0) and (xx <= binX-1) and (yy >= 0) and (yy > binY-1) ):
                   self.hist[-1,xx] += 1
                   
            elif self.outBounds == False:
                 
               if ( (xx >= 0) and (xx <= binX-1) and (yy >= 0) and (yy <= binY-1) ):
                  self.hist[yy,xx] += 1
               else:
                   if cont == 0:
                       print('\033[1;33mWARNING: hist out of bounds.\033[1;37m') 
                       cont = 1  
                          
        return self.hist
    
    ###############################################################################
    ############################################################################### 
    
    def histXYZ (self,XX,A,YY,B,ZZ,C):
 
        # self.coincidence = coincidence
        # self.showStats   = showStats
        #  you pass the vectors XX, YY, ToFx, POPH
                   
        # A = POPH[:,0]
        # B = POPH[:,1]
        # C = POPH[:,2]
               
        # if self.coincidence is True:
        #    print('\t building histograms ... coincidence ON ...')
        # else:
        #    print('\t building histograms ... coincidence OFF ...')
           
        binX   = len(XX) 
        binY   = len(YY)
        binZ   = len(ZZ)
            
        Xmin   = min(XX) 
        Xmax   = max(XX) 
        Ymin   = min(YY) 
        Ymax   = max(YY) 
        Zmin   = min(ZZ)
        Zmax   = max(ZZ) 
               
        self.XY     = np.zeros((binY,binX)) 
        self.XYproj = np.zeros(binX) 
        self.XZ     = np.zeros((binX,binZ)) 
             
        count   = np.zeros((3,2)) #counter for rejected and good evetns
               
        if (not( (len(A) == len(B)) and (len(A) == len(C)))):
            print('\n \t \033[1;31m----> ABORT: X and/or Y and/or T not same length! \033[1;37m\n')
            return self.XY, self.XYproj, self.XZ
    
        xxtemp  =  np.int_(np.around(((binX-1)*((A-Xmin)/(Xmax-Xmin)))))
        yytemp  =  np.int_(np.around(((binY-1)*((B-Ymin)/(Ymax-Ymin)))))
        zztemp  =  np.int_(np.around(((binZ-1)*((C-Zmin)/(Zmax-Zmin)))))
             
        Nev = len(A)
            
        for k in range(0,Nev,1):
                
                xx =  xxtemp[k]
                yy =  yytemp[k]
                zz =  zztemp[k]
                
                #####################################
                if ( (xx >= 0) and (xx <= binX-1) ):
                    
                    # 2D hist X-Y
                    if ( (yy >= 0) and (yy <= binY-1) ):
                        self.XY[yy,xx]  += 1
                        count[0,0] += 1  # if 2D
                    else:
                        count[0,1] += 1  # if 1D
                   
                    # 1D hist X
                    self.XYproj[xx]  += 1
                    count[1,0]  += 1   # if at least 1D
    
                    # # 2D hist X-ToF or X-Lambda
                    # if coincidence == 1:
                    #      if ( (yy >= 0) and (yy <= binY-1) and (zz >= 0) and (zz <= binZ-1) ):
                    #          self.XZ[xx,zz]  += 1
                    #          count[2,0] += 1
                    #      else:
                    #          count[2,1] += 1
                    # elif coincidence == 0:
                    #      if ( (zz >= 0) and (zz <= binZ-1) ):
                    #          self.XZ[xx,zz]  += 1
                    #          count[2,0] += 1
                    #      else:
                    #          count[2,1] += 1
                             
                    if ( (zz >= 0) and (zz <= binZ-1) ):
                             self.XZ[xx,zz]  += 1
                             count[2,0] += 1
                    else:
                             count[2,1] += 1
                             
                else:
                     count[1,1] += 1
                #####################################
                     
        # countn = 100*(count/Nev)
           
        if count[1,1] != 0 :
           print('\n \t \033[1;33mWARNING: %.1f%% out of 1D boundaries \033[1;37m\n' % count[1,1])
    
        # if self.showStats is True:
        #    print(' \t percentage in 2D hist: %.1f%%, out of boundaries or only w %.1f%%' % (countn[0,0],countn[0,1]))
        #    print(' \t percentage in proj hist: %.1f%% (only w), out of boundaries %.1f%%' % (countn[1,0],countn[1,1]))
        #    print(' \t percentage in ToF/Lambda hist: %.1f%%, out of boundaries or only w %.1f%%' % (countn[2,0],countn[2,1]))
        
        return self.XY, self.XYproj, self.XZ
    
###############################################################################
############################################################################### 
    
if __name__ == '__main__' :
    
    configFilePath  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/config/'+"MB300_AMOR_config.json"
    filePathD       = './'+"VMM3a_Freia.pcapng"
   
    parameters  = para.parameters('/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/')
    config = maps.read_json_config(configFilePath)
   
    parameters.loadConfigParameters(config)
    
    parameters.cassettes.cassettes = [1,2]
    
    # parameters.plotting.positionReconstruction = 1
    
    # # parameters.update()
    
    # ax = allAxis()
    
    # ax.axEnergy.start = 0
    # ax.axEnergy.stop  = 65535
    # ax.axEnergy.steps = 128
    
    # ax.axLambda.stop  = 12
    # ax.axLambda.steps = 32
    
    # ax.updateAllAxis()
    
    # ax.axEnergy.update()
    
    # ax.generateAllAxis(config, cassettes)
    
    # ax.axEnergy.stop = 1000
    # ax.updateAllAxis()
    
    # ax.createAllAxis(parameters)
    
    # aa = ax.axStrips_mm.axis
    
    # print(aa)
    
    Nhits = 1e4
    cassettes1 = [1,2,3,4]
      
    bb = sdat.sampleHitsMultipleCassettes(cassettes1)
    bb.generateGlob(Nhits)
    hits = bb.hits

    cc = clu.clusterHits(hits,'globalStat')
    cc.clusterizeManyCassettes(parameters.cassettes.cassettes, parameters.dataReduction.timeWindow)
    events = cc.events
    
    vv = absu.calculateAbsUnits(events, parameters)
    vv.calculatePositionAbsUnit()
    vv.calculateToFandWavelength(0)
    events = vv.events
    
    
    allAxis = hh.allAxis()
    allAxis.createAllAxis(parameters)
    
    selc = events.positionS >= 0
    
    h2D, hProj, hToF = hh.histog().histXYZ(allAxis.axWires.axis, events.positionW, allAxis.axStrips.axis,events.positionS, allAxis.axToF.axis, events.ToF,coincidence=True,showStats=True)
    
    hProjAll = hh.histog().hist1D(allAxis.axWires.axis, events.positionW)
    
    h2D, hProj, hToF = hh.histog().histXYZ(allAxis.axWires.axis, events.positionW[selc], allAxis.axStrips.axis,events.positionS[selc], allAxis.axToF.axis, events.ToF[selc],coincidence=True,showStats=True)
    
    
    # pp = plo.plottingEvents(eventsAT,allAxis)
    # pp.plotXYToF(logScale=False, absUnits = False)
    
    