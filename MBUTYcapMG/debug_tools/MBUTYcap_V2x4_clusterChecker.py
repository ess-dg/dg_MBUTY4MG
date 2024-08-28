#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
###############################################################################
########    V2.4 2023/09/13      francescopiscitelli     ######################
###############################################################################
###############################################################################

import numpy as np
import time
import os
import sys
import matplotlib.pyplot as plt

# import matplotlib
# # matplotlib.use(‘Qt5Agg’)
from PyQt5.QtWidgets import QApplication, QFileDialog, QDialog, QGridLayout, QLabel, QLineEdit
app = QApplication(sys.argv)

### import the library with all specific functions that this code uses 
from lib import libReadPcapngVMM as pcapr
from lib import libSampleData as sdat
from lib import libMapping as maps
from lib import libCluster as clu
from lib import libAbsUnitsAndLambda as absu
from lib import libHistograms as hh
from lib import libFileManagmentUtil as fd
from lib import libParameters as para
from lib import libTerminal as ta
from lib import libPlotting as plo
from lib import libEventsSoftThresholds as thre
from lib import libReducedFileH5 as saveH5

###############################################################################
###############################################################################

# STILL TO IMPLEMENT:
#     - save reduced data - almost done 
#     - monitor lambda
#     - TDC and ADC calibration 
#     - now the mon events  stay in the events array,they need  to be taken out 

###############################################################################
###############################################################################
profiling = para.profiling()
print('----------------------------------------------------------------------')
print('\033[1;32mCiao '+os.environ['USER']+'! Welcome to MBUTY 2.4\033[1;37m')
print('----------------------------------------------------------------------')
plt.close("all")
### check version ###
para.checkPythonVersion()
### get current path ###
currentPath = os.path.abspath(os.path.dirname(__file__))+'/'
parameters  = para.parameters(currentPath)
###############################################################################
###############################################################################
### read json and create parameters for plotting and analisys ###

configFilePath  = currentPath+'config/'

configFileName  = "test.json"

###############################################################################
###############################################################################
### read json and create parameters for plotting and analisys ###
config = maps.read_json_config(configFilePath+configFileName)
parameters.loadConfigAndSetParameters(config)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
### edit parameters for plotting and analisys here ###
###############################################################################
### FILE MANAGMENT  PARAMETERS:
#################################

### ON/OFF if you want to rsync the data     
parameters.fileManagement.sync = False  

### from ... to  ... rsync the data

# test setup
parameters.fileManagement.sourcePath = 'essdaq@172.30.244.233:~/pcaps/'

parameters.fileManagement.sourcePath = 'essdaq@172.30.244.203:~/pcaps/'

# parameters.fileManagement.sourcePath = 'essdaq@192.168.0.1:~/pcaps/'

parameters.fileManagement.destPath   = '/Users/francescopiscitelli/Desktop/dataVMM/'
# parameters.fileManagement.destPath   = currentPath+'data/'
# parameters.fileManagement.destPath   = '/Users/francescopiscitelli/Desktop/dataVMM_old/'


# parameters.fileManagement.destPath  = '/Users/francescopiscitelli/Desktop/ECDC_AMOR/'


# parameters.fileManagement.destPath   = '/Users/francescopiscitelli/Documents/DOC/DATA/202308_Utgard_MBnewAMOR_VMM_muons/DATA/'

###############

parameters.fileManagement.filePath = '/Users/francescopiscitelli/Documents/DOC/DATA/202311_PSI_AMOR_MBnewAMOR_VMM_neutrons/SamplesAndMasks/'

### folder and file to open (file can be a list of files)
# parameters.fileManagement.filePath = parameters.fileManagement.destPath
parameters.fileManagement.fileName = ['freia_1k_pkts_ng.pcapng']
parameters.fileManagement.fileName = ['freiatest.pcapng']
parameters.fileManagement.fileName = ['20231106_142811_duration_s_5_YESneutrons1240K1070Rth280_maskESS_00000.pcapng']

# parameters.fileManagement.fileSerials = np.arange(18,28,1)

### valid otions: 'window','fileName', 'latest', 'secondLast', 'wholeFolder', 'sequence' 
### window opens to selcet file, filename speficified  earlier, last or sencond last file crearted in folder, 
### entire  folder  opend  and analized and cumulated  all togheter 
### sequence opens all filens in     parameters.fileManagement.fileSerials and with fileName
parameters.fileManagement.openMode = 'window'
parameters.fileManagement.openMode = 'fileName'
# parameters.fileManagement.openMode = 'latest'
# parameters.fileManagement.openMode = 'secondLast'
# parameters.fileManagement.openMode = 'wholeFolder'
# parameters.fileManagement.openMode = 'sequence'

###############
### path to threshold  file
parameters.fileManagement.thresholdFilePath = parameters.fileManagement.currentPath+'config/'
parameters.fileManagement.thresholdFileName = 'MB300L_thresholds.xlsx'

###############
### path to  Tshark, in case you open a pcap  it gets converted into pcapng 
parameters.fileManagement.pathToTshark = '/Applications/Wireshark.app/Contents/MacOS/'

###############
### save a hdf file with clusters (reduced file)

### ON/OFF
parameters.fileManagement.saveReducedFileONOFF = False   
parameters.fileManagement.saveReducedPath = '/Users/francescopiscitelli/Desktop/reducedFile/'

parameters.fileManagement.reducedNameMainFolder  = 'entry1'
parameters.fileManagement.reducedCompressionHDFT  = 'gzip'  
parameters.fileManagement.reducedCompressionHDFL  = 9    # gzip compression level 0 - 9

###############################################################################
### ANALISYS PARAMETERS:
#################################

### cassettes to clusterize and to plot, if empty the cassettes in the config file are taken as default
# parameters.cassettes.cassettes = [1,2,3,4,5,6]
# parameters.cassettes.cassettes = np.arange(32,0,-1)

### timeWindow to search for clusters, timeWindow is max time between events in candidate cluster 
### and timeWindow/2 is the recursive time distance between adjacent hits
parameters.dataReduction.timeWindow = 0.5e-6

### 'OFF', 'fromFile' = File With Threhsolds Loaded, 'userDefined' = User defines the Thresholds in an array softTh
parameters.dataReduction.softThresholdType = 'off' 
# parameters.dataReduction.softThresholdType = 'fromFile' 

if parameters.dataReduction.softThresholdType == 'userDefined':
    
    parameters.dataReduction.createThArrays(parameters.cassettes.cassettes, parameters)    
    parameters.dataReduction.softThArray.ThW[:,0] = 15000
    parameters.dataReduction.softThArray.ThS[:,0] = 500   
          
###############################################################################
### WAVELENGTH PARAMETERS:
#################################

### distance in mm from chopper and wires 0 of detector
parameters.wavelength.distance = 8000

##ON/OFF
parameters.wavelength.calculateLambda = False

### ON/OFF plot X vs Lambda 2D plot
parameters.wavelength.plotXLambda   = False
### ON/OFF integrated over single cassettes
parameters.wavelength.plotLambdaDistr = False

parameters.wavelength.lambdaBins  = 128
parameters.wavelength.lambdaRange = [1, 16]   #A

parameters.wavelength.chopperPeriod = 0.15 #s (NOTE: only matters if multipleFramesPerRest > 1)

### if chopper has two openings or more per reset of ToF
parameters.wavelength.multipleFramePerReset = False  #ON/OFF (this only affects the lambda calculation)
parameters.wavelength.numOfBunchesPerPulse  = 2
parameters.wavelength.lambdaMIN             = 2.5     #A

### in seconds, time shift betweeen pickup and chopper edge 
parameters.wavelength.chopperPickUpDelay =  13.5/(2.*180.) * parameters.wavelength.chopperPeriod/parameters.wavelength.numOfBunchesPerPulse 

  
###############################################################################
### MONITOR PARAMETERS:
#################################

### ON/OFF
parameters.MONitor.MONOnOff = True   

### threshold on MON, th is OFF if 0, any other value is ON
parameters.MONitor.MONThreshold = 0   

### ON/OFF plotting (MON ToF and Pulse Height) 
parameters.MONitor.plotMONtofPHS = True  

### in mm, distance of MON from chopper if plotMONtofPH == 1 (needed for lambda calculation if ToF)
parameters.MONitor.MONDistance  = 0   

###############################################################################
### PLOTTING PARAMETERS:
#################################

###############     
### show stat during clustering, option  'globalStat'  stat for all cassettes together, 
### 'individualStat' stat per cassette or None for no stat
parameters.plotting.showStat = 'globalStat'
# parameters.plotting.showStat = 'individualStat'

###############     
### raw plots
parameters.plotting.plotRawReadouts         = True
parameters.plotting.plotReadoutsTimeStamps  = False
parameters.plotting.plotRawHits             = False
parameters.plotting.plotHitsTimeStamps      = False
parameters.plotting.plotHitsTimeStampsVSChannels = False
parameters.plotting.plotChopperResets       = False

# with True disables clustering and mapping for speed reasons, analisys stops at readouts 
bareReadoutsCalc = False

###############
### Instantaneous Rate
parameters.plotting.plotInstRate    = False
parameters.plotting.instRateBin     = 100e-6  # s
 
###############
### ToF plot integrated over individual cassette, one per cassette
parameters.plotting.plotToFDistr    = False

parameters.plotting.ToFrange        = 0.12    # s
parameters.plotting.ToFbinning      = 100e-6 # s
     
parameters.plotting.plotMultiplicity = False 

# parameters.configJsonFile.orientation = 'vertical'

### 'W.max-S.max' is max max,  'W.cog-S.cog' is CoG CoG, 'W.max-S.cog' is wires max and strips CoG 
parameters.plotting.positionReconstruction = 'W.max-S.cog'
parameters.plotting.positionReconstruction = 'W.max-S.max'

### if True plot XY and XtoF plot in absolute unit (mm), if False plot in wire and strip ch no.
parameters.plotting.plotABSunits = False
 
### plot XY and XToF in log scale 
parameters.plotting.plotIMGlog = True

### ON/OFF, if  Tof  and Lambdaplot needs to include only events with strip present (2D) is True otherwise all events also without strip set to False
parameters.plotting.coincidenceWS_ONOFF = True

### ON/OFF, if  invalid ToFs Tofare included in the plots or removed from events 
parameters.plotting.removeInvalidToFs = True

### histogram outBounds param set as True as default (Events out of bounds stored in first and last bin)
parameters.plotting.hitogOutBounds = True

### for absolute  units calculation in X, already loaded in config json file
# parameters.configJsonFile.offset1stWires = 10

##############################      
### PHS

### ON/OFF PHS per channel and global
parameters.pulseHeigthSpect.plotPHS = True

### plot PHS in log scale 
parameters.pulseHeigthSpect.plotPHSlog = True

parameters.pulseHeigthSpect.energyBins = 128
parameters.pulseHeigthSpect.maxEnerg   = 1440

### plot the PHS correaltion wires vs strips
parameters.pulseHeigthSpect.plotPHScorrelation = False

###############################################################################
###############################################################################
########    end of with all the settings you can choose   #####################
########        DO NOT EDIT BELOW THIS LINE!!!!           #####################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
### if you want to save reduced data, it must include lambda, 
### so lambda calculation is turned ON if not yet 
if parameters.fileManagement.saveReducedFileONOFF is True:       
    parameters.wavelength.calculateLambda = True
    print('\n \t Lambda calculation turned ON to save reduced DATA')
    
if parameters.wavelength.plotXLambda or parameters.wavelength.plotLambdaDistr is True:
    parameters.wavelength.calculateLambda = True
    print('\n \t Lambda calculation turned ON to plot Lambda')
###############################################################################
###############################################################################
parameters.update()
###############################################################################
###############################################################################
### create axes for plotting and histograms 

allAxis = hh.allAxis()
allAxis.createAllAxis(parameters)

### here eventually edit axes 
### then update axes with: allAxis.updateAllAxis()

###############################################################################
###############################################################################
### sync the data folder from remote computer to local folder 
if parameters.fileManagement.sync is True:
    transferData = ta.transferDataUtil()
    transferData.syncData(parameters.fileManagement.sourcePath, parameters.fileManagement.destPath)

###############################################################################
###############################################################################
### select data
fileDialogue = fd.fileDialogue(parameters)
fileDialogue.openFile()
###############################################################################
###############################################################################

### initialize readouts cumulated over file list
readouts = pcapr.readouts()
# hits   = maps.hits()
# events = clu.events()


for cont, fileName in enumerate(fileDialogue.fileName):

    print('\033[1;32m-> reading file {} of {}\033[1;37m'.format(cont+1,len(fileDialogue.fileName)))
    
    ### check if a file is pcapng otherwise pcap is converted into pcapng
    conv = ta.pcapConverter(parameters)
    conv.checkExtensionAndConvertPcap(fileDialogue.filePath+fileName)
    if conv.flag is False:
        fileName = conv.fileName_OUT
        
    ### check which Ring, Fen and Hybrid is present in the selected File 
    # pcapr.checkWhich_RingFenHybrid_InFile(fileDialogue.filePath+fileName,parameters.clockTicks.NSperClockTick).check()
    
    ### load data  
    pcap = pcapr.pcapng_reader(fileDialogue.filePath+fileName, parameters.clockTicks.NSperClockTick, config.MONmap.TTLtype, config.MONmap.RingID, timeResolutionType='fine', sortByTimeStampsONOFF = False)
    readouts.append(pcap.readouts)
    
    # md  = maps.mapDetector(pcap.readouts, config)
    # md.mappAllCassAndChannelsGlob()
    # hits.append(md.hits)
     
    # cc = clu.clusterHits(md.hits,parameters.plotting.showStat)
    # cc.clusterizeManyCassettes(parameters.cassettes.cassettes, parameters.dataReduction.timeWindow)
    # eve.append(cc.events)
    


      
####################    
### for debug, generate sample readouts
# aa = sdat.sampleReadouts_2()
# aa.fill()
# readouts = aa.readouts
####################

####################
### for debug, readouts in single array
readoutsArray = readouts.concatenateReadoutsInArrayForDebug()
####################

if bareReadoutsCalc is False:
    
    
    md  = maps.mapDetector(readouts, config)
    md.mappAllCassAndChannelsGlob()
    hits = md.hits
    hitsArray  = hits.concatenateHitsInArrayForDebug()


    ####################
    ### getting the hits for the MONitor
    if parameters.MONitor.MONOnOff is True:
        MON = maps.mapMonitor(readouts, config)
        if MON.flagMONfound is True:
            hitsMON = MON.hits
            
            MONe = clu.hitsMON2events(hitsMON)
            eventsMON = MONe.events
            
            abMON = absu.calculateAbsUnits(eventsMON, parameters, 'MON')
            abMON.calculateToF(parameters.plotting.removeInvalidToFs)
            
            print('\033[1;32m\t MON events: {}\033[1;37m'.format(len(eventsMON.timeStamp)))
            
            
            # CALCULATION OF LAMBDA ON MON NOT YET IMPLEMENTED



    
    ###############################################################################
    ### map data
    # md  = maps.mapDetector(readouts, config)
    # md.mappAllCassAndChannelsGlob()
    # hits = md.hits
    
    # for debug force all hits in a single cassetteno.1 even if from different hybrids
    # hits.Cassette = np.ones(len((hits.Cassette)),dtype='int64')
    
     
    ####################    
    ### for debug, generate sample hits 
    # Nhits = 5e4
    # bb = sdat.sampleHitsMultipleCassettes()
    # bb.generateGlob(Nhits)
    # hits = bb.hits   
    
    # bb = sdat.sampleHitsMultipleCassettes_2()
    # bb.generateGlob()
    # hits = bb.hits  
    ######################
    
    ####################    
    ### for debug, hits in single array 
    # hitsArray = hits.concatenateHitsInArrayForDebug()
    ####################
    
    ###############################################################################
    ### clusterize
    cclu = clu.clusterHits(hits,parameters.plotting.showStat)
    # cclu.clusterizeManyCassettes(parameters.cassettes.cassettes, parameters.dataReduction.timeWindow)
    events = clu.events()
    
    cassettesIDs = parameters.config.DETparameters.cassInConfig
    
    timeWindow = parameters.dataReduction.timeWindow
    
    showStat = parameters.plotting.showStat
    
    intervals = 4
    
    rejCounterAll = np.zeros((5),dtype='int64')
    
    # for cc in cassettesIDs:
        
    for cc in [0]:
        
        # self.clusterize1cassette(cc, timeWindow)
        
        cassette1ID = cc
        
        present = clu.checkCassIDs.checkIfPresentInHits(hits,cassette1ID)
        
        # init an empty event obj
        events1Cass  = clu.events()
        
        
        
        rejCounter = np.zeros((5),dtype='int64')
        
        if present:

            print('\n \t clustering ... Cassette ID ',str(cassette1ID),' ',end='')
           
            # this is a trick to accept also the clusters very close in time otherwise rejected
            # Timewindowrec = np.ceil(timeWindow*1e6/3)/1e6+0.01e-6
            # timeWindow    = timeWindow+0.01e-6
            # changed later to 
           
            # this is a trick to accept also the clusters very close in time otherwise rejected
            # and also timeWindow defines the max time distance from first event and last event in a 
            # candidate cluster.
            # the recursive is the delta time between adjiacdent event in candidate but timeWindow is the max allowed 
            # gap for first and last 
            # and also add 1% to avoid rounding issues *1.01
            # timeWindow    = timeWindow*1.01
            # TimewindowRecursive = timeWindow/2
           
            # Time window is in s but we all do here in ns int64
            # timeWindow converted in ns
            timeWindow_ns = int(round(timeWindow*1e9))
       
            TimeWindowMax          = int(round(timeWindow_ns*1.01))
            TimeWindowRecursive    = int(round(TimeWindowMax/2))
       
            ##########
        
            #select only the data relative to the wanted cassette ID
            selectCassette = hits.Cassette == cassette1ID
            
        
            # for speed the hits are inserted in an array
            # add a line at top [0,0,0,0] not to lose the 1st event
            data = np.zeros((np.shape(hits.timeStamp[selectCassette])[0],6), dtype = 'int64') 
            data[:,0] = hits.timeStamp[selectCassette]
            data[:,1] = hits.WiresStrips[selectCassette]
            data[:,2] = hits.ADC[selectCassette]
            data[:,3] = hits.WorS[selectCassette]
            data[:,4] = hits.PulseT[selectCassette]
            data[:,5] = hits.PrevPT[selectCassette]
            
            data = data[0:100,:]
            

            # add a line at top [0,0,0,0] not to lose the 1st event
            data = np.concatenate( ( np.zeros((1,np.shape(data)[1]), dtype = 'int64'), data ), axis=0)  #add a line at top not to lose the 1st event
            data[0,0] = -2*TimeWindowMax

            # data[:,0] = np.around(data[:,0],decimals=self.resolution) #time rounded at 1us precision is 6 decimals, 7 is 100ns, etc...

            deltaTime = np.diff(data[:,0])                     #1st derivative of time 
            deltaTime = np.concatenate(([0],deltaTime),axis=0) #add a zero at top to restore length of vector
       
            clusterlogic = (np.absolute(deltaTime) <= TimeWindowRecursive) #is zero when a new cluster starts 
            
            
       
            data1 = np.concatenate((data,clusterlogic[:,None]),axis=1) #this is for debugging 
       
            index = np.argwhere(clusterlogic == False) #find the index where a new cluster may start 
           
            #################################
   
            ADCCH = np.zeros((np.shape(data)[0],12),dtype='int64')
    
            # ADCCH2 = np.zeros((np.shape(data)[0],1),dtype='int64')
    
            ADCCH[:,0:3] = (data[:,0:3])  # first 3 columns as data
            ADCCH[:,3]   = (clusterlogic.astype(int)) # col 3 is 0 where a new cluster may start
    
            ADCCH[:,4]   = ((data[:,3] == 0).astype(int))   # wire  
            ADCCH[:,5]   = ((data[:,3] == 1).astype(int))   # strip 
    
            ADCCH[:,6]   = (data[:,1]*ADCCH[:,4]).astype(int)   # wire ch
            ADCCH[:,7]   = (data[:,1]*ADCCH[:,5]).astype(int)   # strip ch

            ADCCH[:,8]   = (data[:,2]*ADCCH[:,4]).astype(int)   # wire ADCs 
            ADCCH[:,9]   = (data[:,2]*ADCCH[:,5]).astype(int)   # strip ADCs 
    
            ADCCH[:,10]  =  (ADCCH[:,4]*ADCCH[:,6]*ADCCH[:,8]).astype(int)    # weighted position on wires
            ADCCH[:,11]  =  (ADCCH[:,5]*ADCCH[:,7]*ADCCH[:,9]).astype(int)   # weighted position on strips
           
            #################################

            NumClusters = np.shape(index)[0]
       
            events1Cass.Nevents = NumClusters
       
            TPHM = np.zeros((NumClusters,9),dtype='int64')  #output data with col0 position wires, col1 poisiton strips, col2 tof, col3 pulse height wires, col4 pulse height strips, col 5 multiplicity w, col 6 muiltiplicity strips
            PO   = np.zeros((NumClusters,9),dtype='float64')
           
            # filling timeStamp column
            TPHM[:,0]  = data[index[:,0],0]   # timeStamp      
            TPHM[:,1]  = data[index[:,0],4]   # PulseT   
            TPHM[:,2]  = data[index[:,0],5]   # PrevPT
    
            #################################
            # add a fake last cluster to make loop up to the very last true cluster
            index = np.concatenate((index,[[np.shape(data)[0]]]),axis=0)
            
            ADCCH = np.concatenate((ADCCH,np.zeros((1,12),dtype='int64')),axis=0)

     
            #################################
            


            #################################
            if  NumClusters >= 0:
               
                for kk in range(0,NumClusters,1):
                                       
                        steps = int(NumClusters/intervals)+1    # equivalent as ceil rounding 
                        if np.mod(kk,steps) == 0 or kk == (NumClusters-1):
                            percents = int(round(100.0 * kk / float(NumClusters), 1))
                            print('['+format(percents,'01d') + '%]',end=' ')
                       
                        clusterq = ADCCH[index[kk,0]:index[kk+1,0],:]
   
   
                        # print('cass: '+str(cassette1ID),' event candi: '+str(kk))
                        # print(np.around(clusterq,1))
                       
                        acceptWindow = ((clusterq[-1,0] - clusterq[0,0]) <= TimeWindowMax)  #max difference in time between first and last in cluster 
                      
                        clusterq = clusterq[clusterq[:,1].argsort(kind='quicksort'),:]  #order cluster by ch number
                  
                        is_wire  = clusterq[:,4] == 1
                        is_strip = clusterq[:,5] == 1
                       
                        # n wires n strips in cluster
                        ww = len(clusterq[is_wire, 4])   #num of wires in cluster
                        ss = len(clusterq[is_strip, 5])  #num of strips in cluster
           
                       
                        if (ww != 0 and ss != 0 and ss <= 63 and ww <= 31 and acceptWindow): #if there is at least 1 wire and 1 strip and no ch number above 32
           
                            mmaxw = clusterq[is_wire, 6][-1]
                            mmaxs = clusterq[is_strip, 7][-1]
                            mminw = clusterq[is_wire, 6][0]
                            mmins = clusterq[is_strip, 7][0]
               
                            neigw = (mmaxw - mminw) == (ww-1) #if event repated is rejected because neigw is 1 even if the same wire is repeated and should be 2 
                            neigs = (mmaxs - mmins) == (ss-1)
                           
                            if (neigw == 1 and neigs == 1):    #if they are neighbour then...
                               
                                rejCounter[0] = rejCounter[0]+1   #counter 2D
                               
                                TPHM[kk,3]   = ww     #multiuplicity wires
                                TPHM[kk,4]   = ss     #multiuplicity strips
                                TPHM[kk,5]   = np.sum(clusterq[:,8],axis=0)   #PH wires
                                TPHM[kk,6]   = np.sum(clusterq[:,9],axis=0)   #PH strips
                                PO[kk,0]     = round((np.sum(clusterq[:,10],axis=0))/(TPHM[kk,5]),2)       #position wires 0 to 31
                                PO[kk,1]     = round((((np.sum(clusterq[:,11],axis=0))/(TPHM[kk,6]))),2)   #position strips from 0 to 31 or up to 63
                                      
                            else:
                                rejCounter[1] = rejCounter[1]+1;                #counter if they are no neighbour 
                               
                        elif (ww >= 1 and ss == 0 and ww <= 31 and acceptWindow): #put in 1D hist only for wires when there is no strip 
                               
                            mmaxw = clusterq[is_wire, 6][-1]
                            mminw = clusterq[is_wire, 6][0] 
               
                            neigw = (mmaxw - mminw) == (ww-1)    #if event repated is rejected because neigw is 1 even if the same wire is repeated and should be 2 
                          
                            if (neigw == 1):    #if they are neighbour then...
                   
                                rejCounter[2] = rejCounter[2]+1                #counter 1D

                                TPHM[kk,3]   = ww     #multiuplicity wires
                                TPHM[kk,5]   = np.sum(clusterq[:,8],axis=0)   #PH wires
                                PO[kk,0]     = round((np.sum(clusterq[:,10],axis=0))/(TPHM[kk,5]),2)     #position wires
                                PO[kk,1]     = -1 #position strips if absent
                                  
                            else:
                                rejCounter[1] = rejCounter[1]+1              #counter if they are no neighbour 
                               
                        elif (ww >= 32 or ss >= 64):
                              rejCounter[3] = rejCounter[3]+1               #counter if cluster above possible limits          
                              print('\n cluster > 32 in either directions w or s -> probably rate too high \n')
                            
                        else:
                            rejCounter[4] = rejCounter[4]+1               #any other case not taken into account previously
                           
               
                # print('\n')        
                  
                rejected = np.logical_and((TPHM[:,3] == 0),(TPHM[:,4] == 0))    #remove rejected from data in rejCoiunter[4] it is when only strips and wire and sgtrip mult is 0, whole row in POPH is 0 actually 
                   
                TPHM     = TPHM[np.logical_not(rejected),:]    #remove rejected from data
                PO       = PO[np.logical_not(rejected),:]      #remove rejected from data
               
                # check if there are some NaN 
                theNanOnes = np.logical_or(np.isnan(PO[:,0]),np.isnan(PO[:,1])) 
               
                if np.sum(theNanOnes) > 0:
                    print('\n --> \033[1;33mWARNING: Some NaN (invalid value) values found in clustering -> removed from events.\033[1;37m',end='')

               
                TPHM     = TPHM[np.logical_not(theNanOnes),:]    #remove nan from data
                PO       = PO[np.logical_not(theNanOnes),:]      #remove nan from data
               
                # self.events.NeventsNotRejAll = self.events.Nevents - (self.rejCounter[1]+self.rejCounter[3]+self.rejCounter[4]);
                # self.events.NeventsNotRej2D   = np.sum(self.POPH[:,1] >= 0)
               
                events1Cass.NeventsNotRejAll = rejCounter[0] + rejCounter[2]
                events1Cass.NeventsNotRej2D  = rejCounter[0]
               
                events1Cass.transform(cassette1ID, TPHM, PO)
       
        ################################
                if showStat == 'individualStat':
                    cclu.someStat(events1Cass,rejCounter)
        ################################
       
                del TPHM
                del PO
               
                if NumClusters == len(hits.timeStamp[selectCassette]):
                    print('--> \033[1;33mWARNING: time window for clustering might be too small!\033[1;37m',end='')
                    time.sleep(1)
                   
                   
        
        events.append(events1Cass)
        
        rejCounterAll += rejCounter
        
        # del events1Cass
        
    if showStat == 'globalStat':
        cclu.someStat(events,rejCounterAll)
    
    ####################    
    ### for debug, events in single array 
    # eventsArray = events.concatenateEventsInArrayForDebug() 
    ####################
    events1CassArray = events1Cass.concatenateEventsInArrayForDebug()
    ####################    
    ### for debug, generate sample events 2
    # dd = sdat.sampleEventsMultipleCassettes(parameters.cassettes.cassettes,'./data/')
    # dd.generateGlob(Nevents)
    # events  = dd.events
    # eventsArray = events.concatenateEventsInArrayForDebug()
    # eventsArray = eventsArray[72:100,:]
    #################### 
        
    
    ###############################################################################
    ### calculate abs units, ToF and lambda
    
    ab = absu.calculateAbsUnits(events, parameters)
    ab.calculatePositionAbsUnit()
    
    ab.calculateToF(parameters.plotting.removeInvalidToFs)
    
    if parameters.wavelength.calculateLambda is True:
        ab.calculateWavelength()
    
    # eventsBTh = ab.events 
    
    events = ab.events 
    
    ####################    
    ### for debug, events in single array 
    # eventsArrayBTh = eventsBTh.concatenateEventsInArrayForDebug() 
    # eventsArrayBTh = events.concatenateEventsInArrayForDebug() 
    ####################    
    
    ###############################################################################
    ###############################################################################
    ### software thresholds
    
    # asth = thre.applyThresholdsToEvents(eventsBTh, parameters.cassettes.cassettes, parameters, parameters.plotting.showStat)
    # asth.thresholdizeAllCassettes()
        
    # events  = asth.events 
    
    # eventsArrayAT = events.concatenateEventsInArrayForDebug()
    
    ###############################################################################
    ###############################################################################
    ### save reduced data to hdf5
      
    if parameters.fileManagement.saveReducedFileONOFF is True: 
        
        fileNameSave = 'temp'
    
        fileNameSave  = fileDialogue.fileName[0]+'_reduced'
        
        sav = saveH5.saveReducedDataToHDF(parameters,parameters.fileManagement.saveReducedPath,fileNameSave)
        # sav.save(eventsBTh)
        sav.save(events)

###############################################################################
###############################################################################
### plot
parameters.HistNotification()

######################
### readouts
if (parameters.plotting.plotRawReadouts  or  parameters.plotting.plotReadoutsTimeStamps) is True:
    plread = plo.plottingReadouts(readouts, config)
    if parameters.plotting.plotRawReadouts is True:
        plread.plotChRaw(parameters.config.DETparameters.cassInConfig)
    if parameters.plotting.plotReadoutsTimeStamps is True:
        plread.plotTimeStamps(parameters.config.DETparameters.cassInConfig)
        
if parameters.plotting.plotChopperResets is True:
    plread1 = plo.plottingReadouts(readouts, config)
    plread1.plotChoppResets()

         
######################

######################
if bareReadoutsCalc is False:
    ### hits
    if (parameters.plotting.plotRawHits or parameters.plotting.plotHitsTimeStamps or parameters.plotting.plotHitsTimeStampsVSChannels) is True:
        plhits = plo.plottingHits(hits, parameters)
        if parameters.plotting.plotRawHits is True:
            plhits.plotChRaw(parameters.config.DETparameters.cassInConfig)
        if parameters.plotting.plotHitsTimeStamps is True:
            plhits.plotTimeStamps(parameters.config.DETparameters.cassInConfig)
        if parameters.plotting.plotHitsTimeStampsVSChannels is True:    
            plhits.plotTimeStampsVSCh(parameters.config.DETparameters.cassInConfig)    
    ######################

######################
### events

if bareReadoutsCalc is False:
    ### XY and XToF
    plev = plo.plottingEvents(events,allAxis,parameters.plotting.coincidenceWS_ONOFF)
    plev.plotXYToF(logScale = parameters.plotting.plotIMGlog, absUnits = parameters.plotting.plotABSunits, orientation = parameters.config.DETparameters.orientation)
    
    # ### ToF per cassette 
    if parameters.plotting.plotToFDistr is True:
        plev.plotToF(parameters.cassettes.cassettes)

    ### lambda
    if parameters.wavelength.plotXLambda is True:
        plev.plotXLambda(logScale = parameters.plotting.plotIMGlog, absUnits = parameters.plotting.plotABSunits)
    ### lambda per cassette
    if parameters.wavelength.plotLambdaDistr is True:
        plev.plotLambda(parameters.cassettes.cassettes)
        
    ### multiplicity 
    if parameters.plotting.plotMultiplicity is True:
        plev.plotMultiplicity(parameters.cassettes.cassettes)

    # ### PHS
    if parameters.pulseHeigthSpect.plotPHS is True:
        plev.plotPHS(parameters.config.DETparameters.cassInConfig, parameters, logScale = parameters.pulseHeigthSpect.plotPHSlog)
    if parameters.pulseHeigthSpect.plotPHScorrelation is True:
        plev.plotPHScorrelation(parameters.cassettes.cassettes, parameters.pulseHeigthSpect.plotPHSlog)
    
    ### instantaneous Rate per cassette
    if parameters.plotting.plotInstRate is True:
        plev.plotInstantaneousRate(parameters.config.DETparameters.cassInConfig)

    ############
    # MON plots
    if parameters.MONitor.MONOnOff is True and parameters.MONitor.plotMONtofPHS is True and MON.flagMONfound is True:
        
        plMON = plo.plottingMON(eventsMON,allAxis)
        plMON.plot_ToF_PHS_MON()

###############################################################################
###############################################################################
#  any other plot that the user wants...

selc = events.positionS >= 0

h2D, _, hToF = hh.histog().histXYZ(allAxis.axWires.axis, events.positionW[selc], allAxis.axStrips.axis, events.positionS[selc], allAxis.axToF.axis, events.ToF[selc]/1e9)

np.sum(h2D)

np.sum(hToF)
# # sel = np.logical_or(hits.WiresStrips == 31,  hits.WiresStrips == 2)
# sel = np.logical_or(hits.WiresStrips == 0,  hits.WiresStrips == 8)

# diffe       = np.diff(hits.timeStamp[sel])

# # print('\n',np.sum(~sel))

# # sortedDiffe = diffe[diffe.argsort()]

# # # print('freq injected in 1 ch '+str(1/(sortedDiffe[-1]*1e-9))+' Hz')

# TD = np.arange(-2000,2000,1)

# hist  = hh.histog().hist1D(TD,diffe) 
# figl, ax = plt.subplots(num=567,figsize=(6,6), nrows=1, ncols=1) 
# ax.step(TD,hist,'b',where='mid')   
# ax.set_yscale('log')
# ax.set_xlabel('delta time (ns)')
# ax.set_ylabel('counts') 
# ax.grid()
# ax.set_xlim((-200,200))
# # figl.suptitle('2FEN_2Hy')
if bareReadoutsCalc is False:     
    figl566, ax566 = plt.subplots(num=56767,figsize=(6,6), nrows=1, ncols=1) 
    ax566.plot(np.arange(0,len(events.ToF)),events.ToF/1e9,'b')   
      
aa = eventsArray[eventsArray[:,3] != -1  , :]

###############################################################################
###############################################################################
plt.show()
profiling.stop()
print('----------------------------------------------------------------------')
###############################################################################
###############################################################################

