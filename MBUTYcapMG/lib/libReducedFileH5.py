#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 14:18:41 2021

@author: francescopiscitelli
"""

import numpy as np
import time
import os
import h5py
import sys


from lib import libReadPcapngVMM as pcapr
from lib import libSampleData as sdat
from lib import libMapping as maps
from lib import libCluster as clu
from lib import libAbsUnitsAndLambda as absu
from lib import libParameters as para

###############################################################################
###############################################################################

class saveReducedDataToHDF():
    def __init__(self, parameters, saveReducedPath='./', fileName='temp'):
        
        self.par = parameters
        
        self.nameMainFolder   = self.par.fileManagement.reducedNameMainFolder

        self.compressionHDFT  = self.par.fileManagement.reducedCompressionHDFT
        self.compressionHDFL  = self.par.fileManagement.reducedCompressionHDFL
        
        ###########################
        
        if not os.path.exists(saveReducedPath):
           os.makedirs(saveReducedPath)

        self.outfile = saveReducedPath+fileName+'.h5'
    
        # check if file already exist and in case yes delete it 
        if os.path.exists(self.outfile):
            print('\033[1;33mWARNING: Reduced DATA file exists, it will be overwritten!\033[1;37m')
            os.system('rm '+self.outfile)
      
        self.fid    = h5py.File(self.outfile, "w")
    
        # create groups in h5 file  
        self.gparam = self.fid.create_group(self.nameMainFolder+'/parameters')
        self.gdet   = self.fid.create_group(self.nameMainFolder+'/detector')
        self.gmon   = self.fid.create_group(self.nameMainFolder+'/monitor')
        
        # self.ginstr = self.fid.create_group(self.nameMainFolder+'/instrument')
        # self.grun   = self.fid.create_group(self.nameMainFolder+'/run')
        
        self.gdet_data = self.gdet.create_group('events')
        self.gmon_data = self.gmon.create_group('hits')
        
        padic =  self.par.__dict__
        for key in  padic.keys():
    
             if key != 'config':
                 
                 gparam_subg = self.gparam.create_group(key)
             
                 temp_dic = padic[key].__dict__
     
                 for key2, value in zip(temp_dic.keys(),temp_dic.values()) :
           
                        # print(key2, value, type(value))
                        
                        if isinstance(value, (bool, float, int)):
                            gparam_subg.create_dataset(key2, data = [value], compression=self.compressionHDFT, compression_opts=self.compressionHDFL)
                        elif isinstance(value, (str)):
                            gparam_subg.attrs.create(key2, value)
                        elif isinstance(value, (list, np.ndarray)):
                            try:
                                gparam_subg.create_dataset(key2, data = value, compression=self.compressionHDFT, compression_opts=self.compressionHDFL)
                            except:
                                
                                # val = value[0]
                                # print(val)
                                gparam_subg.attrs.create(key2,value)
                            
                        else:
                            print('excluded dataset: ->  '+str(key2)+' -> '+str(type(value)))
                            # a=1
                
         
        # 
                 
                 
        # for key, value in {
        #             'clockFreq': self.par.clockTicks.clockFreq,
        #             'NSperClockTick': self.par.clockTicks.NSperClockTick,
        #             'detectorName': self.par.configJsonFile.detectorName,
        #             'cassettesInConfig': self.par.configJsonFile.cassInConfig,
        #             'numOfCassettes': self.par.configJsonFile.numOfCassettes,
        #             'numOfWires': self.par.configJsonFile.numOfWires,
        #             'numOfStrips': self.par.configJsonFile.numOfStrips,
        #             'orientation': self.par.configJsonFile.orientation,
        #             'wirePitch': self.par.configJsonFile.wirePitch,
        #             'stripPitch': self.par.configJsonFile.stripPitch,
        #             'offset': self.par.configJsonFile.offset1stWires,
        #             'inclination': self.par.configJsonFile.bladesInclination,
        #             }.items():
        #     self.gdet_param.create_dataset(key, data=value)
            
        # for key, value in {
        #             'chopperFreq': self.par.wavelength.chopperFreq,
        #             'chopperPeriod': self.par.wavelength.chopperPeriod,
        #             'distance': self.par.wavelength.distance,
        #             }.items():
        #     self.ginstr.create_dataset(key, data=value)
    
        #for key, value in {
        #             'ToF-duration (s)': parameters.ToF,
        #             'DistanceAtWindow (mm)': DistanceAtWindow,
        #             'Distance (mm)': Distance,
        #             'DistanceSampleWindow (mm)': DistanceSampleWindow,
        #             'DistanceSample1stWire mm)' : DistanceSample1stWire,
        #             'PickUpTimeShift (s)': PickUpTimeShift,
        #             'BladeAngularOffset (deg)': BladeAngularOffset,
        #             'OffsetOf1stWires (mm)': OffsetOf1stWires,
        #             }.items():
        #     ginstr.attrs.create(key, value)
    
        
        self.gdet.attrs.create('units: time in ns (int), position in channel number or mm (phys. unit), lambda in Angstrom, Pulse Heigth in a.u. ADC',1)
    
        
        # self.__del__()  
         
    def __del__(self):
        try:
            self.fid.close()
        except:
            pass
 
    def save(self, events, hitsMON=None):
        
        print('-> saving reduced data to h5 file ... ')
        
        self.events  = events
    
        evdic = self.events.__dict__        
        for key, value in  zip(evdic.keys(),evdic.values()) :
            self.gdet_data.create_dataset(key, data = value, compression=self.compressionHDFT, compression_opts=self.compressionHDFL)
            
        if hitsMON is not None:
            
            self.hitsMON = hitsMON
            
            mondic = self.hitsMON.__dict__        
            for key, value in  zip(mondic.keys(),mondic.values()) :
                self.gmon_data.create_dataset(key, data = value, compression=self.compressionHDFT, compression_opts=self.compressionHDFL)
       
       
        # for key, value in {
        #             'Cassette': self.events.Cassette,
        #             'CassetteIDs': self.events.CassetteIDs,
        #             'Duration': self.events.Duration,
        #             'Durations': self.events.Durations,
        #             'Nevents': self.events.Nevents,
        #             'NeventsNotRej2D': self.events.NeventsNotRej2D,
        #             'NeventsNotRejAfterTh': self.events.NeventsNotRejAfterTh,
        #             'NeventsNotRejAll': self.events.NeventsNotRejAll,
        #             'PHS': self.events.PHS,
        #             'PHW': self.events.PHW,
        #             'PrevPT': self.events.PrevPT,
        #             'PulseT': self.events.PulseT,
        #             'ToF': self.events.ToF,
        #             'multS': self.events.multS,
        #             'multW': self.events.multW,
        #             'positionS': self.events.positionS,
        #             'positionSmm': self.events.positionSmm,
        #             'positionW': self.events.positionW,
        #             'positionWmm': self.events.positionWmm,
        #             'positionZmm': self.events.positionZmm,
        #             'timeStamp': self.events.timeStamp,
        #             'wavelength': self.events.wavelength,
                    
        #             }.items():
        #      self.gdet_data.create_dataset(key, data = value, compression=self.compressionHDFT, compression_opts=self.compressionHDFL)
       

 
    
        self.__del__()
        
         
###############################################################################
###############################################################################       
        
class readReducedDataFromHDF():
    def __init__(self, pathAndFileName):
        
        self.parameters = para.parameters()
        self.parameters.init_empty()
        
        self.events = clu.events()
        
        self.hitsMON = maps.hits()
    
        self.pathAndFileName = pathAndFileName
        
        temp2 = os.path.split(self.pathAndFileName)
        self.filePath = temp2[0]+'/'
        self.fileName = temp2[1]
        
        if os.path.exists(self.pathAndFileName) is False:
              print(' \033[1;31m---> File: '+self.fileName+' DOES NOT EXIST \033[1;37m')
              print(' ---> in folder: '+self.filePath)
              # print('  ---> in folder: '+self.filePath+' -> ... it will be skipped!')
               # print(' ---> Exiting ... \n')
               # print('------------------------------------------------------------- \n')
               # sys.exit()
        else:
           
            try:
            
                self.fid    = h5py.File(self.pathAndFileName, "r")
                
            except:
                 
                 print('Unable to opne h5 file')
                 sys.exit()
                 
    def __del__(self):
        try:
            self.fid.close()
        except:
            pass
                 
    def dprint(self, tabs, msg):
        if self.showTree:
            print(tabs+"{}".format(msg))              
      
      
    def read(self, showTree = False):  
        
        self.showTree = showTree
        
        for key_main in self.fid.keys():
            self.dprint('1-',key_main) 
           
            # entry1
            main = self.fid[key_main] 
           
            for key_gr in main.keys():
                self.dprint('2-\t',key_gr)
               
                # detector, monitor, parameters
                group = main[key_gr]
               
                for key_subGr in group.keys():
                    self.dprint('3-\t\t\t',key_subGr)
                    
                    # events
                    subGr = group[key_subGr]
                  
                    if key_gr == 'detector' and key_subGr == 'events':   
                        for keyd in subGr.keys(): 
                            self.dprint('4-\t\t\t\t\t',keyd)
                               
                            self.events.__dict__[keyd] = subGr[keyd][()]
                            
                    if key_gr == 'monitor' and key_subGr == 'hits':   
                        for keyd in subGr.keys(): 
                            self.dprint('4-\t\t\t\t\t',keyd)
                               
                            # here add save into mon hits or events
                            self.hitsMON.__dict__[keyd] = subGr[keyd][()]
                            # print('save mon events still to be implemented')
                            
                    if key_gr == 'parameters':   
                        
                        for keyd in subGr.keys(): 
                            self.dprint('4-\t\t\t\t\t',keyd)
                        
                            self.parameters.__dict__[key_subGr].__dict__[keyd] = subGr[keyd][()]

                        for att, value in zip(subGr.attrs,subGr.attrs.values()):
                            self.dprint('4-\t\t\t\t\t',(att+' -> '+value))  
                            
                            self.parameters.__dict__[key_subGr].__dict__[att] = value

    
        self.__del__()
        
        
  
###############################################################################
###############################################################################

if __name__ == '__main__':

    configFilePath  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/config/'+"MB300_FREIA_config.json"
    filePathD       = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/data/'+"freiatest.pcapng"
   
    tProfilingStart = time.time()
    parameters  = para.parameters('/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/')
    config = maps.read_json_config(configFilePath)
   
    parameters.loadConfigParameters(config)
    
        ### ON/OFF if you want to rsync the data     
    parameters.fileManagement.sync = False  
    
    ### from ... to  ... rsync the data
    parameters.fileManagement.sourcePath = 'essdaq@172.30.244.233:~/pcaps/'
    parameters.fileManagement.destPath   = '/Users/francescopiscitelli/Desktop/dataPcapUtgard/'
    
    # ###############
    
   
    
    ### folder and file to open (file can be a list of files)
    # parameters.fileManagement.filePath = parameters.fileManagement.destPath
    # parameters.fileManagement.fileName = ['freia_1k_pkts_ng.pcapng']
    parameters.fileManagement.fileName = ['freiatest.pcapng']
    
  

    
    # pcap = pcapr.pcapng_reader(filePathD, parameters.clockTicks.NSperClockTick, sortByTimeStampsONOFF = True)
    # readouts = pcap.readouts
    # md  = maps.mapDetector(readouts, config)
    # md.mappAllCassAndChannelsGlob()
    # hits = md.hits
    # cc = clu.clusterHits(hits,parameters.plotting.showStat)
    # cc.clusterizeManyCassettes(parameters.cassettes.cassettes, parameters.dataReduction.timeWindow)
    # events2 = cc.events


    # sav = saveReducedDataToHDF(parameters,'/Users/francescopiscitelli/Desktop/reducedFile/','temp')
    
    # # Nevents = 10
    
    # # dd = sdat.sampleEventsMultipleCassettes([1,2],'/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/data/')
    # # dd.generateGlob(Nevents)
    # # events  = dd.events
    # # eventsArray = events.concatenateEventsInArrayForDebug()


    # # events = 1
    # sav.save(events2)
    
    
    rre = readReducedDataFromHDF('/Users/francescopiscitelli/Desktop/reducedFile/temp.h5')
    rre.read(True)
    eventsOUT  = rre.events
    hitsMONOUT = rre.hitsMON
    paramOUT   = rre.parameters
    
    
    # fid    = h5py.File('/Users/francescopiscitelli/Desktop/reducedFile/temp2.h5', "r")
    
    # events = clu.events()
    # param  = para.parameters()
    # param.init_empty()
    
    # for key1 in fid.keys():
    #         print(key1)
           
    #         ff = fid[key1]
           
    #         for key2 in ff.keys():
               
    #             if key2 == 'parameters':
                   
    #                 print('\t',key2)
               
    #                 fff = ff[key2]
                   
    #                 for key3 in fff.keys():
                        
    #                     if key3 == 'configJsonFile':
                       
    #                         print('\t\t\t',key3)
                       
    #                     # if key3 == 'MONitor' and key2 == 'parameters':
                           
    #                         # self.events.
                           
    #                         ffff = fff[key3]
                            
    #                         # print(ffff.attrs())
                            
    #                         # for a in ffff.attrs:
    #                         #         print(a)
                            
    #                         # param.__dict__[key3] = None
                           
    #                         for key4 in ffff.keys():
                               
    #                             # if  key4 == 'Cassette':
                               
    #                                 print('\t\t\t\t\t',key4)
                               
    #                                 # events.__dict__[key4] = ffff[key4][()]
                                    
    #                                 param.__dict__[key3].__dict__[key4] = ffff[key4][()]
                                    
                                    
                             
    #                         # for key4 in ffff.keys():
    #                         for att, value in zip(ffff.attrs,ffff.attrs.values()):
    #                                 print(' -> ',att,' -> ',value)   
                                    
    #                                 param.__dict__[key3].__dict__[att] = value
                                    
                                
                              