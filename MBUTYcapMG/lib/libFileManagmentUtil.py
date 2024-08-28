#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 13:15:57 2021

@author: francescopiscitelli
"""

import numpy as np

import os
import glob
import sys
from PyQt5.QtWidgets import QFileDialog

# from lib import libMapping as maps
# from lib import libParameters as para

# import libMappingMG as maps
# import libParameters as para

###############################################################################
###############################################################################

class fileDialogue():
    def __init__(self, parameters):
        
        self.openMode     = parameters.fileManagement.openMode
        self.filePath     = parameters.fileManagement.filePath
        self.fileSerials  = parameters.fileManagement.fileSerials
        self.fileName     = parameters.fileManagement.fileName
        
        if isinstance(self.fileName, list) is False:
            self.fileName = [self.fileName]
            
        for kk, ff  in enumerate(self.fileName):
             # add extension if not in fileName
             if len(ff.split('.',1)) == 1:
                 self.fileName[kk] = self.fileName[kk]+'.pcapng'

        
    def openFile(self):
        
        if self.openMode == 'window' :
           self.openWindow()
           
        elif self.openMode == 'fileName' :
           self.openFileName()
           
        elif self.openMode == 'latest' :
           self.openLatestFile()
            
        elif self.openMode == 'secondLast' :
           self.openSecondLastFile() 
            
        elif self.openMode == 'wholeFolder' :
           self.openWholeFolder() 
           
        elif self.openMode == 'sequence' :
           self.openSequence() 
           
        elif self.openMode == None :
            self.doNothing()
           
        else:
            print('\n \033[1;31mPlease select a correct open file mode: window, fileName, latest, secondLast, wholeFolder or sequence \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
            
        if self.openMode == 'wholeFolder' :
            print('\033[1;36mAll Files selected in folder '+self.filePath+' \033[1;37m') 
        elif self.openMode == 'sequence' :
            print('\033[1;36mFile selected: ',self.fileName[0],' with serials: ',end='')
            for ss in self.fileSerials:
                print(ss,', ',end='')
            print('\033[1;37m') 
        elif self.openMode == None :
            pass
        else:
            print('\033[1;36mFile selected: ',end='')
            for ff in self.fileName:
                print(ff,', ',end='')
            print('\033[1;37m') 
            
        if  self.openMode != None :
            self.checkIfExists() 
                     
    def checkIfExists(self):
        
        fileNameThatExist = []
            
        # check if file exists in folder
        for fn in self.fileName:
            
           if os.path.exists(self.filePath+fn) is False:
              print(' \033[1;31m---> File: '+fn+' DOES NOT EXIST \033[1;37m')
              print('  ---> in folder: '+self.filePath+' -> ... it will be skipped!')
              print('NOTE: file name must contain extension, e.g. *.pcapng')
              # print(' ---> Exiting ... \n')
              # print('------------------------------------------------------------- \n')
              # sys.exit()
              
           else:
                
                fileNameThatExist.append(fn)
                  
              
        self.fileName = fileNameThatExist
        
        if  len(self.fileName) == 0:
            print(' \033[1;31m---> NO File EXIST in given list \033[1;37m')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
            
        #####################################
            
    ###############################################

   
    def openWindow(self):
        
        temp1 = QFileDialog.getOpenFileName(None, "Select Files", self.filePath, "pcapng files (*.pcapng *pcap)")
        temp2 = os.path.split(temp1[0])
        self.filePath = temp2[0]+'/'
        self.fileName = [temp2[1]]
        if self.fileName == ['']:
            print('\n \033[1;31mNothing selected! \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
            
    def openFileName(self):
        
        pass
    
    def doNothing(self):
        
        pass 
        # print('nothing to open')
        
    def openLatestFile(self):
        
        listOfFiles = glob.glob(self.filePath+'/*.pcapng') 
        if not len(listOfFiles):
            print('\n \033[1;31mNo file exists in directory \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
            
        latestFile  = max(listOfFiles, key=os.path.getmtime)
        temp        = os.path.split(latestFile)
        self.filePath       = temp[0]+'/'
        self.fileName       = [temp[1]]
        
    def openSecondLastFile(self):
        
        listOfFiles = glob.glob(self.filePath+'/*.pcapng') 
        if not len(listOfFiles):
            print('\n \033[1;31mNo file exists in directory \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
            
        listOfFiles.sort(key=os.path.getmtime)
        temp2        = os.path.split(listOfFiles[-2])
        self.filePath       = temp2[0]+'/'
        self.fileName       = [temp2[1]]
        
    def openWholeFolder(self):
        
        # here add the possibility to specifiy a piece of the filename as a tag
        
        listOfFiles = glob.glob(self.filePath+'/*.pcapng') 

        listOfFiles.sort(key=os.path.getmtime)
        
        fileName = []
        
        for fp in listOfFiles:
  
            temp = fp.rsplit('/',1)[1]
            
            fileName.append(temp)
            
        self.fileName = fileName  
        
    def openSequence(self):
        
        # print(self.fileSerials)
        # print(type(self.fileSerials))
        # self.fileName = []
        
        if isinstance(self.fileName, list) is False:
            self.fileName = [self.fileName]
        
        if len(self.fileName) > 1:
            print('\n \033[1;31mFileName must be a single file for sequence open mode \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
   
        # temp2 = os.path.split(self.fileName)
        
        temp = str(self.fileName[0]).rsplit('.'[-1])
        extension='.'+temp[1]
        
        temp2 = temp[0].rsplit('_',1)
        
        if len(temp2) > 1:
        
            fileName1 = temp2[0]+'_'
            
            self.fileName = []
            
            for ser in self.fileSerials:
                
                self.fileName.append(fileName1+str(format(ser,'05d'))+extension)
                
        elif   len(temp2) == 1:  
            
            print('\n \033[1;31mFileName must be of format blabla_XXXXX.pcapng in sequence mode (XXXXX serial) \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
            
        else:
            
            print('\n \033[1;31mUnknown FileName error in sequence mode \033[1;37m\n')
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
    
###############################################################################
###############################################################################

if __name__ == '__main__' :
    
    path = '/Users/francescopiscitelli/Desktop/dataPcapUtgard/'
    
    configFilePath  = '/Users/francescopiscitelli/Documents/PYTHON/MBUTYcap/config/'
    configFileName  = "MB300_AMOR_config.json"
    config = maps.read_json_config(configFilePath+configFileName)
    
    currentPath = os.path.abspath(os.path.dirname(__file__))+'/'
    parameters  = para.parameters(currentPath)

    parameters.loadConfigParameters(config)

    parameters.fileManagement.filePath = path
    
    # parameters.fileManagement.fileName = ['20211007_114629_duration_s_4_temp_00000.pcapng']

    parameters.fileManagement.fileName = '20211007_114629_duration_s_4_temp_00000.pcapng'
  
    parameters.fileManagement.openMode = 'sequence'
    
    parameters.fileManagement.openMode = 'fileName'
    parameters.fileManagement.openMode = 'window'
    parameters.fileManagement.openMode = 'latest'
    
    
    parameters.fileManagement.fileSerials = np.arange(0,5,1)
    
    parameters.fileManagement.fileSerials = [0,5,6]
    
    fileDetails  = fileDialogue(parameters)
    aa = fileDetails.openFile()
    
    print(fileDetails.fileName)
    
    
    