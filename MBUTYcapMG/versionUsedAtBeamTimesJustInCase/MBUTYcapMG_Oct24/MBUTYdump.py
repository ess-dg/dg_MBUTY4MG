#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 14:24:18 2021

@author: francescopiscitelli
"""

###############################################################################
###############################################################################
########    V2.0 2021/12/14      francescopiscitelli     ######################
###############################################################################
###############################################################################

import argparse
# import os
# import sys
# from datetime import datetime
# import subprocess

from lib import libTerminal as ta

############################################################################### 
###############################################################################

class dumpToFile():
    
     def __init__(self, pathToTshark, interface='en0', destPath='./', fileName='temp', typeOfCapture='packets', quantity=100, numOfFiles=1, delay=0):
         
         rec = ta.dumpToPcapngUtil(pathToTshark, interface, destPath, fileName)
         
         sta = ta.acquisitionStatus(destPath)  
         sta.set_RecStatus()
         
         status = rec.dump(typeOfCapture,quantity,numOfFiles,delay,fileNameOnly=False)
         if status == 0: 
              sta.set_FinStatus()
         else:
              sta.set_RecStatus()
      

###############################################################################
###############################################################################        
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
     
    #  ARGS for FRAPIS LAPTOP
    # parser.add_argument("-i", "--interface", help = "interface from which capture packets", type = str, default = "en0")
    # parser.add_argument("-t", "--tshark", help = "path where tshark is located", type = str, default = "/Applications/Wireshark.app/Contents/MacOS/")
    # parser.add_argument("-e", "--destination", help = "path where to save recorded pcapng files", type = str, default = "/Users/francescopiscitelli/Desktop/reducedFile/")
    
    
    # # #  ARGS for ESSDAQ
    # parser.add_argument("-i", "--interface", help = "interface from which capture packets", type = str, default = "p4p1")
    # parser.add_argument("-t", "--tshark", help = "path where tshark is located", type = str, default = "/usr/sbin/")
    # parser.add_argument("-e", "--destination", help = "path where to save recorded pcapng files", type = str, default = "/home/essdaq/pcaps/")

    # # # ARGS for ESSDAQ EFU
    # parser.add_argument("-i", "--interface", help = "interface from which capture packets", type = str, default = "ens2")
    # parser.add_argument("-t", "--tshark", help = "path where tshark is located", type = str, default = "/usr/sbin/")
    # parser.add_argument("-e", "--destination", help = "path where to save recorded pcapng files", type = str, default = "/home/essdaq/pcaps/")

    # #  ARGS for EFU  JADAQ
    # parser.add_argument("-i", "--interface", help = "interface from which capture packets", type = str, default = "em2")
    # parser.add_argument("-t", "--tshark", help = "path where tshark is located", type = str, default = "/usr/sbin/")
    # parser.add_argument("-e", "--destination", help = "path where to save recorded pcapng files", type = str, default = "/home/efu/data/pcaps/")

    # # # ARGS for MG
    parser.add_argument("-i", "--interface", help = "interface from which capture packets", type = str, default = "enp5s0")
    parser.add_argument("-t", "--tshark", help = "path where tshark is located", type = str, default = "/usr/sbin/")
    parser.add_argument("-e", "--destination", help = "path where to save recorded pcapng files", type = str, default = "/home/mg/data/VMM-Utgard-test-2024/")

    # common  fields
    parser.add_argument("-f", "--file", help = "pcapng filename", type = str, default = "temp")
    parser.add_argument("-n", "--numoffiles", help = "num of files to record in sequence", type = int, default = 1)
    parser.add_argument("-r", "--delay", help = "sleep (s) between consecutive files", type = int, default = 0)
    
    # mutually exclusive fields
    command_group = parser.add_mutually_exclusive_group(required=False)
    command_group.add_argument("-d", "--duration", help = "duration (seconds)", type=int)
    command_group.add_argument("-p", "--packets", help = "num of packets (int)", type=int)
    command_group.add_argument("-s", "--size", help = "file size (kbytes)", type=int)

    ###
    args  = parser.parse_args()

    ###
    if args.duration is None and args.packets is None and args.size is None:
        #  default packets
        args.packets = 100

    ###
    if args.duration is not None:
        typeOfCapture = 'duration'
        quantity = args.duration
        
    if args.packets is not None:
        typeOfCapture = 'packets'
        quantity = args.packets
        
    if args.size is not None:
        typeOfCapture = 'filesize'
        quantity = args.size
        
    rec = dumpToFile(pathToTshark=args.tshark, interface=args.interface, destPath=args.destination, fileName=args.file,typeOfCapture=typeOfCapture,quantity=quantity,numOfFiles=args.numoffiles,delay=args.delay)


# print('--------------------')
# print('interface: '+args.interface)
# print('tshark: '+args.tshark)
# print('path: '+args.destination)
# print('file: '+args.file)
# print('num of files: '+str(args.numoffiles))

# print('duration: '+str(args.duration))
# print('packets: '+str(args.packets))
# print('size: '+str(args.size))
    
   