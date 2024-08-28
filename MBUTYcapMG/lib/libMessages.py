#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 11:13:11 2021

@author: francescopiscitelli
"""

import sys

###############################################################################
###############################################################################

class messages():
    
    def __init__(self):
        
        self.red    = '\033[1;31m'
        self.yellow = '\033[1;33m'
        self.white  = '\033[1;37m'
        self.cyan   = '\033[1;36m'
        self.green  = '\033[1;32m'
        
            # print("{}".format(msg))
            
    def warning(self,msg):     
    
        print(self.yellow,end='')
        print('WARNING: {}'.format(msg))
        print(self.white,end='')
            
    def error(self, msg, exitFlag = False):     
    
        print(self.red,end='')
        print('ERROR: {}'.format(msg))
        print(self.white,end='')
        
        if exitFlag:
            print(' ---> Exiting ... \n')
            print('------------------------------------------------------------- \n')
            sys.exit()
        
    def highlight(self,msg):     
    
        print(self.cyan,end='')
        print('{}'.format(msg))
        print(self.white,end='')
            
            
###############################################################################
###############################################################################

if __name__ == '__main__' :
    
    messages().warning('ciao')
    
    messages().error('ciao \n cdds \n exit',False)
    
    messages().highlight('ciao')
    
