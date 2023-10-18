#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:27:20 2020

@author: mccikpc2
"""

runToDo = [['tinit=280.','tinit=290.'], \
           ['    n_aer1(1:3,1:1)        = 2.27e9, 3.08e9, 1.05e9,' + \
'    d_aer1(1:3,1:1)        = 13e-9   , 32e-9, 123e-9, ' + \
'    sig_aer1(1:3,1:1)      = 0.66   , 0.69, 0.54, ', \
'    n_aer1(1:3,1:1)        = 1.0e7, 5.08e7, 0.55e8,' + \
'    d_aer1(1:3,1:1)        = 13e-9   , 32e-9, 223e-9, ' + \
'    sig_aer1(1:3,1:1)      = 0.66   , 0.69, 0.54, '], \
           ['hm_flag=.false.','hm_flag=.true.'], \
           ['break_flag=0','break_flag=2'], \
           ['mode1_flag=.false.','mode1_flag=.true.'], \
           ['mode2_flag=.false.','mode2_flag=.true.'], \
               ]
columns1=['t','aer','hm','br','m1','m2']
    

outputDir='/tmp'