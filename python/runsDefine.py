#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:27:20 2020

@author: mccikpc2
"""
set1=2
if set1==1:
    runToDo = [['tinit=280.','tinit=290.'], \
               ['n_aer1(1:3,1:1)        = 2.27e9, 3.08e9, 1.05e9, d_aer1(1:3,1:1)        = 13e-9   , 32e-9, 123e-9, sig_aer1(1:3,1:1)      = 0.66   , 0.69, 0.54, ', 
            'n_aer1(1:3,1:1)        = 1.0e7, 5.08e7, 0.55e8,    d_aer1(1:3,1:1)        = 13e-9   , 32e-9, 223e-9,    sig_aer1(1:3,1:1)      = 0.66   , 0.69, 0.54, '], \
               ['hm_flag=.false.','hm_flag=.true.'], \
               ['break_flag=0','break_flag=2'], \
               ['mode1_flag=.false.','mode1_flag=.true.'], \
               ['mode2_flag=.false.','mode2_flag=.true.'], \
                   ]
    columns1=['t','aer','hm','br','m1','m2']
elif set1==2:
    runToDo = [[['tinit=276.869','pinit=57645.1,','rhinit=0.52341,', \
                'n_aer1(1:4,1:1)        = 850e6, 8e6, 210e6, 0.e6, d_aer1(1:4,1:1)        = 0.099e-6  , 0.2e-6, 0.04e-6, 1e-7, sig_aer1(1:4,1:1)      = 0.39   , 0.05, 0.25, 0.25,'],\
                ['tinit=282.069','pinit=65016.9,','rhinit=0.66744,',\
                'n_aer1(1:4,1:1)        = 850e6, 7e6, 650e6, 20.e6, d_aer1(1:4,1:1)        = 0.115e-6  , 0.195e-6, 0.05e-6, 0.26e-6, sig_aer1(1:4,1:1)      = 0.37   , 0.03, 0.3, 0.08,'] ], \
               [['winit=2.5','t_thresh=1600.'], ['winit=5.0','t_thresh=800.']], \
               ['hm_flag=.false.','hm_flag=.true.'], \
               ['break_flag=0','break_flag=2'], \
               ['mode1_flag=.false.','mode1_flag=.true.'], \
               ['mode2_flag=.false.','mode2_flag=.true.'], \
                   ]
    columns1=['case','w','hm','br','m1','m2']
    

outputDir='/tmp'