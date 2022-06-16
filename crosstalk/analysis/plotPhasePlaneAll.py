import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import os
import math
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import copy

from aux_func_States import *
from functions_for_plotting import *

mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams.update({'font.size': 20})
title=''
reglist=['AS','AZ','AU','HS','HU','U3M','U3N','INPUT','UH','IHU']


for filen in os.listdir('data'):
	fsave='figures/'+filen.split(".")[0]
	flist = filen.split(".")[0].split("_")
	if flist[0]=='crosstalk':
		sets=flist[1:]
	else:
		sets=flist[2:]
	for i in range(len(sets)):
		if sets[i].upper()=='IHU':
			sets[i]='HU'
		elif 'UH' in sets[i].upper():
			sets[i]='UH'
		else:
			sets[i]=sets[i].upper()

	if len(sets)==1:
           if ('singles' not in filen) and ('vals' not in filen) and ('no' in filen):# and ('PSF' in filen)
		if ('UH' in sets[0].upper()) or (sets[0].upper() in reglist):
			print filen
			title=sets[0]
			xlabel=sets[0]
			try:
        			plotCoupledPhenotypes_singleLink('data/'+filen,xlabel,title,fsave)
    				plotICSvlambda_singleLink('data/'+filen,xlabel,title,fsave)
				print "success with ",filen
			except:
				print "Issue with", filen
		elif sets[0]=='all':
			print 'a1',filen,sets
		## else ignore as its requires special treatment
	elif len(sets)==2:
           if ('no' in filen):# and  ('PSF'  in filen):
		xlabel=sets[0]
		ylabel=sets[1]
		title = xlabel+' '+ylabel
		try:
			plotCoupledPhenotypes_doubleLink('data/'+filen,xlabel,ylabel,title,fsave)
			print "success with ",filen
		except:
			print "Issue with", filen
	elif False:#len(sets)>0:
		xregList=[]
		if 'UH' in sets:
			xregList=['UH']
		i=0
		while len(xregList)<np.floor(len(sets)/2.):
			xregList+=[sets[i].upper()]
			i+=1
		title = ' '
		try:
		   	plotCoupledPhenotypes_regulatory('data/'+filen,title,fsave,xregList=xregList)
			print "success with ",filen
		except:
			print "Issue with", filen
