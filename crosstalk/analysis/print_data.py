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

from aux_func_States import *

from functions_for_plotting import *
 
'''
for i in [10]:#1,2,3,4,5,6,7,9,10]:
    dirN='../coupledWReg_Ccode/singles/crosstalk_uh/EMT_MR'+str(i)+'/'
    fsave='data/crosstalk_uh'+str(i)+'.txt'
    #try:
    output_results(dirN,fsave,skipUH=False)
    #except:
    #		print "failure with", dirN
'''
for el in ['triples']:#,'higherSets','doubles','PSF','quads','singles',]:
    for dirN in os.listdir('../coupledWReg_Ccode/'+el):
	if "crosstalk_u3n_iHu_input"==dirN:
	    if True:
		if 'PSF' in el:
			fsave = 'data/PSF_'+dirN+'.txt'
		else:
			fsave = 'data/'+dirN+'.txt'
		dirname='../coupledWReg_Ccode/'+el+'/'+dirN+'/'
		if 'crosstalk_all'==dirN:
			output_results(dirname,fsave,skipUH=True)
		else:
			output_results(dirname,fsave)
		if 'uh' in dirN.lower():
			output_results(dirname,fsave,skipUH=True)
			fsave = fsave.replace("data/","data/IND_")
			output_results(dirname,fsave,skipUH=False)

		print 'Succeeded',dirN
	    else:#xcept:
		print "failure with", dirN
'''
for dirN in ['comparison_coupled','crosstalk_comparison']:
	    try:
		fsave = 'data/'+dirN+'.txt'
		dirname='../coupledWReg_Ccode/'+dirN+'/'
		output_results(dirname,fsave)
		print 'Succeeded',dirN
	    except:
		print "failure with", dirN

for el in ['noPartialEMT','bothModified','normal_metabolic']:
    abv={'bothModified':'noHH','noPartialEMT':'noEM','normal_metabolic':'noWO'}
    for dirN in os.listdir('../normalCells_coupled/'+el):
	if True:
	    try:
		fsave = 'data/'+abv[el]+"_"+dirN+'.txt'
		dirname='../normalCells_coupled/'+el+'/'+dirN+'/'
		output_results(dirname,fsave)
		print 'Succeeded',dirN
	    except:
                    print "failure with", dirN
'''
