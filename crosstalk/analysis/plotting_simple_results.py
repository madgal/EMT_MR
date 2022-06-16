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

####plotValuesEMTvMR_simple(filen,overlay=True,save=False,numLinks=1)

for dirn in os.listdir("../coupledWReg_Ccode/"):
    if "crosstalk_uh"==dirn:
        for filen in os.listdir("../coupledWReg_Ccode/"+dirn):
            if  ("EMT_MR7" in filen) and ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle='../coupledWReg_Ccode/'+dirn+"/"+filen
                #plotValuesEMTvMR_simple(filetitle,overlay=True,save=True)
                plotValuesEMTvMR_full(filetitle,save=True)
'''


for dirn in os.listdir("../sensitivity_analysis/"):
    if "params" in dirn:
        for dir2 in ['thresholds','cooperativity']:
            for filen in os.listdir("../sensitivity_analysis/"+dirn+"/"+dir2):
    	        if "res.txt" in filen:
    	    		print filen	
                        filetitle='../sensitivity_analysis/'+dirn+"/"+dir2+'/'+filen
                        #plotValuesEMTvMR_simple(filetitle,overlay=True,save=True)
                        plotValuesEMTvMR_full(filetitle,save=True)
'''
