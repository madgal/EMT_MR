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
'''
for filen in os.listdir("../coupledWReg_Ccode/crosstalk_all"):
            if ("EMT_MR10_HS_1_1_Au_1_1_AZ_0_95_AS_0_95_input_0_u3m_0_2_u3n_0_1_hu_0_1_uh_3" in filen) and  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle='../coupledWReg_Ccode/crosstalk_all/'+filen
                plotValuesEMTvMR_full(filetitle,save=True)
for dirn in os.listdir("../coupledWReg_Ccode/singles"):
    if "crosstalk_u3m"==dirn or "crosstalk_uh"==dirn:
	print dirn
        for filen in os.listdir("../coupledWReg_Ccode/singles/"+dirn):
            if  False and ("res.txt" in filen) and ('uh'not in dirn):# and ('uH' not in filen):
                filetitle='../coupledWReg_Ccode/singles/'+dirn+"/"+filen
                plotValuesEMTvMR_full(filetitle,save=True)
	    elif "EMT_MR10" in filen:
        	for filenA in os.listdir("../coupledWReg_Ccode/singles/"+dirn+"/"+filen):
            		if  ("res.txt" in filenA):
               			filetitle='../coupledWReg_Ccode/singles/'+dirn+"/"+filen+"/"+filenA
 		               	plotValuesEMTvMR_full(filetitle,save=True)
for dirn in os.listdir("../coupledWReg_Ccode/doubles"):
    if "crosstalk_AS_Hu" in dirn:
        for filen in os.listdir("../coupledWReg_Ccode/doubles/"+dirn):
            if  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle='../coupledWReg_Ccode/doubles/'+dirn+"/"+filen
                plotValuesEMTvMR_full(filetitle,save=True)
'''
for dirn in os.listdir("../coupledWReg_Ccode/triples"):
    if "crosstalk_u3n_iHu_input" in dirn:
        for filen in os.listdir("../coupledWReg_Ccode/triples/"+dirn):
            if  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle='../coupledWReg_Ccode/triples/'+dirn+"/"+filen
                plotValuesEMTvMR_full(filetitle,save=True)
'''
for dirn in os.listdir("../coupledWReg_Ccode/quads"):
    if "crosstalk" in dirn:
        for filen in os.listdir("../coupledWReg_Ccode/quads/"+dirn):
            if  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle='../coupledWReg_Ccode/quads/'+dirn+"/"+filen
                plotValuesEMTvMR_full(filetitle,save=True)
for dirn in os.listdir("../coupledWReg_Ccode/cinqs"):
    if "crosstalk" in dirn:
        for filen in os.listdir("../coupledWReg_Ccode/cinqs/"+dirn):
            if  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle='../coupledWReg_Ccode/cinqs/'+dirn+"/"+filen
                plotValuesEMTvMR_full(filetitle,save=True)



for dirn in ['bothModified','noPartialEMT','normal_metabolic']:
    for dir2 in os.listdir("../normalCells_coupled/"+dirn):
        if "crosstalk" in dir2:
            for filen in os.listdir("../normalCells_coupled/"+dirn+"/"+dir2):
                if  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
    		    print filen	
                    filetitle='../normalCells_coupled/'+dirn+'/'+dir2+'/'+filen
                    plotValuesEMTvMR_full(filetitle,save=True)
'''
