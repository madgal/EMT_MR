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
from matplotlib import cm

from aux_func_States import *

__clist=[]
cmap = cm.get_cmap('viridis')
for i in range(200):
	__clist+=[mpl.colors.to_hex(cmap(i))]

__clist=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']

def defineMSE(filen,save=False):
    
    df =pd.read_csv(filen).dropna()
    if len(df)<1000:
        return
    else:


        #res_c=getMSE(df)
        res_c=getMSE2(df)
	total=0
	for k in res_c:
		total+=res_c[k]

        states=getStates(df)
	count=0
	for key in ['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O'] :
		if states[key]>0:
			count+=1

	fileX = filen.split("/")[-1]
	print"%s,%s,%s" %(fileX,count,total)

for dirn in os.listdir("../coupledWReg_Ccode/"):
    if ("crosstalk_uh"==dirn) :
        for filen in os.listdir("../coupledWReg_Ccode/"+dirn):
            if ("EMT_MR" in filen) and ("res" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		#print filen	
                filetitle='../coupledWReg_Ccode/'+dirn+"/"+filen
                #plotValuesEMTvMR_simple(filetitle,overlay=True,save=True)
                defineMSE(filetitle,save=True)
