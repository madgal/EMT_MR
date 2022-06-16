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

def getPval(filen,uval):

    tmp=filen.split("/")[-1].split("_")
    for i in range(len(tmp)):
	tmp_name = tmp[i].upper()
	if 'UH' in tmp_name:
		## uh use P(li,ymi)
		if len(tmp[i+1])==4:
			l1=tmp[i+1][0]
			l2=tmp[i+1][1]
			ym1=tmp[i+1][2]
			ym2=tmp[i+1][3]
		else:
			l1=tmp[i+1]
			l2=tmp[i+2]
			ym1=tmp[i+3]
			ym2=tmp[i+4]
		if ("MR1" in tmp) or ("MR2" in tmp) or ("MR3" in tmp) or ("MR4" in tmp):
			if "MR1" in tmp:
				yustr='[0,0,0]'
			elif "MR2" in tmp:
				yustr='[0,0.01,0.09]'
			elif "MR3" in tmp:
				yustr='[0,0.1,0.2]'
			else:
				yustr='[0,0.5,1.]'
			scaleL = 0.2
			if float(ym1)<3:
				scaleym1= 0.001
				scaleym2= 0.002
				addym1=0.001
				addym2=0.001
			else:
				scaleym1= 0.003
				scaleym2= 0.005
				addym1=0.001
				addym2=0.005
		elif ("MR5" in tmp) or ("MR6" in tmp):
			if "MR5" in tmp:
				yustr='[0,0.01,0.09]'
			else:
				yustr='[0,0.1,0.2]'
			scaleL = 0.1
			scaleym1= 0.002
			scaleym2= 0.004
			addym1=0.001
			addym2=0.001
		elif "MR7" in tmp:
			yustr='[0,0.005,0.05]'
			scaleL = 0.2
			scaleym1= 0.01
			scaleym2= 0.2
			addym1=0.
			addym2=0.
		elif "MR8" in tmp:
			yustr='[0,0.005,0.05]'
			scaleL = 0.2
			if float(ym1)<3:
				scaleym1= 0.001
				scaleym2= 0.002
				addym1=0.006
				addym2=0.051
			else:
				scaleym1= 0.01
				scaleym2= 0.1
				addym1=0.006
				addym2=0.051
		elif "MR9" in tmp:
			yustr='[0,0.001,0.009]'
			scaleL = 0.2
			if float(ym1)<3:
				scaleym1= 0.001
				scaleym2= 0.002
				addym1=0.001
				addym2=0.001
			if float(ym1)<6:
				scaleym1= 0.003
				scaleym2= 0.005
				addym1=0.001
				addym2=0.005
			else:
				scaleym1= 0.01
				scaleym2= 0.2
				addym1=0.
				addym2=0.
			
		li=[1.,float(l1)*scaleL,float(l2)*scaleL]
		ymi=[0.,float(ym1)*scaleym1+addym1,float(ym2)*scaleym2+addym2]
		return P(li,ymi,uval),l1,l2,ym1,ym2
def defineMSE(filen,save=False):
    
    df =pd.read_csv(filen).dropna()
    if len(df)<1000:
        return
    else:


        #res_c=getMSE(df)
        #res_c=getMSE2(df)
	#total=0
	#for k in res_c:
	#	total+=res_c[k]
	res = getPval(filen,np.mean(df['u']))

        states=getStates(df)
	count=0
	for key in ['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O'] :
		if states[key]>0:
			count+=1

	fileX = filen.split("/")[-1]
	print"%s,%s,%s,%s,%s,%s,%s" %(fileX,count,res[0],res[1],res[2],res[3],res[4])

print "filename,states,P,li1,li2,ym1,ym2"
for dirn in os.listdir("../coupledWReg_Ccode/crosstalk_uh"):
    if ("EMT_MR" in dirn) :
        for filen in os.listdir("../coupledWReg_Ccode/crosstalk_uh/"+dirn):
            if ("EMT_MR9" in filen) and ("res" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		#print filen	
                filetitle='../coupledWReg_Ccode/crosstalk_uh/'+dirn+"/"+filen
                #plotValuesEMTvMR_simple(filetitle,overlay=True,save=True)
                defineMSE(filetitle,save=True)
