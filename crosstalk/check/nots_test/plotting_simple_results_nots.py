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


from aux_func_States_nots import *

def plotValuesEMTvMR_full(filen,save=False):
    
    df =pd.read_csv(filen).dropna()
    if len(df.values[:])<1000:
        return
    else:

	fig = plt.figure(figsize=(10,12))
	gs1 = gridspec.GridSpec(3,4,width_ratios=[1,0.1,1,0.2],height_ratios=[1.,0.1,1.])
	#gs1.update(left=0.05,right=8,top=10,bottom=0)
	plt.rcParams['font.size']=20

	ax1 = plt.subplot(gs1[0,0])
	ax2 = plt.subplot(gs1[0,2])
	ax3 = plt.subplot(gs1[2,:3])

        labels={"u":"$\mu_{200}$ (molecules)",'mz':'Zeb mNRA (molecules)','A':'AMPK (nmol/L)','h':'Hif-1 (nmol/L)'}
    ##columns = u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh
        x1L = 0#'u'
        x2L = 6#'A'
        y1L=9#'h'
        y2L=1#'mz'
        x1m=np.max([30000,np.max(df['u'].values)])
        x2m=np.max([1000,np.max(df['A'].values)])
        y1m=np.max([1000,np.max(df['h'].values)])
        y2m=np.max([1200,np.max(df['mz'].values)])
        x1v=np.arange(0,x1m)
        x2v=np.arange(0,x2m)
        y1v=np.arange(0,y1m)
        y2v=np.arange(0,y2m)
        color_list=['k','r','g','b','y','pink','orange','purple','indigo']
	legend_elements=[]

	uncoupled =pd.read_csv("../coupledWReg_Ccode/"+filen).dropna()

	full_setp = np.append(uncoupled.values,df.values,axis=0)
        mean_list = np.mean(full_setp,axis=0)
        std_list = np.std(full_setp,axis=0)
        uncoupledZ = (uncoupled.values-mean_list)/std_list
        uncoupled_fp = np.unique(np.round(uncoupledZ,6),axis=0)
        stateLabels = returnStateLabels(uncoupled_fp)

	for i in range(len(uncoupled_fp)):
		ax1.plot(uncoupled_fp[i,x1L],uncoupled_fp[i,y1L],'o',markersize=10,label=stateLabels[i],color=color_list[i])
		ax2.plot(uncoupled_fp[i,x2L],uncoupled_fp[i,y2L],'o',markersize=10,label=stateLabels[i],color=color_list[i])

	for i in range(len(uncoupled_fp)):
		legend_elements+=[ Patch(facecolor=color_list[i], edgecolor=color_list[i], label=stateLabels[i])]

	tmp=[]
	for i in range(len(df.values[:])):
		tmp+=[df.values[i]]


	reduced = np.unique(np.round(tmp,0),axis=0)
        for i in range(len(reduced)):
		resZ = (reduced[i]-mean_list)/std_list
		
		#distances=0
		#for j in [0,1,2,4,5,6,7,8,9,10]:
		#	distances+=((resZ[j]-uncoupled_fp[:,j])**2)

		##distances = np.sum((df_res.values[i]-uncoupled_fp)**2,axis=1)
		distances = np.sum((resZ-uncoupled_fp)**2,axis=1)
                location= np.argwhere(np.min(distances)==distances)[:,0][0]
		ax1.plot(resZ[x1L],resZ[y1L],'X',markersize=22,color=color_list[location])
		ax2.plot(resZ[x2L],resZ[y2L],'X',markersize=22,color=color_list[location])

    	ax2.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')

	ax1.set_xlabel(labels['u'],fontsize=20)
	ax2.set_xlabel(labels['A'],fontsize=20)
	ax1.set_ylabel(labels['h'],fontsize=20)
	ax2.set_ylabel(labels['mz'],fontsize=20)


        res_uc=getStates(uncoupled)
        res_c=getStates(df)

	count=0
	xlabs=[]
	
	for key in ['E','EM','M','W','WO','O','E/W','E/WO','E/O','EM/W','EM/WO','EM/O','M/W','M/WO','M/O']:#res_uc:
		if count==0 and res_c[key]!=0:
			ax3.bar(count+0.2,np.log10(res_c[key]/10.),color='k',width=0.2,label='coupled')
			ax3.bar(count,np.log10(res_uc[key]/10.),color='r',width=0.2,label='uncoupled')
		elif res_c[key]!=0:
			ax3.bar(count+0.2,np.log10(res_c[key]/10.),color='k',width=0.2)
			ax3.bar(count,np.log10(res_uc[key]/10.),color='r',width=0.2)
		else:
			ax3.bar(count,np.log10(res_uc[key]/10.),color='r',width=0.2)
		xlabs+=[key]
		count+=1

	ax3.set_xticks(np.arange(0,count))
	ax3.set_xticklabels(xlabs,rotation=50)
    	ax3.legend( bbox_to_anchor=(1.05, 1), loc='upper left')
	ax3.set_ylabel("Initial conditions\nleading to SS (%)")
	ax3.set_xlabel("State label")

        if save:
		title = filen.replace("res.txt",".png")
                fig.savefig(title,bbox_inches='tight')
       	else:
              	plt.show()
	plt.close()
###########################################

for dirn in os.listdir("."):
    if "crosstalk" in dirn:
        for filen in os.listdir(dirn):
            if  ("res.txt" in filen):# and ('uh' not in filen) and ('uH' not in filen):
		print filen	
                filetitle=dirn+"/"+filen
		try:
                	plotValuesEMTvMR_full(filetitle,save=True)
		except:
			print "Issue " ,filetitle
