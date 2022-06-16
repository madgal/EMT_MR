import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import os
import math
import time
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import ticker


from aux_func_States import *


hlw = 3

__clist=[]
cmap1= cm.get_cmap('tab20b')
cmap2= cm.get_cmap('tab20c')

count,base=0,0
for i in range(20):
    __clist+=[mpl.colors.to_hex(cmap1(base+count*4))]
    __clist+=[mpl.colors.to_hex(cmap2(base+count*4))]
    count+=1
    if count==5:
            count=0
            base+=1

#__clist=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']



##########
##########
##########
####def plotValuesEMTvMR_full(filen,save=False,compI=False):
####
####def plotICSvlambda_singleLink(dirN,xlabel,title,fsave):
####def plotICSvlambda_singleLinkComp(dirN,xlabel,title,fsave):
####def plotICSvlambda_doubleLink(dirN,xlabel,title,fsave,compI=False):
#### plotICSvlambda_RegComp(filen,xaxis,title,fsave,scaleList=[]):
####def    __plotSet1(xdir,ydir,col,xlabel,title,fsave,comp=False,labels=False,xlim=[-1],ylim=[-1],compI=False):
####def    __plotSet2(xdir,ydir,col,xlabel,title,fsave,comp=False,labels=False,compI=False):
####
####def plotCoupledPhenotypes_singleLink(dirN,xlabel,title,fsave,logScale=False,checkRanges=False):
####def plotCoupledPhenotypes_doubleLink(fileN,xaxis,yaxis,title,fsave,xlim=[-1],ylim=[-1],checkRanges=False):
####def plotCoupledPhenotypes_regulatory(fileN,title,fsave,xlog=False,ylog=False,xlim=[-1],ylim=[-1],checkRanges=False):
####def plotCoupledPhenotypes_singleLinkComp(fileN,yaxis,title,fsave,logScale=False,checkRanges=False):
####
####def output_results(dirN,fsave=''):
####
##########
##########
##########

def    __plotSet1(xdir,ydir,xlabel,title,fsave,col='k',comp=False,labels=False,xlim=[-1],ylim=[-1],compI=False):
    ## 9 plot figures for coupled steady states
    fig = plt.figure()
    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(left=0.05, right=2, wspace=0.5,hspace=0.3,top=3,bottom=0)

    axC1 = plt.subplot(gs1[0, 0])
    axC2 = plt.subplot(gs1[0, 1])
    axC3 = plt.subplot(gs1[0, 2])
    axC4 = plt.subplot(gs1[1, 0])
    axC5 = plt.subplot(gs1[1, 1])
    axC6 = plt.subplot(gs1[1, 2])
    axC7 = plt.subplot(gs1[2, 0])
    axC8 = plt.subplot(gs1[2, 1])
    axC9 = plt.subplot(gs1[2, 2])

    if col=='k':
	col={}
	for key in xdir:
		col[key]=[]
		for i in range(len(xdir[key])):
			col[key]+=['k']
	
    for key in xdir:
	xdir[key]= np.array(xdir[key])
	ydir[key]= np.array(ydir[key])
    
    if comp:
	if not compI:
	        df =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    		tempStates=getStates(df)
	else:
		df = pd.read_csv("data/crosstalk_input.txt",header=0)
		tempStates={}
		row = np.argwhere(df['INPUT'].values==compI)[:,0][0]
		for key in ['E','EM','M','W','O','WO','E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']:
			tempStates[key] = df[key][row]
        nics=1000.
        axC1.plot(xdir['E/O'], xdir['E/O']*0.+tempStates['E/O']/nics*100.,'r-')
        axC2.plot(xdir['E/WO'], xdir['E/WO']*0.+tempStates['E/WO']/nics*100.,'r-')
        axC3.plot(xdir['E/W'], xdir['E/W']*0.+tempStates['E/W']/nics*100.,'r-')
        axC4.plot(xdir['EM/O'],xdir['EM/O']*0.+tempStates['EM/O']/nics*100.,'r-')
        axC5.plot(xdir['EM/WO'],xdir['EM/WO']*0.+tempStates['EM/WO']/nics*100.,'r-')
        axC6.plot(xdir['EM/W'],xdir['EM/W']*0.+tempStates['EM/W']/nics*100.,'r-')
        axC7.plot(xdir['M/O'], xdir['M/O']*0.+tempStates['M/O']/nics*100.,'r-')
        axC8.plot(xdir['M/WO'], xdir['M/WO']*0.+tempStates['M/WO']/nics*100.,'r-')
        axC9.plot(xdir['M/W'], xdir['M/W']*0.+tempStates['M/W']/nics*100.,'r-')
        
    usedL=[]
    for i in range(len(xdir['E/O'])):
        if labels and (labels['E'][i] not in usedL):
	    usedL+=[labels['E'][i]]
            axC1.plot(xdir['E/O'][i],ydir['E/O'][i],marker='o',color=col['E/O'][i],label=labels['E/O'][i])
            axC2.plot(xdir['E/WO'][i],ydir['E/WO'][i],marker='o',color=col['E/WO'][i],label=labels['E/WO'][i])
            axC3.plot(xdir['E/W'][i],ydir['E/W'][i],marker='o',color=col['E/W'][i],label=labels['E/W'][i])
            axC4.plot(xdir['EM/O'][i],ydir['EM/O'][i],marker='o',color=col['EM/O'][i],label=labels['EM/O'][i])
            axC5.plot(xdir['EM/WO'][i],ydir['EM/WO'][i],marker='o',color=col['EM/WO'][i],label=labels['EM/WO'][i])
            axC6.plot(xdir['EM/W'][i],ydir['EM/W'][i],marker='o',color=col['EM/W'][i],label=labels['EM/W'][i])
            axC7.plot(xdir['M/O'][i],ydir['M/O'][i],marker='o',color=col['M/O'][i],label=labels['M/O'][i])
            axC8.plot(xdir['M/WO'][i],ydir['M/WO'][i],marker='o',color=col['M/WO'][i],label=labels['M/WO'][i])
            axC9.plot(xdir['M/W'][i],ydir['M/W'][i],marker='o',color=col['M/W'][i],label=labels['M/W'][i])
        else:
            axC1.plot(xdir['E/O'][i],ydir['E/O'][i],marker='o',color=col['E/O'][i])
            axC2.plot(xdir['E/WO'][i],ydir['E/WO'][i],marker='o',color=col['E/WO'][i])
            axC3.plot(xdir['E/W'][i],ydir['E/W'][i],marker='o',color=col['E/W'][i])
            axC4.plot(xdir['EM/O'][i],ydir['EM/O'][i],marker='o',color=col['EM/O'][i])
            axC5.plot(xdir['EM/WO'][i],ydir['EM/WO'][i],marker='o',color=col['EM/WO'][i])
            axC6.plot(xdir['EM/W'][i],ydir['EM/W'][i],marker='o',color=col['EM/W'][i])
            axC7.plot(xdir['M/O'][i],ydir['M/O'][i],marker='o',color=col['M/O'][i])
            axC8.plot(xdir['M/WO'][i],ydir['M/WO'][i],marker='o',color=col['M/WO'][i])
            axC9.plot(xdir['M/W'][i],ydir['M/W'][i],marker='o',color=col['M/W'][i])

    axC1.set_ylabel("E/O states (%)")
    axC2.set_ylabel("E/WO states (%)")
    axC3.set_ylabel("E/W states (%)")
    axC4.set_ylabel("EM/O states (%)")
    axC5.set_ylabel("EM/WO states (%)")
    axC6.set_ylabel("EM/W states (%)")
    axC7.set_ylabel("M/O states (%)")
    axC8.set_ylabel("M/WO states (%)")
    axC9.set_ylabel("M/W states (%)")
    
    axC1.set_xlabel(xlabel)
    axC2.set_xlabel(xlabel)
    axC3.set_xlabel(xlabel)
    axC4.set_xlabel(xlabel)
    axC5.set_xlabel(xlabel)
    axC6.set_xlabel(xlabel)
    axC7.set_xlabel(xlabel)
    axC8.set_xlabel(xlabel)
    axC9.set_xlabel(xlabel)
    
    axC9.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    fig.savefig(fsave+"_LvCS.png",bbox_inches='tight')
    plt.close()
    #plt.show()

###################################
###################################
###################################
#########################
#########################

def    __plotSet2(xdir,ydir,xlabel,title,fsave,col='k',comp=False,labels=False,compI=False):
    ## 6 plot figures for individual steady states
    fig = plt.figure()
    gs2 = gridspec.GridSpec(2, 3)
    gs2.update(left=0.05, right=2, wspace=0.5,hspace=0.3,top=2,bottom=0)
    
    axS1 = plt.subplot(gs2[0, 0])
    axS2 = plt.subplot(gs2[0, 1])
    axS3 = plt.subplot(gs2[0, 2])
    axS4 = plt.subplot(gs2[1, 0])
    axS5 = plt.subplot(gs2[1, 1])
    axS6 = plt.subplot(gs2[1, 2])

    if col=='k':
	col={}
	for key in xdir:
		col[key]=[]
		for i in range(len(xdir[key])):
			col[key]+=['k']
    for key in xdir:
	xdir[key]= np.array(xdir[key])
	ydir[key]= np.array(ydir[key])
	
    if comp:
    	if not compI:
	        df =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    		tempStates=getStates(df)
	else:
		df = pd.read_csv("data/crosstalk_input.txt",header=0)
		tempStates={}
		row = np.argwhere(df['INPUT'].values==compI)[:,0][0]
		for key in ['E','EM','M','W','O','WO','E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']:
			tempStates[key] = df[key][row]
        nics=1000.
        axS1.plot(xdir['E'], xdir['E']*0.+tempStates['E']/nics*100.,'r-')
        axS2.plot(xdir['EM'],xdir['EM']*0.+tempStates['EM']/nics*100.,'r-')
        axS3.plot(xdir['M'], xdir['M']*0.+tempStates['M']/nics*100.,'r-')
        axS4.plot(xdir['O'], xdir['O']*0.+tempStates['O']/nics*100.,'r-')
        axS5.plot(xdir['WO'],xdir['WO']*0.+tempStates['WO']/nics*100.,'r-')
        axS6.plot(xdir['W'], xdir['W']*0.+tempStates['W']/nics*100.,'r-')
    
    usedL = []
    for i in range(len(xdir['E/O'])):
        if labels and (labels['E'][i] not in usedL):
	    usedL+=[labels['E'][i]]
            axS1.plot(xdir['E'][i],ydir['E'][i],marker='o',color=col['E'][i],label=labels['E'][i])
            axS2.plot(xdir['EM'][i],ydir['EM'][i],marker='o',color=col['EM'][i],label=labels['EM'][i])
            axS3.plot(xdir['M'][i],ydir['M'][i],marker='o',color=col['M'][i],label=labels['M'][i])
            axS4.plot(xdir['O'][i],ydir['O'][i],marker='o',color=col['O'][i],label=labels['O'][i])
            axS5.plot(xdir['WO'][i],ydir['WO'][i],marker='o',color=col['WO'][i],label=labels['WO'][i])
            axS6.plot(xdir['W'][i],ydir['W'][i],marker='o',color=col['W'][i],label=labels['W'][i])
        else:
            axS1.plot(xdir['E'][i],ydir['E'][i],marker='o',color=col['E'][i])
            axS2.plot(xdir['EM'][i],ydir['EM'][i],marker='o',color=col['EM'][i])
            axS3.plot(xdir['M'][i],ydir['M'][i],marker='o',color=col['M'][i])
            axS4.plot(xdir['O'][i],ydir['O'][i],marker='o',color=col['O'][i])
            axS5.plot(xdir['WO'][i],ydir['WO'][i],marker='o',color=col['WO'][i])
            axS6.plot(xdir['W'][i],ydir['W'][i],marker='o',color=col['W'][i])
     
    axS1.set_ylabel("E states (%)")
    axS2.set_ylabel("EM states (%)")
    axS3.set_ylabel("M states (%)")
    axS4.set_ylabel("O states (%)")
    axS5.set_ylabel("WO states (%)")
    axS6.set_ylabel("W states (%)")
    
    axS1.set_xlabel(xlabel)
    axS2.set_xlabel(xlabel)
    axS3.set_xlabel(xlabel)
    axS4.set_xlabel(xlabel)
    axS5.set_xlabel(xlabel)
    axS6.set_xlabel(xlabel)
    
    axS6.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    fig.savefig(fsave+"_LvOS.png",bbox_inches='tight')
    plt.close()
    #plt.show()




####################3
###########################
###########################
def plotValuesEMTvMR_full(filen,save=None,compI=None,comp=True,PSF=None):
    
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

	uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
	full_setp = np.append(uncoupled.values,df.values,axis=0)
        mean_list = np.mean(full_setp,axis=0)
        std_list = np.std(full_setp,axis=0)
        uncoupledZ = (uncoupled.values-mean_list)/std_list
        uncoupled_fp = np.unique(np.round(uncoupledZ,6),axis=0)

	for i in range(len(uncoupled_fp)):
		ax1.plot(uncoupled_fp[i,x1L],uncoupled_fp[i,y1L],'o',markersize=10)#,label=stateLabels[i],color=color_list[i])
		ax2.plot(uncoupled_fp[i,x2L],uncoupled_fp[i,y2L],'o',markersize=10),#label=stateLabels[i],color=color_list[i])

	#for i in range(len(uncoupled_fp)):
	#	legend_elements+=[ Patch(facecolor=color_list[i], edgecolor=color_list[i], label=stateLabels[i])]

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
	ax1.set_ylabel(labels['h'],fontsize=20)
	ax2.set_xlabel(labels['A'],fontsize=20)
	#####3ax2.set_ylabel(labels['mz'],fontsize=20)
	ax2.set_ylabel(labels['h'],fontsize=20)


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
        ax1.text(1, 1.1, getTitle(filen,np.mean(df['u'])), horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,fontsize=30)
        if save:
		title = filen.replace("res.txt",".png")
                fig.savefig(title,bbox_inches='tight')
       	else:
              	plt.show()
	plt.close()
###########################################
###################################
###################################
###################################
def plotCoupledPhenotypes_singleLinkComp(fileN,yaxis,title,fsave,scaleList=[],xlim=[],checkRanges=False):

    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    xlabel='Increasing regulation'
    ylabel=[]
    labelL=[]
    hybridVL=[]
    x1L,x2L,y1L,y2L=[],[],[],[]
    for jj in range(len(fileN)):
        hybridV=[]
        xval,yval,label=[],[],[]
        df = pd.read_csv(fileN[jj],header=0)
        nics = float(df['nics'][0])
	if len(scaleList)>0:
		scaleX=scaleList[jj]
	else:
		if 'uh'.upper() in yaxis[jj].upper():
	        	scaleX = np.max(df['UH'])
		else:
	        	scaleX = np.max(df[yaxis[jj].upper()])
        for i in range(len(df.values[:])):
		if 'uh'.upper() in yaxis[jj].upper():
			tmpx = df['UH'][i]
		else:	
			tmpx = df[yaxis[jj].upper()][i]
		if tmpx>1:
			tmpx =(tmpx-1)/scaleX
		else:
			tmpx = 1-tmpx
    		xval+=[tmpx]
    		yval+=[jj]
		hybridV+=[df['EM/WO'][i]/nics*100.]
        ylabel+=[yaxis[jj]]


        labelTmp=getStateListfromFile(df)
        compHybrid=getCompHybridVal(compI=False)

        lens = []
        for i in range(len(labelTmp)):
            lens+=[len(labelTmp[i])]
        if 9 in lens:
              if 0 not in xval and ('input' not in fileN) and ('input' not in fsave):
                    xval+=[0]
                    yval+=[jj]
                    ind = np.argwhere(np.array(lens)==9)[:,0][0]
                    labelTmp+=[np.array(labelTmp)[ind]]
	    	    hybridV+=[compHybrid]
        elif 0 not in xval and ('input' not in fileN) and ('input' not in fsave):
            xval+=[0]
            yval+=[jj]
            labelTmp+=[['E/O','E/WO','E/W','M/O','M/WO','M/W','EM/O','EM/WO','EM/W']]
	    hybridV+=[compHybrid]
        inds = np.argsort(xval)
        xval = np.array(xval)[inds]
        yval = np.array(yval)[inds]
        label= np.array(labelTmp)[inds]
        hybridV= np.array(hybridV)[inds]

        x1,x2,y1,y2,label,hybridV=getDataForPlot(xval,yval,label,hybridV)
	x1L+=x1
	x2L+=x2
	y1L+=y1
	y2L+=y2
	labelL+=label
	hybridVL+=hybridV

    color,star_colors,colList = getPlotData(labelL)
    hatch,edgcolors,regType=getHatchForPlot(star_colors,color,hybridVL,compHybrid,colList)

    fig = plt.figure(figsize=(20,14))
    mpl.rc('hatch',linewidth=hlw)
    for i in range(len(x1L)):
            	plt.fill_between([x1L[i],x2L[i]],y1L[i],y2L[i],facecolor=color[i],edgecolor=edgcolors[i],hatch=hatch[i],linewidth=0.0)
    legend_elements = getLegend(regType)

    plt.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(title)
    plt.yticks([])
    plt.xlabel(xlabel,fontsize=40)
    plt.yticks(np.arange(0,len(ylabel)),ylabel)
    fig.savefig(fsave+"stateProg.png",bbox_inches='tight')
    plt.close()

###################################
###################################
############################################
###################################
###################################
###################################
def plotCoupledPhenotypes_singleLink(fileN,xaxis,title,fsave,checkRanges=False,compI=False,xlim=[],ylim=[]):

    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    xlabel = '$\lambda_{'+xaxis+'}$'
    xval,yval,label=[],[],[]
    hybridV = []
    df = pd.read_csv(fileN,header=0)
    #if 'G' in df.columns:
    #    df=df.drop(columns=['G','O','mg','mo','t'])

    nics = float(df['nics'][0])
    for i in range(len(df.values[:])):
	xval+=[df[xaxis.upper()][i]]
	yval+=[0]
	hybridV+=[df['EM/WO'][i]/nics*100.]

    compHybrid=getCompHybridVal(compI=False)

    label=getStateListfromFile(df)

    lens = []
    for i in range(len(label)):
        lens+=[len(label[i])]
    if 9 in lens:
          if 1 not in xval and ('input' not in fileN) and ('input' not in fsave) and ('PSF' not in fileN) and ('PSF' not in fsave) :
                xval+=[1.]
                yval+=[0.]
                ind = np.argwhere(np.array(lens)==9)[:,0][0]
                label+=[np.array(label)[ind]]
		hybridV+=[compHybrid]
    elif 1 not in xval and ('input' not in fileN) and ('input' not in fsave) and ('PSF' not in fileN) and ('PSF' not in fsave) :
        xval+=[1.]
        yval+=[0.]
        label+=[['E/O','E/WO','E/W','M/O','M/WO','M/W','EM/O','EM/WO','EM/W']]
	hybridV+=[compHybrid]


    if "input".upper() not in fileN.upper():
        plotPhasePlane(xval,yval,label,hybridV,checkRanges,title,xlabel,'',compHybrid,fsave,xlim,ylim)
    else:
        plotPhasePlane(xval,yval,label,hybridV,checkRanges,title,xlabel,'',compHybrid,fsave,xlim,ylim,hatching=False)
###################################
###################################
def plotICSvlambda_RegComp(filen,xaxis,title,fsave,scaleList=[]):

    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams.update({'font.size': 20})

    xdir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    ydir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    labels={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    col={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    

    for jj in range(len(filen)):
        df = pd.read_csv(filen[jj],header=0)
        nics = df['nics'][0]
	if len(scaleList)>0:
		scaleX=scaleList[jj]
	else:
		if 'uh'.upper() in xaxis[jj].upper():
	        	scaleX = np.max(df['UH'])
		else:
	        	scaleX = np.max(df[xaxis[jj].upper()])

        for i in range(len(df)):
        	for key in xdir:
			if 'uh'.upper() in xaxis[jj].upper():
        			tmpx=df['UH'][i]
			else:
 	       			tmpx=df[xaxis[jj].upper()][i]
			if tmpx>1:
				tmpx =(tmpx-1)/scaleX
			else:
				tmpx = 1-tmpx
 	       		xdir[key]+=[tmpx]
        		ydir[key]+=[df[key][i]/nics*100.]
			col[key]+=[__clist[jj]]
			labels[key]+=[xaxis[jj]]
                
    if 'input' in fsave.lower():
	compI=xdir['E'][0]## get the value of input
	comp=True
    else:
 	comp=True
	compI=False

    xlabel="Increasing regulation"
    __plotSet1(xdir,ydir,xlabel,title,fsave,col=col,labels=labels,comp=comp,compI=compI)
    __plotSet2(xdir,ydir,xlabel,title,fsave,col=col,labels=labels,comp=comp,compI=compI)
###################################
###################################
###################################
####################################
def plotICSvlambda_singleLinkComp(fileList,xaxis,title,fsave,scaleList):

    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams.update({'font.size': 20})

    xdir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    ydir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    col={}
    for key in xdir:
	col[key]=[]

    for jj in range(len(fileList)):
        df = pd.read_csv(fileList[jj],header=0)
        nics = float(df['nics'][0])
	if len(scaleList)>0:
		scaleX=scaleList[jj]
	else:
		if 'uh'.upper() in xaxis[jj].upper():
	        	scaleX = np.max(df['UH'])
		else:
	        	scaleX = np.max(df[xaxis[jj].upper()])
        for i in range(len(df)):
		if 'uh'.upper() in xaxis[jj].upper():
			tmpx = df['UH'][i]
		else:	
			tmpx = df[xaxis[jj].upper()][i]
		if tmpx>1:
			tmpx =(tmpx-1)/scaleX
		else:
			tmpx = 1-tmpx
    	        for key in xdir:
    			xdir[key]+=[tmpx]
    			ydir[key]+=[df[key][i]/nics*100.]
			col[key]+=[__clist[jj]]
            

    xlabel="$\lambda_{"+str(xaxis)+"}$"
    comp=True
    compI=False
    __plotSet1(xdir,ydir,xlabel,title,fsave,comp=comp,compI=compI,col=col)
    __plotSet2(xdir,ydir,xlabel,title,fsave,comp=comp,compI=compI,col=col)
###################################
###################################
####################################
def plotICSvlambda_singleLink(filen,xaxis,title,fsave):

    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams.update({'font.size': 20})

    xdir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    ydir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    

    df = pd.read_csv(filen,header=0)
    #if 'G' in df.columns:
    #    df=df.drop(columns=['G','O','mg','mo','t'])
    nics = df['nics'][0]

    for i in range(len(df)):
	for key in xdir:
		xdir[key]+=[df[xaxis.upper()][i]]
		ydir[key]+=[df[key][i]/nics*100.]
        
    if 'input' in fsave.lower():
	compI=xdir['E'][0]## get the value of input
	comp=True
    else:
 	comp=True
	compI=False

    xlabel="$\lambda_{"+str(xaxis)+"}$"
    __plotSet1(xdir,ydir,xlabel,title,fsave,comp=comp,compI=compI)
    __plotSet2(xdir,ydir,xlabel,title,fsave,comp=comp,compI=compI)
###################################
###################################
###################################
###################################
###################################
###################################
def plotICSvlambda_doubleLink(filen,title,fsave,compI=False,comp=True):

    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams.update({'font.size': 20})

    xdir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    ydir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    zdir={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    labels={'E/O':[],'E/WO':[],'E/W':[],
              'M/O':[],'M/WO':[],'M/W':[],
             'EM/O':[],'EM/WO':[],'EM/W':[],
             'E':[],'EM':[],'M':[],
             'O':[],'WO':[],'W':[]}
    
    if comp:
	if PSF:
    		res_9state=getStates(df,PSF=True)
    	if not compI:
	        df =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    		res_9state=getStates(df)
	else:
		df = pd.read_csv("data/crosstalk_input.txt",header=0)
		tempStates={}
		row = np.argwhere(df['INPUT'].values==compI)[:,0][0]
		for key in ['E','EM','M','W','O','WO','E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']:
			res_9state[key] = df[key][row]

    df = pd.read_csv(filen,header=0)
    nics = df['nics'][0]

    ## header is AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH
    keyList = fsave.split("_")[1:]
    for i in range(len(keyList)):
	if keyList[i].upper()=='IHU':
		keyList[i]='Hu'

    constants={}
    fsaveDir={}
    el = keyList[0].upper()
    el2 = keyList[1].upper()
    constants[el] = np.unique(df[el.upper()])
    fsaveDir[el] = fsave.split("_")[0]+"_"+el2+'_'+el
    constants[el2] = np.unique(df[el2.upper()])
    fsaveDir[el2] = fsave.split("_")[0]+"_"+el+'_'+el2
    xaxis={el2:el,el:el2}

    for elKey in constants:
	for jj in range(len(constants[elKey])):
		inds = list(np.argwhere(df[elKey.upper()].values==constants[elKey][jj])[:,0])
    		xlabel = '$\lambda_{'+str(xaxis[elKey])+'}$'

                if len(inds)>0:
                    for key in xdir:
    	    	        inds = np.unique(inds)
			xdir[key],ydir[key]=[],[]
    	    	        for i in range(len(df.values[inds])):
    	    			xdir[key]+=[df[xaxis[elKey].upper()].values[inds][i]]
    	    			ydir[key]+=[df[key].values[inds][i]/nics*100.]
                else:
                    for key in xdir:
			xdir[key],ydir[key]=[],[]
                        for i in range(len(df.values[:])):
    	    			xdir[key]+=[df[xaxis[elKey].upper()][i]]
    	    			ydir[key]+=[df[key][i]/nics*100.]

                for key in xdir:
                    xdir[key]=np.array(xdir[key])
                    ydir[key]=np.array(ydir[key])
                
		fsave2 = fsaveDir[elKey]+'_'+str(constants[elKey][jj])
                __plotSet1(xdir,ydir,xlabel,title,fsave2,comp=comp,compI=compI)
                __plotSet2(xdir,ydir,xlabel,title,fsave2,comp=comp,compI=compI)

#################
#################
#################
#################333
####################3
def plotCoupledPhenotypes_doubleLink(fileN,xaxis,yaxis,title,fsave,constants={},xlim=[-1],ylim=[-1],checkRanges=False):

    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    xval,yval,lab=[],[],[]
    hybridV = []
    df = pd.read_csv(fileN,header=0)
    nics = float(df['nics'][0])
    inds =[]
    if len(constants.keys())>0:
	for key in constants:
		inds += list(np.argwhere(df[key.upper()].values==constants[key])[:,0])

    if len(inds)>0:
	    inds = np.unique(inds)
	    for i in range(len(df.values[inds])):
		xval+=[df[xaxis.upper()].values[inds][i]]
		yval+=[df[yaxis.upper()].values[inds][i]]
		hybridV+=[df['EM/WO'].values[inds][i]/nics*100.]
    else:
        for i in range(len(df.values[:])):
    		xval+=[df[xaxis.upper()][i]]
	    	yval+=[df[yaxis.upper()][i]]
		hybridV+=[df['EM/WO'][i]/nics*100.]

    lab=getStateListfromFile(df)

    xlabel='$\lambda_{'+xaxis+'}$'
    ylabel='$\lambda_{'+yaxis+'}$'
    compHybrid=getCompHybridVal(compI=False)

    plotPhasePlane(xval,yval,lab,hybridV,checkRanges,title,xlabel,ylabel,compHybrid,fsave,xlim,ylim)

####################3
####################3
####################3
####################3
####################3
###########################
###########################

def plotCoupledPhenotypes_regulatory(fileN,title,fsave,xregList=['AS','AZ','AU','HS','HU','INPUT'],ylim=[-1,-1],xlim=[-1,-1]):

    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    ###regs = ['AS','AZ','Au','HS','HU','input','u3m','u3n','uh']
    #xregList = ['AS','AZ','AU','HS','HU','INPUT']
    #yregList=['U3M','U3N','UH']
    
    tickList={'x':{},'y':{}}
    label_dir={'x':{},'y':{}}
    opp_label={'x':{},'y':{}}
    tickLabels={'x':{},'y':{}}

    df = pd.read_csv(fileN,header=0)

    # get all active links
    reg_list=fileN.split(".")[0].split("_")[1:]
    if reg_list[0]=='crosstalk':
	reg_list=reg_list[1:]

    ## split the list of links based on x and y (try to get half and half):input manually or EMT on x
    reg_dict={'x':[],'y':[]}
    for i in range(len(reg_list)):
        if reg_list[i].upper() in xregList:
            reg_dict['x']+=[reg_list[i].upper()]
        elif reg_list[i].upper()=='IHU':
            reg_dict['x']+=['HU']
        else:
            reg_dict['y']+=[reg_list[i].upper()]

    ### get the values of crosstalk associated with the reg_dict
    ### label_dir[key][i] are foldchange lists for reg_dict[key][i]
    for key in ['x','y']:
        for i in range(len(reg_dict[key])):
            label_dir[key][i]= list(np.sort(np.unique(df[reg_dict[key][i].upper()])))

    ####now rearrange to get xvals,yvals,labels, and hybridV
    xval,yval,label,hybridV,xlabels,ylabels=getPlotData_regulatory(df,label_dir,reg_dict)
    compHybrid=getCompHybridVal(compI=False)

    plotPhasePlane(xval,yval,label,hybridV,False,title,'EMT reg','MR reg',compHybrid,fsave,[],[])

####################3
###########################
###########################
def output_dataVals(dirN,fsave=''):

    if not fsave:
	dirn_red = dirN.split("/")[-1]
	fsave = "data/"+dirn_red+".txt"
    nics= dirN.split("_")[-2]
    fileo=open(fsave,"w")

    fileo.write("AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV,u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,E,EM,M,W,WO,O\n")

    labels=['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','LAMDAUH','LAMBDAUH']
    for filen in os.listdir(dirN):
        if 'res.txt' in filen:
            tmp= filen.split("_")
            df =pd.read_csv(dirN+filen).dropna()
	    cxs = {'AS':1.,'AZ':1.,'AU':1.,'HS':1.,'HU':1.,'INPUT':50000,'U3M':1.,'U3N':1.,'UH':1.,'UHV':0}
	    cxs,lamda_counted=add_lamdas(cxs,tmp,df)
	    if not lamda_counted:
		tmp = [dirN.replace("/",'').split("_")[-1],tmp[3],tmp[4],tmp[-2],'0']
		cxs,lamda_counted = add_lamdas(cxs,tmp,df)

	    try:
                #if 'G' in df.columns:
                # 	res=getStates_ind(df.drop(columns=['G','O','mg','mo','t']))
                if True:#else:
                 	res=getStates_ind(df)

		###u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh
		tmp=np.unique(np.round(df.values,2),axis=0)
		for i in range(len(tmp)):
			inds = np.argwhere(tmp[i]==np.round(df.values,3))[:,0]
			for j in range(len(inds)):
			        tmpRES={'E':0,'EM':0,'M':0,'W':0,'WO':0,'O':0}
				emv,mrv = res[inds[j]].split("/")
				tmpRES[emv]=1
				tmpRES[mrv]=1
	         		fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(cxs['AS'],cxs['AZ'],cxs['AU'],cxs['HS'],cxs['HU'],cxs['INPUT'],cxs['U3M'],cxs['U3N'],cxs['UH'],cxs['UHV'],df['u'].values[inds[j]],df['mz'].values[inds[j]],df['Z'].values[inds[j]],df['ms'].values[inds[j]],df['u3'].values[inds[j]],df['S'].values[inds[j]],df['A'].values[inds[j]],df['Rmt'].values[inds[j]],df['Rnox'].values[inds[j]],df['h'].values[inds[j]],df['mh'].values[inds[j]],tmpRES['E'],tmpRES['EM'],tmpRES['M'],tmpRES['W'],tmpRES['WO'],tmpRES['O']))
	    except:
		print "Issue with ", dirN,filen

    fileo.close()


#########################
#########################
def output_results(dirN,fsave='',skipUH=None):

    if not fsave:
	dirn_red = dirN.split("/")[-1]
	if (not skipUH) or ('UH' in dirN.upper()):
		fsave = "data/IND_"+dirn_red+".txt"
	else:
		fsave = "data/"+dirn_red+".txt"
    nics= dirN.split("_")[-2]
    fileo=open(fsave,"w")

    fileo.write("AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV,E,EM,M,W,WO,O,M/W,M/WO,M/O,EM/W,EM/WO,EM/O,E/W,E/WO,E/O,nics\n")

    labels=['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','LAMDAUH','LAMBDAUH']
    for filen in os.listdir(dirN):
        if 'res.txt' in filen:
            tmp= filen.split("_")
            df =pd.read_csv(dirN+filen).dropna()
	    cxs = {'AS':1.,'AZ':1.,'AU':1.,'HS':1.,'HU':1.,'INPUT':50000,'U3M':1.,'U3N':1.,'UH':1.,'UHV':0}
	    cxs,lamda_counted=add_lamdas(cxs,tmp,df,skipUH)
	    if not lamda_counted:
		tmp = [dirN.replace("/",'').split("_")[-1],tmp[3],tmp[4],tmp[-2],'0']
		cxs,lamda_counted = add_lamdas(cxs,tmp,df,skipUH)

            nics=float(tmp[-2])
            if len(df.values[:])!=nics:
		print len(df.values[:]),'!=',nics
	    #print dirN, filen

            if 'noHH_' in fsave:
		noHH='noHH'
            elif 'noEM_' in fsave:
		noHH='noEM'
            elif 'noWO_' in fsave:
		noHH='noWO'
	    if True:#try:
		 if 'UH' in filen.upper() and not skipUH:
                      if 't' in df.columns:
                      	res=getStates_uh(df.drop(columns=['G','O','mg','mo','t']),PSF=True)
	              elif 'G' in df.columns:
                      #	res=getStates_uh(df.drop(columns=['G','O','mg','mo']))
                      	res=getStates_uh(df,PSF=True)
                      else:
                      	res=getStates_uh(df)#,noHH=noHH) using noHH doesn't increase smoothness
	              ##fileo.write("AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV,E,EM,M,W,WO,O,M/W,M/WO,M/O,EM/W,EM/WO,EM/O,E/W,E/WO,E/O,nics\n")
                      for i in range(len(res)):
	                  fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(cxs['AS'][i],cxs['AZ'][i],cxs['AU'][i],cxs['HS'][i],cxs['HU'][i],cxs['INPUT'][i],cxs['U3M'][i],cxs['U3N'][i],cxs['UH'][i],cxs['UHV'][i],res[i]['E'],res[i]['EM'],res[i]['M'],res[i]['W'],res[i]['WO'],res[i]['O'],res[i]['M/W'],res[i]['M/WO'],res[i]['M/O'],res[i]['EM/W'],res[i]['EM/WO'],res[i]['EM/O'],res[i]['E/W'],res[i]['E/WO'],res[i]['E/O'],nics))

		 else:
                      if 't' in df.columns:
                      	res=getStates(df.drop(columns=['G','O','mg','mo','t']),PSF=True)
	              elif 'G' in df.columns:
                      	res=getStates(df,PSF=True)
                      else:
                      	res=getStates(df)#,noHH=noHH) using noHH doesn't increase smoothness
	              fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(cxs['AS'],cxs['AZ'],cxs['AU'],cxs['HS'],cxs['HU'],cxs['INPUT'],cxs['U3M'],cxs['U3N'],cxs['UH'],cxs['UHV'],res['E'],res['EM'],res['M'],res['W'],res['WO'],res['O'],res['M/W'],res['M/WO'],res['M/O'],res['EM/W'],res['EM/WO'],res['EM/O'],res['E/W'],res['E/WO'],res['E/O'],nics))
	    else:#xcept:
		print "Issue with ", dirN,filen

    fileo.close()


#########################
#########################
def plotPhasePlane(xval,yval,label,hybridV,checkRanges,title,xlabel,ylabel,compHybrid,fsave,xlim=[],ylim=[],hatching=True,xticks=[],yticks=[],tickLabels=[]):
    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    label= np.array(label)[inds]
    hybridV= np.array(hybridV)[inds]

    x1,x2,y1,y2,label,hybridV=getDataForPlot(xval,yval,label,hybridV)
    color,star_colors,colList = getPlotData(label)
    hatch,edgcolors,regType=getHatchForPlot(star_colors,color,hybridV,compHybrid,colList)

    fig = plt.figure(figsize=(20,14))
    ax1 = plt.subplot(111)

    if checkRanges:
        for i in range(len(x1)):
            ax1.fill_between([x1[i],x2[i]],y1[i],y2[i],facecolor=color[i])    
	for i in range(len(xval)):
            ax1.plot(xval[i],yval[i],'o',markersize=20,markeredgecolor='w',markerfacecolor='k',markeredgewidth=4)
    elif hatching:
    	mpl.rc('hatch',linewidth=hlw)
        for i in range(len(x1)):
            	ax1.fill_between([x1[i],x2[i]],y1[i],y2[i],facecolor=color[i],edgecolor=edgcolors[i],hatch=hatch[i],linewidth=0.0)
    else:
        for i in range(len(x1)):
            ax1.fill_between([x1[i],x2[i]],y1[i],y2[i],facecolor=color[i])    

    legend_elements = getLegend(regType)

    ax1.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.set_title(title)
    if ylabel=='':
    	ax1.set_yticks([])

    if xticks==[]:
        if len(xlim)>1:
        	ax1.set_xlim(xlim[0],xlim[1])
        else:
        	ax1.set_xlim(np.min(x1),np.max(x2))
    else:
	axes = []
	if len(xticks.keys())>2:
		axes +=[ax1.twiny()]
	for ax in axes:
    		mpl.rc('hatch',linewidth=hlw)
		ax.fill_between([x1[0],x2[0]],y1[0],y2[0],facecolor=color[0],edgecolor=edgcolors[0],hatch=hatch[0],linewidth=0.0,alpha=0.0)
		ax.xaxis.set_ticks_position('bottom')
	
	added0,added1=False,False
	count=0
	for key in xticks:
		if key%2==0:
			if count==0:
				ax1.xaxis.set_minor_locator(ticker.FixedLocator(np.array(xticks[key])))
				ax1.tick_params(which='minor',length=15,color='k')
    				ax1.set_xticklabels([])
			else:
				axes[count-1].xaxis.set_minor_locator(ticker.FixedLocator(np.array(xticks[key])))
				axes[count-1].tick_params(which='minor',length=15,color='k')
    				axes[count-1].set_xticklabels([])
			added0=True
		else:##key%2==1
			if count==0:
				ax1.xaxis.set_major_locator(ticker.FixedLocator(np.array(xticks[key])))
				ax1.tick_params(which='major',length=45,color='k')
			else:
				axes[count-1].xaxis.set_major_locator(ticker.FixedLocator(np.array(xticks[key])))
				axes[count-1].tick_params(which='major',length=25,color='k')
    				axes[count-1].set_xticklabels([])
			added1=True
		if added0 and added1:
			count+=1
			added0,added1=False,False
	if len(xticks.keys())>1 and len(tickLabels['x'])<50:
    		ax1.set_xticklabels(tickLabels['x'],minor=True,rotation=75)
	else:
    		ax1.set_xticklabels(tickLabels['x'],minor=False,rotation=75)


    if yticks==[]:
        if len(ylim)>1:
        	ax1.set_ylim(ylim[0],ylim[1])
        else:
        	ax1.set_ylim(np.min(y1),np.max(y2))
    else:
	axes = []
	if len(yticks.keys())>2:
		axes +=[ax1.twinx()]
	for ax in axes:
    		mpl.rc('hatch',linewidth=hlw)
		ax.fill_between([x1[0],x2[0]],y1[0],y2[0],facecolor=color[0],edgecolor=edgcolors[0],hatch=hatch[0],linewidth=0.0,alpha=0.0)
		ax.yaxis.set_ticks_position('left')

	count=0
	added0,added1=False,False
	for key in yticks:
		if key%2==0:
			if count==0:
				ax1.yaxis.set_minor_locator(ticker.FixedLocator(np.array(yticks[key])))
				ax1.tick_params(which='minor',length=15,color='k')
    				ax1.set_yticklabels([])
			else:
				axes[count-1].yaxis.set_minor_locator(ticker.FixedLocator(np.array(yticks[key])))
				axes[count-1].tick_params(which='minor',length=15,color='k')
    				axes[count-1].set_yticklabels([])
			added0=True
		else:##key%2==1
			if count==0:
				ax1.yaxis.set_major_locator(ticker.FixedLocator(np.array(yticks[key])))
				ax1.tick_params(which='major',length=45,color='k')
			else:
				axes[count-1].yaxis.set_major_locator(ticker.FixedLocator(np.array(yticks[key])))
				axes[count-1].tick_params(which='major',length=25,color='k')
    				axes[count-1].set_yticklabels(tickLabels['y'][key])
			added1=True
		if added0 and added1:
			count+=1
			added0,added1=False,False
	if len(xticks.keys())>1 and len(tickLabels['x'])<50:
    		ax1.set_yticklabels(tickLabels['y'],minor=True,rotation=25)
	else:
    		ax1.set_yticklabels(tickLabels['y'],minor=False,rotation=25)


    ax1.set_xlabel(xlabel,fontsize=40)
    if type(ylabel)==list:
        ax1.set_yticks(np.arange(0,len(ylabel)),ylabel)
    else:
    	ax1.set_ylabel(ylabel,fontsize=40)
    if checkRanges:
    	fig.savefig(fsave+"stateProg_ranges.png",bbox_inches='tight')
    else:
    	fig.savefig(fsave+"stateProg.png",bbox_inches='tight')
    plt.close()
