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



##########
##########
##########
####def plotValuesEMTvMR_full(filen,save=False,compI=False):
####
####def plotICSvlambda_singleLink(dirN,xlabel,title,fsave):
####def plotICSvlambda_doubleLink(dirN,xlabel,title,fsave,compI=False):
####def    __plotSet1(xdir,ydir,col,xlabel,title,fsave,comp=False,labels=False,xlim=[-1],ylim=[-1],compI=False):
####def    __plotSet2(xdir,ydir,col,xlabel,title,fsave,comp=False,labels=False,compI=False):
####
####def plotCoupledPhenotypes_singleLink(dirN,xlabel,title,fsave,logScale=False,onlyShowHybrid=False):
####def plotCoupledPhenotypes_doubleLink(fileN,xaxis,yaxis,title,fsave,xlim=[-1],ylim=[-1],onlyShowHybrid=False):
####def plotCoupledPhenotypes_regulatory(fileN,title,fsave,xlog=False,ylog=False,xlim=[-1],ylim=[-1],onlyShowHybrid=False):
####
####def output_results(dirN,fsave=''):
####
##########
##########
##########

def    __plotSet1(xdir,ydir,col,xlabel,title,fsave,comp=False,labels=False,xlim=[-1],ylim=[-1],compI=False):
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
        
    compV = np.min(xdir['E/O'])
    for i in range(len(xdir['E/O'])):
        if labels and xdir['E/O'][i]==compV:
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
    
    axC9.legend()
    fig.savefig(fsave+"_LvCS.png",bbox_inches='tight')
    plt.show()

###################################
###################################
###################################
#########################
#########################

def    __plotSet2(xdir,ydir,col,xlabel,title,fsave,comp=False,labels=False,compI=False):
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
    
    compV = np.min(xdir['E'])
    for i in range(len(xdir['E/O'])):
        if labels and xdir['E'][i]==compV:
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
    
    axS6.legend()
    fig.savefig(fsave+"_LvOS.png",bbox_inches='tight')
    plt.show()



def plotValuesEMTvMR_full(filen,save=False,compI=False):
    
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

        if comp:
    	    	if not compI:
		        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
		else:
		        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_input/EMT_MR_comp_input_"+str(compI)+"_1000_res.txt").dropna()

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
###################################
###################################
###################################
def plotCoupledPhenotypes_singleLink(fileN,xaxis,title,fsave,logScale=False,onlyShowHybrid=False):

    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    xval,yval,label=[],[],[]
    df = pd.read_csv(fileN,header=0)
    print len(df.values[:]),df.columns
    for i in range(len(df.values[:])):
	xval+=[df[xaxis.upper()][i]]
	yval+=[0]

    xlabel='$\lambda_{'+xaxis+'}$'

    label=getStateListfromFile(df,onlyShowHybrid) 

    lens = []
    for i in range(len(label)):
        lens+=[len(label[i])]
    if 9 in lens:
          if 1 not in xval:
                xval+=[1.]
                yval+=[0.]
                ind = np.argwhere(np.array(lens)==9)[:,0][0]
                label+=[np.array(label)[ind]]
    elif 1 not in xval:
        xval+=[1.]
        yval+=[0.]
        label+=[['E/O','E/WO','E/W','M/O','M/WO','M/W','EM/O','EM/WO','EM/W']]

    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    label= np.array(label)[inds]


    x1,x2,y1,y2,label=getDataForPlot(xval,yval,label)
    color,star_colors,colList = getPlotData(label)
    star_x,star_y=getStarsForPlot(x1,x2,y1,y2,star_colors,color)

    fig = plt.figure(figsize=(20,14))

    for i in range(len(x1)):
        plt.fill_between([x1[i],x2[i]],y1[i],y2[i],facecolor=color[i])    
    for i in range(len(star_x)):
    	plt.plot(star_x[i],star_y[i],'*',markersize=70,markeredgecolor='w',markerfacecolor='k',markeredgewidth=4)
	
    legend_elements = getLegend(colList)

    plt.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(title)
    plt.yticks([])
    plt.ylim(min(y1),max(y2))
    plt.xlim(min(x1),max(x2))
    if logScale:
        plt.xscale('log')
	xlab = np.arange(1.,np.max(xval)+1)
	plt.xticks([1,2,4,7,10],[1,2,4,7,10])
    	plt.xlabel(xlabel+" (log10)",fontsize=40)
    else:
    	plt.xlabel(xlabel,fontsize=40)
    if onlyShowHybrid:
    	fig.savefig(fsave+"stateProg_grouped.png",bbox_inches='tight')
    else:
    	fig.savefig(fsave+"stateProg.png",bbox_inches='tight')
    #plt.show()
###################################
###################################
###################################
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
    nics = df['nics'][0]

    for i in range(len(df)):
	for key in xdir:
		xdir[key]+=[df[xaxis.upper()][i]]
		ydir[key]+=[df[key][i]/nics*100.]
        
    col={}
    for key in xdir:
        xdir[key]=np.array(xdir[key])
        ydir[key]=np.array(ydir[key])
	col[key]=[]
	for i in range(len(xdir[key])):
	    col[key]+=['k']

    if 'input' in fsave.lower():
	compI=xdir['E'][0]## get the value of input
	comp=True
    else:
 	comp=True
	compI=False

    xlabel="$\lambda_{"+str(xaxis)+"}$"
    __plotSet1(xdir,ydir,col,xlabel,title,fsave,comp=comp,compI=compI)
    __plotSet2(xdir,ydir,col,xlabel,title,fsave,comp=comp,compI=compI)
###################################
###################################
###################################
###################################
###################################
###################################
def plotICSvlambda_doubleLink(dirN,xlabel,title,fsave,compI=False):

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
    	if not compI:
	        df =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    		res_9state=getStates(df)
	else:
		df = pd.read_csv("data/crosstalk_input.txt",header=0)
		tempStates={}
		row = np.argwhere(df['INPUT'].values==compI)[:,0][0]
		for key in ['E','EM','M','W','O','WO','E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']:
			res_9state[key] = df[key][row]

    for filen in os.listdir(dirN):
        if 'res.txt' in filen:
            tmp= filen.split("_")
	    start=3
	    if "comp" not in filen:
		start = 2
            lamdaX = float(tmp[start+1])+float(tmp[start+2])/10**(len(tmp[start+2]))
	    if 'uH' in dirN or 'uh' in dirN:
                lamdaY = float(tmp[start+4])
                nics=float(tmp[start+5])
	    else:
                lamdaY = float(tmp[start+4])+float(tmp[start+5])/10**(len(tmp[start+5]))
                nics=float(tmp[start+6])
            df =pd.read_csv(dirN+filen).dropna()
            if len(df.values[:])==nics:
                res=getStates(df)
                for key in res:
                    ydir[key]+=[res[key]/nics*100.]
                    xdir[key]+=[lamdaX]
                    zdir[key]+=[lamdaY]
                    labels[key]+=[lamdaY]
                    
            if 1. not in xdir['E']:
                for key in res_9state:
                    ydir[key]+=[res_9state[key]/nics*100.]
                    xdir[key]+=[1.]
                    zdir[key]+=[1.]
                    labels[key]+=[1.]
        
    cols = np.unique(zdir['E'])
    for key in zdir:
        for i in range(len(zdir[key])):
            el = zdir[key][i]
            ind = np.argwhere(cols==el)[:,0][0]
            zdir[key][i]=__clist[ind]
    
    for key in xdir:
        xdir[key]=np.array(xdir[key])
        ydir[key]=np.array(ydir[key])
        zdir[key]=np.array(zdir[key])
        labels[key]=np.array(labels[key])
    __plotSet1(xdir,ydir,zdir,xlabel,title,fsave,comp=True,labels=labels,compI=compI)
    __plotSet2(xdir,ydir,zdir,xlabel,title,fsave,comp=True,labels=labels,compI=compI)

#################
#################
#################
#################333
####################3
def plotCoupledPhenotypes_doubleLink(fileN,xaxis,yaxis,title,fsave,constants={},xlim=[-1],ylim=[-1],onlyShowHybrid=False):

    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    xval,yval,lab=[],[],[]
    df = pd.read_csv(fileN,header=0)
    inds =[]
    if len(constants.keys())>0:
	for key in constants:
		inds += list(np.argwhere(df[key.upper()].values==constants[key])[:,0])

    if len(inds)>0:
	    inds = np.unique(inds)
	    for i in range(len(df.values[inds])):
		xval+=[df[xaxis.upper()].values[inds][i]]
		yval+=[df[yaxis.upper()].values[inds][i]]
    else:
        for i in range(len(df.values[:])):
    		xval+=[df[xaxis.upper()][i]]
	    	yval+=[df[yaxis.upper()][i]]

    xlabel='$\lambda_{'+xaxis+'}$'
    ylabel='$\lambda_{'+yaxis+'}$'


    lab=getStateListfromFile(df,onlyShowHybrid) 

    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    lab= np.array(lab)[inds]


    label=np.array(lab)
    x1,x2,y1,y2,label=getDataForPlot(xval,yval,label)
    color,star_colors,colList = getPlotData(label)
    star_x,star_y=getStarsForPlot(x1,x2,y1,y2,star_colors,color)


    fig = plt.figure(figsize=(20,14))
    for i in range(len(x1)):
        plt.fill_between([x1[i],x2[i]],y1[i],y2[i],facecolor=color[i])    

    for i in range(len(star_x)):
    	plt.plot(star_x[i],star_y[i],'*',markersize=70,markeredgecolor='w',markerfacecolor='k',markeredgewidth=4)
            
    legend_elements = getLegend(colList)
    print getColorMatch(colList)
    print np.unique(color)
    for i in range(len(color)):
	print x1[i],x2[i],y1[i],y2[i],color[i]
            
    plt.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(title)
    plt.xlabel(xlabel,fontsize=40)
    plt.ylabel(ylabel,fontsize=40)
    if len(xlim)>1:
    	plt.xlim(xlim[0],xlim[1])
    else:
    	plt.xlim(np.min(x1),np.max(x2))
    if len(ylim)>1:
    	plt.ylim(ylim[0],ylim[1])
    else:
    	plt.ylim(np.min(y1),np.max(y2))

    if onlyShowHybrid:
    	fig.savefig(fsave+"stateProg_grouped.png",bbox_inches='tight')
    else:
    	fig.savefig(fsave+"stateProg.png",bbox_inches='tight')
    #plt.show()

####################3
####################3
####################3
####################3
####################3
###########################
###########################

def plotCoupledPhenotypes_regulatory(fileN,title,fsave,xlog=False,ylog=False,xlim=[-1],ylim=[-1],onlyShowHybrid=False):

    fileo=open("tmp","w")
    mpl.rcParams['xtick.labelsize'] = 40
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams.update({'font.size': 30})

    xval,yval,lab=[],[],[]
    df = pd.read_csv(fileN,header=0)

    ### FIX THIS
    for i in range(len(df.values[:])):
	xval+=[df['HS'][i]]
	yval+=['uh']

    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    lab= np.array(lab)[inds]

    print consistent_states(xval,yval,lab)

    label=np.array(lab)

    x1,x2,y1,y2,label=getDataForPlot(xval,yval,label)
    color,star_colors,colList = getPlotData(label)
    star_x,star_y=getStarsForPlot(x1,x2,y1,y2,star_colors,color)
    
    fig = plt.figure(figsize=(20,14))
    for i in range(len(x1)):
        plt.fill_between([x1[i],x2[i]],y1[i],y2[i],facecolor=color[i])    
    for i in range(len(star_x)):
    	plt.plot(star_x[i],star_y[i],'*',markersize=70,markeredgecolor='w',markerfacecolor='k',markeredgewidth=4)
            
    legend_elements = getLegend(colList)

    plt.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(title)
    plt.xlabel('EMT_reg',fontsize=40)
    plt.ylabel('MR_reg',fontsize=40)
    if ylog:
	plt.yscale('log')
    if xlog:
	plt.xscale('log')
    if len(xlim)>1:
    	plt.xlim(xlim[0],xlim[1])
    else:
    	plt.xlim(np.min(x1),np.max(x2))
    if len(ylim)>1:
    	plt.ylim(ylim[0],ylim[1])
    else:
    	plt.ylim(np.min(y1),np.max(y2))
    if onlyShowHybrid:
    	fig.savefig(fsave+"stateReg_grouped.png",bbox_inches='tight')
    else:
    	fig.savefig(fsave+"stateReg.png",bbox_inches='tight')
    #plt.show()
    fileo.close()

def output_results(dirN,fsave=''):

    if not fsave:
	dirn_red = dirN.split("/")[-1]
	fsave = "data/"+dirn_red+".txt"
    nics= dirN.split("_")[-2]
    fileo=open(fsave,"w")

    fileo.write("AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,E,EM,M,W,WO,O,M/W,M/WO,M/O,EM/W,EM/WO,EM/O,E/W,E/WO,E/O,nics\n")

    labels=['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','LAMDAUH','LAMBDAUH']
    for filen in os.listdir(dirN):
        if 'res.txt' in filen:
            tmp= filen.split("_")
            df =pd.read_csv(dirN+filen).dropna()
	    cxs = {'AS':1.,'AZ':1.,'AU':1.,'HS':1.,'HU':1.,'INPUT':50000,'U3M':1.,'U3N':1.,'UH':1.}
	    cxs,lamda_counted=add_lamdas(cxs,tmp,df)
	    if not lamda_counted:
		tmp = [dirN.replace("/",'').split("_")[-1],tmp[3],tmp[4],tmp[-2],'0']
		cxs,lamda_counted = add_lamdas(cxs,tmp,df)

            nics=float(tmp[-2])
            if len(df.values[:])!=nics:
		print len(df.values[:]),'!=',nics
            res=getStates(df)
	    ##fileo.write("AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,E,EM,M,W,WO,O,M/W,M/WO,M/O,EM/W,EM/WO,EM/O,E/W,E/WO,E/O,nics\n")
	    fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(cxs['AS'],cxs['AZ'],cxs['AU'],cxs['HS'],cxs['HU'],cxs['INPUT'],cxs['U3M'],cxs['U3N'],cxs['UH'],res['E'],res['EM'],res['M'],res['W'],res['WO'],res['O'],res['M/W'],res['M/WO'],res['M/O'],res['EM/W'],res['EM/WO'],res['EM/O'],res['E/W'],res['E/WO'],res['E/O'],nics))

    fileo.close()
