import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.lines import Line2D
import os
import math

#mpl.rcParams['xtick.labelsize'] = 20
#mpl.rcParams['ytick.labelsize'] = 20
#mpl.rcParams.update({'font.size': 20})
def plotHeatmap(x,y,nbin,xlabel='',ylabel='',title='',hold=False):
    [HH,xh,yh] = np.histogram2d(x,y,bins=(nbin,nbin))
    HH = HH.T/(1.)
    xmin,xmax=xh[0],xh[-1]
    ymin,ymax=yh[0],yh[-1]
    plt.imshow(HH,interpolation='nearest',origin='low',extent=[xmin,xmax,ymin,ymax],aspect='auto',cmap='gnuplot2_r')#,cmap='jet',vmin=0.,vmax=100,aspect='auto')
    plt.colorbar(label="Number of Initial Conditions")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel) 
    plt.title(title)   
    if not hold:
        plt.show()
    plt.close()
    return [HH,xh,yh]

def stateThresholds():
    return {'Zeb mRNA':[ 62+(301-62)/2.,301+(990-301)/2.],
            '$\\mu_{200}$':[1265+(12389-1265)/2.,12389+(19098-12389)/2.],
            'Hif-1':[26+(292-26)/2.,292+(490-292)/2.],
            'AMPK':[69+(327-69)/2.,327+(436-327)/2.]}
def getStates_fromMap(df_res,hold=False):

    x2 = df_res['A']
    y2 = df_res['h']
    x1 = df_res['u']
    y1 = df_res['mz']
    hm1= plotHeatmap(x1,y1,10,'$\mu_{200}$','Zeb mRNA','Steady States',hold)
    hm2 =plotHeatmap(x2,y2,10,'AMPK','Hif-1','Steady States',hold)
    res2,res1={'n':[],'u':[],'mz':[]},{'n':[],'A':[],'h':[]}
    
    thresh = stateThresholds()
            
    E = (df_res['u']>thresh['$\mu_{200}$'][1])&(df_res['mz']<=thresh['Zeb mRNA'][0])
    M = (df_res['u']<=thresh['$\mu_{200}$'][0])&(df_res['mz']>thresh['Zeb mRNA'][1])
    EM =(df_res['u']>thresh['$\mu_{200}$'][0])&(df_res['u']<=thresh['$\mu_{200}$'][1])&(df_res['mz']>thresh['Zeb mRNA'][0])&(df_res['mz']<=thresh['Zeb mRNA'][1])
    
    
    O = (df_res['A']>thresh['AMPK'][1])&(df_res['h']<=thresh['Hif-1'][0])
    W = (df_res['A']<=thresh['AMPK'][0])&(df_res['h']>thresh['Hif-1'][1])
    WO =(df_res['A']>thresh['AMPK'][0])&(df_res['A']<=thresh['AMPK'][1])&(df_res['h']>thresh['Hif-1'][0])&(df_res['h']<=thresh['Hif-1'][1])
    
    
    #print E,M,EM,E+EM+M,np.sum(df_res['u']>=-11)
    #print O,W,WO,O+W+WO
    
    tmp = np.sum(E)+np.sum(EM)+np.sum(M)+np.sum(O)+np.sum(WO)+np.sum(W)
    resultsS = {'name':['E','EM','M','O','WO','W','unknown'],
                'amount':[np.sum(E),np.sum(EM),np.sum(M),np.sum(O),np.sum(WO),np.sum(W),len(df_res['Z'])-tmp]}
    
    tmp = np.sum(E&O)+np.sum(E&WO)+np.sum(E&W)+np.sum(EM&O)+np.sum(EM&WO)+np.sum(EM&W)+np.sum(M&O)+np.sum(M&WO)+np.sum(M&W)
    resultsE = {}
    resultsE = {'name':['E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W','unknown'],
                'amount':[np.sum(E&O),np.sum(E&WO),np.sum(E&W),np.sum(EM&O),np.sum(EM&WO),np.sum(EM&W),np.sum(M&O),np.sum(M&WO),np.sum(M&W),len(df_res['Z'])-tmp]}
    maxV = np.max(resultsE['amount'])
    #print np.sum(resultsE['amount']),len(df_res['Z'])
    return resultsE,resultsS,maxV
def getFileSet2(direct,sets):
    lists={'x':[],'y':[],'c':[],'l':[]}
    count=0
    for filen in os.listdir(direct):
        if "res.txt" in filen:
            #print filen
            #plotRes(direct+filen)
            name= filen.split(".")[0].split("_")[2]
            lamda=''
            while name[-1].isdigit():
                lamda=name[-1]+lamda
                name = name[:-1]
            if ("AZ" in filen) or ("AS" in filen) or ("IHu" in filen)  :
                scale=0.1
            else:
                scale=1.
            lists['x']+=[count]
            lists['y']+=[filen]
            lists['c']+=[0]
            lists['l']+=[[int(lamda[sets[0]:sets[1]]),int(lamda[sets[1]:])]]
            count+=1
    return lists
def getMeshGrid(final):
    finalX={}
    for key in final:
        xval = np.array(final[key]['x'])
        yval = np.array(final[key]['y'])
        zval = np.array(final[key]['c'])

        tmpX = np.unique(xval)
        tmpY = np.unique(yval)
        x2,y2,z2=[],[],[]

        for val in tmpX:
            tmp1,tmp2,tmp3=[],[],[]
            for val2 in tmpY:
                inds = np.argwhere(xval==val)[:,0]
                ind2 = np.argwhere(yval==val2)[:,0]
                ind = np.intersect1d(inds,ind2)
                tmp1+=[float(val)]
                tmp2+=[float(val2)]
                tmp3+=[float(zval[ind])]
            x2+=[tmp1]    
            y2+=[tmp2]
            z2+=[tmp3]    

	finalX[key]={}
        finalX[key]['x'] = np.array(x2)
        finalX[key]['y'] = np.array(y2)
        finalX[key]['c'] = np.array(z2).reshape(finalX[key]['y'].shape)

    return finalX

def plotLoops(dirc,results,xlab,ylab,title):#,legL,labelAll=False):
    
    ogFile="coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt"
    nics=1000.

    df_res =pd.read_csv(ogFile).dropna()
    nocross,noCs,tmp  =getStates_fromMap(df_res,hold=True)
    

    noCRes={}
    for i in range(len(nocross['amount'])):
        key = nocross['name'][i]
        val = nocross['amount'][i]
	noCRes[key]=val
    for i in range(len(noCs['amount'])):
        key = noCs['name'][i]
        val = noCs['amount'][i]
	noCRes[key]=val

    final={'E':{'x':[],'y':[],'c':[]},'EM':{'x':[],'y':[],'c':[]},'M':{'x':[],'y':[],'c':[]},'O':{'x':[],'y':[],'c':[]},'WO':{'x':[],'y':[],'c':[]},'W':{'x':[],'y':[],'c':[]},'E/O':{'x':[],'y':[],'c':[]},'E/WO':{'x':[],'y':[],'c':[]},'E/W':{'x':[],'y':[],'c':[]},'EM/O':{'x':[],'y':[],'c':[]},'EM/WO':{'x':[],'y':[],'c':[]},'EM/W':{'x':[],'y':[],'c':[]},'M/O':{'x':[],'y':[],'c':[]},'M/WO':{'x':[],'y':[],'c':[]},'M/W':{'x':[],'y':[],'c':[]}}
        
    for i in range(len(results['y'])):
        try:
            df = pd.read_csv(dirc+results['y'][i]).dropna()
            mapResE,resSum,maxC =getStates_fromMap(df,hold=True)
            tmp={}
            for j in range(len(mapResE['name'])):
                tmp[mapResE['name'][j]]=mapResE['amount'][j]/nics*100.
            for j in range(len(resSum['name'])):
                tmp[resSum['name'][j]]=resSum['amount'][j]/nics*100.

            results['y'][i] = tmp
        except:
            results['y'][i]={'E':-100,'EM':-100,'M':-100,'O':-100,'WO':-100,'W':-100,'E/O':-100,'E/WO':-100,'E/W':-100,'EM/O':-100,'EM/WO':-100,'EM/W':-100,'M/O':-100,'M/WO':-100,'M/W':-100}

	print results['y'][i]
	print noCRes

	for key in results['y'][i]:
	    if key!='unknown':
		final[key]['x']+=[results['l'][i][0]]
		final[key]['y']+=[results['l'][i][1]]
		if results['y'][i][key]<noCRes[key]:
			final[key]['c']+=[0]#['r']
		elif results['y'][i][key]>noCRes[key]:
			final[key]['c']+=[100]#['r']
		else:
			final[key]['c']+=[50]#['k']


    finalX = getMeshGrid(final)
    plotCoupledLoop(finalX,final,xlab,ylab,title)

def plotCoupledLoop(results2,results,xlab,ylab,title):#,legL,labelAll):

    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(30,30))
    plt.subplots_adjust(hspace=0.5,right=0.85)
    my_cmap= cm.get_cmap('cividis')
    my_norm=Normalize(vmin=0,vmax=np.max(results['E/O']['c']))

    colors=my_cmap(my_norm(results['E/O']['c']))

    for i in range(len(results['E/O']['x'])):
        ax1.plot(results['E/O']  ['x'][i],results['E/O']  ['y'][i],'.',markersize=20,color=colors[i])
        ax2.plot(results['E/WO'] ['x'][i],results['E/WO'] ['y'][i],'.',markersize=20,color=colors[i])
        ax3.plot(results['E/W']  ['x'][i],results['E/W']  ['y'][i],'.',markersize=20,color=colors[i])
        ax4.plot(results['EM/O'] ['x'][i],results['EM/O'] ['y'][i],'.',markersize=20,color=colors[i])
        ax5.plot(results['EM/WO']['x'][i],results['EM/WO']['y'][i],'.',markersize=20,color=colors[i])
        ax6.plot(results['EM/W'] ['x'][i],results['EM/W'] ['y'][i],'.',markersize=20,color=colors[i])
        ax7.plot(results['M/O']  ['x'][i],results['M/O']  ['y'][i],'.',markersize=20,color=colors[i])
        ax8.plot(results['M/WO'] ['x'][i],results['M/WO'] ['y'][i],'.',markersize=20,color=colors[i])
        ax9.plot(results['M/W']  ['x'][i],results['M/W']  ['y'][i],'.',markersize=20,color=colors[i])

    #ax1.contourf(results['E/O'] ['x'],results['E/O'] ['y'],results['E/O'] ['c'])#,cmap=my_cmap)
    #ax2.contourf(results['E/WO']['x'],results['E/WO']['y'],results['E/WO']['c'])#,cmap=my_cmap)
    #ax3.contourf(results['E/W'] ['x'],results['E/W'] ['y'],results['E/W'] ['c'])#,cmap=my_cmap)
    #ax4.contourf(results['EM/O']['x'],results['EM/O']['y'],results['EM/O']['c'])#,cmap=my_cmap)
    #ax5.contourf(results['EM/WO']['x'],results['EM/WO']['y'],results['EM/WO']['c'])#,cmap=my_cmap)
    #ax6.contourf(results['EM/W']['x'],results['EM/W']['y'],results['EM/W']['c'])#,cmap=my_cmap)
    #ax7.contourf(results['M/O'] ['x'],results['M/O'] ['y'],results['M/O'] ['c'])#,cmap=my_cmap)
    #ax8.contourf(results['M/WO']['x'],results['M/WO']['y'],results['M/WO']['c'])#,cmap=my_cmap)
    #ax9.contourf(results['M/W'] ['x'],results['M/W'] ['y'],results['M/W'] ['c'])#,cmap=my_cmap)

    ax1.set_xlabel(xlab)
    ax2.set_xlabel(xlab)
    ax3.set_xlabel(xlab)
    ax4.set_xlabel(xlab)
    ax5.set_xlabel(xlab)
    ax6.set_xlabel(xlab)
    ax7.set_xlabel(xlab)
    ax8.set_xlabel(xlab)
    ax9.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax2.set_ylabel(ylab)
    ax3.set_ylabel(ylab)
    ax4.set_ylabel(ylab)
    ax5.set_ylabel(ylab)
    ax6.set_ylabel(ylab)
    ax7.set_ylabel(ylab)
    ax8.set_ylabel(ylab)
    ax9.set_ylabel(ylab)
    ax1.set_title("E/O")
    ax2.set_title("E/WO")
    ax3.set_title("E/W")
    ax4.set_title("EM/O")
    ax5.set_title("EM/WO")
    ax6.set_title("EM/W")
    ax7.set_title("M/O")
    ax8.set_title("M/WO")
    ax9.set_title("M/W")


    legend_elements = [ Line2D([0], [0], marker='.', color='w', label='Higher', markerfacecolor=my_cmap(my_norm(100)), markersize=20),
			 Line2D([0], [0], marker='.', color='w', label='Equal', markerfacecolor=my_cmap(my_norm(50)), markersize=20),
			 Line2D([0], [0], marker='.', color='w', label='Lower', markerfacecolor=my_cmap(my_norm(0)), markersize=20) ]                   

    fig.legend(handles=legend_elements,loc="center right",borderaxespad=0.1,title="Legend")
    fig.savefig(title+"_coupledStates.png",bbox_inches='tight')
    plt.show()
    plt.close()
#####################
#####################
#####################


res_hS_aS = getFileSet2("coupledWReg_Ccode/crosstalk_hS_aS/",[0,1])
plotLoops("coupledWReg_Ccode/crosstalk_hS_aS/",res_hS_aS,xlab='Foldchange HS',ylab='Foldchange AS',title='HS_AS')#,legL='$\lambda_{AS}$=')

res_mrNr = getFileSet2("coupledWReg_Ccode/crosstalk_mrNr/",[0,1])
print res_mrNr
plotLoops("coupledWReg_Ccode/crosstalk_mrNr/",res_mrNr,xlab='Foldchange $\mu_{34}$Rmt',ylab='Foldchange $\mu_{34}$Rnox',title='mrNR')

res_uh_hu = getFileSet2("coupledWReg_Ccode/crosstalk_uh_hu/",[0,2])
print res_uh_hu
plotLoops("coupledWReg_Ccode/crosstalk_uh_hu/",res_uh_hu,xlab='Par. set (L,Y)',ylab='Foldchange Hu',title='uh_hu')

res_uh_hu_mrNr = getFileSet2("coupledWReg_Ccode/crosstalk_uh_hu_mrNr/",[2,3])
print res_uh_hu_mrNr
plotLoops("coupledWReg_Ccode/crosstalk_uh_hu_mrNr/",res_uh_hu_mrNr,xlab='Par. set',ylab='Foldchange HS',title='uh_hu_mrnr')

res_uh_hu_au = getFileSet2("coupledWReg_Ccode/crosstalk_uh_hu_au/",[2,3])
plotLoops("coupledWReg_Ccode/crosstalk_uh_hu_au/",res_uh_hu_au,xlab='Foldchange Hu',ylab='Foldchange Au',title='uh_hu_au')

res_hS_aS_AZ = getFileSet2("coupledWReg_Ccode/crosstalk_hS_aS_AZ/",[0,1])
plotLoops("coupledWReg_Ccode/crosstalk_hS_aS_AZ/",res_hS_aS_AZ,xlab='Par. set',ylab='Foldchange HS',title='HS_AS_AZ')
