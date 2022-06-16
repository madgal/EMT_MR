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

mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams.update({'font.size': 20})

xlabel,title,fsave='','','tmp'

def plotRESvlambda_singleLink(fileN,xlabel,title,fsave,typeV):
    df = pd.read_csv(fileN)
    print 'E',np.unique(df['E'])
    print 'EM',np.unique(df['EM'])
    print 'M',np.unique(df['M'])
    xval=df[xlabel.upper()].values
    
    if typeV=='EMT':
        cols={'E':'r','EM':'k','M':'b'}
        compV={'E':1,'EM':2,'M':3}
	keyList=['u','mz','Z','ms','u3','S']
	sub='EMT'
    else:
        cols={'O':'r','WO':'k','W':'b'}
        compV={'O':1,'WO':2,'W':3}
	keyList=['A','Rmt','Rnox','h','mh']
	sub='MR'
    for el in compV.keys():
        check=[]
        fig,axs = plt.subplots(2,3,figsize=(15,6))
        row,col=0,0
	print el,compV[el]
	inds = np.argwhere(df[el].values==1)[:,0]
	print len(inds),len(df.values)
	yval = compV[el]*np.ones(len(inds))
        for key in keyList:
            if el not in check and row==0 and col==0:
                axs[row,col].plot(xval[inds],df[key].values[inds],color=cols[el],marker='*',label=el,linestyle='')
                check+=[el]
            else:
                axs[row,col].plot(xval[inds],df[key].values[inds],color=cols[el],marker='*',linestyle='')
            axs[row,col].set_xlabel("$\lambda_{"+xlabel+"}$")
            axs[row,col].set_ylabel(key)
            col+=1
            if col==3:
                row+=1
                col=0
        plt.suptitle(xlabel)
        axs[0,0].legend()
        fsave='figures/Expressions_'+xlabel+'_'+el+'_'+sub+'.png'
        fig.savefig(fsave,bbox_inches='tight')
	plt.close()

    ### put on top now
    fig,axs = plt.subplots(2,3,figsize=(15,6))
    check=[]
    for el in compV.keys():
        row,col=0,0
	inds = np.argwhere(df[el].values==1)[:,0]
	print len(inds),len(df.values)
	yval = compV[el]*np.ones(len(inds))
        for key in keyList:
            if el not in check and row==0 and col==0:
                axs[row,col].plot(xval[inds],df[key].values[inds],color=cols[el],marker='*',label=el,linestyle='')
                check+=[el]
            else:
                axs[row,col].plot(xval[inds],df[key].values[inds],color=cols[el],marker='*',linestyle='')
            axs[row,col].set_xlabel("$\lambda_{"+xlabel+"}$")
            axs[row,col].set_ylabel(key)
            col+=1
            if col==3:
                row+=1
                col=0
    plt.suptitle(xlabel)
    axs[0,0].legend()
    fsave='figures/Expressions_'+xlabel+'_'+sub+'_comp.png'
    fig.savefig(fsave,bbox_inches='tight')
    plt.close()



dirN='data/vals_crosstalk_Au.txt'
xlabel='AU'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'EMT')

dirN='data/vals_crosstalk_iHu.txt'
xlabel='HU'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'EMT')

dirN='data/vals_crosstalk_AZ.txt'
xlabel='AZ'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'EMT')

dirN='data/vals_crosstalk_AS.txt'
xlabel='AS'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'EMT')

dirN='data/vals_crosstalk_HS.txt'
xlabel='HS'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'EMT')

dirN='data/vals_crosstalk_u3m.txt'
xlabel='U3M'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'MR')
    
dirN='data/vals_crosstalk_u3n.txt'
xlabel='U3N'
plotRESvlambda_singleLink(dirN,xlabel,title,fsave,'MR')

