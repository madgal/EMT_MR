import matplotlib
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib import cm

from aux_func_States import getLegend
from aux_func_States import equals
from aux_func_States import getStates
from aux_func_States import getCompHybridVal


labelD={'W':'W','O':'O','WO':'W/O',
        'E':'E','EM':'E/M','M':'M',
        'E/O':'E-O','E/W':'E-W','E/WO':'E-W/O',
        'EM/O':'E/M-O','EM/W':'E/M-W','EM/WO':'E/M-W/O',
        'M/O':'M-O','M/W':'M-W','M/WO':'M-W/O'}

alpha={'W':0.5,'E/W':1.,'EM/W':1.,'M/W':1.,'WO':0.5,'E/WO':1.,'EM/WO':1.,'M/WO':1.,'O':0.5,'E/O':1.,'EM/O':1.,'M/O':1.,'E':0.5,'M':0.5,'EM':0.5}


color={'W':'y','O':'b','WO':'g','EM':'y','M':'b','E':'g'}
mark = {'E':'dotted','M':'--','EM':'-.','WO':'dotted','W':'--','O':'-.'}


def getColor(label):
    
    ### Reduced EMT states
    if equals(label,['E']):
        return 'tomato'
    if equals(label,['M']):
        return 'cornflowerblue'
    if equals(label,['EM']):
        return 'cyan'
    if equals(label,['E','M']):
        return 'plum'
    if equals(label,['E','EM']):
        return 'orange'
    if equals(label,['M','EM']):
        return 'lightgreen'
    if equals(label,['E','M','EM']):
        return 'yellow'
    
    ### REduced to Metabolic states
    if equals(label,['W']):
        return 'tomato'
    if equals(label,['O']):
        return 'cornflowerblue'
    if equals(label,['WO']):
        return 'cyan'
    if equals(label,['W','O']):
        return 'plum'
    if equals(label,['W','WO']):
        return 'orange'
    if equals(label,['O','WO']):
        return 'lightgreen'
    if equals(label,['O','W','WO']):
        return 'yellow'
    
    ### make the all 9 coupled states always be black
    if equals(label,['E/O','E/WO','EM/O','M/O','EM/WO','M/WO','M/W','EM/W','E/W']):
        return 'yellow'
    
    
    ##generate the color dictionary based on the keys so that you can easily add results (511 total possible combos)
    keys=[  ['E/W'],['M/W'],['EM/W'],
            ['E/O'],['M/O'],['EM/O'],
            ['E/WO'],['M/WO'],['EM/WO'],

            ['E/O','EM/WO'],
            ['EM/O','EM/WO'],
            ['EM/O','M/O'],
            ['E/O','E/WO'],
            ['E/O','E/W'],
            ['M/O','M/W'],
            ['M/WO','M/W'],
            ['M/WO','M/O'],
            ['EM/O','M/WO'],
            ['E/O','EM/O'],
            ['E/O','M/O'],
            ['M/O','EM/WO'],
            ['E/WO','EM/WO'],

            ['E/O','E/W','E/WO'],
            ['E/O','E/WO','EM/W'],
            ['E/O','E/WO','EM/WO'],
            ['E/O','E/WO','M/W'],
            ['E/O','EM/W','EM/WO'],
            ['E/O','M/W','EM/WO'],
            ['E/O','M/W','M/WO'],
            ['E/O','M/O','M/W'],
            ['M/W','M/WO','M/O'],
            ['E/O','EM/O','M/O'],
            ['EM/O','M/O','M/W'],
            ['EM/O','EM/W','EM/WO'],
            ['EM/WO','M/O','M/W'],
            ['E/WO','M/O','M/W'],
            ['EM/O','M/WO','EM/WO'],
            ['EM/O','M/O','EM/WO'],

            ['E/O','EM/WO','M/W','M/WO'],
            ['E/O','EM/WO','M/W','M/O'],
            ['E/W','EM/WO','M/O','M/WO'],
            ['E/WO','EM/WO','M/O','M/WO'],
            ['E/WO','EM/WO','M/O','M/W'],
            ['E/W','EM/W','M/O','M/WO'],
            ['E/W','E/WO','E/O','M/W'],
            ['EM/O','M/W','M/WO','M/O'],
            ['E/O','M/W','M/WO','M/O'],
            ['E/O','EM/O','M/O','M/W'],
            ['E/O','E/WO','EM/W','M/W'],
            ['E/O','E/WO','EM/WO','M/W'],
            ['E/WO','EM/W','M/O','M/WO'],
            ['E/O','E/W','M/O','M/W'],
            ['E/O','E/WO','M/O','M/W'],
            ['E/O','E/WO','M/WO','M/W'],
            ['E/O','EM/W','M/O','M/W'],
            ['E/O','E/WO','EM/O','M/O'],
            ['E/WO','EM/O','EM/WO','M/O'],

            ['E/W','EM/W','M/O','M/W','M/WO'],
            ['E/WO','EM/O','EM/WO','M/O','M/WO'],
            ['E/W','EM/WO','M/W','M/O','M/WO'],
            ['E/WO','EM/WO','M/W','M/O','M/WO'],
            ['E/O','E/W','E/WO','EM/W','M/W'],
            ['E/O','E/W','E/WO','M/WO','M/W'],
            ['E/O','E/W','EM/W','M/O','M/W'],
            ['E/O','EM/O','M/O','M/WO','M/W'],
            ['EM/WO','EM/O','M/O','M/WO','M/W'],
            ['E/WO','EM/W','M/O','M/W','M/WO'],
            ['E/WO','E/O','M/O','M/W','M/WO'],
            ['EM/WO','E/O','M/O','M/W','M/WO'],
            ['E/O','E/WO','EM/WO','M/O','M/W'],
            ['E/O','E/WO','EM/WO','M/WO','M/W'],
            ['E/O','E/WO','EM/W','M/WO','M/W'],

            ['E/O','E/WO','EM/O','EM/WO','M/O','M/WO'],
            ['E/O','E/W','EM/O','EM/W','M/O','M/W'],
            ['E/WO','EM/O','EM/WO','M/O','M/W','M/WO'],
            ['E/O','E/W','E/WO','EM/W','M/WO','M/W'],
            ['E/O','E/W','E/WO','M/O','M/WO','M/W'],
            ['E/O','EM/W','E/WO','M/O','M/WO','M/W'],
            ['E/O','EM/O','EM/WO','M/O','M/W','M/WO'],
            ['E/O','E/WO','EM/WO','M/O','M/W','M/WO'],
            ['EM/W','EM/O','EM/WO','M/O','M/W','M/WO'],
            ['E/O','EM/O','EM/W','M/O','M/W','M/WO'],

            ['E/O','E/WO','EM/O','EM/WO','M/WO','M/O','M/W'],
            ['E/O','EM/WO','EM/O','EM/W','M/WO','M/O','M/W'],
            ['E/O','E/WO','EM/W','EM/WO','M/WO','M/O','M/W'],
            ['E/O','E/WO','E/W','EM/W','M/WO','M/O','M/W'],

            ['E/O','E/WO','E/W','EM/W','EM/WO','M/W','M/O','M/WO'],
            ['E/O','E/WO','EM/O','EM/W','EM/WO','M/O','M/WO','M/W'],
             ]
    
    __clist=[]
    cmap1= cm.get_cmap('tab20b')
    cmap2= cm.get_cmap('tab20c')
    count,base=0,0
    for i in range(20):
        if matplotlib.colors.to_hex(cmap1(base+count*4))=='yellow':
            count+=1
        else:
            __clist+=[matplotlib.colors.to_hex(cmap1(base+count*4))]
            __clist+=[matplotlib.colors.to_hex(cmap2(base+count*4))]
            count+=1
        if count==5:
                count=0
                base+=1
    __clist+=['darkorange','beige','chartreuse','cyan','lime','magenta','darkviolet',
              'blue','red','saddlebrown','tomato','pink','salmon',
              'olive','wheat','mediumslateblue','hotpink',
              'seagreen','darkkhaki','thistle',
              'lightgreen','lightsteelblue','pink','mistyrose',
              'silver','rosybrown','coral','lightyellow','darkseagreen','paleturquoise',
              'mediumvioletred','darkviolet','crimson','cornflowerblue','dodgerblue','springgreen','lime','pink']
    for i in range(20):
        if matplotlib.colors.to_hex(cmap1(base+count*4))=='yellow':
            count+=1
        else:
            __clist+=[matplotlib.colors.to_hex(cmap1(base+count*4))]
            __clist+=[matplotlib.colors.to_hex(cmap2(base+count*4))]
            count+=1
        if count==5:
                count=0
                base+=1
                
    cmap = matplotlib.cm.get_cmap('hsv')
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=len(keys))

    fileout = open('colorsUsed.txt','a')
    for i in range(len(keys)):
        if equals(label,keys[i]):
            if i>=len(__clist):
                print "NEED MORE COLORS",len(keys),len(__clist)
                
	    fileout.write("%s\n" %i)
            return __clist[i%len(__clist)]
            #return cmap(norm(i))


##############3
##############3
##############3
def colorMap(fileN,reduced=None):
    df = pd.read_csv(fileN)
    lab={}
    if 'col' in df.columns:#reduced:
        keys=['E','EM','M','W','O','WO']
        mainKey='col'
    else:
        keys=['E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']
        mainKey='t3'
    for i in range(len(df.values)):
        tmp=[]
        for key in keys:
            if key in df.keys():
                if df[key].values[i]==1:
                    tmp+=[key]
        lab[df[mainKey].values[i]] =tmp
        
    return lab
##############3
##############3
##############3
def plotPhases(axa,subF,titleStart,xlabel=None,ylabel=None,xlim=None,ylim=None,legend=None,switchX=None,reduced=None,fs=20,noEMWO=None,legendLoc=None):
    axa.text(-0.2, 1.07, subF,transform=axa.transAxes,verticalalignment='top', horizontalalignment='right',color='black', fontsize=fs+5)
    labs = colorMap(titleStart+"_legend.txt",reduced)
    if not switchX:
        xk = 'X'
        yk='Y'
    else:
        xk='Y'
        yk='X'
    df = pd.read_csv(titleStart+".txt")
    for i in range(len(df)):
            axa.fill_between([df[xk+'1'].values[i],df[xk+'2'].values[i]],df[yk+'1'].values[i],df[yk+'2'].values[i],facecolor=getColor(labs[df['color'].values[i]]),linewidth=0.0)

    try:
        if not noEMWO:
            dfHE = pd.read_csv(titleStart+"_HHexists.txt")
            for i in range(len(dfHE)):
                axa.fill_between([dfHE[xk+'1'].values[i],dfHE[xk+'2'].values[i]],dfHE[yk+'1'].values[i],dfHE[yk+'2'].values[i],facecolor='none',hatch='.',edgecolor='k',linewidth=0.0)
            dfHO = pd.read_csv(titleStart+"_HHOnlhy.txt")
            for i in range(len(dfHO)):
                axa.fill_between([dfHO[xk+'1'].values[i],dfHO[xk+'2'].values[i]],dfHO[yk+'1'].values[i],dfHO[yk+'2'].values[i],facecolor='none',hatch='.',edgecolor='r',linewidth=0.0)
    except:
          print "noEM"

    if legend:
	retL,retC=[],[]
        dfl = pd.read_csv(titleStart+"_legend.txt")
        Lflag=False
        try:
            Lflag=True
            legend_elements=[]
            for i in range(len(dfl)):
                tmp=[]
                for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
                    if dfl[key].values[i]==1:
                        tmp+=[key]
                legend_elements+=[ Patch(facecolor=getColor(labs[dfl['t3'].values[i]]),label=tmp)]
		retL+=[tmp]
		retC+=[getColor(labs[dfl['t3'].values[i]])]
        except:
            Lflag=False
            print "not coupled"
        try:
            if not Lflag:
                Lflag=True
                legend_elements=[]
                for i in range(len(dfl)):
                    tmp=[]
                    for key in ['E','M','EM']:
                        if dfl[key].values[i]==1:
                            tmp+=[key]
                    legend_elements+=[ Patch(facecolor=getColor(labs[dfl['col'].values[i]]),label=tmp)]
                    retL+=[tmp]
                    retC+=[getColor(labs[dfl['col'].values[i]])]
        except:
            Lflag=False
            print "not coupled"
        try:
            if not Lflag:
                Lflag=True
                legend_elements=[]
                for i in range(len(dfl)):
                    tmp=[]
                    for key in ['W','WO','O']:
                        if dfl[key].values[i]==1:
                            tmp+=[key]
                    legend_elements+=[ Patch(facecolor=getColor(labs[dfl['col'].values[i]]),label=tmp)]
                    retL+=[tmp]
                    retC+=[getColor(labs[dfl['col'].values[i]])]
        except:
            print "not coupled"
        if not legendLoc and legend!='return':
        	axa.legend(handles=legend_elements, bbox_to_anchor=(5.1, 1.))
	elif legend!='return':
        	axa.legend(handles=legend_elements, bbox_to_anchor=legendLoc)
    if xlabel:
        if switchX:
            axa.set_xlabel(ylabel)
        else:
            axa.set_xlabel(xlabel)
    if ylabel:
        if switchX:
            axa.set_ylabel(xlabel)
        else:
            axa.set_ylabel(ylabel)
    if xlim:
        if switchX:
            axa.set_xlim(ylim)
        else:
            axa.set_xlim(xlim)
    if ylim:
        if switchX:
            axa.set_ylim(xlim)
        else:
            axa.set_ylim(ylim)
    if legend=='return':
        return retC,retL

##############3
##############3
##############3
def plotEMWOheatmap(axa,subF,titleStart,xlabel=None,ylabel=None,xlim=None,ylim=None,switchX=None,fs=None):
    if not fs:
	fs=20
    axa.text(-0.1, 1.05, subF,transform=axa.transAxes,verticalalignment='top', horizontalalignment='right',color='black', fontsize=fs+5)
    df = pd.read_csv(titleStart+"_EMWO.txt")
    if not switchX:
        xk = 'X'
        yk='Y'
    else:
        xk='Y'
        yk='X'

    compHybrid=getCompHybridVal(compI=False)
    plt.scatter(df[xk],df[yk],c=df['EMWO'],cmap='bwr',norm=MidpointNormalize(midpoint=compHybrid,vmin=0,vmax=100))
    plt.colorbar(label='E/M-W/O (%)')
    if xlabel:
        if switchX:
            axa.set_xlabel(ylabel)
        else:
            axa.set_xlabel(xlabel)
    if ylabel:
        if switchX:
            axa.set_ylabel(xlabel)
        else:
            axa.set_ylabel(ylabel)
    if xlim:
        if switchX:
            axa.set_xlim(ylim)
        else:
            axa.set_xlim(xlim)
    if ylim:
        if switchX:
            axa.set_ylim(xlim)
        else:
            axa.set_ylim(ylim)

##############3
##############3
##############3
def plotPhaseICS(ax,subF,titleStart,typeN,xlabel=None,xlim=None,ylabel=None,noICS=False,fs=None):
    ##########
    if not fs:
        fs=20
        
    ax.text(-0.15, 1.15, subF,transform=ax.transAxes,verticalalignment='top', horizontalalignment='right',color='black', fontsize=fs+5)

    labs = colorMap(titleStart+"_legend.txt",reduced=True)
    df = pd.read_csv(titleStart+".txt")
    for i in range(len(df)):
            ax.fill_between([df['X1'].values[i],df['X2'].values[i]],-10,110,facecolor=getColor(labs[df['color'].values[i]]),linewidth=0.0)

    
    if typeN=='MR':
        keyList=['O','WO','W']
    else:
        keyList=['M','EM','E']
    ax.set_ylim(-1,101)
    if not noICS:
        mark = {'O':'dotted','W':'--','WO':'-.','M':'dotted','E':'--','EM':'-.'}
        df2 = pd.read_csv(titleStart+"_ics.txt")
        for key in keyList:
            inds = np.argsort(df2['x'].values)
            ax.plot(df2['x'].values[inds],df2[key].values[inds],linestyle=mark[key],linewidth=5,color='k',label=labelD[key])
        ax.legend(loc='center left',frameon=False,bbox_to_anchor=(1.05,0.5),fontsize=fs)
        ax.set_ylabel("Initial\nConditions (%)")
        if xlim:
            if xlim==-1:
                ax.set_xlim(0,np.max(df2['x'].values))
            else:
                ax.set_xlim(xlim)
    else:
        ax.set_yticks([])
        if xlim==-1:
                ax.set_xlim(0,np.max(df['X2'].values))
        else:
                ax.set_xlim(xlim)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("Initial\nConditions (%)")

##############3
##############3
##############3
def getCompV(el):
    df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    mapRes = getStates(df)
    return mapRes[el]/10.
##############3
##############3
##############3
def plotICS(ax,df,key,key1,key2,fs=None):
    if not fs:
	fs=20
    xvals=df['icsNumber'].values
    inds = np.argsort(xvals)
    yv=df[key+"Avg"][inds]+df[key1+"Avg"][inds]+df[key2+"Avg"][inds]
    yv1=df[key2+"Avg"][inds]+df[key1+"Avg"][inds]
    yv2=df[key2+"Avg"][inds]
    ye=df[key+"Std"][inds]
    ye1=df[key1+"Std"][inds]
    ye2=df[key2+"Std"][inds]
    ax.bar(np.arange(0,len(df)),yv,label=labelD[key])#,yerr=ye,label=key)
    ax.bar(np.arange(0,len(df)),yv1,label=labelD[key1])#,yerr=ye1,label=key1)
    ax.bar(np.arange(0,len(df)),yv2,label=labelD[key2])#,yerr=ye2,label=key2)
    ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left',frameon=False,fontsize=fs)
    ax.set_xlabel("Initial conditions")
    ax.set_ylabel("Steady states (%)")
    xvals=xvals[inds]
    ax.set_xticks(np.arange(0,len(df)))
    ax.set_xticklabels(xvals,rotation=70)

##############3
##############3
##############3
def plotICS_coupled(gs1,df,colLength):
    row,col=0,0
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
        ax = plt.subplot(gs1[row,col])
        inds = np.argsort(df['icsNumber'].values)
        ax.bar(np.arange(len(df['icsNumber'])),df[key+"Avg"].values[inds],yerr=df[key+"Std"].values[inds])
        ax.set_title(key)
        ax.set_ylabel("Steady States (%)")
        ax.set_xlabel("Initial conditions")
        ax.set_xticks(np.arange(len(df['icsNumber'])))
        ax.set_xticklabels(df['icsNumber'].values[inds],rotation=30)
        col+=2
        if col>=colLength+1:
            row+=2
            col=0 

##############3
##############3
##############3
def plotBreakdown(ax,df,masterKey,fs=None):
    if not fs:
	fs=20
    names = {'E':'Epithelial','EM':'E/M','M':'Mesenchymal','W':"Warburg",'WO':'W/O','O':'OXPHOS'}
    if masterKey in ['W','WO','O']:
        keyList=['E/'+masterKey,'EM/'+masterKey,'M/'+masterKey]
        title = "Breakdown of coupled states with\n"+names[masterKey]+" metabolic phenotype" 
        keyLoc=1
    else:
        keyList=[masterKey+"/O",masterKey+"/WO",masterKey+"/W"]
        title = "Breakdown of coupled states with\n"+names[masterKey]+" phenotype"
        keyLoc=1
    for key in keyList:#['E/O','EM/O','M/O','O']:#['E/O','E/WO','E/W','M/O','M/WO','M/W','EM/O','EM/WO','EM/W']:
            x=df['x'].values
            y=df[key].values
            inds = np.argsort(x)
            x=x[inds]
            y=y[inds]
            ax.plot(x,np.array(y)/getCompV(key),color=color[masterKey],linestyle=mark[key.split("/")[keyLoc]],label=labelD[key],lw=6,markersize=20,alpha=alpha[key])
    ax.set_xlabel("$\lambda_{}$")
    ax.set_ylabel("Initial\nconditions (%)")
    ax.set_title(title,fontsize=fs-3,pad=0.6)
    ax.legend(fontsize=fs-10,frameon=False)#bbox_to_anchor=(1.1, 0.5), loc='center left',frameon=False,fontsize=30)

##############3
##############3
##############3
class MidpointNormalize(matplotlib.colors.Normalize):
    ### Copied from Matplotlib documentation
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

##############3
##############3
##############3
def plotNumberBreakdown(ax,df,maxRegType,pType,fs=None,titleSub=None):
    if not fs:
	fs=20
    emwo,emo,emw = 'EM/WO','EM/O','EM/W'
    ewo,eo,ew = 'E/WO','E/O','E/W'
    mwo,mo,mw = 'M/WO','M/O','M/W'
    if maxRegType=='min':
        maxReg = np.min(df['x'].values)
    else:
        maxReg = np.max(df['x'].values)
        
    xloc = np.argwhere(df['x'].values==maxReg)[:,0][0]
    
    if pType=='MR':
        ax.bar(0.7,df[emwo].values[xloc],color='tab:green',label='E/M',width=0.2)
        ax.bar(1.7,df[emo].values[xloc],color='tab:green',width=0.2)
        ax.bar(2.7,df[emw].values[xloc],color='tab:green',width=0.2)
        ax.bar(0.9,df[mwo].values[xloc],color='tab:orange',label='M',width=0.2)
        ax.bar(1.9,df[mo].values[xloc],color='tab:orange',width=0.2)
        ax.bar(2.9,df[mw].values[xloc],color='tab:orange',width=0.2)
        ax.bar(1.1,df[ewo].values[xloc],color='tab:blue',label='E',width=0.2)
        ax.bar(2.1,df[eo].values[xloc],color='tab:blue',width=0.2)
        ax.bar(3.1,df[ew].values[xloc],color='tab:blue',width=0.2)
        ax.set_xticks([1,2,3])
        ax.set_xticklabels(['W/O','O','W'])
    else:
        ax.bar(0.7,df[emwo].values[xloc],color='tab:green',label='W/O',width=0.2)
        ax.bar(0.9,df[emo].values[xloc],color='tab:blue',label='W',width=0.2)
        ax.bar(1.1,df[emw].values[xloc],color='tab:orange',label='O',width=0.2)
        ax.bar(1.7,df[mwo].values[xloc],color='tab:green',width=0.2)
        ax.bar(1.9,df[mo].values[xloc],color='tab:blue',width=0.2)
        ax.bar(2.1,df[mw].values[xloc],color='tab:orange',width=0.2)
        ax.bar(2.7,df[ewo].values[xloc],color='tab:green',width=0.2)
        ax.bar(2.9,df[eo].values[xloc],color='tab:blue',width=0.2)
        ax.bar(3.1,df[ew].values[xloc],color='tab:orange',width=0.2)
        ax.set_xticks([1,2,3])
        ax.set_xticklabels(['E/M','M','E'])
        
    ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left',frameon=False,fontsize=fs)
    ax.set_ylabel("Steady states (%)")
    if not titleSub:
        ax.set_title('$\lambda=$'+str(maxReg))
    else:
        ax.set_title('$\lambda_{'+titleSub+'}=$'+str(maxReg))

