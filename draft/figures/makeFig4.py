import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib import cm
from shapely import geometry 
from descartes import PolygonPatch
import alphashape

from PIL import Image

from aux_func_States import getLegend
from aux_func_States import equals

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'Times New Roman'

labelD={'W':'W','O':'O','WO':'W/O',
        'E':'E','EM':'E/M','M':'M',
        'E/O':'E-O','E/W':'E-W','E/WO':'E-W/O',
        'EM/O':'E/M-O','EM/W':'E/M-W','EM/WO':'E/M-W/O',
        'M/O':'M-O','M/W':'M-W','M/WO':'M-W/O'}

def getColor(label):
    
    ### Reduced EMT states
    if equals(label,['E']):
        return 'r'
    if equals(label,['M']):
        return 'b'
    if equals(label,['EM']):
        return 'y'
    if equals(label,['E','M']):
        return 'purple'
    if equals(label,['E','EM']):
        return 'orange'
    if equals(label,['M','EM']):
        return 'g'
    if equals(label,['E','M','EM']):
        return 'k'
    
    ### REduced to Metabolic states
    if equals(label,['W']):
        return 'r'
    if equals(label,['O']):
        return 'b'
    if equals(label,['WO']):
        return 'y'
    if equals(label,['W','O']):
        return 'purple'
    if equals(label,['W','WO']):
        return 'orange'
    if equals(label,['O','WO']):
        return 'g'
    if equals(label,['O','W','WO']):
        return 'k'
    
    ### make the all 9 coupled states always be black
    if equals(label,['E/O','E/WO','EM/O','M/O','EM/WO','M/WO','M/W','EM/W','E/W']):
        return 'k'
    
    
    ##generate the color dictionary based on the keys so that you can easily add results (511 total possible combos)
    keys=[  ['E/W'],['M/W'],['EM/W'],
            ['E/O'],['M/O'],['EM/O'],
            ['E/WO'],['M/WO'],['EM/WO'],

            ['E/O','EM/WO'],
            ['E/O','E/WO'],
            ['E/O','E/W'],

            ['E/O','E/W','E/WO'],
            ['E/O','E/WO','EM/W'],
            ['E/O','E/WO','M/W'],
            ['E/O','EM/W','EM/WO'],
            ['E/O','M/W','EM/WO'],
            ['E/O','M/W','M/WO'],
            ['M/W','M/WO','M/O'],
            ['E/O','EM/O','M/O'],

            ['E/O','EM/WO','M/W','M/WO'],
            ['E/W','EM/WO','M/O','M/WO'],
            ['E/WO','EM/WO','M/O','M/WO'],
            ['E/W','EM/W','M/O','M/WO'],
            ['E/W','E/WO','E/O','M/W'],
            ['EM/O','M/W','M/WO','M/O'],
            ['E/O','EM/O','M/O','M/W'],

            ['E/W','EM/W','M/O','M/W','M/WO'],
            ['E/WO','EM/O','EM/WO','M/O','M/WO'],
            ['E/W','EM/WO','M/W','M/O','M/WO'],
            ['E/WO','EM/WO','M/W','M/O','M/WO'],
            ['E/O','E/W','E/WO','EM/W','M/W'],
            ['E/O','E/W','E/WO','M/WO','M/W'],
            ['E/O','EM/O','M/O','M/WO','M/W'],
            ['EM/WO','EM/O','M/O','M/WO','M/W'],

            ['E/O','E/WO','EM/O','EM/WO','M/O','M/WO'],
            ['E/WO','EM/O','EM/WO','M/O','M/W','M/WO'],
            ['E/O','E/W','E/WO','EM/W','M/WO','M/W'],
            ['E/O','E/W','E/WO','M/O','M/WO','M/W'],
            ['E/O','EM/O','EM/WO','M/O','M/W','M/WO'],
            ['EM/W','EM/O','EM/WO','M/O','M/W','M/WO'],
            ['E/O','EM/O','EM/W','M/O','M/W','M/WO'],

            ['E/O','E/WO','EM/O','EM/WO','M/WO','M/O','M/W'],
            ['E/O','E/WO','E/W','EM/W','M/WO','M/O','M/W'],
            ['E/O','EM/WO','EM/O','EM/W','M/WO','M/O','M/W'],

            ['E/O','E/WO','EM/O','EM/W','EM/WO','M/O','M/WO','M/W'],
            ['E/O','E/WO','E/W','EM/W','EM/WO','M/W','M/O','M/WO'],
             ]
    
    __clist=[]
    cmap1= cm.get_cmap('tab20b')
    cmap2= cm.get_cmap('tab20c')
    count,base=0,0
    for i in range(20):
        __clist+=[matplotlib.colors.to_hex(cmap1(base+count*4))]
        __clist+=[matplotlib.colors.to_hex(cmap2(base+count*4))]
        count+=1
        if count==5:
                count=0
                base+=1
    __clist+=['yellow','beige','chartreuse','cyan','lime','magenta','darkviolet','blue','red','saddlebrown']
                
    cmap = matplotlib.cm.get_cmap('hsv')
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=len(keys))

    for i in range(len(keys)):
        if equals(label,keys[i]):
            if i>=len(__clist):
                print("NEED MORE COLORS")
            return __clist[i%len(__clist)]
            #return cmap(norm(i))

def colorMap(fileN,reduced=None):
    df = pd.read_csv(fileN)
    lab={}
    if reduced:
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

fig = plt.figure(figsize=(5,5))
gs1 = gridspec.GridSpec(3,3, height_ratios=[1.,.2,1], width_ratios=[1,0.3,1])
gs1.update(left=0.05, right=2., wspace=0.05,hspace=0.05,top=2.,bottom=0)
fs =20
matplotlib.rcParams.update({'font.size':fs})
matplotlib.rcParams.update({'hatch.linewidth':3})


__clist=[]
cmap1= cm.get_cmap('tab20b')
cmap2= cm.get_cmap('tab20c')

count,base=0,0
for i in range(20):
    __clist+=[matplotlib.colors.to_hex(cmap1(base+count*4))]
    __clist+=[matplotlib.colors.to_hex(cmap2(base+count*4))]
    count+=1
    if count==5:
            count=0
            base+=1

fillColors=__clist#['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']

##########
##########
ax3a = plt.subplot(gs1[0,0])
ax3a.text(-0.1, 1.05, 'A',transform=ax3a.transAxes,verticalalignment='top', horizontalalignment='right',color='black', fontsize=fs+5)
 
labs = colorMap("data_3a_option1_legend.txt")
df = pd.read_csv("data_3a_option1.txt")
for i in range(len(df)):
    ax3a.fill_between([df['X1'].values[i],df['X2'].values[i]],df['Y1'].values[i],df['Y2'].values[i],facecolor=getColor(labs[df['color'].values[i]]),linewidth=0.0)

dfHE = pd.read_csv("data_3a_option1_HHexists.txt")
for i in range(len(dfHE)):
    ax3a.fill_between([dfHE['X1'].values[i],dfHE['X2'].values[i]],dfHE['Y1'].values[i],dfHE['Y2'].values[i],facecolor='none',hatch='.',edgecolor='k',linewidth=0.0)
dfHO = pd.read_csv("data_3a_option1_HHOnlhy.txt")
for i in range(len(dfHO)):
    ax3a.fill_between([dfHO['X1'].values[i],dfHO['X2'].values[i]],dfHO['Y1'].values[i],dfHO['Y2'].values[i],facecolor='none',hatch='.',edgecolor='r',linewidth=0.0)

dfl = pd.read_csv("data_3a_option1_legend.txt")
legend_elements=[]
for i in range(len(dfl)):
    tmp=[]
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
        if dfl[key].values[i]==1:
            tmp+=[key]
    legend_elements+=[ Patch(facecolor=fillColors[dfl['t3'].values[i]], 
                             edgecolor=dfl['t1'].values[i],label=tmp)]
#ax3a.legend(handles=legend_elements, bbox_to_anchor=(2.1, -.2))
ax3a.set_xlabel("$\lambda_{\mu_{34}->noxROS}$")
ax3a.set_ylabel("$\lambda_{\mu_{34}->mtROS}$")
ax3a.set_title("Regulation of metabolism by EMT can only \nresult in upregulated H/H")

##########
ax3b = plt.subplot(gs1[0,2])
ax3b.text(-0.1, 1.05, 'B',transform=ax3b.transAxes,verticalalignment='top', horizontalalignment='right',color='black', fontsize=fs+5)
 
labs = colorMap("data_3b_option2_legend.txt")
df = pd.read_csv("data_3b_option2.txt")
for i in range(len(df)):
        ax3b.fill_between([df['X1'].values[i],df['X2'].values[i]],df['Y1'].values[i],df['Y2'].values[i],facecolor=getColor(labs[df['color'].values[i]]),linewidth=0.0)

dfHE = pd.read_csv("data_3b_option2_HHexists.txt")
for i in range(len(dfHE)):
    ax3b.fill_between([dfHE['X1'].values[i],dfHE['X2'].values[i]],dfHE['Y1'].values[i],dfHE['Y2'].values[i],facecolor='none',hatch='.',edgecolor='k',linewidth=0.0)
dfHO = pd.read_csv("data_3b_option2_HHOnlhy.txt")
for i in range(len(dfHO)):
    ax3b.fill_between([dfHO['X1'].values[i],dfHO['X2'].values[i]],dfHO['Y1'].values[i],dfHO['Y2'].values[i],facecolor='none',hatch='.',edgecolor='r',linewidth=0.0)

    
dfl = pd.read_csv("data_3b_option2_legend.txt")
legend_elements=[]
for i in range(len(dfl)):
    tmp=[]
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
        if dfl[key].values[i]==1:
            tmp+=[key]
    legend_elements+=[ Patch(facecolor=fillColors[dfl['t3'].values[i]], 
                             edgecolor=dfl['t1'].values[i],label=tmp)]
#ax3b.legend(handles=legend_elements, bbox_to_anchor=(1.1, -0.5))
ax3b.set_xlabel("$\lambda_{AMPK-|Snail}$")#"$\lambda_{Hif1-|\mu_{200}}$")
ax3b.set_ylabel("$\lambda_{AMPK-|\mu_{200}}$")#"Input to Snail")
ax3b.set_title("Regulation of EMT by metabolism\n can only result in upregulated H/H")



##########
ax3c = plt.subplot(gs1[2,0])
ax3c.text(-0.1, 1.05, 'C',transform=ax3b.transAxes,verticalalignment='top', horizontalalignment='right',color='black', fontsize=fs+5)
labs = colorMap("data_3c0_legend.txt")
df = pd.read_csv("data_3c0.txt")
for i in range(len(df)):
        ax3c.fill_between([df['X1'].values[i],df['X2'].values[i]],df['Y1'].values[i],df['Y2'].values[i],facecolor=getColor(labs[df['color'].values[i]]),linewidth=0.0)

#dfHE = pd.read_csv("data_3c0_HHexists.txt")
#for i in range(len(dfHE)):
#    ax3c.fill_between([dfHE['X1'].values[i],dfHE['X2'].values[i]],dfHE['Y1'].values[i],dfHE['Y2'].values[i],facecolor='none',hatch='.',edgecolor='k',linewidth=0.0)
#dfHO = pd.read_csv("data_3c0_HHOnlhy.txt")
#for i in range(len(dfHO)):
#    ax3c.fill_between([dfHO['X1'].values[i],dfHO['X2'].values[i]],dfHO['Y1'].values[i],dfHO['Y2'].values[i],facecolor='none',hatch='.',edgecolor='r',linewidth=0.0)

dfl = pd.read_csv("data_3c0_legend.txt")
legend_elements=[]
for i in range(len(dfl)):
    tmp=[]
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
        if dfl[key].values[i]==1:
            tmp+=[key]
    legend_elements+=[ Patch(facecolor=fillColors[dfl['t3'].values[i]], 
                             edgecolor=dfl['t1'].values[i],label=tmp)]
#ax3c.legend(handles=legend_elements, bbox_to_anchor=(2.5, 0.8))
ax3c.set_xlabel("$\lambda_{Hif1-|\mu_{200}}$")
ax3c.set_ylabel("$\lambda_{\mu_{34}->mtROS}$")
legend_elements =[Patch(facecolor=fillColors[0],label='downregulated E/M-W/O'),
                    Patch(facecolor=fillColors[1],label='upregulated E/M-W/O'),
                    Patch(facecolor=fillColors[2],label='no change E/M-W/O'),
                    Patch(facecolor=fillColors[3],label='only E/M-W/O exists')]

dfHE = pd.read_csv("data_3c0_HHOnlhy.txt")
tmpx=list(dfHE['X1'].values)
tmpx+=list(dfHE['X2'].values)
tmpy=list(dfHE['Y1'].values)
tmpy+=list(dfHE['Y2'].values)
    
points = np.zeros((len(tmpx),2))
points[:,0]=tmpx
points[:,1]=tmpy
ashape = alphashape.alphashape(points,8)
print(ashape)
ax3c.add_patch(PolygonPatch(ashape,alpha=0.9,hatch='.',edgecolor='r',facecolor='none'))
ax3c.scatter(tmpx,tmpy)
    
plt.show()
fig.savefig("Figure4.png",bbox_inches='tight')#,dpi=300)
