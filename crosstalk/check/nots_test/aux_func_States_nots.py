import numpy as np
import pandas as pd
import os
import math
from math import factorial
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib as mpl

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



########## Metabolic Phenotypes
######{'A': 436.47213160474945, 'Rnox': 1.6461845679152642, 'Rmt': 128.97826781691597, 'h': 26.69259410719074}
######{'A': 385.53355124455993, 'Rnox': 2.857316799224224, 'Rmt': 137.55212530487026, 'h': 126.85433451255018}
######{'A': 327.50827624867554, 'Rnox': 8.970170738431756, 'Rmt': 138.2393863317781, 'h': 292.59779290727334}
######{'A': 163.0219643626051, 'Rnox': 21.429749060524145, 'Rmt': 54.95997630766442, 'h': 381.53445469342324} 
######{'A': 68.67768314578954, 'Rnox': 30.525042671185414, 'Rmt': 8.353361802666083, 'h': 480.4124463778733} 
########## EMT phenotypes
######'u': 1265.0145554683656, 'ms': 451.1540645013371, 'u3': 8177.938943062448, 'mz': 990.4010240761947},
######'u': 4322.718712975423, 'ms': 457.2059660505549, 'u3': 14831.841141418037, 'mz': 592.8640392330298},
######'u': 12389.217608300634, 'ms': 458.4642670196209, 'u3': 16908.617669803563, 'mz': 301.84367817798204}, 
######'u': 15365.410894097318, 'ms': 458.5016342609288, 'u3': 16975.501915824585, 'mz': 179.92788317173316})
######'u': 19098.53802971536, 'ms': 458.5115606841859, 'u3': 16993.323702636688, 'mz': 62.206071606052106},

#######################
#######################
def stateThresholds():
    return {'E':{'u':15365.41,'mz':179.928},'M':{'u':4322.72,'mz':592.86},
            'EM':{'mz':[179.928,592.86],'u':[4322.72,15365.41]},
            'W':{'h':381.53,'A':163.02},'O':{'h':126.85,'A':385.53},
            'WO':{'h':[126.85,381.53],'A':[163.02,385.53]}
           }



#########################
#########################
#########################
def returnStateLabels(res):
#['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O']

	final=[]
	uvals=np.unique(res[:,0])
	hvals=np.unique(res[:,9])
	for i in range(len(res)):
		if res[i][0]==np.min(uvals):
			tmp='M/'
		elif res[i][0]==np.max(uvals):
			tmp='E/'
		else:
			tmp='EM/'
		if res[i][9]==np.min(hvals):
			tmp+='O'
		elif res[i][9]==np.max(hvals):
			tmp+='W'
		else:
			tmp+='WO'
		
		final+=[tmp]
	return final

def getStates(df_res,compI=False):
    ##columns = u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh

    results={'E/O':0,'E/WO':0,'E/W':0,
             'M/O':0,'M/WO':0,'M/W':0,
             'EM/O':0,'EM/WO':0,'EM/W':0,
             'E':0,'EM':0,'M':0,
             'O':0,'WO':0,'W':0}

    if not compI:
	        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    else:
	        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_input/EMT_MR_comp_input_"+str(compI)+"_1000_res.txt").dropna()
    full_setp = np.append(uncoupled.values,df_res.values,axis=0)
    #mean_list = np.mean(uncoupled.values,axis=0)
    mean_list = np.mean(full_setp,axis=0)
    std_list = np.std(full_setp,axis=0)
    uncoupledZ = (uncoupled.values-mean_list)/std_list
    uncoupled_fp = np.unique(np.round(uncoupledZ,6),axis=0)

    stateLabels = returnStateLabels(uncoupled_fp)


    for i in range(len(df_res)):
	resZ = (df_res.values[i]-mean_list)/std_list
	
	distances = np.sum((resZ-uncoupled_fp)**2,axis=1)
        #print 'A',np.argwhere(np.min(distances)==distances)[:,0],np.sum(distances>0),distances
        location= np.argwhere(np.min(distances)==distances)[:,0][0]
	results[stateLabels[location]]+=1


    results['E'] = results['E/W']+ results['E/WO']+ results['E/O']
    results['EM'] = results['EM/W']+ results['EM/WO']+ results['EM/O']
    results['M'] = results['M/W']+ results['M/WO']+ results['M/O']
    results['W'] = results['E/W']+ results['EM/W']+ results['M/W']
    results['WO'] = results['E/WO']+ results['EM/WO']+ results['M/WO']
    results['O'] = results['E/O']+ results['EM/O']+ results['M/O']

    return results


###############3
###############3
###############3
def equals(list1,list2):
    if len(list1)!=len(list2):
        return False

    for el in list1:# since they are equal lengths
        if el not in list2:
            return False
        
    # all elements are in both lists and are of equal length
    return True

#########################
#########################
#########################
#########################
#########################
def retMeshGrid(x,y,c):
    cx = np.unique(x)
    cy = np.unique(y)
    res={}
    for el1 in cx:
        res[el1]={}    
        for el2 in cy:
            res[el1][el2]=-1

    for i in range(len(x)):
        res[x[i]][y[i]]=c[i]

    xval,yval,zval=[],[],[]
    for k in res:
        tmp1,tmp2,tmp3=[],[],[]
        for k2 in res[k]:
            tmp1+=[mapUH_vals(k,1)]
            tmp2+=[mapUH_vals(k2,2)]
            tmp3+=[res[k][k2]]

        xval+=[tmp1]
        yval+=[tmp2]
        zval+=[tmp3]
        
    return xval,yval,zval
###########
