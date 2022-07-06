import numpy as np
import pandas as pd
import os
import math
import time
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


##__clist=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']



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
def returnStateLabels(compI,PSF,noHH=None):
#['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O']

	if (not compI) and (not PSF) and (not noHH):
		uncoupled =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
		uncoupled_fp = np.unique(np.round(uncoupled.values,6),axis=0)
		res = np.unique(np.round(uncoupled.values,6),axis=0)
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
	elif (not PSF) and (not noHH):
		uncoupled =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_input/EMT_MR_comp_input_"+str(compI)+"_1000_res.txt").dropna()
		uncoupled_fp = np.unique(np.round(uncoupled.values,6),axis=0)
		res = np.unique(np.round(uncoupled.values,6),axis=0)
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
	elif (not PSF):
		if noHH=='noEM':
			uncoupled =pd.read_csv("~/Research/EMT_MR/crosstalk/normalCells_coupled/noPartialEMT/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
			uncoupled_fp = np.unique(np.round(uncoupled.values,6),axis=0)
			res = np.unique(np.round(uncoupled_fp,6),axis=0)
			tmpAdd=[]
			for el in res[:,9]:
				inds = np.argwhere(res[:,9]==el)[:,0]
				tmp = []
				for i in range(11):
					tmp+=[np.median(res[inds,i])]
				tmpAdd +=[tmp]

			tmpAdd = np.array(tmpAdd)
			uncoupled_fp = np.concatenate((uncoupled_fp,tmpAdd),axis=0)
			
		elif noHH=='noWO':
			uncoupled =pd.read_csv("~/Research/EMT_MR/crosstalk/normalCells_coupled/normal_metabolic/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
			uncoupled_fp = np.unique(np.round(uncoupled.values,6),axis=0)
			res = np.unique(np.round(uncoupled_fp,6),axis=0)
			tmpAdd=[]
			for el in res[:,0]:
				inds = np.argwhere(res[:,0]==el)[:,0]
				tmp = []
				for i in range(11):
					tmp+=[np.median(res[inds,i])]
				tmpAdd +=[tmp]

			tmpAdd = np.array(tmpAdd)
			uncoupled_fp = np.concatenate((uncoupled_fp,tmpAdd),axis=0)
			
		else:## noHH
			uncoupled =pd.read_csv("~/Research/EMT_MR/crosstalk/normalCells_coupled/bothModified/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
			uncoupled_fp = np.unique(np.round(uncoupled.values,6),axis=0)
			res = np.unique(np.round(uncoupled_fp,6),axis=0)
		
			tmpAdd=[]
			for el in res[:,9]:
				inds = np.argwhere(res[:,9]==el)[:,0]
				tmp = []
				for i in range(11):
					tmp+=[np.median(res[inds,i])]
				tmpAdd +=[tmp]

			tmpAdd = np.array(tmpAdd)
			uncoupled_fp = np.concatenate((uncoupled_fp,tmpAdd),axis=0)
	
			tmpAdd=[]
			for el in res[:,0]:
				inds = np.argwhere(res[:,0]==el)[:,0]
				tmp = []
				for i in range(11):
					tmp+=[np.median(res[inds,i])]
				tmpAdd +=[tmp]

			tmpAdd = np.array(tmpAdd)
			uncoupled_fp = np.concatenate((uncoupled_fp,tmpAdd),axis=0)
			
		final=[]
		res = uncoupled_fp.copy()
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

	else:
		em_values = pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/PSF/getEMSets/state_EM_1000_res.txt").dropna()
		m_values = pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/PSF/getEMSets/state_M_1000_res.txt")
		e_values = pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/PSF/getEMSets/state_E_1000_res.txt")

		em_res = np.unique(em_values.values,axis=0)
		m_res = np.unique(m_values.values,axis=0)
		e_res = np.unique(e_values.values,axis=0)
		uncoupled_fp = np.concatenate((e_res,em_res,m_res),axis=0)
		final=[]
		hvals=np.unique(e_res[:,9])
		for i in range(len(e_res)):
			tmp= 'E/'
			if e_res[i][9]==np.min(hvals):
				tmp+='O'
			elif e_res[i][9]==np.max(hvals):
				tmp+='W'
			else:
				tmp+='WO'
			final+=[tmp]
			
		for i in range(len(em_res)):
			tmp= 'EM/'
			if em_res[i][9]==np.min(hvals):
				tmp+='O'
			elif em_res[i][9]==np.max(hvals):
				tmp+='W'
			else:
				tmp+='WO'
			final+=[tmp]
			
		for i in range(len(m_res)):
			tmp= 'M/'
			if m_res[i][9]==np.min(hvals):
				tmp+='O'
			elif m_res[i][9]==np.max(hvals):
				tmp+='W'
			else:
				tmp+='WO'
			final+=[tmp]
			
	return uncoupled_fp,final

def getStates(df_res,compI=None,PSF=None,noHH=None):
    if compI==None:
         compI=False
    ##columns = u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh

    results={'E/O':0,'E/WO':0,'E/W':0,
             'M/O':0,'M/WO':0,'M/W':0,
             'EM/O':0,'EM/WO':0,'EM/W':0,
             'E':0,'EM':0,'M':0,
             'O':0,'WO':0,'W':0}

    uncoupled_fp,stateLabels = returnStateLabels(compI,PSF,noHH)
    for i in range(len(df_res)):
        prob=df_res.values[i]/uncoupled_fp
        ####tmpS=np.abs(np.sum(prob*np.log10(prob),axis=1))
        tmpS=np.abs(np.sum(np.log10(prob)**2,axis=1))
        loc=np.argwhere(np.min(tmpS)==tmpS)[:,0][0]
        results[stateLabels[loc]]+=1


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
def mapUH_vals(mystr,setV):
    tmp = mystr.split("_")
    #if setV==1:
    #    return float(tmp[0])*0.2+float(tmp[1])*0.02
    #else:
    #    if int(tmp[0])<3:
    #        return float(tmp[0])*0.004+0.003+float(tmp[1])*0.02
    #    else:
    #        return float(tmp[0])*0.04+float(tmp[1])*0.1
    return float(tmp[0])+0.1*float(tmp[1])       
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

def M(i,n,u):
	
	return ((u/10000.)**i)/((1+u/10000.)**n)

def P(li,ymi,u):
	## P=1 no silencing, P=0 full silencing
        km = 0.143 ## for mRNA of Hif-1 as that is the one we are varying
        return (L2(li,u))/(1.+Ym2(ymi,u)/km)
def Ym2(ymi,u):
        total=0
        for i in range(len(ymi)):
                total+= ymi[i]*1.*factorial(len(ymi))/factorial(i)/factorial(len(ymi)-i)*M(i,len(ymi),u)
        return total
def L2(li,u):
        total=0
        for i in range(len(li)):
                total+= li[i]*1.*factorial(len(li))/factorial(i)/factorial(len(li)-i)*M(i,len(li),u)
        return total

def consistent_states(xvalues,yvalues,states):

	consistent=True
	list_err=[]
	for el1 in np.unique(np.array(xvalues)):
		ind_x = np.argwhere(xvalues==el1)[:,0]
		values = yvalues[ind_x]
		for el in np.unique(np.array(values)):
			inds = np.argwhere(values==el)[:,0]
			if len(inds)>1:
				for i in range(1,len(inds)):
					if states[inds[0]]!=states[inds[i]]:
						consistent=False
						list_err +=[el]
						break
	return consistent,list_err

###########
def getColorMatch(colList):
    color_list={}
    if type(colList[0])==list or type(colList[0])==np.ndarray:
        for i in range(len(colList)):
            if len(colList[i])!=9:
                if len(colList[i])==1 and  colList[i]==['N/A']:
                    color_list['w']=colList[i]
                else:
                    color_list[__clist[i]]=colList[i]
            else:
                color_list['k']=colList[i]
    else:
        if len(colList)!=9:
                if colList[0]=='N/A':
                       color_list['w']=colList[0]
                else:
                       color_list[__clist[0]]=colList
        else:
            color_list['k']=colList

    return color_list

def getLegend(regT):
    legend_elements=[]
    if type(regT[0][0])==list or type(regT[0][0])==np.ndarray:
        for i in range(len(regT)):
            colList=regT[i][0]
            ecr=regT[i][1]
            hatchVal=regT[i][2]
            fcr=regT[i][3]
            if len(colList)!=9:
                if len(colList)==1 and  colList==['N/A']:
                	legend_elements += [Patch(facecolor='w',   edgecolor='w',label=colList,hatch=hatchVal)]
                else:
                	legend_elements += [Patch(facecolor=fcr,  edgecolor=ecr, label=colList,hatch=hatchVal)]
            else:
                legend_elements += [Patch(facecolor='k', edgecolor='k',
                                         label=colList,hatch=hatchVal)]
    else:
        colList=regT[0][0]
        ecr=regT[0][1]
        hatchVal=regT[0][2]
        fcr=regT[i][3]
        if len(colList)!=9:
                if colList=='N/A':
                        legend_elements += [Patch(facecolor='w', edgecolor='w',  label=colList,hatch=hatchVal)]
                else:
                        legend_elements += [Patch(facecolor=fcr, edgecolor=ecr,
                                     label=colList,hatch=hatchVal)]
        else:
              legend_elements += [Patch(facecolor='k', edgecolor='k',
                                     label=colList,hatch=hatchVal)]

    legend_elements += [Patch(facecolor='w',label='Black hatch= No Hybrid/Hybrid',hatch='/|')]
    legend_elements += [Patch(facecolor='k',edgecolor='w',label='White hatch=Hybrid/Hybrid',hatch='-')]


    return legend_elements

def getPlotData(label):
    try:
    	colList=np.unique(label,axis=0)
    except:
    	colList=np.unique(label)
    color=[]

    if (type(colList[0])==list or type(colList[0])==np.ndarray) and len(colList)==1:
         colList=colList[0]
    star_colors=[]## colors where EM/WO exists
    for i in range(len(label)):
        if type(colList[0])==list or type(colList[0])==np.ndarray:
            for j in range(len(colList)):
                if equals(label[i],colList[j]):
                    if len(colList[j])!=9:
                         if len(colList[j])==1 and colList[j]==['N/A']:
                                color+=['w']
                                if 'EM/WO' in label[i]:
                                     star_colors+=['w']
                         else:
                                color+=[__clist[j]]
                                if 'EM/WO' in label[i]:
                                      star_colors+=[__clist[j]]
                    else:
                        color+=['k']
                        if 'EM/WO' in label[i]:
                             star_colors+=['k']
                    break
        else:	
            if equals(label[i],colList):
                if len(colList)!=9:
                         if colList[0]=='N/A':
                             color+=['w']
                             if 'EM/WO' in label[i]:
                                star_colors+=['w']
                         else:
                              color+=[__clist[0]]
                              if 'EM/WO' in label[i]:
                                   star_colors+=[__clist[0]]
                else:
                    color+=['k']
                    if 'EM/WO' in label[i]:
                         star_colors+=['k']

    star_colors = list(np.unique(star_colors))
 
    return color,star_colors,colList

def getStateList(res):
        tmp =[]
        for key in res:
            if "/" in key and res[key]>0:
                     tmp+=[key]
        
        if tmp==[]:
                 tmp=['N/A']

        return tmp
def getHatchForPlot(star_colors,color,hybridV,compHybrid,colList):
    matchCols = getColorMatch(colList)## colList[i],plotcol
    regType=[]
    for el in np.unique(color):
       if 'EM/WO' not in matchCols[el]:
            regType+=[[matchCols[el],'k','/|',el]]
       else:
            regType+=[[matchCols[el],'w','-',el]]#el,'',el]]

    hatches=[]
    ecolors=[]

    for j in range(len(color)):
        if len(star_colors)>0:
            added=False
            for k in range(len(star_colors)):
                if star_colors[k]==color[j]:
                        listV = matchCols[color[j]]
                        if 'EM/WO' not in listV:
                                added=True
                                ecolors+=['k']
                                hatches+=['/|']
                        else:
                             if hybridV[j]<compHybrid:
                                    ecolors+=['w']
                                    hatches+=['-']
                             else:# hybridV[j]>= compHybrid:
                                   ecolors+=['w']#star_colors[k]]
                                   hatches+=['-']#'']
                             added=True
            if not added:
                 hatches+=['/|']
                 ecolors+=['k']
        else:## no star colors
                 hatches+=['/|']
                 ecolors+=['k']

    return hatches,ecolors,regType

def getDataForPlot(xval,yval,lab,extras=None):
  

    if type(extras)!=list and type(extras)!=np.ndarray and extras==None:
        extras=[]

    x1,x2,y1,y2=[],[],[],[]
    label=[]
    extraNew=[]
    tmp=np.unique(np.array(xval))
    for i in range(len(tmp)):
        el=tmp[i]
        ind = np.argwhere(el==np.array(xval))[:,0]
        tmpz=np.array(lab)[ind]
        if len(extras)>0:
             tmpextras = np.array(extras)[ind]

        ind2= np.argsort(np.array(yval)[ind])
        tmpz=tmpz[ind2]
        if len(extras)>0:
             tmpextras = tmpextras[ind2]
        tmp2=np.array(yval)[ind][ind2]
        for j in range(len(tmp2)):
                el2=tmp2[j]
                if len(tmp)==1:
                    tmpx1=el-0.5
                    tmpx2=el+0.5
                elif i==0:
                    tmpx1=el
                    tmpx2=(tmp[i+1]-tmp[i])/2.+tmp[i]
                elif el==np.max(tmp):
                    tmpx1=(tmp[i]-tmp[i-1])/2.+tmp[i-1]
                    tmpx2=el
                else:
                    tmpx1=(tmp[i]-tmp[i-1])/2.+tmp[i-1]
                    tmpx2=(tmp[i+1]-tmp[i])/2.+tmp[i]
                if len(tmp2)==1:
                    tmpy1=el2-0.5
                    tmpy2=el2+0.5
                elif j==0:
                    tmpy1=el2
                    tmpy2=(tmp2[j+1]-tmp2[j])/2.+tmp2[j]
                elif el2==np.max(tmp2):
                    tmpy1=(tmp2[j]-tmp2[j-1])/2.+tmp2[j-1]
                    tmpy2=el2
                else:
                    tmpy1=(tmp2[j]-tmp2[j-1])/2.+tmp2[j-1]
                    tmpy2=(tmp2[j+1]-tmp2[j])/2.+tmp2[j]
                label+=[tmpz[j]]
                x1+=[tmpx1]
                x2+=[tmpx2]
                y1+=[tmpy1]
                y2+=[tmpy2]
                if len(extras)>0:
                    extraNew+=[tmpextras[j]]

    if len(extras)>0:
    	return x1,x2,y1,y2,label,extraNew
    return x1,x2,y1,y2,label

def getStateListfromFile(df):

	lab=[]
	for i in range(len(df.values[:])):
                tmp =[]
                for key in ['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O']:
                    if df[key][i]>0:
                         tmp+=[key]
                if tmp==[]:
                      tmp=['N/A']

                lab+=[tmp]
	return lab


def add_lamdas(cxs,tmp,df):
    lamda_counted=False
    labels=['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','LAMDAUH','LAMBDAUH']
    if len(df)>0:
        for i in range(len(tmp)):
                tmp_name = tmp[i].upper()
                if tmp_name=='INPUT':
                    if float(tmp[i+1])<100:
                            cxs[tmp_name] = float(tmp[i+1])*10000.
                    else:
                          cxs[tmp_name] = float(tmp[i+1])
                    lamda_counted=True
                elif 'UH' in tmp_name:
                   li,ymi,yui,strLi,strYm,yustr,val=getUHSet(tmp,i)
                   uavg = np.mean(df['u'].values)
                   cxs['UH'] =P(li,ymi,uavg)+yui

                   cxs['UHV'] =val
                   lamda_counted=True


                elif tmp_name in labels or tmp_name=='IHU':
                   if tmp_name=='IHU':
                        tmp_name='HU'
                   cxs[tmp_name] = float(tmp[i+1])+float(tmp[i+2])/10.**(len(tmp[i+2]))
                   lamda_counted=True 
    return   cxs,lamda_counted

def getPhaseBoundaries(xval,yval,colors):
	# boundList==np.unique(colors)
	#labels=colors

    boundList = np.unique(colors)
    xmin=np.ones(len(boundList))*1000.
    xmax=np.ones(len(boundList))*-1.
    ymin=np.ones(len(boundList))*1000.
    ymax=np.ones(len(boundList))*-1.
    col_list=[]
    for i in range(len(ymax)):
         col_list+=['']

    x1,x2,y1,y2=[],[],[],[]
    tmp=np.unique(np.array(xval))
    for i in range(len(tmp)):
        el=tmp[i]
        ind = np.argwhere(el==np.array(xval))[:,0]
        tmp_list = np.array(colors)[ind]

        ind2= np.argsort(np.array(yval)[ind])
        tmp_list = tmp_list[ind2]
        tmp2=np.array(yval)[ind][ind2]
        for j in range(len(tmp2)):
                el2=tmp2[j]
                if len(tmp)==1:
                    tmpx1=el-0.5
                    tmpx2=el+0.5
                elif i==0:
                    tmpx1=el
                    tmpx2=(tmp[i+1]-tmp[i])/2.+tmp[i]
                elif el==np.max(tmp):
                    tmpx1=(tmp[i]-tmp[i-1])/2.+tmp[i-1]
                    tmpx2=el
                else:
                    tmpx1=(tmp[i]-tmp[i-1])/2.+tmp[i-1]
                    tmpx2=(tmp[i+1]-tmp[i])/2.+tmp[i]
                if len(tmp2)==1:
                       tmpy1=el2-0.5
                       tmpy2=el2+0.5
                elif j==0:
                    tmpy1=el2
                    tmpy2=(tmp2[j+1]-tmp2[j])/2.+tmp2[j]
                elif el2==np.max(tmp2):
                    tmpy1=(tmp2[j]-tmp2[j-1])/2.+tmp2[j-1]
                    tmpy2=el2
                else:
                    tmpy1=(tmp2[j]-tmp2[j-1])/2.+tmp2[j-1]
                    tmpy2=(tmp2[j+1]-tmp2[j])/2.+tmp2[j]

                if len(boundList)>0 and (tmp_list[j] in boundList):
                   for k in range(len(boundList)):
                         if boundList[k]==tmp_list[j]:
                            if xmin[k]>tmpx1:
                                xmin[k]=tmpx1
                            if xmax[k]<tmpx2:
                                xmax[k]=tmpx2
                            if ymin[k]>tmpy1:
                                ymin[k]=tmpy1
                            if ymax[k]<tmpy2:
                                ymax[k]=tmpy2
                            col_list[k] = tmp_list[j]

    bound_x,bound_y=[],[]
    for i in range(len(xmax)):
        bound_x+=[[xmin[i],xmax[i]]]
        bound_y+=[[ymin[i],ymax[i]]]
    

    return bound_x,bound_y,col_list

def __makeLamdaString(my_str,val):

    labels={'AS':'AMPK-|Snail','AZ':'AMPK-|Zeb','AU':'AMPK -> \mu_{200}','HS':'Hif1->Snail','HU':'Hif1-|\mu_{200}','INPUT':'Input->Snail','U3M':'\mu_{34}->mtROS','U3N':'\mu_{34}->noxROS','UH':'\mu_{200}->Hif1'}
    return '$\lambda_{'+str(labels[my_str])+'}$='+str(val)

def getUHSet(tmp,i):
	if len(tmp[i+1])==3:
		l1=tmp[i+1][0]
		ym1=tmp[i+1][1]
		yu1=tmp[i+1][2]
		val=tmp[i+1]
	elif len(tmp[i+1])==4:
		l1=tmp[i+1][0]
		l2=tmp[i+1][1]
		ym1=tmp[i+1][2]
		ym2=tmp[i+1][3]
		val=tmp[i+1]
	else:
		l1=tmp[i+1]
		l2=tmp[i+2]
		ym1=tmp[i+3]
		ym2=tmp[i+4]
		val=tmp[i+1]+tmp[i+2]+tmp[i+3]+tmp[i+4]
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
	elif "MR10" in tmp or "MR"==tmp[1]:### this should* be the one we are using
		l1 = int(l1)
		ym1=int(ym1)
		yu1=int(yu1)
		list_l=[[1.,0.,0.],[1.,0.2,0.2],[1.,0.5,0.3],[1.,0.9,0.8]]
		list_ym=[[0.,0.,0.],[0.,0.002,0.01],[0.,0.01,0.5],[0.,2.,4.]]
		list_yu=[[0.,0.001,0.009],[0.,0.01,0.09],[0.,0.1,0.2],[0.,0.5,1.]]
		li = list_l[l1]
		ymi = list_l[ym1]
		strLi= str(int(l1))
		strYm=str(int(ym1))
		yustr=str(yu1)
		li=list_l[l1]
		ymi=list_ym[ym1]
		yui=yu1
	if "MR10" not in tmp and "MR"!=tmp[1]:
		li=[1.,float(l1)*scaleL,float(l2)*scaleL]
		ymi=[0.,float(ym1)*scaleym1+addym1,float(ym2)*scaleym2+addym2]
		yui=0
               
		strLi= str(int(l1))+','+str(int(l2))
		strYm=str(int(ym1))+','+str(int(ym2))
		li=[1.,float(l1)*scaleL,float(l2)*scaleL]
		ymi=[0.,float(ym1)*scaleym1+addym1,float(ym2)*scaleym2+addym2]

	return li,ymi,yui,strLi,strYm,yustr,val

def getTitle(filen,uval):

    title=''
    labels=['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','LAMDAUH','LAMBDAUH']
    tmp=filen.split("/")[-1].split("_")
    for i in range(len(tmp)):
        tmp_name = tmp[i].upper()
        if tmp_name=='INPUT':
           title+=__makeLamdaString('INPUT',float(tmp[i+1]))
        elif 'UH' in tmp_name:
           li,ymi,yui,strLi,strYm,yustr,val=getUHSet(tmp,i)
           mystr='P='+str(P(li,ymi,uval))+'\tli=[1,'+strLi+']\t ymi=[0,'+strYm+']\t yui='+yustr
           title+=__makeLamdaString('UH',mystr)

        elif tmp_name in labels or tmp_name=='IHU':
           if tmp_name=='IHU':
             tmp_name='HU'
           title+=__makeLamdaString(tmp_name,str(float(tmp[i+1])+float(tmp[i+2])/10.**(len(tmp[i+2]))))

        title+='\t'
	
    return title

#########################

def getCompHybridVal(compI=None):
    if compI==None:
       compI=False

    if not compI:
         df_comp =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
         tempStates=getStates(df_comp)
         compHybrid=tempStates['EM/WO']/10.
    else:
      df_cmp = pd.read_csv("data/crosstalk_input.txt",header=0)
      row = np.argwhere(df_cmp['INPUT'].values==compI)[:,0][0]
      compHybrid = df_cmp['EM/WO'][row]/10.
 
    return compHybrid

def sortDescending(yticks,labels,regList):
	count=0
	len_list=np.zeros(len(yticks.keys()))
	for i in range(len(len_list)):
		len_list[i] = len(yticks[i])
	inds = np.argsort(len_list)[::-1]

	new_ticks={}
	count=0
	reg_label={}
	for i in inds:
		new_ticks[count]=yticks[i]
		reg_label[count]=regList[i]
		count+=1

	new_tickLabels={}
	for key in new_ticks:
		for el in new_ticks[key]:
			new_tickLabels[el]=''

	for key in new_ticks:
		checkList=[]
		for el in new_ticks[0]:
			for i in range(len(labels[el])):
				if labels[el][i].keys()[0]==reg_label[key]:
					if labels[el][i].values()[0] not in checkList:
						##checkList+=[labels[el][i].values()[0]]
						if new_tickLabels[el]=='':
							new_tickLabels[el]=str(reg_label[key])+'='+str(labels[el][i].values()[0])
						else:
							new_tickLabels[el]+='   '+str(reg_label[key])+'='+str(labels[el][i].values()[0])
					else:
						new_tickLabels[el]+=''
					break

	keylist=np.sort(new_tickLabels.keys())
	ticklabs=[]
	for i in keylist:
		ticklabs+=[new_tickLabels[i]]
		
	return new_ticks,ticklabs



#########################
#########################
#########################
def getFullReg(df,reg_list,label_dir,tickList,opp_label,count=None,i=None,tmp=None,els=None):

    if count==None:
        count=0
    if i==None:
         i=0
    if tmp==None:
        tmp={}
    if els==None:
        els={}
    ### reg_list = {'x':['A','B','...],'y':['A2','B2','...]}
    ### label_dir = {'x':{'A':{f1:[],f2:[]}}}, etc where f1 and f2 are foldchang values
    
    if len(label_dir)==1:
        vals = np.sort(df[reg_list[i].upper()])
        for el in vals:#
            tickList[0]+=[count]
            opp_label[count]=[{reg_list[i]:el}]
            count+=1
    elif i==len(label_dir)-1:
        for el in np.sort(label_dir[i].keys()):
            inds={}
            for j in els: 
                inds[j] = np.argwhere(df[reg_list[j].upper()]==els[j])[:,0]
                
            if len(inds.keys())>1:
                indsf = np.intersect1d(inds[0],inds[1])
            elif len(inds.keys())==1:
                indsf = inds[0]
            else:
                indsf=[]

            if len(inds.keys())>2:
                for j in range(2,len(inds.keys())):
                    indsf = np.intersect1d(indsf,inds[j])
            
            if len(inds.keys())>0 and len(indsf)>0:
                label_dir[i][el]+=[count]
                if count not in opp_label:
                    opp_label[count]=[{reg_list[i]:el}]
                tickList[i]+=[count]
                
                is_el = np.argsort(els.values())
                for j in is_el:#els:
                    label_dir[j][els[j]]+=[count]
                    opp_label[count]+=[{reg_list[j]:els[j]}]
                for j in els:
                    tmp[j]+=[count]
                count+=1
    else:
        if i==0:
            els={}
            tmp={}
            for j in label_dir:
                tmp[j]=[]
        for el in np.sort(label_dir[i].keys()):
            els[i]=el
            label_dir,tickList,opp_label,count,tmp=getFullReg(df,reg_list,label_dir,tickList,opp_label,count,i+1,tmp,els)
            tickList[i]+=[np.max(tmp[i])]
            
    return label_dir,tickList,opp_label,count,tmp
##################3
#########################
#########################
def getLamdaSets(lvx,lvy,xk=None,yk=None,count=0,tmpx=None,tmpy=None,finx=None,finy=None):

    if finx is None:
       finx=[]
    if finy is None:
       finy=[]
    if tmpx is None:
        tmpx=[]
    if tmpy is None:
       tmpy=[]

    if xk is None and yk is None:
        xk = lvx.keys()
        yk = lvy.keys()
    if count==(len(xk)+len(yk)-1):
        ntx,nty = [],[]
        for i in range(len(lvy[yk[-1]])):
            ty = tmpy+[lvy[yk[-1]][i]]
            ntx+=[tmpx]
            nty+=[ty]
        finx+=ntx
        finy+=nty
        return finx,finy
    else:
        if count <len(xk):
            for i in range(len(lvx[xk[count]])):
                t2 = tmpx+[lvx[xk[count]][i]]
                finx,finy=getLamdaSets(lvx,lvy,xk,yk,count+1,t2,tmpy,finx,finy)
        elif count <len(xk)+len(yk):
            for i in range(len(lvy[yk[count%len(xk)]])):
                t2 = tmpy+[lvy[yk[count%len(xk)]][i]]
                finx,finy=getLamdaSets(lvx,lvy,xk,yk,count+1,tmpx,t2,finx,finy)
        else:
            return finx,finy
    return finx,finy
        
##
def getPlotData_regulatory(df,lamdaList,regList):
    allLabels = getStateListfromFile(df)
    nics = df['nics'][0]*1.

    xlabels,ylabels= getLamdaSets(lamdaList['x'],lamdaList['y'])

    xvals=np.zeros(len(xlabels))
    xunique=np.unique(xlabels,axis=0)
    for i in range(len(xunique)):
        xvals = xvals+(np.mean(xunique[i]==xlabels,axis=1)==1)*i
    yvals=np.zeros(len(ylabels))
    yunique=np.unique(ylabels,axis=0)
    for i in range(len(yunique)):
       yvals = yvals+(np.mean(yunique[i]==ylabels,axis=1)==1)*i

    label,hybridV=[],[]
    for i in range(len(xlabels)):
	  ##start with the x
          inds = np.argwhere(xlabels[i][0]==df[regList['x'][0].upper()])[:,0]
          for j in range(1,len(xlabels[i])):
             tmp = np.argwhere(xlabels[i][j]==df[regList['x'][j].upper()])[:,0]
             inds = np.intersect1d(tmp,inds)
          for k in range(len(ylabels[i])):
             tmp = np.argwhere(ylabels[i][k]==df[regList['y'][k].upper()])[:,0]
             inds = np.intersect1d(tmp,inds)

          if len(inds)==1:
              label+=[allLabels[inds[0]]]
              hybridV+=[df['EM/WO'].values[inds[0]]/nics*100.]
          else:
              label+=[[-1]]
              hybridV+=[-1]
              if len(inds)>0:
                  print("Too many inds" )

    indsF= np.argwhere(np.array(hybridV)!=-1)[:,0]
    xvalF = xvals[indsF]
    yvalF = yvals[indsF]
    labelF = list(np.array(label)[indsF])
    hybridVF = list(np.array(hybridV)[indsF])
    xlabelsF = list(np.array(xlabels)[indsF])
    ylabelsF = list(np.array(ylabels)[indsF])

    return xvalF,yvalF,labelF,hybridVF,xlabelsF,ylabelsF

##################
def determineRegRange(filen,comp,compI):
    
    df =pd.read_csv(filen).dropna()
    if len(df)>0:

        file0 = open("upreg_files.txt","a")
        fileF = open("onlyHH_files.txt","a")
        file1 = open("eq_files.txt","a")
        file2 = open("downreg_files.txt","a")
        file3 = open("none_files.txt","a")

        for i in range(len(df)):
                tmp = int(df['INPUT'].values[i])
                if tmp!=50000:
                      if tmp not in compI.keys():
                              print(tmp," not in input comparison dictionary")
			
                      compV = compI[tmp]
                else:
                    compV=comp

                if compV<=0:
                       compV=0

                key = 'EM/WO'
                ## header is AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH
                listV = df['EM/WO'].values[i]
                pdif = np.abs(listV-compV)*1./compV*100.
                if listV==df['nics'].values[i]:
                      fileF.write("%s,%s,%s" %(filen,listV/10.,pdif))
                      for k in ['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','UHV']:
                          try:
                                fileF.write(",%s" %df[k].values[i])
                          except:
                                fileF.write(",")
                      fileF.write("\n")
                elif listV==0:
                      file3.write("%s,%s,%s" %(filen,listV/10.,pdif))
                      for k in ['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','UHV']:
                          try:
                              file3.write(",%s" %df[k].values[i])
                          except:
                                file3.write(",")
                      file3.write("\n")
                elif listV>compV:
                      file0.write("%s,%s,%s" %(filen,listV/10.,pdif))
                      for k in ['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','UHV']:
                          try:
                              file0.write(",%s" %df[k].values[i])
                          except:
                                file0.write(",")
                      file0.write("\n")
                elif listV==compV:
                      file1.write("%s,%s,%s" %(filen,listV/10.,pdif))
                      for k in ['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','UHV']:
                          try:
                              file1.write(",%s" %df[k].values[i])
                          except:
                                file1.write(",")
                      file1.write("\n")
                else:
                      file2.write("%s,%s,%s" %(filen,listV/10.,pdif))
                      for k in ['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','UHV']:
                          try:
                              file2.write(",%s" %df[k].values[i])
                          except:
                                file2.write(",")
                      file2.write("\n")

        file0.close()
        file1.close()
        file2.close()
        fileF.close()
        file3.close()
##################
def determineGroups(filen):
    
	df =pd.read_csv(filen).dropna()
	if len(df)>0:

		compV={0:['E/O','M/W'],1:['E/W','M/O'],2:['E/O','M/W','EM/WO'],3:['E/W','M/O','EM/WO'],4:['E/WO','M/O','EM/W'],5:['E/WO','M/W','EM/O'],6:['E/O','M/WO','EM/W'],7:['E/W','M/WO','EM/O'],8:['EM/W','E/O'],9:['EM/W','E/WO'],10:['EM/O','E/W'] ,11:['EM/O','E/WO'] ,12:['EM/WO','E/W'] ,13:['EM/WO','E/O'] ,14:['EM/W','M/O'] ,15:['EM/W','M/WO'] ,16:['EM/O','M/W'] ,17:['EM/O','M/WO'] ,18:['EM/WO','M/W'] ,19:['EM/WO','M/O']}
		comps={0:0,1:0,2:0, 3:0,4:0, 5:0, 6:0, 7:0 ,8:0 ,9:0 ,10:0 ,11:0 ,12:0 ,13:0 ,14:0 ,15:0 ,16:0 ,17:0 ,18:0 ,19:0}
		file0 = open("groups.txt","a")

		for i in range(len(df)):
                        ## header is AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV,E,EM,M,W,WO,O,M/W,M/WO,M/O,EM/W,EM/WO,EM/O,E/W,E/WO,E/O,nics

			tmp=0
			labels=[]
			for key in ['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O']:
				if df[key].values[i]>0:
					tmp+=1
					labels+=[key]


			if len(labels)<=3:
				#print labels
				for el in comps.keys():
					#print el,compV[el],equals(compV[el],labels)
					if equals(compV[el],labels): 
						#print compV[el]
						comps[el]+=1



			if tmp<=3:
				tmpEMT = []
				tmpMR =[]
				for i in range(len(labels)):
					tmp = labels[i].split("/")
					if tmp[0] not in tmpEMT:
						tmpEMT+=[tmp[0]]
					if tmp[1] not in tmpMR:
						tmpMR+=[tmp[1]]
				
				if len(tmpEMT)==len(tmpMR) and len(tmpMR)==len(labels):
					file0.write("%s" %filen)
					for i in range(3):
						if i>=len(labels):
							file0.write(",")
						else:
							file0.write(",%s" %labels[i])

					file0.write("\n")

		file0.close()


		file1 = open("groups_list.txt","a")
		for el in comps:
					if comps[el]>0:
						file1.write("%s" %filen)
						for i in range(3):
							if i>= len(compV[el]):
								file1.write(",")
							else:
								file1.write(",%s" %compV[el][i])
						file1.write("\n")
		file1.close()


