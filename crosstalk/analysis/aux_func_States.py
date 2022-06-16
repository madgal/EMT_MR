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
cmap = cm.get_cmap('viridis')
for i in range(200):
	__clist+=[mpl.colors.to_hex(cmap(i))]

__clist=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']



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

def getMSE(df_res,compI=False):
    ##columns = u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh

    if not compI:
	        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    else:
	        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_input/EMT_MR_comp_input_"+str(compI)+"_1000_res.txt").dropna()

    full_setp = np.append(uncoupled.values,df_res.values,axis=0)
    mean_list = np.mean(full_setp,axis=0)
    std_list = np.std(full_setp,axis=0)
    uncoupledZ = (uncoupled.values-mean_list)/std_list
    uncoupled_fp = np.unique(np.round(uncoupledZ,6),axis=0)

    stateLabels = returnStateLabels(uncoupled_fp)

    mse_res={'E/O':0,'E/WO':0,'E/W':0,
             'M/O':0,'M/WO':0,'M/W':0,
             'EM/O':0,'EM/WO':0,'EM/W':0}

    for i in range(len(df_res)):
	resZ = (df_res.values[i]-mean_list)/std_list
	
	distances = np.sum((resZ-uncoupled_fp)**2,axis=1)
        location= np.argwhere(np.min(distances)==distances)[:,0][0]
	
	tmp_mse = np.sum((resZ-uncoupled_fp[location])**2,axis=0)
	mse_res[stateLabels[location]]+=tmp_mse

    for key in mse_res:
	mse_res[key] = mse_res[key]*1./(len(df_res))

    return mse_res


def getMSE2(df_res,compI=False):
    ##columns = u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh

    if not compI:
	        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
    else:
	        uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_input/EMT_MR_comp_input_"+str(compI)+"_1000_res.txt").dropna()
    full_setp = np.append(uncoupled.values,df_res.values,axis=0)
    mean_list = np.mean(full_setp,axis=0)
    std_list = np.std(full_setp,axis=0)
    uncoupledZ = (uncoupled.values-mean_list)/std_list
    uncoupled_fp = np.unique(np.round(uncoupledZ,6),axis=0)

    stateLabels = returnStateLabels(uncoupled_fp)

    mse_res={'E/O':0,'E/WO':0,'E/W':0,
             'M/O':0,'M/WO':0,'M/W':0,
             'EM/O':0,'EM/WO':0,'EM/W':0}

    for i in range(len(df_res)):
	resZ = (df_res.values[i]-mean_list)/std_list
	
	distances = np.sum((resZ-uncoupled_fp)**2,axis=1)
        location= np.argwhere(np.min(distances)==distances)[:,0][0]
	
	tmp_mse = 0.5*np.sum((resZ+uncoupled_fp[location])/uncoupled_fp[location],axis=0)/len(resZ)
	mse_res[stateLabels[location]]+=tmp_mse

    for key in mse_res:
	mse_res[key] = mse_res[key]*1./(len(df_res))

    return mse_res



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
	
	#distances=0
	#for j in [0,1,2,4,5,6,7,8,9,10]:
	#	distances+=((resZ[j]-uncoupled_fp[:,j])**2)

	##distances = np.sum((df_res.values[i]-uncoupled_fp)**2,axis=1)
	distances = np.sum((resZ-uncoupled_fp)**2,axis=1)
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
    color_list=[]
    if type(colList[0])==list or type(colList[0])==np.ndarray:
        for i in range(len(colList)):
            if len(colList[i])!=9:
		if colList[i]==['N/A']:
			color_list+=[[colList[i],'w']]
		else:
                	color_list+=[[colList[i],__clist[i]]]
            else:
                color_list+=[[colList[i],'k']]
    else:
        if len(colList)!=9:
		if colList[0]=='N/A':
                	color_list+=[[colList[0],'w']]
		else:
                	color_list+=[[colList,__clist[0]]]
        else:
            color_list+=[[colList,'k']]

    return color_list

def getLegend(colList):
    legend_elements=[]
    if type(colList[0])==list or type(colList[0])==np.ndarray:
        for i in range(len(colList)):
            if len(colList[i])!=9:
		if colList[i]==['N/A']:
                	legend_elements += [Patch(facecolor='w', edgecolor='w', label=colList[i])]
		else:
                	legend_elements += [Patch(facecolor=__clist[i], edgecolor=__clist[i], label=colList[i])]
            else:
                legend_elements += [Patch(facecolor='k', edgecolor='k',
                                         label=colList[i])]
    else:
        if len(colList)!=9:
		if colList[0]=='N/A':
                	legend_elements += [Patch(facecolor='w', edgecolor='w', label=colList[0])]
		else:
	            	legend_elements += [Patch(facecolor=__clist[0], edgecolor=__clist[0],
                                     label=colList)]
        else:
            legend_elements += [Patch(facecolor='k', edgecolor='k',
                                     label=colList)]                


    return legend_elements

def getPlotData(label):
    try:
    	colList=np.unique(label,axis=0)
    except:
    	colList=np.unique(label)
    color=[]

    star_colors=[]## colors where EM/WO exists
    for i in range(len(label)):
        if type(colList[0])==list or type(colList[0])==np.ndarray:
            for j in range(len(colList)):
                if equals(label[i],colList[j]):
                    if len(colList[j])!=9:
			if colList[j]==['N/A']:
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
			if colList[j]=='N/A':
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

def getStateList(res,onlyShowHybrid): 
	tmp =[]
        for key in res:
            if "/" in key and res[key]>0:
    		    tmp+=[key]
        
    	if onlyShowHybrid:
    		reduced=[]
    		for i in range(len(tmp)):
    			reduced += [tmp[i].split("/")[0],tmp[i].split("/")[1]]
    		tmp=[]
        
    	if onlyShowHybrid=='bothI' and (('EM' in reduced) or ('WO' in reduced)):
    		if 'EM' in reduced:
               		tmp=['EM']
    		if 'WO' in reduced:
               		tmp+=['WO']
    	elif onlyShowHybrid=='bothX' and (('EM' in reduced) and ('WO' in reduced)):
               tmp=['EM/WO']
    	elif onlyShowHybrid=='EM' and (('EM' in reduced)):
               tmp=['EM']
    	elif onlyShowHybrid=='WO' and (('WO' in reduced)):
               tmp=['WO']
        
    	if tmp==[]:
    		    tmp=['N/A']

	return tmp


def getStarsForPlot(x1,x2,y1,y2,star_colors,color):

    xmin=np.ones(len(star_colors))*1000.
    xmax=np.ones(len(star_colors))*-1.
    ymin=np.ones(len(star_colors))*1000.
    ymax=np.ones(len(star_colors))*-1.
    for j in range(len(color)):
        if len(star_colors)>0:
            for k in range(len(star_colors)):
                if star_colors[k]==color[j]:
			if xmin[k]>x1[j]:
				xmin[k]=x1[j]
			if xmax[k]<x2[j]:
				xmax[k]=x2[j]
			if ymin[k]>y1[j]:
				ymin[k]=y1[j]
			if ymax[k]<y2[j]:
				ymax[k]=y2[j]
                

    star_x,star_y=[],[]
    for i in range(len(xmax)):
        star_x+=[(xmax[i]-xmin[i])/2.+xmin[i]]
        star_y+=[(ymax[i]-ymin[i])/2.+ymin[i]]
    

    return star_x,star_y

def getDataForPlot(xval,yval,lab):

    x1,x2,y1,y2=[],[],[],[]
    label=[]
    tmp=np.unique(np.array(xval))
    for i in range(len(tmp)):
        el=tmp[i]
        ind = np.argwhere(el==np.array(xval))[:,0]
        tmpz=np.array(lab)[ind]

        ind2= np.argsort(np.array(yval)[ind])
        tmpz=tmpz[ind2]
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

    return x1,x2,y1,y2,label

def getStateListfromFile(df,onlyShowHybrid):
	lab=[]
	for i in range(len(df.values[:])):
		tmp =[]
		for key in ['M/W','M/WO','M/O','EM/W','EM/WO','EM/O','E/W','E/WO','E/O']:
                    if df[key][i]>0:
    	    		    tmp+=[key]
                
    	    	if onlyShowHybrid:
    	    		reduced=[]
    	    		for i in range(len(tmp)):
    	    			reduced += [tmp[i].split("/")[0],tmp[i].split("/")[1]]
    	    		tmp=[]
                
    	    	if onlyShowHybrid=='bothI' and (('EM' in reduced) or ('WO' in reduced)):
    	    		if 'EM' in reduced:
                       		tmp=['EM']
    	    		if 'WO' in reduced:
                       		tmp+=['WO']
    	    	elif onlyShowHybrid=='bothX' and (('EM' in reduced) and ('WO' in reduced)):
                       tmp=['EM/WO']
    	    	elif onlyShowHybrid=='EM' and (('EM' in reduced)):
                       tmp=['EM']
    	    	elif onlyShowHybrid=='WO' and (('WO' in reduced)):
                       tmp=['WO']
                
    	    	if tmp==[]:
    	    		    tmp=['N/A']

		lab+=[tmp]
	return lab


def add_lamdas(cxs,tmp,df):
    lamda_counted=False
    labels=['AS','AZ','AU','HS','HU','INPUT','U3M','U3N','UH','LAMDAUH','LAMBDAUH']
    for i in range(len(tmp)):
	tmp_name = tmp[i].upper()
	if tmp_name=='INPUT':
		cxs[tmp_name] = float(tmp[i+1])
		lamda_counted=True
	elif 'UH' in tmp_name:
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
			scaleL = 0.1
			scaleym1= 0.002
			scaleym2= 0.004
			addym1=0.001
			addym2=0.001
		elif "MR7" in tmp:
			scaleL = 0.2
			scaleym1= 0.01
			scaleym2= 0.2
			addym1=0.
			addym2=0.

		li=[1.,float(l1)*scaleL,float(l2)*scaleL]
		ymi=[0.,float(ym1)*scaleym1+addym1,float(ym2)*scaleym2+addym2]
		uavg = np.mean(df['u'].values)

		cxs['UH'] =P(li,ymi,uavg)
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
