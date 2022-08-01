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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import ticker
import sys
from descartes import PolygonPatch
import alphashape


from aux_func_States import *

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



###########################
###########################
############################
def CompData(fileOut,dirN):
    fileo = open("data/"+fileOut+".txt","w")
    fileo.write("icsNumber")
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
        fileo.write(",%sAvg,%sStd" %(key,key))
    fileo.write("\n")


    for fileN in os.listdir("/home/madeline/Research/EMT_MR/crosstalk/coupledWReg_Ccode/"+dirN+"/"):
        if 'EMT' in fileN:
            fileS = fileN.split("_")[:-3]
            fileS= '_'.join(fileS)+"_"
            break

    keys=[]
    res={}
    for x in [100,500,1000,2000]:#,5000,10000]: 
        if x not in res.keys():
            res[x]={'E/O':[],'E/W':[],'E/WO':[],'EM/O':[],'EM/W':[],'EM/WO':[],'M/O':[],'M/W':[],'M/WO':[],'E':[],'EM':[],'M':[],'O':[],'WO':[],'W':[]}

        for i in range(10):
                df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/"+dirN+"/"+fileS+str(i)+"_"+str(x)+"_res.txt").dropna()

                mapRes = getStates(df)
                for key in mapRes:
                    res[x][key] +=[mapRes[key]/float(x)*100.]

        fileo.write("%s" %x)
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
             fileo.write(",%s,%s" %(np.mean(res[x][key]),np.std(res[x][key])))
        fileo.write("\n")
    fileo.close()

###############################
###############################
def getDataSets_ics(name,Ftitle,xaxis,start=None,constants=None,title=None):
    if constants!=None:
        for key in constants:
            name+="_"+key
    if start==None:
        start=''
    
    if title==None:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+start+"crosstalk_"+name+".txt",header=0)
    else:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+title,header=0)
    nics = float(df['nics'][0])

    if constants!=None:
        keys =list(constants.keys())
        key = keys[0]
        if type(constants[key])!=list:
            inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        else:
            inds = list(np.argwhere(df[key.upper()].values==constants[key][0])[:,0])
            for i in range(1,len(constants[key])):
                inds +=list(np.argwhere(df[key.upper()].values==constants[key][i])[:,0])
            inds = np.array(inds)
        for i in range(1,len(keys)):
           if type(constants[keys[i]])!=list:
                tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
                inds = np.intersect1d(inds,tmpI)
           else:
                tmpI=[]
                for j in range(len(constants[keys[i]])):
                    tmpI += list(np.argwhere(df[keys[i].upper()].values==constants[keys[i]][j])[:,0])
                inds = np.intersect1d(inds,np.array(tmpI))
    else:
        inds=np.argwhere(df['nics'].values>-1)[:,0]


    fileo = open("data/"+Ftitle+"_ics.txt",'w')
    fileo.write('x,E/O,E/W,E/WO,EM/O,EM/W,EM/WO,M/O,M/W,M/WO,E,EM,M,O,WO,W\n')
    for i in range(len(inds)):
        fileo.write("%s" %(df[xaxis.upper()][inds[i]]))
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
            fileo.write(",%s" %(df[key][inds[i]]/nics*100.))
        fileo.write("\n")


    if 1 not in df[xaxis.upper()].values and start==None:
       	df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
        mapRes = getStates(df)
        fileo.write("1.")
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
            fileo.write(",%s" %(mapRes[key]/10.))
        fileo.write("\n")


    fileo.close()

 
############
############
def getDataSets_uhSinglesnoWO(title):
	df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/singles_crosstalk_uh10.txt")
	elV = {'O':1.,'WO':2.,'W':3.}
	y = df['O'].values*elV['O']+df['WO'].values*elV['WO']+df['W'].values*elV['W']
	xun_ok=[]
	for i in range(4):
		for j in range(1,4):
			for k in range(4):
				if (j==1 and k<=0) or (j==2 and k<=1) or (j==3):
					xun_ok+=[i*100+j*10+k]
	tmp = np.round(df['UH'].values,3)
	for i in [3,2,1]:
		inds = np.argwhere(tmp>i)[:,0]
		tmp[inds]=tmp[inds]-i
	xf,yf=[],[]
	for el in xun_ok:
		inds = np.argwhere(df['UHV'].values==el)[:,0]
		for i in inds:
			xf+=[tmp[i]]
			yf+=[y[i]]

	xf=np.array(xf)
	yf=np.array(yf)

	xfics,yfics=[],{'W':[],'O':[],'WO':[]}
	for el in np.unique(xf):
		inds = np.argwhere(xf==el)[:,0]	
		#print(el,len(inds),np.sum(yf[inds]==1),np.sum(yf[inds]==2),np.sum(yf[inds]==3))
		xfics+=[el]
		yfics['W']  +=[np.sum(yf[inds]==elV['W'])*100./len(inds) ]
		yfics['O']  +=[np.sum(yf[inds]==elV['O']) *100./len(inds) ]
		yfics['WO'] +=[np.sum(yf[inds]==elV['WO']) *100./len(inds) ]

	fileo=open("data/"+title+"_ics.txt",'w')
	fileo.write("x,W,O,WO\n")
	for i in range(len(xfics)):
		fileo.write("%s,%s,%s,%s\n" %(xfics[i],yfics['W'][i],yfics['O'][i],yfics['WO'][i]))
	fileo.close()

	res={'O':0,'WO':0,'W':0}
	breaks=[]
	for key in elV:
		inds = np.argwhere(yf==elV[key])[:,0]
		res[key] = [np.min(xf[inds]),np.max(xf[inds])]
		breaks+=[res[key][0]]
		breaks+=[res[key][1]]

	breaks = np.unique(breaks)
	x1,x2=[],[]
	prx = {'O':[],'W':[],'WO':[]}
	for i in range(1,len(breaks)):
		x1 += [breaks[i-1]]
		x2 += [breaks[i]]
		prx['O']+=[0]
		prx['WO']+=[0]
		prx['W']+=[0]
	for key in elV:
		res[key][0] = res[key][0]+0.0001
		res[key][1] = res[key][1]-0.0001
	for i in range(len(x1)):
		for key in elV:
			if (res[key][0]>=x1[i] and res[key][0]<x2[i]) or  (res[key][1]>=x1[i] and res[key][1]<x2[i]) or (res[key][0]<=x1[i] and res[key][1]>=x2[i]):
				##print("Y",x1[i],x2[i],res[key],key)
				prx[key][i]=1
				
	fileo=open("data/"+title+".txt",'w')	
	fileo.write("X1,X2,Y1,Y2,color\n")
	for i in range(len(x1)):
		tmpL=[]
		for key in prx:
			if prx[key][i]==1:
				tmpL+=[key] 
		if equals(tmpL,['O']):
			colD=0
		elif equals(tmpL,['W']):
			colD=1
		elif equals(tmpL,['WO']):
			colD=2
		elif equals(tmpL,['W','O']):
			colD=3
		elif equals(tmpL,['W','WO']):
			colD=4
		elif equals(tmpL,['O','WO']):
			colD=5
		elif equals(tmpL,['O','W','WO']):
			colD=6
		else:
			colD=-1
		fileo.write("%s,%s,-0.5,0.5,%s\n" %(x1[i],x2[i],colD))
	fileo.close()
 
	fileo = open("data/"+title+"_legend.txt",'w')
	fileo.write('O,W,WO,col\n')
	fileo.write("1,0,0,0\n")
	fileo.write("0,1,0,1\n")
	fileo.write("0,0,1,2\n")
	fileo.write("1,1,0,3\n")
	fileo.write("0,1,1,4\n")
	fileo.write("1,0,1,5\n")
	fileo.write("1,1,1,6\n")
	fileo.write("0,0,0,-1\n")
	fileo.close()

############
############
def getDataSets_uhSingles(title):
	df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/singles_crosstalk_uh10.txt")
	typeList=['E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']
	elV={}
	for i in range(len(typeList)):
		elV[typeList[i]]=i*1.
	y=0
	for i in range(len(typeList)):
		y+=df[typeList[i]].values*elV[typeList[i]]

	xun_ok=[]
	for i in range(4):
		for j in range(1,4):
			for k in range(4):
				if (j==1 and k<=0) or (j==2 and k<=1) or (j==3):
					xun_ok+=[i*100+j*10+k]
	tmp = np.round(df['UH'].values,3)
	for i in [3,2,1]:
		inds = np.argwhere(tmp>i)[:,0]
		tmp[inds]=tmp[inds]-i
	xf,yf=[],[]
	for el in xun_ok:
		inds = np.argwhere(df['UHV'].values==el)[:,0]
		for i in inds:
			xf+=[tmp[i]]
			yf+=[y[i]]

	xf=np.array(xf)
	yf=np.array(yf)

	xfics,yfics=[],{}
	for el in typeList:
		yfics[el]=[]
	for el in np.unique(xf):
		inds = np.argwhere(xf==el)[:,0]	
		#print(el,len(inds),np.sum(yf[inds]==1),np.sum(yf[inds]==2),np.sum(yf[inds]==3))
		xfics+=[el]
		for key in typeList:
			yfics[key]  +=[np.sum(yf[inds]==elV[key])*100./len(inds) ]

	fileo=open("data/"+title+"_ics.txt",'w')
	fileo.write("x")
	for el in typeList:
		fileo.write(",%s" %el)
	fileo.write("\n")

	for i in range(len(xfics)):
		fileo.write("%s" %(xfics[i]))
		for el in typeList:
			fileo.write(",%s" %(yfics[el][i]))
	fileo.write("\n")
	fileo.close()

	res={}
	#for el in typeList:
	#	res[el]=[]
	breaks=[]
	for key in elV:
		inds = np.argwhere(yf==elV[key])[:,0]
		if len(inds)>0:
			res[key] = [np.min(xf[inds]),np.max(xf[inds])]
			breaks+=[res[key][0]]
			breaks+=[res[key][1]]

	breaks = np.unique(breaks)
	x1,x2=[],[]
	prx = {}
	for el in typeList:
		prx[el]=[]
	for i in range(1,len(breaks)):
		x1 += [breaks[i-1]]
		x2 += [breaks[i]]
		for el in typeList:
			prx[el]+=[0]
	for key in res:
		res[key][0] = res[key][0]+0.001
		res[key][1] = res[key][1]-0.001
	for i in range(len(x1)):
		for key in res:
			if (res[key][0]>=x1[i] and res[key][0]<x2[i]) or  (res[key][1]>=x1[i] and res[key][1]<x2[i]) or (res[key][0]<=x1[i] and res[key][1]>=x2[i]):
				##print("Y",x1[i],x2[i],res[key],key)
				prx[key][i]=1
				

	groups=[]
	fileo=open("data/"+title+".txt",'w')	
	fileo.write("X1,X2,Y1,Y2,color\n")
	for i in range(len(x1)):
		tmpL=[]
		for key in prx:
			if prx[key][i]==1:
				tmpL+=[key] 
		if len(groups)==0:
			colD=0
			groups+=[tmpL]
		else:
			flag=False
			for j in range(len(groups)):
				if equals(tmpL,groups[j]):
					colD=j
					flag=True
					break

			if not flag:
				colD=len(groups)
				groups+=[tmpL]

		fileo.write("%s,%s,-0.5,0.5,%s\n" %(x1[i],x2[i],colD))
	fileo.close()
 
	fileo = open("data/"+title+"_legend.txt",'w')
	for el in typeList:
		fileo.write("%s," %el)
	fileo.write('t3\n')
	for i in range(len(groups)):
		for el in typeList:
			if el in groups[i]:
				fileo.write("1,")
			else:
				fileo.write("0,")
		fileo.write("%s\n" %i)
	fileo.close()


############
############
def getDataSets_uhSingles_red(title):
	df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/singles_crosstalk_uh10.txt")
	elV = {'O':1.,'WO':2.,'W':3.}
	y = df['O'].values*elV['O']+df['WO'].values*elV['WO']+df['W'].values*elV['W']
	xun_ok=[]
	for i in range(4):
		for j in range(1,4):
			for k in range(4):
				if (j==1 and k<=0) or (j==2 and k<=1) or (j==3):
					xun_ok+=[i*100+j*10+k]
	tmp = np.round(df['UH'].values,3)
	for i in [3,2,1]:
		inds = np.argwhere(tmp>i)[:,0]
		tmp[inds]=tmp[inds]-i
	xf,yf=[],[]
	for el in xun_ok:
		inds = np.argwhere(df['UHV'].values==el)[:,0]
		for i in inds:
			xf+=[tmp[i]]
			yf+=[y[i]]

	xf=np.array(xf)
	yf=np.array(yf)

	xfics,yfics=[],{'W':[],'O':[],'WO':[]}
	for el in np.unique(xf):
		inds = np.argwhere(xf==el)[:,0]	
		#print(el,len(inds),np.sum(yf[inds]==1),np.sum(yf[inds]==2),np.sum(yf[inds]==3))
		xfics+=[el]
		yfics['W']  +=[np.sum(yf[inds]==elV['W'])*100./len(inds) ]
		yfics['O']  +=[np.sum(yf[inds]==elV['O']) *100./len(inds) ]
		yfics['WO'] +=[np.sum(yf[inds]==elV['WO']) *100./len(inds) ]

	fileo=open("data/"+title+"_ics.txt",'w')
	fileo.write("x,W,O,WO\n")
	for i in range(len(xfics)):
		fileo.write("%s,%s,%s,%s\n" %(xfics[i],yfics['W'][i],yfics['O'][i],yfics['WO'][i]))
	fileo.close()

	res={'O':0,'WO':0,'W':0}
	breaks=[]
	for key in elV:
		inds = np.argwhere(yf==elV[key])[:,0]
		res[key] = [np.min(xf[inds]),np.max(xf[inds])]
		breaks+=[res[key][0]]
		breaks+=[res[key][1]]

	breaks = np.unique(breaks)
	x1,x2=[],[]
	prx = {'O':[],'W':[],'WO':[]}
	for i in range(1,len(breaks)):
		x1 += [breaks[i-1]]
		x2 += [breaks[i]]
		prx['O']+=[0]
		prx['WO']+=[0]
		prx['W']+=[0]
	for key in elV:
		res[key][0] = res[key][0]+0.001
		res[key][1] = res[key][1]-0.001
	for i in range(len(x1)):
		for key in elV:
			if (res[key][0]>=x1[i] and res[key][0]<x2[i]) or  (res[key][1]>=x1[i] and res[key][1]<x2[i]) or (res[key][0]<=x1[i] and res[key][1]>=x2[i]):
				##print("Y",x1[i],x2[i],res[key],key)
				prx[key][i]=1
				
	fileo=open("data/"+title+".txt",'w')	
	fileo.write("X1,X2,Y1,Y2,color\n")
	for i in range(len(x1)):
		tmpL=[]
		for key in prx:
			if prx[key][i]==1:
				tmpL+=[key] 
		if equals(tmpL,['O']):
			colD=0
		elif equals(tmpL,['W']):
			colD=1
		elif equals(tmpL,['WO']):
			colD=2
		elif equals(tmpL,['W','O']):
			colD=3
		elif equals(tmpL,['W','WO']):
			colD=4
		elif equals(tmpL,['O','WO']):
			colD=5
		elif equals(tmpL,['O','W','WO']):
			colD=6
		else:
			colD=-1
		fileo.write("%s,%s,-0.5,0.5,%s\n" %(x1[i],x2[i],colD))
	fileo.close()
 
	fileo = open("data/"+title+"_legend.txt",'w')
	fileo.write('O,W,WO,col\n')
	fileo.write("1,0,0,0\n")
	fileo.write("0,1,0,1\n")
	fileo.write("0,0,1,2\n")
	fileo.write("1,1,0,3\n")
	fileo.write("0,1,1,4\n")
	fileo.write("1,0,1,5\n")
	fileo.write("1,1,1,6\n")
	fileo.write("0,0,0,-1\n")
	fileo.close()

def getDataSets_reduced(namex,title,typeR,namey=None,constants=None,start=None,ftitle=None):
    xval,yval,label=[],[],[]
    hybridV = []
    if namey!=None and constants!=None:
        name=namex+"_"+namey
        for key in constants:
            name+="_"+key
    elif namey!=None:
        name=namex+"_"+namey
    else:
        name=namex

    if start==None:
        start=''
    if ftitle!=None:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+ftitle,header=0)
    else:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+start+"crosstalk_"+name+".txt",header=0)
    nics = float(df['nics'][0])

    if constants!=None:
        keys =list(constants.keys())
        key = keys[0]
        if type(constants[key])!=list:
            inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        else:
            inds = list(np.argwhere(df[key.upper()].values==constants[key][0])[:,0])
            for i in range(1,len(constants[key])):
                inds +=list(np.argwhere(df[key.upper()].values==constants[key][i])[:,0])
            inds = np.array(inds)
        for i in range(1,len(keys)):
           if type(constants[keys[i]])!=list:
                tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
                inds = np.intersect1d(inds,tmpI)
           else:
                tmpI=[]
                for j in range(len(constants[keys[i]])):
                    tmpI += list(np.argwhere(df[keys[i].upper()].values==constants[keys[i]][j])[:,0])
                inds = np.intersect1d(inds,np.array(tmpI))
    else:
        inds=np.argwhere(df['nics'].values>-1)[:,0]


    for i in range(len(df.values[inds])):
        if 'UH' in namex.upper():
            xval+=[df['UH'][inds[i]]]
        elif 'HU' in namex.upper():
            xval+=[df['HU'][inds[i]]]
        else:
            xval+=[df[namex.upper()][inds[i]]]
        if namey!=None:
            if 'UH' in namey.upper():
                yval+=[df['UH'][inds[i]]]
            elif 'HU' in namey.upper():
                yval+=[df['HU'][inds[i]]]
            else:
                yval+=[df[namey.upper()][inds[i]]]
        else:
                yval+=[0]
        hybridV+=[df['EM/WO'][inds[i]]/nics*100.]

    compHybrid=getCompHybridVal(compI=False)
    labelTmp=getStateListfromFile(df)

    label=[]
    for i in range(len(inds)):
       label+=[labelTmp[inds[i]]]

    if 1 not in xval and namey==None and start==None:
        xval+=[1]
        yval+=[0]
        label+=[['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']]
        hybridV+=[compHybrid]

    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    label= np.array(label)[inds]
    hybridV= np.array(hybridV)[inds]

    x1,x2,y1,y2,label,hybridV=getDataForPlot(xval,yval,label,hybridV)



    fileo = open("data/"+title+".txt",'w')
    fileo.write("X1,X2,Y1,Y2,color\n")
    for i in range(len(x1)):
         if hybridV[i]==100:
            tmpC=3
         elif hybridV[i]==compHybrid:
            tmpC=2
         elif hybridV[i]>compHybrid:
            tmpC=1
         else:#hybridV[i]<compHybrid:
            tmpC=0
         tmpL=[]
         for j in range(len(label[i])):
              if typeR=='MR':
                  tmp2 = label[i][j].split("/")[-1]
              else:
                  tmp2 = label[i][j].split("/")[0]
              if tmp2 not in tmpL:
                  tmpL+=[tmp2] 

         if typeR=='MR':
                if equals(tmpL,['O']):
    	       	     colD=0
                elif equals(tmpL,['W']):
    	       	     colD=1
                elif equals(tmpL,['WO']):
   	       	     colD=2
                elif equals(tmpL,['W','O']):
    	       	     colD=3
                elif equals(tmpL,['W','WO']):
    	       	     colD=4
                elif equals(tmpL,['O','WO']):
                     colD=5
                elif equals(tmpL,['O','W','WO']):
                   colD=6
                else:
                   colD=-1
         else:
            if equals(tmpL,['E']):
                colD=0
            elif equals(tmpL,['M']):
                colD=1
            elif equals(tmpL,['EM']):
                colD=2
            elif equals(tmpL,['M','E']):
                colD=3
            elif equals(tmpL,['M','EM']):
                colD=4
            elif equals(tmpL,['E','EM']):
                colD=5
            elif equals(tmpL,['E','M','EM']):
                colD=6
            else:
                colD=-1
         fileo.write("%s,%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i],colD))
    fileo.close()

    if typeR=='MR':
        fileo = open("data/"+title+"_legend.txt",'w')
        fileo.write('O,W,WO,col\n')
        fileo.write("1,0,0,0\n")
        fileo.write("0,1,0,1\n")
        fileo.write("0,0,1,2\n")
        fileo.write("1,1,0,3\n")
        fileo.write("0,1,1,4\n")
        fileo.write("1,0,1,5\n")
        fileo.write("1,1,1,6\n")
        fileo.write("0,0,0,-1\n")
        fileo.close()
    else:
        fileo = open("data/"+title+"_legend.txt",'w')
        fileo.write('E,M,EM,col\n')
        fileo.write("1,0,0,0\n")
        fileo.write("0,1,0,1\n")
        fileo.write("0,0,1,2\n")
        fileo.write("1,1,0,3\n")
        fileo.write("0,1,1,4\n")
        fileo.write("1,0,1,5\n")
        fileo.write("1,1,1,6\n")
        fileo.write("0,0,0,-1\n")
        fileo.close()

############
############
def getDataSets(namex,title,namey=None,constants=None,start=None,Ftitle=None):
    xval,yval,label=[],[],[]
    hybridV = []
    if namey!=None and constants!=None:
        name=namex+"_"+namey
        for key in constants:
            name+="_"+key
    elif namey!=None:
        name=namex+"_"+namey
    else:
        name=namex

    if start==None:
       start=''
    
    if Ftitle==None:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+start+"crosstalk_"+name+".txt",header=0)
    else:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+Ftitle,header=0)
    nics = float(df['nics'][0])

    if constants!=None:
        keys =list(constants.keys())
        key = keys[0]
        if type(constants[key])!=list:
            inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        else:
            inds = list(np.argwhere(df[key.upper()].values==constants[key][0])[:,0])
            for i in range(1,len(constants[key])):
                inds +=list(np.argwhere(df[key.upper()].values==constants[key][i])[:,0])
            inds = np.array(inds)
        for i in range(1,len(keys)):
           if type(constants[keys[i]])!=list:
                tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
                inds = np.intersect1d(inds,tmpI)
           else:
                tmpI=[]
                for j in range(len(constants[keys[i]])):
                    tmpI += list(np.argwhere(df[keys[i].upper()].values==constants[keys[i]][j])[:,0])
                inds = np.intersect1d(inds,np.array(tmpI))
    else:
        inds=np.argwhere(df['nics'].values>-1)[:,0]

    ## acutal values are corred
    for i in range(len(df.values[inds])):
        if 'UH' in namex.upper():
                xval+=[df['UH'][inds[i]]]
        elif 'HU' in namex.upper():
                xval+=[df['HU'][inds[i]]]
        else:
                xval+=[df[namex.upper()][inds[i]]]
        if namey!=None:
            if 'UH' in namey.upper():
                yval+=[df['UH'][inds[i]]]
            elif 'HU' in namey.upper():
                yval+=[df['HU'][inds[i]]]
            else:
                yval+=[df[namey.upper()][inds[i]]]
        else:
            yval+=[0]
        hybridV+=[df['EM/WO'][inds[i]]*1./nics*100.]


    compHybrid=getCompHybridVal(compI=False)
    labelTmp=getStateListfromFile(df)

    label=[]
    for i in range(len(inds)):
         label+=[labelTmp[inds[i]]]


    if 1 not in xval and namey==None and start==None:
        xval+=[1]
        yval+=[0]
        label+=[['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']]
        hybridV+=[compHybrid]

    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    label= np.array(label)[inds]
    hybridV= np.array(hybridV)[inds]


    x1,x2,y1,y2,label,hybridV=getDataForPlot(xval,yval,label,hybridV)
    color,star_colors,colList = getPlotData(label)
    hatch,edgcolors,regType=getHatchForPlot(star_colors,color,hybridV,compHybrid,colList)

    tmp = np.unique(color)
    colD={}
    count=0
    for el in tmp:
        colD[el]=count
        count+=1

    fileo = open("data/"+title+".txt",'w')
    fileo.write("X1,X2,Y1,Y2,color,ecolor,hatch,colorHHReg\n")
    for i in range(len(x1)):
         if hybridV[i]==100:
              tmpC=3
         elif hybridV[i]==compHybrid:
              tmpC=2
         elif hybridV[i]>compHybrid:
              tmpC=1
         else:#hybridV[i]<compHybrid:
             tmpC=0
         fileo.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i],colD[color[i]],edgcolors[i],hatch[i],tmpC))
    fileo.close()

    fileo = open("data/"+title+"_legend.txt",'w')
    fileo.write('E/O,E/W,E/WO,EM/O,EM/W,EM/WO,M/O,M/W,M/WO,E,EM,M,O,WO,W,t1,t2,t3\n')
    for i in range(len(regType)):
        tmpD={'E':0,'M':0,'EM':0,'W':0,'WO':0,'O':0}
        if 'E/O' in regType[i][0]:
             fileo.write("1")
             tmpD['E']=1
             tmpD['O']=1
        else:
             fileo.write("0")
        for key in ['E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
            if key in regType[i][0]:
                 fileo.write(",1")
                 for k2 in tmpD:
                      if k2==key.split("/")[0] or k2==key.split("/")[1]:
                          tmpD[k2]=1
            else:
                 fileo.write(",0")
	
        for key in ['E','EM','M','O','WO','W']:
            fileo.write(",%s" %tmpD[key])

        fileo.write(",%s,%s,%s\n" %(regType[i][1],regType[i][2],colD[regType[i][3]]))
    fileo.close()
############
def getDataSets_EMWO(namex,title,alphaVal,namey=None,constants=None,start=None,Ftitle=None):
    xval,yval,label=[],[],[]
    hybridV = []
    if namey!=None and constants!=None:
        name=namex+"_"+namey
        for key in constants:
           name+="_"+key
    elif namey!=None:
        name=namex+"_"+namey
    else:
         name=namex

    if start==None:
        start=''
    
    if Ftitle==None:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+start+"crosstalk_"+name+".txt",header=0)
    else:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+Ftitle,header=0)
    nics = float(df['nics'][0])

    if constants!=None:
        keys =list(constants.keys())
        key = keys[0]
        if type(constants[key])!=list:
            inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        else:
            inds = list(np.argwhere(df[key.upper()].values==constants[key][0])[:,0])
            for i in range(1,len(constants[key])):
                inds +=list(np.argwhere(df[key.upper()].values==constants[key][i])[:,0])
            inds = np.array(inds)
        for i in range(1,len(keys)):
           if type(constants[keys[i]])!=list:
                tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
                inds = np.intersect1d(inds,tmpI)
           else:
                tmpI=[]
                for j in range(len(constants[keys[i]])):
                    tmpI += list(np.argwhere(df[keys[i].upper()].values==constants[keys[i]][j])[:,0])
                inds = np.intersect1d(inds,np.array(tmpI))
    else:
        inds=np.argwhere(df['nics'].values>-1)[:,0]

    ## acutal values are corred
    for i in range(len(df.values[inds])):
        if 'UH' in namex.upper():
                xval+=[df['UH'][inds[i]]]
        elif 'HU' in namex.upper():
                xval+=[df['HU'][inds[i]]]
        else:
                xval+=[df[namex.upper()][inds[i]]]
        if namey!=None:
            if 'UH' in namey.upper():
               yval+=[df['UH'][inds[i]]]
            elif 'HU' in namey.upper():
               yval+=[df['HU'][inds[i]]]
            else:
               yval+=[df[namey.upper()][inds[i]]]
        else:
               yval+=[0]
        hybridV+=[df['EM/WO'][inds[i]]*1./nics*100.]


    compHybrid=getCompHybridVal(compI=False)

    if 1 not in xval and namey==None and start==None:
         xval+=[1]
         yval+=[0]
         label+=[['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']]
         hybridV+=[compHybrid]

    fileoE = open("data/"+title+"_EMWO.txt",'w')
    fileoE.write("X,Y,EMWO,Change\n")
    for i in range(len(xval)):
         tmp = (hybridV[i]/compHybrid)
         fileoE.write("%s,%s,%s,%s\n" %(xval[i],yval[i],hybridV[i],tmp))
    fileoE.close()
    
############
############
############
############
def getDataSets_HHLoc(namex,title,alphaVal,namey=None,constants=None,start=None,Ftitle=None):
    xval,yval,label=[],[],[]
    hybridV = []
    if namey!=None and constants!=None:
        name=namex+"_"+namey
        for key in constants:
           name+="_"+key
    elif namey!=None:
        name=namex+"_"+namey
    else:
         name=namex

    if start==None:
        start=''
    
    if Ftitle==None:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+start+"crosstalk_"+name+".txt",header=0)
    else:
        df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+Ftitle,header=0)
    nics = float(df['nics'][0])

    if constants!=None:
        keys =list(constants.keys())
        key = keys[0]
        if type(constants[key])!=list:
            inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        else:
            inds = list(np.argwhere(df[key.upper()].values==constants[key][0])[:,0])
            for i in range(1,len(constants[key])):
                inds +=list(np.argwhere(df[key.upper()].values==constants[key][i])[:,0])
            inds = np.array(inds)
        for i in range(1,len(keys)):
           if type(constants[keys[i]])!=list:
                tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
                inds = np.intersect1d(inds,tmpI)
           else:
                tmpI=[]
                for j in range(len(constants[keys[i]])):
                    tmpI += list(np.argwhere(df[keys[i].upper()].values==constants[keys[i]][j])[:,0])
                inds = np.intersect1d(inds,np.array(tmpI))
    else:
        inds=np.argwhere(df['nics'].values>-1)[:,0]

    ## acutal values are corred
    for i in range(len(df.values[inds])):
        if 'UH' in namex.upper():
                xval+=[df['UH'][inds[i]]]
        elif 'HU' in namex.upper():
                xval+=[df['HU'][inds[i]]]
        else:
                xval+=[df[namex.upper()][inds[i]]]
        if namey!=None:
            if 'UH' in namey.upper():
               yval+=[df['UH'][inds[i]]]
            elif 'HU' in namey.upper():
               yval+=[df['HU'][inds[i]]]
            else:
               yval+=[df[namey.upper()][inds[i]]]
        else:
               yval+=[0]
        hybridV+=[df['EM/WO'][inds[i]]*1./nics*100.]


    compHybrid=getCompHybridVal(compI=False)
    labelTmp=getStateListfromFile(df)

    label=[]
    for i in range(len(inds)):
       label+=[labelTmp[inds[i]]]

    if 1 not in xval and namey==None and start==None:
         xval+=[1]
         yval+=[0]
         label+=[['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']]
         hybridV+=[compHybrid]

    inds = np.argsort(xval)
    xval = np.array(xval)[inds]
    yval = np.array(yval)[inds]
    label= np.array(label)[inds]
    hybridV= np.array(hybridV)[inds]

    x1,x2,y1,y2,label,hybridV=getDataForPlot(xval,yval,label,hybridV)


    xmin,xmax,ymin,ymax=[],[],[],[]
    xminO,xmaxO,yminO,ymaxO=[],[],[],[]
    xt,yt=[],[]
    xt0,yt0=[],[]
    fileoE = open("data/"+title+"_HHexists.txt",'w')
    fileoO = open("data/"+title+"_HHOnlhy.txt",'w')
    fileoE.write("X1,X2,Y1,Y2\n")
    fileoO.write("X1,X2,Y1,Y2\n")
    for i in range(len(x1)):
        if hybridV[i]!=100 and hybridV[i]>0:
                fileoE.write("%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i]))
                xmin+=[x1[i]]
                ymin+=[y1[i]]
                xmax+=[x2[i]]
                ymax+=[y2[i]]
                xt+=[x1[i],x2[i]]#,xval[i]]
                yt+=[y1[i],y2[i]]#,yval[i]]
		
        elif hybridV[i]==100:
                 fileoO.write("%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i]))
                 xminO+=[x1[i]]
                 yminO+=[y1[i]]
                 xmaxO+=[x2[i]]
                 ymaxO+=[y2[i]]
                 xt0+=[x1[i],x2[i]]#,xval[i]]
                 yt0+=[y1[i],y2[i]]#,yval[i]]
    fileoE.close()
    fileoO.close()
    
    if len(xt)>0:
        points = np.zeros((len(xt),2))
        points[:,0]=xt
        points[:,1]=yt
        ashape = alphashape.alphashape(points,alphaVal)
        try:
            exterior = list(ashape.exterior.coords)
        
            fileoH = open("data/"+title+"_hull.txt",'w')
            fileoH.write("x,y\n")
            for i in range(len(exterior)):
                    fileoH.write("%s,%s\n" %(exterior[i][0],exterior[i][1]))
            fileoH.close()
        except:
           print("Error in getting hull",len(xt),ashape)


    if len(xt0)>0:
        points = np.zeros((len(xt0),2))
        points[:,0]=xt0
        points[:,1]=yt0
        ashape = alphashape.alphashape(points,9.)
        try:
            exterior = list(ashape.exterior.coords)
        
            fileoH = open("data/"+title+"_hullO.txt",'w')
            fileoH.write("x,y\n")
            for i in range(len(exterior)):
                    fileoH.write("%s,%s\n" %(exterior[i][0],exterior[i][1]))
            fileoH.close()
        except:
           print("Error in getting hull",len(xt),ashape)

############
############
############
def getData_6():
	getDataSets('Au','data_6c3',namey='HS',start='PSF_')
	getDataSets_uhSingles('data_3_uh')
	getDataSets_ics('u3m','data_3_u3m','u3m')
	getDataSets_reduced('u3n','data_3_MR1_red','MR',namey='u3m',constants={'UHV':300},ftitle='crosstalk_uh_u3n_u3m.txt')
	getDataSets_HHLoc('u3n','data_3_MR1_red',6,namey='u3m',constants={'UHV':300},Ftitle='crosstalk_uh_u3n_u3m.txt')

##################
##################
##################
def getData_1c():
    fileo = open("data/data_1c.txt","w")
    fileo.write("icsNumber")
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
        fileo.write(",%sAvg,%sStd" %(key,key))
    fileo.write("\n")

    keys=[]
    res={}
    for x in [100,500,1000,2000]:#,5000,10000]: 
        if x not in res.keys():
            res[x]={'E/O':[],'E/W':[],'E/WO':[],'EM/O':[],'EM/W':[],'EM/WO':[],'M/O':[],'M/W':[],'M/WO':[],'E':[],'EM':[],'M':[],'O':[],'WO':[],'W':[]}

        for i in range(10):
                df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_"+str(i)+"_"+str(x)+"_res.txt").dropna()

                mapRes = getStates(df)
                for key in mapRes:
                    res[x][key] +=[mapRes[key]/float(x)*100.]

        fileo.write("%s" %x)
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
             fileo.write(",%s,%s" %(np.mean(res[x][key]),np.std(res[x][key])))
        fileo.write("\n")
    fileo.close()

    ####################3
    fileo = open("data/data_1d.txt","w")
    fileo.write("icsNumber")
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
        fileo.write(",%sAvg,%sStd" %(key,key))
    fileo.write("\n")

    keys=[]
    res={}
    for x in [100,500,1000,2000]:#,5000,10000]: 
        if x not in res.keys():
            res[x]={'E/O':[],'E/W':[],'E/WO':[],'EM/O':[],'EM/W':[],'EM/WO':[],'M/O':[],'M/W':[],'M/WO':[],'E':[],'EM':[],'M':[],'O':[],'WO':[],'W':[]}

        for i in range(10):
                df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/comparison_coupled/EMT_MR10_AS_0_95_AZ_0_95_Hu_0_9_u3m_0_7_u3n_0_7_HS_1_1_Au_1_01_"+str(i)+"_"+str(x)+"_res.txt").dropna()

                mapRes = getStates(df)
                for key in mapRes:
                    res[x][key] +=[mapRes[key]/float(x)*100.]

        fileo.write("%s" %x)
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
             fileo.write(",%s,%s" %(np.mean(res[x][key]),np.std(res[x][key])))
        fileo.write("\n")
    fileo.close()



##################
def getData_NOHHICS():
    fileo = open("data/data_noHHics.txt","w")
    fileo.write("icsNumber")
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
        fileo.write(",%sAvg,%sStd" %(key,key))
    fileo.write("\n")

    keys=[]
    res={}
    for x in [100,500,1000,2000]:#,5000,10000]: 
        if x not in res.keys():
            res[x]={'E/O':[],'E/W':[],'E/WO':[],'EM/O':[],'EM/W':[],'EM/WO':[],'M/O':[],'M/W':[],'M/WO':[],'E':[],'EM':[],'M':[],'O':[],'WO':[],'W':[]}

        for i in range(4):
                df =pd.read_csv("~/Research/EMT_MR/crosstalk/normalCells_coupled/bothModified/crosstalk_comparison/EMT_MR_comp_"+str(i)+"_"+str(x)+"_res.txt").dropna()

                mapRes = getStates(df)
                for key in mapRes:
                    res[x][key] +=[mapRes[key]/float(x)*100.]

        fileo.write("%s" %x)
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
             fileo.write(",%s,%s" %(np.mean(res[x][key]),np.std(res[x][key])))
        fileo.write("\n")
    fileo.close()
##########################
#
##################
def getData_PSF_ICS():
    fileo = open("data/data_PSF_ics.txt","w")
    fileo.write("icsNumber")
    for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
        fileo.write(",%sAvg,%sStd" %(key,key))
    fileo.write("\n")

    keys=[]
    res={}
    for x in [200,500,1000,2000]:#,5000,10000]: 
        if x not in res.keys():
            res[x]={'E/O':[],'E/W':[],'E/WO':[],'EM/O':[],'EM/W':[],'EM/WO':[],'M/O':[],'M/W':[],'M/WO':[],'E':[],'EM':[],'M':[],'O':[],'WO':[],'W':[]}

        for i in range(10):
                df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/PSF/crosstalk_comparison/EMT_MR_"+str(i)+"_"+str(x)+"_res.txt").dropna()
                mapRes = getStates(df,PSF=True)
                for key in mapRes:
                    res[x][key] +=[mapRes[key]/float(x)*100.]

        fileo.write("%s" %x)
        for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']:
             fileo.write(",%s,%s" %(np.mean(res[x][key]),np.std(res[x][key])))
        fileo.write("\n")
    fileo.close()
##########################
##########################
##########################
##########################
def getParamData():

	###  code in ~/Research/EMT_MR/understanding/checks/metabolism_HMRNA/final_par
	### and ~/Research/EMT_MR/understanding/checks/metabolism
	print("x")


###########################################
###########################################
###########################################
###########################################
def getUH_correct(filename,name,xlab):
        ### y is UH
        ### x is other
        df = pd.read_csv("/home/madeline/Research/EMT_MR/crosstalk/analysis/data/"+filename)
        xun_ok=[]
        for i in range(4):
               for j in range(1,4):
                       for k in range(4):
                             if (j==1 and k<=0) or (j==2 and k<=1) or (j==3):
                                      xun_ok+=[i*100+j*10+k]
        inds = []
        for el in xun_ok:
            inds += list(np.argwhere(df['UHV'].values==el)[:,0])

        usefulData = np.round(df.values[inds],4)
        states = np.unique(usefulData[:,16:-1],axis=0)

        tags = { 0:'M/W' ,1:'M/WO', 2:'M/O',
                3:'EM/W', 4:'EM/WO', 5:'EM/O',
                6:'E/W', 7:'E/WO', 8:'E/O'}
        results = { 'M/W':{'x':[],'y':[]} ,'M/WO':{'x':[],'y':[]}, 'M/O':{'x':[],'y':[]},
           'EM/W':{'x':[],'y':[]}, 'EM/WO':{'x':[],'y':[]}, 'EM/O':{'x':[],'y':[]},
           'E/W':{'x':[],'y':[]}, 'E/WO':{'x':[],'y':[]}, 'E/O':{'x':[],'y':[]}}


        eps = 0.00001
        xbreaks,ybreaks=[],[]
        for el in states:
            for i in range(len(el)):
                tmp = usefulData[:,16:-1]
                if i==0:
                    locs = np.argwhere(el[i]==tmp[:,i])[:,0]
                else:
                    locs2 = np.argwhere(el[i]==tmp[:,i])[:,0]
                    locs = np.intersect1d(locs,locs2)

            labels={'AS':0, 'AZ':1, 'AU':2, 'HS':3, 'HU':4, 'INPUT':5, 'U3M':6, 'U3N':7, 'UH':8, 'UHV':9}

            if xlab:
                        xv = usefulData[locs,labels[xlab]]
                        yv = usefulData[locs,8]        
            else:
                        xv = usefulData[locs,8]
                        yv = (usefulData[locs,8])*0.
    
            sloc = np.argwhere(el==1)[:,0]
    
            for i in sloc:
                ss = tags[i]
                results[ss]['x'] +=[np.min(xv),np.max(xv)]
                xbreaks +=[np.min(xv),np.max(xv)]
                if xlab:
                     results[ss]['y'] +=[np.min(yv),np.max(yv)]
                     ybreaks +=[np.min(yv),np.max(yv)]
                else:
                     results[ss]['y'] +=[-0.5,0.5]
                     ybreaks +=[-0.5,0.5]
 

        xbreaks = np.sort(np.unique(xbreaks))
        ybreaks = np.sort(np.unique(ybreaks))
        
        x1,x2,y1,y2=[],[],[],[]
        finalRes=[]
        for i in range(len(xbreaks)-1):
            for j in range(len(ybreaks)-1):

                tmpR = []
                for el in results:
                    xlow,xhigh,ylow,yhigh = False,False,False,False
                    if i<len(xbreaks)-1 and len(results[el]['x'])>0:
                        xlow = results[el]['x'][0]<xbreaks[i]+eps
                        xhigh = results[el]['x'][1]>xbreaks[i+1]-eps
                    elif   len(results[el]['x'])>0:
                        xlow = results[el]['x'][0]<xbreaks[i]+eps
                        xhigh=True
                    if j<len(ybreaks)-1 and len(results[el]['y'])>0:
                        ylow = results[el]['y'][0]<ybreaks[j]+eps
                        yhigh = results[el]['y'][1]>ybreaks[j+1]-eps
                    elif len(results[el]['y'])>0:
                        ylow = results[el]['y'][0]<ybreaks[j]+eps
                        yhigh=True
                

                    #print xlow,xhigh,results[el]['x']
                    #print ylow,yhigh,results[el]['y']     ,'\n'     
                    if xlow and xhigh and ylow and yhigh:
                        tmpR+=[el]
                if len(tmpR)>0:# and xlow and xhigh and ylow and yhigh:
                        x1+=[xbreaks[i]]
                        x2+=[xbreaks[i+1]]                
                        y1+=[ybreaks[j]]                
                        y2+=[ybreaks[j+1]]
                        finalRes+=[tmpR]
        colLst=list(np.unique(finalRes))

        fileo = open("data/"+name+".txt",'w')

        fileo.write("X1,X2,Y1,Y2,color\n")

        for i in range(len(x1)):

            for j in range(len(colLst)):
                if equals(finalRes[i],colLst[j]):
                    colN = j+1
            fileo.write("%s,%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i],colN))

        fileo.close()
        #####################33
        fileo = open("data/"+name+"_HHOnlhy.txt",'w')

        fileo.write("X1,X2,Y1,Y2\n")

        for i in range(len(x1)):
            if equals(finalRes[i],['EM/WO']):
                fileo.write("%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i]))

        fileo.close()
        #############################
        fileo = open("data/"+name+"_HHexists.txt",'w')

        fileo.write("X1,X2,Y1,Y2\n")

        for i in range(len(x1)):
            if 'EM/WO' in finalRes[i]:
                fileo.write("%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i]))

        fileo.close()
        #############################
        fileo = open("data/"+name+"_legend.txt",'w')
        fileo.write('E/O,E/W,E/WO,EM/O,EM/W,EM/WO,M/O,M/W,M/WO,E,EM,M,O,WO,W,t3\n')
        
        for i in range(len(colLst)):
            tmp=[]
            for el in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
                if el in colLst[i]:
                    fileo.write("1,")
                    tmp+=el.split("/")
                else:
                    fileo.write("0,")            
            for el in ['E','EM','M','O','WO','W']:
                if el in tmp:
                    fileo.write("1,")
                else:
                    fileo.write("0,")
            fileo.write("%s\n" %(i+1))
        fileo.close()
############
def getData_Fig4_uh_u3n_u3m():
    df = pd.read_csv("/home/madeline/Research/EMT_MR/crosstalk/analysis/data/IND_crosstalk_uh_u3n_u3m.txt")
    typeList=['E/O','E/WO','E/W','EM/O','EM/WO','EM/W','M/O','M/WO','M/W']
    elV={}
    for i in range(len(typeList)):
                    elV[typeList[i]]=i*1.
    z=0
    for i in range(len(typeList)):
                    z+=df[typeList[i]].values*elV[typeList[i]]

    ### first make sure you only use parameters confirming to biology
    xun_ok=[]
    for i in range(7):
        for j in range(5):
            for k in range(6):
                if (j==0 and k==4) or (j==1 and (k==0 or k==4)) or (j==2 and (k<=2 or k==4)) or (j==3 and k<=5) or (j==4):
                          xun_ok+=[i*100+j*10+k]

    print(xun_ok)
    ##### Get rid of silencing shifts you previously generated
    tmp = np.round(df['UH'].values,3)
    tmpy = df['U3N'].values
    for i in [9,8,7,6,5,4,3,2,1]:
         inds = np.argwhere(tmp>i)[:,0]
         if len(inds)>0:
               tmp[inds]=tmp[inds]-i

    #### Only keep if the value of U3M=0 and the UH parameters are correct
    ### x = uh, y = u3n, z = coupled state
    xf,yf,zf=[],[],[]
    inds2 = np.argwhere(df['U3M'].values==0)[:,0]
    for i in inds2:#range(len(df['UHV'].values)):
            if (str(df['UHV'].values[i])=='xxx'):
                xf+=[tmp[i]]
                yf+=[tmpy[i]]
                zf+=[z[i]]
            else:
                 df['UHV'].values[i]=int(df['UHV'].values[i])
    for el in xun_ok:
                    inds = np.argwhere(df['UHV'].values==el)[:,0]
                    inds = np.intersect1d(inds,inds2)
                    for i in inds:
                            xf+=[tmp[i]]
                            yf+=[tmpy[i]]
                            zf+=[z[i]]
    xf=np.array(xf)
    xf = 1.-xf
    yf=np.array(yf)
    yf = 1.-yf
    zf=np.array(zf)

    ## for each value of the x and y figure out what the states are possible
    xfics,yfics,zfics=[],[],{}
    for el in typeList:
                    zfics[el]=[]
    for el in np.unique(xf):
                    inds = np.argwhere(xf==el)[:,0]
                    ytmp = yf[inds]
                    ztmp = zf[inds]
                    for el2 in np.unique(ytmp):
                        inds2 = np.argwhere(ytmp==el2)[:,0]
                        xfics+=[el]
                        yfics+=[el2]
                        for key in typeList:
                                zfics[key]  +=[np.sum(ztmp[inds2]==elV[key])*100./len(inds2) ]


    testZ =[]
    for i in range(len(zfics['E/W'])):
        tmp=[]
        for el in zfics:
            if zfics[el][i]>0:
                tmp+=[el]
        testZ+=[tmp]
    #################################################3
    #################################################3

    count=0
    eps=0.01
    unqXlist = np.unique(np.round(xfics,3) )
    flag=True
    hashDict={}
    while flag and count<10:
            flag=False
            newset=[]
            count+=1
            for i in range(len(unqXlist)):
                added=False
                for j in  range(len(unqXlist)):
                    if i!=j:
                        if np.abs(unqXlist[i]-unqXlist[j])<eps:
                            tmpV = np.round(np.mean([unqXlist[i],unqXlist[j]]),2)
                            if tmpV not in hashDict.keys():
                                hashDict[tmpV]=[unqXlist[j]]
                                hashDict[tmpV]=[unqXlist[i]]
                            else:
                                hashDict[tmpV]+=[unqXlist[i]]
                                hashDict[tmpV]+=[unqXlist[j]]
                            newset+=[tmpV]
                            flag=True
                            added=True
                            break
                if not added:
                    if unqXlist[i] not in hashDict:
                        hashDict[unqXlist[i]]=[unqXlist[i]]
                    else:
                        hashDict[unqXlist[i]]+=[unqXlist[i]]
                    newset+=[unqXlist[i]]
            unqXlist =list(np.unique(newset))

    xfics=np.array(xfics)
    yfics = np.array(yfics)
    for el in hashDict:
        for key in np.unique(hashDict[el]):
            inds = np.argwhere(xfics==key)[:,0]
            xfics[inds] = el

    finX,finY,finZ=[],[],[]
    for el in np.unique(xfics):
        inds = np.argwhere(xfics==el)[:,0]
        tmpy = yfics[inds]
        tmpz = np.array(testZ)[inds]
        for el2 in np.unique(tmpy):
            inds2 = np.argwhere(tmpy==el2)[:,0]
            tmpV=[]
            for i in inds2:
                tmpV+=tmpz[i]
            finX+=[el]
            finY+=[el2]
            finZ+=[tmpV]

    #################################################3
    #################################################3
    #################################################3
    x1,x2,y1,y2,label=getDataForPlot(finX,finY,finZ)
    color,star_colors,colList = getPlotData(label)
    compHybrid=30
    hybridV = np.ones(len(x1))
    hatch,edgcolors,regType=getHatchForPlot(star_colors,color,hybridV,compHybrid,colList)

    tmp = np.unique(color)
    colD={}
    count=0
    for el in tmp:
        colD[el]=count
        count+=1

    title = "data/data_3_MR3_single"
    fileo = open(title+".txt",'w')
    fileo.write("X1,X2,Y1,Y2,color\n")
    for i in range(len(x1)):
         fileo.write("%s,%s,%s,%s,%s\n" %(y1[i],y2[i],x1[i],x2[i],colD[color[i]]))
    fileo.close()

    fileo = open(title+"_legend.txt",'w')
    fileo.write('E/O,E/W,E/WO,EM/O,EM/W,EM/WO,M/O,M/W,M/WO,E,EM,M,O,WO,W,t3\n')
    for i in range(len(regType)):
        tmpD={'E':0,'M':0,'EM':0,'W':0,'WO':0,'O':0}
        if 'E/O' in regType[i][0]:
             fileo.write("1")
             tmpD['E']=1
             tmpD['O']=1
        else:
             fileo.write("0")
        for key in ['E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO']:
            if key in regType[i][0]:
                 fileo.write(",1")
                 for k2 in tmpD:
                      if k2==key.split("/")[0] or k2==key.split("/")[1]:
                          tmpD[k2]=1
            else:
                 fileo.write(",0")
	
        for key in ['E','EM','M','O','WO','W']:
            fileo.write(",%s" %tmpD[key])

        fileo.write(",%s\n" %(colD[regType[i][3]]))
    fileo.close()
    ####################################
    ####################################

    fileoE = open(title+"_HHexists.txt",'w')
    fileoE.write("X1,X2,Y1,Y2\n")
    for i in range(len(x1)):
        if 'EM/WO' in label[i]:
            fileoE.write("%s,%s,%s,%s\n" %(y1[i],y2[i],x1[i],x2[i]))
    fileoE.close()
    fileoO = open(title+"_HHOnlhy.txt",'w')
    fileoO.write("X1,X2,Y1,Y2\n")
    fileoO.close()

############ 
############
#
############ 
# 
##########################
##########################
##########################
##########################



#############################
#### Fig S1
#############################
## all in jupyter notebook

#############################
#### Fig S2
#############################
getData_1c()

#############################
#### Fig S3
#############################
getData_1c()

#############################
#### Fig S4
#############################
CompData("data_s1","comparison_coupled")

#############################
#### Fig S5
#############################
getDataSets('u3m','data_s2a')
getDataSets_ics('u3m','data_s2a','u3m',title='crosstalk_u3m.txt')

#############################
#### Fig S6
#############################
### use data that was in the figure (used to be figure) need to tranfer and fix

#############################
#### Fig S7
#############################
getUH_correct("crosstalk_uh_u3m.txt","data_uh_u3m",'U3M')##S7
getUH_correct("crosstalk_uh_u3n.txt","data_uh_u3n",'U3N')## S7

#############################
#### Fig S8
#############################
getData_Fig4_uh_u3n_u3m()

#############################
#### Fig S9
#############################
getDataSets('HS','data_s2g')
getDataSets_ics('HS','data_s2g','HS',title='crosstalk_HS.txt')

#############################
#### Fig S10
#############################
getDataSets('Hu','data_s2f',Ftitle='crosstalk_iHu.txt')
getDataSets_ics('Hu','data_s2f','Hu',title='crosstalk_iHu.txt')

#############################
#### Fig S11
#############################
getDataSets('input','data_s2h')
getDataSets_ics('input','data_s2h','input',title='crosstalk_input.txt')

#############################
#### Fig S12
#############################
getDataSets('Au','data_s2e')
getDataSets_ics('Au','data_s2e','Au',title='crosstalk_Au.txt')

#############################
#### Fig S13
#############################
getDataSets('AS','data_s2c')
getDataSets_ics('AS','data_s2c','AS',title='crosstalk_AS.txt')

#############################
#### Fig S14
#############################
getDataSets('AZ','data_s2d')
getDataSets_ics('AZ','data_s2d','AZ',title='crosstalk_AZ.txt')

#############################
#### Fig S15
#############################
getDataSets('AS','data_AS_HS',namey='HS')
getDataSets_HHLoc('AS','data_AS_HS',8,namey='HS')
getDataSets_EMWO('AS','data_AS_HS',8,namey='HS')

getDataSets('Au','data_Au_HS',namey='HS')
getDataSets_HHLoc('Au','data_Au_HS',8,namey='HS')
getDataSets_EMWO('Au','data_Au_HS',8,namey='HS')

getDataSets('AS','data_AS_Hu',namey='Hu')
getDataSets_HHLoc('AS','data_AS_Hu',8,namey='Hu')
getDataSets_EMWO('AS','data_AS_Hu',8,namey='Hu')

getDataSets('Au','data_Au_Hu',namey='Hu')
getDataSets_HHLoc('Au','data_Au_Hu',8,namey='Hu')
getDataSets_EMWO('Au','data_Au_Hu',8,namey='Hu')

#############################
#### Fig S16
#############################
getDataSets('u3m','data_4a_si',namey='AS',constants={'input':60000.},Ftitle='crosstalk_AS_u3m_input.txt')
getDataSets_HHLoc('u3m','data_4a_si',8,namey='AS',constants={'input':60000.},Ftitle='crosstalk_AS_u3m_input.txt')
getDataSets('Hu','data_4b_si',namey='u3m',constants={'AS':0.2},Ftitle='crosstalk_AS_Hu_u3m.txt')
getDataSets_HHLoc('Hu','data_4b_si',8,namey='u3m',constants={'AS':0.2},Ftitle='crosstalk_AS_Hu_u3m.txt')
getDataSets('u3m','data_s4',namey='iHu',constants={'input':20000.})
getDataSets_HHLoc('u3m','data_s4',8,namey='iHu',constants={'input':20000.})
getDataSets_EMWO('u3m','data_s4',8,namey='iHu',constants={'input':20000.})
getDataSets('u3m','data_s4d',namey='AZ',constants={'input':20000.,'AS':0.95,'Au':1.1,'HS':1.1,'Hu':0.1,'UHV':310,'U3N':0.2},Ftitle='crosstalk_all.txt')
getDataSets_HHLoc('u3m','data_s4d',8,namey='AZ',constants={'input':20000.,'AS':0.95,'Au':1.1,'HS':1.1,'Hu':0.1,'UHV':310,'U3N':0.2},Ftitle='crosstalk_all.txt')
#############################
#### Fig S17
#############################
getData_NOHHICS()

#############################
#### Fig S18
#############################
getUH_correct("noWO_crosstalk_uh.txt","data_nowo_uh",None)

#############################
#### Fig S19
#############################
getDataSets('AZ','data_noem_AZ',start='noEM_')
getDataSets_ics('AZ','data_noem_AZ','AZ',start='noEM_',title="noEM_crosstalk_AZ.txt")
getDataSets('Au','data_noem_Au',start='noEM_')
getDataSets_ics('Au','data_noem_Au','Au',start='noEM_',title="noEM_crosstalk_Au.txt")
getDataSets('AS','data_noem_AS',start='noEM_')
getDataSets_ics('AS','data_noem_AS','AS',start='noEM_',title="noEM_crosstalk_AS.txt")
getDataSets('HS','data_noem_HS',start='noEM_')
getDataSets_ics('HS','data_noem_HS','HS',start='noEM_',title="noEM_crosstalk_HS.txt")

#############################
#### Fig S20
#############################
getDataSets('AS','data_noem_AS_HS',namey='HS',start='noEM_')
getDataSets_HHLoc('AS','data_noem_AS_HS',8,namey='HS',start='noEM_')
getDataSets_EMWO('AS','data_noem_AS_HS',8,namey='HS',start='noEM_')

getDataSets('Au','data_noem_Au_HS',namey='HS',start='noEM_')
getDataSets_HHLoc('Au','data_noem_Au_HS',8,namey='HS',start='noEM_')
getDataSets_EMWO('Au','data_noem_Au_HS',8,namey='HS',start='noEM_')

getDataSets('AS','data_noem_AS_Hu',namey='Hu',start='noEM_')
getDataSets_HHLoc('AS','data_noem_AS_Hu',8,namey='Hu',start='noEM_')
getDataSets_EMWO('AS','data_noem_AS_Hu',8,namey='Hu',start='noEM_')

getDataSets('Au','data_noem_Au_Hu',namey='Hu',start='noEM_')
getDataSets_HHLoc('Au','data_noem_Au_Hu',8,namey='Hu',start='noEM_')
getDataSets_EMWO('Au','data_noem_Au_Hu',8,namey='Hu',start='noEM_')

#############################
#### Fig S21
#############################
getDataSets('Hu','data_s6cnoHH',namey='u3m',constants={'input':10000.},start='noHH_')
getDataSets_HHLoc('Hu','data_s6cnoHH',8,namey='u3m',constants={'input':10000.},start='noHH_')

#############################
#### Fig S22
#############################

getDataSets('AS','data_PSF_AS',start='PSF_')
getDataSets_ics('AS','data_PSF_AS','AS',start='PSF_',title="PSF_crosstalk_AS.txt")
getDataSets('Au','data_PSF_Au',start='PSF_')
getDataSets_ics('Au','data_PSF_Au','Au',start='PSF_',title="PSF_crosstalk_Au.txt")

getDataSets('HS','data_PSF_HS',start='PSF_')
getDataSets_ics('HS','data_PSF_HS','HS',start='PSF_',title="PSF_crosstalk_HS.txt")
getDataSets('Hu','data_PSF_Hu',start='PSF_')
getDataSets_ics('Hu','data_PSF_Hu','Hu',start='PSF_',title="PSF_crosstalk_Hu.txt")

getDataSets('u3m','data_PSF_u3m',start='PSF_')
getDataSets_ics('u3m','data_PSF_u3m','u3m',start='PSF_',title="PSF_crosstalk_u3m.txt")
getDataSets('u3n','data_PSF_u3n',start='PSF_')
getDataSets_ics('u3n','data_PSF_u3n','u3n',start='PSF_',title="PSF_crosstalk_u3n.txt")

#############################
#### Fig S23
#############################
getDataSets('Hu','data_7cPSF',namey='u3m',constants={'input':10000.},start='PSF_')
getDataSets_HHLoc('Hu','data_7cPSF',8,namey='u3m',constants={'input':10000.},start='PSF_')
getDataSets_EMWO('Hu','data_7cPSF',8,namey='u3m',constants={'input':10000.},start='PSF_')
getDataSets('u3m','data_4cComp',namey='iHu',constants={'input':10000.})
getDataSets_HHLoc('u3m','data_4cComp',8,namey='iHu',constants={'input':10000.})
#############################
#############################
exit()



getDataSets('u3n','data_u3nHuInput',namey='iHu',constants={'input':60000.})
getDataSets_HHLoc('u3n','data_u3nHuInput',8,namey='iHu',constants={'input':60000.})
exit()

## comparison data

## full for 4A

getDataSets('u3n','data_s2',namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_s2',7,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_EMWO('u3n','data_s2',7,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')

## Full for 4B
getDataSets('AS','data_s3',namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},Ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')
getDataSets_HHLoc('AS','data_s3',4,namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},Ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')
getDataSets_EMWO('AS','data_s3',4,namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},Ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')

## Full for 4C


### Full for 5C
getDataSets('Hu','data_s5',namey='u3m',constants={'input':10000.},start='noHH_')
getDataSets_HHLoc('Hu','data_s5',8,namey='u3m',constants={'input':10000.},start='noHH_')
getDataSets_EMWO('Hu','data_s5',8,namey='u3m',constants={'input':10000.},start='noHH_')

getDataSets('Hu','data_s6',namey='u3m',constants={'AZ':0.9},start='noHH_',Ftitle='noHH_crosstalk_AZ_Hu_u3m.txt')
getDataSets_HHLoc('Hu','data_s6',8,namey='u3m',constants={'AZ':0.9},start='noHH_',Ftitle='noHH_crosstalk_AZ_Hu_u3m.txt')
getDataSets_EMWO('Hu','data_s6',8,namey='u3m',constants={'AZ':0.9},start='noHH_',Ftitle='noHH_crosstalk_AZ_Hu_u3m.txt')
getDataSets_EMWO('AS','data_s6b',8,namey='HS')


## noWO 
## noEM 
#getDataSets('uh','data_nowo_uh',start='noWO_')
#getDataSets_ics('uh','data_nowo_uh','uh',start='noWO_',title="noWO_crosstalk_uh.txt")
###################################################
### REDO WITH getDataSets_uhSinglesnoWO ######################################################
###################################################
getDataSets('u3m','data_nowo_u3m',start='noWO_')
getDataSets_ics('u3m','data_nowo_u3m','u3m',start='noWO_',title="noWO_crosstalk_u3m.txt")

getDataSets('u3n','data_nowo_u3n',start='noWO_')
getDataSets_ics('u3n','data_nowo_u3n','u3n',start='noWO_',title="noWO_crosstalk_u3n.txt")


getDataSets('AZ','data_s8',namey='Hu',start='noEM_')
getDataSets_HHLoc('AZ','data_s8',8,namey='Hu',start='noEM_')
getDataSets_EMWO('AZ','data_s8',8,namey='Hu',start='noEM_')

####
#getDataSets_uhSingles('data_s2b_uh')




getData_1c()
getData_NOHHICS()
getData_PSF_ICS()
##getDataSets('u3n','data_PSF_u3n_u3m',namey='u3m',start='PSF_',Ftitle='PSF_crosstalk_u3m_u3n.txt')
getDataSets_HHLoc('u3n','data_PSF_u3n_u3m',8,namey='u3m',start='PSF_',Ftitle='PSF_crosstalk_u3m_u3n.txt')
getDataSets_EMWO('u3n','data_PSF_u3n_u3m',8,namey='u3m',start='PSF_',Ftitle='PSF_crosstalk_u3m_u3n.txt')

getDataSets('AS','data_PSF_AS_Hu',namey='Hu',start='PSF_')
getDataSets_HHLoc('AS','data_PSF_AS_Hu',8,namey='Hu',start='PSF_')
getDataSets_EMWO('AS','data_PSF_AS_Hu',8,namey='Hu',start='PSF_')

getDataSets('Au','data_PSF_Au_HS',namey='HS',start='PSF_')
getDataSets_HHLoc('Au','data_PSF_Au_HS',8,namey='HS',start='PSF_')
getDataSets_EMWO('Au','data_PSF_Au_HS',8,namey='HS',start='PSF_')

### u3m, u3n, uh
getDataSets_EMWO('u3n','data_3_MR1_red',6,namey='u3m',constants={'UHV':300},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_EMWO('u3n','data_3_MR2_red',6,namey='u3m',constants={'UHV':301},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_EMWO('u3n','data_3_MR3_red',6,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_EMWO('u3n','data_3_MR4_red',6,namey='u3m',constants={'UHV':311},Ftitle='crosstalk_uh_u3n_u3m.txt')

getDataSets_HHLoc('u3n','data_3_MR1_red',6,namey='u3m',constants={'UHV':300},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR2_red',6,namey='u3m',constants={'UHV':301},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR3_red',6,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR4_red',6,namey='u3m',constants={'UHV':311},Ftitle='crosstalk_uh_u3n_u3m.txt')

getDataSets_reduced('u3n','data_3_MR1_red','MR',namey='u3m',constants={'UHV':300},ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_reduced('u3n','data_3_MR2_red','MR',namey='u3m',constants={'UHV':301},ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_reduced('u3n','data_3_MR3_red','MR',namey='u3m',constants={'UHV':310},ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_reduced('u3n','data_3_MR4_red','MR',namey='u3m',constants={'UHV':311},ftitle='crosstalk_uh_u3n_u3m.txt')

getDataSets('u3n','data_3_MR1',namey='u3m',constants={'UHV':300},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR1',6,namey='u3m',constants={'UHV':300},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets('u3n','data_3_MR2',namey='u3m',constants={'UHV':301},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR2',6,namey='u3m',constants={'UHV':301},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets('u3n','data_3_MR3',namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR3',6,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets('u3n','data_3_MR4',namey='u3m',constants={'UHV':311},Ftitle='crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR4',6,namey='u3m',constants={'UHV':311},Ftitle='crosstalk_uh_u3n_u3m.txt')

#getDataSets_uhSingles('data_s_uh')


getDataSets('u3n','data_u3n_u3m',namey='u3m')
getDataSets_HHLoc('u3n','data_u3n_u3m',8,namey='u3m')
getDataSets_EMWO('u3n','data_u3n_u3m',8,namey='u3m')








############





