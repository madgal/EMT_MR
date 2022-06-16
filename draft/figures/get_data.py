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
def getData_1b():
    
        df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
        x1L = 0#'u'
        x2L = 6#'A'
        y1L=9#'h'
        y2L=1#'mz'

        color_list=['k','r','g','b','y','pink','orange','purple','indigo']
        marker=['p','o','s']
        colorD={'E':'k','M':'b','EM':'m'}
        markerD={'WO':'P','W':'o','O':'s'}

        fileo = open("data_1b_zval.txt","w")
        file2 = open("data_1b.txt","w")
        columns = df.columns
        fileo.write("%s,%s,%s,%s,label,color,marker\n"%(columns[x1L],columns[y1L],columns[x2L],columns[y2L]))
        file2.write("%s,%s,%s,%s,label,color,marker\n"%(columns[x1L],columns[y1L],columns[x2L],columns[y2L]))

        full_setp=df.values
        fs_inds= np.unique(full_setp,axis=0)
        uncoupled_fp,stateLabels = returnStateLabels(False,False)

        for i in range(len(uncoupled_fp)):
            [emt,mr] =stateLabels[i].split("/")
            col=colorD[emt]
            mark=markerD[mr]
            fileo.write("%s,%s,%s,%s,%s,%s,%s\n" %(uncoupled_fp[i,x1L],uncoupled_fp[i,y1L],uncoupled_fp[i,x2L],uncoupled_fp[i,y2L],stateLabels[i],col,mark))
            file2.write("%s,%s,%s,%s,%s,%s,%s\n" %(fs_inds[i,x1L],fs_inds[i,y1L],fs_inds[i,x2L],fs_inds[i,y2L],stateLabels[i],col,mark))

        fileo.close()
        file2.close()
############################
############################
############################
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


    fileo = open(Ftitle+"_ics.txt",'w')
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
	for i in range(7):
		for j in range(5):
			for k in range(6):
				if (j==0 and k==4) or (j==1 and (k==0 or k==4)) or (j==2 and (k<=2 or k==4)) or (j==3 and k<=5) or (j==4):
					xun_ok+=[i*100+j*10+k]
	tmp = np.round(df['UH'].values,3)
	for i in [9,8,7,6,5,4,3,2,1]:
		inds = np.argwhere(tmp>i)[:,0]
		if len(inds)>0:
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

	fileo=open(title+"_ics.txt",'w')
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
	fileo=open(title+".txt",'w')	
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
	x1v = np.max(x2)
	x2v = 1.
	colD=len(groups)+1
	fileo.write("%s,%s,-0.5,0.5,%s\n" %(x1v,x2v,colD))
	fileo.close()
 
	fileo = open(title+"_legend.txt",'w')
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
	for el in typeList:## add in the full set
		fileo.write("1,")
	fileo.write("%s\n" %colD)
	fileo.close()


############
#############
############
def getDataSets_uhSingles_red(title):
	df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/singles_crosstalk_uh10.txt")
	elV = {'O':1.,'WO':2.,'W':3.}
	y = df['O'].values*elV['O']+df['WO'].values*elV['WO']+df['W'].values*elV['W']
	xun_ok=[]
	for i in range(7):
		for j in range(5):
			for k in range(6):
				if (j==0 and k==4) or (j==1 and (k==0 or k==4)) or (j==2 and (k<=2 or k==4)) or (j==3 and k<=5) or (j==4):
					xun_ok+=[i*100+j*10+k]
	tmp = np.round(df['UH'].values,3)
	for i in [9,8,7,6,5,4,3,2,1]:
		inds = np.argwhere(tmp>i)[:,0]
		if len(inds)>0:
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
		xfics+=[el]
		yfics['W']  +=[np.sum(yf[inds]==elV['W'])*100./len(inds) ]
		yfics['O']  +=[np.sum(yf[inds]==elV['O']) *100./len(inds) ]
		yfics['WO'] +=[np.sum(yf[inds]==elV['WO']) *100./len(inds) ]

	fileo=open(title+"_ics.txt",'w')
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
				prx[key][i]=1
				
	fileo=open(title+".txt",'w')	
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
 
	fileo = open(title+"_legend.txt",'w')
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

    if not ftitle:
       ftitle=''
    if not name:
          name=''
    if ('UH' in ftitle.upper()) or  ('UH' in name.upper()):
        xun_ok=[]
        for i in range(7):
            for j in range(5):
                        for k in range(6):
                              if (j==0 and k==4) or (j==1 and (k==0 or k==4)) or (j==2 and (k<=2 or k==4)) or (j==3 and k<=5) or (j==4):
                                     xun_ok+=[i*100+j*10+k]

        tmp = np.round(df['UH'].values,3)
        for i in [9,8,7,6,5,4,3,2,1]:
             inds = np.argwhere(tmp>i)[:,0]
             if len(inds)>0:
                    tmp[inds]=tmp[inds]-i

        finalInds=[]
        for el in xun_ok:
            inds2 = np.argwhere(df['UHV'].values==el)[:,0]
            finalInds+=list(inds2)

        #print(len(inds),len(finalInds))
        inds = np.intersect1d(inds,np.array(finalInds))
        #print(len(inds))
        

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
            else                   :
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



    fileo = open(title+".txt",'w')
    fileo.write("X1,X2,Y1,Y2,color,HHV\n")
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
         fileo.write("%s,%s,%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i],colD,hybridV[i]))
    fileo.close()

    if typeR=='MR':
        fileo = open(title+"_legend.txt",'w')
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
        fileo = open(title+"_legend.txt",'w')
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

    fileo = open(title+".txt",'w')
    fileo.write("X1,X2,Y1,Y2,color,ecolor,hatch,colorHHReg,HHV\n")
    for i in range(len(x1)):
         if hybridV[i]==100:
              tmpC=3
         elif hybridV[i]==compHybrid:
              tmpC=2
         elif hybridV[i]>compHybrid:
              tmpC=1
         else:#hybridV[i]<compHybrid:
             tmpC=0
         fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(x1[i],x2[i],y1[i],y2[i],colD[color[i]],edgcolors[i],hatch[i],tmpC,hybridV[i]))
    fileo.close()

    fileo = open(title+"_legend.txt",'w')
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
    fileoE = open(title+"_HHexists.txt",'w')
    fileoO = open(title+"_HHOnlhy.txt",'w')
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
        
            fileoH = open(title+"_hull.txt",'w')
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
        
            fileoH = open(title+"_hullO.txt",'w')
            fileoH.write("x,y\n")
            for i in range(len(exterior)):
                    fileoH.write("%s,%s\n" %(exterior[i][0],exterior[i][1]))
            fileoH.close()
        except:
           print("Error in getting hull",len(xt),ashape)
############
############
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

    fileoE = open(title+"_EMWO.txt",'w')
    fileoE.write("X,Y,EMWO,Change\n")
    for i in range(len(xval)):
         tmp = (hybridV[i]/compHybrid)
         fileoE.write("%s,%s,%s,%s\n" %(xval[i],yval[i],hybridV[i],tmp))
    fileoE.close()
    
############
############
############
############
############
def getDataSets_EMWO_diff(namex,title,namey,Ftitle,CompTitle,constants=None):
    ### get current data
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

    xval,yval=[],[]
    hybridV = []
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



    ### get comparison data
    df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+CompTitle,header=0)
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

    xvalC,yvalC,zvalC=[],[],[]
    ## acutal values are corred
    for i in range(len(df.values[inds])):
        if 'UH' in namex.upper():
                xvalC+=[df['UH'][inds[i]]]
        elif 'HU' in namex.upper():
                xvalC+=[df['HU'][inds[i]]]
        else:
                xvalC+=[df[namex.upper()][inds[i]]]
        if namey!=None:
            if 'UH' in namey.upper():
               yvalC+=[df['UH'][inds[i]]]
            elif 'HU' in namey.upper():
               yvalC+=[df['HU'][inds[i]]]
            else:
               yvalC+=[df[namey.upper()][inds[i]]]
        else:
               yvalC+=[0]
        zvalC+=[df['EM/WO'][inds[i]]*1./nics*100.]


    ## calculate difference
    xf,yf,zf=[],[],[]
    for i in range(len(xval)):
        indx = np.argwhere(xval[i]==xvalC)[:,0]
        indy = np.argwhere(yval[i]==yvalC)[:,0]
        inds = np.intersect1d(indx,indy)
        if len(inds)>0:
            diff = hybridV[i]-zvalC[inds[0]]

            xf+=[xval[i]]
            yf+=[yval[i]]
            zf+=[diff]
	

    ####
    fileoE = open(title+"_EMWO_diff.txt",'w')
    fileoE.write("X,Y,Diff\n")
    for i in range(len(xf)):
         fileoE.write("%s,%s,%s\n" %(xf[i],yf[i],zf[i]))
    fileoE.close()
   
############
############
def getData_Fig4_uh_u3n_u3mNOX():
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
    tmpy = df['U3M'].values
    for i in [9,8,7,6,5,4,3,2,1]:
         inds = np.argwhere(tmp>i)[:,0]
         if len(inds)>0:
               tmp[inds]=tmp[inds]-i

    #### Only keep if the value of U3M=0 and the UH parameters are correct
    ### x = uh, y = u3n, z = coupled state
    xf,yf,zf=[],[],[]
    inds2 = np.argwhere(df['U3N'].values==0)[:,0]
    inds3 = np.argwhere(df['U3M'].values>0)[:,0]
    inds2 = np.intersect1d(inds2,inds3)
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
    yf=np.array(yf)
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

    title = "data_3_MR3_single_NOX0"
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
############
############
##
def getData_2():
	getDataSets('u3n','data_2a')
	getDataSets_ics('u3n','data_2a','u3n',title='crosstalk_u3n.txt')
############
def getData_4():
	getDataSets('u3n','data_4a',namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
	getDataSets_HHLoc('u3n','data_4a',7,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')

	getDataSets('AS','data_4b',namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},Ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')
	getDataSets_HHLoc('AS','data_4b',4,namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},Ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')

	getDataSets('u3m','data_4c',namey='iHu',constants={'input':20000.})
	getDataSets_HHLoc('u3m','data_4c',8,namey='iHu',constants={'input':20000.})

	getDataSets('u3m','data_4d',namey='AZ',constants={'input':20000.,'AS':0.95,'Au':1.1,'HS':1.1,'Hu':0.1,'UHV':310,'U3N':0.2},Ftitle='crosstalk_all.txt')
	getDataSets_HHLoc('u3m','data_4d',8,namey='AZ',constants={'input':20000.,'AS':0.95,'Au':1.1,'HS':1.1,'Hu':0.1,'UHV':310,'U3N':0.2},Ftitle='crosstalk_all.txt')

def getData_5():
	getDataSets_reduced('u3m','data_5a','MR',start='noWO_')
	getDataSets_ics('u3m','data_5a','u3m',start='noWO_',title="noWO_crosstalk_u3m.txt")
	getDataSets('u3m','data_5a_full',start='noWO_',Ftitle="noWO_crosstalk_u3m.txt")
	getDataSets_HHLoc('u3m','data_5a_full',8,start='noWO_',Ftitle="noWO_crosstalk_u3m.txt")

	getDataSets_reduced('Hu','data_5b',"EMT",start='noEM_')
	getDataSets_ics('Hu','data_5b','Hu',start='noEM_')
	getDataSets('Hu','data_5b_full',start='noEM_',Ftitle="noEM_crosstalk_Hu.txt")
	getDataSets_HHLoc('Hu','data_5b_full',8,start='noEM_',Ftitle="noEM_crosstalk_Hu.txt")

	getDataSets('Hu','data_5c',namey='u3m',constants={'input':10000.},start='noHH_')
	getDataSets_HHLoc('Hu','data_5c',8,namey='u3m',constants={'input':10000.},start='noHH_')


def getData_6():
	getDataSets('Hu','data_6c1',namey='u3m',constants={'input':10000.},start='PSF_')
	getDataSets_HHLoc('Hu','data_6c1',8,namey='u3m',constants={'input':10000.},start='PSF_')
	getDataSets_EMWO('Hu','data_6c1',8,namey='u3m',constants={'input':10000.},start='PSF_')
	getDataSets('AS','data_6c2',namey='Hu',start='PSF_')
	getDataSets_HHLoc('AS','data_6c2',8,namey='Hu',start='PSF_')
	getDataSets_EMWO('AS','data_6c2',8,namey='Hu',start='PSF_')
	getDataSets('Au','data_6c3',namey='HS',start='PSF_')
	getDataSets_HHLoc('Au','data_6c3',8,namey='HS',start='PSF_')
	getDataSets_EMWO('Au','data_6c3',8,namey='HS',start='PSF_')

def getData_3():
	getDataSets_reduced('Au','data_3_Au','EMT')
	getDataSets_reduced('HS','data_3_HS','EMT')
	getDataSets_reduced('u3m','data_3_u3m','MR')
	getDataSets_uhSingles('data_3_uh')
	##getDataSets_reduced('uh10','data_3_uh','MR',constants={'UHV':tmp},ftitle='singles_crosstalk_uh10.txt')

	getDataSets_ics('Au','data_3_Au','Au')
	getDataSets_ics('HS','data_3_HS','HS')
	getDataSets_ics('u3m','data_3_u3m','u3m')
	###getDataSets_ics('uh10','data_3_uh','UH',constants={'UHV':tmp},title='singles_crosstalk_uh10.txt')


	getDataSets_reduced('u3n','data_3_MR1_red','MR',namey='u3m',constants={'UHV':310},ftitle='crosstalk_uh_u3n_u3m.txt')
	getDataSets_HHLoc('u3n','data_3_MR1_red',6,namey='u3m',constants={'UHV':310},Ftitle='crosstalk_uh_u3n_u3m.txt')
	getDataSets_reduced('u3n','data_3_MR2_red','MR',namey='u3m',constants={'UHV':333},ftitle='crosstalk_uh_u3n_u3m.txt')
	getDataSets_HHLoc('u3n','data_3_MR2_red',6,namey='u3m',constants={'UHV':333},Ftitle='crosstalk_uh_u3n_u3m.txt')
	getDataSets_reduced('u3n','data_3_MR3_red','MR',namey='uh',constants={'U3M':0.},ftitle='IND_crosstalk_uh_u3n_u3m.txt')
	getDataSets_HHLoc('u3n','data_3_MR3_red',6,namey='uh',constants={'U3M':0.},Ftitle='IND_crosstalk_uh_u3n_u3m.txt')

	getDataSets_reduced('AS','data_3_EMT_red','EMT',namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')
	getDataSets_HHLoc('AS','data_3_EMT_red',8,namey='Au',constants={'input':50000.,'AZ':0.4,'HS':1.1,'Hu':0.1},Ftitle='crosstalk_AS_AZ_Au_HS_Hu_input.txt')


########################
########################
########################
########################
########################
########################
getData_Fig4_uh_u3n_u3mNOX()
getDataSets_uhSingles('data_uh_coupled')
exit()

getData_1b()
getData_2()
getData_3()
getData_4()
getData_5()
getData_6()

getDataSets_EMWO_diff('Hu','data_6c1','u3m',"PSF_crosstalk_Hu_u3m_input.txt","crosstalk_u3m_iHu_input.txt",constants={'input':20000})
getDataSets_EMWO_diff('Hu','data_6c1_1','u3m',"PSF_crosstalk_Hu_u3m_input.txt","crosstalk_u3m_iHu_input.txt",constants={'input':10000})
getDataSets_EMWO('Hu','data_testTris',8,namey='u3m',constants={'input':10000.},Ftitle='crosstalk_u3m_iHu_input.txt')
getDataSets_EMWO('Hu','data_testTris2',8,namey='u3m',constants={'input':20000.},Ftitle='crosstalk_u3m_iHu_input.txt')
getDataSets_EMWO('Hu','data_testPSF',8,namey='u3m',constants={'input':20000.},Ftitle='PSF_crosstalk_Hu_u3m_input.txt')



####

getDataSets('u3m','data_AS_u3m_input',namey='AS',constants={'input':60000.},Ftitle='crosstalk_AS_u3m_input.txt')
getDataSets_HHLoc('u3m','data_AS_u3m_input',8,namey='AS',constants={'input':60000.},Ftitle='crosstalk_AS_u3m_input.txt')

getDataSets('Hu','data_AS_Hu_u3m',namey='u3m',constants={'AS':0.2},Ftitle='crosstalk_AS_Hu_u3m.txt')
getDataSets_HHLoc('Hu','data_AS_Hu_u3m',8,namey='u3m',constants={'AS':0.2},Ftitle='crosstalk_AS_Hu_u3m.txt')

getDataSets('u3n','data_3_MR3_full',namey='uh',constants={'U3M':0.},Ftitle='IND_crosstalk_uh_u3n_u3m.txt')
getDataSets_HHLoc('u3n','data_3_MR3_full',6,namey='uh',constants={'U3M':0.},Ftitle='IND_crosstalk_uh_u3n_u3m.txt')

