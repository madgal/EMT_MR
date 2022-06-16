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
        mean_list = np.mean(full_setp,axis=0)
        std_list = np.std(full_setp,axis=0)
        uncoupledZ = (df.values-mean_list)/std_list
        uncoupled_fp = np.unique(np.round(uncoupledZ,6),axis=0)
        stateLabels = returnStateLabels(uncoupled_fp)

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
def getData_1c():
    fileo = open("data_1c.txt","w")
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
        keys =constants.keys()
        key = keys[0]
        inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        for i in range(1,len(keys)):
            tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
            inds = np.intersect1d(inds,tmpI)
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
############
def getDataSets_reduced(namex,title,typeR,namey=None,constants=None,start=None):
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
    df = pd.read_csv("~/Research/EMT_MR/crosstalk/analysis/data/"+start+"crosstalk_"+name+".txt",header=0)
    nics = float(df['nics'][0])

    if constants!=None:
       keys =constants.keys()
       key = keys[0]
       inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
       for i in range(1,len(keys)):
           tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
           inds = np.intersect1d(inds,tmpI)
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



    fileo = open(title+".txt",'w')
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
        inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
        for i in range(1,len(keys)):
            tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
            inds = np.intersect1d(inds,tmpI)
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
          inds = np.argwhere(df[key.upper()].values==constants[key])[:,0]
          for i in range(1,len(keys)):
               tmpI = np.argwhere(df[keys[i].upper()].values==constants[keys[i]])[:,0]
               inds = np.intersect1d(inds,tmpI)
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

#######3
#######3
#######3
def getData_2d():
    fileo = open("data_2d.txt","w")
    fileo.write("xlabel")
    keyList=['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO','W']
    for key in keyList:
        fileo.write(",%s" %(key))
    fileo.write("\n")

    count=0
    for dirn in os.listdir("../../crosstalk/coupledWReg_Ccode/singles/."):
       if "crosstalk_uh"!=dirn and "crosstalk_input"!=dirn and "crosstalk_u3m"!=dirn:
                     df =pd.read_csv("../../crosstalk/analysis/data/"+dirn+".txt")
                     tmpL=dirn.split("_")[1:][0].upper()
                     if tmpL=='IHU':
                          tmpL='HU'
                     mystr ="$\lambda_{"+tmpL+"}=$"
                     nics =df['nics'].values[0]*1.

                     set_reg= df[tmpL].values
                     minr,maxr=np.min(set_reg),np.max(set_reg)
                     if minr<1 and maxr<=1:
                          ## inhibitory
                          regV = minr
                     elif maxr>1 and minr>=1:
                          ## activation
                          regV = maxr
                     else:
                          print("issue with", dirn)
                          regV=1

                     ind = np.argwhere(regV==set_reg)[:,0]
                     if len(ind)!=1:
                          print(len(ind)," issue for ",dirn)
                     else:
                          mystr=mystr+str(regV)
                          fileo.write("%s" %(mystr))
                          for key in keyList:
                               fileo.write(",%s" %(float(df[key].values[ind])/float(nics)*100.))
                          fileo.write("\n")
                          count+=1

       elif "crosstalk_uh"==dirn or "crosstalk_u3m"==dirn:
           if dirn=="crosstalk_uh":
                 df =pd.read_csv("../../crosstalk/analysis/data/"+dirn+"10.txt")
           else:
                 df =pd.read_csv("../../crosstalk/analysis/data/"+dirn+".txt")
           tmpL=dirn.split("_")[1:][0].upper()
           if tmpL=='UH':
                 regV= 0.817
                 compV=2.817
                 roundV=3
           if tmpL=='U3M':
                 regV=0.1
                 compV=0.1
                 roundV=6
           mystr ="$\lambda_{"+tmpL+"}=$"
           nics =df['nics'].values[0]*1.

           set_reg= np.round(df[tmpL].values,roundV)
           ind = np.argwhere(compV==set_reg)[:,0]
           if len(ind)!=1:
                 print(len(ind)," issue for ",dirn)
           else:
                 mystr=mystr+str(regV)
                 fileo.write("%s" %(mystr))
                 for key in keyList:
                        fileo.write(",%s" %(float(df[key].values[ind])/float(nics)*100.))
                 fileo.write("\n")
                 count+=1
    fileo.close()
##########3
##########3
##########3
def getData_comp():
    
	df =pd.read_csv("~/Research/EMT_MR/crosstalk/coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()

	fileo = open("data_comp.txt","w")
	for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO']:#,'W']:
		fileo.write("%s," %key)
	fileo.write("W\n")

	mapRes = getStates(df)
	for key in ['E/O','E/W','E/WO','EM/O','EM/W','EM/WO','M/O','M/W','M/WO','E','EM','M','O','WO']:#,'W']:
		fileo.write("%s," %(mapRes[key]/10.))
	fileo.write("%s\n" %(mapRes['W']/10.))
	fileo.close()
############3
############3
############3
############3
getDataSets_reduced('AZ','data/d_no_EM_AZ_red','EMT',start='noEM_')
getDataSets_reduced('HS','data/d_no_EM_HS_red','EMT',start='noEM_')
getDataSets_reduced('Hu','data/d_no_EM_Hu_red','EMT',start='noEM_')
getDataSets_reduced('AZ','data/d_no_EM_AZ_Hu_red','EMT',start='noEM_',namey='Hu')
getDataSets_reduced('Hu','data/d_no_EM_Hu_input_red','EMT',start='noEM_',namey='input')


getDataSets_reduced('u3m','data/d_noWO_u3m_red','MR',start='noWO_')
getDataSets_reduced('u3n','data/d_noWO_u3n_red','MR',start='noWO_')
