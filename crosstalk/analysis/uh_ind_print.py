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
 
def main():
      for i in [10]:#1,2,3,4,5,6,7,9,10]:
         dirN='../coupledWReg_Ccode/singles/crosstalk_uh/EMT_MR'+str(i)+'/'
         fsave='data/singles_crosstalk_uh'+str(i)+'.txt'
         if True:#try:
         	output_results_singles(dirN,fsave)
         else:#xcept:
		print "failure with", dirN



#######################
def output_results_singles(dirN,fsave=''):

    if not fsave:
	dirn_red = dirN.split("/")[-1]
	fsave = "data/"+dirn_red+".txt"
    nics= dirN.split("_")[-2]
    fileo=open(fsave,"w")

    fileo.write("UH,UHV,E,EM,M,W,WO,O,M/W,M/WO,M/O,EM/W,EM/WO,EM/O,E/W,E/WO,E/O,nics\n")

    for filen in os.listdir(dirN):
        if 'res.txt' in filen:
            tmp= filen.split("_")
            df =pd.read_csv(dirN+filen).dropna()
	    for i in range(len(df)):
		    df_tmp = df.iloc[[i]]
		    cxs = {'UH':1.,'UHV':0}
		    cxs,lamda_counted=add_lamdas(cxs,tmp,df_tmp)
		    if not lamda_counted:
			tmp = [dirN.replace("/",'').split("_")[-1],tmp[3],tmp[4],tmp[-2],'0']
			cxs,lamda_counted = add_lamdas(cxs,tmp,df_tmp)

		    if True:#try:
                         if 'G' in df_tmp.columns:
                         	res=getStates(df_tmp.drop(columns=['G','O','mg','mo','t']))
                         else:
                         	res=getStates(df_tmp)
		         fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(cxs['UH'],cxs['UHV'],res['E'],res['EM'],res['M'],res['W'],res['WO'],res['O'],res['M/W'],res['M/WO'],res['M/O'],res['EM/W'],res['EM/WO'],res['EM/O'],res['E/W'],res['E/WO'],res['E/O']))
		    else:#xcept:
			print "Issue with ", dirN,filen

    fileo.close()


main()
