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
import time

from aux_func_States import *
from functions_for_plotting import *

file0 = open("upreg_files.txt","w")
fileF = open("onlyHH_files.txt","w")
file1 = open("eq_files.txt","w")
file2 = open("downreg_files.txt","w")
file3 = open("none_files.txt","w")
file0.write("filen,amount,percentDiff,AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV\n")
file1.write("filen,amount,percentDiff,AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV\n")
file2.write("filen,amount,percentDiff,AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV\n")
fileF.write("filen,amount,percentDiff,AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV\n")
file3.write("filen,amount,percentDiff,AS,AZ,AU,HS,HU,INPUT,U3M,U3N,UH,UHV\n")
file0.close()
file1.close()
file2.close()
fileF.close()
file3.close()

uncoupled =pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_res.txt").dropna()
comp=getStates(uncoupled)['EM/WO']

compI={}
for filen in os.listdir("../coupledWReg_Ccode/singles/crosstalk_input/"):
	if "res.txt" in filen:
		key = int(str(filen).split("_")[3])
		uncoupledI =pd.read_csv("../coupledWReg_Ccode/singles/crosstalk_input/EMT_MR_input_"+str(key)+"_1000_res.txt").dropna()
		compI[key]=getStates(uncoupledI)['EM/WO']
		if key>1000:
			key = key/10000.
		compI[key]=getStates(uncoupledI)['EM/WO']

print compI
for filen in os.listdir("data"):
            if  ("crosstalk_comparison.txt"!=filen) and ("crosstalk_input.txt"!=filen) and ('crosstalk' in filen) and ('singles' not in filen) and ('vals' not in filen) and ('IND' not in filen) and ('bkp' not in filen):
		print filen
                filetitle='data/'+filen
                determineRegRange(filetitle,comp,compI)



####
### determine how they're coupled if there are 3 states, 3 EM , and 3WO types

file3 = open("groups_list.txt","w")
file3.write("filen,A,B,C\n")
file3.close()
file3 = open("groups.txt","w")
file3.write("filen,A,B,C\n")
file3.close()
for filen in os.listdir("data"):
            if  ("crosstalk_comparison.txt"!=filen) and ("crosstalk_input.txt"!=filen) and ('crosstalk' in filen) and ('singles' not in filen) and ('vals' not in filen) and ('IND' not in filen) and ('bkp' not in filen):
                filetitle='data/'+filen
		determineGroups(filetitle)
