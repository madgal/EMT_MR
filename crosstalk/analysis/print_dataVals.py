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
 
for i in [10]:#1,2,3,4,5,6,7,9,10]:
    dirN='../coupledWReg_Ccode/singles/crosstalk_uh/EMT_MR'+str(i)+'/'
    fsave='data/vals_crosstalk_uh'+str(i)+'.txt'
    try:
    	output_dataVals(dirN,fsave)
	print 'Succeeded',dirN
    except:
		print "failure with", dirN

for dirN in os.listdir('../coupledWReg_Ccode/singles'):
	if "crosstalk_" in dirN and 'crosstalk_uh' not in dirN:
	    try:
		fsave = 'data/vals_'+dirN+'.txt'
		dirname='../coupledWReg_Ccode/singles/'+dirN+'/'
		output_dataVals(dirname,fsave)
		print 'Succeeded',dirN
	    except:
		print "failure with", dirN
