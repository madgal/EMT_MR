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

for i in range(1,8):
    dirN='../coupledWReg_Ccode/crosstalk_uh/EMT_MR'+str(i)+'/'
    fsave='data/crosstalk_uh'+str(i)+'.txt'
    output_results(dirN,fsave)

for dirN in os.listdir('../coupledWReg_Ccode/'):
	if 'crosstalk' in dirN:
	    try:
		fsave = 'data/'+dirN+'.txt'
		dirname='../coupledWReg_Ccode/'+dirN+'/'
		output_results(dirname,fsave)
		print 'Succeede',dirN
	    except:
		print "failure with", dirN
