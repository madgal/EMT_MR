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

#from aux_func_States import *
#from functions_for_plotting import *

'''
file0 = open("upreg_files.txt","w")
fileF = open("onlyHH_files.txt","w")
file1 = open("eq_files.txt","w")
file2 = open("downreg_files.txt","w")
file3 = open("none_files.txt","w")
'''


df = pd.read_csv("onlyHH_files.txt")
print np.unique(df['filen'])
data/crosstalk_u3m_iHu_input.txt
