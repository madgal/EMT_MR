import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib import cm
import os
import math

mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams.update({'font.size': 20})

#####################
#####################
#####################
#####################
#####################
#####################

def compICS(fileDir):

	icsOG = pd.read_csv("../coupledWReg_Ccode/crosstalk_comparison/EMT_MR_comp_0_1000_ics.txt")

	for file in  os.listdir(fileDir):
		if 'ics.txt' in file:
			icsN = pd.read_csv(fileDir+file)

			for key in icsN.columns:

				tmpa = np.sort(icsOG[key].values)
				tmpb = np.sort(icsN[key].values)
				print len(tmpa),len(tmpb),key


				plt.plot(tmpa,marker='o',markersize=10)
				plt.plot(tmpb,marker='*',markersize=5)

				for i in range(len(tmpb)):
					print tmpa[i],tmpb[i]

				plt.show()

			print file
			exit()

regs_uh = compICS("../coupledWReg_Ccode/crosstalk_uh/")
regs_hS_aS = compICS("../coupledWReg_Ccode/crosstalk_hS_aS/")
regs_mrNr = compICS("../coupledWReg_Ccode/crosstalk_mrNr/")
regs_uh_hu = compICS("../coupledWReg_Ccode/crosstalk_uh_hu/")
regs_uh_hu_mrNr = compICS("../coupledWReg_Ccode/crosstalk_uh_hu_mrNr/")
regs_uh_hu_au = compICS("../coupledWReg_Ccode/crosstalk_uh_hu_au/")
regs_hS_aS_AZ = compICS("../coupledWReg_Ccode/crosstalk_hS_aS_AZ/")
