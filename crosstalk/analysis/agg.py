import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


dfO = pd.read_csv("onlyHH_files.txt")
dfE = pd.read_csv("eq_files.txt")
dfN = pd.read_csv("none_files.txt")
dfU = pd.read_csv("upreg_files.txt")
dfD = pd.read_csv("downreg_files.txt")


print "Only"
for el in np.unique(dfO['filen'].values):
	if "data/crosstalk" in el:
		print el

print "\n\n"
print "Upreg"
for el in np.unique(dfU['filen'].values):
	if "data/crosstalk" in el:
		inds = np.argwhere(dfU['filen'].values==el)[:,0]
		tmp = dfU['amount'].values[inds]
		print np.sum(tmp>90),'out of',len(tmp),'\t',el

#print "Equal", np.unique(dfE['filen'].values)
#print "None", np.unique(dfN['filen'].values)
#print "Downreg", np.unique(dfD['filen'].values)

