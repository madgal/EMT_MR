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
import shutil



for filen in os.listdir("data/"):
	if ('IND' in filen) and (not 'bkp' in filen):#(('uh' in filen) or ('uH' in filen) or ('all' in filen)) and (('bkp' not in filen) and ('singles' not in filen) and ('values' not in filen)):
		print filen
		df = pd.read_csv("data/"+filen)
		if len(df)>10000:
			shutil.copy("data/"+filen, "data/bkp_"+filen)
			### df.columns[:10]=['AS', 'AZ', 'AU', 'HS', 'HU', 'INPUT', 'U3M', 'U3N', 'UH','UHV']
			uniques= np.unique(df.values[:,:10],axis=0)
			uniques_F= np.unique(df.values[:,:],axis=0)
			if len(uniques)==len(uniques_F):## then we know there isn't inconsistencies in which P(u) values correspond to which sets
				df_un = pd.DataFrame(uniques_F, columns = df.columns,index=None)

				df_un.to_csv("data/"+filen,sep=',',index=None)
	
