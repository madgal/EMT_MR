import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_all = pd.read_csv("../../crosstalk_comparison/EMT_MR_comp_0_1000_res.txt")

df_EM = pd.read_csv("state_EM_1000_res.txt")
#df_El = pd.read_csv("state_E_lou_0_5_1000_res.txt")
#df_Ml = pd.read_csv("state_M_loz_0_5_1000_res.txt")
df_E  = pd.read_csv("state_E_Au_15_0_AS_0_0_AZ_0_0_1000_res.txt")
df_E2 = pd.read_csv("state_E_1000_res.txt")#pd.read_csv("state_E_Au_3_1_AS_0_5_AZ_0_5_1000_res.txt")
##df_E2 = pd.read_csv("state_E_Au_1_1_AS_0_9_AZ_0_9_1000_res.txt")
#df_E3 = pd.read_csv("state_E_lou_0_5_Au_1_1_AS_0_9_AZ_0_9_1000_res.txt")
#df_E4 = pd.read_csv("state_E_lou_0_8_Au_1_1_AS_0_9_AZ_0_9_1000_res.txt")
#df_E5 = pd.read_csv("state_E_lou_0_8_Au_10_1_AS_0_0_AZ_0_0_1000_res.txt")
#df_E6 = pd.read_csv("state_E_lou_0_5_Au_10_1_AS_0_0_AZ_0_0_1000_res.txt")
df_M  = pd.read_csv("state_M_Hu_0_0_HS_20_0_1000_res.txt")
df_M2 = pd.read_csv("state_M_1000_res.txt")#pd.read_csv("state_M_Hu_0_5_HS_5_1_1000_res.txt")
##df_M2 = pd.read_csv("state_M_Hu_0_9_HS_1_1_1000_res.txt")
#df_M3 = pd.read_csv("state_M_loz_0_5_Hu_0_9_HS_1_1_1000_res.txt")
#df_M4 = pd.read_csv("state_M_loz_0_8_Hu_0_9_HS_1_1_1000_res.txt")
#df_M5 = pd.read_csv("state_M_loz_0_8_Hu_0_0_HS_10_1_1000_res.txt")
#df_M6 = pd.read_csv("state_M_loz_0_5_Hu_0_0_HS_10_1_1000_res.txt")



for key in ['u','u3','ms','S','mz','Z','G','O','mg','mo']:#df_all.columns:
	if key in df_all.columns:
		print key, np.unique(df_E[key]),np.unique(df_EM[key]),np.unique(df_M[key]),np.unique(df_all[key])
		print key, np.unique(df_E2[key]),np.unique(df_M2[key])
	else:
		print key,np.unique(df_E[key]),np.unique(df_EM[key]),np.unique(df_M[key])
		print key, np.unique(df_E2[key]),np.unique(df_M2[key])

print ''
for key in ['u','u3','ms','S','mz','Z','G','O','mg','mo']:#df_all.columns:
	if key in df_all.columns:
		print key, np.mean(df_E[key]),np.mean(df_EM[key]),np.mean(df_M[key]),np.mean(df_all[key])
		print key, np.mean(df_E2[key]),np.mean(df_EM[key]),np.mean(df_M2[key]),np.mean(df_all[key])
	else:
		print key,np.mean(df_E[key]),np.mean(df_EM[key]),np.mean(df_M[key])
		print key, np.mean(df_E2[key]),np.mean(df_EM[key]),np.mean(df_M2[key])


exit()
for key in ['u','u3','ms','S','mz','Z','G','O','mg','mo']:#df_all.columns:
	maxv = np.max([np.max(df_E[key]),np.max(df_EM[key]),np.max(df_M[key])])
	minv = np.min([np.min(df_E[key]),np.min(df_EM[key]),np.min(df_M[key])])
	#bins = np.linspace(minv,maxv,100)
	plt.hist(df_E[key],color='r',label='E')
	plt.hist(df_EM[key],color='b',label='EM')
	plt.hist(df_M[key],color='g',label='M')
	plt.xlim(minv,maxv)
	plt.legend()
	plt.show()



