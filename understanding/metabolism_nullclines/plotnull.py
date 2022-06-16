import numpy as np
import matplotlib.pyplot as plt
import json

with open("data_metabolism_nullcline_og.json") as dataF:
	res = json.load(dataF)

fp = res['fixed_points']
null = res['nullA_A']
plt.plot(res['nullA_A'],res['nullA_H'],label='A')
plt.plot(res['nullH_A'],res['nullH_H'],label='A')
plt.title("og")
plt.show()


with open("data_metabolism_nullcline_miRNA.json") as dataF:
	res = json.load(dataF)

fp = res['fixed_points']
null = res['nullA_A']
plt.plot(res['nullA_A'],res['nullA_H'],label='A')
plt.plot(res['nullH_A'],res['nullH_H'],label='A')
plt.title("miRNA")
plt.show()
