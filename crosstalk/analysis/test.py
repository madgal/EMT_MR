import numpy as np
import pandas as pd
import os
import math
from math import factorial
import matplotlib.pyplot as plt


def P(li,ymi):
	## P=1 no silencing, P=0 full silencing
	km = 0.143 ## for mRNA of Hif-1 as that is the one we are varying
        return L2(li)/(1.+Ym2(ymi)/km)
def Ym2(ymi):
        total=0
        for i in range(len(ymi)):
                total+= ymi[i]*1.*factorial(len(ymi))/factorial(i)/factorial(len(ymi)-i)
        return total
def L2(li):
        total=0
        for i in range(len(li)):
                total+= li[i]*1.*factorial(len(li))/factorial(i)/factorial(len(li)-i)
        return total


count=0
for i in range(11):
	for j in range(11):
		plt.plot(count,P([1.,i,j],[0,0,0]),'*')
		count+=1

plt.show()
for i in range(11):
	for j in range(11):
		plt.plot(count,P([1.,1.,1.],[0,i,j]),'*')
		count+=1
plt.show()
