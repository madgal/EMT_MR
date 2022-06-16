import numpy as np
import pandas as pd
import os
import math
from math import factorial
import matplotlib.pyplot as plt


def getLamdaSets(lvx,lvy,xk=[],yk=[],count=0,tmpx=[],tmpy=[],finx=[],finy=[]):

    print count,finx
    if count==0:
        print 'A',xk,yk,count,tmpx,tmpy,finx#,finy
        print 'B',lvx#,lvy
    else:
        print 'C'

    if xk==[] and yk==[]:
        xk = lvx.keys()
        yk = lvy.keys()
    if count==(len(xk)+len(yk)-1):
        ntx,nty = [],[]
        for i in range(len(lvy[yk[-1]])):
            ty = tmpy+[lvy[yk[-1]][i]]
            ntx+=[tmpx]
            nty+=[ty]
        finx+=ntx
        finy+=nty
        return finx,finy
    else:
        if count <len(xk):
            for i in range(len(lvx[xk[count]])):
                t2 = tmpx+[lvx[xk[count]][i]]
                finx,finy=getLamdaSets(lvx,lvy,xk,yk,count+1,t2,tmpy,finx,finy)
        elif count <len(xk)+len(yk):
            for i in range(len(lvy[yk[count%len(xk)]])):
                t2 = tmpy+[lvy[yk[count%len(xk)]][i]]
                finx,finy=getLamdaSets(lvx,lvy,xk,yk,count+1,tmpx,t2,finx,finy)
        else:
            return finx,finy
    return finx,finy
    

lx = {0:[0.1,0.2,0.3]}
ly={0:[1,2,3],1:[3,4,5]}
lx1={0:[3,4,5]}
ly2={0:[9,8,7]}

r1x,r1y=getLamdaSets(lx,ly)


r2x,r2y=getLamdaSets(lx1,ly2)


print "1",r1x
print "2",r2x
