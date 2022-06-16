import numpy as np
import matplotlib
import math
matplotlib.rcParams.update({'font.size': 30})
from matplotlib.lines import Line2D
from parameters_both import *

from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()

##################################3
##################################3
def M(u,u0,nX,i):
	return ((u/u0)**i)/((1+u/u0)**nX)
def combination(n,i):
	return math.factorial(n)/math.factorial(i)/math.factorial(n-i)
def L(u,u0,nX):
	li = [1.,0.6,0.3,0.1,0.05,0.05,0.05]
	total=0
	for i in range(nX+1):
		total +=li[i]*combination(nX,i)*M(u,u0,nX,i)
	return total
def Yu(u,u0,nX):
	yui = [0,0.005,0.05,0.5,0.5,0.5,0.5]
	total=0
	for i in range(nX+1):
		total +=i*yui[i]*combination(nX,i)*M(u,u0,nX,i)
	return total
def Ym(u,u0,nX):
	ymi = [0,0.04,0.2,1.,1.,1.,1.]
	total=0
	for i in range(nX+1):
		total +=ymi[i]*combination(nX,i)*M(u,u0,nX,i)
	return total
def CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm):
	return y*(gn+(A/A0rm)**narm)/(1.+(h/h0hrm)**nhrm+(A/A0rm)**narm)
def CompRn(g0,h,h0rn,nhrn,ghn,A,A0rn,gan,narn):
	return (g0+ghn*(h/h0rn)**nhrn+gan*(A/A0rn)**narn)/(1.+(h/h0rn)**nhrn+(A/A0rn)**narn)

def H(X,X0,l,nX):
	return (l+(1.-l)/(1.+(X/X0)**nX))
#################
#################
#################
def fZ(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return gz*mz*L(u,u0,nu)-kz*Z
def fS(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return gs*ms*L(u3,u30,nu3)-ks*S
def fu3(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return gu3*H(Z,Z0u3,lamdazu3,nzu3)*H(S,S0u3,lamdaSu3,nsu3)-ms*Yu(u3,u30,nu3)-ku3*u3
def fms(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return gms*H(h,1,1,1)*H(S,S0ms,lamdaSms,nsms)*H(A,A0s,lamdaAs,nAs)-ms*Ym(u3,u30,nu3)-kms*ms
def fmz(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return gmz*H(A,A0m,lamdaAm,nAm)*H(Z,Z0m,lamdaZm,nzm)*H(S,S0m,lamdaSm,nsm)-mz*Ym(u,u0,nu)-kmz*mz
def fu(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return gu*H(A,A0Au,lamdaAu,nAu)*H(Z,Z0u,lamdazu,nzu)*H(S,S0u,lamdaSu,nsu)*H(h,1,1,1)*H(A,1,1,1)-mz*Yu(u,u0,nu)-ku*u
def fA(mz,Z,S,ms,u,u3,h,A,rm,rn):
	R = rm+rn
	return ga*H(R,R0ra,lambdara,nra)*H(h,h0ha,lambdaha,nha)*H(A,A0aa,lambdaaa,naa)-ka*A
def fRMT(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return grm*H(A,A0aR,lambdaar,nar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*rm
def fRNOX(mz,Z,S,ms,u,u3,h,A,rm,rn):
	return grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*rn
def fH(mz,Z,S,ms,u,u3,h,A,rm,rn):
	R = rm+rn
	return gh*H(A,A0ah,lambdaah,nah)-kh*h*H(R,R0rh,lambdarh,nrh)*H(h,h0hh,lambdahh,nhh)
#############################################
#############################################
#DSargsF.ics = {'u':1000,'mz':1000,'S':200000,'ms':2000,'u3':4000}
#DSargsF.xdomain = {'u':[0,
#DSargsF.tdomain = [0,5000]
#DSargsF.algparams = {'init_step':0.1}
#DSargsF.pars['A']=300
#el= [0,0.2,0.3,0.35,0.4,0.6,1.]
#DSargsF.pars['lamdaAs']=el[i]
#############################################
#############################################


def getICS():
	## mz,Z,S,ms,u,u3,H,A,rm,rn
	return np.random.uniform(0,100000),np.random.uniform(0,100000),np.random.uniform(0,400000),np.random.uniform(0,100000),np.random.uniform(0,100000),np.random.uniform(0,100000),np.random.uniform(0,2000),np.random.uniform(0,2000),np.random.uniform(0,2000),np.random.uniform(0,2000)


def runSim(mz,Z,S,ms,u,u3,h,A,rm,rn,tmax=1000,dt=0.1):

	for  t in range(int(tmax/dt)):

		tmpa = fZ(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpb = fS(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpc = fu3(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpd = fms(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpe = fmz(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpf = fu(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpg = fA(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmph = fH(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpi = fRNOX(mz,Z,S,ms,u,u3,h,A,rm,rn)
		tmpj = fRMT(mz,Z,S,ms,u,u3,h,A,rm,rn)

		mz = mz + tmpe*dt
		Z = Z + tmpa*dt
		S = S + tmpb*dt
		ms = ms + tmpd*dt
		u = u + tmpf*dt
		u3 = u3 + tmpc*dt
		h = h + tmpg*dt
		A = A + tmpg*dt
		rn = rn + tmpi*dt
		rm = rm + tmpj*dt


		if (tmpa**2+tmpb**2+tmpc**2+tmpd**2+tmpe**2+tmpf**2+tmpg**2+tmph**2+tmpi**2+tmpj**2)<1:		
			return mz,Z,S,ms,u,u3,h,A,rm,rn
	return mz,Z,S,ms,u,u3,h,A,rm,rn
#######################
lamdaAs=1.
nAs=2
A0s=100.
nICS=5000
lamdaAu=1.
nAu=1
A0Au=100.
lamdaAm=1
nAm=2
A0m=100.

#DSargsF.pars['lamdaAs']=el[i]

l_rank ={}
for i in range(nprocs):
	l_rank[i]=[]
count=0
for el in [0,0.1,0.3,0.5,0.7,1.,2.]:
	l_rank[count%nprocs]+=[el]
	count+=1

for lamdaAm in l_rank[myrank]:
	file1 = open("ics_Full_AZ"+str(lamdaAm)+".txt","w")
	file2 = open("res_Full_AZ"+str(lamdaAm)+".txt","w")
	file1.write("mz,Z,S,ms,u,u3,h,A,rm,rn\n")
	file2.write("mz,Z,S,ms,u,u3,h,A,rm,rn\n")
	np.random.seed(3142931923)
	for i in range(nICS):
		mz,Z,S,ms,u,u3,h,A,rm,rn = getICS()
		file1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(mz,Z,S,ms,u,u3,h,A,rm,rn))
		mz,Z,S,ms,u,u3,h,A,rm,rn = runSim(mz,Z,S,ms,u,u3,h,A,rm,rn,5000,0.1)
		file2.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(mz,Z,S,ms,u,u3,h,A,rm,rn))

	file1.close()
	file2.close()
