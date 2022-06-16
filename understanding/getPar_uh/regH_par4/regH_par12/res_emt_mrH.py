import matplotlib
import numpy as np
import math
matplotlib.rcParams.update({'font.size': 30})

from mpi4py import MPI
comm = MPI.COMM_WORLD
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()

from parameters_cross import *
gmh=10.
kmh=.1667
gh=1.5
kh=1.5

def H(X,X0,nX,lamdaX):
	return lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)

def M(i,n,x,x0):
	return (x/x0)**i/(1+(x/x0))**n
def combination(k,n):
	return math.factorial(n)/math.factorial(k)/math.factorial(n-k)
def L(X,X0,nX):
	total=0
	for i in range(nX+1):
		total+= li[i]*combination(i,nX)*M(i,nX,X,X0)
	return total

def Ym(X,X0,nX):
	total=0
	for i in range(nX+1):
		total+= ymi[i]*combination(i,nX)*M(i,nX,X,X0)
	return total
def Yu(X,X0,nX):
	total=0
	for i in range(nX+1):
		total+= i*yui[i]*combination(i,nX)*M(i,nX,X,X0)
	return total
def CompRmt(yyX,gnX,hX,h0rmtX,nhmtX,AX,A0rmtX,namtX):
	return yyX*(gnX+(AX/A0rmtX)**namtX)/(1.+(hX/h0rmtX)**nhmtX+(AX/A0rmtX)**namtX)
def CompRn(g0X,hX,h0rnX,nhnX,ghnX,AX,A0rnX,ganX,nanxX):
	return (g0X+ghnX*(hX/h0rnX)**nhnX+ganX*(AX/A0rnX)**nanxX)/(1.+(hX/h0rnX)**nhnX+(AX/A0rnX)**nanxX)
def R(R1,R2):
	return R1+R2


def uF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gu*H(Z,Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)-mz*Yu(u,u0,nu)-ku*u
def mzF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gmz*H(Z,Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)-mz*Ym(u,u0,nu)-kmz*mz
def ZF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gz*mz*L(u,u0,nu)-kz*Z
def msF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	I = 50000.
	return gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)-ms*Ym(u3,u30,nu3)-kms*ms
def u3F(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gu3*H(Z,Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu(u3,u30,nu3)-ku3*u3
def SF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gs*ms*L(u3,u30,nu3)-ks*S

def AF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A
def RmtF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt
def RnoxF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox
def hF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gh*H(A,A0ah,nah,lambdaah)-kh*h*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)
def hF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gh*mh*L(0.,u0,nuh)-kh*h
def mhF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):
	return gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ym(0.,u0,nuh)

def getICS():
	mz_range=2000
	S_range=250000
	ms_range=1000
	u3_range=20000
	Z_range=700000
	u_range=25000
	A_range=1000
	h_range=1000
	Rmt_range=1000
	Rnox_range=1000
	mh_range=1000
	#mz_range=100000 S_range=400000 ms_range=100000 u3_range=100000 Z_range=100000 u_range=100000 A_range=2000 h_range=2000 Rmt_range=2000 Rnox_range=2000
	return np.random.randint(0,u_range),np.random.randint(0,mz_range),np.random.randint(0,Z_range),np.random.randint(0,ms_range),np.random.randint(0,u3_range),np.random.randint(0,S_range),np.random.randint(0,A_range),np.random.randint(0,Rmt_range),np.random.randint(0,Rnox_range),np.random.randint(0,h_range),np.random.randint(0,mh_range)

def run(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):

	total=5000
	dt=0.1
	for i in range(int(total/dt)):

		utmp =     uF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		mztmp =    mzF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		ztmp =     ZF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		mstmp =    msF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		u3tmp =    u3F(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		Stmp =     SF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		Atmp =     AF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		rmttmp =   RmtF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		rnoxtmp =  RnoxF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		htmp =     hF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		mhtmp =     mhF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
		

		u =    dt*utmp +u
		mz =   dt*mztmp +mz
		Z =    dt*ztmp +Z
		ms =   dt*mstmp +ms
		u3 =   dt*u3tmp +u3
		S =    dt*Stmp +S
		A =    dt*Atmp +A
		Rmt =  dt*rmttmp +Rmt
		Rnox = dt*rnoxtmp +Rnox
		h =    dt*htmp +h
		mh =    dt*mhtmp +mh
	
	return [u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh]

##################
##################
##################

for totalSim in [1000]:#[100,1000,10000,50000]:
	timeL={}
	tempTotal=0
	for rank in range(1,nprocs-1):
		timeL[rank]=int(totalSim/(nprocs-1))
		tempTotal+=int(totalSim/(nprocs-1))

	timeL[nprocs-1]=totalSim-tempTotal
	
	if myrank!=0:
		data={'ics':[],'res':[]}
		for i in  range(timeL[myrank]):
			u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh = getICS()
			res = run(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
			data['ics']+=[[u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh]]
			data['res']+=[res]
		comm.send(data,dest=0)
		cont=comm.recv(source=0)

	if myrank==0:
		filei = open("emtmrH_ics_"+str(totalSim)+".txt","w")
		filei.write("u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh\n")
		fileo = open("emtmrH_res_"+str(totalSim)+".txt","w")
		fileo.write("u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh\n")

		for i in range(1,nprocs):
			data = comm.recv(source=i)
			for i in range(len(data['ics'])):
				[u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh] = data['ics'][i]
				res = data['res'][i]
				filei.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh))
				fileo.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(res[0],res[1],res[2],res[3],res[4],res[5],res[6],res[7],res[8],res[9],res[10]))
	
		fileo.close()
		filei.close()

		for i in range(1,nprocs):
			comm.send(True,dest=i)
