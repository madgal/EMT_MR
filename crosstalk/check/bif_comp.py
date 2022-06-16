import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

gz = 100. 	
gu = 2100. 	
grn = 40.
grm =150.
gms =90. 
gu3 = 1350. 
gs = 100.
ga = 30. 
gn = 0.2
g1 = 5. 
g2 =  0.2
gmh=10.
gh=1.5
gmz = 11. 	
y =8.

ks = 0.125
ka = 0.2 
kz=0.1 		
ku = 0.05 	
kmz=0.5		
kms =0.5 
ku3 =0.05
krm = 5. 
krn = 5.
kmh=.143
kh=1.75

ymi=[0,0.04, 0.2,1.0,1., 1.,1.]
yui=[0,0.005,0.05,0.5,0.5,0.5,0.5]
li = [1.,0.6,0.3,0.1,0.05,0.05,0.05]
ymih=[0.,0.,0.]#[0,0.04, 0.2}
yuih=[0,0.005,0.05]
lih = [1.,0.,0.]#[1.,0.6,0.3]

Z0u = 220000.   
I0m = 50000. 
S0ms = 200000.
u30 = 10000. 
S0u3 = 300000. 
h0ha= 250.
A0aa= 350.
h0hrm = 200. 
h0hrn = 250. 
A0rn = 150. 
A0rm = 150. 
A0aR = 350. 
R0ra =  100. 
Z0u3 = 600000. 
Z0m = 25000.	
S0u=180000.	
S0m=180000.	
u0=10000.	
h0u = 200.
A0u=300.
A0m =300.
A0ms= 300.
h0ms = 200.
nu30rn = 10000.
nu30rm = 10000.
A0ah = 250.
R0rh = 300.
h0hh = 80. 

lamdahu = 1.  
lamdaAu = 1.
lamda3n = 1.
lamda3m= 1.
lamdaAms = 1.  
lamdahms = 1. 
lamdaAm=1.

lambdarh = 0.2
lambdahh = 0.1
lambdaah = 0.1
lambdaar = 0.25
lambdara = 8.
lambdaha= 0.1
lambdaaa = 0.2

lamdazu=0.1	
lamdaSu=0.1	
lamdaIm =10.
lamdaSms = 0.1
lamdazu3 = 0.2
lamdaSu3 = 0.1
lamdaZm=7.5	
lamdaSm=10.	

naa = 2
nu3 = 2
nzu3 = 2
nsms = 1
nIm = 2
nsu3 = 1
nra= 4
nha= 1
nar = 2
narm = 4
nhrm = 2
nhrn = 2
nzu =3
nsu=2
nzm=2
nsm=2
nu=6
narn = 2
nah = 1
nrh = 4
nhh = 4
nuh=2
nhu = 1
nAu = 1
nAm=2
nAms =2
nhms = 2 
n3m=3
n3n=2

def H( X, X0, nX, lamdaX):
	return lamdaX+(1.-lamdaX)/(1.+pow((X/X0),nX))

def M( i, n, x, x0):
	return pow((x/x0),i)/pow((1+(x/x0)),n)

def L( X, X0, nX, liX):
	total=0
	for i in range(nX+1):
		total= total+ liX[i]*combination(i,nX)*M(i,nX,X,X0)
	return total
	
def Ym( X, X0, nX, ymiX):
	total=0
	for i in range(nX+1):
		total=total+ ymiX[i]*combination(i,nX)*M(i,nX,X,X0)
	return total
	
def Yu( X, X0, nX, yuiX):
	total=0
	for i in range(nX+1):
		total= total+i*yuiX[i]*combination(i,nX)*M(i,nX,X,X0)
	return total
def CompRmt( yyX, gnX, hX, h0rmtX, nhmtX, AX, A0rmtX,  namtX):
	return yyX*(gnX+pow((AX/A0rmtX),namtX))/(1.+pow((hX/h0rmtX),nhmtX)+pow((AX/A0rmtX),namtX))
def CompRn( g0X, hX, h0rnX, nhnX, ghnX, AX, A0rnX, ganX,  nanxX):
	return (g0X+ghnX*pow((hX/h0rnX),nhnX)+ganX*pow((AX/A0rnX),nanxX))/(1.+pow((hX/h0rnX),nhnX)+pow((AX/A0rnX),nanxX))

def R( R1, R2):
	return R1+R2
def combination( k,  n):
	return math.factorial(n)/math.factorial(k)/math.factorial(n-k)

def uF( u, mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh):
	return gu*H(Z,Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)-mz*Yu(u,u0,nu,yui)-mh*Yu(0,u0,nuh,yuih)-ku*u
def mzF( u, mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh):
	return gmz*H(Z,Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)*H(A,A0m,nAm,lamdaAm)-mz*Ym(u,u0,nu,ymi)-kmz*mz
def AF( u, mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh):
	return ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A
def hF( u, mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh):
	return gh*mh*L(0.,u0,nuh,lih)-kh*h
def mhF( u, mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh):
	return gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ym(0.,u0,nuh,ymih)


def getICS(size):
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
	mh_range=150#1000

	ics =[np.random.randint(0,u_range,size),np.random.randint(0,mz_range,size),np.random.randint(0,Z_range,size),np.random.randint(0,ms_range,size),np.random.randint(0,u3_range,size),np.random.randint(0,S_range,size),np.random.randint(0,A_range,size),np.random.randint(0,Rmt_range,size),np.random.randint(0,Rnox_range,size),np.random.randint(0,h_range,size),np.random.randint(0,mh_range,size)]

	return ics

def run(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh):

	u =     uF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
	mz =    mzF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
	A=     AF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
	h =     hF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
	mh =     mhF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)
	res = {'u':u,'mz':mz,'A':A,'h':h,'mh':mh}
	return res

seeds = pd.read_csv("seedsForSims.txt")
[u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh] = getICS(100000)
res = run(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh)

plt.scatter(A,h,c=res['A'],cmap='seismic',marker='.')
plt.colorbar()
plt.show()
plt.scatter(mh,h,c=res['h'],cmap='seismic',marker='.')
plt.colorbar()
plt.show()

plt.scatter(A,mh,c=res['A'],cmap='seismic',marker='.')
plt.colorbar()
plt.show()
plt.scatter(A,mh,c=res['mh'],cmap='seismic',marker='.')
plt.colorbar()
plt.show()

plt.scatter(u,mz,c=res['u'],cmap='seismic',marker='.')
plt.colorbar()
plt.show()
plt.scatter(u,mz,c=res['mz'],cmap='seismic',marker='.')
plt.colorbar()
plt.show()


