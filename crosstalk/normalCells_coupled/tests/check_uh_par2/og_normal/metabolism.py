import numpy as np
import matplotlib.pyplot as plt
import math

grn = 40.
grm =150.
ga = 30. 
gn = 0.2
g1 = 5. 
g2 =  0.2
gmh=10.
gh=15.
y =3.### =3 for normal  8. for cancer 

ka=0.2 
krm = 5. 
krn = 5.
kh=0.3#  0.3 for normal 0.25 for cancer

h0ha= 250.
A0aa= 350.
h0hrm = 200. 
h0hrn = 250. 
A0rn = 150. 
A0rm = 150. 
A0aR = 350. 
R0ra =  100. 
A0ah = 250.
R0rh = 300.
h0hh = 80. 

lambdarh = 0.2
lambdahh = 0.1
lambdaah = 0.1
lambdaar = 0.25
lambdara = 8.
lambdaha= 0.1
lambdaaa = 0.2

naa = 2
nra= 4
nha= 1
nar = 2
narm = 4
nhrm = 2
nhrn = 2
narn = 2
nah = 1
nrh = 4
nhh = 4

#############################
#############################
#############################
#############################


##  units

## [dA/dt] = nM/hr
## [dH/dt] = nM/hr
## [dRmt/dt] = uM/min
## [dRnox/dt] = uM/min

def H(X,X0,nX,lamdaX):
	return (lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX))
#######################################
def CompRmt(yyX,gnX,hX,h0rmtX,nhmtX,AX,A0rmtX,namtX):
	return yyX*(gnX+((AX/A0rmtX)**namtX))/(1.+((hX/h0rmtX)**nhmtX)+((AX/A0rmtX)**namtX))
def CompRn(g0X,hX,h0rnX,nhnX,ghnX,AX,A0rnX,ganX, nanxX):
	return (g0X+ghnX*((hX/h0rnX)**nhnX)+ganX*((AX/A0rnX)**nanxX))/(1.+((hX/h0rnX)**nhnX)+((AX/A0rnX)**nanxX))
def fA(R,h,A):
	return ga*H(R,R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A
def fRMT(A,h,Rmt):
	return grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt
def fRNOX(A,h,Rnox):
	return grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox
def fH(R,h,A):
	return gh*H(A,A0ah,nah,lambdaah)-kh*h*H(h,h0hh,nhh,lambdahh)*H(R,R0rh,nrh,lambdarh)

def run(Rmt,Rnox,A,h):
	dt=0.1
	totalT=10000

	Al,Rmtl,Rnoxl,hl = [],[],[],[]
	for i in range(0,int(totalT/dt)):
		R = Rmt+Rnox
		tmpA = fA(R,h,A)
		tmph = fH(R,h,A)
		tmprm = fRMT(A,h,Rmt)
		tmprn = fRNOX(A,h,Rnox)

		A = A+ tmpA*dt
		Rmt = Rmt + tmprm*dt
		Rnox = Rnox + tmprn*dt
		h = h+ tmph*dt

		#Al+=[A]
		#Rmtl+=[Rmt]
		#hl+=[h]
		#Rnoxl+=[Rnox]

		

	return [Rmt,Rnox,A,h]
	#return [Rmtl,Rnoxl,Al,hl]


fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
Rm=np.random.uniform(0,1000,1000)
Rn=np.random.uniform(0,1000,1000)
A=np.random.uniform(0,1000,1000)
h=np.random.uniform(0,1000,1000)

[Rmtl,Rnoxl,Al,hl] = run(Rm,Rn,A,h)
ax1.plot(Rmtl,'*')
ax2.plot(Rnoxl,'*')
ax3.plot(Al,'*')
ax4.plot(hl,'*')
ax1.set_title("Rmt")
ax2.set_title("Rnox")
ax3.set_title("A")
ax4.set_title("h")
plt.show()

plt.plot(Al,hl,'*')
print 'A',np.unique(Al)
print 'h', np.unique(hl)
plt.show()
