#######################################################
#######################################################
############## EMT circuit params ##############
############## for L>Y            #####################
#######################################################
gz = 100. 	##  molecules/hr
gu = 2100. 	##  molecules/hr
gmz = 11.#12.5  	##  molecules/hr

kz=0.1 		## 1/hr
ku = 0.05 	## 1/hr
kmz=0.5		## 1/hr

Z0u = 220000.# 280000.    ## molecules/hr
Z0m = 25000.	## molecules/hr
h0=-1.
A0=-1.
S0u=180000.	## molecules/hr
S0m=180000.	## molecules/hr
u0=10000.	## molecules/hr


nzu =3		##
nhu=1#fix
nAu=1#fix
nsu=2		##
nzm=2		##
nhm=1#fix
nAm=1#fix
nsm=2		##
nu=6		##

lamdazu=0.1		##
lamdaSu=0.1		##
lamdaZm=7.5		##
lamdaSm=10.		##
lamdahu=1.
lamdaAu=1.
lamdahm=1.
lamdaAm=1.

yui=[0,0.005,0.05,0.5,0.5,0.5,0.5]
li = [1.,0.6,0.3,0.1,0.05,0.05,0.05]
li0  = 1.
li1  = 0.6
li2  = 0.3
li3  = 0.1
li4  = 0.05
li5  = 0.05
li6  = 0.05

yui0=0
yui1=0.005
yui2=0.05
yui3=0.5
yui4=0.5
yui5=0.5
yui6=0.5
ymi0= 0
ymi1= 0.04
ymi2= 0.2
ymi3= 1.0
ymi4= 1.
ymi5= 1.
ymi6= 1.

ymi=[0,0.04, 0.2,1.0,1., 1.,1.]
'''
gu = #1500 for L, 1500 for Y, 2900 for L~Y, 2100 for L>Y, 1550 for L<Y
gmz =  #12.5 for L, 15 for Y, 30 for L~Y, 11 for L>Y, 20 for L<Y
Z0u=#200K for L, 200K for Y, 300K for L~Y, 220K for L>Y, 250K for L<Y
Z0m=#50K for L, 50K for Y, 30K for L~Y, 25K for L>Y, 50K for L<Y
#LOnly
#li = [1.,0.5,0.2,0.02,0.02,0.02,0.02]
#Yonly
#ymi = [0,0.3,1.5,7.5,7.5,7.5,7.5]
#both equal
#li = [1.,0.7,0.5,0.1,0.05,0.05,0.05]
#ymi = [0,0.1,0.5,2.5,2.5,2.5,2.5]
'''

gms = 90. # molecule/hr
I0m = 50000. #Molecule
nIm = 2
lamdaIm =10.
S0ms = 200000.# mol
nsms = 1
lamdaSms = 0.1
u30 = 10000. # mol
nu3 = 2
kms =0.5 #/molecule/hr
gu3 = 1350. #molecule/hr
Z0u3 = 600000. # molecule
nzu3 = 2
lamdazu3 = 0.2
S0u3 = 300000. # molecule
nsu3 = 1
lamdaSu3 = 0.1
ku3 =0.05
gs = 9000.
ks = 0.0625

S0S =200000 
nss=1
lamdass=0.1


h0hs=200#  this is just a guess not sure how to get the exact value
lambdahs = 7. # fix
nhs = 2 #  2 binding sites based on 10.1074/jbc.M115.636944
