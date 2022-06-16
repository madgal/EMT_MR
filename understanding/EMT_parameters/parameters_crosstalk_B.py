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
gs = 100.
ks = 0.125
#######################################################
############## Metabolism circuit params ##############
#######################################################
## for AMPK
ga = 30. # n</h
R0ra =  100. #uM
lambdara = 8.
nra= 4
h0ha= 250.# nM
lambdaha= 0.1
nha= 1
A0aa= 350. #nM
lambdaaa = 0.2
naa = 2
ka = 0.2 #\hr

# for Rmt
grm =150. # uM/min
lambdaar = 0.25
A0aR = 350. #nM
nar = 2
y =3.
gn = 0.2
h0hrm = 200. #nM
nhrm = 2
A0rm = 150. #nM
narm = 4
krm = 5. #\min

#for rnox
grn = 40. # uM/min
gn = 0.2
h0hrn = 250. #nM
nhrn = 2 ## maybe this is the issue
g1 = 5. 
A0rn = 150. #nM
g2 =  0.2
narn = 2
krn = 5. #\min
narn = 2

## for H
gh = 15. #nM/h
A0ah = 250. #nM
lambdaah = 0.1
nah = 1
kh = 0.25 #\hour
R0rh = 300. #uM
lambdarh = 0.2 ## opposite because acting on degredation
nrh = 4
h0hh = 80. # nM
lambdahh = 0.1 ## opposite because acting on degredation
nhh = 4


#### 
h0hs=200#  this is just a guess not sure how to get the exact value
lamdahs = 7. # fix
nhs = 2 #  2 binding sites based on 10.1074/jbc.M115.636944
h0=-1.
A0=-1.
I0Is=50000.
lambdaIs=10.
nIs=2


h0u=1.
lamdahu=1.
nhu=2
A0u=1.
lamdaAu=1.
nAu=2
A0m=1.
lamdaAm=1.
nAm=2
gmh= 10 # nM/h?? 1.1
kmh=.5 #/hr


nuh=1
u0uh=200 # fix
lambdauh=1. # fix

A0s=1.
lamdaAs=1.
nAs=8.
