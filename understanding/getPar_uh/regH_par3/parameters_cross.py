#######################################################
#######################################################
############## EMT circuit params ##############
############## for L>Y            #####################
#######################################################
gz = 100. 	##  molecules/hr
gu = 2100. 	##  molecules/hr
gmz = 11. 	##  molecules/hr

kz=0.1 		## 1/hr
ku = 0.05 	## 1/hr
kmz=0.5		## 1/hr

Z0u = 220000.    ## molecules/hr
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

ymi=[0,0.04, 0.2,1.0,1., 1.,1.]
yui=[0,0.005,0.05,0.5,0.5,0.5,0.5]
li = [1.,0.6,0.3,0.1,0.05,0.05,0.05]

ymih=[0.,0.,0.]#[0,0.04, 0.2]
lih = [1.,0.,0.]#[1.,0.6,0.3]

gms =90. # molecule/hr
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
nra= 4.
h0ha= 250.# nM
lambdaha= 0.1
nha= 1.
A0aa= 350. #nM
lambdaaa = 0.2
naa = 2.
ka = 0.2 #\hr

# for Rmt
grm =150. # uM/min
lambdaar = 0.25
A0aR = 350. #nM
nar = 2.
y =8.
gn = 0.2
h0hrm = 200. #nM
nhrm = 2.
A0rm = 150. #nM
narm = 4.
krm = 5. #\min

#for rnox
grn = 40. # uM/min
gn = 0.2
h0hrn = 250. #nM
nhrn = 2. ## maybe this is the issue
g1 = 5. 
A0rn = 150. #nM
g2 =  0.2
narn = 2.
krn = 5. #\min
narn = 2.

## for H
gh = 15. #nM/h
A0ah = 250. #nM
lambdaah = 0.1
nah = 1.
kh = 0.25 #\hour
R0rh = 300. #uM
lambdarh = 0.2 ## opposite because acting on degredation
nrh = 4.
h0hh = 80. # nM
lambdahh = 0.1 ## opposite because acting on degredation
nhh = 4.

nuh=2
