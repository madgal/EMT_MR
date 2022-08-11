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
y =3.
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

