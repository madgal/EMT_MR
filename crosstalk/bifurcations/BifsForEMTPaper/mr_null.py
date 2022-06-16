import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
from pylab import *
import copy
import math


from bif_functions import *


def reduce_fp(tmp_Ffp):
    #tmpb = copy.deepcopy(tmp_Ffp)
    finalV=[]
    i=0
    while i<len(tmp_Ffp) and i<50:
        dat= np.sum((tmp_Ffp[i]-tmp_Ffp[:])**2,axis=1)<100
        finalV+=[tmp_Ffp[i]]
        inds = np.argwhere(dat==1)[:,0]
        tmp_Ffp=np.delete(tmp_Ffp,inds,axis=0)
        i+=1

    finalV=np.array(finalV)
    
    finalArr = []
    for i in range(len(finalV)):
        finalArr+=[{'A':finalV[i][iA],'h':finalV[i][ih],'Rmt':finalV[i][imr],'Rnox':finalV[i][inr],'mh':finalV[i][imh]}]
    return finalArr


DSargsF = PyDSTool.args(name='MRFull',checklevel=2)
DSargsF.fnspecs = {'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
		  'M' : (['i','n','x','x0'],'(x/x0)**i/(1+(x/x0))**n'),
                  'indLh': (['i'],'if(i==0,lih0,if(i==1,lih1,if(i==2,lih2,lih2)))'),
                  'indYmh': (['i'],'if(i==0,ymih0,if(i==1,ymih1,if(i==2,ymih2,ymih2)))'),
                  'combination':(['k','n'],'special_gamma(n+1)/special_gamma(k+1)/special_gamma(n-k+1)'),
		  'Lh' : (['X','X0','nX'],'sum(i,0,2,lih[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Ymh' : (['X','X0','nX'],'sum(i,0,2,ymi[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'R':(['R1','R2'],'R1+R2'),
	 	  'CompRmt':(['yyX', 'gnX', 'hX', 'h0rmtX', 'nhmtX', 'AX', 'A0rmtX',  'namtX'],'yyX*(gnX+(AX/A0rmtX)**namtX)/(1.+(hX/h0rmtX)**nhmtX+(AX/A0rmtX)**namtX)'),
		  'CompRn':(['g0X', 'hX', 'h0rnX', 'nhnX', 'ghnX', 'AX', 'A0rnX', 'ganX',  'nanxX'],'(g0X+ghnX*(hX/h0rnX)**nhnX+ganX*(AX/A0rnX)**nanxX)/(1.+(hX/h0rnX)**nhnX+(AX/A0rnX)**nanxX)')
			}

DSargsF.varspecs = { 'A':'ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A',
		   'h':'gh*mh*Lh(0.,u0,nuh)-kh*h',
		   'Rmt':'grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt*H(u3,nu30rm,n3m,lamda3m)',
		   'Rnox':'grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox*H(u3,nu30rn,n3n,lamda3n)',
 		   'mh':'gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ymh(0.,u0,nuh)',
		}
			
DSargsF.pars ={'gz':100.,'gu':2100.,'grn':40.,'grm':150.,'gms':90.,'gu3':1350.,'gs':100.,'ga':30.,'gn':0.2,'g1':5.,'g2':0.2,'gmh':10.,'gh':1.5,'gmz':11.,'y':8.,
		'ks':0.125,'ka':0.2,'kz':0.1,'ku':0.05,'kmz':0.5,'kms':0.5,'ku3':0.05,'krm':5.,'krn':5.,'kmh':.143,'kh':1.75,
		'Z0u':220000.,'I0m':50000.,'S0ms':200000.,'u30':10000.,'S0u3':300000.,'h0ha':250.,'A0aa':350.,'h0hrm':200.,'h0hrn':250.,'A0rn':150.,'A0rm':150.,
		'A0aR':350.,'R0ra':100.,'Z0u3':600000.,'Z0m':25000,'S0u':180000.,'S0m':180000.,'u0':10000.,'h0u':200.,'A0u':300.,'A0m':300.,'A0ms':300.,
		'h0ms':200.,'nu30rn':10000.,'nu30rm':10000.,'A0ah':250.,'R0rh':300.,'h0hh':80.,
		'lamdahu':1.,'lamdaAu':1.,'lamda3n':1.,'lamda3m':1.,'lamdaAms':1.,'lamdahms':1.,'lamdaAm':1.,
		'lambdarh':0.2,'lambdahh':0.1,'lambdaah':0.1,'lambdaar':0.25,'lambdara':8.,'lambdaha':0.1,'lambdaaa':0.2,
		'lamdazu':0.1,'lamdaSu':0.1,'lamdaIm':10.,'lamdaSms':0.1,'lamdazu3':0.2,'lamdaSu3':0.1,'lamdaZm':7.5,'lamdaSm':10.,
		'naa':2,'nu3':2,'nzu3':2,'nsms':1,'nIm':2,'nsu3':1,'nra':4,'nha':1,'nar':2,'narm':4,'nhrm':2,'nhrn':2,
		'nzu':3,'nsu':2,'nzm':2,'nsm':2,'nu':6,'narn':2,'nah':1,'nrh':4,'nhh':4,'nuh':2,'nhu':1,'nAu':1,'nAm':2,'nAms':2,'nhms':2,'n3m':3,'n3n':2,
		'I':50000,
		'ymi0':0.,'ymi1':0.04,'ymi2':0.2,'ymi3':1.,'ymi4':1.,'ymi5':1.,'ymi6':1.,
		'yui0':0.,'yui1':0.005,'yui2':0.05,'yui3':0.5,'yui4':0.5,'yui5':0.5,'yui6':0.5,
		'li0':1.,'li1':0.6,'li2':0.3,'li3':0.05,'li4':0.05,'li5':0.05,'li6':0.05,
		'ymih0':0.,'ymih1':0.04,'ymih2':0.2,
		'yuih0':0.,'yuih1':0.005,'yuih2':0.05,
		'lih0':1.,'lih1':0.6,'lih2':0.3 ,
		'u':0,'u3':0}

##############
##############
#####
[{'A': 384.99831311027367, 'h': 19.91610769590981, 'Rmt': 104.15145578851006, 'Rnox': 1.6320911843192776, 'mh': 29.210291287334385},
{'A': 150.73497696093543, 'h': 266.7233481040565, 'Rmt': 51.020662041702394, 'Rnox': 15.484389738296143, 'mh': 391.1942438859495},
{'A': 79.07984883904831, 'h': 327.83410134719895, 'Rmt': 12.775080377329802, 'Rnox': 23.62897002100897, 'mh': 480.8233486425585}]

DSargsF.ics = {'h':100,'A':100,'Rnox':30,'Rmt':100,'mh':100}
DSargsF.xdomain = {'A':[0,500],'h':[0,500],'Rmt':[0,150],'Rnox':[0,50],'mh':[0,600]}
DSargsF.tdomain = [0,9000]
DSargsF.algparams = {'init_step':0.1}
new_domains = {'A':[0,500],'h':[0,500],'Rmt':[0,150],'Rnox':[0,50],'mh':[0,600]}
	
ode = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+8,eps=1e-10)

tmp_Ffp=[]
for i in range(len(Ffixed_points)):
    tmp_Ffp+=[list(Ffixed_points[i].values())]
tmp_Ffp = np.array(tmp_Ffp)
columns=np.array(list(Ffixed_points[0].keys()))
iA = np.argwhere(columns=='A')[0][0]
ih = np.argwhere(columns=='h')[0][0]
imh = np.argwhere(columns=='mh')[0][0]
imr = np.argwhere(columns=='Rmt')[0][0]
inr = np.argwhere(columns=='Rnox')[0][0]

fixed_points = reduce_fp(tmp_Ffp)

fileO = open("data/MR_fp.txt",'w')
fileO.write('A,h,Rmt,Rnox,mh\n')
for i in range(len(fixed_points)):
	fileO.write("%s,%s,%s,%s,%s\n" %(  fixed_points[i]['A'],fixed_points[i]['h'], fixed_points[i]['Rmt'], fixed_points[i]['Rnox'], fixed_points[i]['mh']))

fileO.close()

### nullclines
solC_u = nullclines('A',DSargsF,fixed_points,new_domains,2000,step=1e+1)
try:
	fileO = open("data/nullA_mr.txt",'w')
	fileO.write('A,h,Rmt,Rnox,mh\n')
	for i in range(len(solC_u['mh'])):
		fileO.write("%s,%s,%s,%s,%s\n" %( solC_u['A'][i], solC_u['h'][i], solC_u['Rmt'][i], solC_u['Rnox'][i], solC_u['mh'][i]))
	fileO.close()
except:
	fileO = open("data/nullA_mr.txt",'w')
	fileO.write('h,Rmt,Rnox,mh\n')
	for i in range(len(solC_u['mh'])):
		fileO.write("%s,%s,%s,%s\n" %(solC_u['h'][i], solC_u['Rmt'][i], solC_u['Rnox'][i], solC_u['mh'][i]))
	fileO.close()

### nullclines
solC_u = nullclines('h',DSargsF,fixed_points,new_domains,2000,step=1e+1)
try:
	fileO = open("data/nullh_mr.txt",'w')
	fileO.write('A,h,Rmt,Rnox,mh\n')
	for i in range(len(solC_u['mh'])):
		fileO.write("%s,%s,%s,%s,%s\n" %( solC_u['A'][i], solC_u['h'][i], solC_u['Rmt'][i], solC_u['Rnox'][i], solC_u['mh'][i]))
	fileO.close()
except:
	fileO = open("data/nullh_mr.txt",'w')
	fileO.write('A,Rmt,Rnox,mh\n')
	for i in range(len(solC_u['mh'])):
		fileO.write("%s,%s,%s,%s\n" %( solC_u['A'][i], solC_u['Rmt'][i], solC_u['Rnox'][i], solC_u['mh'][i]))
	fileO.close()
