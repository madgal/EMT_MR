import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
from pylab import *
import copy
import math


from bif_functions import *


def new_reducefp(tmp_Ffp):
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
        finalArr+=[{'u':finalV[i][iu],'mz':finalV[i][imz],'u3':finalV[i][iu3],'S':finalV[i][iS],'ms':finalV[i][ims]}]
    return finalArr

DSargsF = PyDSTool.args(name='emtFull',checklevel=2)
DSargsF.fnspecs = {'Z':(['uX','mzX'],'gz*mzX*L(uX,u0,nu)/kz'),
		  'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
		  'M' : (['i','n','x','x0'],'(x/x0)**i/(1+(x/x0))**n'),
                  'indL': (['i'],'if(i==0,li0,if(i==1,li1,if(i==2,li2,if(i==3,li3,if(i==4,li4,if(i==5,li5,li6))))))'),
                  'indYm': (['i'],'if(i==0,ymi0,if(i==1,ymi1,if(i==2,ymi2,if(i==3,ymi3,if(i==4,ymi4,if(i==5,ymi5,ymi6))))))'),
                  'indYu': (['i'],'if(i==0,yui0,if(i==1,yui1,if(i==2,yui2,if(i==3,yui3,if(i==4,yui4,if(i==5,yui5,yui6))))))'),
                  'combination':(['k','n'],'special_gamma(n+1)/special_gamma(k+1)/special_gamma(n-k+1)'),
                  'combinationU':(['k','n'],'k*special_gamma(n+1)/special_gamma(k+1)/special_gamma(n-k+1)'),
		  'L' : (['X','X0','nX'],'sum(i,0,6,li[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Ym' : (['X','X0','nX'],'sum(i,0,6,ymi[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Yu' : (['X','X0','nX'],'sum(i,0,6,yui[i]*combinationU([i],nX)*M([i],nX,X,X0))'),
		  'L3' : (['X','X0','nX'],'sum(i,0,2,li[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Ym3' : (['X','X0','nX'],'sum(i,0,2,ymi[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Yu3' : (['X','X0','nX'],'sum(i,0,2,yui[i]*combinationU([i],nX)*M([i],nX,X,X0))'),
			}

DSargsF.varspecs = {'u':'gu*H(Z(u,mz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)-mz*Yu(u,u0,nu)-ku*u',
		   'mz':'gmz*H(Z(u,mz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)-ms*Ym(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(u,mz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                   'S':'gs*ms*L3(u3,u30,nu3)-ks*S'
		}
			
DSargsF.pars ={'gz':100.,'gu':2100.,'gms':90.,'gu3':1350.,'gs':100.,'gmz':11.,
		'ks':0.125,'ka':0.2,'kz':0.1,'ku':0.05,'kmz':0.5,'kms':0.5,'ku3':0.05,
		'Z0u':220000.,'I0m':50000.,'S0ms':200000.,'u30':10000.,'S0u3':300000.,
		'Z0u3':600000.,'Z0m':25000,'S0u':180000.,'S0m':180000.,'u0':10000.,
		'lamdazu':0.1,'lamdaSu':0.1,'lamdaIm':10.,'lamdaSms':0.1,'lamdazu3':0.2,'lamdaSu3':0.1,'lamdaZm':7.5,'lamdaSm':10.,
		'nu3':2,'nzu3':2,'nsms':1,'nIm':2,'nsu3':1,
		'nzu':3,'nsu':2,'nzm':2,'nsm':2,'nu':6,
		'I':50000.,
		'ymi0':0.,'ymi1':0.04,'ymi2':0.2,'ymi3':1.,'ymi4':1.,'ymi5':1.,'ymi6':1.,
		'yui0':0.,'yui1':0.005,'yui2':0.05,'yui3':0.5,'yui4':0.5,'yui5':0.5,'yui6':0.5,
		'li0':1.,'li1':0.6,'li2':0.3,'li3':0.1,'li4':0.05,'li5':0.05,'li6':0.05}
		

##############
##############
#####
#####
DSargsF.ics = {'mz':500,'S':200000,'ms':500,'u3':5000,'u':1000}
DSargsF.xdomain = {'u':[0,20000],'mz':[0,1000],'S':[0,250000],'ms':[0,1000],'u3':[0,20000]}
new_domains = {'u':[0,20000],'mz':[0,1000],'S':[0,250000],'ms':[0,1000],'u3':[0,20000]}
DSargsF.tdomain = [0,5000]
DSargsF.algparams = {'init_step':0.4}

######{'A': 436.47213160474945, 'h': 26.69259410719074}
######{'A': 327.50827624867554,'h': 292.59779290727334}
######{'A': 68.67768314578954, 'h': 480.4124463778733}
ode = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
Ffixed_points = pp.find_fixedpoints(ode,n=25,maxsearch=1e+5,eps=1e-10)

tmp_Ffp=[]
for i in range(len(Ffixed_points)):
    tmp_Ffp+=[list(Ffixed_points[i].values())]
tmp_Ffp = np.array(tmp_Ffp)
columns=np.array(list(Ffixed_points[0].keys()))
iu = np.argwhere(columns=='u')[0][0]
imz = np.argwhere(columns=='mz')[0][0]
iu3 = np.argwhere(columns=='u3')[0][0]
iS = np.argwhere(columns=='S')[0][0]
ims = np.argwhere(columns=='ms')[0][0]

fixed_points = new_reducefp(tmp_Ffp)


fileO = open("data/EMW_fp.txt",'w')
fileO.write('S,ms,mz,u,u3\n')
for i in range(len(fixed_points)):
	fileO.write("%s,%s,%s,%s,%s\n" %(  fixed_points[i]['S'],fixed_points[i]['ms'], fixed_points[i]['mz'], fixed_points[i]['u'], fixed_points[i]['u3']))

fileO.close()

### nullclines
solC_u = nullclines('u',DSargsF,fixed_points,new_domains,1000)
try:
	fileO = open("data/nullU_emt.txt",'w')
	fileO.write('S,ms,mz,u,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i]))
	fileO.close()
except:
	fileO = open("data/nullU_emt.txt",'w')
	fileO.write('S,ms,mz,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u3'][i]))
	fileO.close()

### nullclines
solC_u = nullclines('mz',DSargsF,fixed_points,new_domains,1000)
try:
	fileO = open("data/nullMZ_emt.txt",'w')
	fileO.write('S,ms,mz,u,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i]))
	fileO.close()
except:
	fileO = open("data/nullMZ_emt.txt",'w')
	fileO.write('S,ms,u,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['u'][i], solC_u['u3'][i]))
	fileO.close()
