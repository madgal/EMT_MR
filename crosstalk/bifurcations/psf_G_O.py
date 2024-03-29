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
        finalArr+=[{'u':finalV[i][iu],'mz':finalV[i][imz],'u3':finalV[i][iu3],'S':finalV[i][iS],'ms':finalV[i][ims],'mg':finalV[i][img],'mo':finalV[i][imo],'G':finalV[i][iG],'O':finalV[i][iO]}]
    return finalArr



DSargsF = PyDSTool.args(name='psfFull',checklevel=2)
DSargsF.fnspecs = {'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
                  'Z':(['mzX','uX','u0X','nuX','kzX','gzX'],'gzX*mzX*L(uX,u0X,nuX)/kzX'),
		  'M' : (['i','n','x','x0'],'(x/x0)**i/(1+(x/x0))**n'),
                  'indL': (['i'],'if(i==0,li0,if(i==1,li1,if(i==2,li2,if(i==3,li3,if(i==4,li4,if(i==5,li5,li6))))))'),
                  'indYm': (['i'],'if(i==0,ymi0,if(i==1,ymi1,if(i==2,ymi2,if(i==3,ymi3,if(i==4,ymi4,if(i==5,ymi5,ymi6))))))'),
                  'indYu': (['i'],'if(i==0,yui0,if(i==1,yui1,if(i==2,yui2,if(i==3,yui3,if(i==4,yui4,if(i==5,yui5,yui6))))))'),
                  'indYuh': (['i'],'if(i==0,yuih0,if(i==1,yuih1,if(i==2,yuih2,yuih2)))'),
                  'combination':(['k','n'],'special_gamma(n+1)/special_gamma(k+1)/special_gamma(n-k+1)'),
                  'combinationU':(['k','n'],'k*special_gamma(n+1)/special_gamma(k+1)/special_gamma(n-k+1)'),
		  'L' : (['X','X0','nX'],'sum(i,0,6,li[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Ym' : (['X','X0','nX'],'sum(i,0,6,ymi[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Yu' : (['X','X0','nX'],'sum(i,0,6,yui[i]*combinationU([i],nX)*M([i],nX,X,X0))'),
		  'L3' : (['X','X0','nX'],'sum(i,0,2,li[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Ym3' : (['X','X0','nX'],'sum(i,0,2,ymi[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Yu3' : (['X','X0','nX'],'sum(i,0,2,yui[i]*combinationU([i],nX)*M([i],nX,X,X0))'),
		  'Yuh' : (['X','X0','nX'],'sum(i,0,2,yuih[i]*combinationU([i],nX)*M([i],nX,X,X0))'),
			}

DSargsF.varspecs = {'u':'gu*H(Z(mz,u,u0,nu,kz,gz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)*H(O,o0u,nou,lou)-mz*Yu(u,u0,nu)-mh*Yuh(0,u0,nuh)-ku*u',
		   'mz':'gmz*H(Z(mz,u,u0,nu,kz,gz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)*H(A,A0m,nAm,lamdaAm)*H(O,o0z,noz,loz)*H(G,g0mz,ngmz,lgmz)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)*H(h,h0ms,nhms,lamdahms)*H(A,A0ms,nAms,lamdaAms)-ms*Ym(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(mz,u,u0,nu,kz,gz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                   'S':'gs*ms*L3(u3,u30,nu3)-ks*S',
                   'mg':'gmg*H(Z(mz,u,u0,nu,kz,gz),Z0mg,nzmg,lzmg)-kmg*mg',
                   'mo':'gmo*H(O,o0mo,nomo,lomo)*H(Z(mz,u,u0,nu,kz,gz),Z0mo,nzmo,lzmo)*H(G,G0mo,nGmo,lGmo)-kmo*mo',
                   'G':'gg*mg-kgg*G',
                   'O':'gov*mo-kov*O',
		}
			
DSargsF.pars ={'gz':100.,'gu':2100.,'gms':90.,'gu3':1350.,'gs':100.,'gmz':11.,
		'ks':0.125,'kz':0.1,'ku':0.05,'kmz':0.5,'kms':0.5,'ku3':0.05,
		'Z0u':220000.,'I0m':50000.,'S0ms':200000.,'u30':10000.,'S0u3':300000.,
		'Z0u3':600000.,'Z0m':25000,'S0u':180000.,'S0m':180000.,'u0':10000.,'h0u':200.,'A0u':300.,'A0m':300.,'A0ms':300.,
		'h0ms':200.,
		'lamdahu':1.,'lamdaAu':1.,'lamdaAms':1.,'lamdahms':1.,'lamdaAm':1.,
		'lamdazu':0.1,'lamdaSu':0.1,'lamdaSms':0.1,'lamdazu3':0.2,'lamdaSu3':0.1,'lamdaZm':7.5,'lamdaSm':10.,
		'nu3':2,'nzu3':2,'nsms':1,'nIm':2,'nsu3':1,
		'nzu':3,'nsu':2,'nzm':2,'nsm':2,'nu':6,'nuh':2,'nhu':1,'nAu':1,'nAm':2,'nAms':2,'nhms':2,
		'I':50000.,
		'ymi0':0.,'ymi1':0.04,'ymi2':0.2,'ymi3':1.,'ymi4':1.,'ymi5':1.,'ymi6':1.,
		'yui0':0.,'yui1':0.005,'yui2':0.05,'yui3':0.5,'yui4':0.5,'yui5':0.5,'yui6':0.5,
		'li0':1.,'li1':0.6,'li2':0.3,'li3':0.1,'li4':0.05,'li5':0.05,'li6':0.05,
		'ymih0':0.,'ymih1':0.04,'ymih2':0.2,
		'yuih0':0.,'yuih1':0.005,'yuih2':0.05,
		'lih0':1.,'lih1':0.6,'lih2':0.3 ,
		'h':0.,'A':0. ,'mh':0.,
                'gt':90.,'gg':200.,'gmg':22.,'gov':200.,'gmo':22.,
                'kt':0.4,'kgg':0.1,'kmg':0.5,'kov':0.1,'kmo':0.5,
		'nou':1,'noz':1,'nomo':2,'nzmo':1,'nzmg':3,'ngmz':1,
		'lou':0.1,'loz':0.1,'lomo':0.1,'lzmo':0.5,'lzmg':0.5,'lgmz':0.1,
                'o0u':250000.,'o0z':25000.,'o0mo':25000.,'Z0mo':10000.,'Z0mg':10000.,'g0mz':25000.,
                'lamdaIm':16.,
                'G0mo':25000., 'nGmo':2, 'lGmo':0.7 } 
##############
##############
#####
DSargsF.ics={'u':1000,'mz':400,'S':200000,'ms':500,'u3':4000,'G':60000.,'O':22000.,'mg':20.,'mo':10.}#'Z':1000}
DSargsF.xdomain = {'u':[0,15000],'mz':[0,500],'S':[0,500000],'ms':[0,1000],'u3':[0,20000],'G':[0,100000],'O':[0,30000],'mg':[0,50],'mo':[0,20]}
DSargsF.tdomain = [0,8000]
DSargsF.algparams = {'init_step':0.4}
new_domains = {'u':[0,15000],'mz':[0,500],'S':[0,500000],'ms':[0,1000],'u3':[0,20000],'G':[0,100000],'O':[0,30000],'mg':[0,50],'mo':[0,20]}


##### SET 1
ode = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+5,eps=1e-5)

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
img = np.argwhere(columns=='mg')[0][0]
imo = np.argwhere(columns=='mo')[0][0]
iG = np.argwhere(columns=='G')[0][0]
iO = np.argwhere(columns=='O')[0][0]

fixed_points = reduce_fp(tmp_Ffp)

fileO = open("data/PSF_fp.txt",'w')
fileO.write('S,ms,mz,u,u3,G,O,mg,mo\n')
for i in range(len(fixed_points)):
	fileO.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(  fixed_points[i]['S'],fixed_points[i]['ms'], fixed_points[i]['mz'], fixed_points[i]['u'], fixed_points[i]['u3'], fixed_points[i]['G'], fixed_points[i]['O'], fixed_points[i]['mg'], fixed_points[i]['mo']))

fileO.close()

### nullclines
solC_u = nullclines('u',DSargsF,fixed_points,new_domains,1000)
try:
	fileO = open("data/nullU_psf.txt",'w')
	fileO.write('S,ms,mz,u,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i], solC_u[i]['G'], solC_u[i]['O'], solC_u[i]['mg'], solC_u[i]['mo']))
	fileO.close()
except:
	fileO = open("data/nullU_psf.txt",'w')
	fileO.write('S,ms,mz,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u3'][i], solC_u[i]['G'], solC_u[i]['O'], solC_u[i]['mg'], solC_u[i]['mo']))
	fileO.close()

### nullclines
solC_u = nullclines('mz',DSargsF,fixed_points,new_domains,1000)
try:
	fileO = open("data/nullMZ_psf.txt",'w')
	fileO.write('S,ms,mz,u,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i], solC_u[i]['G'], solC_u[i]['O'], solC_u[i]['mg'], solC_u[i]['mo']))
	fileO.close()
except:
	fileO = open("data/nullMZ_psf.txt",'w')
	fileO.write('S,ms,u,u3\n')
	for i in range(len(solC_u['S'])):
		fileO.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['u'][i], solC_u['u3'][i], solC_u[i]['G'], solC_u[i]['O'], solC_u[i]['mg'], solC_u[i]['mo']))
	fileO.close()
