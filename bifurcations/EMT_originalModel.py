import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
from pylab import *
import copy
import math
import matplotlib
matplotlib.rcParams.update({'font.size': 30})
import random as rand
import json


from parameters_metabolism_miRNA import *

def stability(fp,DSarg,eps=0.1):
	DSargs = copy.deepcopy(DSarg)
	ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
	out=[]
	for i in range(len(fp)):
		X={}
		stable=True
		for k in fp[i].keys():
			X[k] = fp[i][k]*(1+eps*rand.sample(list([-1,1]),1)[0])
			if X[k]<0:
				X[k]=0
			elif X[k]>DSargs.xdomain[k][1]:
				X[k] = DSargs.xdomain[k][1]
	
		ode.set(ics=X)
		traj=ode.compute('traj')
		X = traj.sample()[-1]
		for k in fp[0].keys():
			if np.abs(X[k]-fp[i][k]) >eps*fp[i][k]:
				stable=False
		if stable:
			out+=['S']
		else:
			out+=['U']
	return out
def reduce_fp(fixed_points):
        final = []
        for el in fixed_points:
                if len(final)>0:
                        addI=True
                        for item in final:
                                if np.abs(el['ms']-item['ms'])<0.001:
                                        if np.abs(el['mz']-item['mz'])<0.001:
                                                if np.abs(el['S']-item['S'])<0.001:
                                                	if np.abs(el['u']-item['u'])<0.001:
                                                		if np.abs(el['u3']-item['u3'])<0.001:
                                                       			addI=False
                        if addI:
                                final+=[el]
                else:
                        final+=[el]
        return tuple(final)

def nullclines(freepar,DSarg,fixed_points,new_domains):
    #### Setup Nullcline for N2
    DSargs_X = copy.deepcopy(DSarg)
    del DSargs_X.varspecs[freepar]
    DSargs_X.pars[freepar] = fixed_points[0][freepar]
    DSargs_X.ics = copy.deepcopy(fixed_points[0])
    del DSargs_X.ics[freepar]
    DSargs_X.pdomain = {freepar:new_domains[freepar]}
    DSargs_X.xdomain = copy.deepcopy(new_domains)
    del DSargs_X.xdomain[freepar]
    ode_X = PyDSTool.Generator.Vode_ODEsystem(DSargs_X)

    ### set arguments for numerical solution to ODE
    PCargs_X = PyDSTool.args(name='nullcline',type='EP-C') ## equilibrium point curve labeled ND-Test
    PCargs_X.freepars = [freepar]
    ## following 3 parameters set through trial and error
    PCargs_X.MaxNumPoints = 1000## max number for each call of forward and backward
    PCargs_X.MaxStepSize = 1e+2
    PCargs_X.MinStepSize = 1e-1
    PCargs_X.StepSize = 1e+1
    PCargs_X.StopAtPoints = ['B'] ##  stops searching once this point is hit
    PCargs_X.LocBifPoints = 'all' ##
    PCargs_X.SaveEigen = True

    PyCont_X = PyDSTool.ContClass(ode_X) # Set up continuation class
    PyCont_X.newCurve(PCargs_X)
    PyCont_X['nullcline'].forward()
    PyCont_X['nullcline'].backward()
    return PyCont_X['nullcline'].sol

##  units
DSargs = PyDSTool.args(name='emtFull',checklevel=2)
DSargs.fnspecs = {'Z':(['uX','mzX'],'gz*mzX*L(uX,u0,nu)/kz'),
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

DSargs.varspecs = {'u':'gu*H(Z(u,mz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)-mz*Yu(u,u0,nu)-ku*u',
		   'mz':'gmz*H(Z(u,mz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)-ms*Ym(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(u,mz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                   'S':'gs*ms*L3(u3,u30,nu3)-ks*S'
		}
			
DSargs.pars ={'gz':100.,'gu':2100.,'gms':90.,'gu3':1350.,'gs':100.,'gmz':11.,
		'ks':0.125,'ka':0.2,'kz':0.1,'ku':0.05,'kmz':0.5,'kms':0.5,'ku3':0.05,
		'Z0u':220000.,'I0m':50000.,'S0ms':200000.,'u30':10000.,'S0u3':300000.,
		'Z0u3':600000.,'Z0m':25000,'S0u':180000.,'S0m':180000.,'u0':10000.,
		'lamdazu':0.1,'lamdaSu':0.85,'lamdaIm':10.,'lamdaSms':0.1,'lamdazu3':0.2,'lamdaSu3':0.1,'lamdaZm':7.5,'lamdaSm':17.,
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
DSargs.ics = {'mz':1000,'S':200000,'ms':2000,'u3':4000,'u':1000}
DSargs.xdomain = {'u':[0,100000],'mz':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000]}
DSargs.tdomain = [0,5000]
DSargs.algparams = {'init_step':0.4}

ode=PyDSTool.Generator.Vode_ODEsystem(DSargs)

Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
fixed_points = reduce_fp(Ffixed_points)
print(fixed_points)
stabilities = stability(fixed_points,DSargs)
new_domains = {'u':[0,100000],'mz':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000]}


solC_U = nullclines('u',DSargs,fixed_points,new_domains)
solC_Z = nullclines('mz',DSargs,fixed_points,new_domains)


results={'fixed_points':fixed_points,'nullZ_Z':list(solC_a['mz']),'nullZ_U':list(solC_a['u']),'nullU_Z':list(solC_h['mz']),'nullU_U':list(solC_h['u']),'stability':stabilities}

json_data = json.dumps(results,indent=4)
with open("data_EMT_nullcline_original.json",'w') as outfile:
	outfile.write(json_data)
