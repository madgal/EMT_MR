import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
import matplotlib
from pylab import *
import copy
import math
matplotlib.rcParams.update({'font.size': 30})
from matplotlib.lines import Line2D


from mpi4py import MPI
myrank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()

if nprocs <5:
	print("Need  at least 5")
	exit()

from parameters_emt_B import *

def getCont(DSarg,fixed_points,maxnum,step,fp='I'):
        DSargs = copy.deepcopy(DSarg)
        freepar=fp
        DSargs.pdomain={fp:[0,2500]}
        odex = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        odex.set(ics=fixed_points[0])
        PCargs = PyDSTool.args(name='EMTest',type='EP-C')
        PCargs.freepars =[fp]
        PCargs.MaxNumPoints = maxnum ## max number for each call of forward and backward
        PCargs.MaxStepSize = 1e+2
        PCargs.MinStepSize = 1e-2
        PCargs.StepSize = step
        PCargs.StopAtPoints = ['B'] ##  stops searching once this point is hit
        PCargs.LocBifPoints = 'all' ##
        PCargs.SaveEigen = True

        PyCont = PyDSTool.ContClass(odex) # Set up continuation class
        PyCont.newCurve(PCargs)
        PyCont['EMTest'].forward()
        PyCont['EMTest'].backward()

        return PyCont

def reduce_fp(fixed_points):

	final = []
	for el in fixed_points:
		if len(final)>0:
			addI=True
			for item in final:
				if np.abs(el['u']-item['u'])<0.001:
					if np.abs(el['mz']-item['mz'])<0.001:
						addI=False
			if addI:	
				final+=[el]
				
		else:
			final+=[el]

	return tuple(final)


DSargsF = PyDSTool.args(name='emtFull',checklevel=2)
DSargsF.fnspecs = {'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
		  'Z':(['mzX','uX','u0X','nuX','kzX','gzX'],'gzX*mzX*L(uX,u0X,nuX)/kzX'),	
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
		  'Yu3' : (['X','X0','nX'],'sum(i,0,2,yui[i]*combinationU([i],nX)*M([i],nX,X,X0))')} 

DSargsF.pars = {	'gz':gz, 'u0':u0,'nu':nu,'kz':kz,
                'li0':li0, 'li1':li1, 'li2':li2, 'li3':li3, 'li4':li4, 'li5':li5, 'li6':li6,
                'ymi0':ymi0, 'ymi1':ymi1, 'ymi2':ymi2, 'ymi3':ymi3, 'ymi4':ymi4, 'ymi5':ymi5, 'ymi6':ymi6,
                'yui0':yui0, 'yui1':yui1, 'yui2':yui2, 'yui3':yui3, 'yui4':yui4, 'yui5':yui5, 'yui6':yui6,
		'gu':gu, 'Z0u':Z0u,'nzu':nzu, 'lamdazu':lamdazu,'S0u':S0u,'nsu':nsu,'lamdaSu':lamdaSu,'u0':u0,'nu':nu,'ku':ku,
		'gmz':gmz, 'Z0m':Z0m,'nzm':nzm, 'lamdaZm':lamdaZm,'S0m':S0m,'nsm':nsm,'lamdaSm':lamdaSm,'u0':u0,'nu':nu,'kmz':kmz,
		'gms':gms, 'I0m':I0m, 'nIm':nIm, 'lamdaIm':lamdaIm, 'S0ms':S0ms, 'nsms':nsms, 'lamdaSms':lamdaSms, 'u30':u30,
                'nu3':nu3, 'kms':kms, 'gu3':gu3, 'Z0u3':Z0u3, 'nzu3':nzu3, 'lamdazu3':lamdazu3, 'S0u3':S0u3, 'nsu3':nsu3,
		'lamdaSu3':lamdaSu3, 'ku3':ku3, 'gs':gs, 'ks':ks, 'lamdahs':7., 'nhs':2., 'h0hs':200. ,'h':300 }
	
DSargsF.varspecs = { 'u':'gu*H(Z(mz,u,u0,nu,kz,gz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)-mz*Yu(u,u0,nu)-ku*u',
		   'mz':'gmz*H(Z(mz,u,u0,nu,kz,gz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(h,h0hs,nhs,lamdahs)*H(S,S0ms,nsms,lamdaSms)-ms*Ym3(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(mz,u,u0,nu,kz,gz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                    'S':'gs*ms*L3(u3,u30,nu3)-ks*S'

		}
			
DSargsF.ics = {'u':1000,'mz':1000,'S':200000,'ms':2000,'u3':4000}
DSargsF.xdomain = {'u':[0,100000],'mz':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000]}
DSargsF.tdomain = [0,5000]
DSargsF.algparams = {'init_step':0.1}
DSargsF.pars['h']=300

##############
##############
##############
def plotCC(x,y):
	cl= ['b','r','k','g','m','y']
	el= [6.5,6.8,7.,7.2,7.5]
	ll = []
	labs=[]
	minv,maxv=1000000,0
	for i in range(len(el)):
		DSargsF.pars['lamdahs']=el[i]
		odeF = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
		Ffixed_points = pp.find_fixedpoints(odeF,n=5,maxsearch=1e+4,eps=1e-10)
		fixed_pointsF = reduce_fp(Ffixed_points)
		PyContF=getCont(DSargsF,fixed_pointsF,2000,1e+0,fp=x)
		PyContF.display([x,y],stability=True,color=cl[i],linewidth=2)
		PyContF.plot.toggleLabels("off")

		solV = PyContF['EMTest'].sol
		if np.min(solV[y])<minv:
			minv = np.min(solV[y])
		if np.max(solV[y])>maxv:
			maxv = np.max(solV[y])

		labs+=[el[i]]
		ll +=[Line2D([0],[0],color=cl[i],lw=2)]

	return [ll,labs,maxv,minv]


#####
#####
if myrank==0:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('h','mz')
	plt.title('')
	plt.xlim(0,2500)
	plt.ylim(0,maxv)
	plt.xlabel("Hif1 (molecules)",fontsize=20)
	plt.ylabel("Zeb mRNA (molecules)",fontsize=20)
	plt.legend(ll,labs,fontsize=20,frameon=False,bbox_to_anchor=(1.04,1))
	plt.savefig("EM_MR_hmz.png",bbox_inches='tight')
	#plt.show()
	plt.close()

if myrank==1:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('h','S')
	plt.title('')
	plt.xlim(0,2500)
	plt.ylim(0,maxv)
	plt.xlabel("Hif1 (molecules)",fontsize=20)
	plt.ylabel("Snail (molecules)",fontsize=20)
	plt.legend(ll,labs,fontsize=20,frameon=False,bbox_to_anchor=(1.04,1))
	plt.savefig("EM_MR_hs.png",bbox_inches='tight')
	#plt.show()
	plt.close()

if myrank==2:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('h','u')
	plt.title('')
	plt.xlim(0,2500)
	plt.ylim(0,maxv)
	plt.xlabel("Hif1 (molecules)",fontsize=20)
	plt.ylabel("u200 (molecules)",fontsize=20)
	plt.legend(ll,labs,fontsize=20,frameon=False,bbox_to_anchor=(1.04,1))
	plt.savefig("EM_MR_hu.png",bbox_inches='tight')
	#plt.show()
	plt.close()

if myrank==3:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('h','ms')
	plt.title('')
	plt.xlim(0,2500)
	plt.ylim(0,maxv)
	plt.xlabel("Hif1 (molecules)",fontsize=20)
	plt.ylabel("Snail mRNA (molecules)",fontsize=20)
	plt.legend(ll,labs,fontsize=20,frameon=False,bbox_to_anchor=(1.04,1))
	plt.savefig("EM_MR_hms.png",bbox_inches='tight')
	#plt.show()
	plt.close()

if myrank==4:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('h','u3')
	plt.title('')
	plt.xlim(0,2500)
	plt.ylim(0,maxv)
	plt.xlabel("Hif1 (molecules)",fontsize=20)
	plt.ylabel("u34 (molecules)",fontsize=20)
	plt.legend(ll,labs,fontsize=20,frameon=False,bbox_to_anchor=(1.04,1))
	plt.savefig("EM_MR_hu3.png",bbox_inches='tight')
	#plt.show()
	plt.close()
