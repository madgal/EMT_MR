import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
import matplotlib
from pylab import *
import copy
import math
matplotlib.rcParams.update({'font.size': 30})


from parameters_emt_B import *


def plotContPartial(DSarg,fixed_points,x,y,title,maxnum,step,col='k'):
        DSargs = copy.deepcopy(DSarg)
        freepar=x
        DSargs.pdomain={x:[0,100000]}
        odex = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        odex.set(ics=fixed_points[0])
        PCargs = PyDSTool.args(name='EMTest',type='EP-C')
        PCargs.freepars =[x]
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
        PyCont.display([x,y],stability=True,color=col)
        PyCont.plot.toggleLabels("off")
        plt.title('')
        plt.xlabel(x+" (molecules)",fontsize=20)
        plt.ylabel(y+"(molecules)",fontsize=20)
        #plt.savefig(title+".png",bbox_inches='tight')
        #plt.show()
        #plt.close()

def plotNullcline(sol_u,sol_mz,fixed_points,title):
        fig = plt.figure(figsize=(12,9))
        fig.subplots_adjust(hspace=0.4,wspace=0.4)
        plt.subplot(1,1,1)
        plt.plot(sol_u['u']/1000.,sol_u['mz'],color='g',label='dU/dt=0 and dZ/dt=0',linewidth=3)
        plt.plot(sol_mz['u']/1000.,sol_mz['mz'],color='b',label='dmz/dt=0 and dZ/dt=0',linewidth=3)

        colors=['k','k','k','w','w']
        for el in range(len(fixed_points)):
                plt.plot(fixed_points[el]['u']/1000.,fixed_points[el]['mz'],'o',color=colors[el],markersize=15,markeredgecolor='k',markeredgewidth=1.5)
        plt.xlabel('$\mu_{200}$ (K molecules)')
        plt.ylabel('Zeb mRNA (molecules)')
        plt.legend(frameon=False)
        plt.title("EMT")
        plt.savefig(title+".png",bbox_inches='tight')
        #plt.show()
        plt.close()

def plotCont(DSarg,fixed_points,title,maxnum,step):
        DSargs = copy.deepcopy(DSarg)
        freepar='I'
        DSargs.pdomain={'I':[0,100000]}
        odex = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        odex.set(ics=fixed_points[0])
        PCargs = PyDSTool.args(name='EMTest',type='EP-C')
        PCargs.freepars =['I']
        PCargs.MaxNumPoints = maxnum## max number for each call of forward and backward
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
        PyCont.display(['S','mz'],stability=True)
        PyCont.plot.toggleLabels("off")
        plt.savefig(title+".png",bbox_inches='tight')
        #plt.show()
        plt.close()


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


def nullclines(freepar,DSarg,fixed_points,new_domains,maxnum):
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
    PCargs_X.MaxNumPoints = maxnum## max number for each call of forward and backward
    PCargs_X.MaxStepSize = 1e+2
    PCargs_X.MinStepSize = 1e-0
    PCargs_X.StepSize = 1e+2
    PCargs_X.StopAtPoints = ['B'] ##  stops searching once this point is hit
    PCargs_X.LocBifPoints = 'all' ## 
    PCargs_X.SaveEigen = True

    PyCont_X = PyDSTool.ContClass(ode_X) # Set up continuation class
    PyCont_X.newCurve(PCargs_X)
    PyCont_X['nullcline'].forward()
    PyCont_X['nullcline'].backward()
    return PyCont_X['nullcline'].sol


DSargs = PyDSTool.args(name='emt',checklevel=2)
DSargs.fnspecs = {'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
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

DSargs.pars = {	'gz':gz, 'u0':u0,'nu':nu,'kz':kz,
                'li0':li0, 'li1':li1, 'li2':li2, 'li3':li3, 'li4':li4, 'li5':li5, 'li6':li6,
                'ymi0':ymi0, 'ymi1':ymi1, 'ymi2':ymi2, 'ymi3':ymi3, 'ymi4':ymi4, 'ymi5':ymi5, 'ymi6':ymi6,
                'yui0':yui0, 'yui1':yui1, 'yui2':yui2, 'yui3':yui3, 'yui4':yui4, 'yui5':yui5, 'yui6':yui6,
		'gu':gu, 'Z0u':Z0u,'nzu':nzu, 'lamdazu':lamdazu,'S0u':S0u,'nsu':nsu,'lamdaSu':lamdaSu,'u0':u0,'nu':nu,'ku':ku,
		'gmz':gmz, 'Z0m':Z0m,'nzm':nzm, 'lamdaZm':lamdaZm,'S0m':S0m,'nsm':nsm,'lamdaSm':lamdaSm,'u0':u0,'nu':nu,'kmz':kmz,
		'gms':gms, 'I0m':I0m, 'nIm':nIm, 'lamdaIm':lamdaIm, 'S0ms':S0ms, 'nsms':nsms, 'lamdaSms':lamdaSms, 'u30':u30,
                'nu3':nu3, 'kms':kms, 'gu3':gu3, 'Z0u3':Z0u3, 'nzu3':nzu3, 'lamdazu3':lamdazu3, 'S0u3':S0u3, 'nsu3':nsu3,
		'lamdaSu3':lamdaSu3, 'ku3':ku3, 'gs':gs, 'ks':ks}
	
DSargs.varspecs = { 'u':'gu*H(Z(mz,u,u0,nu,kz,gz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)-mz*Yu(u,u0,nu)-ku*u',
		   'mz':'gmz*H(Z(mz,u,u0,nu,kz,gz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(I,I0m,nIm,lamdaIm)*H(S,S0ms,nsms,lamdaSms)-ms*Ym3(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(mz,u,u0,nu,kz,gz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                    'S':'gs*ms*L3(u3,u30,nu3)-ks*S'

		}
			
DSargs.ics = {'u':1000,'mz':1000,'S':200000,'ms':2000,'u3':4000}
DSargs.xdomain = {'u':[0,100000],'mz':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000]}
DSargs.tdomain = [0,5000]
DSargs.algparams = {'init_step':0.1}
DSargs.pars['I']=50000

ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)

Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)

fixed_points = reduce_fp(Ffixed_points)
print( fixed_points)

#######3nullcines

new_domains = {'u':[0,20000],'mz':[0,1500],'S':[0,400000],'ms':[0,20000],'u3':[0,20000]}

#solmz = nullclines('mz',DSargs,fixed_points,new_domains,1000)
#solu = nullclines('u',DSargs,fixed_points,new_domains,1000)

#plotNullcline(solu,solmz,fixed_points,'EMT_full')
#plotCont(DSargs,fixed_points,'EMT_fullCont',2000,1e+0)



'''
###################################
count=0
for i in [8.,9.,10.,11.,12.]:
    try:
        DSargs.pars['lamdaIm']=i
        DSargs.xdomain = {'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}
        ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
        fixed_points = reduce_fp(Ffixed_points)
        new_domains ={'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}

        solmz = nullclines('mz',DSargs,fixed_points,new_domains,1000)
        solu = nullclines('u',DSargs,fixed_points,new_domains,1000)
        plotNullcline(solu,solmz,fixed_points,'EMTf_lamdaI_'+str(count))
        count+=1
    except:
        print( "Error in lambdaI")
###################################
DSargs.pars['lamdaIm']=10.
for i in range(30000,70000,1000):
    try:
        DSargs.pars['I']=i
        DSargs.xdomain = {'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}
        ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
        fixed_points = reduce_fp(Ffixed_points)
        new_domains ={'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}

        solmz = nullclines('mz',DSargs,fixed_points,new_domains,1000)
        solu = nullclines('u',DSargs,fixed_points,new_domains,1000)
        plotNullcline(solu,solmz,fixed_points,'EMTf_I_'+str(i))
    except:
        print( "Error in I")
###################################
DSargs.pars['I']=50000.
count=0
for i in [0,0.05,0.1,0.2,0.3]:
    try:
        DSargs.pars['lamdaSu']=i
        DSargs.xdomain = {'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}
        ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
        fixed_points = reduce_fp(Ffixed_points)
        new_domains ={'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}
    
        solmz = nullclines('mz',DSargs,fixed_points,new_domains,1000)
        solu = nullclines('u',DSargs,fixed_points,new_domains,1000)
        plotNullcline(solu,solmz,fixed_points,'EMTf_lamdaSu_'+str(i))
        count+=0
    except:
        print( "Error in lambdaSu")
###################################
DSargs.pars['lamdaSu']=lamdaSu
count=0
for i in [8.,9.,10.,11.,12.]:
    try:
        DSargs.pars['lamdaSm']=i
        DSargs.xdomain = {'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}
        ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
        fixed_points = reduce_fp(Ffixed_points)
        new_domains ={'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}

        solmz = nullclines('mz',DSargs,fixed_points,new_domains,1000)
        solu = nullclines('u',DSargs,fixed_points,new_domains,1000)
        plotNullcline(solu,solmz,fixed_points,'EMTf_lamdaSu_'+str(i))
        count+=0
    except:
        print( "Error in lambdaSm")
'''

##############
##############
##############
DSargs.xdomain = {'u':[0,100000],'mz':[0,5000],'S':[0,400000],'ms':[0,50000],'u3':[0,50000]}

el=[0,1,2,3,4]
col=['r','m','k','g','b']
fig = plt.figure()
for i in range(len(el)):
	try:
                DSargs.pars['lamdaIm']=el[i]
                ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
                Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
                fixed_points = reduce_fp(Ffixed_points)
                plotContPartial(DSargs,fixed_points,'I','u','EMT_fullContU_lamdaIm',2000,1e+0,col[i])
                #plotContPartial(DSargs,fixed_points,'I','mz','EMT_fullCont_lamdaIm',2000,1e+0,col[i])
	except:
		print("Error in lamdaIm"+str(i))
plt.savefig("EMT_fullContU_lamdaIm.png",bbox_inches='tight')
plt.close()

##############
DSargs.pars['lamdaIm']=10.
el=[50000]
col=['k']
fig = plt.figure()
for i in range(len(el)):
	try:
                DSargs.pars['I']=el[i]
                ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
                Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
                fixed_points = reduce_fp(Ffixed_points)
                plotContPartial(DSargs,fixed_points,'I','u','EMT_fullContU_I',2000,1e+0,col[i])
                #plotContPartial(DSargs,fixed_points,'I','mz','EMT_fullCont_I',2000,1e+0,col[i])
	except:
		print("Error in I"+str(i))
plt.savefig("EMT_fullContU_I.png",bbox_inches='tight')
plt.close()

##############
DSargs.pars['I']=50000
el=[0,0.1,0.2,0.5]
col=['m','k','g','b']
fig = plt.figure()
for i in range(len(el)):
	try:
                DSargs.pars['lamdaSu']=el[i]
                ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
                Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
                fixed_points = reduce_fp(Ffixed_points)
                plotContPartial(DSargs,fixed_points,'I','u','EMT_fullContU_lamdaSu',2000,1e+0,col[i])
                #plotContPartial(DSargs,fixed_points,'I','mz','EMT_fullCont_lamdaSu',2000,1e+0,col[i])
	except:
		print("Error in lamdaSu"+str(i))
plt.savefig("EMT_fullContU_lamdaSu.png",bbox_inches='tight')
plt.close()

##############
DSargs.pars['lamdaSu']=lamdaSu
el=[ 8,9,10,12,15]
col=['r','m','k','g','b']
fig = plt.figure()
for i in range(len(el)):
	try:
                DSargs.pars['lamdaSm']=el[i]
                ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
                Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
                fixed_points = reduce_fp(Ffixed_points)
                plotContPartial(DSargs,fixed_points,'I','u','EMT_fullContU_lamdaSm',2000,1e+0,col[i])
                #plotContPartial(DSargs,fixed_points,'I','mz','EMT_fullCont_lamdaSm',2000,1e+0,col[i])
	except:
		print("Error in lamdaSm"+str(i))
plt.savefig("EMT_fullContU_lamdaSm.png",bbox_inches='tight')
plt.close()
