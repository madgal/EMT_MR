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

if nprocs <3:
	print("Need  at least 2")
	exit()

from parameters_metabolism import *
def getCont(DSarg,fixed_points,maxnum,step,fp='I',fpdom=[0,2500],name="EMTest"):
        DSargs = copy.deepcopy(DSarg)
        del DSargs.varspecs[fp]
        DSargs.pars[fp] = fpdom[0]
        del DSargs.xdomain[fp]
        DSargs.pdomain={fp:fpdom}
        del DSargs.ics[fp]
        odex = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        #odex.set(ics=fixed_points[0])
        PCargs = PyDSTool.args(name=name,type='EP-C')
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


def reduce_fp(fixed_points):

	final = []
	for el in fixed_points:
		if len(final)>0:
			addI=True
			for item in final:
				if np.abs(el['h']-item['h'])<0.001:
					if np.abs(el['A']-item['A'])<0.001:
						addI=False
			if addI:	
				final+=[el]
				
		else:
			final+=[el]

	return tuple(final)


DSargsF = PyDSTool.args(name='emtFull',checklevel=2)
DSargsF.fnspecs = {'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
		  'CompRmt':(['yyX','gnX','hX','h0rmtX','nhmtX','AX','A0rmtX','namtX'],'yyX*(gnX+(AX/A0rmtX)**namtX)/(1.+(hX/h0rmtX)**nhmtX+(AX/A0rmtX)**namtX)'),
		  'CompRn':(['g0X','hX','h0rnX','nhnX','ghnX','AX','A0rnX','ganX','nanxX'],'(g0X+ghnX*(hX/h0rnX)**nhnX+ganX*(AX/A0rnX)**nanxX)/(1.+(hX/h0rnX)**nhnX+(AX/A0rnX)**nanxX)'),   
		 'R':(['R1','R2'],'R1+R2')}

DSargsF.pars = { 'ga':ga,'ka':ka,'A0aa':A0aa,'A0ah':A0ah,'A0aR':A0aR,'A0rn':A0rn,'A0rm':A0rm,'lambdaaa':lambdaaa,'lambdaah':lambdaah,'lambdaar':lambdaar,'y':y,
		'g2':g2,'naa':naa,'nah':nah,'nar':nar,'narm':narm,'narn':narn,'gh':gh,'kh':kh,'h0hh':h0hh,'h0ha':h0ha,'h0hrm':h0hrm,
		'h0hrn':h0hrn,'lambdahh':lambdahh,'lambdaha':lambdaha,'g1':g1,'nhh':nhh,'nha':nha,'nhrm':nhrm,'nhrn':nhrn,'grm':grm,'grn':grn,'gn':gn,
		'krm':krm,'krn':krn,'R0ra':R0ra,'R0rh':R0rh, 'lambdara':lambdara,'lambdarh':lambdarh,'nra':nra,'nrh':nrh,
		'y':8,'kh':0.25,'u0uh':10000.,'nuh':1,'lamdauh':0.8,'u0u3R':10000,'nkrn':4,'lamdau3mR':2.,'u3':1000,'u':1000
		}
	
DSargsF.varspecs = { 'A':'ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A',
		   'Rmt':'grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt',
		   'Rnox':'grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox',
		   'h':'gh*H(A,A0ah,nah,lambdaah)*H(u,1,1,1.)-kh*h*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)'
		}
			
DSargsF.ics = {'A':300,'h':300,'Rmt':300,'Rnox':300}
DSargsF.xdomain = {'A':[0,1000],'h':[0,1000],'Rmt':[0,1000],'Rnox':[0,1000]}
DSargsF.tdomain = [0,5000]
DSargsF.algparams = {'init_step':0.1}

##############
##############
##############
def plotCC(x,y,fp='u3'):
	cl= ['b','m','k','g','r','y','orange']
	el= [3.,4.,5.,6.,7.]
	ll,labs=[],[]
	minv,maxv=1000000,0
	for i in range(len(el)):
		DSargsF.pars['krn']=el[i]
		odeF = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
		Ffixed_points = pp.find_fixedpoints(odeF,n=5,maxsearch=1e+4,eps=1e-10)
		fixed_pointsF = reduce_fp(Ffixed_points)
		PyContF=getCont(DSargsF,fixed_pointsF,2000,1e+0,fp=fp,fpdom=[1,2500])
		PyContF.display([x,y],stability=True,color=cl[i],linewidth=2,label=el)
		PyContF.plot.toggleLabels("off")

		sol = PyContF['EMTest'].sol
		if np.min(sol[y])<minv:
			minv = np.min(sol[y])
		if np.max(sol[y])>maxv:
			maxv = np.max(sol[y])

		labs +=[el[i]]
		ll +=[Line2D([0],[0],color=cl[i],lw=2)]
	return [ll,labs,minv,maxv]
#####
#####
if myrank==0:
	try:
		new_domains =  {'A':[0,1000],'h':[0,1000],'Rmt':[0,1000],'Rnox':[0,1000]}
		maxnum=2000
		odeF = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
		Ffixed_points = pp.find_fixedpoints(odeF,n=5,maxsearch=1e+4,eps=1e-10)
		fixed_pointsF = reduce_fp(Ffixed_points)
		DSargsF.pars['krn']=krn

		solh= nullclines('h',DSargsF,fixed_pointsF,new_domains,maxnum)
		sola= nullclines('A',DSargsF,fixed_pointsF,new_domains,maxnum)
		fig = plt.figure()
		plotNullcline(solh,sola,fixed_pointsF,'MR')
		plt.title('')
		plt.xlabel("AMPK",fontsize=20)
		plt.ylabel("Hif1 ",fontsize=20)
		plt.savefig("EMnull_krn_AH.png",bbox_inches='tight')
		#plt.show()
		plt.close()
	except:
		print("E0")

if myrank==1:
	try:
		fig = plt.figure()
		[ll,labs,maxv,minv]=plotCC('A','h','A')
		plt.title('')
		#plt.xlim(0,1500)
		plt.ylim(ymin=0)#minv,maxv)
		plt.xlabel("ampk ",fontsize=20)
		plt.ylabel("hif1 ",fontsize=20)
		plt.legend(ll,labs)
		plt.savefig("EM_krn_AH.png",bbox_inches='tight')
		#plt.show()
		plt.close()
	except:
		print("E1")

if myrank==2:
	try:
		fig = plt.figure()
		[ll,labs,maxv,minv]=plotCC('A','Rnox','A')
		plt.title('')
		#plt.xlim(0,1500)
		plt.ylim(ymin=0)#minv,maxv)
		plt.xlabel("ampk ",fontsize=20)
		plt.ylabel("hif1 ",fontsize=20)
		plt.legend(ll,labs)
		plt.savefig("EM_krn_AN.png",bbox_inches='tight')
		#plt.show()
		plt.close()
	except:
		print("E2")
if myrank==3:
	try:
		fig = plt.figure()
		[ll,labs,maxv,minv]=plotCC('A','Rmt','A')
		plt.title('')
		#plt.xlim(0,1500)
		plt.ylim(ymin=0)#minv,maxv)
		plt.xlabel("ampk ",fontsize=20)
		plt.ylabel("hif1 ",fontsize=20)
		plt.legend(ll,labs)
		plt.savefig("EM_krn_AM.png",bbox_inches='tight')
		#plt.show()
		plt.close()
	except:
		print("E3")
