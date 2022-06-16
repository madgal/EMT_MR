import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
import matplotlib
from pylab import *
import copy
import random as rand


def nullclines(freepar,DSarg,fixed_points,new_domains,maxnum=100,maxstep=1e+2,minstep=1e-0,step=1e+2):
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
    PCargs_X.MaxStepSize = maxstep
    PCargs_X.MinStepSize = minstep
    PCargs_X.StepSize = step
    PCargs_X.StopAtPoints = ['B'] ##  stops searching once this point is hit
    PCargs_X.LocBifPoints = 'all' ## 
    PCargs_X.SaveEigen = True

    PyCont_X = PyDSTool.ContClass(ode_X) # Set up continuation class
    PyCont_X.newCurve(PCargs_X)
    PyCont_X['nullcline'].forward()
    PyCont_X['nullcline'].backward()
    return PyCont_X['nullcline'].sol



def getCont(DSarg,fixed_points,maxnum,step,fp,fp_domain):
        DSargs = copy.deepcopy(DSarg)
        freepar=fp
        DSargs.pdomain={fp:fp_domain}
        odex = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        odex.set(ics=fixed_points[0])
        PCargs = PyDSTool.args(name='EMT-MR',type='EP-C')
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
        PyCont['EMT-MR'].forward()
        PyCont['EMT-MR'].backward()

        return PyCont

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
