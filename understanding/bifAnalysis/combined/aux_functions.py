import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
import matplotlib
from pylab import *
import copy
import math
matplotlib.rcParams.update({'font.size': 30})
from matplotlib.lines import Line2D

def getCont(DSarg,fixed_points,maxnum,step,fp='I',fpdom=[0,2500],name="EMTest"):
        DSargs = copy.deepcopy(DSarg)
        del DSargs.varspecs[fp]
        DSargs.pars[fp] = fpdom[0]
        del DSargs.ics[fp]
        del DSargs.xdomain[fp]
        DSargs.pdomain={fp:fpdom}
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


def plotCC(x,y,param,DSargsF,name,fpdom,cl,el):
	ll = []
	labs=[]
	minv,maxv=1000000,0
	for i in range(len(el)):
		DSargsF.pars[param]=el[i]
		odeF = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
		Ffixed_points = pp.find_fixedpoints(odeF,n=10,maxsearch=1e+6,eps=1e-10)
		fixed_pointsF = reduce_fp(Ffixed_points)
		DSargsF.ics = fixed_pointsF
		PyContF=getCont(DSargsF,fixed_pointsF,2000,1e+0,x,fpdom,name)
		PyContF.display([x,y],stability=True,color=cl[i],linewidth=2)
		PyContF.plot.toggleLabels("off")

		solV = PyContF[name].sol
		if np.min(solV[y])<minv:
			minv = np.min(solV[y])
		if np.max(solV[y])>maxv:
			maxv = np.max(solV[y])

		labs+=[el[i]]
		ll +=[Line2D([0],[0],color=cl[i],lw=2)]

	return [ll,labs,maxv,minv]

