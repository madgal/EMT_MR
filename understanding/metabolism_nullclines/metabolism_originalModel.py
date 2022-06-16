import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
from pylab import *
import copy
import math
import matplotlib
matplotlib.rcParams.update({'font.size': 30})
import random as rand
import json


from parameters_metabolism_og import *

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
                                if np.abs(el['h']-item['h'])<0.001:
                                        if np.abs(el['A']-item['A'])<0.001:
                                                if np.abs(el['Rnox']-item['Rnox'])<0.001:
                                                	if np.abs(el['Rmt']-item['Rmt'])<0.001:
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
    PCargs_X.MaxNumPoints = 100000## max number for each call of forward and backward
    PCargs_X.MaxStepSize = 1e+2
    PCargs_X.MinStepSize = 1e-1
    PCargs_X.StepSize = 1e-0
    PCargs_X.StopAtPoints = ['B'] ##  stops searching once this point is hit
    PCargs_X.LocBifPoints = 'all' ##
    PCargs_X.SaveEigen = True

    PyCont_X = PyDSTool.ContClass(ode_X) # Set up continuation class
    PyCont_X.newCurve(PCargs_X)
    PyCont_X['nullcline'].forward()
    PyCont_X['nullcline'].backward()
    return PyCont_X['nullcline'].sol

##  units
## [dA/dt] = nM/hr
## [dH/dt] = nM/hr
## [dRmt/dt] = uM/min
## [dRnox/dt] = uM/min

## metabolism for cancer cells
DSargs = PyDSTool.args(name='metabolism-cc',checklevel=2)

DSargs.fnspecs = {'H':(['X','X0','lamdaX','nX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
		  'CompRmt':(['yyX','gnX','hX','h0rmtX','nhmtX','AX','A0rmtX','namtX'],'yyX*(gnX+(AX/A0rmtX)**namtX)/(1.+(hX/h0rmtX)**nhmtX+(AX/A0rmtX)**namtX)'),
		  'CompRn':(['g0X','hX','h0rnX','nhnX','ghnX','AX','A0rnX','ganX','nanxX'],'(g0X+ghnX*(hX/h0rnX)**nhnX+ganX*(AX/A0rnX)**nanxX)/(1.+(hX/h0rnX)**nhnX+(AX/A0rnX)**nanxX)'),   
		 'R':(['R1','R2'],'R1+R2')}

DSargs.pars = {'ga'         :ga, 'ka'        :ka, 'A0aa'      :A0aa, 'A0ah'      :A0ah, 'A0aR'      :A0aR, 'A0rn'      :A0rn, 'A0rm'      :A0rm, 'lambdaaa'  :lambdaaa,
		'lambdaah'  :lambdaah, 'lambdaar'  :lambdaar, 'y'         :y, 'g2'        :g2, 'naa'       :naa, 'nah'       :nah, 'nar'       :nar, 'narm'      :narm,
		'narn'      :narn, 'gh'        :gh, 'kh'        :kh, 'h0hh'      :h0hh, 'h0ha'      :h0ha, 'h0hrm'     :h0hrm, 'h0hrn'     :h0hrn, 'lambdahh'  :lambdahh,
		'lambdaha'  :lambdaha, 'g1'        :g1, 'nhh'       :nhh, 'nha'       :nha, 'nhrm'      :nhrm, 'nhrn'      :nhrn, 'grm'       :grm, 'grn'       :grn,
		'gn'        :gn, 'krm'       :krm, 'krn'       :krn, 'R0ra'      :R0ra, 'R0rh'      :R0rh, 'lambdara'  :lambdara, 'lambdarh'  :lambdarh, 'nra'       :nra,
		'nrh'       :nrh}

DSargs.varspecs = {'A':'ga*H(R(Rmt,Rnox),R0ra,lambdara,nra)*H(h,h0ha,lambdaha,nha)*H(A,A0aa,lambdaaa,naa)-ka*A',
		   'Rmt':'grm*H(A,A0aR,lambdaar,nar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt',
		   'Rnox':'grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox',
		   'h':'gh*H(A,A0ah,lambdaah,nah)-kh*h*H(h,h0hh,lambdahh,nhh)*H(R(Rmt,Rnox),R0rh,lambdarh,nrh)'}

## cancer
DSargs.pars['y']= 8
DSargs.pars['kh']= 0.25
DSargs.ics={'A':300,'h':300,'Rmt':300,'Rnox':300}
DSargs.xdomain={'A':[0,1000],'h':[0,1000],'Rmt':[0,1000],'Rnox':[0,1000]}
DSargs.tdomain = [0,5000]
DSargs.algparams = {'init_step':0.1}

ode=PyDSTool.Generator.Vode_ODEsystem(DSargs)

Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+4,eps=1e-10)
fixed_points = reduce_fp(Ffixed_points)
stabilities = stability(fixed_points,DSargs)
new_domains={'A':[0,1000],'h':[0,1000],'Rmt':[0,1000],'Rnox':[0,1000]}


solC_h = nullclines('h',DSargs,fixed_points,new_domains)
solC_a = nullclines('A',DSargs,fixed_points,new_domains)


results={'fixed_points':fixed_points,'nullA_A':list(solC_a['A']),'nullA_H':list(solC_a['h']),'nullH_A':list(solC_h['A']),'nullH_H':list(solC_h['h']),'stability':stabilities}

json_data = json.dumps(results,indent=4)
with open("data_metabolism_nullcline_og.json",'w') as outfile:
	outfile.write(json_data)
