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

def getCont(DSarg,fixed_points,maxnum,step,fp,fp_domain):
        DSargs = copy.deepcopy(DSarg)
        if fp in DSargs.varspecs.keys():
            del DSargs.varspecs[fp]
        if fp in DSargs.ics.keys():
      	    del DSargs.ics[fp]
        if fp in DSargs.xdomain.keys():
            del DSargs.xdomain[fp]
        #print fixed_points[0][fp],fp
        if fp in fixed_points[0].keys():
            DSargs.pars[fp]=fixed_points[0][fp]
        else:
            DSargs.pars[fp]=1000
        #print DSargs.pars[fp]
        freepar=fp
        DSargs.pdomain={fp:fp_domain}
        up_fixed_points = copy.deepcopy(fixed_points)
        for i in range(len(up_fixed_points)):
            if fp in up_fixed_points[i].keys():
                del up_fixed_points[i][fp]
        odex = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        odex.set(ics=up_fixed_points[0])
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



DSargsF = PyDSTool.args(name='emtFull',checklevel=2)
DSargsF.fnspecs = {'Z':(['u','mz'],'gz*mz*L(u,u0,nu)/kz'),
		  'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
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

DSargsF.varspecs = {'u':'gu*H(Z(u,mz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)-mz*Yu(u,u0,nu)-mh*Yuh(0,u0,nuh)-ku*u',
		   'mz':'gmz*H(Z(u,mz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)*H(A,A0m,nAm,lamdaAm)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)*H(h,h0ms,nhms,lamdahms)*H(A,A0ms,nAms,lamdaAms)-ms*Ym(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(u,mz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                   'S':'gs*ms*L3(u3,u30,nu3)-ks*S'
		}
			
DSargsF.pars ={'gz':100.,'gu':2100.,'grn':40.,'grm':50.,'gms':90.,'gu3':1350.,'gs':100.,'ga':30.,'gn':0.2,'g1':5.,'g2':0.2,'gmh':10.,'gh':1.5,'gmz':11.,'y':8.,
		'ks':0.125,'ka':0.2,'kz':0.1,'ku':0.05,'kmz':0.5,'kms':0.5,'ku3':0.05,'krm':5.,'krn':5.,'kmh':.143,'kh':1.75,
		'Z0u':220000.,'I0m':50000.,'S0ms':200000.,'u30':10000.,'S0u3':300000.,'h0ha':250.,'A0aa':350.,'h0hrm':200.,'h0hrn':250.,'A0rn':150.,'A0rm':150.,
		'A0aR':350.,'R0ra':100.,'Z0u3':600000.,'Z0m':25000,'S0u':180000.,'S0m':180000.,'u0':10000.,'h0u':200.,'A0u':300.,'A0m':300.,'A0ms':300.,
		'h0ms':200.,'nu30rn':10000.,'nu30rm':10000.,'A0ah':250.,'R0rh':300.,'h0hh':80.,
		'lamdahu':1.,'lamdaAu':1.,'lamda3n':1.,'lamda3m':1.,'lamdaAms':1.,'lamdahms':1.,'lamdaAm':1.,
		'lambdarh':0.2,'lambdahh':0.1,'lambdaah':0.1,'lambdaar':0.25,'lambdara':8.,'lambdaha':0.1,'lambdaaa':0.2,
		'lamdazu':0.1,'lamdaSu':0.1,'lamdaIm':10.,'lamdaSms':0.1,'lamdazu3':0.2,'lamdaSu3':0.1,'lamdaZm':7.5,'lamdaSm':10.,
		'naa':2,'nu3':2,'nzu3':2,'nsms':1,'nIm':2,'nsu3':1,'nra':4,'nha':1,'nar':2,'narm':4,'nhrm':2,'nhrn':2,
		'nzu':3,'nsu':2,'nzm':2,'nsm':2,'nu':6,'narn':2,'nah':1,'nrh':4,'nhh':4,'nuh':2,'nhu':1,'nAu':1,'nAm':2,'nAms':2,'nhms':2,'n3m':3,'n3n':2,
		'I':50000.,
		'ymi0':0.,'ymi1':0.04,'ymi2':0.2,'ymi3':1.,'ymi4':1.,'ymi5':1.,'ymi6':1.,
		'yui0':0.,'yui1':0.005,'yui2':0.05,'yui3':0.5,'yui4':0.5,'yui5':0.5,'yui6':0.5,
		'li0':1.,'li1':0.6,'li2':0.3,'li3':0.1,'li4':0.05,'li5':0.05,'li6':0.05,
		'ymih0':0.,'ymih1':0.04,'ymih2':0.2,
		'yuih0':0.,'yuih1':0.005,'yuih2':0.05,
		'lih0':1.,'lih1':0.6,'lih2':0.3 ,
		'h':0.,'A':0. ,'mh':0.}

##############
##############
#####
#####
DSargsF.ics = {'mz':1000,'S':200000,'ms':2000,'u3':4000,'u':1000}
DSargsF.xdomain = {'u':[0,100000],'mz':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000]}
new_domains = {'u':[0,100000],'mz':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000]}
DSargsF.tdomain = [0,5000]
DSargsF.algparams = {'init_step':0.4}

lambdaVals = { 'lamdahu':[1.,0.85,0.8,0.7,0.6,0.2,0.0],
		'lamdahms':[1.,1.1,1.19,1.24,1.3,9.,15.],
		'lamdaAms':[1.,0.96,0.9,0.73,0.6,0.0],
		'lamdaAm':[1.,0.95,0.8,0.72,0.5,0.0],
		'lamdaAu':[1.,1.05,1.4,1.8,3,4]}


if nprocs!=5:
	print("Nprocs==5")
	exit()

#lambdaVals = { 'lamdahu':[0.85,0.8,0.7,0.6,0.2,0.0],
#		'lamdaAu':[4.,3.,1.8,1.5,1.1,1.05],
#		'lamdaAms':[1., 0.96,0.9,0.7, 0.6,0.53],
#		'lamdahms':[1.1,1.19,1.24,1.3,9.,15.],
#		'lamdaAm':[0.95,0.9,0.72,0.5,0.36,0.0]}
#if nprocs!=5:
#	print("Nprocs==5")
#	exit()


######{'A': 436.47213160474945, 'h': 26.69259410719074}
######{'A': 327.50827624867554,'h': 292.59779290727334}
######{'A': 68.67768314578954, 'h': 480.4124463778733}

val_w={'A':68.7,'h':480.4}
val_wo={'A':327.5,'h':292.6}
val_o={'A':436.5,'h':26.7}

print(list(lambdaVals.keys()))
lamda = list(lambdaVals.keys())[myrank]

for lamI in range(len(lambdaVals[lamda])):	

	DSargsF.pars[lamda] = lambdaVals[lamda][lamI]

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


	fileO = open("data_bif/EMW_fp_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
	fileO.write('S,ms,mz,u,u3\n')
	for i in range(len(fixed_points)):
		fileO.write("%s,%s,%s,%s,%s\n" %(  fixed_points[i]['S'],fixed_points[i]['ms'], fixed_points[i]['mz'], fixed_points[i]['u'], fixed_points[i]['u3']))

	fileO.close()


	if 'h' in lamda:
		free_key='h'
	else:
		free_key='A'
	PyContF=getCont(DSargsF,fixed_points,8000,1e-1,free_key,[0,2000])
	sol_cont = PyContF['EMT-MR'].sol
	lp_cont=True
	i=1
	fileO = open("data_bif/EMW_LP_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
	fileO.write('LPN,S,ms,mz,u,u3,'+free_key+'\n')
	while lp_cont and i<1000:
		try:
			mystr = 'LP'+str(i)
			tmpLP = PyContF['EMT-MR'].getSpecialPoint(mystr).labels['LP']['data'].X
			fileO.write("%s,%s,%s,%s,%s,%s,%s\n" %(mystr, tmpLP['S'], tmpLP['ms'], tmpLP['mz'], tmpLP['u'], tmpLP['u3'], tmpLP[free_key]))
		except:
			lp_cont=False
		i+=1
	fileO.close()


	fileO = open("data_bif/EMW_bif_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
	fileO.write('stability,S,ms,mz,u,u3,'+free_key+'\n')
	for i in range(len(sol_cont)):
		fileO.write("%s,%s,%s,%s,%s,%s,%s\n" %(sol_cont[i].labels['EP']['stab'], sol_cont['S'][i], sol_cont['ms'][i], sol_cont['mz'][i], sol_cont['u'][i], sol_cont['u3'][i], sol_cont[free_key][i]))

	fileO.close()



	### nullclines
	try:
		DSargsF.pars[free_key]=val_w[free_key]
		solC_u = nullclines('u',DSargsF,fixed_points,new_domains,1000)
		fileO = open("data/nullU_W_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
		fileO.write('S,ms,mz,u,u3,'+free_key+'\n')
		for i in range(len(solC_u['u'])):
			fileO.write("%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i], solC_u[free_key][i]))
		fileO.close()
	except:
		print('Issue with nullcline uw')


	try:
		DSargsF.pars[free_key]=val_w[free_key]
		solC_mz = nullclines('mz',DSargsF,fixed_points,new_domains,1000)
		fileO = open("data/nullmz_W_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
		fileO.write('S,ms,mz,u,u3,'+free_key+'\n')
		for i in range(len(solC_mz['u'])):
			fileO.write("%s,%s,%s,%s,%s,%s\n" %( solC_mz['S'][i], solC_mz['ms'][i], solC_mz['mz'][i], solC_mz['u'][i], solC_mz['u3'][i], solC_mz[free_key][i]))
		fileO.close()
	except:
		print('Issue with nullcline zw')


	try:
		DSargsF.pars[free_key]=val_wo[free_key]
		solC_u = nullclines('u',DSargsF,fixed_points,new_domains,1000)
		fileO = open("data/nullU_WO_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
		fileO.write('S,ms,mz,u,u3,'+free_key+'\n')
		for i in range(len(solC_u['u'])):
			fileO.write("%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i], solC_u[free_key][i]))
		fileO.close()
	except:
		print('Issue with nullcline uwo')


	try:
		DSargsF.pars[free_key]=val_wo[free_key]
		solC_mz = nullclines('mz',DSargsF,fixed_points,new_domains,1000)
		fileO = open("data/nullmz_WO_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
		fileO.write('S,ms,mz,u,u3,'+free_key+'\n')
		for i in range(len(solC_mz['u'])):
			fileO.write("%s,%s,%s,%s,%s,%s\n" %( solC_mz['S'][i], solC_mz['ms'][i], solC_mz['mz'][i], solC_mz['u'][i], solC_mz['u3'][i], solC_mz[free_key][i]))
		fileO.close()
	except:
		print('Issue with nullcline zwo')


	try:
		DSargsF.pars[free_key]=val_o[free_key]
		solC_u = nullclines('u',DSargsF,fixed_points,new_domains,1000)
		fileO = open("data/nullU_O_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
		fileO.write('S,ms,mz,u,u3,'+free_key+'\n')
		for i in range(len(solC_u['u'])):
			fileO.write("%s,%s,%s,%s,%s,%s\n" %( solC_u['S'][i], solC_u['ms'][i], solC_u['mz'][i], solC_u['u'][i], solC_u['u3'][i], solC_u[free_key][i]))
		fileO.close()
	except:
		print('Issue with nullcline uo')


	try:
		DSargsF.pars[free_key]=val_o[free_key]
		solC_mz = nullclines('mz',DSargsF,fixed_points,new_domains,1000)
		fileO = open("data/nullmz_O_"+lamda+"_"+str(lambdaVals[lamda][lamI])+".txt",'w')
		fileO.write('S,ms,mz,u,u3,'+free_key+'\n')
		for i in range(len(solC_mz['u'])):
			fileO.write("%s,%s,%s,%s,%s,%s\n" %( solC_mz['S'][i], solC_mz['ms'][i], solC_mz['mz'][i], solC_mz['u'][i], solC_mz['u3'][i], solC_mz[free_key][i]))
		fileO.close()
	except:
		print('Issue with nullcline zo')


DSargsF.pars[lamda] = 1.




