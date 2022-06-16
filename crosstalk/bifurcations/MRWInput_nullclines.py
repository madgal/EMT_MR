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


DSargsF = PyDSTool.args(name='MRFull',checklevel=2)
DSargsF.fnspecs = {'H':(['X','X0','nX','lamdaX'],'lamdaX+(1.-lamdaX)/(1.+(X/X0)**nX)'),
		  'M' : (['i','n','x','x0'],'(x/x0)**i/(1+(x/x0))**n'),
                  'indLh': (['i'],'if(i==0,lih0,if(i==1,lih1,if(i==2,lih2,lih2)))'),
                  'indYmh': (['i'],'if(i==0,ymih0,if(i==1,ymih1,if(i==2,ymih2,ymih2)))'),
                  'combination':(['k','n'],'special_gamma(n+1)/special_gamma(k+1)/special_gamma(n-k+1)'),
		  'Lh' : (['X','X0','nX'],'sum(i,0,2,lih[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'Ymh' : (['X','X0','nX'],'sum(i,0,2,ymi[i]*combination([i],nX)*M([i],nX,X,X0))'),
		  'R':(['R1','R2'],'R1+R2'),
	 	  'CompRmt':(['yyX', 'gnX', 'hX', 'h0rmtX', 'nhmtX', 'AX', 'A0rmtX',  'namtX'],'yyX*(gnX+(AX/A0rmtX)**namtX)/(1.+(hX/h0rmtX)**nhmtX+(AX/A0rmtX)**namtX)'),
		  'CompRn':(['g0X', 'hX', 'h0rnX', 'nhnX', 'ghnX', 'AX', 'A0rnX', 'ganX',  'nanxX'],'(g0X+ghnX*(hX/h0rnX)**nhnX+ganX*(AX/A0rnX)**nanxX)/(1.+(hX/h0rnX)**nhnX+(AX/A0rnX)**nanxX)')
			}

DSargsF.varspecs = { 'A':'ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A',
		   'h':'gh*mh*Lh(0.,u0,nuh)-kh*h',
		   'Rmt':'grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt*H(u3,nu30rm,n3m,lamda3m)',
		   'Rnox':'grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox*H(u3,nu30rn,n3n,lamda3n)',
 		   'mh':'gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ymh(0.,u0,nuh)',
		}
			
DSargsF.pars ={'gz':100.,'gu':2100.,'grn':40.,'grm':150.,'gms':90.,'gu3':1350.,'gs':100.,'ga':30.,'gn':0.2,'g1':5.,'g2':0.2,'gmh':10.,'gh':1.5,'gmz':11.,'y':8.,
		'ks':0.125,'ka':0.2,'kz':0.1,'ku':0.05,'kmz':0.5,'kms':0.5,'ku3':0.05,'krm':5.,'krn':5.,'kmh':.143,'kh':1.75,
		'Z0u':220000.,'I0m':50000.,'S0ms':200000.,'u30':10000.,'S0u3':300000.,'h0ha':250.,'A0aa':350.,'h0hrm':200.,'h0hrn':250.,'A0rn':150.,'A0rm':150.,
		'A0aR':350.,'R0ra':100.,'Z0u3':600000.,'Z0m':25000,'S0u':180000.,'S0m':180000.,'u0':10000.,'h0u':200.,'A0u':300.,'A0m':300.,'A0ms':300.,
		'h0ms':200.,'nu30rn':10000.,'nu30rm':10000.,'A0ah':250.,'R0rh':300.,'h0hh':80.,
		'lamdahu':1.,'lamdaAu':1.,'lamda3n':1.,'lamda3m':1.,'lamdaAms':1.,'lamdahms':1.,'lamdaAm':1.,
		'lambdarh':0.2,'lambdahh':0.1,'lambdaah':0.1,'lambdaar':0.25,'lambdara':8.,'lambdaha':0.1,'lambdaaa':0.2,
		'lamdazu':0.1,'lamdaSu':0.1,'lamdaIm':10.,'lamdaSms':0.1,'lamdazu3':0.2,'lamdaSu3':0.1,'lamdaZm':7.5,'lamdaSm':10.,
		'naa':2,'nu3':2,'nzu3':2,'nsms':1,'nIm':2,'nsu3':1,'nra':4,'nha':1,'nar':2,'narm':4,'nhrm':2,'nhrn':2,
		'nzu':3,'nsu':2,'nzm':2,'nsm':2,'nu':6,'narn':2,'nah':1,'nrh':4,'nhh':4,'nuh':2,'nhu':1,'nAu':1,'nAm':2,'nAms':2,'nhms':2,'n3m':3,'n3n':2,
		'I':50000,
		'ymi0':0.,'ymi1':0.04,'ymi2':0.2,'ymi3':1.,'ymi4':1.,'ymi5':1.,'ymi6':1.,
		'yui0':0.,'yui1':0.005,'yui2':0.05,'yui3':0.5,'yui4':0.5,'yui5':0.5,'yui6':0.5,
		'li0':1.,'li1':0.6,'li2':0.3,'li3':0.05,'li4':0.05,'li5':0.05,'li6':0.05,
		'ymih0':0.,'ymih1':0.04,'ymih2':0.2,
		'yuih0':0.,'yuih1':0.005,'yuih2':0.05,
		'lih0':1.,'lih1':0.6,'lih2':0.3 ,
		'u':0,'u3':0}

##############
##############
#####
DSargsF.ics = {'h':100,'A':100,'Rnox':100,'Rmt':100,'mh':100}
DSargsF.xdomain = {'A':[0,1000],'h':[0,1000],'Rmt':[0,1000],'Rnox':[0,1000],'mh':[0,1000]}
DSargsF.tdomain = [0,5000]
DSargsF.algparams = {'init_step':0.4}
new_domains = {'A':[0,1000],'h':[0,1000],'Rmt':[0,1000],'Rnox':[0,1000],'mh':[0,1000]}
	
	
if myrank==0:
	ode = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
	Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+6,eps=1e-10)
	#print Ffixed_points
	fixed_points = reduce_fp(Ffixed_points)
	print(fixed_points)

	fig = plt.figure()
	solC_h = nullclines('h',DSargsF,fixed_points,new_domains,1000)
	solC_a = nullclines('A',DSargsF,fixed_points,new_domains,1000)
	plt.plot(solC_h['A'],solC_h['h'],color='b',label='dA/dt=0',linewidth=3)
	plt.plot(solC_a['A'],solC_a['h'],color='g',label='dH/dt=0',linewidth=3)
	for i in range(len(fixed_points)):
		plt.plot(fixed_points[i]['A'],fixed_points[i]['h'])
	fig.savefig("nullclines_figures/mr_null.png",bbox_inches='tight')
	#plt.show()
	plt.close()


###############################

lamda_value_dict={}
ss_value_dict={}
lamda_label_dict={}
ss_label_dict={}
title_dict={}


for i in range(nprocs):
	lamda_value_dict[i]=[]
	ss_value_dict[i]=[]
	lamda_label_dict[i]=[]
	ss_label_dict[i]=[]
	title_dict[i]=[]

## HS, Hu, AS, AZ, Au
lamList = ['lamda3n','lamda3m']
lamda_value = [ [[0,1,0,0],[0,2,4,4],[0,3,0,0],[0,3,9,7],[0,4,0,0]],
		[[0,0,0,0],[0,3,0,0],[0,5,0,5],[0,7,6,0]]]
ss_label = ['u3','u3']
ss_value = [[8177.9,16908.6,16993.3],[8177.9,16908.6,16993.3]]
state_Names = [['M','EM','E'],['M','EM','E']]
crossLab = ['u3Rnox','u3Rmt']

count=0
for i in range(len(lamList)):
	for el1 in lamda_value[i]:
		for j in range(len(ss_value[i])):
			title_dict[count%nprocs] += ["nullclines_figures/nullcline_MR_"+state_Names[i][j]+"_"+crossLab[i]+"_"+str(el1[0])+"_"+str(el1[1])+str(el1[2])+".png"]
			lamda_label_dict[count%nprocs]+=[lamList[i]]
			lamda_value_dict[count%nprocs]+=[el1]
			ss_label_dict[count%nprocs]+=[ss_label[i]]
			ss_value_dict[count%nprocs]+=[ss_value[i][j]]
			count+=1

### HS 
for i in range(len(lamda_value_dict[myrank])):
	el = lamda_value_dict[myrank][i] 
	ssv = ss_value_dict[myrank][i] 
	lam= lamda_label_dict[myrank][i] 
	ssl = ss_label_dict[myrank][i] 
	title = title_dict[myrank][i]

	print(el,ssv,lam,ssl)
	print(title,'\n')

	try:
		DSargsF.pars[ssl]=ssv
		DSargsF.pars[lam]=el[0]*1.+el[1]*0.1+el[2]*0.01
		ode = PyDSTool.Generator.Vode_ODEsystem(DSargsF)
		Ffixed_points = pp.find_fixedpoints(ode,n=5,maxsearch=1e+6,eps=1e-10)
		#print Ffixed_points
		fixed_points = reduce_fp(Ffixed_points)
		print(fixed_points)

		fig = plt.figure()
		solC_h = nullclines('h',DSargsF,fixed_points,new_domains,1000)
		solC_a = nullclines('A',DSargsF,fixed_points,new_domains,1000)
		plt.plot(solC_h['A'],solC_h['h'],color='b',label='dA/dt=0',linewidth=3)
		plt.plot(solC_a['A'],solC_a['h'],color='g',label='dH/dt=0',linewidth=3)
		for i in range(len(fixed_points)):
			plt.plot(fixed_points[i]['A'],fixed_points[i]['h'])
		#plt.show()
		fig.savefig(title,bbox_inches='tight')
		plt.close()
	except:
		print("Failed for "+title)
