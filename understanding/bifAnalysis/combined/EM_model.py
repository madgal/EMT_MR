import PyDSTool
from PyDSTool.Toolbox import phaseplane as pp
import matplotlib
from pylab import *
import copy
import math
matplotlib.rcParams.update({'font.size': 30})
from matplotlib.lines import Line2D
from  parameters_crosstalk import *

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
		  'Yu3' : (['X','X0','nX'],'sum(i,0,2,yui[i]*combinationU([i],nX)*M([i],nX,X,X0))'),
		  'CompRmt':(['yyX','gnX','hX','h0rmtX','nhmtX','AX','A0rmtX','namtX'],'yyX*(gnX+(AX/A0rmtX)**namtX)/(1.+(hX/h0rmtX)**nhmtX+(AX/A0rmtX)**namtX)'),
		  'CompRn':(['g0X','hX','h0rnX','nhnX','ghnX','AX','A0rnX','ganX','nanxX'],'(g0X+ghnX*(hX/h0rnX)**nhnX+ganX*(AX/A0rnX)**nanxX)/(1.+(hX/h0rnX)**nhnX+(AX/A0rnX)**nanxX)'),   
		  'R':(['R1','R2'],'R1+R2')}

DSargsF.pars = {	'gz':gz, 'u0':u0,'nu':nu,'kz':kz,
                'li0':li0, 'li1':li1, 'li2':li2, 'li3':li3, 'li4':li4, 'li5':li5, 'li6':li6,
                'ymi0':ymi0, 'ymi1':ymi1, 'ymi2':ymi2, 'ymi3':ymi3, 'ymi4':ymi4, 'ymi5':ymi5, 'ymi6':ymi6,
                'yui0':yui0, 'yui1':yui1, 'yui2':yui2, 'yui3':yui3, 'yui4':yui4, 'yui5':yui5, 'yui6':yui6,
		'gu':gu, 'Z0u':Z0u,'nzu':nzu, 'lamdazu':lamdazu,'S0u':S0u,'nsu':nsu,'lamdaSu':lamdaSu,'u0':u0,'nu':nu,'ku':ku,
		'gmz':gmz, 'Z0m':Z0m,'nzm':nzm, 'lamdaZm':lamdaZm,'S0m':S0m,'nsm':nsm,'lamdaSm':lamdaSm,'u0':u0,'nu':nu,'kmz':kmz,
		'gms':gms, 'I0m':I0m, 'nIm':nIm, 'lamdaIm':lamdaIm, 'S0ms':S0ms, 'nsms':nsms, 'lamdaSms':lamdaSms, 'u30':u30,
                'nu3':nu3, 'kms':kms, 'gu3':gu3, 'Z0u3':Z0u3, 'nzu3':nzu3, 'lamdazu3':lamdazu3, 'S0u3':S0u3, 'nsu3':nsu3,
		'lamdaSu3':lamdaSu3, 'ku3':ku3, 'gs':gs, 'ks':ks, 'lamdahs':7., 'nhs':2., 'h0hs':200. ,
                'ga':ga,'ka':ka,'A0aa':A0aa,'A0ah':A0ah,'A0aR':A0aR,'A0rn':A0rn,'A0rm':A0rm,'lambdaaa':lambdaaa,'lambdaah':lambdaah,'lambdaar':lambdaar,'y':y,
		'g2':g2,'naa':naa,'nah':nah,'nar':nar,'narm':narm,'narn':narn,'gh':gh,'kh':kh,'h0hh':h0hh,'h0ha':h0ha,'h0hrm':h0hrm,
		'h0hrn':h0hrn,'lambdahh':lambdahh,'lambdaha':lambdaha,'g1':g1,'nhh':nhh,'nha':nha,'nhrm':nhrm,'nhrn':nhrn,'grm':grm,'grn':grn,'gn':gn,
		'krm':krm,'krn':krn,'R0ra':R0ra,'R0rh':R0rh, 'lambdara':lambdara,'lambdarh':lambdarh,'nra':nra,'nrh':nrh,
		'y':8,'kh':0.25,'u0uh':10000.,'nuh':1,'lamdauh':0.8,'u0u3R':10000,'nu3nR':4,'lamdau3nR':2.}
	
DSargsF.varspecs = { 'u':'gu*H(Z(mz,u,u0,nu,kz,gz),Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)-mz*Yu(u,u0,nu)-ku*u',
		   'mz':'gmz*H(Z(mz,u,u0,nu,kz,gz),Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)-mz*Ym(u,u0,nu)-kmz*mz',
		   'ms':'gms*H(h,h0hs,nhs,1.)*H(S,S0ms,nsms,lamdaSms)-ms*Ym3(u3,u30,nu3)-kms*ms',
                   'u3':'gu3*H(Z(mz,u,u0,nu,kz,gz),Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu3(u3,u30,nu3)-ku3*u3',
                   'S':'gs*ms*L3(u3,u30,nu3)-ks*S',
                   'A':'ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A',
		   'Rmt':'grm*H(A,A0aR,nar,lambdaar)*H(u3,1,1,1.)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt',
		   'Rnox':'grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox*H(u3,1.,1.,1.)',
		   'h':'gh*H(A,A0ah,nah,lambdaah)*H(u,u0uh,nuh,lamdauh)-kh*h*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)' 
		}
			
DSargsF.ics = {'u':1000,'mz':1000,'S':200000,'ms':2000,'u3':4000,'A':300,'h':300,'Rmt':300,'Rnox':300}
DSargsF.xdomain = {'u':[0,100000],'mz':[0,100000],'Z':[0,100000],'S':[0,400000],'ms':[0,100000],'u3':[0,100000],'A':[0,5000],'h':[0,5000],'Rmt':[0,5000],'Rnox':[0,5000]}
DSargsF.tdomain = [0,5000]
DSargsF.algparams = {'init_step':0.1}
