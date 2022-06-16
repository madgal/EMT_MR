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
	print("Need  at least 3")
	exit()

from aux_functions import *
from EM_model import *

#####
#####
if myrank==0:
	fig = plt.figure()

	[ll,labs,maxv,minv]=plotCC('A','h','lamdauh',DSargsF,'EMTest',[0,1000],['k'],[1.])
	plt.title('')
	plt.xlim(0,2500)
	plt.ylim(minv,maxv)
	plt.xlabel("ampk",fontsize=20)
	plt.ylabel("hif1",fontsize=20)
	plt.legend(ll,labs,fontsize=20,frameon=False,bbox_to_anchor=(1.04,1))
	plt.savefig("test_MR.png",bbox_inches='tight')
	#plt.show()
	plt.close()

#####
if myrank==1:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('u','mz','lamdauh',DSargsF,'EMTest',[0,1000],['k'],[1.])
	plt.title('')
	#plt.xlim(0,1500)
	plt.ylim(ymin=0)#minv,maxv)
	plt.xlabel("u200 ",fontsize=20)
	plt.ylabel("zeb mrna ",fontsize=20)
	plt.legend(ll,labs)
	plt.savefig("test_EMT.png",bbox_inches='tight')
	#plt.show()
	plt.close()

#####
if myrank==2:
	fig = plt.figure()
	[ll,labs,maxv,minv]=plotCC('A','mz','lamdauh',DSargsF,'EMTest',[0,1000],['k'],[1.])
	plt.title('')
	#plt.xlim(0,1500)
	plt.ylim(ymin=0)#minv,maxv)
	plt.xlabel("ampk ",fontsize=20)
	plt.ylabel("zeb mrna ",fontsize=20)
	plt.legend(ll,labs)
	plt.savefig("test_EMTMR.png",bbox_inches='tight')
	#plt.show()
	plt.close()
