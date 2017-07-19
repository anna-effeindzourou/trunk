#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, glob, fileinput, sys
from numpy import *


## fedinition for my range NOTE: works only for min step size 0.001
def myrange(begin, end, increment):
	return (arange(int(begin*1000), int(end*1000+1), int(increment*1000))/1000.)

#AnalysisType = 'const'
#lambdaFs = myrange(0.70,1.00,0.01)
#lambdaus = myrange(0.05,0.15,0.01)

AnalysisType = 'rnd'
#lambdaFs = myrange(0.91,0.99,0.02)
lambdaFs = myrange(0.5,1,0.01)
lambdaus = myrange(0.,0.2,0.01)
#lambdaFs = [0.75,0.85]
#lambdaus = [0.20]

if AnalysisType == 'rnd':
	newFile = open("rnd.table", "w")
	newFile.write("description lambdaF lambdau seed ## parameter names\n")
	a_range = range(1,101,1) # values must be integer
	counter=0
	for lambdau in lambdaus:
		for lambdaF in lambdaFs:
			for seed in a_range:
				if seed<10:
					newFile.write( '%s_lF%.2f_lu%.2f_seed00%i %.2f %.2f %i\n'%(AnalysisType,lambdaF,lambdau,seed,lambdaF,lambdau,seed) )
				elif seed<100:
					newFile.write( '%s_lF%.2f_lu%.2f_seed0%i %.2f %.2f %i\n'%(AnalysisType,lambdaF,lambdau,seed,lambdaF,lambdau,seed) )
				else:
					newFile.write( '%s_lF%.2f_lu%.2f_seed%i %.2f %.2f %i\n'%(AnalysisType,lambdaF,lambdau,seed,lambdaF,lambdau,seed) )
				counter+=1
else:
	newFile = open("const.table", "w")
	newFile.write("description lambdaF lambdau seed ## parameter names\n")
	a_range = range(1,101,1) # values must be integer
	counter=0
	for lambdau in lambdaus:
		for lambdaF in lambdaFs:
			newFile.write( '%s_lF%.2f_lu%.2f %.2f %.2f -1\n'%(AnalysisType,lambdaF,lambdau,lambdaF,lambdau) )
			counter+=1

newFile.close()

print 'Number of simulations: ',counter
print 'Estimated time for simulation with -j8 [h]: ', counter*10./8./3600
print 'run e.g. "yade-batch -j8 rnd.table itaTest.py > rnd.log &"'
print 'Check AnalysysType in main script!'