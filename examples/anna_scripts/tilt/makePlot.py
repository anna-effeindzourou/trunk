#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, glob, fileinput, sys
from pylab import *

#### for LaTeX like style
#params = {'backend': 'ps','axes.labelsize': 12,'text.fontsize': 12,'legend.fontsize': 12,'xtick.labelsize': 12,'ytick.labelsize': 12,'text.usetex': True,'font.family': 'serif'}
params = {'lines.linewidth': 3,'backend': 'ps','axes.labelsize': 16,'text.fontsize': 12, 'font.size': 12, 'legend.fontsize': 8,'xtick.labelsize': 14,'ytick.labelsize': 14,'text.usetex': True,'font.family': 'serif'}
rcParams.update(params)

#### Comparison for various Young's moduli
#figure(1)
#### kwarg
#E1={'color':'k','linestyle':'-','label':'$E_S$=5e6Pa (trunk)'}
#E2={'color':'r','linestyle':'--','label':'$E_S$=5e6Pa'}
#kwargs=[E1,E2]

#files=['5n','5e6']

#### kwarg
E0={'color':'r','linestyle':'--','label':r'pty, $\beta_r$:0'}
E1={'color':'r','linestyle':'-','label':r'ptz, $\beta_r$:0'}
E2={'color':'r','linestyle':'-.','label':r'dy, $\beta_r$:0'}
E21={'color':'r','linestyle':'dotted','label':r'dy, $\beta_r$:0'}

E3={'color':'g','linestyle':'--','label':r'pty, $\beta_r$:0.01'}
E4={'color':'g','linestyle':'-','label':r'ptz, $\beta_r$:0.01'}
E5={'color':'g','linestyle':'-.','label':r'dy, $\beta_r$:0.01'}
E51={'color':'g','linestyle':'-','label':r'dy, $\beta_r$:0.01'}

E6={'color':'k','linestyle':'dotted','label':r'pty, $\beta_r$:0.1'}
E7={'color':'k','linestyle':'-','label':r'ptz, $\beta_r$:0.1'}
E8={'color':'k','linestyle':'-.','label':r'dy, $\beta_r$:0.1'}
E81={'color':'k','linestyle':'dotted','label':r'dy, $\beta_r$:0.1'}

E9={'color':'y','linestyle':'dotted','label':'pty,b=0.1,stop, indicator'}
E10={'color':'g','linestyle':'dotted','label':'ptz,b=0.1,stop, indicator'}
E11={'color':'k','linestyle':'dotted','label':'dy,b=0.1,stop, indicator'}




fig1, ax1 = plt.subplots()



# t		alpha		unbF		pbY		pbZ		pby		pbz		ptY		ptZ		pty		ptz		dy		vty		vtz		indicator		v0x		v0y		v0z		v1x		v1y		v1z		v2x		v2y		v2z		v3x		v3y		v3z		v4x		v4y		0  
#v4z		v5x		v5y		v5z		v6x		v6y		v6z		v7x		v7y		v7z		vboxx		vboxy		vboxz		Ec		F0x		F0y		F0z		F1x		F1y		F1z		F2x		F2y		F2z		F3x		F3y		F3z		F4x		F4y		F4z		

#F5x		F5y		F5z		F6x		F6y		F6z		F7x		F7y		F7z		Fboxx		Fboxy		Fboxz		Fnorm		angle_F		angle_box



# t		alpha		unbF		pbY		pbZ		pby		pbz		ptY		ptZ		pty		ptz		dy		vty		vtz		indicator		v0x		v0y		v0z		v1x		v1y		v1z		v2x		v2y		v2z		v3x		v3y		v3z		v4x		v4y		v4z		v5x		v5y		v5z		v6x		v6y		v6z		v7x		v7y		v7z		vboxx		vboxy		vboxz		Ec		F0x		F0y		F0z		F1x		F1y		F1z		F2x		F2y		F2z		F3x		F3y		F3z		F4x		F4y		F4z		F5x		F5y		F5z		F6x		F6y		F6z		F7x		F7y		F7z		Fboxx		Fboxy		Fboxz		Fnorm		angle_F		angle_box
E9={'color':'y','linestyle':'-','label':'v7y'}
E10={'color':'g','linestyle':'-','label':'v7x'}
E11={'color':'k','linestyle':'-','label':'a slope'}

datafile='tgravity_fac_1_fs_40_betaR_0.bz2'
data1 = genfromtxt(datafile,comments='#')
print datafile
#ax1.plot(data1[:,1],data1[:,38],**E9)
ax1.plot(data1[:,1],-data1[:,40],**E10)

a=[]
angl=[]
print len(data1[:,20])
n=1
n2=n
x=[]
y=[]

for i in range(n,len(data1[:,20]),n2):
	if ((data1[i-n,1]-data1[i,1])!=0):
		al=(data1[i-n,2]-data1[i,2])/(data1[i-n,1]-data1[i,1])
		a.append(al)
		angl.append(data1[i,1])
		#y.append(data1[i,20])
	else:
		al=0
		#print 'ouch'
nl=len(a[:])
val=max(a[:])
for i in range(0,nl):
	y.append (val*angl[i])
print len(a[:]),len(y)
#cost=a[v]
print val
ax1.plot(angl[:],a[:],**E11)
#ax1.plot(data1[:,0],data1[:,2],**E6)
#ax1.plot(data1[:,0],data1[:,2],**E51)
#ax1.plot(data1[:,0],data1[:,37],**E5)
#ax1.plot(data1[:,0],data1[:,38],**E6)
#ax1.plot(data1[:,0],data1[:,39],**E9)
ax1.grid()
#ax2.set_xlim([0,1]) 
#ax2.set_title('betaR=0.1')
ax1.legend(ncol=3,loc=0)
ax1.set_xlabel(r'Time [s]')
#ax1.set_xlabel(r'$\alpha$ ($^{\circ}$)')
#ax1.set_ylabel(r'Slope a=\frac{}{}')
##m=.05
##ax2.set_ylim([-m,m]) 
#ax1.set_ylim([-3,0.5]) 
ax1.set_xlim([30,40]) 
#ax1.set_ylim([-0.1,0.5]) 
fig1.savefig('plot_fig1.pdf')



close('all')
