import matplotlib.pyplot as plt
import pylab as pl
import pickle
import healpy as hp
import sys

from myclass import *
from mcmc import *
from data_model import *

# ====================================================================
# selection = ['3225','3225h','3235','5114','12232','32232','myhaze']
# model_list = [ 'myhaze' ]
# print model_list

tag = 'southern'	

dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/'
# file = dir + '../../../sed_dm_myhaze_2step_'+tag+'_fixCMB_noDISC.txt'
file = dir + '../../../sed_dm_myhaze_1step_'+tag+'.txt'

cfreq, haze_flux, haze_flux_er, haze_flux_er_scaled = np.loadtxt( file, usecols=(0,1,2,3), unpack=True )
# print cfreq, haze_flux, haze_flux_er
haze_flux_er = haze_flux_er_scaled

# model_file = dir + 'all-old_haze.txt'
model_file = dir + 'all-new_haze.txt'
content = open( model_file, 'r' )
models = {}
models_h = {}
for line in content:
	if line[0] != ';':
		line = line.strip('\n')
		fields = line.strip('\n').split()
		# print fields
		# print len( fields )
		key = fields[0].strip().split('_')[3]
		if model_file.find('new') >=0:
			value = np.array( fields[1:], dtype=float )
		else:
			value = np.array( fields[1:-1], dtype=float )
		if key[-1] == 'h':
			models_h[key] = value
		else:
			models[key] = value

# --- Log plot ---
# ax = plt.subplot(122)
# plt.plot( cfreq, haze_flux, 's-', label='Haze flux', color='m' )
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_ylim([1.e-4,1])
# ax.set_xlim([8,100])
# plt.xlabel(r'$\nu$ [GHz]')
# plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')

fig = plt.figure(1,(13,8), dpi=80)
ax = plt.subplot(121)
plt.errorbar( cfreq, haze_flux, yerr=haze_flux_er, label='Haze flux', color='m', fmt='s-' )
plt.plot( cfreq, cfreq*0., '--', color='k')

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim([1.e-4,1])
ax.set_xlim([8,100])
plt.xlabel(r'$\nu$ [GHz]')
plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')


first = True
for key in models_h.keys():
	value = models_h[key]
	# print value
	if first:
		plt.plot( cfreq, value*1.e20, ':', color='b', label='Theory h')
		first = False
	else:
		plt.plot( cfreq, value*1.e20, ':', color='b')

plt.legend( loc='upper left',prop={'size':12})

ax = plt.subplot(122)
plt.errorbar( cfreq, haze_flux, yerr=haze_flux_er, label='Haze flux', color='m', fmt='s-' )
plt.plot( cfreq, cfreq*0., '--', color='k')

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim([1.e-4,1])
ax.set_xlim([8,100])
plt.xlabel(r'$\nu$ [GHz]')
plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')

first = True
for key in models.keys():
	value = models[key]
	if first:
		plt.plot( cfreq, value*1.e20, ':', color='c', label='Theory')
		first = False
	else:
		plt.plot( cfreq, value*1.e20, ':', color='c')

plt.legend( loc='upper left',prop={'size':12})
fig.savefig( '../dm-model_vs_data_'+tag+'.png' )

# --- Amplitude (cross-section) fit ---
# I think I can use themap_regression, since it's the same algorithm, but to check it,
# I can use the np.linalg routine.
# Actually it's much faster to write the solution

A = {}
A_h = {}
# --- Log plot ---
fig = plt.figure(2,(10,8), dpi=80)
ax = plt.subplot(111)
# plt.plot( cfreq, haze_flux, 's-', label='Haze flux', color='m' )
plt.errorbar( cfreq, haze_flux, yerr=haze_flux_er, fmt='s-', label='Haze flux', color='m' )
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim([1.e-2,10])
ax.set_xlim([20,80])
plt.xlabel(r'$\nu$ [GHz]')
plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')

first = True
min_c2 = ( '',1.e30, 1. )
for key in models.keys():
	value = models[key] * 1.e20
	a = np.sum( haze_flux*value/haze_flux_er**2 ) / np.sum( value**2/haze_flux_er**2 )
	aer = np.sqrt( 1. / np.sum( value**2/haze_flux_er**2 ) )
	c2 = np.sum( (haze_flux-a*value) / haze_flux_er**2 )
	A[key] = (a, aer, c2)
	if c2 < min_c2[1]:
		min_c2 = ( key, c2, a, aer )
	if first:
		plt.plot( cfreq, value*a, '--', color='b', label='Theory')
		first = False
	else:
		plt.plot( cfreq, value*a, '--', color='b')

bf = models[min_c2[0]] * 1.e20 * min_c2[2]
plt.plot( cfreq, bf, 'd-', color='r', label=r'Best fit: '+min_c2[0]+', $\sigma$='+str('%.3f' %min_c2[2])+r'$\pm$'+str('%.3f' %min_c2[3] ) )
plt.legend( loc='upper left',prop={'size':12})

first = True
min_c2_h = ( '',1.e30, 1. )
for key in models_h.keys():
	value = models_h[key] * 1.e20
	a = np.sum( haze_flux*value/haze_flux_er**2 ) / np.sum( value**2/haze_flux_er**2 )
	aer = np.sqrt( 1. / np.sum( value**2/haze_flux_er**2 ) )
	c2 = np.sum( (haze_flux-a*value) / haze_flux_er**2 )
	A_h[key] = (a, aer, c2)
	if c2 < min_c2_h[1]:
		min_c2_h = ( key, c2, a, aer )
	if first:
		plt.plot( cfreq, value*a, '--', color='c', label='Theory h')
		first = False
	else:
		plt.plot( cfreq, value*a, '--', color='c')

bf = models_h[min_c2_h[0]] * 1.e20 * min_c2_h[2]
plt.plot( cfreq, bf, 'd-', color='k', label=r'Best fit: '+min_c2_h[0]+', $\sigma$='+str('%.3f' %min_c2_h[2])+r'$\pm$'+str('%.3f' %min_c2_h[3] ) )
plt.legend( loc='upper left',prop={'size':12})

fig.savefig( '../dm-model_fit_'+tag+'.png' )

# for key in A.keys():
# 	print key, A[key]

print min_c2
print min_c2_h

pl.show()



