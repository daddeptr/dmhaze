import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import pylab as pl

dir = '/global/scratch2/sd/dpietrob/DM-haze/'

mass_val = [ 6.99, 11.6, 19.1, 31.7, 52.4, 73.3, 103, 144, 201, 281, 393, 550, 769, 1080]

tag = 'bubbles'
one_plot = False
six_plot = True

tc2 = []

if six_plot:
	fig = plt.figure(1, (14,8), dpi=80)
	cnt = 0
	for channel in [1, 2]:
		for slope in range(1,4):
			cnt += 1
			pt = int( str(2)+str(3)+str(cnt) )
			ax = plt.subplot(pt)
			m = []
			c2 = []
			sig = []
			low = []
			up = []
			for mass in range(1,6):
				for mag in [13, 21, 32, 4, 5]:
					model = str(mass) + str(channel) + str(slope) + str(mag) + "h"
					print model
					bf = open( dir+'results/mcmc_chain_'+tag+'_'+model+'_2stp_08x.likestats', 'r')
					lines = bf.readlines()
					bf.close()
					m.append(mass_val[mass-1])
					tc2.append( float( lines[0].split('=')[1] ) )
					c2.append( float( lines[0].split('=')[1] ) )
					sig.append( float( lines[3].split()[1] ) )
					low.append( float( lines[3].split()[2] ) )
					up.append( float( lines[3].split()[3] ) )

			ler = np.array(sig)-np.array(low)
			uer = np.array(up)-np.array(sig)
			color = np.array( ( np.array(c2) - np.min(c2) ) )  #/ 100.
			print np.min(c2), np.max(c2)
			print color, np.min(color), np.max(color)

	# for i, c in enumerate(color):
	# 	plt.errorbar( [m[i]], [sig[i]], yerr=[[ler[i]],[uer[i]]], fmt='-', color='k' )
	# plt.xlabel('Mass [GeV]')
	# plt.ylabel(r'$\sigma$')

	# ax = plt.subplot(122)
			plt.scatter(m, sig, c=color, marker='d')
			if cnt > 3:
				plt.xlabel('Mass [GeV]')
			if cnt == 1 or cnt == 4:
				plt.ylabel(r'$\sigma$')
			if channel == 1:
				ax.set_ylim([0.9,11])
			elif channel == 2:
				ax.set_ylim([0.1,2])
			ax.set_xlim([5,60])
			# ax.set_yscale('log')
			ax.set_xscale('log')
			plt.title('Decaying channel '+str(channel)+', slope '+str(slope))
			#from matplotlib.text import OffsetFrom	
			#offset_from = OffsetFrom(an1, (0.5, 0))
			# if cnt == 1:
			#	ax.annotate("Channel "+str(channel), xy=(10, 2))

# 			cdict = {
# 			  'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
# 			  'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
# 			  'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
# 			}
# 
# 			cm = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
# 			cb = plt.colorbar( cmap=cm )
 			cb = plt.colorbar()
			cb.ax.xaxis.set_label_text(r'$\Delta\chi^2$')


if one_plot:
	for channel in [1, 2]:
		for slope in range(1,4):
			for mass in range(1,6):
				for mag in [13, 21, 32, 4, 5]:
					model = str(mass) + str(channel) + str(slope) + str(mag)
					print model
					bf = open( dir+'results/mcmc_chain_'+tag+'_'+model+'_2stp_08x.likestats', 'r')
					lines = bf.readlines()
					bf.close()
					m.append(mass_val[mass-1])
					c2.append( float( lines[0].split('=')[1] ) )
					sig.append( float( lines[3].split()[1] ) )
					low.append( float( lines[3].split()[2] ) )
					up.append( float( lines[3].split()[3] ) )

	fig = plt.figure(1, (10,10))
	ax = plt.subplot(111)
	ler = np.array(sig)-np.array(low)
	uer = np.array(up)-np.array(sig)
	color = np.array( ( np.array(c2) - np.min(c2) ) )  #/ 100.
	print np.min(c2), np.max(c2)
	print color, np.min(color), np.max(color)

	# for i, c in enumerate(color):
	# 	plt.errorbar( [m[i]], [sig[i]], yerr=[[ler[i]],[uer[i]]], fmt='-', color='k' )
	# plt.xlabel('Mass [GeV]')
	# plt.ylabel(r'$\sigma$')

	# ax = plt.subplot(122)
	plt.scatter(m, sig, c=color, marker='d')
	plt.xlabel('Mass [GeV]')
	plt.ylabel(r'$\sigma$')
	ax.set_ylim([0.9,11])
	ax.set_xlim([5,60])
	# ax.set_yscale('log')
	ax.set_xscale('log')


	cdict = {
	  'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
	  'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
	  'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
	}

	cm = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
	cb = plt.colorbar( cmap=cm )
	cb.ax.xaxis.set_label_text(r'$\Delta\chi^2$')

# --- from setup_matplotlib
# 
# import numpy as np
# from matplotlib import rcParams, rc
# 
# # common setup for matplotlib
# params = {'backend': 'pdf',
#           'savefig.dpi': 300, # save figures to 300 dpi
#           'axes.labelsize': 10,
#           'text.fontsize': 10,
#           'legend.fontsize': 10,
#           'xtick.labelsize': 10,
#           'ytick.major.pad': 6,
#           'xtick.major.pad': 6,
#           'ytick.labelsize': 10,
#           #'text.usetex': True,
#           #'font.family':'sans-serif',
#           # free font similar to Helvetica
#           'font.sans-serif':'FreeSans'}
# 
# # use of Sans Serif also in math mode
# rc('text.latex', preamble='\usepackage{sfmath}')
# 
# rcParams.update(params)
# 
# def cm2inch(cm):
#     """Centimeters to inches"""
#     return cm *0.393701
# 
# # ---
# from matplotlib import ticker
# def format_func(x, pos):
#     if np.abs(x) < 90:
#         formatted = "%d" % x
#     else:
#         formatted = "%d^{%d}" % (10*np.sign(x), int(np.log10(np.abs(x))))
#     out = r"$%s$" % formatted 
#     return out 
# formatter = ticker.FuncFormatter(format_func)
# colorbar_ticks = np.arange(0,100,10)
# colorbar_boundaries = np.concatenate([-1 * np.logspace(0, 3, 80)[::-1], np.linspace(-1, 1, 10), np.logspace(0, 7, 150)])
# 
# cdict = {
#   'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
#   'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
#   'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
# }
# 
# cm = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
# 
# def do_plot(n, f, title):
#     #plt.clf()
#     plt.subplot(1, 1, n)
#     plt.pcolor(m, sig, f(c2), cmap=cm, vmin=-4, vmax=4)
#     plt.title(title)
#     plt.colorbar()
# 
# do_plot(1, lambda x:x, "all")
# 
# cb = mpl.colorbar.ColorbarBase( ax, cmap=cm,
# 	orientation='horizontal',
# 	boundaries=colorbar_boundaries,
# 	ticks=colorbar_ticks,
# 	format=formatter)
# 	
#cb.ax.xaxis.set_label_text('Mass [GeV]')

pl.show()


