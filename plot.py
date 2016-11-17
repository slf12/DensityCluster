import logging
import numpy as np
import matplotlib.pyplot as plt
from cluster import *
from sklearn import manifold

rho_threshold = -1
delta_threshold = -1

def plot_rho_delta(rho, delta):
	'''
	plot rho delta to choose threshold
	argments:
		rho: local density
		delta: min distance of higher local density
	returns:
		rho_threshold, delta_threshold
	'''
	logger.info('plot: rho delta to choose threshold for rho and delta')
	styles = ['k.', 'g.', 'r.', 'c.', 'm.', 'y.', 'b.']
	assert len(rho) == len(delta)
	plt.figure(0)
	plt.clf()
	plt.plot(rho[1:], delta[1:], styles[0])
	plt.title('rho-delta')
	plt.xlabel('rho')
	plt.ylabel('delta')
	plt.ylim(bottom = 0)
	plt.show()
