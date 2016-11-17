import logging
import numpy as np
import matplotlib.pyplot as plt
from sklearn import manifold
from cluster import *
logger = logging.getLogger('Plot')

rho_threshold = -1
delta_threshold = -1

def onclick(event):
	global rho_threshold, delta_threshold
	rho_threshold = event.xdata
	delta_threshold = event.ydata
	logger.info('plot: rho threshold : ' + str(rho_threshold) + ' delta threshold: ' + str(delta_threshold))

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
	#plt.clf()
	fig = plt.figure(0)
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	plt.plot(rho[1:], delta[1:], styles[0])
	plt.title('rho-delta')
	plt.xlabel('rho')
	plt.ylabel('delta')
	plt.ylim(bottom = 0)
	plt.show()
	return rho_threshold, delta_threshold

def plot_cluster(dcluster):
	'''
	plot cluster
	argments:
		dcluster: cluster

	'''
	logger.info('plot: cluster result, start multi-dimensional scaling')
	dp = np.zeros((dcluster.max_id, dcluster.max_id), dtype = np.float32)
	cls = []
	for i in xrange(1, dcluster.max_id):
		for j in xrange(i + 1, dcluster.max_id + 1):
			dp[i - 1, j - 1] = dcluster.distance[(i, j)]
			dp[j - 1, i - 1] = dcluster.distance[(i, j)]
		cls.append(dcluster.cluster[i])
	cls.append(dcluster.cluster[dcluster.max_id])
	cls = np.array(cls, dtype = np.float32)
	fo = open(r'./result.txt', 'w')
	center_index = 1
	for index  in xrange(1, dcluster.max_id + 1):
		if dcluster.cluster[index] == -1:
			fo.write(str(index) + ' ' + str(center_index) + ' ' + str(dcluster.halo[index]) + '\n')
			center_index = center_index + 1
			continue
		fo.write(str(index) + ' ' + str(dcluster.cluster[index]) + ' ' + str(dcluster.halo[index]) + '\n')
	fo.close()
	styles = ['k.', 'g.', 'r.', 'b.', 'm.', 'y.', 'c.', 'or', 'ob', 'og', 'ok', '^r', '+r', 'sr', 'dr', '<r', 'pr']

	seed = np.random.RandomState(seed=3)
	mds = manifold.MDS(max_iter=200, eps=1e-4, n_init=1,dissimilarity="precomputed")
	dp_mds = mds.fit_transform(dp)
	#dp_mds = mds.fit_transform(dp)
	logger.info("plot: end mds, start plot")
	xdata = dp_mds[:, 0]
	ydata = dp_mds[:, 1]
	style_list = cls
	clses = set(style_list)
	xs, ys = {}, {}
	for i in xrange(len(xdata)):
		try:
			xs[style_list[i]].append(xdata[i])
			ys[style_list[i]].append(ydata[i])
		except KeyError:
			xs[style_list[i]] = [xdata[i]]
			ys[style_list[i]] = [ydata[i]]
	added = 1
	for idx, cls in enumerate(clses):
		if cls == -1:
			style = styles[0]
			added = 0
		else:
			style = styles[idx + added]
		plt.plot(xs[cls], ys[cls], style)
	plt.title('cluster')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.show()
	logger.info("plot: end plot")
