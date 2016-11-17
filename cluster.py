import logging
import sys
import math
import numpy as np
from plot import *
logger = logging.getLogger('Density Cluster')


def load_data(dis_file):
	'''
	load data from distance file
	arguments: 
		dis_file : file name, the format is 1-column index1, 2-column index 2, 3-column distance
	return: 
		distance dict, max distance, min distance, max index
	'''
	logger.info('progress: load data start')
	distance = {}
	max_dis, min_dis = 0.0, sys.float_info.max
	max_id = 0
	with open(dis_file, 'r') as df:
		for line in df:
			index1, index2, dis = line.split(' ')
			index1 = int(index1)
			index2 = int(index2)
			max_id = max(max_id, index1, index2)
			dis = float(dis)
			min_dis, max_dis = min(min_dis, dis), max(max_dis, dis)
			distance[(index1, index2)] = dis
			distance[(index2, index1)] = dis

	for i in xrange(max_id):
		distance[(i, i)] = 0.0
	logger.info('progress: load data end')
	return distance, max_dis, min_dis, max_id

def select_dc(distance, max_id, max_dis, min_dis, auto_select_dc):
	'''
	select local density threshold, default is used in the paper, is 2.0, if you do not want use this, you can set auto_select_dc = True
	arguments: 
		distance: distance dict
		max_id : max index
		max_dis: max distance in distance
		min_dis: min distance in distance\
		auto: auto select dc or not
	returns:
		dc : local density threshold
	'''
	if auto_select_dc:
		return auto_select_dc_fuc(distance, max_id, max_dis, min_dis)
	precent = 2.0
	position = int(max_id * (max_id - 1) / 2 * precent / 100)
	dc = sorted(distance.values())[position * 2 + max_id] # because there is double of distance, and distance(i,i) is 0.0
	logger.info('progress: dc is ' + str(dc))
	return dc

def auto_select_dc_fuc(distance, max_id, max_dis, min_dis):
	'''
	auto select local density threshold
	'''
	dc = (min_dis + max_dis) / 2
	while  True:
		precent = float(sum([1 for v in distance.values() if v < dc])) / max_id ** 2
		if precent >= 0.01 and precent <= 0.02:
			break
		if precent < 0.01:
			min_dis = dc
		else :
			max_dis = dc
		dc = (max_dis + min_dis) / 2
		if max_dis - min_dis < 0.0001:
			break
	return dc

def compute_local_density(distance, max_id, dc, guass):
	'''
	compute local density
	arguments:
		distance: distance dict
		max_id: max index
		dc: local density threshold
		guass: use gauss func or not
		cutoff: use cutoff func or not
	return:
		local density vector that index is point index that starts from 1
	'''
	logger.info('progress: local density start')
	guass_func = lambda dij, dc : math.exp(- (dij / dc) ** 2)
	cutoff_func = lambda dij, dc: 1 if dij < dc else 0
	func = guass and guass_func or cutoff_func
	rho = [-1] + [0] * max_id # index start from 1
	for i in xrange(1, max_id):
		for j in xrange(i + 1, max_id+1):
			rho[i] += func(distance[(i, j)], dc)
			rho[j] += func(distance[(i, j)], dc)
		if i % (max_id / 10) == 0:
			logger.info('at index ' + str(i))
	logger.info('progress: local density end')	
	return np.array(rho, np.float32)

def min_distance_for_higher_point(distance, rho, max_id, max_dis):
	'''
	compute all points's min distance from the points that have higher local density 
	Arguments:
		distance: distance dict
		rho: local density
		max_id: max index
		max_dis: max distance
	returns:
		min_distance vector, nearest neighbor vector
	'''
	logger.info("progess: compute min distance to nearest higher density neighbor")
	rho_sort_index = np.argsort(-rho)
	delta = [0.0] + [float(max_dis)] * (len(rho) - 1)
	nearest_neighbor = [0] * len(rho)
	delta[rho_sort_index[0]] = -1
	for i in xrange(1, max_id):
		for j in xrange(0, i):
			original_i = rho_sort_index[i]
			original_j = rho_sort_index[j]
			if distance[(original_i, original_j)] < delta[original_i]:
			 	delta[original_i] = distance[(original_i, original_j)]
			 	nearest_neighbor[original_i] = original_j

		if i % (max_id / 10) == 0:
			logger.info('progress: at index ' + str(i))
	delta[rho_sort_index[0]] = max(delta)
	return np.array(delta, np.float32), np.array(nearest_neighbor, np.float32)

def  choose_cluster_center(delta, rho, nearest_neighbor, rho_threshold, delta_threshold):
	'''
	get cluster, ccenter
	arguments:
		delta: distance
		rho: local density
		nearest_neighbor: index that has higher density and min distance
		rho_threshold: local density threshold 
		delta_threshold: distance threshold
	returns:
		cluster center id
	'''
	logger.info('progress: choose cluster center start')
	cluster = [0] * len(rho)
	ccenter = [-1]
	flag = 0
	for index, (density_temp, distance_temp, neighbor) in enumerate(zip(rho, delta, nearest_neighbor)):
		if index == 0:
			continue
		if density_temp >= rho_threshold and distance_temp >= delta_threshold:
			flag = flag + 1
			ccenter.append(index)
			cluster[index] = flag
	for value in ccenter:
		if value == -1:
			continue
		logger.info('id of cluster center ' + str(value))
	logger.info('progress: choose cluster center end')
	return cluster, ccenter

def do_cluster(cluster, ccenter, delta, rho, nearest_neighbor):
	'''
	progress the points that is not center
	arguments:
		cluster: cluster
		ccenter: index of center
		delta: min distance
		rho: local density
		nearest_neighbor: index that point's nearest neighbor
	'''
	logger.info('progress: cluster start')
	rho_sort_index = np.argsort(-rho)
	for index in rho_sort_index:
		if index == 0 or index in ccenter:
			continue
		cluster[index] = cluster[int(nearest_neighbor[index])]
	for value in ccenter:
		if value == -1:
			continue
		cluster[value] = -1

	logger.info('progress: cluster end')
	return cluster, ccenter

def halo(cluster, ccenter, distance, dc, max_id, rho):
	'''
	compute halo, 0 is center, 1 is halo
	arguments:
		cluster: cluster
		ccenter: index of center
		distance: distance dict
		max_id: max index
		dc: local density threshold
	return:
		halo vector
	'''
	if len(ccenter) < 2:
		halo = [0] * (max_id + 1)
		return halo
	halo = [1] * (max_id + 1)
	bord_rho = [0.0] * len(ccenter)
	for i in xrange(1, max_id):
		for j in xrange(1 + i, max_id + 1):
			if (cluster[i] != cluster[j]) and (distance[(i, j)] <= dc):
				rho_aver = (rho[i] + rho[j]) / 2.0
				if rho_aver > bord_rho[cluster[i]]:
					bord_rho[cluster[i]] = rho_aver
				if rho_aver > bord_rho[cluster[j]]:
					bord_rho[cluster[j]] = rho_aver
	print bord_rho
	for i in xrange(1, max_id + 1):
		if rho[i] < bord_rho[cluster[i]]:
			halo[i] = 0
	return halo

class DensityCluster( object ):
	def get_local_density(self, load_data_func, distance_file, auto_select_dc, guass):
		'''
		use function defined before together to get the local density
		Arguments:
			load_data_func: the function that load the data from file, we can define different function to load different distance
			distance_file: the distance file
			dc: local density threshold 
			auto_select_dc: auto select dc or not
		return: 
			distance: distance dict
			max_dis: max distance
			min_dis: min distance
			max_id: max index
			rho: loacl density
		'''
		distance, max_dis, min_dis, max_id = load_data_func(distance_file)
		dc = select_dc(distance, max_id, max_dis, min_dis, auto_select_dc)
		rho = compute_local_density(distance, max_id, dc, guass)

		return distance, max_dis, min_dis, max_id, rho, dc


	def cluster(self, distance_file,  load_data_func, auto_select_dc, guass):
		'''
		cluster 
		arguments:
			distance_file: distance file
			load_data_func: load data function 
			density_threshold: density threshold that  use when choose cluster center
			distance_threshold: distance threshold that use when choose cluster center
			dc : local density threshold
			auto_select_dc : auto select dc or not
		returns:

		'''
		distance, max_dis, min_dis, max_id, rho, dc = self.get_local_density(load_data_func, distance_file, auto_select_dc, guass)
		delta, nearest_neighbor = min_distance_for_higher_point(distance, rho, max_id, max_dis)
		rho_threshold, delta_threshold = plot_rho_delta(rho, delta)
		if rho_threshold == -1:
			rho_threshold = 20
			delta_threshold = 0.1
		cluster, ccenter = choose_cluster_center(delta, rho, nearest_neighbor, rho_threshold, delta_threshold)
		cluster, ccenter = do_cluster(cluster, ccenter, delta, rho, nearest_neighbor)
		halo_v = halo(cluster, ccenter, distance, dc, max_id, rho)
		self.cluster = cluster
		self.ccenter = ccenter
		self.max_id = max_id
		self.distance = distance
		self.halo = halo_v

