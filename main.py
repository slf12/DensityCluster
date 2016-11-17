import logging
from plot import *
from cluster import *

if __name__ == '__main__':
	logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
	dpcluster = DensityPeakCluster()
	dpcluster.cluster('data\example_distance.dat', load_data, auto_select_dc = False)