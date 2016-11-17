import logging
from cluster import *
from plot import *
import argparse
def parse_args(): 
    """
    Parse input arguments
    """
    parser = argparse.ArgumentParser(description='density cluster')
    parser.add_argument('--file_name', help='distance file',
                        default='data\example_distances.dat', type=str)
    parser.add_argument('--guass',
                        help='function that compute rho',
                        default=True, type=bool)
    parser.add_argument('--auto_select_dc',
                        help='select dc or not',
                        default=False, type=bool)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

if __name__ == '__main__':
	args = parse_args()
	logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
	dcluster = DensityCluster()
	dcluster.cluster(args.file_name, load_data, args.auto_select_dc, args.guass)
	plot_cluster(dcluster)