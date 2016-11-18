# DensityCluster
An implement for 'Clustering by fast search and find of density peaks' in science 2014.

implement details:  
	1. now only has a functin that loads data from disatnce file paper providing. You can create a new function to load your data.   
	2. dc (local density threshold) can be choosed auto or not. If not, I choose precent = 2.0 to compute dc. If yes, I will choose a dc that makes precent in 1.0 ~ 2.0.  
	3. When computing rho, you can choose guass or cutoff.  
	4. You can choose density threshold and distance threshold by click on the figure of rho-delta.  
	5. cluster result is in result.txt.  

how to run:  
	python main.py --guass=False --file_name=data\example_distances.dat --auto_select_dc=False  

	guass: use guass or not, if --guass=False, will use cutoff. default is True.
	file_name : distance_file, default is data\example_distances.dat
	auto_select_dc: auto select dc or not. default is False.
