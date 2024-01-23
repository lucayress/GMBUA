# A Generalized Multiscale Bundle-Based Hyperspectral Sparse Unmixing Algorithm #

This is the authors' implementation of the paper [1]. If you use this software please cite the following in any resulting
publication:

    [1] A Generalized Multiscale Bundle-Based Hyperspectral Sparse Unmixing Algorithm
        L. C. Ayres, R.A. Borsoi, J.C.M. Bermudez, S.J.M. de Almeida.
        IEEE Geoscience and Remote Sensing Letters, 2024.

The code is implemented in MATLAB R2022a and includes the main files:
-  demo_synthetic_compare.m  
-  demo_cuprite_compare.m	 
-  find_most_consistent_run_abundances.m 
-  genMUA.m
-  genMUA_social.m
-  match_endmembers.m
-  plot_bundles.m
-  plot_hist.m
-  plot_results_cuprite.m
-  plot_results_cuprite_rainclouds.py
-  plot_results_synthetic.m
-  plot_results_synthetic_rainclouds.py
-  ./test_data/               
-  ./real_data/                
-  ./vlfeat-0.9.21/            
-  README.md                 
It also includes code associated to the papers [2] and [3].

## INSTALLING & RUNNING:
Just start MATLAB and run the script demo_synthetic_compare.m to compare the algorithms for the synthetic data (SD) and the script demo_cuprite_compare.m to compare the algorithms for the real data (Cuprite).

The output files will be saved in the ./test_data/results and ./real_data/results directories and the following scripts are used to plot the results:
-  plot_results_cuprite.m
-  plot_results_synthetic.m
-  plot_results_cuprite_rainclouds.py
-  plot_results_synthetic_rainclouds.py  

### Important Notice About VLFeat
If you encounter problems with the "vl_slic" or "vl_setup" functions, try to download the latest version of the toolbox at http://www.vlfeat.org/install-matlab.html. 


## NOTES:
1.  The MUA_SLIC algorithm is provided by Ricardo Borsoi [2]
    at https://github.com/BehnoodRasti/SUnCNN

2.  The SUnCNN algorithm is provided by Behnood Rasti [4]
    at https://github.com/BehnoodRasti/SUnCNN

3.  The raincloud plots are generated in Python from 
	https://github.com/RainCloudPlots/RainCloudPlots

## References:

	[2] A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing
        R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard.
        IEEE Geoscience and Remote Sensing Letters, 2018.
		
	[3] Hyperspectral Image Unmixing with Endmember Bundles and Group Sparsity Inducing Mixed Norms.
		L. Drumetz, T.R. Meyer, J. Chanussot, L.A. Bertozzi and C. Jutten.
		IEEE Transactions on Image Processing, 28(7), 3435-3450, 2019
			
	[4] SUnCNN: Sparse Unmixing Using Unsupervised Convolutional Neural Network
        B. Rasti and B. Koirala
        IEEE Geoscience and Remote Sensing Letters, 2021.		