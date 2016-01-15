# EMD_HSI
Earth Movers Distance-Based Simultaneous Comparison of Hyperspectral Endmembers and Proportions

Matlab implementation

****************************************************************
****************************************************************
NOTE: If the EMD Algorithm is used in any publication or presentation, the following reference must be cited:

A. Zare and D. T. Anderson, "Earth Movers Distance-Based Simultaneous Comparison of Hyperspectral Endmembers and Proportions," IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing , In Press.

****************************************************************
****************************************************************
Note: This code uses MEX files.  Your MEX compiler will need to be set up accordingly to be able to run this code.  The first time your code is run, set the buildFlag parameter to 1.  After the initial compilation, this buildFlag can be set to 0. 

The command to run the EMD code: 

EMDtot = compareUnmixing(P1, P2, E1, E2, gdtype, buildFlag) 
Inputs:
   P1 - NxM matrix - each row is the proportion vector for one data point
        given endmembers E1
   P1 - NxL matrix - each row is the proportion vector for one data point
        given endmembers E2
   E1 - DxM matrix - each column is an endmember of dimension D
   E2 - DxL matrix - each column is an endmember of dimension D
   gdtype - string - Indicates the ground distance to use to compare the
        endmember sets.  Options are: 'seuclidean' (squared Euclidean
        distance), 'euclidean' (Euclidean distance), 'sad' (Spectral Angle
        Distance), or 'sid' (Spectral Information Divergence).
   buildFlag - Set to 1 if mex files need to be compiled, 0 otherwise.

Outputs:
    EMDtot - Aggregated EMD value (Aggregation is the sum in this implementation)

****************************************************************
Run the Simulated 2-Dimensional Data Demo: (Simulated data experiment (c) Case 3 of varying number of endmembers from the EMD paper)

To Run the EMD Algorithm, with the example simulated 2-dimensional data set run the script: DEMO_RunSimulatedData.m

This script will estimate EMD_tot with SED, SAM and SID with 3,4 and 5 endmembers estimated from the ICE algorithm with random initialization 
   
****************************************************************
****************************************************************
Code Dependencies: 

EMD Dependencies: The EMD code used here, relies on the pdist function from the MATLAB
 statistics toolbox.  This code also uses the EMD implementation by Y. Rubner found at:
   http://ai.stanford.edu/~rubner/emd/default.htm
 This code also uses the MATLAB mex interface that can be found at:
   By M. Alipour
   http://www.mathworks.com/matlabcentral/fileexchange/12936-emd--earth-movers-distance--mex-interface

DEMO Dependencies: The 2-D simulated demo uses SPICE Sparsity Promoting Iterated Constrained Endmembers Algorithm Implementation to run the ICE algorithm.  
 If the SPICE Algorithm is used in any publication or presentation, the following reference must be cited:
 Zare, A.; Gader, P.; , "Sparsity Promoting Iterated Constrained Endmember Detection in Hyperspectral Imagery," IEEE Geoscience and Remote Sensing Letters, vol.4, no.3, pp.446-450, July 2007.

 SPICE code uses QPC: Quadratic Programming in C (http://sigpromu.org/quadprog/index.html).    
 The Quadratic Programming implementation can be downloaded here: http://sigpromu.org/quadprog/register.php?target=/quadprog/download.php
 and is free for academic/non-commercial use.  If this code is used in any publications, the software must be referenced. 
****************************************************************
****************************************************************

If you have any questions, please contact:

% Author:  Alina Zare
% University of Missouri, Electrical and Computer Engineering
% Email Address: zarea@missouri.edu
% Created: December 2013
% Latest Revision: December 19, 2013
% This product is Copyright (c) 2013 University of Missouri.
% All rights reserved.



