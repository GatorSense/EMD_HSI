%% Script to run EMD on simulated data
% Simulated 2-Dimensional Data Demo: (Simulated data experiment (c) Case 3 of varying number of endmembers from the EMD paper)
% This script will estimate EMD_tot with SED, SAM and SID for 3,4 and 5 endmembers estimated from the ICE algorithm with random initialization.
% Resconstructed data points(computed using the estimated endmembers and estimated proportions values) and estimated endmembers are plotted in each case
%%
%  If this code is used in any publication or presentation, the following
%  reference must be included:
%
% A. Zare, D.T. Anderson, "Earth Movers Distance-Based Simultaneous Comparison of Hyperspectral 
% Endmembers and Proportions," IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, In Press.
%
%%
% Dependencies:  This demo uses SPICE Sparsity Promoting Iterated
% Constrained Endmembers Algorithm Implementation to run ICE.
% If the SPICE Algorithm is referenced in any publication or presentation, the following reference must be cited:
% Zare, A.; Gader, P.; , “Sparsity Promoting Iterated Constrained Endmember Detection in Hyperspectral Imagery,” 
% IEEE Geoscience and Remote Sensing Letters, vol.4, no.3, pp.446-450, July 2007.
%
% SPICE code uses QPC: Quadratic Programming in C (http://sigpromu.org/quadprog/index.html).    
% The Quadratic Programming implementation can be downloaded here: http://sigpromu.org/quadprog/register.php?target=/quadprog/download.php
% and is free for academic/non-commercial use.  If this code is used in any publications, the software must be referenced.  
%%
% Dependencies: The EMD code used here, relies on the pdist function from the MATLAB
% statistics toolbox.  This code also uses the EMD implementation by Y. Rubner found at:
%   http://ai.stanford.edu/~rubner/emd/default.htm
% This code also uses the MATLAB mex interface that can be found at:
%   By M. Alipour
%   http://www.mathworks.com/matlabcentral/fileexchange/12936-emd-earth-movers-distance-mex-interface%
%%
% Author:  Alina Zare
% University of Missouri, Electrical and Computer Engineering
% Email Address: zarea@missouri.edu
% Created: December 2013
% Latest Revision: December 19, 2013
% This product is Copyright (c) 2013 University of Missouri.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF MISSOURI AND
% CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%% 
clear all;close all;clc;
addpath('./qpc');
addpath('/qpc');

% Build mex files for EMD
mex  emd_mex.c emd.c   -O  -v
clc;
%% Generate Proportions and Data
fprintf('(EMD)Earth Movers Distance-Based Simultaneous Comparison of Hyperspectral Endmembers and Proportions \n');
pause(1);
fprintf('Simulated 2-D Data Demo:(Simulated data experiment (c) Case 3 of varying number of endmembers from the EMD paper\n');
fprintf('------------------------------------------------\n');
fprintf('Generating Simulated 2-D Data and True Endmembers\n');
[X,Etrue,Ptrue]=generateSimulatedData();
%% Run ICE with 2 Endmembers and compare results using EMD
fprintf('Running with 2 Endmembers \n');

% ICE Parameters
inputData = X';
parameters.u = 0.0001; %Trade-off parameter between RSS and V term
parameters.gamma = 0; %Sparsity parameter
parameters.M = 2; %Initial number of endmembers
parameters.endmemberPruneThreshold = 0;
parameters.changeThresh = 1e-3; %Used as the stopping criterion
parameters.iterationCap = 5000; %Alternate stopping criterion
parameters.produceDisplay = 0;
initEM = horzcat(Etrue, mean(Etrue, 2), mean(Etrue,2)+.01);
parameters.initEM = initEM(:, 1:parameters.M); 
parameters.options = optimset('Display', 'off');
[Eest_2, Pest_2] = SPICE(inputData, parameters);

% EMD
fprintf('Calculating EMD using SED, SAM, and SID \n');

EMDtot_SED_2 = compareUnmixing(Ptrue, Pest_2, Etrue, Eest_2, 'seuclidean', 0);
EMDtot_SAM_2 = compareUnmixing(Ptrue, Pest_2, Etrue, Eest_2, 'sad', 0);
EMDtot_SID_2 = compareUnmixing(Ptrue, Pest_2, Etrue, Eest_2, 'sid', 0);
fprintf('------------------------------------------------\n');
pause(.3);
%% Run ICE with 3 Endmembers and compare results using EMD
fprintf('Running with 3 Endmembers \n');

% ICE
parameters.M = 3; %Initial number of endmembers
parameters.initEM = initEM(:, 1:parameters.M); 
[Eest_3, Pest_3] = SPICE(inputData, parameters);

% EMD
fprintf('Calculating EMD using SED, SAM, and SID \n');
EMDtot_SED_3 = compareUnmixing(Ptrue, Pest_3, Etrue, Eest_3, 'seuclidean', 0);
EMDtot_SAM_3 = compareUnmixing(Ptrue, Pest_3, Etrue, Eest_3, 'sad', 0);
EMDtot_SID_3 = compareUnmixing(Ptrue, Pest_3, Etrue, Eest_3, 'sid', 0);
fprintf('------------------------------------------------\n');
pause(1);
%% Run ICE with 4 Endmembers and compare results using EMD
fprintf('Running with 4 Endmembers \n');

% ICE
parameters.M = 4; %Initial number of endmembers
parameters.initEM = initEM(:, 1:parameters.M); 
[Eest_4, Pest_4] = SPICE(inputData, parameters);

% EMD
fprintf('Calculating EMD using SED, SAM, and SID \n');
EMDtot_SED_4 = compareUnmixing(Ptrue, Pest_4, Etrue, Eest_4, 'seuclidean', 0);
EMDtot_SAM_4 = compareUnmixing(Ptrue, Pest_4, Etrue, Eest_4, 'sad', 0);
EMDtot_SID_4 = compareUnmixing(Ptrue, Pest_4, Etrue, Eest_4, 'sid', 0);
fprintf('------------------------------------------------\n');
pause(1);
%% Run ICE with 5 Endmembers and compare results using EMD
fprintf('Running with 5 Endmembers \n');

% SPICE
parameters.M = 5; %Initial number of endmembers
parameters.initEM = initEM(:, 1:parameters.M); 
[Eest_5, Pest_5] = SPICE(inputData, parameters);

% EMD
fprintf('Calculating EMD using SED, SAM, and SID \n');
EMDtot_SED_5 = compareUnmixing(Ptrue, Pest_5, Etrue, Eest_5, 'seuclidean', 0);
EMDtot_SAM_5 = compareUnmixing(Ptrue, Pest_5, Etrue, Eest_5, 'sad', 0);
EMDtot_SID_5 = compareUnmixing(Ptrue, Pest_5, Etrue, Eest_5, 'sid', 0);
fprintf('------------------------------------------------\n');
%% Plot Results

% Plot true data and true endmembers
fprintf('Plotting Results\n')
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
subplot(2,4,[ 2 3 ]);
scatter(X(:,1),X(:,2));
hold on;
tmpE = Etrue';
scatter(tmpE(:,1),tmpE(:,2),100,'r','fill');
hold off; 
title('Two Dimensional Simualted Data with Three Endmembers');

% 2 Endmembers- Plot reconstructed data and estimated endmembers
Xrecons_2 = Pest_2*Eest_2';
subplot(2,4,5);
scatter(Xrecons_2(:,1),Xrecons_2(:,2));
hold on;
tmpE = Eest_2';
scatter(tmpE(:,1),tmpE(:,2),100,'r','fill');
hold off;
title({'Result with 2 Endmembers';'Estimated Endmembers and Reconstructed Data Points';['EMD-SED =', (num2str(EMDtot_SED_2)), '; EMD-SAM =', (num2str(EMDtot_SAM_2)),'; EMD-SID =', (num2str(EMDtot_SID_2))]});

% 3 Endmembers- Plot reconstructed data and estimated endmembers
Xrecons_3 = Pest_3*Eest_3';
subplot(2,4,6);
scatter(Xrecons_3(:,1),Xrecons_3(:,2));
hold on;
tmpE = Eest_3';
scatter(tmpE(:,1),tmpE(:,2),100,'r','fill');
hold off;
title({'Result with 3 Endmembers';'Estimated Endmembers and Reconstructed Data Points';['EMD-SED =', (num2str(EMDtot_SED_3)), '; EMD-SAM =', (num2str(EMDtot_SAM_3)),'; EMD-SID =', (num2str(EMDtot_SID_3))]});

% 4 Endmembers- Plot reconstructed data and estimated endmembers
Xrecons_4 = Pest_4*Eest_4';
subplot(2,4,7);
scatter(Xrecons_4(:,1),Xrecons_4(:,2));
hold on;
tmpE = Eest_4';
scatter(tmpE(:,1),tmpE(:,2),100,'r','fill');
hold off;
title({'Result with 4 Endmembers';'Estimated Endmembers and Reconstructed Data Points';['EMD-SED =', (num2str(EMDtot_SED_4)), '; EMD-SAM =', (num2str(EMDtot_SAM_4)),'; EMD-SID =', (num2str(EMDtot_SID_4))]});

% 5 Endmembers- Plot reconstructed data and estimated endmembers
Xrecons_5 = Pest_5*Eest_5';
subplot(2,4,8);
scatter(Xrecons_5(:,1),Xrecons_5(:,2));
hold on;
tmpE = Eest_5';
scatter(tmpE(:,1),tmpE(:,2),100,'r','fill');
hold off;
title({'Result with 5 Endmembers';'Estimated Endmembers and Reconstructed Data Points';['EMD-SED =', (num2str(EMDtot_SED_5)), '; EMD-SAM =', (num2str(EMDtot_SAM_5)),'; EMD-SID =', (num2str(EMDtot_SID_5))]});

%%
clearvars -EXCEPT X Etrue Ptrue EMDtot_SAM_2 EMDtot_SAM_3 EMDtot_SAM_4 EMDtot_SAM_5 EMDtot_SED_2 EMDtot_SED_3 EMDtot_SED_4 EMDtot_SED_5 EMDtot_SID_2 ...
    .. EMDtot_SID_3 EMDtot_SID_3 EMDtot_SID_4 EMDtot_SID_5 Eest_2 Eest_3 Eest_4 Eest_5 Pest_2 Pest_3 Pest_4 Pest_5 scrsz;
%% end of script

