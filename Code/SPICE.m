function [endmembers, P] = SPICE(inputData, parameters)

% SPICE Sparsity Promoting Iterated Constrained Endmembers Algorithm
%       Finds Endmembers and Unmixes Input Data
%
% Syntax: [endmembers, P] = SPICE(inputData, parameters)
%
% Inputs:
%   inputData - double Mat - NxM matrix of M data points of
%       dimensionality N (i.e.  M pixels with N spectral bands, each pixel is
%       a column vector)
%   parameters - struct - The struct contains the following fields:
%                   1. u : Regularization Parameter for RSS and V terms
%                   2. gamma: Gamma Constant for SPT term
%                   3. changeThresh: Stopping Criteria, Change threshold
%                       for Objective Function.
%                   4. M: Initial Number of endmembers
%                   5. iterationCap: Maximum number of iterations
%                   6. endmemberPruneThreshold: Proportion threshold used
%                      to prune endmembers
%                   7. produceDisplay : Set to 1 if a progress display is
%                       wanted
%                   8. initEM: Set to nan to randomly select endmembers,
%                       otherwise NxM matrix of M endmembers with N spectral
%                       bands, Number of endmembers must equal parameters.M
% Outputs:
%   endmembers - double Mat - NxM matrix of M endmembers with N spectral
%       bands
%   P - double Mat - NxM matrix of abundances corresponding to M input
%       pixels and N endmembers
% Other m-files required: unmix, Matlab Optimization Toolbox
%
% Author: Alina Zare
% University of Missouri, Electrical and Computer Engineering
% Email Address: azare@ufl.edu
% Created: August 2006
% Latest Revision: November 22, 2011
% This product is Copyright (c) 2013 University of Missouri and University
% of Florida
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

addpath('.\qpc');
parameters.pruningIteration = 1;
parameters.options = optimset('Display', 'off');
M = parameters.M;
X = inputData;

if(isnan(parameters.initEM))
    %Find Random Initial Endmembers
    randIndices = randperm(size(inputData,2));
    randIndices = randIndices(1:parameters.M);
    endmembers = inputData(:,randIndices);
    parameters.initEM = endmembers;
else
    %Use endmembers provided
    M = size(parameters.initEM, 2);
    endmembers = parameters.initEM;
end


%N is the number of pixels, RSSreg is the current objective function total.
N = size(X,2);
RSSreg = inf;
change = inf;

iteration = 0;
P = ones(N,M)*(1/M);
lambda = N*parameters.u/((M-1)*(1-parameters.u));
Im = eye(M);
I1 = ones([M,1]);
while( change > parameters.changeThresh && iteration < parameters.iterationCap)
    
    iteration = iteration + 1;
    
    %Given Endmembers, minimize P -- Quadratic Programming Problem
    [P] = unmix2(X, endmembers, parameters.gamma, P);
    
    %Given P minimize Endmembers
    endmembersPrev = endmembers;
    endmembers = ((P'*P + lambda*(Im - (I1*I1')/M))\(P'*X'))';
    
    %Prune Endmembers below pruning threshold
    pruneFlag = 0;
    pruneIndex = max(P)<parameters.endmemberPruneThreshold;
    minmaxP = min(max(P));
    if(sum(pruneIndex) > 0)
        pruneFlag = 1;
        endmembers = endmembers(:,logical(1-pruneIndex));
        P = P(:, logical(1-pruneIndex));
        M = M - sum(pruneIndex);
        lambda = N*parameters.u/((M-1)*(1-parameters.u));
        Im = eye(M);
        I1 = ones([M,1]);
    end
    
    %Calculate RSSreg (the current objective function value)
    sqerr = (X-(endmembers*P')).^2; 
    RSS = sum(sum(sqerr));
    V = sum(sum(endmembers.*endmembers,2) - (1/M)*sum(endmembers,2).*2)/(M-1);
    SPT = M*parameters.gamma;
    RSSprev = RSSreg;
    RSSreg = (1-parameters.u)*(RSS/N) + parameters.u*V + SPT;
    
    %Determine if Change Threshold has been reached
    change = (abs(RSSreg - RSSprev));
    
    if(parameters.produceDisplay)
        disp(' ');
        disp(strcat('Change in Objective Function Value: ', num2str(change)))
        disp(strcat('Minimum of Maximum Proportions: ', num2str(minmaxP)))
        disp(strcat('Number of Endmembers: ', num2str(M)))
        disp(strcat('Iteration: ', num2str(iteration)))
        disp(' ');
    end
    
end
