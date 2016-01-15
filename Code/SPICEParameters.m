
function [parameters] = SPICEParameters()

%   This function sets the parameters to be used during the SPICE algorithm
%   parameters - struct - The struct contains the following fields:
%                   1. u : Regularization Parameter to trade off between the RSS and V terms
%                   2. gamma: Gamma Constant for SPT term, Controls degree
%                       of sparsity among endmembers
%                   3. changeThresh: Stopping Criteria, Change threshold
%                       for Objective Function.
%                   4. M: Initial Number of endmembers
%                   5. iterationCap: Maximum number of iterations
%                   6. endmemberPruneThreshold: Proportion threshold used
%                       to prune endmembers
%                   7. produceDisplay : Set to 1 if a progress display is
%                       wanted
%                   9. initEM: Set to nan to randomly select endmembers,
%                       otherwise NxM matrix of M endmembers with N spectral
%                       bands, Number of endmembers must equal parameters.M
%
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
% THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF FLORIDA AND
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

parameters.u = 0.001; %Trade-off parameter between RSS and V term
parameters.gamma = 10; %Sparsity parameter
parameters.M = 20; %Initial number of endmembers
parameters.endmemberPruneThreshold = 1e-9;

parameters.changeThresh = 1e-4; %Used as the stopping criterion
parameters.iterationCap = 5000; %Alternate stopping criterion

parameters.produceDisplay = 1;
parameters.initEM = nan; %This randomly selects parameters.M initial endmembers from the input data
parameters.options = optimset('Display', 'off');
