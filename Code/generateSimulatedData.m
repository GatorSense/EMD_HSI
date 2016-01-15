%% function[X,Etrue,Ptrue]=generateSimulatedData()
% Algorithm-I (EMD paper) Simualted 2-D Data Generation
% Output- X(N x D)
%       - Etrue(D x M)
%       - Ptrue(N x M)
%%
%  If this code is used or referenced in any publication or presentation, the following
%  reference must be included:
%
%  A. Zare, D. Anderson, "Earth Movers Distance-Based Simultaneous
%  Comparison of Hyperspectral Endmembers and Proportions," IEEE Journal of
%  Selected Topics in Applied Earth Observations and Remote Sensing, In
%  Press.
%
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
function[X,Etrue,Ptrue]=generateSimulatedData()
    Nx = 500;                                                                  % Number of data points
    Dimensions = 2;                                                            % Number of dimensions
    load 'Etrue.mat';                                                          % Load true endmembers
    noEndmembers = size(Etrue,2);
    X = zeros(Nx,Dimensions);                                                  % Initialize data matrix
    Ptrue=zeros(Nx,noEndmembers);                                              % Initialize proportion matrix
    smallConst=0.000001;
    for i = 1:Nx
        tmpE = zeros(Dimensions,noEndmembers);
        tmpMj = mnrnd(1,[1/noEndmembers,1/noEndmembers,1/noEndmembers]);       % Sample mj from multinomial distribution
        mj = find(tmpMj==1);                                                   % Determine number of endmeber/s in pixel i
        index = randperm(noEndmembers,mj);                                     % Sample mj endmembers from full set of M endmembers
        tmpE(:,index) = Etrue(:,index);                                        

        Alpha    = exprnd(Dimensions,[1 mj]);                                  % Sample alpha from exponential distribution with mean=2
        A        = repmat(Alpha,[1 1]);
        Y        = randg(A); 
        v        = sum(Y,2);
        Ptrue(i,index)    = Y./repmat(v, [1, size(A,2)]);                      % Generate proportion using Dirchlet Distribution
        aa = isnan(Ptrue(i,:));                                                % Avoid NaN's(if any)
        if(find(aa==1))
            Ptrue(i,find(aa==1))=smallConst;
        end
        X(i,:) = Ptrue(i,:)*tmpE';                                             % Generate ith data point
    end
    Ptrue = Ptrue + smallConst;                                                % Add small constant to avoid zero proportion values
    X = X + rand(Nx, size(Etrue,1))*.0001;                                     % Add zero-mean Gaussian random noise
    clearvars -EXCEPT Ptrue Etrue X;
end
% end