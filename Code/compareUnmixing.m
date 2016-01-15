function EMDtot = compareUnmixing(P1, P2, E1, E2, gdtype, buildFlag)
%%
% EMDtot = compareUnmixing(P1, P2, E1, E2, gdtype, buildFlag)
%
% Inputs:
%   P1 - NxM matrix - each row is the proportion vector for one data point
%        given endmembers E1
%   P1 - NxL matrix - each row is the proportion vector for one data point
%        given endmembers E2
%   E1 - DxM matrix - each column is an endmember of dimension D
%   E2 - DxL matrix - each column is an endmember of dimension D
%   gdtype - string - Indicates the ground distance to use to compare the
%        endmember sets.  Options are: 'seuclidean' (squared Euclidean
%        distance), 'euclidean' (Euclidean distance), 'sad' (Spectral Angle
%        Distance), or 'sid' (Spectral Information Divergence).
%   buildFlag - Set to 1 if mex files need to be compiled, 0 otherwise.
%
% Outputs:
%    EMDtot - Aggregated  EMD value (Aggregation is the sum in this
%    implementation)
%
% Note: This code uses MEX files.  Your MEX compiler will need to be set up accordingly to be able to run this code.  The first time your code is run, set the buildFlag parameter to 1.  After the initial compilation, this buildFlag can be set to 0. 
%
%%
% Dependencies: This code relies on the pdist function from the MATLAB
% statistics toolbox.  This code also uses the EMD implementation by Y. Rubner found at:
%   http://ai.stanford.edu/~rubner/emd/default.htm
% This code also uses the MATLAB mex interface that can be found at:
%   By M. Alipour
%   http://www.mathworks.com/matlabcentral/fileexchange/12936-emd-earth-movers-distance-mex-interface%
%
%
%%
%  If this code is used in any publication or presentation, the following
%  reference must be included:
%
% A. Zare, D.T. Anderson, "Earth Movers Distance-Based Simultaneous Comparison of Hyperspectral 
% Endmembers and Proportions," IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, In Press.
%
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

if(buildFlag)
    mex  emd_mex.c emd.c   -O  -v
end

if(strcmp(gdtype, 'seuclidean'))
    %Compute squared Euclidean ground distance between every pair of endmembers
    GD = bsxfun(@plus,dot(E1,E1,1)',dot(E2,E2,1))-2*(E1'*E2);
    
elseif(strcmp(gdtype, 'euclidean'))
    %Compute  Euclidean ground distance between every pair of endmembers
    GD = sqrt(bsxfun(@plus,dot(E1,E1,1)',dot(E2,E2,1))-2*(E1'*E2));
    
elseif(strcmp(gdtype, 'sad'))
    %Compute Spectral Angle Distance between every pair of endmembers 
    N1 = sqrt(sum(E1.*E1));
    N2 = sqrt(sum(E2.*E2));
    N = N1'*N2;
    GD= acos(round(((E1'*E2)./N)*1000)/1000);
    
elseif(strcmp(gdtype, 'sid'))
    %Compute Spectral Information Divergence between every pair of
    %endmembers; Endmembers are assumed to be appropriately normalized for
    %the SID computation (non-negative and sum-to-one across wavelengths)
    for i = 1:size(E1,2)
        for j = 1:size(E2,2)
            GD(i,j) = sum(E1(:,i).*log(E1(:,i)./(E2(:,j)+eps)+eps)) + sum(E2(:,j).*log(E2(:,j)./(E1(:,i)+eps)+eps));
        end
    end
else
    error('Incorrect gdtype error: gdtype should be set to either seuclidean, euclidean, sad or sid');
end

if(sum(sum(isnan(GD))) || sum(sum(isinf(GD))) || sum(sum(~isreal(GD))))
    error('Error in ground distance computation. Please use another ground distance or correct/normalize the endmembers appropriately');
end

%Compute each emd for each data point
[N] = size(P1,1);
for i = 1:N
    emd(i) = emd_mex(P1(i,:), P2(i,:), GD);
end

%Compute total EMD value 
EMDtot = sum(emd);
end