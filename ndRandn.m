function [Y] = ndRandn(mu, CY, num)

% mu = desired mean of sample
% CX = desired covariance of sample
% num = # of samples (num x N)

% generates a set of
% samples drawn from an N-dimensional Gaussian distribution with the specified mean
% (an N-vector) and covariance (an NxN matrix). The parameter num should be optional
% (defaulting to 1) and should specify the number of samples to return. The returned 
% value should be a matrix with num rows each containing a sample of N elements.
if nargin < 2
    disp('Not enough input arguments')
end 
if nargin == 2
    num = 1;
end 

N = length(mu);
X = randn(num,N);

% Y'Y = XM' XM
% M'M = C
% M'M = VrootD'rootDV'
% M = rootDV'

[u,s,v] = svd(CY);
M = sqrt(s)*v';
Y = X*M;

Y = Y+ mu;

end 

