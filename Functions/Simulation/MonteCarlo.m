function [EX,VarX] = MonteCarlo(bcrn,N,T)
%% Monte Carlo
% Computes unbiased estimates of the expectation and Variance
% of the given biochemical reaction network at time T
%
% Inputs:
%    bcrn - a biochemical reaction network struct
%    N    - number of simulations to compute the estimates
%    T    - the time to estimate E[X(T)] and Var[X(T)]
% Outputs:
%   EX    - estimate of E[X(T)]
%   VarX  - estimate of Var[X(T)]
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% to store realisations
R = zeros(length(bcrn.X0),N);

% Generate stochastic samples with exact method
for i=1:N
    [X,t] = GillespieDirectMethod(bcrn,T);
    R(:,i) = X(:,end); % only need the final timestep
end

EX = sum(R,2)/N;
EX2 = sum(R.^2,2)/N;
% unbiased estimate of variance
VarX = (N/(N-1))*(EX2 - EX.^2);
