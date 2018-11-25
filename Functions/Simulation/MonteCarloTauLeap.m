function [EZ,VarZ] = MonteCarloTauLeap(bcrn,N,T,tau)
%% Monte Carlo
% Computes biased estimates of the expectation and Variance
% of the given biochemical reaction network at time T
%
% Inputs:
%    bcrn - a biochemical reaction network struct
%    N    - number of simulations to compute the estimates
%    T    - the time to estimate E[Z(T)] and Var[Z(T)]
% Outputs:
%   EZ    - estimate of E[Z(T)]
%   VarZ  - estimate of Var[Z(T)]
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% to store realisations
R = zeros(length(bcrn.X0),N);

% Generate stochastic samples with approximate method
for i=1:N
    [Z,t] = TauLeapingMethod(bcrn,T,tau);
    R(:,i) = Z(:,end); % only need the final timestep
end

EZ = sum(R,2)/N;
EZ2 = sum(R.^2,2)/N;
% unbiased estimate of variance
VarZ = (N/(N-1))*(EZ2 - EZ.^2);
