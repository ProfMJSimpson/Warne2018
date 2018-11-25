function [EBC,VarBC] = MonteCarloBiasCorrection(bcrn,N,T,tau_f,M)
%% Monte Carlo
% Computes estimates of the difference in Tau leaping realisations
% E[Z_f - Z_c] with low variance based on Correlated Tau-leaping pairs
%
% Inputs:
%    bcrn - a biochemical reaction network struct
%    N    - number of simulations to compute the estimates
%    T    - the time to estimate E[Z_f(T) - Z_c(T)] and Var[Z_f(T) - Z_c(T)]
%    tau_f - tau timestep for computing Z_f realisations
%    M     - scale factor for tau_c = M*tau_f used for Z_c
% Outputs:
%   EBC    - estimate of E[Z_f(T) - Z_c(T)]
%   VarBC  - estimate of Var[Z_f(T) - Z_c(T)]
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% to store realisations
R = zeros(length(bcrn.X0),N);

% Generate coupled pairs (Z_f,Z_c) 
for i=1:N
    [Z_f,Z_c,t] = CorTauLeapingMethod(bcrn,T,tau_f,M);
    R(:,i) = Z_f(:,end) - Z_c(:,end);
end

EBC = sum(R,2)/N;
EBC2 = sum(R.^2,2)/N;
% unbiased variance estimator
VarBC = (N/(N-1))*(EBC2 - EBC.^2);
