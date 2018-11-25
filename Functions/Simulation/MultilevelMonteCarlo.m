function [EZ,VarZ] = MultilevelMonteCarlo(bcrn,Nl,L,M,tau0,T)
%% Multilevel Monte Carlo
% Computes biased estimates of the expectation of a given biochemical 
% reaction network at time T using Multilevel Monte Carlo (MLMC) methods. 
%
% Inputs:
%    bcrn - a biochemical reaction network
%    Nl   - L+1 vector of simulation numbers for each level
%    L    - number of bias correction levels
%    M    - Tau-leaping step scale factor between levels
%    tau0 - Tau-leaping step size of the coarsest estimator
%    T    - the time to estimate E[Z_L(T)]
% Outputs:
%    EZ   - MLMC estimate of E[Z_L(T)]
%    VarZ - MLMC estimator variance, Var[EZ]
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% compute coarse estimates E[Z_0(T)]
[EZ0,VarZ0] = MonteCarloTauLeap(bcrn,Nl(1),T,tau0);

% compute bias corrections, estimate E[Z_l(T) - Z_{l-1}(T)]
BC = zeros(length(bcrn.X0),L);
VarBC = zeros(length(bcrn.X0),L);
for l=1:L
    tau_f = tau0*M^-l;
    [BC(:,l),VarBC(:,l)] = MonteCarloBiasCorrection(bcrn,Nl(l+1),T,tau_f,M);
end
% compute MLMC telescoping sum for final estimator
EZ = EZ0 + sum(BC,2);
VarZ = sum([VarZ0,VarBC] .* (repmat(Nl.^-1,[length(bcrn.X0),1])),2);

