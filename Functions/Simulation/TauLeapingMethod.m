function [Z,t] = TauLeapingMethod(bcrn,T,tau)
%% The Tau Leaping Method
% An approximate stochastic simulation algorithm
%
% Inputs:
%    bcrn - a biochemical reaction network struct
%    T    - the end time of the simulation
%    tau   - the timestep
% Outputs:
%    Z    -  time series of copy number vectors
%    t    -  vector of times
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise
Nt = floor(T/tau);
Z = zeros(length(bcrn.X0),Nt+1);
t = zeros(1,Nt+1);
Z(:,1) = bcrn.X0;

for i=1:Nt
    % compute propensities
    a = bcrn.a(Z(:,i),bcrn.k); a(a < 0 ) = 0;
    % generate poisson variates
    Y = poissrnd(a*tau);
    % update copy numbers
    Z(:,i+1) = Z(:,i) + (bcrn.nu') * Y;
    t(i+1) = t(i) + tau;
end
