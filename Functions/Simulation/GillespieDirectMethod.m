function [X,t] = GillespieDirectMethod(bcrn,T)
%% Gillespie Direct Method
% An exact stochastic simulation algorithm
% 
% Inputs:
%    bcrn - a biochemical reaction network struct
%    T    - the end time of the simulation
% Outputs:
%    X    -  time series of copy number vectors
%    t    -  vector of reaction times
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise
X = [bcrn.X0];
t = [0];
while true
    % compute propensities
    a = bcrn.a(X(:,end),bcrn.k);
    % sample exponential waiting time
    dt = exprnd(1/sum(a));
    % check if the simulation is finished
    if t(end) + dt <= T
        % sample the next reaction event
        j = randsample(bcrn.M,1,true,a);
        % update copy numbers  
        X = [X,X(:,end) + bcrn.nu(j,:)'];
        t = [t, t(end) + dt];
    else
        return;
    end
end
