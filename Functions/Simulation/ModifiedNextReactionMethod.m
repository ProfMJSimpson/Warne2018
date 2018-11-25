function [X,t] = ModifiedNextReactionMethod(bcrn,T)
%% Modified Next Reaction Method
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

%initialise
X = [bcrn.X0];
t = [0];
T_r = zeros(bcrn.M,1);
% generate M unit-time exponential variates
P = exprnd(1,[bcrn.M,1]);
while true
    % compute propensities
    a = bcrn.a(X(:,end),bcrn.k);
    % determine which reaction channel fires next
    dt = (P - T_r) ./ a;
    dt(a <= 0) = Inf;
    [delta,mu] = min(dt);
    if t(end) + delta <= T
        %update copy numbers
        X = [X,X(:,end) + bcrn.nu(mu,:)'];
        t = [t,t(end) + delta];
        T_r = T_r + a*delta;
        % update next reaction time for the firing channel
        P(mu) = P(mu) + exprnd(1);
    else
        return;
    end
end

