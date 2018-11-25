function [theta,prior,pps] = ABCRejectionSamplerV2(N,p,f,rho,epsilon)
%% Rejection Sampler for approximate Bayesian computation.
%
%
% Inputs:
%    N - the number of ABC posterior samples
%    p - function that generates iid samples from the parameter joint
%        prior distribution
%    f - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - the discrepancy acceptance threshold
%
% Outputs:
%    theta - a matrix of ABC posterior samples
%    prior - a matrix of all prior samples
%    pps   - a matrix of prior predicitive samples
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

theta = [];
prior = [];
pps = [];
while size(theta,2) < N
    % generate trial from the prior
    theta_trial = p();
    % generate simulated data using these parameters
    D_s = f(theta_trial);
    % accept/reject
    if rho(D_s) <= epsilon
        theta = [theta,theta_trial];
    end
    prior = [prior,theta_trial];
    pps = [pps;D_s];
end
