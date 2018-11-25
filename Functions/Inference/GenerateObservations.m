function [Y_obs] = GenerateObservations(bcrn,k,X0,N,I,t,sigma)
%% Generate Observations
% Generates N Observations of state dimensions I of stochastic process X(t)
% with rate parameters k and initial condition X(0) = X0 = X0. Observations are 
% taken at time points t(1),..,t(Nt) with observation error 
% Y_obs(t(i)) = X(t(i)) + zeta, where zeta ~ N(0,sigma^2)
%
% Inputs:
%    bcrn - a BCRN structure
%    k    - true kinetic rate parameters for the model
%    X0   - true initial condition of the BCRN
%    N    - number of sample paths in dataset
%    I    - vector of observable chemical species indices
%    t    - vector of observation times
%    sigma - standard deviation of observation error
% 
% Outputs:
%     a Table of observations Y_obs
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

Nt = length(t);
Y_obs = zeros(length(I),Nt,N);
bcrn.k = k;
for i=1:N
    % stochastic simulation
    [X_r,t_r] = GillespieDirectMethod(bcrn,t(Nt));
    % just store the state of selected species and times
    for j=1:Nt
        [J] = find(t_r <= t(j));
        Y_obs(:,j,i) = X_r(I,J(end));
    end
end
% perturb with noise
zeta = normrnd(0,sigma,size(Y_obs));
Y_obs = Y_obs + zeta;
