function [theta,W] = ABCSMCSampler(N,p,p_pdf,K,K_pdf,f,rho,epsilon)
%% Sequential Monte Carlo Sampler for approximate Bayesian computaion
%
% Inputs:
%    N - the number of particles
%    p - prior distribution sampler
%    p_pdf - prior PDF 
%    K - proposal kernel process, generates new trials
%    K_pdf - the PDF of the proposal kernel
%    f - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - a sequence of discrepancy acceptance thresholds
%
% Outputs:
%    theta - the sequence of particles
%    W     - the sequence of weights
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise
T = length(epsilon);
theta0 = p();
theta = zeros(length(theta0),N,T);
W = zeros(N,T);
[theta(:,:,1)] = ABCRejectionSampler(N,p,f,rho,epsilon(1));
W(:,1) = 1/N;
% sequential sampling
for t=2:T
    % generate new generation of particles
    for i=1:N
        % rejections steps
        r = inf;
        while r > epsilon(t)
            % sample from weighted particles
            j = randsample(N,1,true,W(:,t-1));
            % generate new particle based on proposal kernel
            theta(:,i,t) = K(theta(:,j,t-1));
            theta(theta(:,i,t) < 0,i,t) = 0;
            D_s = f(theta(:,i,t));
            r = rho(D_s);
        end
        % recompute particle weight using optimal backward kernel
        back_K = 0;
        for j=1:N
            back_K = back_K + W(j,t-1)*K_pdf(theta(:,i,t),theta(:,j,t-1));
        end
        W(i,t) = p_pdf(theta(:,i,t))/back_K;
    end
    % resample 
    theta_rs  = theta(:,:,t);
    J = randsample(N,N,true,W(:,t));
    theta(:,:,t) = theta_rs(:,J);
    % re-set weights
    W(:,t) = 1/N;
end
