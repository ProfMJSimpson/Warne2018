function [theta] = ABCMCMCSampler(T,p_pdf,K,K_pdf,f,rho,epsilon,theta0)
%% Markov Chain Monte Carlo Sampler for approximate Bayesian computaion
%
% Inputs:
%    T - the number of steps to simulate the Markov Chain
%    p_pdf - the PDF function of the parameter joint prior distribution
%    K - proposal kernel process, generates new trial
%    K_pdf - the PDF function of the proposal kernel
%    f - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - the discrepancy acceptance threshold
%    theta0 - initial condition of the Markov Chain
%
% Outputs:
%    theta - The Markov Chain trajectory 
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise the chain
theta = zeros(length(theta0),T);
theta(:,1) = theta0;

for i=2:T
    % generate new proposal based on previous state
    theta_trial = K(theta(:,i-1));
    theta_trial(theta_trial < 0) = 0;
    % generate simulated data using these parameters
    D_s = f(theta_trial);
    % accept/reject proposal
    if rho(D_s) <= epsilon
        % compute Metropolis-Hastings acceptance probability
        h = min(1,(p_pdf(theta_trial)/p_pdf(theta(:,i-1))) * ...
            (K_pdf(theta(:,i),theta_trial)/K_pdf(theta_trial,theta(:,i))));
        % accept/reject transition    
        u = unifrnd(0,1);
        if u <= h
            % accept proposed move
            theta(:,i) = theta_trial;
        else
            % stay at current state
            theta(:,i) = theta(:,i-1);
        end
    else
        % stay at current state
        theta(:,i) = theta(:,i-1);
    end
end
