function [Z_f,Z_c,t] = CorTauLeapingMethod(bcrn,T,tau_f,M)
%% Correlated pairs of Tau Leaping simulations
% Generates tau leaping simulations for nested tau's.
% 
% Inputs:
%    bcrn - a biochemical reaction network struct
%    T    - the end time of the simulation
%    tau_f   - the timestep of the fine-grain simulation
%    M    - The integer scale factor tau_c = tau_f*M
% Outputs:
%    Z_f  -  time series of copy number vectors using tau_f
%    Z_c  -  time series of copy number vectors using tau_c
%    t    -  vector of times
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% obtain coarse-grain timestep
tau_c = M*tau_f;
% initialise
Nt = floor(T/tau_f);

Z_f = zeros(length(bcrn.X0),Nt+1);
Z_c = zeros(length(bcrn.X0),Nt+1);
t = zeros(1,Nt+1);
Z_f(:,1) = bcrn.X0;
Z_c(:,1) = bcrn.X0;

for i=1:Nt
    % compute fine-grain propensities
    a_f = bcrn.a(Z_f(:,i),bcrn.k); a_f(a_f < 0 ) = 0;
    % compute coarse-grain propensities every M steps
    if mod(i-1,M) == 0
        a_c = bcrn.a(Z_c(:,i),bcrn.k); a_c(a_c < 0) = 0;
        Z_tmp = Z_c(:,i);
    end
    % update virtual propensities
    b = [zeros(size(a_c)),a_c,a_f] +  min(a_f,a_c)*[1,-1,-1];
    Y = poissrnd(b*tau_f);
    % update fine-grain copy numbers
    Z_f(:,i+1) = Z_f(:,i) + (bcrn.nu') * (Y(:,1) + Y(:,3));
    t(i+1) = t(i) + tau_f;
    % update intermediate course grain
    Z_tmp = Z_tmp + (bcrn.nu') * (Y(:,1) + Y(:,2));
    % update coarse-grain trajectory every M steps
    if mod(i,M) == 0
        Z_c(:,i+1) = Z_tmp;
    else
        Z_c(:,i+1) = Z_c(:,i);
    end
end
