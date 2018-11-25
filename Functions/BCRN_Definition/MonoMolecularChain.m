function [bcrn] = MonoMolecularChain(k,X0)
%% MonoMolecularChain 
% Construction of a mono-molecular chain bio-chemical reaction network.
%   k(1)    k(2)           k(N)    k(N+1)
% 0 -> X(1) -> X(2) -> ... -> X(N) -> 0
%
% Inputs:
%    k - vector of kinetic rate parameters
%    X0 - initial state vector
%
% Outputs:
%     a BCRN struct
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

bcrn = struct();
% kinetic rate parameters
bcrn.k = k;              
% number of reactions
bcrn.M = length(k);                   
% number of chemical species
bcrn.N = bcrn.M - 1;                       
% reactant stoichiometries
bcrn.nu_minus = [zeros(1,bcrn.N);diag([ones(1,bcrn.N)])];  
% product stoichiometries
bcrn.nu_plus = [diag([ones(1,bcrn.N)]);zeros(1,bcrn.N)];         
% stoichiometric matrix
bcrn.nu = bcrn.nu_plus - bcrn.nu_minus;         
% initial copy numbers
bcrn.X0 = X0;                
% propensity function
bcrn.a = @(X,k) k.*[1;X];                 
                                    
