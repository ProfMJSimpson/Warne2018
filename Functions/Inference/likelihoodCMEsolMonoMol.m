function [lh] = likelihoodCMEsolMonoMol(k,Y_obs,t);
%% The likelihood of a realisation of the Monomolecular chain
% based on the explicit analytic solution to the CME.
%
% Inputs:
%    k - given values for kinetic rate parameters 
%    Y_obs - observations of a mono-molecular chain sample path
%    t     - observation times.
%
% Outputs:
%     The likelihood of Y_obs given k, that is p(Y_obs | k)
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology
lh = 1;
for j=2:size(Y_obs,2)
    lh = lh*CMEsolMonoMol(k,Y_obs(1,j),Y_obs(2,j),t,Y_obs(1,j-1),Y_obs,2,j-1);
end
