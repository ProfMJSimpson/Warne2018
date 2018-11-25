%% Demonstration of Direct application of Bayes' theorem using
%  the solution to the Chemical Master Equation. 
%
% NOTE: This approach is not accurate without a very find grid
%       to enable computation of the evidence.
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
rng(513,'twister');

% generate data from discrete sampling of a single realisation, 
% no observation error
[monomol] = MonoMolecularChain([1.0;0.1;0.05],[100;0]);
t = [0;25;50;75;100];
Y_obs = GenerateObservations(monomol,1,[1;2],t,0);

% create uniform joint prior
kmax = [2;0.2;0.1];
kmin = [eps;eps;eps];
p = @(k) prod(unifpdf(k,kmin,kmax));

% build a grid over the prior support
d = 20;
ks = [linspace(kmin(1),kmax(1),d);
      linspace(kmin(2),kmax(2),d);
      linspace(kmin(3),kmax(3),d)]
[K1,K2,K3] = ndgrid(ks(1,:),ks(2,:),ks(3,:))
dk = (kmax - kmin)/(d-1);

% Compute Posterior
PK = zeros(size(K1));
PY_obs = 0;
for i = 1:(d^3)
    % compute likelihood for this parameter set using CME solution
    likelihood = 1;
    for j = 2:size(Y_obs,2)
        Pt = CMEsolMonoMol([K1(i);K2(i);K3(i)],...
                                              Y_obs(1,j),Y_obs(2,j),t(j)-t(j-1), ...
                                              Y_obs(1,j-1),Y_obs(2,j-1))
        likelihood = likelihood*Pt;
    end
    % the posterior
    PK(i) = likelihood*p([K1(i);K2(i);K3(i)]);
end

% compute evidence term
PY_obs = trapz(ks(1,:),trapz(ks(2,:),trapz(ks(3,:),PK,3),2),1);

% normalise to obtain the posterior
PK = PK /PY_obs;

% compute Marginals
Pm1 = trapz(ks(2,:),trapz(ks(3,:),PK,3),2);
Pm2 = trapz(ks(1,:),trapz(ks(3,:),PK,3),1);
Pm3 = trapz(ks(1,:),trapz(ks(2,:),PK,2),1);

% plot marginals
figure;
hold on;
plot(ks(1,:),squeeze(Pm1),'--b');
plot(ks(2,:),squeeze(Pm2),'--r');
plot(ks(3,:),squeeze(Pm3),'--k');

