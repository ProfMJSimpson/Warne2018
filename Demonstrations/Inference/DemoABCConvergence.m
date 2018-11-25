%% Demonstration of approximate Bayesian computation convergence
% as epsilon -> 0
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
rng(513,'twister');

% generate data from discrete sampling of a single realisation, 
% no observation error
k_true = [1.0;0.1;0.05]; 
X0 = [100;0];
t = [0;25;50;75;100];
[monomol] = MonoMolecularChain(k_true,X0);

% assume no observation error for now.
Y_obs = GenerateObservations(monomol,k_true,X0,1,[1;2],t,0);

% discrepancy thresholds to demonstrate convergence
epsilon = [100,50,25,20,15];

% discrepancy function as a function of simulated data
rho = @(X_s) sqrt(sum((X_s(:) - Y_obs(:)).^2));

% create uniform joint prior
kmax = [2;0.2;0.1];
kmin = [eps;eps;eps];
p = @() unifrnd(kmin,kmax);

% Simulation as a function of k only
f = @(k) GenerateObservations(monomol,k,X0,1,[1;2],t,0);

%% collect ABC samples at each discrepancy threshold
theta_epsilon  = cell(length(epsilon),1);
C_epsilon = zeros(size(epsilon));
N = 100; % small sample number of performance sake
for j = 1:length(epsilon)
    tic;
    theta_epsilon{j} = ABCRejectionSampler(N,p,f,rho,epsilon(j));
    C_epsilon(j) = toc
end

%% now compute with exact likelihood rejection samples for comparison
theta_exact = [];
c_exact = 0;
c = 9e-15; % We a priori know the upper bound on the likelihood, 
           % it is necessary to increase acceptance rates enough to 
           % be feasible.
tic;
while size(theta_exact,2) < N
    % generate trial from prior
    k_trial = p();
    %compute acceptance probability based on the likelihood
    L = 1;
    for i = 2:size(Y_obs,2)
        L = L*CMEsolMonoMol(k_trial,Y_obs(1,i),Y_obs(2,i),t(i) - t(i-1),...
                            Y_obs(1,i-1),Y_obs(2,i-1));
    end
    % accept with probability L/c
    u = unifrnd(0,1);
    if u <= L/c
        theta_exact = [theta_exact,k_trial];
    end
end
C_exact = toc;

%% plot results
for i=1:3
    figure;
    hold on;
    labels = cell(length(epsilon)+1,1);
    for j=1:length(epsilon)
        ksdensity(theta_epsilon{j}(i,:),'Support','positive','BoundaryCorrection','reflection'); 
        xlabel(['k_',num2str(i)]);ylabel(['p_\epsilon (k_',num2str(i),' | Y_{obs})']);
        labels{j} = ['\epsilon = ',num2str(epsilon(j))];
    end
    labels{end} = '\epsilon = 0'; 
    [F,x] = ksdensity(theta_exact(i,:),'Support','positive','BoundaryCorrection','reflection');
    plot(x,F,'--k');
    legend(labels);
end

