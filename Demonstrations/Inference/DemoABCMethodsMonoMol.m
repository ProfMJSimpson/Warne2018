%% Demonstration of Monte Carlo methods for
%  approximate Bayesian computation 
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

%
% discrepancy function as a function of simulated data
rho = @(X_s) sqrt(sum((X_s(:) - Y_obs(:)).^2));

% Simulation as a function of k only
f = @(k) GenerateObservations(monomol,k,X0,1,[1;2],t,0);
% prior support (uniform)
kmax = [2;0.2;0.1];
kmin = [eps;eps;eps];

%% Set up ABC MLMC 
% discrepancy threshold sequence
epsilon = [100;50;25;20;15];
% create uniform joint prior
supp0.l = kmin;
supp0.u = kmax;
p = @(l,u) unifrnd(l,u);
% sequence of sample numbers
N = [800;400;200;100;50];

%% Run and time ABC MLMC
fprintf('Running ABC MLMC...\n');
tic;
[E_mlmc,V_mlmc,F_mlmc] = ABCMLMC(N,p,supp0,f,rho,epsilon)
C_mlmc = toc;
fprintf('ABC MLMC Completed in %f Sec\n',C_mlmc);

%% Set up ABC Rejection
% discrepancy threshold
epsilon = 15;
% Prior sampler
p = @() unifrnd(kmin,kmax);
% number of samples
N = 100;

%% Run and time ABC Rejection
fprintf('Running ABC Rejection...\n');
tic;
theta_rej = ABCRejectionSampler(N,p,f,rho,epsilon);
E_rej = mean(theta_rej,2);
V_rej = (1/(N-1))*(mean(theta_rej.^2,2) - E_rej.^2);
C_rej = toc;
fprintf('ABC Rejection Completed in %f Sec\n',C_rej);

%% Set up ABC MCMC 
% discrepancy threshold
epsilon = 15;
% Prior sampler and PDF
p = @() unifrnd(kmin,kmax);
p_pdf = @(k) prod(unifpdf(k,kmin,kmax));
% create proposal Kernel sampler and PDF
Sigma = diag((0.05*(kmax-kmin)).^2);
K = @(k) mvnrnd(k,Sigma)';
K_pdf = @(k_n,k_p) mvnpdf(k_n,k_p,Sigma);
% number of steps in the Markov Chain
T = 500000;
burnin = 100000;
thin = 10000;

%% Run and ABC MCMC
tic;
fprintf('Running ABC MCMC...\n');
theta0 = ABCRejectionSampler(1,p,f,rho,epsilon);
theta_mcmc = ABCMCMCSampler(T,p_pdf,K,K_pdf,f,rho,epsilon,theta0);
E_mcmc = mean(theta_mcmc(:,burnin:thin:T),2);
V_mcmc = (1/(N-1))*(mean(theta_mcmc(:,burnin:thin:T).^2,2) - E_mcmc.^2);
C_mcmc = toc;
fprintf('ABC MCMC Completed in %f Sec\n',C_mcmc);

%% Set up ABC SMC
% discrepancy threshold sequence
epsilon = [100;50;25;20;15];
% Prior sampler and PDF
p = @() unifrnd(kmin,kmax);
p_pdf = @(k) prod(unifpdf(k,kmin,kmax));
% create proposal Kernel sampler and PDF
Sigma = diag((0.05*(kmax-kmin)).^2);
K = @(k) mvnrnd(k,Sigma)';
K_pdf = @(k_n,k_p) mvnpdf(k_n,k_p,Sigma);
% number of particles
N = 100;

%% Run and time ABC SMC
fprintf('Running ABC SMC...\n');
tic;
[theta_smc,W] = ABCSMCSampler(N,p,p_pdf,K,K_pdf,f,rho,epsilon);
E_smc = theta_smc(:,:,end)*W(:,end);
V_smc = (1/(N-1))*(mean(theta_smc(:,:,end).^2,2) - E_smc.^2);
C_smc = toc;
fprintf('ABC SMC Completed in %f Sec\n',C_smc);

%% collect results
comp = [E_mlmc',V_mlmc',C_mlmc;
        E_rej',V_rej',C_rej;
        E_mcmc',V_mcmc',C_mcmc;
        E_smc',V_smc',C_smc]
comp(:,4:6) = 1.95*sqrt(comp(:,4:6));

%% plot Markov Chain transient behaviour
for i=1:3
    figure;
    plot(theta_mcmc(i,:),'b','LineWidth',2);xlabel('t');ylabel(['k_',num2str(i)])
    hold on;
    plot([0,length(theta_mcmc(i,:))],[k_true(i),k_true(i)],'--k','LineWidth',2)
end

%% plot marginal densities 
ks = [linspace(0,2,1000);linspace(0,0.2,1000);linspace(0,0.1,1000)];
for i=1:3
    figure;
    hold on;
    plot(ks(i,2:end),diff(F_mlmc{i}(ks(i,:)))./diff(ks(i,:)));
    ksdensity(theta_rej(i,:),'Support','positive','BoundaryCorrection','reflection'); 
    ksdensity(theta_mcmc(i,burnin:thin:T),'Support','positive','BoundaryCorrection','reflection');
    ksdensity(theta_smc(i,:,end),'Support','positive','BoundaryCorrection','reflection'); 
    ylabel(['p_\epsilon (k_',num2str(i),')']);xlabel(['k_',num2str(i)])
    xlim([ks(i,1),ks(i,end)]);
end
