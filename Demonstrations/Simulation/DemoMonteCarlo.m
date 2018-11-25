%% Comparison of Monte Carlo approaches
% for the forwards problem
% Task to compute E[B(T)] for the monomolecular chain
% with 95% confidence interval of +-h
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

T = 100;
% Build mono-molecular chain
[monomol] = MonoMolecularChain([10.0;0.1;0.05],[1000;0]);

% The target error  
h = [1,0.5,0.25,0.125,0.06125];

% MLMC parameters
L = [1,2,3,3,4]; M = 2;
tau0 = [1,1,1,0.5,0.5]

%% Monte Carlo using MLMC method
fprintf('Performance test, Multilevel Monte Carlo  ...\n')
C_mlmc = zeros(size(h))';
E_mlmc = zeros(size(h))';
CI_mlmc = zeros(size(h))';
for i=1:length(h)
    % determine Nl using 100 trial samples at each level
    tic;
    V_test = zeros(2,L(i)+1);
    C_test = zeros(1,L(i)+1);
    [~,V_test(:,1)] = MonteCarloTauLeap(monomol,100,T,tau0(i));
    C_test(1) = toc;
    for l=1:L(i)
        tic;
        [~,V_test(:,l+1)] = MonteCarloBiasCorrection(monomol,100,T,tau0(i)*M^-l,M);
        C_test(l+1) = toc;
    end
    % optimally choose Nl for 95% CI +-h 
    % (based on Lagrange Multiplier)
    c = (1.95^2)*(h(i)^-2)*sum(sqrt(V_test(2,:).*C_test),2);
    Nl = ceil(c*sqrt(V_test(2,:) ./ C_test))
    % compute the estimate using MLMC
    tic;
    [EZ,VarZ] = MultilevelMonteCarlo(monomol,Nl,L(i),M,tau0(i),T);
    C_mlmc(i) = toc + sum(C_test)
    % the estimate   C_mc 
    E_mlmc(i) = EZ(2);
    CI_mlmc(i) = (1.95)*sqrt(VarZ(2));
end

%% Monte Carlo using approximate Tau-Leaping Method
fprintf('Performance test, Monte Carlo with Tau-Leaping ...\n')
C_tau = zeros(size(h))';
E_tau = zeros(size(h))';
CI_tau = zeros(size(h))';
for i=1:length(h)
    % choose tau to match MLMC bias
    tau = tau0(i)*M^-L(i) 
    % determine N using 100 trial samples
    tic;
    [~,VarZ] = MonteCarloTauLeap(monomol,100,T,tau);
    N = ceil((1.95^2)*VarZ(2)/((h(i))^2))
    % compute the estimate using standard Monte Carlo 
    [EZ,VarZ] = MonteCarloTauLeap(monomol,N,T,tau);
    C_tau(i) = toc
    % the estimate   
    E_tau(i) = EZ(2);
    CI_tau(i) = (1.95)*sqrt(VarZ(2)/N);
end

%% Monte Carlo using exact Gillespie method
fprintf('Performance test, Monte Carlo with Gillespie ...\n')
C_mc = zeros(size(h))';
E_mc = zeros(size(h))';
CI_mc = zeros(size(h))';
for i=1:length(h)
    % determine N using 100 trial samples
    tic;
    [~,VarX] = MonteCarlo(monomol,100,T);
    N = ceil((1.95^2)*VarX(2)/(h(i)^2))
    % compute the estimate using standard Monte Carlo 
    [EX,VarX] = MonteCarlo(monomol,N,T);
    C_mc(i) = toc
    % the estimate   
    E_mc(i) = EX(2);
    CI_mc(i) = (1.95)*sqrt(VarX(2)/N);
end



%% plot results
loglog(CI_mc,C_mc,'ob--','LineWidth',2);
hold on
loglog(CI_tau,C_tau,'^r--','LineWidth',2);
loglog(CI_mlmc,C_mlmc,'sk--','LineWidth',2);
grid on;
ylabel('Compute Time (sec)');
xlabel('Error Tolerance');
legend({'Monte Carlo (GDM)','Monte Carlo (TLM)','Multilevel Monte Carlo'});

