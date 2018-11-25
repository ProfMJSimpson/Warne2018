%% Demonstration of the key processes in ABC methods
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
T = t(end);
Nt = length(t);
[monomol] = MonoMolecularChain(k_true,X0);
% assume no observation error for now.
Y_obs = GenerateObservations(monomol,k_true,X0,1,[1;2],t,0);
rng(513,'twister');
[X_r,t_r] = GillespieDirectMethod(monomol,T);
figure
hold on;
    Xn = reshape([X_r;X_r],size(X_r).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r;t_r],[1,2*length(t_r)]); tn(1) = [];
    ha = plot(tn,Xn(2,:),'r','LineWidth',2); ha.Color(4) = 0.5; 

plot(t,Y_obs(2,:),'+k','LineWidth',2);
xlabel('time (sec)');
ylabel('copy numbers (molecules)');
% discrepancy threshold
epsilon = [15]

% discrepancy function 
rho = @(X_s) sqrt(sum((X_s(:) - Y_obs(:)).^2));
prior = [];
posterior = [];
klist=[0.03,0.04,0.05,0.06,0.07]
for k3=klist
    % a prior sample with low posterior probability
    theta_low = [1.0;0.1;k3];
    %theta_low = [0.5;0.08;0.03];
    % a prior sample with high posterior probability
    theta_high = k_true
    
    N = 1000;
    
    %% generate prior samples  for low probability
    monomol.k = theta_low
    X_r = cell(N,1); t_r = cell(N,1);
    rho_val = zeros(N,1);
    for i=1:N
        [X_r{i},t_r{i}] = GillespieDirectMethod(monomol,T);
        % get simulated data from BCRN sample path
        X_s = zeros(size(Y_obs));
        for j=1:Nt
            [J] = find(t_r{i} <= t(j));
            X_s(:,j) = X_r{i}(:,J(end)); 
        end
        rho_val(i) = rho(X_s);
        if rho_val(i) <= epsilon
        X_s - Y_obs
        rho_val(i)
        end
    end
    
    [~,I] = sort(rho_val, 'descend');
    
    colourmap = [102,194,165; % [rho > 100]
                 252,141,98;  % [100 >= rho > 50]
                 141,160,203; % [50 >= rho > 25]
                 231,138,195; % [25 >= rho > 20]
                 166,216,84;  % [ 20 >= rho > 15]
                 255,217,47]/255; % [15 >= rho]
    
    figure;
    hold on;
    for i=1:length(I)
        prior = [prior,k3];
        if rho_val(I(i)) <= epsilon
            col = colourmap(5,:);
            alpha = 1.0;
            posterior = [posterior,k3];
        else
            col = [0.5,0.5,0.5];
            alpha = 0.05;
        end
        % plot in colour
        Xn = reshape([X_r{I(i)};X_r{I(i)}],size(X_r{I(i)}).*[1,2]); Xn(:,end) = [];
        tn = reshape([t_r{I(i)};t_r{I(i)}],[1,2*length(t_r{I(i)})]); tn(1) = [];
        ha = plot(tn,Xn(2,:),'Color',col,'LineWidth',2); ha.Color(4) = alpha; 
        ylim([0,100]);
    end
    % plot data
    plot(t,Y_obs(2,:),'+k','LineWidth',2);
    xlabel('time (sec)');
    ylabel('copy numbers (molecules)');
    
    end
edges = klist - (klist(2)-klist(1))/2
edges = [edges,(klist(end) + (klist(2)-klist(1))/2)]
figure;
histogram(prior,edges);
ylabel('Prior sample count')
xlabel('k_3');
figure;
histogram(posterior,edges);
xlabel('k_3');
ylabel('Posterior sample count')

%% now do bivariate ABC posterior

% create uniform joint prior
kmax = [2;0.2;0.1];
kmin = [eps;eps;eps];
p = @() unifrnd(kmin,kmax);
% Simulation as a function of k only
f = @(k) GenerateObservations(monomol,k,X0,1,[1;2],t,0);

N = 100; % small sample number of performance sake
[theta_epsilon,theta_prior] = ABCRejectionSamplerV2(N,p,f,rho,epsilon);
%% plot bivariate histograms
figure;
subplot(3,1,1);
histogram2(theta_epsilon(1,:)',theta_epsilon(2,:)',[10,10]);
zlabel('Posterior sample count')
xlabel('k_1');
ylabel('k_2');
subplot(3,1,2);
histogram2(theta_epsilon(2,:)',theta_epsilon(3,:)',[10,10]);
zlabel('Posterior sample count')
xlabel('k_2');
ylabel('k_3');
subplot(3,1,3);
histogram2(theta_epsilon(3,:)',theta_epsilon(1,:)',[10,10]);
zlabel('Posterior sample count')
xlabel('k_3');
ylabel('k_1');
figure;
subplot(3,1,1);
histogram2(theta_prior(1,:)',theta_prior(2,:)',[10,10]);
xlabel('k_1');
ylabel('k_2');
zlabel('Prior sample count')

subplot(3,1,2);
histogram2(theta_prior(2,:)',theta_prior(3,:)',[10,10]);
xlabel('k_2');
ylabel('k_3');
zlabel('Prior sample count')

subplot(3,1,3);
histogram2(theta_prior(3,:)',theta_prior(1,:)',[10,10]);
xlabel('k_3');
ylabel('k_1');
zlabel('Prior sample count')
