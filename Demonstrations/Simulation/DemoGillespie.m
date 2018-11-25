%% Demonstration of the Gillespie Direct Method
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

% Build Michaelis-Menten model
[michment] = MichaelisMenten([0.001;0.005;0.01],100,100);
% simulate
[X,t] = GillespieDirectMethod(michment,80);
% Plot
subplot(1,2,1);
Xn = reshape([X;X],size(X).*[1,2]); Xn(:,end) = [];
tn = reshape([t;t],[1,2*length(t)]);tn(1) = []; 
plot(tn,Xn,'LineWidth',2); xlim([0,80]); ylim([0,100]);
xlabel('t (sec)'); ylabel('copy numbers (molecules)');
legend({'S','E','C','P'});

% Build mono-molecular chain
[monomol] = MonoMolecularChain([1.0;0.1;0.05],[100;0]);
%simulate
tic;
[X,t] = GillespieDirectMethod(monomol,100);
toc;
% Plot
subplot(1,2,2);
Xn = reshape([X;X],size(X).*[1,2]); Xn(:,end) = [];
tn = reshape([t;t],[1,2*length(t)]); tn(1) = [];
plot(tn,Xn,'LineWidth',2); xlim([0,100]); ylim([0,100]);
xlabel('t (sec)'); ylabel('copy numbers (molecules)');
legend({'A','B'});
