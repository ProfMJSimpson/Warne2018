%% Demonstration of correlated Tau-Leap 
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;

% Build mono-molecular chain
[monomol] = MonoMolecularChain([1.0;0.1;0.05],[100;0]);

% generate N realisations with same RNG seed to ensure
% the tau = 1 sample path is the same
rng(502,'twister');
[Zf_r,Zc_r2,t_r] = CorTauLeapingMethod(monomol,100,1,2)
rng(502,'twister');
[Zf_r,Zc_r4,t_r] = CorTauLeapingMethod(monomol,100,1,4)
rng(502,'twister');
[Zf_r,Zc_r8,t_r] = CorTauLeapingMethod(monomol,100,1,8)
hold on;
% plot correlated sample paths to visualise the strong convergence
Zfn = reshape([Zf_r;Zf_r],size(Zf_r).*[1,2]); Zfn(:,end) = [];
Zcn2 = reshape([Zc_r2;Zc_r2],size(Zc_r2).*[1,2]); Zcn2(:,end) = [];
Zcn4 = reshape([Zc_r4;Zc_r4],size(Zc_r4).*[1,2]); Zcn4(:,end) = [];
Zcn8 = reshape([Zc_r8;Zc_r8],size(Zc_r8).*[1,2]); Zcn8(:,end) = [];
tn = reshape([t_r;t_r],[1,2*length(t_r)]); tn(1) = [];
plot(tn,Zfn(2,:),'-','LineWidth',2); 
plot(tn,Zcn2(2,:),'--','LineWidth',2); 
plot(tn,Zcn4(2,:),'-.','LineWidth',2); 
plot(tn,Zcn8(2,:),':','LineWidth',2); 
xlim([0,100]); ylim([0,80]); legend({'\tau = 1','\tau = 2','\tau = 4','\tau = 8'});
xlabel('t (sec)'); ylabel('Copy Numbers (Molecules)');

