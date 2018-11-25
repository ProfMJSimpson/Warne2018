%% Demonstration of multiple stochastic realisations
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

% generate N realisations
N = 4;
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    [X_r{i},t_r{i}] = GillespieDirectMethod(michment,80);
end

% plot samples with trasparant overlay (hint: for large N, make alpha smaller)
alpha = 0.5;
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    ha = plot(tn,Xn(1,:),'b','LineWidth',2); ha.Color(4) = alpha; 
    hb = plot(tn,Xn(2,:),'r','LineWidth',2); hb.Color(4) = alpha;
    hc = plot(tn,Xn(3,:),'Color',[237,177,32]/255,'LineWidth',2); hc.Color(4) = alpha; 
    hd = plot(tn,Xn(4,:),'Color',[127,47,142]/255,'LineWidth',2); hd.Color(4) = alpha;
end
xlim([0,80]); ylim([0,100]); legend({'S','E','C','P'});
xlabel('time (sec)'); ylabel('copy numbers (molecules)');
