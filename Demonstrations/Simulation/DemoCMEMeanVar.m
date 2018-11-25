%% Demonstration of CME derived Means and Variances  for 
% the mono-molecular chain
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

% generate N realisations
N = 100;
X_r = cell(N,1); t_r = cell(N,1);
for i=1:N
    [X_r{i},t_r{i}] = GillespieDirectMethod(monomol,100);
end

% plot samples with trasparant overlay (hint: for large N, make alpha smaller)
alpha = 0.05;
hold on;
for i=1:N
    Xn = reshape([X_r{i};X_r{i}],size(X_r{i}).*[1,2]); Xn(:,end) = [];
    tn = reshape([t_r{i};t_r{i}],[1,2*length(t_r{i})]); tn(1) = [];
    ha = plot(tn,Xn(1,:),'b','LineWidth',2); ha.Color(4) = alpha; 
    hb = plot(tn,Xn(2,:),'r','LineWidth',2); hb.Color(4) = alpha;
end
xlim([0,100]); ylim([0,100]); 
xlabel('time (sec)'); ylabel('copy numbers (molecules)');

hold on;

% CME Moment evolution
k = monomol.k;
a0 = monomol.X0(1);
b0 = monomol.X0(2);

A = [k(1)/k(2), a0 - k(1)/k(2), 0,0,0,0;
     k(1)/k(3), (k(2)*a0 -k(1))/(k(3)-k(2)), (b0 - (k(2)*a0 -k(1))/(k(3)-k(2)) -k(1)/k(3)),0,0,0;
     k(1)/k(2), a0 - k(1)/k(2),0,-a0,0,0;
     k(1)/k(3), (k(2)*a0 -k(1))/(k(3)-k(2)), (b0 - (k(2)*a0 -k(1))/(k(3)-k(2)) -k(1)/k(3)),...
    -a0*k(2)/(k(3)-k(2)), (1 - 2*k(2)/(k(3) + k(2)))*(a0*k(2)/(k(3)-k(2))) - b0, (2*a0*k(2)^2)/(k(3)^2 - k(2)^2);
     0,0,0, a0*k(2)/(k(3) - k(2)),0,-a0*k(2)/(k(3)-k(2))]

t = linspace(0,100,1000);
F = [ones(1,length(t));
     exp(-k(2)*t);
     exp(-k(3)*t);
     exp(-2*k(2)*t);
     exp(-2*k(3)*t);
     exp(-(k(3)+k(2))*t)];

% plot solution
M = A*F;
plot(t,M(1,:),'--k','LineWidth',2);
plot(t,M(2,:),'--k','LineWidth',2);
plot(t,M(1,:) + 2*sqrt(M(3,:)),':k','LineWidth',2);
plot(t,M(1,:) - 2*sqrt(M(3,:)),':k','LineWidth',2);
plot(t,M(2,:) + 2*sqrt(M(4,:)),':k','LineWidth',2);
plot(t,M(2,:) - 2*sqrt(M(4,:)),':k','LineWidth',2);

