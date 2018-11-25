%% Demonstration of stationary distribution computation
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology


% initialise random number generator for reproducibility
rng(502,'twister');
h = figure;
tic
% Build mono-molecular chain
[monomol] = MonoMolecularChain([1.0;0.1;0.05],[100;0]);

% generate N long time realisations
N = 1000;
X_r = cell(N,1); t_r = cell(N,1);
A = zeros(N,1);
B = zeros(N,1);
for i=1:N
    [X_r{i},t_r{i}] = GillespieDirectMethod(monomol,1000);
    A(i) = X_r{i}(1,end);
    B(i) = X_r{i}(2,end);
end

% plot histograms of A and B 
hold on;
edges = [0,((0:max(B)) +0.5),max(B)+1]
histogram(A,edges,'Normalization','probability');
histogram(B,edges,'Normalization','probability');

% Actual Mean and Variances from CME derivation
Ma = monomol.k(1)/monomol.k(2);
Mb = monomol.k(1)/monomol.k(3);
Va = Ma; Vb = Mb;
a = 0:(max(B)+1);
b = a;

% plot Gaussian approximations
plot(a,normpdf(a,Ma,sqrt(Va)),'--b','LineWidth',2);
plot(b,normpdf(b,Mb,sqrt(Vb)),'--r','LineWidth',2);
xlabel('copy numbers (molecules)');
ylabel('stationary probability');

%% compute exact marginals in A and B from the full CME
a =0:75
b = 0:75
[A,B] = ndgrid(a,b);
P = zeros(size(A));
Pa = zeros(size(a));
Pb = zeros(size(b));
for i=1:length(a)
    for j=1:length(b);
        P(i,j) = CMEsolMonoMol(monomol.k,A(i,j),B(i,j),inf,100,0);
    end
end
Pa = sum(P,2)';
Pb =  sum(P,1);

% plot marginals
plot(a,Pa,':b','LineWidth',2);
plot(b,Pb,':r','LineWidth',2);

% compute transient and stationary full CME solution
h2 = cell(4,1)
hi = 1;
for T = [5,10,20,60,inf]
    a =0:75
    b = 0:75
    [A,B] = ndgrid(a,b);
    P = zeros(size(A));
    Pa = zeros(size(a));
    Pb = zeros(size(b));
    for i=1:length(a)
        for j=1:length(b);
            P(i,j) = CMEsolMonoMol(monomol.k,A(i,j),B(i,j),T,100,0);
        end
    end
   
    % plot full CME stationary distribution
    h2{hi}=figure;
    contourf(A,B,P,4);
    xlabel('A');
    ylabel('B');
    txt = {'P(A,B,t)'};
    text(80,-04,txt)
    view(2);
    xlim([0,max(a)]);
    ylim([0,max(b)])
    h3 = colorbar
    hi = hi+1
end
toc
