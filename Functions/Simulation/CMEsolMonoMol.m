function [P] = CMEsolMonoMol(k,a,b,t,a0,b0)
%% Explicit solution to the small mono-molecular chain
%   k1   k2   k3
% 0 -> A -> B -> 0
% Computes Pr(A(t) = a, B(t) = b | A(0) = a0, B(0) = b0)
%
% Inputs:
%     k - kinetic rate parameters
%     a - population of A molecules
%     b - population of B molecules
%     t - time since t0
%     a0 - initial population of A molecules
%     b0 - initial population of B molecules
% 
% Output:
%    the CME solution P(a,b,t | a0,b0)
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Queensland University of Technology

% pre-allocate arrays 
lp_poiss = zeros(length(0:a)*length(0:b),1);
lp_m_sum = zeros(length(0:a)*length(0:b),1);
% compute rate equation solutions lambda(t) and pa(t), pb(t)
lambda_a = (k(1)/k(2))*(1-exp(-k(2)*t));
lambda_b = (k(1)/(k(3)-k(2)))*((k(2)-2*k(3))*exp(-k(3)*t)/k(3) + exp(-k(2)*t)) +k(1)/k(3);
pa_a = exp(-k(2)*t);
pa_b = (exp(-k(2)*t) - exp(-k(3)*t))*(k(2)/(k(3) -k(2)));
pb_b = exp(-k(3)*t);
k = 1;
% apply the nested convolutions using logs

% compute convolution of product Poisson distribution with nested convolution
% of multinomial distributions

% special caase for stationary distribution
if isinf(t)
    P = exp(logProductPoisson(a,b,lambda_a,lambda_b));
    return;
end

for a_w = 0:a
    for b_w = 0:b
        b_z = max([0;b_w - b0]):min([b_w;max([0;a0-a_w])]);
        lp_mA = zeros(size(b_z))';
        lp_mB = zeros(size(b_z))';
        % compute convolution of multinomial distributions
        for i=1:length(b_z)
           lp_mA(i) = logMultinomialA(a_w,b_z(i),a0,pa_a,pa_b);
           lp_mB(i) = logMultinomialB(0,b_w-b_z(i),b0,pb_b);
        end
        lp_m_sum(k) = logexpsum([lp_mA+lp_mB]);
        lp_poiss(k) = logProductPoisson(a-a_w,b-b_w,lambda_a,lambda_b);
        k = k + 1;
    end
end
lP = logexpsum([lp_m_sum+lp_poiss]);
P = exp(lP);
if isnan(P)
    P = 0;
end


function [l_sum] = logexpsum(X)
%% computes log(sum_{i=0}^N exp(x_i)) while avoiding numerical
% overflow/underflow
if isempty(X)
    l_sum = -inf;
else
    %% numerically stable log sum of exponents
    M = max(X);
    if isinf(M)
        l_sum = -inf;
    else
        X = X - M;
        l_sum = log(sum(exp(X))) + M;
    end
end

function [lp] = logProductPoisson(a,b,lambda_a,lambda_b)
%% Computes the log probability of the product Poisson distribution
if a >=0 && b >= 0
    lp = -(abs(lambda_a) + abs(lambda_b)) + a*log(lambda_a) + b*log(lambda_b) ...
         - sum(log(1:a)) - sum(log(1:b));
else
    lp = -inf;
end

function [lp] = logMultinomialA(a,b,a0,pa_a,pa_b)
%% Computes the log probability of the Multinomial distribution for A
if a + b <= a0
    lp = sum(log((a0-a-b+1):a0)) + (a0 - a - b)*log(1 -abs(pa_a) - abs(pa_b)) ...
        + a*log(pa_a) + b*log(pa_b) - sum(log(1:a)) - sum(log(1:b));
else
    lp = -inf;
end

function [lp] = logMultinomialB(a,b,b0,pb_b)
%% Computes the log probability of the Multinomial distribution for B
if a == 0 && b <= b0
    lp = sum(log((b0-b+1):b0)) + (b0 - b)*log(1 - abs(pb_b)) ...
         + b*log(pb_b) - sum(log(1:b));
else
    lp = -inf;
end
