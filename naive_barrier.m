%=============================================================
% This program prices a knock in barrier call option 
% by naive Monte Carlo simulation and importance sampling.
%
% I'm trying to find which of the two cases: (i) fixed time; or
% (ii) state dependence gives a better IS scheme.
%
% Shawn Ban
% Jan 28, 2007
%
%=============================================================

clear all
close all
tic

% Declaration of variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 5000; % number of iterations
s0 = 50; % Initial stock price.
H = 52; 
K = 51; % Needs to be high.
r = 0.05; % Interest rate
sigma = 0.1; % Volatility. 
T = 1; %For simplicity sake, assume 1 year.
timesteps = 500; % Also for simplicity, 1000 time steps.


naive_payoff = zeros(N,1); % Our payoff matrix



%true_value = s0*normcdf(d1,0,1) - strike*exp(-r*time)*normcdf(d2,0,1); %Black-Scholes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive Monte Carlo Case                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N 

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = zeros(timesteps,1);
    W = [0; cumsum(randn(timesteps,1).*sqrt(T/timesteps))];
    for j = 2:timesteps+1
        S(j) = s0 * exp ((r - 0.5 * (sigma^2))*T + sigma * W(j));
    end
    if max(S) >= H
        W_final = W(timesteps+1); % Final value of Brownian motion
        S_final = S(timesteps+1);
    
        if S_final >= K
             naive_payoff(i,1) = S_final - K;
        else naive_payoff(i,1) = 0;
        end
    else
        naive_payoff(i,1) = 0;
    end

    
end

naive_price = exp(-r*T)*(sum(naive_payoff)/N)
naive_error = std(naive_payoff)/sqrt(N)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importance sampling setup                 % 
%                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = 200; %Number of points in matrix
t=linspace(0.01,T-0.01, J); %time, this shows up on the row index.
k=linspace(0.01,K-0.01, J); %target stock price, this shows up on the column index.

Y = zeros(J,J);
for i = 1:J;
    for j = 1:J;
        Y(i,j) = (t(i)/(2*sigma^2))*((log(H)-log(s0))/t(i) - r + 0.5*sigma^2)^2 + ((T-t(i))/(2*sigma^2))*((log(k(j))-log(H))/(T-t(i)) - r + 0.5*sigma^2)^2 - log(K-k(j));
    end
end
[a,b] = min(Y);
[c,d] = min(a);
k_star = k(d)
t_star = t(b(d))
toc
