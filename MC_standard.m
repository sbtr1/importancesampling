%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program prices a standard European call option 
% by naive Monte Carlo simulation. 
% We'll compare results with our importance
% sampling algorithm.
%
% Shawn Ban
% Sep 21, 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
tic

% Declaration of variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10000; % number of iterations

s0 = 50; % Initial stock price

r = 0.05;

sigma = 0.01; % Volatility. Needs to be small. how small?

time = 1; %For simplicity sake, assume 1 year.

timesteps = 1000; % Also for simplicity, 1000 time steps.

strike = 52; % Needs to be high. How high?

payoff = zeros(N,1); % Our payoff matrix

stock_matrix = zeros(N,1); % Our stock matrix

d1 = (log(s0/strike) + ((r+sigma*sigma)/2)*time)/sigma*sqrt(time);

d2 = d1 - sigma*sqrt(time);

true_value = s0*normcdf(d1,0,1) - strike*exp(-r*time)*normcdf(d2,0,1); %Black-Scholes

for i = 1:N %Simulate 10000 draws

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps))];

    W_final = W(timesteps+1); % Final value of Brownian motion

    % Final stock price:
    S_final = s0 * exp ((r - 0.5 * (sigma^2))*time + sigma * W_final);

    if S_final >= strike
         payoff(i,1) = S_final - strike;
    else payoff(i,1) = 0;
    end

stock_matrix(i,1) = S_final;
    
end

sum_payoff = sum(payoff);
price = sum(payoff)/N
error = std(payoff)/sqrt(N)
true_value
toc
