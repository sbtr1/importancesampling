%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program prices a vanilla option by naive Monte
% Carlo simulation. 
% We'll compare results with our importance
% sampling algorithm.
%
% Shawn Ban
% Sep 11, 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
tic

% Declaration of variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 200000; % number of iterations

s0 = 48; % Initial stock price

r = 0.05;

sigma = 0.01; % Volatility. Needs to be small. how small?

time = 1; %For simplicity sake, assume 1 year.

timesteps = 1000; % Also for simplicity, 1000 time steps.

strike = 52.7; % Needs to be high. How high?

payoff = zeros(N,1); % Our payoff matrix

stock_matrix = zeros(N,1); % Our stock matrix

true_value = exp(-r*time)*(1-normcdf((log(strike/s0)-((r-0.5*sigma*sigma)*time))/(sigma*sqrt(time))  ,0,1));

for i = 1:N %Simulate 10000 draws

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps))];

    W_final = W(timesteps+1); % Final value of Brownian motion

    % Final stock price:
    S_final = s0 * exp ((r - 0.5 * (sigma^2))*time + sigma * W_final);

    if S_final >= strike
         payoff(i,1) = 1;
    else payoff(i,1) = 0;
    end

stock_matrix(i,1) = S_final;
    
end

sum_payoff = sum(payoff);
price = sum(payoff)/N
error = std(payoff)/sqrt(N)
true_value
toc
