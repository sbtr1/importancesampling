%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program prices a vanilla option by Monte Carlo simulation
% after using an importance sampling change of measure.
% For comparison with naive Monte Carlo case.
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

interest_rate = 0.05;

sigma = 0.01; % Volatility. Needs to be small. how small?

time = 1; %For simplicity sake, assume 1 year.

timesteps = 1000; % Also for simplicity, 1000 time steps.

strike = 52.7; % Needs to be high. How high?

payoff = zeros(N,1); % Our payoff matrix

stock_matrix = zeros(N,1); % Our stock matrix

u_star = (1/sigma*time)*(log(strike/s0)) + 0.5*sigma - (interest_rate/sigma); 

for i = 1:N %Simulate 10000 draws

% Time variable in case you want to plot:
% t = (0:1:number_timesteps)'/number_timesteps;

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    W = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps) + u_star*(time/timesteps))];

    W_final = W(1001); % Final value of Brownian motion

    % Final stock price:
    S_final = s0 * exp ((interest_rate - 0.5 * (sigma^2))*time + sigma * W_final);

    if S_final >= strike
         payoff(i,1) = 1 * exp (-u_star*W_final + 0.5*u_star*u_star*time);
    else payoff(i,1) = 0;
    end

stock_matrix(i,1) = S_final;
    
end

sum_payoff = sum(payoff);
price = sum(payoff)/N
error = std(payoff)/sqrt(N)
toc
