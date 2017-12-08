%=============================================================
% This program prices a standard European call option 
% by both naive Monte Carlo simulation
% and importance sampling.
%
% NB: Remember to change values in f_prime function for now.
%
% Shawn Ban
% Sep 21, 2007
%
%=============================================================

clear all
close all
tic

% Declaration of variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10000; % number of iterations
s0 = 51; % Initial stock price
strike = 53; % Needs to be high.
r = 0.05; % Interest rate
sigma = 0.1; % Volatility. 
time = 1; %For simplicity sake, assume 1 year.
timesteps = 1000; % Also for simplicity, 1000 time steps.


naive_payoff = zeros(N,1); % Our payoff matrix
naive_stock_matrix = zeros(N,1); % Our stock matrix
IS_payoff = zeros(N,1); % IS payoff matrix
IS_stock_matrix = zeros(N,1); % IS stock matrix

d1 = (log(s0/strike) + (r+0.5*sigma*sigma)*time)/(sigma*sqrt(time));
d2 = d1 - (sigma*sqrt(time));

%true_value = s0*normcdf(d1,0,1) - strike*exp(-r*time)*normcdf(d2,0,1); %Black-Scholes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive Monte Carlo Case                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N 

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps))];

    W_final = W(timesteps+1); % Final value of Brownian motion

    % Final stock price:
    S_final = s0 * exp ((r - 0.5 * (sigma^2))*time + sigma * W_final);

    if S_final >= strike
         naive_payoff(i,1) = S_final - strike;
    else naive_payoff(i,1) = 0;
    end

naive_stock_matrix(i,1) = S_final;
    
end

naive_price = sum(naive_payoff)/N;
naive_error = std(naive_payoff)/sqrt(N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importance Sampling Case                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(0,2);
z= (log(strike/s0)-(r-0.5*sigma*sigma)*time)/(sigma*time);
x = x+z;
y = 0.5*(x.^2)*time - log(s0*exp((r-0.5*sigma*sigma)*time + x * time * sigma) - strike); 
plot(x,y)
figure
y_prime = x*time - ((s0*exp((r-0.5*sigma*sigma) + x*sigma*time)*time*sigma)/ (s0*exp((r-0.5*sigma*sigma) + x*sigma*time) - strike));
plot(x,y_prime)
% We have to minimize this function.



u_star = Bisection('f_prime',z,z+1);

for i = 1:N 

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ISW = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps) + u_star*(time/timesteps))];

    ISW_final = ISW(timesteps+1); % Final value of Brownian motion

    % Final stock price:
    IS_final = s0 * exp ((r - 0.5 * (sigma^2))*time + sigma * ISW_final);

    if IS_final >= strike
         IS_payoff(i,1) = (IS_final - strike) * exp (-u_star*ISW_final + 0.5*u_star*u_star*time);
    else IS_payoff(i,1) = 0;
    end

IS_stock_matrix(i,1) = IS_final;
end

IS_price = sum(IS_payoff)/N;
IS_error = std(IS_payoff)/sqrt(N);
disp(sprintf('u* is %2.9f',u_star));
%disp(sprintf('True value is                %12.9f',true_value));
disp(sprintf('Naive Monte Carlo price is   %12.9f',naive_price));
disp(sprintf('Importance sampling price is %12.9f',IS_price));
disp(sprintf('Naive Monte Carlo error is   %12.9f',naive_error));
disp(sprintf('Importance sampling error is %12.9f',IS_error));    
toc
