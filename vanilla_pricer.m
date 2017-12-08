x%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program prices a vanilla call option 
% by both naive Monte Carlo simulation
% and Monte Carlo simulation with importance sampling.
%
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

N = input('Enter number of iterations: '); 
s0 = input('Enter initial stock price: '); 
strike = input('Enter strike price: ');
r = input('Enter interest rate: ');
sigma = input('Enter volatility: ');
time = 1; %For simplicity sake, assume 1 year.
timesteps = 1000; % Also for simplicity, 1000 time steps.



naive_payoff = zeros(N,1); % Our payoff matrix
naive_stock_matrix = zeros(N,1); % Our stock matrix
IS_payoff = zeros(N,1); % IS payoff matrix
IS_stock_matrix = zeros(N,1); % IS stock matrix
u_star = (1/(sigma*time))*(log(strike/s0)) + 0.5*sigma - (r/sigma); 

%true_value = exp(-r*time)*(1-normcdf((log(strike/s0)-((r-0.5*sigma*sigma)*time))/(sigma*sqrt(time)),0,1));

% Naive Monte Carlo pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:N %Simulate 10000 draws

    % Simulation of Brownian motion from N(0, time/timesteps):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        W = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps))];

        W_final = W(timesteps+1); % Final value of Brownian motion

        % Final stock price:
        S_final = s0 * exp ((r - 0.5 * (sigma^2))*time + sigma * W_final);

        if S_final >= strike
             naive_payoff(i,1) = 1*exp(-r*time);
        else naive_payoff(i,1) = 0;
        end

    naive_stock_matrix(i,1) = S_final;
    
    end

    naive_price = sum(naive_payoff)/N;
    naive_error = std(naive_payoff)/sqrt(N);
    
    
% Importance Sampling pricer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:N %Simulate 10000 draws

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ISW = [0; cumsum(randn(timesteps,1).*sqrt(time/timesteps) + u_star*(time/timesteps))];

    ISW_final = ISW(timesteps+1); % Final value of Brownian motion

    % Final stock price:
    IS_final = s0 * exp ((r - 0.5 * (sigma^2))*time + sigma * ISW_final);

    if IS_final >= strike
         IS_payoff(i,1) = 1 * exp (-u_star*ISW_final + 0.5*u_star*u_star*time)*exp(-r*time);
    else IS_payoff(i,1) = 0;
    end

    IS_stock_matrix(i,1) = IS_final;
    
    end
    
    IS_price = sum(IS_payoff)/N;
    IS_error = std(IS_payoff)/sqrt(N);
    
    
%disp(sprintf('True value is                %12.9f',true_value));
disp(sprintf('Naive Monte Carlo price is   %12.9f',naive_price));
disp(sprintf('Importance sampling price is %12.9f',IS_price));
disp(sprintf('Naive Monte Carlo error is   %12.9f',naive_error));
disp(sprintf('Importance sampling error is %12.9f',IS_error));
toc