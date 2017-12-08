%=============================================================
% This program prices a knock in barrier put option 
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

N = 10000; % number of iterations
s0 = 50; % Initial stock price.
H = 54; 
K = 52; % Needs to be high.
r = 0.05; % Interest rate
sigma = 0.1; % Volatility. 
T = 1; %For simplicity sake, assume 1 year.
timesteps = 200; % Also for simplicity, 1000 time steps.

naive_payoff = zeros(N,1); % Our payoff matrix
IS1_payoff = zeros(N,1);
IS2_payoff = zeros(N,1);% Our payoff matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive Monte Carlo Case                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N 

% Simulation of Brownian motion from N(0, time/timesteps):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = zeros(timesteps+1,1);
    S(1) = s0;
    W = [0; cumsum(randn(timesteps,1).*sqrt(T/timesteps))];
    for j = 2:timesteps+1
        S(j) = s0 * exp ((r - 0.5 * (sigma^2))*((j-1)*(T/timesteps)) + sigma * W(j));
    end
    S_final = S(timesteps+1);
    
    if max(S) >= H    
        if K >= S_final
             naive_payoff(i,1) = K-S_final;
        else naive_payoff(i,1) = 0;
        end
    else
        naive_payoff(i,1) = 0;
    end 
end

naive_price = exp(-r*T)*(sum(naive_payoff)/N);
naive_error = std(naive_payoff)/sqrt(N);


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
k_star = k(d); %These are the values that minimize Y
t_star = t(b(d)); %These are the values that minimize Y
t_bar = T - t_star;


u_star1 = (1/(sigma*t_star))*(log(H/s0)) + 0.5*sigma - (r/sigma); 


% Importance Sampling scheme 1: fixed time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:N 

% Simulation of Brownian motion from time 0 to t_star:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ISW1 = [0; cumsum(randn(timesteps,1).*sqrt(t_star/timesteps) + u_star1*(t_star/timesteps))];
    S1 = zeros(2*timesteps+1,1);
    S1(1) = s0;
    for j = 2:timesteps+1
        S1(j) = s0 * exp ((r - 0.5 * (sigma^2))*(j-1)*(t_star/timesteps) + sigma * ISW1(j));
    end
    ISW1_final = ISW1(timesteps+1);
    S_star = S1(timesteps+1); %stock price at t_star
    
    %Now we simulate from t_star to T (call this t_bar for ease):
    
    u_star2 = (1/(sigma*t_bar))*(log(k_star/S_star)) + 0.5*sigma - (r/sigma);
    
    ISW2 = [0; cumsum(randn(timesteps,1).*sqrt(t_bar/timesteps) + u_star2*(t_bar/timesteps))];
    ISW2_diff = ISW2(timesteps+1) - ISW1_final;
    for j = 2:timesteps+1
        S1(timesteps+j) = S_star * exp ((r - 0.5 * (sigma^2))*(j-1)*(t_bar/timesteps) + sigma * ISW2(j));
    end
    IS1_final = S1(2*timesteps+1); %final stock price

    if max(S1) >= H
        if K >= IS1_final
             IS1_payoff(i,1) = (K - IS1_final)*exp (-u_star1*ISW1_final + 0.5*u_star1*u_star1*t_star)*exp (-u_star2*ISW2_diff + 0.5*u_star2*u_star2*t_bar)*exp(-r*T);
        else IS1_payoff(i,1) = 0;
        end
    else
        IS1_payoff(i,1) = 0;
    end 
    
    end
    
IS1_price = sum(IS1_payoff)/N;
IS1_error = std(IS1_payoff)/sqrt(N);

% Importance Sampling scheme 2: state dependence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for i = 1:N 

% Simulation of Brownian motion from time 0 to t_star:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S2 = s0;
    t_count = 0;
    while ((S2 < H) & (t_count < timesteps))
        ISW3 = randn(1).*sqrt(T/timesteps) + u_star1*(T/timesteps);
        S2 = S2 * exp ((r - 0.5 * (sigma^2))*(T/timesteps) + sigma * ISW3);
        t_count = t_count + 1;
    end
    if t_count >= timesteps
        IS2_final = S2;
        IS2_payoff(i,1) = 0;
    else
       t_hat = T - (t_count/timesteps) * T;
       u_star3 = (1/(sigma*t_hat))*(log(k_star/S2)) + 0.5*sigma - (r/sigma);
       ISW4 = [0; cumsum(randn(timesteps,1).*sqrt(t_hat/timesteps) + u_star3*(t_hat/timesteps))];
       IS2_final =  S2 * exp ((r - 0.5 * (sigma^2))*(t_hat) + sigma * ISW4(timesteps+1));
       if K >= IS2_final
            IS2_payoff(i,1) = K - IS2_final;
       else
            IS2_payoff(i,1) = 0;
       end
    end
end

IS2_price = sum(IS2_payoff)/N;
IS2_error = std(IS2_payoff)/sqrt(N);



disp('                                                  ');
disp('-----------------------------------------------');
disp(sprintf('k_star is   %12.9f',k_star));
disp(sprintf('t_star is %12.9f',t_star));
disp(sprintf('Naive Monte Carlo price is   %12.9f',naive_price));
disp(sprintf('Importance sampling 1 price is %12.9f',IS1_price));
disp(sprintf('Importance sampling 2 price is %12.9f',IS2_price));
disp(sprintf('Naive Monte Carlo error is   %12.9f',naive_error));
disp(sprintf('Importance sampling 1 error is %12.9f',IS1_error));
disp(sprintf('Importance sampling 2 error is %12.9f',IS2_error));
toc
