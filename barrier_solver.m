clear all
S = 50; %initial stock price
H = 53; %the barrier
T = 1; %time, set to 1 for now.
sigma = 0.2; %volatility
r = 0.05; %interest rate
K = 52; %strike price.
N = 200; %Number of points
t=linspace(0.01,0.99, N); %time, this shows up on the row index.
k=linspace(0.01,51.99, N); %target stock price, this shows up on the column index.

Y = zeros(N,N);
for i = 1:N;
    for j = 1:N;
        Y(i,j) = (t(i)/(2*sigma^2))*((log(H)-log(S))/t(i) - r + 0.5*sigma^2)^2 + ((T-t(i))/(2*sigma^2))*((log(k(j))-log(H))/(T-t(i)) - r + 0.5*sigma^2)^2 - log(K-k(j));
    end
end
[a,b] = min(Y);
[c,d] = min(a);
k_star = k(d)
t_star = t(b(d))