function c = f_prime_bar(u)

time = 0.339949749;
r = 0.05;
sigma = 0.1;
s0 = 50;
strike = 51;
c = u*time - ((s0*exp((r-0.5*sigma*sigma) + u *sigma*time)*time*sigma)/ (strike - s0*exp((r-0.5*sigma*sigma) + u*sigma*time));
