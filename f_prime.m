function c = f_prime(u)

time = 1;

r = 0.05;

sigma = 0.1;

s0 = 50;

strike = 53;


c = u*time - ((s0*exp((r-0.5*sigma*sigma) + u *sigma*time)*time*sigma)/ (s0*exp((r-0.5*sigma*sigma) + u*sigma*time) - strike));
