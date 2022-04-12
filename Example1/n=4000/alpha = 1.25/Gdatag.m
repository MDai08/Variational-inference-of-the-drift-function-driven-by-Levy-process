%% Generate sample paths for SDE with Brownian motion
% dX_t = g(X_t)dt +sigma*dB_t
% dX_t = f(X_t)dt +sigma*dB_t
% g(X_t) = 4*x-4*x.^3-beta*(x.^2+2*x+1)
% f(X_t) = 4*x-4*x.^3
% beta = 2;
tic;
clear;
randn('state',100);
rand('state',100);
format long;

tic;
alpha = 1.25;
beta = 0;
gam = 1;
delta = 0;

T = 200; 
b = 1;
[n,N] = size(beta);
sigma = @(x)                    1; 
f = @(b,x)                   x-x.^3; 
h1 = @(x)                       1;
% f = @(beta,x)           beta*(x.^2+2*x+1); 
% f = @(x)    x-x^3;
M = 2000000;                                                                   %the number of sample path
dt=T/M;                                                                     %the time step
t=0:dt:T;

X = zeros(M+1,N);
X(1,:) = -1;                                                                   %initial point  
for j = 1 : N
    for i = 1 : M
        Z1 = sqrt(dt)*randn;
        Z2 = stblrnd(alpha,beta,gam,delta);
        X(i+1,j) = X(i,j) + f(beta(j),X(i,j))*dt + sigma(X(i,j))*Z1+ h1(X(i))*(dt)^(1/alpha)*Z2; %alpha !=1
    end
end

N = 4000;
st = M/N;
for i = 1 : N+1
    Y(i) = X(1+(i-1)*st);
end
s = 0 : T/N : T;
plot(s,Y,'r*')

toc
save  data  Y

