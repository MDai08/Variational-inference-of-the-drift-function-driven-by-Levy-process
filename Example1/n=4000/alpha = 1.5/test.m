%% Explicit form of the estimator
%  kernel matrix G, block-diagonal weight matrix A
%  Q=-inv(I+CGA)*C*F (no Levy)
%  Q=-inv(I+CGA)*C*(F+Z) (with Levy) 
%  F=1/2*sum(tr[D(z)K(z,z')'']|z'=z_j)

%%  initial data
clear,clc;

V = load('data.mat');           %% samples with levy

y = V.Y';                    
[m,n] = size(y);
C = 1;                          %% parameter of weight                           
alpha = 1.5;                   %% parameter of Levy process
l = 1.2;                        %% superparameter of kernel
%r = @(y)    4*y-4*y.^3;         %% drift term of known
d = @(y)     1;                 %% diffusion term
A = eye(m);                     %% block-diagonal weight matrix
I = eye(m);

%% calculate vector Q
G = zeros(m);                   %% kernel matrix
for i = 1 : m
    for j = 1 : m
        G(i,j) = exp(-(y(i)-y(j))^2/(2*l^2))/(l^2)-(exp(-(y(i)-y(j))^2/(2*l^2))*((y(i)-y(j))^2))/(l^4);
    end
end

F = zeros(m,1);                 %% 求梯度F(公式中的梯度y)
for i = 1 : m
    for j = 1 : m
        F(i) = F(i) + exp(-(y(i)-y(j))^2/(2*l^2))*((y(i)-y(j))/(l^4))+exp(-(y(i)-y(j))^2/(2*l^2))*(2*(y(i)-y(j))/(l^4))-exp(-(y(i)-y(j))^2/(2*l^2))*((y(i)-y(j))^3/(l^6));
    end
end

Z = zeros(m,1);
cons = alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));  %the constant of alpha-stable symmetric jump measure
% f = @(y) cons*(abs(y).^(-1-alpha)).*(exp(-(y(i)-(y(j)+y)).^2/(2*l.^2))-exp(-(X(i)-X(j)).^2/(2*l.^2)));
for i=1:m
    for j=1:m       
        g = @(z) cons*(abs(z).^(-1-alpha)).*(-((y(i)-(y(j)+z))/(l.^2)).*exp(-(y(i)-(y(j)+z)).^2/(2*l.^2))+((y(i)-y(j))/(l.^2)).*exp(-(y(i)-y(j)).^2/(2*l.^2)));     
        Z(i) = Z(i)+integral(g,-Inf,-0.001)+integral(g,0.001,Inf);
    end
end

% R = zeros(m,1);
% for i = 1 : m
%     R(i) = r(y(i));
% end

% Q = -inv(I+C*G*A)*C*(G*R+F/2+Z);

Q = -inv(I+C*G*A)*C*(F/2+Z);
%% Calculate the function P
step = 3/(m-1);
x = -1.5:step:1.5;  

L = zeros(m,1);         %求L(公式中的y)
for i=1:m
    for j=1:m
        L(i) = L(i)+(exp(-(x(i) - y(j))^2/(2*l^2))*(x(i) - y(j))^2)/(l^4) - exp(-(x(i) - y(j))^2/(2*l^2))/l^2;       
    end
end

J = zeros(m,1);         %求积分项
for i=1:m
    for j=1:m
       f = @(z) cons*(abs(z).^(-1-alpha)).*(exp(-(x(i)-(y(j)+z)).^2/(2*l.^2))-exp(-(x(i)-y(j)).^2/(2*l.^2)));
       J(i) = J(i)+integral(f,-Inf,-0.001)+integral(f,0.001,Inf);
    end
end

K = zeros(m);           %梯度kernel
for i=1:m
    for j=1:m
        K(i,j) = ((x(i)-y(j))/(l^2))*(exp(-(x(i)-y(j)).^2/(2*l.^2)));
    end
end

P = zeros(m,1);   %求\psi
for i=1:m
%     P(i) = -C*(K(i,:)*A*Q+(K(i,:)*R+L(i)/2+J(i)));
     P(i) = -C*(K(i,:)*A*Q+(L(i)/2+J(i)));
end

% PF = zeros(m-2,1);
% RPF = zeros(m-2,1);
% RRPF = zeros(m-2,1);
% for i=1:m-2
%     PF(i)=(P(i+2)-P(i))/(2*step);
%     RPF(i)=(x(i+1)-x(i+1)^3); 
%     RPF(i)=(-1)*(x(i+1)^2+2*x(i+1)+1); 
% end

% for i = 1 : m
%     RRPF(i) = (-1)*(y(i)^2+2*y(i)+1); 
% end

PF = zeros(m-2,1);
for i=1:m-2
    PF(i)=(P(i+2)-P(i))/(2*step);
    RPF(i)=x(i+1)-x(i+1)^3;
end






% d = (1./(2.*(m))).*(sum(RRPF.^2))
% d1 = (1./(2.*(m))).*(sum(Q.^2))

%% 误差值
E(1)=0;
for i=1:m-2
    E(i+1)=E(i)+(PF(i)-RPF(i)).^2;
end
E = E(m-1)/(m-2);


plot(x(2:end-1),PF,'r-',x(2:end-1),RPF,'b-')


save drift120 m x y G F Q L K P PF RPF E
