close all;
clear all;
clc;


%%1
N = 500;

n = [0:N-1]; % vetor de amostras 
s = randn(1, N)*sqrt(.01); %ruído branco

fi = rand*2*pi;

x = sin((2*pi*n/10) + (pi/6) + fi);
u = 5*sin((2*pi*n/10) + fi);
d = s + x;

%%%%%%%%%%%%%%%%%

xcorr(u, 1, 'biased'); % verificação matriz R (em -1, 0 e em 1)
xcorr(d, u, 1, 'biased'); % verificaão matriz p  -> d(n-1)u(n), d(n)u(n), d(n+1)u(n)

%%% 

R = [ 12.5, 10.1127; 10.1127, 12.5];
p = [2.1651; 1.0168];

wO = R\p;

sigma_s = sqrt(var(s));

J_min = sigma_s - transpose(p)*wO;

H = freqz(wO);


%%2
lambda = eig(R);

mi_max = 2/max(lambda);

%%3
[W,erro]=lms(u, d, 2, 500, 0.03);

figure(1);
title('Comparação')
plot(n, erro, 'red');
hold on ;
plot(u, 'blue');
hold on;
plot(s, 'green');

legend('erro', 'u', 's');


figure(2);
plot(W);
hold on;
plot(ones(N)*wO(1));
hold on;
plot(ones(N)*wO(2));

legend('w0 iterado','w1 iterado');

w = -5:0.01:5;
J = zeros(N,N);
for i = 1:N
    for j = 1:N
        wij = [w(i); w(j)];
        J(i, j) = var(d) + 2*transpose(p)*w + transpose(w)*R*w;
    end
end

plot(J)




