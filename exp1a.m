close all;;
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

J_min = sigma_s^2 - transpose(p)*wO;

H = freqz(wO);


%%2
lambda = eig(R);

mi_max = 2/max(lambda);

%%3
[W,erro]= lms(u, d, 2, 500, 0.03);

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTANDO CURVA DE NIVEL 
w0_passo = linspace(-.5, .5, 100);
w1_passo = linspace(-.5, .5, 100);
[W0, W1] = meshgrid(w0_passo, w1_passo);
J = zeros(100,100);

for i = 1:100
    for j = 1:100
        wij = [W0(i,j); W1(i,j)];
        J(i, j) = var(d) + 2*transpose(p)*wij + transpose(wij)*R*wij;
    end
end

%PLOTANDO
figure(3)
curvas = 0.5:.3:3; %curvas que serão plotadas
contour(W0, W1, J, curvas)
hold on 
plot(wO(1), wO(2), 'p', 'MarkerSize', 12) %wO calculado analiticamente


%Encontrando os valores de iteração
J_iterado = zeros(1,501);
for i = 1:1:501
    w = transpose(W(i,:));
    J_iterado(i) = var(d) + 2*transpose(p)*w + transpose(w)*R*w;
end 

%desenhando nas curvas de nivel
for i=1:50:501
    plot(W(i, 1), W(i, 2), '.', 'MarkerSize', 12);
end

%%PLOTANDO e²(n) em dB
figure(4);
plot(20*log10(J_iterado));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%e
MSE = mean(J_iterado);
EMSE = MSE - min(J_iterado);

%%
mi003 = J_iterado;
[a,mi001]= lms(u, d, 2, 500, 0.01);
[mi005]= lms(u, d, 2, 500, 0.05);

figure(5);
plot(mi001)
hold on
plot(mi003)
hold on
plot(mi005)
legend('\mu = 0.01', '\mu = 0.03', '\mu = 0.05');

%%
%não sei se faz sentido, quando eu pego d-y não dá x











