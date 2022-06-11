close all;
clear;
clc;

%% Definições
N = 500;

n = [0:N-1]; % Vetor de amostras 
s = randn(1, N)*sqrt(.01); % Ruído branco de média nula e variância 0.01

fi = rand*2*pi; % Variável aleatória distribuida uniformemente entre 0 e 2pi

x = sin((2*pi*n/10) + (pi/6) + fi); % Interferência
u = 5*sin((2*pi*n/10) + fi); % Sinal correlacionado com a interferência
d = s + x; % Resposta desejada para o filtro adaptativo

%% A

xcorr(u, 1, 'biased'); % verificação matriz R
xcorr(d, u, 1, 'biased'); % verificaão matriz p

%%% 

R = [ 12.5, 10.1127; 10.1127, 12.5];
p = [2.1651; 1.0168];

wO = R\p;

J_min = var(d) - transpose(wO)*p;

figure(1);
freqz(wO);
title('Resposta em frequência do filtro adaptativo');

%% B
lambda = eig(R);

mi_max = 2/max(lambda);

%% C
[W,erro]= lms(u, d, 2, N, 0.03);

figure(2);
plot(n, erro, 'red');
hold on ;
plot(n, u, 'blue');
hold on;
plot(n, s, 'green');
legend('e(n)', 'u(n)', 's(n)');
title('Sinais e(n), u(n) e s(n)');
xlabel('n');

figure(3);
plot(W);
hold on;
plot(ones(N)*wO(1));
hold on;
plot(ones(N)*wO(2));
legend('w0 iterado','w1 iterado');
title('Iteração dos coeficientes');
xlabel('Número de iterações');
ylabel('Valor dos coeficientes');

% Curvas de nível
w0_passo = linspace(-.5, .5, 100);
w1_passo = linspace(-.5, .5, 100);
[W0, W1] = meshgrid(w0_passo, w1_passo);
J = zeros(100,100);

for i = 1:100
    for j = 1:100
        wij = [W0(i,j); W1(i,j)];
        J(i, j) = var(d) - 2*transpose(p)*wij + transpose(wij)*R*wij;
    end
end

figure(4)
contour(W0, W1, J);
hold on ;
plot(wO(1), wO(2), 'p', 'MarkerSize', 12) % Coeficientes ótimos para comparação

% Trajetória dos coeficientes

plot(W(:, 1), W(:, 2), '.-', 'MarkerSize', 5);
title('Curva de nível do LMS');
xlabel('W0');
ylabel('W1');

% Gráfico de e²(n) em dB

figure(5);
plot(10*log10(erro.^2));
title('Erro quadrático');
xlabel('n');
ylabel('Erro (dB)');

%% D

mi_max_exp = 0;
limite = 0.05;
for mi = 0.03:0.0001:0.08
   if(mi_max_exp == 0)
       [W_temp, erro_temp] = lms(u, d, 2, N, mi);
       if ((W_temp(N+1, 1)-wO(1)>limite) || (W_temp(N+1, 1)-wO(1)<-limite) || (W_temp(N+1, 2)-wO(2)>limite) || (W_temp(N+1, 2)-wO(2)<-limite))
           mi_max_exp = mi;
       end
   end
end

%% E
J003 = zeros(1, N);
for i = 1:N
    %Para cada realização, um novo valor de fi e um novo valor de s
    s_temp = randn(1, N)*sqrt(.01); %ruído branco
    fi_temp = rand*2*pi;
    x_temp = sin((2*pi*n/10) + (pi/6) + fi_temp);
    u_temp = 5*sin((2*pi*n/10) + fi_temp);
    d_temp = s_temp + x_temp;

    [W_temp, erro_temp] = lms(u_temp, d_temp, 2, N, 0.03);
    
    J003 = J003 + erro_temp.^2;
end
J003 = J003/N;

figure(6);
plot(10*log10(J003));
MSE = mean(J003((N-100):N));
MSE_vec = ones(1,500)*10*log10(MSE);
hold on;
plot(MSE_vec);
legend('J(n)','MSE');
title('Determinação gráfica do MSE')
xlabel('n');
ylabel('J(n) (dB) ');

%%
J001 = zeros(1, N);
J005 = zeros(1, N);
for i = 1:N
    %Para cada realização, um novo valor de fi e um novo valor de s
    s_temp = randn(1, N)*sqrt(.01); %ruído branco
    fi_temp = rand*2*pi;
    x_temp = sin((2*pi*n/10) + (pi/6) + fi_temp);
    u_temp = 5*sin((2*pi*n/10) + fi_temp);
    d_temp = s_temp + x_temp;

    [W_temp001, erro_temp001] = lms(u_temp, d_temp, 2, N, 0.01);
    [W_temp005, erro_temp005] = lms(u_temp, d_temp, 2, N, 0.05);
    
    J001 = J001 + erro_temp001.^2;
    J005 = J005 + erro_temp005.^2;
end
J001 = J001/N;
J005 = J005/N;

figure(7);
plot(10*log10(J001))
hold on
plot(10*log10(J003))
hold on
plot(10*log10(J005))
legend('\mu = 0.01', '\mu = 0.03', '\mu = 0.05');
title('Comparação entre diferentes valores de \mu')
xlabel('n');
ylabel('J(n) (dB) ');