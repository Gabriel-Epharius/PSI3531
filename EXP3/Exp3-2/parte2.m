close all;
clear all;
clc;
%% A
load h0_QMF.mat

figure(1);
freqz(h0);
title('Resposta em frequência do filtro h0(n)');


figure(2);
h1 = h0.*real(exp(1i*pi*(0:length(h0)-1)));
freqz(h1);
title('Resposta em frequência do filtro h1(n)');

figure(3);
g0 = 2*h0;
freqz(g0);
title('Resposta em frequência do filtro g0(n)');

figure(4);
g1 = -2*h1;
freqz(g1);
title('Resposta em frequência do filtro g1(n)');

%% B
f_max = 4000; % A frequência deve variar de 0 a 4 kHz
t_c_max = 4; % Tempo contínuo com 4 segundos de duração
f_a = 8000; % A frequência de amostragem é de 8 kHz

A_0 = 2*pi*f_max/(2*t_c_max);
alpha_0 = A_0*(1/f_a)^2;

N = 4*f_a; % Número de amostras 
n = 0:1:(N-1);
x = zeros(1, N);
x = cos(alpha_0*(n.^2));

x_filtrado = filter(h0, 1, x);

x_M_filtrado = x_filtrado(1:2:end);

x_L_reconstruindo = zeros(1, N);
x_L_reconstruindo(1:2:end) = x_M_filtrado(1:1:end); %Multiplicando o sinal superamostrado por 2

x_sintese_PB = filter(g0, 1, x_L_reconstruindo);

x_filtrado_PA = filter(h1, 1, x);

x_M_filtrado_PA = x_filtrado_PA(1:2:end);

x_L_reconstruindo_PA = zeros(1, N);
x_L_reconstruindo_PA(1:2:end) = x_M_filtrado_PA(1:1:end);

x_sintese_PA = filter(g1, 1, x_L_reconstruindo_PA);

y = x_sintese_PA + x_sintese_PB;

n = length(h0)-1;
diff = y(n+1:end) - x(1:end-n); 

figure(5);
spectrogram(y, 512, 256, 1024);
title('Espectrograma do sinal reconstruído');

figure(6);
plot(diff);
title('Diferença entre o sinal reconstruído e o sinal original')



