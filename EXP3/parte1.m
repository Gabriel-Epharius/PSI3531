clear all;
clc;
close all;

%% A

% Definição do sinal chirp
f_max = 4000; % A frequência deve variar de 0 a 4 kHz
t_c_max = 4; % Tempo contínuo com 4 segundos de duração
f_a = 8000; % A frequência de amostragem é de 8 kHz

A_0 = 2*pi*f_max/(2*t_c_max);
alpha_0 = A_0*(1/f_a)^2;

N = 4*f_a; % Número de amostras 
n = 0:1:(N-1);
x = zeros(1, N);
x = cos(alpha_0*(n.^2));

figure(1);
spectrogram(x, 512, 256, 1024, f_a);
title('Espectrograma do sinal chirp linear');

%% B

% Redução de taxa de amostragem por um fator M = 2
x_M = x(1:2:end);

% Aumentando a taxa de amostragem por um fator L = 2
% Inserindo zeros entre as amostras do sinal subamostrado
x_L = zeros(1, N);
x_L(1:2:end) = 2*x_M(1:1:end); %Multiplicando o sinal superamostrado por 2

% for i= 1:1:N/2
%     x_L(i*2-1) = 2*x_M(i); %Multiplicando o sinal superamostrado por 2
% end

figure(2);
spectrogram(x_M, 512, 256, 1024, f_a);
title('Espectrograma do sinal com taxa de amostragem reduzida');

figure(3);
spectrogram(x_L, 512, 256, 1024, f_a);
title('Espectrograma do sinal com taxa de amostragem aumentada');

% Comandos para ouvir os sinais
%sound(x, f_a);
%sound(x_M, f_a/2);
%sound(x_L, f_a);

%% C

% Projeto de filtro passa-baixas FIR com fase linear por trechos
f = [1960 2040]; % Frequência de corte aproximadamente igual a 2000 Hz
ripple = [2*10^(-4) 10^(-4)];
ganho = [1 0];
fs = f_a; % Frequência de amostragem do sinal superamostrado

[n,fo,ao,w] = firpmord(f,ganho,ripple,fs);
b = firpm(n,fo,ao,w);

% Gráficos do módulo e da fase do filtro projetado
figure(4);
freqz(b,1,1024,fs);
title('Resposta em frequência do filtro passa-baixas projetado');

%% D

% Filtrando o sinal superamostrado
x_L_filtrado = filter(b, 1, x_L);

figure(5);
spectrogram(x_L_filtrado, 512, 256, 1024);
title('Espectrograma do sinal filtrado e com taxa de amostragem aumentada');

figure(6);
plot(x_L_filtrado);
title('Gráfico do sinal filtrado no tempo');
xlabel('n');
ylabel('Sinal filtrado');

% Comando para ouvir o sinal filtrado
%sound(x_L_filtrado, f_a);

%% E

% É necessário filtrar o sinal x antes da redução da taxa de amostragem
% para evitar aliasign

x_filtrado = filter(b, 1, x);

figure(7);
spectrogram(x_filtrado, 512, 256, 1024);
title('Espectrograma do sinal filtrado');

% Redução de taxa de amostragem por um fator M = 2
x_M_filtrado = x_filtrado(1:2:end);

% Aumentando a taxa de amostragem por um fator L = 2
% Inserindo zeros entre as amostras do sinal subamostrado

x_L_reconstruindo = zeros(1, N);
x_L_reconstruindo(1:2:end) = 2*x_M_filtrado(1:1:end); %Multiplicando o sinal superamostrado por 2

x_sintese = filter(b, 1, x_L_reconstruindo);

figure(8);
spectrogram(x_sintese, 512, 256, 1024);
title('Espectrograma do sinal reconstruído');

