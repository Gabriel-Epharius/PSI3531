close all;
clear all;
clc;
%% PARTE 3 

load Analise_Sintese32.mat
%% A
% figure(1);
% for i = 1:32
%     freqz(PQMF32_Hfilters(i,:))
%     hold on
% end
% title('Resposta em frequ�ncia dos 32 filtros de an�lise');

%% B
f_max = 22050; % A frequência deve variar de 0 a 4 kHz
t_c_max = 20; % Tempo contínuo com 4 segundos de duração
f_a = 44100; % A frequência de amostragem é de 8 kHz

A_0 = 2*pi*f_max/(2*t_c_max);
alpha_0 = A_0*(1/f_a)^2;

N = t_c_max*f_a; % Número de amostras 
n = 0:1:(N-1);
x = cos(alpha_0*(n.^2));

figure(2);
spectrogram(x, 512, 256, 1024, f_a);
title('Espectrograma do sinal chirp linear');

N_canais = 32;

x_filtrado = zeros(N_canais,N);
x_M_filtrado = zeros(N_canais,ceil(N/32));
x_L_reconstruindo = zeros(N_canais,N);
x_sintese = zeros(N_canais,N);
y = zeros(1,N);

for i = 1:1:N_canais;
    x_filtrado(i,:) = filter(PQMF32_Hfilters(i,:), 1, x);

    x_M_filtrado(i,:) = x_filtrado(i, 1:32:end);
    
    x_L_reconstruindo(i, 1:32:end) = x_M_filtrado(i, 1:1:end);
    
    x_sintese(i,:) = filter(PQMF32_Gfilters(i,:), 1, x_L_reconstruindo(i,:));
end

y = sum(x_sintese);

figure(3);
spectrogram(y, 512, 256, 1024, f_a);
title('Espectograma do sinal chirp reconstru�do');

ordem = 512-1;
diff = y(ordem+1:end) - x(1:end-ordem); 

figure(4);
plot(diff);
title('Diferen�a entre o sinal reconstru�do e o sinal original')

%% C

NR_PQMF=snr(x, y,0); % Qual � desse atraso que orienta o roteiro ? Precisa mudar alguma coisa ? 

%% D

violino = transpose(audioread("violino.wav"));
size_violino = length(violino);

%Primeiramente n�o vou quantizar, depois � s� aplicar aquela fun��o do
%Vitor.

violino_filtrado = zeros(N_canais,size_violino);
violino_M_filtrado = zeros(N_canais,ceil(size_violino/32));
violino_L_reconstruindo = zeros(N_canais,size_violino);
violino_sintese = zeros(N_canais,size_violino);

for i = 1:1:N_canais;
    violino_filtrado(i,:) = filter(PQMF32_Hfilters(i,:), 1, violino);

    violino_M_filtrado(i,:) = violino_filtrado(i, 1:32:end);
    
    violino_L_reconstruindo(i, 1:32:end) = violino_M_filtrado(i, 1:1:end);
    
    violino_sintese(i,:) = filter(PQMF32_Gfilters(i,:), 1, violino_L_reconstruindo(i,:));
end


SNRs = zeros(32,1);
passo = ceil(size_violino/32);

for i =1:32
    SNRs(i) = snr(violino, violino_sintese(i,:), 0);
end

figure(4)
plot(SNRs)
title("An�lise dos SNRs em cada um dos 32 canais")

figure(5)
plot(violino_sintese(1,:));
title("Gr�fico do canal 0");

figure(6)
plot(violino_sintese(6,:));
title("Gr�fico do canal 5");

figure(7)
plot(violino_sintese(11,:));
title("Gr�fico do canal 10");


violino_final = sum(violino_sintese);
SNR_violino = snr(violino, violino_final, 0);

