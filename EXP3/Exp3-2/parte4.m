clc;
clear;
close all;

load Analise_Sintese32.mat;
[violino, fa] = audioread("violino.wav");
violino = transpose(violino);
%%
n_bits = 7;


Qviolino = midtreadQ(violino, n_bits, 1); % Não tenho certeza se o 3 argumento é esse, mas está funcionando.

size_violino = length(violino);
N_canais = 32;

Qviolino_filtrado = zeros(N_canais,size_violino);
Qviolino_M_filtrado = zeros(N_canais,ceil(size_violino/32));
Qviolino_L_reconstruindo = zeros(N_canais,size_violino);
Qviolino_sintese = zeros(N_canais,size_violino);


for i = 1:1:N_canais;
    Qviolino_filtrado(i,:) = filter(PQMF32_Hfilters(i,:), 1, Qviolino);

    Qviolino_M_filtrado(i,:) = Qviolino_filtrado(i, 1:32:end);
    
    Qviolino_L_reconstruindo(i, 1:32:end) = Qviolino_M_filtrado(i, 1:1:end);
    
    Qviolino_sintese(i,:) = filter(PQMF32_Gfilters(i,:), 1, Qviolino_L_reconstruindo(i,:));
end


SNRs = zeros(32,1);
passo = ceil(size_violino/32);

for i =1:32
    SNRs(i) = snr(violino, Qviolino_sintese(i,:), 0); %Aqui comparamos com o violino original.
end

figure(4)
plot(SNRs)
title("Análise dos SNRs em cada um dos 32 canais")


Qviolino_final = sum(Qviolino_sintese);
SNR_QViolino = snr(violino, Qviolino_final, 0);

%sound(Qviolino_final);
%Não esqueça de ouvir. **Fica horroroso com 4 bits.**

%% B
Nh = length(PQMF32_Hfilters);

L_violino = lappedQ(violino, fa, PQMF32_Hfilters, PQMF32_Gfilters, 4);
L_violino = L_violino(Nh:end-Nh);

SNR_LViolino = snr(violino, L_violino,0);

%% C
Ladap_violino = lappedQadap(violino, fa, PQMF32_Hfilters, PQMF32_Gfilters, 4);
Ladap_violino = Ladap_violino(Nh:end-Nh);

SNR_LadapViolino = snr(violino, Ladap_violino,0); %Do nada ele fala de um fator. Que fator é esse ?



