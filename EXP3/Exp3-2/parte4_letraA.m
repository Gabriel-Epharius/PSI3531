clc;
clear;
close all;

load Analise_Sintese32.mat;
ordem = 512-1;
[violino, fa] = audioread("violino.wav");
violino = transpose(violino);

%% A
n_bits = 4;

size_violino = length(violino);
N_canais = 32;

Qviolino_filtrado = zeros(N_canais,size_violino);
Qviolino_M_filtrado = zeros(N_canais,ceil(size_violino/32));
Qviolino_L_reconstruindo = zeros(N_canais,size_violino);
Qviolino_sintese = zeros(N_canais,size_violino);


for i = 1:1:N_canais;
    Qviolino_filtrado(i,:) = filter(PQMF32_Hfilters(i,:), 1, violino);

    Qviolino_M_filtrado(i,:) = Qviolino_filtrado(i, 1:32:end);
    
    Qviolino_L_reconstruindo(i, 1:32:end) = Qviolino_M_filtrado(i, 1:1:end);

    Qviolino_L_reconstruindo_quatizado = midtreadQ(Qviolino_L_reconstruindo, n_bits, 1);
    
    Qviolino_sintese(i,:) = filter(PQMF32_Gfilters(i,:), 1, Qviolino_L_reconstruindo(i,:));
end

SNRs = zeros(32,1);

figure(1);
for i =1:32
    SNRs(i) = snr(violino, Qviolino_sintese(i,:), ordem, 0); %Aqui comparamos com o violino original.
    hold on
end
title("Análise dos SNRs em cada um dos 32 canais")

figure(2);
Qviolino_final = sum(Qviolino_sintese);
SNR_QViolino = snr(violino, Qviolino_final, ordem, 0);

%sound(Qviolino_final);
%Não esqueça de ouvir.

%% B
Nh = length(PQMF32_Hfilters);

L_violino = lappedQ(violino, fa, PQMF32_Hfilters, PQMF32_Gfilters, 32);

SNR_LViolino = snr(violino, L_violino, ordem-1, 0); %checa pra mim se o atraso é isso mesmo, eu sempre me confundo.

%% C
Ladap_violino = lappedQadap(violino, fa, PQMF32_Hfilters, PQMF32_Gfilters, 4);

SNR_LadapViolino = snr(violino, Ladap_violino, ordem-1, 0); %aqui também.

%D

[Lpsico_violino, SMRs, vetor_Nbits] = lappedPsico(violino, fa, PQMF32_Hfilters, PQMF32_Gfilters, 192000);

SNR_Lpsico = snr(violino, Ladap_violino, ordem-1, 0); %aqui também ao quadrado.

%% E

%Agora precisa plotar os gráficos com SMRs e vetor_NBits mas ainda não sei
%como fazer.

figure(3)
surf(SMRs);
title("SMRs")
xlabel()

figure(4)
surf(vetor_Nbits);
title("Nbits")

