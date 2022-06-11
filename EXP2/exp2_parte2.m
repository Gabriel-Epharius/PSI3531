close all;
clear;
clc;

% Definições
L = 120;
N = L/2;
Q = 512;
K = 2;
gama = 1;
[sinal, fs] = audioread('antarctica.wav'); 

%% 1
% Cortando sinal em trechos de 240 amostras (30 ms)
% trechos = zeros(L,46);
% 
% for i = 1:46
%         %trechos(:,46) = sinal(240*(i-1)+1:240*i);
% end

%% 2
% Base com 512 funções aleatórias
fnc_base = randn(N, Q);

%% 3
% Condição inicial do filtro de trato vocal
zs = zeros(1, 10); 

%% 4

gerado = zeros(size(sinal));

for i = 1:138
    % A
    % Determinando os coeficientes LPC do quadro atual
    trecho = sinal(L*(ceil(i/2)-1)+1:L*ceil(i/2));
    [aq, sig]= lpc(trecho.*hamming(L), 10);
    
    subquadro = trecho(1+ L/2:2*L/2);

    % C
    % Filtrando as 512 sequ"Encias da base pelo filtro
    fnc_base_filt = zeros(N, Q);
    fnc_base_filt = filter(1, aq, fnc_base);

    % D
    [y0, zs] = filter(1, aq, zeros(N,1));

    % E
    e0 = subquadro - y0;

    % F
    [ganhos, indices] = find_Nbest_components(e0, fnc_base_filt, K);
    
    % G
    d = fnc_base_filt(:,indices)*ganhos;
    
    % H
    [gerado_trecho, zs] = filter(1, aq, d, zs);   
    
    gerado(N*(i-1)+1:N*i) = d;

end

figure(1);
plot(sinal);
hold on
plot(gerado);
legend("Sinal original", "CELP");
title("Comparação do sinal original com o reconstruído (codificador CELP)");

audiowrite("antarctida_CELP.wav", gerado, 400);







